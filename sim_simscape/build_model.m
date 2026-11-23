function build_model(S)
%BUILD_MODEL Create three_leg_switched.slx: Simscape switched circuit plus Simulink controller.
%   S = jsondecode of scenarios.json. Circuit values are fixed at build time; the run-specific
%   values (initial state, AC phase, controller data) come from workspace variables C, Ts,
%   Tsw, Tstop and fine_dec set by run_case.
if nargin < 1
    S = jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')), 'scenarios.json')));
end
p = S.params; sw = S.switch;
L = lib_paths();
load_system('ee_lib'); load_system('fl_lib'); load_system('nesl_utility');
mdl = 'three_leg_switched';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
set_param(mdl, 'Solver', 'FixedStepDiscrete', 'FixedStep', 'Ts', 'StopTime', 'Tstop', ...
    'SolverMode', 'SingleTasking', 'SingleTaskRateTransMsg', 'none', ...
    'ReturnWorkspaceOutputs', 'on', 'SignalLogging', 'off', 'SaveTime', 'off', ...
    'SaveOutput', 'off');
P = @(h) get_param(h, 'PortHandles');
num = @(x) num2str(x, 17);
blk = @(lib, name) add_block(lib, [mdl '/' name]);

% ---------------------------------------------------------------- Simscape circuit
sc = blk(L.solver, 'Solver');
set_param(sc, 'UseLocalSolver', 'on', 'LocalSolverChoice', 'NE_PARTITIONING_ADVANCER', ...
    'LocalSolverSampleTime', 'Ts', 'DoFixedCost', 'on');
gnd = blk(L.ref, 'Gnd');
cdc = blk(L.cap, 'Cdc');
set_param(cdc, 'c', num(p.Cdc), 'r', '0', 'vc_specify', 'on', 'vc_priority', 'High', 'vc', 'C.Vdc0');
ccs = blk(L.ccs, 'DC1_DC2');
cb = blk(L.cap, 'Cb');
set_param(cb, 'c', num(p.Cb), 'r', '0', 'vc_specify', 'on', 'vc_priority', 'High', 'vc', 'C.vb0');
lb = blk(L.ind, 'Lb');
set_param(lb, 'l', num(p.Lb), 'r', num(p.Rb), 'i_L_specify', 'on', 'i_L_priority', 'High', 'i_L', 'C.ib0');
l1 = blk(L.ind, 'L1');
set_param(l1, 'l', num(p.L1), 'r', num(p.R1), 'i_L_specify', 'on', 'i_L_priority', 'High', 'i_L', 'C.i10');
l2 = blk(L.ind, 'L2');
set_param(l2, 'l', num(p.L2), 'r', num(p.R2), 'i_L_specify', 'on', 'i_L_priority', 'High', 'i_L', 'C.i20');
e1 = blk(L.acv, 'E1');
set_param(e1, 'amp', num(p.E1), 'frequency', num(p.f), 'shift', 'C.theta0_deg');
e2 = blk(L.acv, 'E2');
set_param(e2, 'amp', num(p.E2), 'frequency', num(p.f), 'shift', ['C.theta0_deg + ' num(rad2deg(p.phi2))]);
for b = [cdc cb lb l1 l2 e1 e2]                % let run-specific values change without recompiling
    for f = {'vc_conf', 'i_L_conf', 'shift_conf'}
        try, set_param(b, f{1}, 'runtime'); catch, end
    end
end

legs = {'A', 'B', 'N'}; hsw = zeros(3, 2);
for k = 1:3
    for j = 1:2
        h = blk(L.isw, sprintf('SW_%s_%s', legs{k}, char('h' + (j == 2) * ('l' - 'h'))));
        set_param(h, 'Ron', num(sw.Ron), 'Vth', '0.5', 'Goff', '1e-6', ...
            'diode_param', 'ee.enum.converters.protectiondiode.nodynamics', ...
            'diode_Vf', num(sw.Vf), 'diode_Ron', num(sw.Ron));
        hsw(k, j) = h;
    end
end
% Bus connections (switch RConn1 = drain/positive, RConn2 = source/negative)
add_line(mdl, P(cdc).LConn(1), P(ccs).LConn(1));                      % DC+
add_line(mdl, P(cdc).RConn(1), P(ccs).RConn(2));                      % DC-
add_line(mdl, P(gnd).LConn(1), P(cdc).RConn(1));
add_line(mdl, P(sc).RConn(1), P(cdc).RConn(1));
for k = 1:3
    add_line(mdl, P(hsw(k, 1)).RConn(1), P(cdc).LConn(1));            % upper drain to DC+
    add_line(mdl, P(hsw(k, 1)).RConn(2), P(hsw(k, 2)).RConn(1));      % pole node
    add_line(mdl, P(hsw(k, 2)).RConn(2), P(cdc).RConn(1));            % lower source to DC-
end
poleA = P(hsw(1, 2)).RConn(1); poleB = P(hsw(2, 2)).RConn(1); poleN = P(hsw(3, 2)).RConn(1);
add_line(mdl, poleA, P(l1).LConn(1)); add_line(mdl, P(l1).RConn(1), P(e1).LConn(1));
add_line(mdl, P(e1).RConn(1), poleN);
add_line(mdl, poleB, P(l2).LConn(1)); add_line(mdl, P(l2).RConn(1), P(e2).LConn(1));
add_line(mdl, P(e2).RConn(1), poleN);
add_line(mdl, poleN, P(lb).LConn(1)); add_line(mdl, P(lb).RConn(1), P(cb).LConn(1));
add_line(mdl, P(cb).RConn(1), P(cdc).RConn(1));

% Probes: [i1 i2 ib vb Vdc] and the pole-N voltage
probes = {l1, 'i_L'; l2, 'i_L'; lb, 'i_L'; cb, 'vc'; cdc, 'vc'; hsw(3, 2), 'v'};
mx = blk('simulink/Signal Routing/Mux', 'MeasMux'); set_param(mx, 'Inputs', '5');
for k = 1:6
    pr = blk(L.probe, sprintf('Probe%d', k));
    simscape.probe.setBoundBlock(pr, probes{k, 1}); simscape.probe.setVariables(pr, string(probes{k, 2}));
    if k <= 5, add_line(mdl, sprintf('Probe%d/1', k), sprintf('MeasMux/%d', k)); end
end
ud = blk('simulink/Discrete/Unit Delay', 'MeasDelay');      % breaks algebraic loops
set_param(ud, 'SampleTime', 'Ts', 'InitialCondition', 'C.meas0');
add_line(mdl, 'MeasMux/1', 'MeasDelay/1');
udv = blk('simulink/Discrete/Unit Delay', 'VNDelay');
set_param(udv, 'SampleTime', 'Ts', 'InitialCondition', 'C.vb0');
add_line(mdl, 'Probe6/1', 'VNDelay/1');

% ---------------------------------------------------------------- Simulink control
blk('simulink/Sources/Clock', 'Clock');
mlfn(mdl, 'Monitor', 'function [ints, inst] = Monitor(meas, vN, tau, C)', ...
    '[ints, inst] = fine_monitor(meas, vN, tau, C);');
add_line(mdl, 'MeasDelay/1', 'Monitor/1'); add_line(mdl, 'VNDelay/1', 'Monitor/2');
add_line(mdl, 'Clock/1', 'Monitor/3');

zoh = {'ZOH_meas', 'MeasDelay/1'; 'ZOH_ints', 'Monitor/1'; 'ZOH_tau', 'Clock/1'};
for k = 1:3
    z = blk('simulink/Discrete/Zero-Order Hold', zoh{k, 1}); set_param(z, 'SampleTime', 'Tsw');
    add_line(mdl, zoh{k, 2}, [zoh{k, 1} '/1']);
end
mlfn(mdl, 'Controller', 'function [duty, pd, lg] = Controller(meas, ints, tau, C)', ...
    '[duty, pd, lg] = ctrl_step(meas, ints, tau, C);');
for k = 1:3, add_line(mdl, [zoh{k, 1} '/1'], sprintf('Controller/%d', k)); end

udd = blk('simulink/Discrete/Unit Delay', 'DutyDelay');     % one-period computation delay
set_param(udd, 'SampleTime', 'Tsw', 'InitialCondition', 'C.duty0');
add_line(mdl, 'Controller/1', 'DutyDelay/1');
mlfn(mdl, 'PWM', 'function g = PWM(duty, tau, C)', 'g = pwm_gates(duty, tau, C);');
add_line(mdl, 'DutyDelay/1', 'PWM/1'); add_line(mdl, 'Clock/1', 'PWM/2');
dm = blk('simulink/Signal Routing/Demux', 'GateDemux'); set_param(dm, 'Outputs', '6');
add_line(mdl, 'PWM/1', 'GateDemux/1');
order = [1 1; 1 2; 2 1; 2 2; 3 1; 3 2];                     % [A_hi A_lo B_hi B_lo N_hi N_lo]
for k = 1:6
    s = blk(L.sps, sprintf('GateSPS%d', k));
    add_line(mdl, sprintf('GateDemux/%d', k), sprintf('GateSPS%d/1', k));
    add_line(mdl, P(s).RConn(1), P(hsw(order(k, 1), order(k, 2))).LConn(1));
end

% DC1 minus DC2 as a current into DC+: (pd - Pdc2) / max(Vdc, 1)
sel = blk('simulink/Signal Routing/Selector', 'VdcSel');
set_param(sel, 'InputPortWidth', '5', 'Indices', '5');
add_line(mdl, 'MeasDelay/1', 'VdcSel/1');
mm = blk('simulink/Math Operations/MinMax', 'VdcFloor'); set_param(mm, 'Function', 'max', 'Inputs', '2');
one = blk('simulink/Sources/Constant', 'One'); set_param(one, 'Value', '1');
add_line(mdl, 'VdcSel/1', 'VdcFloor/1'); add_line(mdl, 'One/1', 'VdcFloor/2');
pdc2 = blk('simulink/Sources/Constant', 'Pdc2'); set_param(pdc2, 'Value', num(p.Pdc2));
sm = blk('simulink/Math Operations/Sum', 'NetPower'); set_param(sm, 'Inputs', '+-');
add_line(mdl, 'Controller/2', 'NetPower/1'); add_line(mdl, 'Pdc2/1', 'NetPower/2');
dv = blk('simulink/Math Operations/Product', 'Current'); set_param(dv, 'Inputs', '*/');
add_line(mdl, 'NetPower/1', 'Current/1'); add_line(mdl, 'VdcFloor/1', 'Current/2');
sps = blk(L.sps, 'CurrentSPS');
add_line(mdl, 'Current/1', 'CurrentSPS/1'); add_line(mdl, P(sps).RConn(1), P(ccs).RConn(1));

% Logging and protective trip
tw = blk('simulink/Sinks/To Workspace', 'LogCtl');
set_param(tw, 'VariableName', 'ctl_log', 'SaveFormat', 'Array');
add_line(mdl, 'Controller/3', 'LogCtl/1');
ti = blk('simulink/Sinks/To Workspace', 'LogInst');
set_param(ti, 'VariableName', 'inst_log', 'SaveFormat', 'Array', 'MaxDataPoints', '1');
add_line(mdl, 'Monitor/2', 'LogInst/1');
fm = blk('simulink/Signal Routing/Mux', 'FineMux'); set_param(fm, 'Inputs', '3');
add_line(mdl, 'Clock/1', 'FineMux/1'); add_line(mdl, 'MeasDelay/1', 'FineMux/2');
add_line(mdl, 'VNDelay/1', 'FineMux/3');
tf = blk('simulink/Sinks/To Workspace', 'LogFine');
set_param(tf, 'VariableName', 'fine_log', 'SaveFormat', 'Array', 'Decimation', 'fine_dec');
add_line(mdl, 'FineMux/1', 'LogFine/1');
ss = blk('simulink/Signal Routing/Selector', 'StopSel');
set_param(ss, 'InputPortWidth', '27', 'Indices', '27');
add_line(mdl, 'Controller/3', 'StopSel/1');
st = blk('simulink/Sinks/Stop Simulation', 'Trip');
add_line(mdl, 'StopSel/1', 'Trip/1');

Simulink.BlockDiagram.arrangeSystem(mdl);
save_system(mdl, fullfile(fileparts(mfilename('fullpath')), [mdl '.slx']));
fprintf('built %s\n', mdl);
end

function mlfn(mdl, name, sig, body)
% MATLAB Function block whose last argument C is a (tunable) parameter structure.
add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/' name]);
ch = find(sfroot, '-isa', 'Stateflow.EMChart', 'Path', [mdl '/' name]);
ch.Script = sprintf('%s\n%%#codegen\n%s\n', sig, body);
d = ch.find('-isa', 'Stateflow.Data', 'Name', 'C');
d.Scope = 'Parameter';
end
