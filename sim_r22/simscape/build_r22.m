function build_r22()
%BUILD_R22 Create r22_switched.slx: R22's three-leg decoupling converter in Simscape Electrical.
%   Grid (with L_ac) between legs a and b, C_ac with L_c between legs c and b, C_dc with a
%   resistive load. Run-specific values come from the workspace structure C (run_r22.m).
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(here, '..', '..', 'sim_simscape'));
L = lib_paths();
load_system('ee_lib'); load_system('fl_lib'); load_system('nesl_utility');
mdl = 'r22_switched';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
set_param(mdl, 'Solver', 'FixedStepDiscrete', 'FixedStep', 'Ts', 'StopTime', 'Tstop', ...
    'SolverMode', 'SingleTasking', 'SingleTaskRateTransMsg', 'none', ...
    'ReturnWorkspaceOutputs', 'on', 'SignalLogging', 'off', 'SaveTime', 'off', 'SaveOutput', 'off');
P = @(h) get_param(h, 'PortHandles');
blk = @(lib, name) add_block(lib, [mdl '/' name]);

sc = blk(L.solver, 'Solver');
set_param(sc, 'UseLocalSolver', 'on', 'LocalSolverChoice', 'NE_PARTITIONING_ADVANCER', ...
    'LocalSolverSampleTime', 'Ts', 'DoFixedCost', 'on');
gnd = blk(L.ref, 'Gnd');
cdc = blk(L.cap, 'Cdc');
set_param(cdc, 'c', 'C.C_dc', 'r', '0', 'vc_specify', 'on', 'vc_priority', 'High', 'vc', 'C.U0');
rl = blk(L.res, 'Rload'); set_param(rl, 'R', 'C.R_load');
lac = blk(L.ind, 'Lac');
set_param(lac, 'l', 'C.L_ac', 'r', 'C.R_ac', 'i_L_specify', 'on', 'i_L_priority', 'High', 'i_L', 'C.ia0');
lc = blk(L.ind, 'Lc');
set_param(lc, 'l', 'C.L_c', 'r', 'C.R_c', 'i_L_specify', 'on', 'i_L_priority', 'High', 'i_L', 'C.ic0');
cac = blk(L.cap, 'Cac');
set_param(cac, 'c', 'C.C_ac', 'r', '0', 'vc_specify', 'on', 'vc_priority', 'High', 'vc', 'C.uc0');
src = blk(L.acv, 'Grid');
set_param(src, 'amp', 'C.U_ac', 'frequency', 'C.f', 'shift', '0');

legs = {'a', 'b', 'c'}; hsw = zeros(3, 2);
for k = 1:3
    for j = 1:2
        h = blk(L.isw, sprintf('S_%s_%s', legs{k}, char('h' + (j == 2) * ('l' - 'h'))));
        set_param(h, 'Ron', 'C.Ron', 'Vth', '0.5', 'Goff', '1e-6', ...
            'diode_param', 'ee.enum.converters.protectiondiode.nodynamics', ...
            'diode_Vf', 'C.Vf', 'diode_Ron', 'C.Ron');
        hsw(k, j) = h;
    end
end
add_line(mdl, P(gnd).LConn(1), P(cdc).RConn(1));
add_line(mdl, P(sc).RConn(1), P(cdc).RConn(1));
add_line(mdl, P(rl).LConn(1), P(cdc).LConn(1)); add_line(mdl, P(rl).RConn(1), P(cdc).RConn(1));
for k = 1:3
    add_line(mdl, P(hsw(k, 1)).RConn(1), P(cdc).LConn(1));
    add_line(mdl, P(hsw(k, 1)).RConn(2), P(hsw(k, 2)).RConn(1));
    add_line(mdl, P(hsw(k, 2)).RConn(2), P(cdc).RConn(1));
end
pa = P(hsw(1, 2)).RConn(1); pb = P(hsw(2, 2)).RConn(1); pc = P(hsw(3, 2)).RConn(1);
add_line(mdl, P(src).LConn(1), P(lac).LConn(1));      % grid + -> L_ac -> node a
add_line(mdl, P(lac).RConn(1), pa);
add_line(mdl, P(src).RConn(1), pb);                    % grid - -> node b
add_line(mdl, P(cac).LConn(1), P(lc).LConn(1));        % C_ac + -> L_c -> node c
add_line(mdl, P(lc).RConn(1), pc);
add_line(mdl, P(cac).RConn(1), pb);                    % C_ac - -> node b

% Measurements [i_a i_c u_cac U], delayed one step to break algebraic loops
probes = {lac, 'i_L'; lc, 'i_L'; cac, 'vc'; cdc, 'vc'};
mx = blk('simulink/Signal Routing/Mux', 'MeasMux'); set_param(mx, 'Inputs', '4');
for k = 1:4
    pr = blk(L.probe, sprintf('Probe%d', k));
    simscape.probe.setBoundBlock(pr, probes{k, 1}); simscape.probe.setVariables(pr, string(probes{k, 2}));
    add_line(mdl, sprintf('Probe%d/1', k), sprintf('MeasMux/%d', k));
end
ud = blk('simulink/Discrete/Unit Delay', 'MeasDelay');
set_param(ud, 'SampleTime', 'Ts', 'InitialCondition', 'C.meas0');
add_line(mdl, 'MeasMux/1', 'MeasDelay/1');
sel = blk('simulink/Signal Routing/Selector', 'USel'); set_param(sel, 'InputPortWidth', '4', 'Indices', '4');
add_line(mdl, 'MeasDelay/1', 'USel/1');
blk('simulink/Sources/Clock', 'Clock');

% Fine-rate PWM and integral of U
mlfn(mdl, 'Fine', 'function [g, intU] = Fine(duty, U, tau, C)', '[g, intU] = fine_r22(duty, U, tau, C);');
% Controller at the switching rate
zoh = {'ZOH_meas', 'MeasDelay/1'; 'ZOH_int', 'Fine/2'; 'ZOH_tau', 'Clock/1'};
for k = 1:3
    z = blk('simulink/Discrete/Zero-Order Hold', zoh{k, 1}); set_param(z, 'SampleTime', 'Tsw');
    add_line(mdl, zoh{k, 2}, [zoh{k, 1} '/1']);
end
mlfn(mdl, 'Controller', 'function [duty, lg] = Controller(meas, intU, tau, C)', ...
    '[duty, lg] = ctrl_r22(meas, intU, tau, C);');
for k = 1:3, add_line(mdl, [zoh{k, 1} '/1'], sprintf('Controller/%d', k)); end
udd = blk('simulink/Discrete/Unit Delay', 'DutyDelay');
set_param(udd, 'SampleTime', 'Tsw', 'InitialCondition', 'C.duty0');
add_line(mdl, 'Controller/1', 'DutyDelay/1');
add_line(mdl, 'DutyDelay/1', 'Fine/1'); add_line(mdl, 'USel/1', 'Fine/2'); add_line(mdl, 'Clock/1', 'Fine/3');
dm = blk('simulink/Signal Routing/Demux', 'GateDemux'); set_param(dm, 'Outputs', '6');
add_line(mdl, 'Fine/1', 'GateDemux/1');
order = [1 1; 1 2; 2 1; 2 2; 3 1; 3 2];
for k = 1:6
    s = blk(L.sps, sprintf('GateSPS%d', k));
    add_line(mdl, sprintf('GateDemux/%d', k), sprintf('GateSPS%d/1', k));
    add_line(mdl, P(s).RConn(1), P(hsw(order(k, 1), order(k, 2))).LConn(1));
end

% Logs: controller rows plus the applied duties, and fine-rate i_a for the THD
tw = blk('simulink/Sinks/To Workspace', 'LogCtl');
set_param(tw, 'VariableName', 'ctl_log', 'SaveFormat', 'Array');
add_line(mdl, 'Controller/2', 'LogCtl/1');
fm = blk('simulink/Signal Routing/Mux', 'FineMux'); set_param(fm, 'Inputs', '2');
add_line(mdl, 'Clock/1', 'FineMux/1'); add_line(mdl, 'MeasDelay/1', 'FineMux/2');
tf = blk('simulink/Sinks/To Workspace', 'LogFine');
set_param(tf, 'VariableName', 'fine_log', 'SaveFormat', 'Array', 'Decimation', 'fine_dec');
add_line(mdl, 'FineMux/1', 'LogFine/1');

Simulink.BlockDiagram.arrangeSystem(mdl);
save_system(mdl, fullfile(here, [mdl '.slx']));
fprintf('built %s\n', mdl);
end

function mlfn(mdl, name, sig, body)
add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/' name]);
ch = find(sfroot, '-isa', 'Stateflow.EMChart', 'Path', [mdl '/' name]);
ch.Script = sprintf('%s\n%%#codegen\n%s\n', sig, body);
d = ch.find('-isa', 'Stateflow.Data', 'Name', 'C');
d.Scope = 'Parameter';
end
