function test_ports()
%TEST_PORTS Establish terminal orientation of the switch block and current-source polarity.
%   1) Ideal Semiconductor Switch, gate off, 10 V through 10 ohm applied from RConn2 to
%      RConn1: current flows only if the protection diode points RConn2 -> RConn1, which
%      means RConn1 is the positive (drain/collector) terminal.
%   2) Gate on: the switch must conduct (about 1 A).
%   3) Controlled Current Source with +1 A, LConn1 -> capacitor +, RConn2 -> reference:
%      the capacitor charges positive if current leaves LConn1 into the circuit.
load_system('ee_lib'); load_system('fl_lib'); load_system('nesl_utility');
L = lib_paths();
for gate = [0 1]
    o = run_switch_test(L, gate);
    fprintf('switch gate=%d, 10 V from RConn2 to RConn1 via 10 ohm: current = %+.3f A\n', gate, o);
end
fprintf('current source +1 A for 2 us into 1 uF (LConn1 -> C+): vC = %+.3f V\n', run_ccs_test(L));
end

function iR = run_switch_test(L, gate)
mdl = 'tmp_sw'; if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl); set_param(mdl, 'Solver', 'FixedStepDiscrete', 'FixedStep', '1e-7', 'StopTime', '2e-6');
sc = add_block(L.solver, [mdl '/Solver']);
set_param(sc, 'UseLocalSolver', 'on', 'LocalSolverChoice', 'NE_PARTITIONING_ADVANCER', 'LocalSolverSampleTime', '1e-7');
gnd = add_block(L.ref, [mdl '/Gnd']);
v = add_block(L.dcv, [mdl '/V']); set_param(v, 'v0', '10');
r = add_block(L.res, [mdl '/R']); set_param(r, 'R', '10');
sw = add_block(L.isw, [mdl '/SW']);
set_param(sw, 'Ron', '0.05', 'Vth', '0.5', 'diode_param', 'ee.enum.converters.protectiondiode.nodynamics', ...
    'diode_Vf', '0.8', 'diode_Ron', '0.05');
k = add_block('simulink/Sources/Constant', [mdl '/K']); set_param(k, 'Value', num2str(gate));
s = add_block(L.sps, [mdl '/S']);
P = @(h) get_param(h, 'PortHandles');
add_line(mdl, 'K/1', 'S/1');
add_line(mdl, P(s).RConn(1), P(sw).LConn(1));
add_line(mdl, P(v).LConn(1), P(r).LConn(1));          % source + -> resistor
add_line(mdl, P(r).RConn(1), P(sw).RConn(2));         % resistor -> switch RConn2
add_line(mdl, P(sw).RConn(1), P(v).RConn(1));         % switch RConn1 -> source -
add_line(mdl, P(gnd).LConn(1), P(v).RConn(1));
add_line(mdl, P(sc).RConn(1), P(v).RConn(1));
iR = probe_out(L, mdl, r, "i");
close_system(mdl, 0);
end

function vC = run_ccs_test(L)
mdl = 'tmp_cs'; if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl); set_param(mdl, 'Solver', 'FixedStepDiscrete', 'FixedStep', '1e-7', 'StopTime', '2e-6');
sc = add_block(L.solver, [mdl '/Solver']);
set_param(sc, 'UseLocalSolver', 'on', 'LocalSolverChoice', 'NE_PARTITIONING_ADVANCER', 'LocalSolverSampleTime', '1e-7');
gnd = add_block(L.ref, [mdl '/Gnd']);
cs = add_block(L.ccs, [mdl '/CS']); c = add_block(L.cap, [mdl '/C']); set_param(c, 'c', '1e-6');
k = add_block('simulink/Sources/Constant', [mdl '/K']); s = add_block(L.sps, [mdl '/S']);
P = @(h) get_param(h, 'PortHandles');
add_line(mdl, 'K/1', 'S/1');
add_line(mdl, P(s).RConn(1), P(cs).RConn(1));         % PS command
add_line(mdl, P(cs).LConn(1), P(c).LConn(1));         % LConn1 -> C+
add_line(mdl, P(cs).RConn(2), P(c).RConn(1));         % RConn2 -> C-
add_line(mdl, P(gnd).LConn(1), P(c).RConn(1));
add_line(mdl, P(sc).RConn(1), P(c).RConn(1));
vC = probe_out(L, mdl, c, "vc");
close_system(mdl, 0);
end

function y = probe_out(L, mdl, target, var)
pr = add_block(L.probe, [mdl '/Probe']);
simscape.probe.setBoundBlock(pr, target); simscape.probe.setVariables(pr, var);
o = add_block('simulink/Sinks/To Workspace', [mdl '/Out']);
set_param(o, 'VariableName', 'y', 'SaveFormat', 'Array');
add_line(mdl, 'Probe/1', 'Out/1');
out = sim(mdl);
y = out.y(end);
end
