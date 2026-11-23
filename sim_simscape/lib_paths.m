function L = lib_paths()
%LIB_PATHS Library paths of the Simscape blocks used (several contain line breaks).
nl = newline;
L.isw    = ['ee_lib/Semiconductors &' nl 'Converters/Ideal Semiconductor' nl 'Switch'];
L.cap    = 'fl_lib/Electrical/Electrical Elements/Capacitor';
L.ind    = 'fl_lib/Electrical/Electrical Elements/Inductor';
L.res    = 'fl_lib/Electrical/Electrical Elements/Resistor';
L.ref    = 'fl_lib/Electrical/Electrical Elements/Electrical Reference';
L.acv    = 'fl_lib/Electrical/Electrical Sources/AC Voltage Source';
L.dcv    = 'fl_lib/Electrical/Electrical Sources/DC Voltage Source';
L.ccs    = ['fl_lib/Electrical/Electrical Sources/Controlled Current' nl 'Source'];
L.solver = ['nesl_utility/Solver' nl 'Configuration'];
L.sps    = ['nesl_utility/Simulink-PS' nl 'Converter'];
L.probe  = 'nesl_utility/Probe';
end
