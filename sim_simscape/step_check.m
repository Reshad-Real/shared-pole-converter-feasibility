function step_check()
%STEP_CHECK Solver-step convergence at the tightest 2 us dead-time case (90 deg).
%   Reruns the case at 0.1 us and 0.05 us and saves results/step_check.json.
here = fileparts(mfilename('fullpath')); addpath(here);
S = jsondecode(fileread(fullfile(here, 'scenarios.json')));
load_system(fullfile(here, 'three_leg_switched.slx'));
R = struct('phase_deg', 90, 'td_us', 2, 'runs', struct('Ts_us', {}, 'metrics', {}));
for Ts = [1e-7 5e-8]
    r = run_case(S, 90, 'conventional', 2, false, Ts);
    R.runs(end + 1) = struct('Ts_us', Ts * 1e6, 'metrics', r.metrics);
    fprintf('Ts %.2f us: Vavg min %.2f, violation=%d at %.2f ms\n', Ts * 1e6, ...
        r.metrics.Vdc_min, r.metrics.violation, r.metrics.t_violation_ms);
end
fid = fopen(fullfile(here, 'results', 'step_check.json'), 'w');
fwrite(fid, jsonencode(R)); fclose(fid);
end
