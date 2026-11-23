function run_check()
%RUN_CHECK Run the WP2 campaign in the Simscape model and save results for comparison.
%   Pair runs (conventional at theta_a and theta_b, LP plan at theta_a) at dead time 1 us with
%   fine logging, the pair at 0 and 2 us, and the 24-phase conventional sweep at 1 and 2 us.
%   Output: results/simscape_results.json, read by compare_simscape.py.
here = fileparts(mfilename('fullpath'));
addpath(here);
S = jsondecode(fileread(fullfile(here, 'scenarios.json')));
if ~exist(fullfile(here, 'three_leg_switched.slx'), 'file'), build_model(S); end
load_system(fullfile(here, 'three_leg_switched.slx'));
outdir = fullfile(here, 'results'); if ~exist(outdir, 'dir'), mkdir(outdir); end
th_a = S.theta_a_deg; th_b = S.theta_b_deg;

R = struct('pair', struct(), 'pair_td', struct(), 'sweep', struct());
cases = {'conventional_theta_a', th_a, 'conventional'; 'conventional_theta_b', th_b, 'conventional'; ...
         'plan_theta_a', th_a, 'plan'};
for k = 1:size(cases, 1)
    r = run_case(S, cases{k, 2}, cases{k, 3}, 1, true);
    R.pair.(cases{k, 1}) = slim(r, true);
    fprintf('%-22s td 1 us: violation=%d t_viol=%.2f ms Vdc %.1f-%.1f (%.0f s)\n', cases{k, 1}, ...
        r.metrics.violation, r.metrics.t_violation_ms, r.metrics.Vdc_min, r.metrics.Vdc_max, r.wall_s);
    for td = [0 2]
        r = run_case(S, cases{k, 2}, cases{k, 3}, td, false);
        R.pair_td.(sprintf('%s_td%dus', cases{k, 1}, td)) = r.metrics;
        fprintf('%-22s td %d us: violation=%d t_viol=%.2f ms\n', cases{k, 1}, td, ...
            r.metrics.violation, r.metrics.t_violation_ms);
    end
end
for td = [1 2]
    rows = repmat(struct('phase_deg', 0, 'metrics', struct()), numel(S.sweep_deg), 1);
    for q = 1:numel(S.sweep_deg)
        r = run_case(S, S.sweep_deg(q), 'conventional', td, false);
        rows(q).phase_deg = S.sweep_deg(q); rows(q).metrics = r.metrics;
        fprintf('sweep td %d us, %5.0f deg: violation=%d (%.0f s)\n', td, S.sweep_deg(q), ...
            r.metrics.violation, r.wall_s);
    end
    R.sweep.(sprintf('td%dus', td)) = rows;
end
R.solver = struct('Ts', 1e-7, 'local_solver', 'partitioning', 'switch', 'Ideal Semiconductor Switch', ...
                  'matlab', version);
fid = fopen(fullfile(outdir, 'simscape_results.json'), 'w');
fwrite(fid, jsonencode(R)); fclose(fid);
fprintf('saved %s\n', fullfile(outdir, 'simscape_results.json'));
end

function s = slim(r, keep_ts)
% Metrics plus per-period series and a fine-rate window around the step (-0.5 to 3 ms).
s.metrics = r.metrics;
if keep_ts
    f = {'t', 'vb', 'ib', 'Vdc', 'pd', 'vN_cmd', 'lo', 'hi', 'Vdc_avg', 'vb_avg', 'vN_avg'};
    for k = 1:numel(f), s.ts.(f{k}) = r.ts.(f{k}); end
    if isfield(r, 'fine')
        w = r.fine.t >= -0.5e-3 & r.fine.t <= 3e-3;
        g = fieldnames(r.fine);
        for k = 1:numel(g), s.fine.(g{k}) = r.fine.(g{k})(w); end
    end
end
end
