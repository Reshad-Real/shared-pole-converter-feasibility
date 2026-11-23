function run_r22(keys)
%RUN_R22 Run the R22 benchmark cases in Simscape and save the last 0.1 s of each.
%   run_r22() runs all rectifier cases and the decoupling-off case; run_r22({'rect_phi0_GDPWM'})
%   runs a subset. Output: sim_r22/results/simscape/<case>.mat with per-period rows
%   [tau i_a i_c u_cac U U_avg dA dB dC] and fine-rate [t i_a i_c u_cac U] at 1 us.
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(here, '..', '..', 'sim_simscape'));
outdir = fullfile(here, '..', 'results', 'simscape'); if ~exist(outdir, 'dir'), mkdir(outdir); end
if ~exist(fullfile(here, 'r22_switched.slx'), 'file'), build_r22(); end
load_system(fullfile(here, 'r22_switched.slx'));
cases = case_list();
if nargin > 0, cases = cases(ismember({cases.key}, keys)); end
for k = 1:numel(cases)
    c = cases(k);
    C = params(c);
    in = Simulink.SimulationInput('r22_switched');
    in = in.setVariable('C', C).setVariable('Ts', C.Ts).setVariable('Tsw', C.Tsw) ...
           .setVariable('Tstop', 0.25).setVariable('fine_dec', 20);
    t0 = tic; out = sim(in); wall = toc(t0);
    lg = as_rows(out.ctl_log); fl = out.fine_log;
    keep = lg(:, 1) >= 0.15; keepf = fl(:, 1) >= 0.15;
    per = lg(keep, :); fine = fl(keepf, :); meta = c; meta.wall_s = wall; meta.Tsw = C.Tsw; meta.f = C.f;
    save(fullfile(outdir, [c.key '.mat']), 'per', 'fine', 'meta');
    fprintf('%-22s %5.0f s  U mean %.1f V\n', c.key, wall, mean(per(:, 6)));
end
end

function cs = case_list()
cs = struct('key', {}, 'scheme', {}, 'phi_deg', {}, 'R_load', {}, 'apd', {});
cs(end + 1) = struct('key', 'apd_off', 'scheme', 3, 'phi_deg', 0, 'R_load', 80, 'apd', 0);
names = {'SVPWM', 'DPWM1', 'GDPWM'}; loads = [0 80; 15 83; 45 113];
for i = 1:3
    for s = 1:3
        cs(end + 1) = struct('key', sprintf('rect_phi%d_%s', loads(i, 1), names{s}), 'scheme', s, ...
                             'phi_deg', loads(i, 1), 'R_load', loads(i, 2), 'apd', 1); %#ok<AGROW>
    end
end
end

function C = params(c)
% Circuit, control and initial state of MODEL_SPEC.md for one case.
C.Ts = 50e-9; C.Tsw = 1 / 40e3; C.Nper = round(C.Tsw / C.Ts); C.td = 200e-9; C.Nd = round(C.td / C.Ts);
C.f = 50; C.w = 2 * pi * C.f; C.U_ac = 220 * sqrt(2);
C.L_ac = 1.44e-3; C.R_ac = 0.1; C.L_c = 1.15e-3; C.R_c = 0.1; C.C_ac = 130e-6; C.C_dc = 135e-6;  % C_ac from R22 eq. (9)
C.U_ref = 400; C.Ron = 0.15; C.Vf = 3.0;
C.w_ci = 2 * pi * 2000; C.w_cc = 2 * pi * 2000; C.w_cv = 2 * pi * 300; w_dc = 2 * pi * 10;
g = C.U_ac / (2 * C.C_dc * C.U_ref); C.Kp = w_dc / g; C.Ki = C.Kp * w_dc / 4;
C.n_maf = round(0.01 / C.Tsw);
C.statcom = 0; C.I_q = 0; C.phi = deg2rad(c.phi_deg); C.apd = c.apd; C.scheme = c.scheme;
C.R_load = c.R_load;
C.I_d0 = 2 * (C.U_ref^2 / C.R_load) * 1.02 / (C.U_ac * cos(C.phi));
I = C.I_d0; beta = C.phi;
C.U0 = C.U_ref; C.ia0 = I * sin(beta);
if C.apd
    X = C.U_ac - (C.R_ac + 1i * C.w * C.L_ac) * I * exp(1i * beta); Y = I * exp(1i * beta);
    A = abs(X * Y) / 2; th = (angle(X * Y) - pi / 2) / 2;
    V = sqrt(2 * A / (C.w * C.C_ac * (1 - C.w^2 * C.L_c * C.C_ac)));
    C.uc0 = V * sin(th); C.ic0 = -C.C_ac * C.w * V * cos(th);
else
    C.uc0 = 0; C.ic0 = 0;
end
C.meas0 = [C.ia0; C.ic0; C.uc0; C.U0];
% Duties for the first period (SVPWM-style split; the run settles within the first cycles)
t = 0.5 * C.Tsw;
u_ab = C.U_ac * sin(C.w * t) - C.L_ac * I * C.w * cos(C.w * t + beta);
u_cb = C.apd * C.uc0;
u = [(2 * u_ab - u_cb) / 3; (-u_ab - u_cb) / 3; (-u_ab + 2 * u_cb) / 3];
C.duty0 = min(max(0.5 + (u - (max(u) + min(u)) / 2) / C.U0, 0), 1);
end

function x = as_rows(x)
if ndims(x) == 3, x = reshape(x, [], size(x, 3)).'; end
end
