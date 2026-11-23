function res = run_case(S, phase_deg, policy, td_us, zoom, Ts)
%RUN_CASE Run one step test in three_leg_switched.slx.
%   S = jsondecode(scenarios.json); policy = 'conventional' or 'plan'; td_us = dead time;
%   zoom = true logs the fine-rate signals every 1 us; Ts = solver step (default 0.1 us).
%   Returns metrics and time series.
if nargin < 5, zoom = false; end
if nargin < 6, Ts = 1e-7; end
mdl = 'three_leg_switched';
if ~bdIsLoaded(mdl), load_system(mdl); end
in = Simulink.SimulationInput(mdl);
C = case_params(S, phase_deg, strcmp(policy, 'plan'), td_us, Ts);
in = in.setVariable('C', C);
in = in.setVariable('Ts', C.Ts);
in = in.setVariable('Tsw', C.Tsw);
in = in.setVariable('Tstop', C.Tpre + S.params.n_cycles / S.params.f);
in = in.setVariable('fine_dec', ifelse(zoom, 10, 1e9));
t0 = tic;
out = sim(in);
res.wall_s = toc(t0);
lg = as_rows(out.ctl_log);
F = log_fields();
for k = 1:numel(F), res.ts.(F{k}) = lg(:, k); end
res.ts.t = res.ts.tau - C.Tpre;                       % time after the step
last = lg(end, :);
inst = reshape(out.inst_log(:, :, end), 1, []);    % only the last sample is kept
res.metrics = struct( ...
    'Vdc_min', last(24), 'Vdc_max', last(25), 'Vdc_inst_min', inst(1), 'Vdc_inst_max', inst(2), ...
    'I_pole_peak', last(26), 'I_pole_peak_inst', inst(3), 'sat_time_ms', 1e3 * last(21), ...
    'synth_time_ms', 1e3 * last(22), 'min_vN_headroom_V', last(23), 'violation', last(19) > 0, ...
    't_violation_ms', ifelse(last(20) >= 0, 1e3 * last(20), NaN), ...
    'current_violation_inst', inst(3) > S.params.I_pole_max);
if zoom
    fl = out.fine_log;
    res.fine = struct('t', fl(:, 1) - C.Tpre, 'i1', fl(:, 2), 'i2', fl(:, 3), 'ib', fl(:, 4), ...
                      'vb', fl(:, 5), 'Vdc', fl(:, 6), 'vN', fl(:, 7));
end
res.phase_deg = phase_deg; res.policy = policy; res.td_us = td_us;
end

function F = log_fields()
F = {'tau', 'i1', 'i2', 'ib', 'vb', 'Vdc', 'pd', 'vN_cmd', 'lo', 'hi', 'dA', 'dB', 'dN', ...
     'Vdc_avg', 'vb_avg', 'vN_avg', 'sat', 'synth', 'viol', 't_viol', 'sat_t', 'synth_t', ...
     'headroom', 'Vavg_min', 'Vavg_max', 'Iavg_pk', 'stop'};
end

function C = case_params(S, phase_deg, plan, td_us, Ts)
p = S.params; sw = S.switch;
sc = S.scenarios.(matlab.lang.makeValidName(sprintf('%g', phase_deg)));
C.Ts = Ts;                                     % solver step (0.1 us: 500 steps per period)
C.Tsw = 1 / sw.fsw;
C.Nper = round(C.Tsw / C.Ts);
C.td = td_us * 1e-6;
C.Nd = round(C.td / C.Ts);
C.f = p.f; C.w = 2 * pi * p.f; C.T = 1 / p.f;
C.t0 = sc.t0; C.theta0_deg = sc.theta0_deg; C.Tpre = S.pre_cycles / p.f;
C.E1 = p.E1; C.E2 = p.E2; C.phi2 = p.phi2;
C.L1 = p.L1; C.R1 = p.R1; C.L2 = p.L2; C.R2 = p.R2; C.Lb = p.Lb; C.Rb = p.Rb;
C.Cb = p.Cb; C.Cdc = p.Cdc;
C.P1 = p.P1; C.P2 = p.P2_pre; C.dP = S.dP; C.Pdc2 = p.Pdc2; C.Tr2 = p.Tr2;
C.Vnom = p.Vdc_nom; C.Vmin = p.Vdc_min; C.Vmax = p.Vdc_max; C.Imax = p.I_pole_max;
C.slew = p.slew_dc1; C.Pmin = p.Pdc1_min; C.Pmax = p.Pdc1_max;
C.wc = S.gains.wc; C.k_bus = S.gains.k_bus; C.k_rec = S.gains.k_rec; C.k_tot = S.gains.k_tot;
C.wc_ac = sw.wc_ac;
C.P_pre = S.P_pre; C.P_post = S.P_post; C.Ebref_pre = S.Eb_ref_pre; C.Ebref_post = S.Eb_ref_post;
C.M = S.M; C.hist0 = sc.hist(:); C.pd0 = sc.pd;
C.policy = double(plan);
C.pTs = S.plan.Ts; C.pt = S.plan.t(:); C.pEb = S.plan.Eb(:); C.pEd = S.plan.Ed(:);
C.pdEd = gradient(C.pEd, C.pTs); C.ppd = S.plan.pd(:);
C.trip_after = 1e-3;
C.Vdc0 = sc.Vdc; C.vb0 = sc.vb; C.ib0 = sc.ib; C.i10 = sc.i1; C.i20 = sc.i2;
C.meas0 = [sc.i1; sc.i2; sc.ib; sc.vb; sc.Vdc];
% Duties for the first period, from the same law with the initial state
t = C.t0 + 0.5 * C.Tsw;
w = C.w; I1 = 2 * C.P1 / C.E1; A2 = 2 * C.P2 / C.E2;
u1 = C.E1 * sin(w * t) + C.L1 * I1 * w * cos(w * t) + C.R1 * I1 * sin(w * t);
u2 = C.E2 * sin(w * t + C.phi2) + C.L2 * A2 * w * cos(w * t + C.phi2) + C.R2 * A2 * sin(w * t + C.phi2);
vN = C.vb0 + C.Rb * C.ib0;
C.duty0 = min(max([(u1 + vN); (u2 + vN); vN] / C.Vdc0, 0), 1);
end

function x = as_rows(x)
% To Workspace saves matrix-shaped signals as rows x cols x samples; return samples x values.
if ndims(x) == 3
    x = reshape(x, [], size(x, 3)).';
end
end

function y = ifelse(c, a, b)
if c, y = a; else, y = b; end
end
