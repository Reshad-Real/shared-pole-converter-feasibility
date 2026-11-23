function [duty, pd_out, lg] = ctrl_step(meas, ints, tau, C)
%CTRL_STEP Controller of CONTROLLER_SPEC.md section 4, run once per switching period.
%   meas = [i1 i2 ib vb Vdc] sampled at the period start; ints = running integrals of
%   [Vdc vb iA iB iN vN]; tau = simulation time (absolute time is C.t0 + tau).
%   duty = [dA dB dN] for the next period; pd_out = DC1 power; lg = log row (see LOG_FIELDS).
%#codegen
persistent init pd ibref_prev hist hpos hsum ints_prev t_viol sat_t synth_t headroom ...
    Vavg_min Vavg_max Iavg_pk viol_any
Tsw = C.Tsw;
if isempty(init) || tau < 0.5 * C.Ts
    init = true; pd = C.pd0; ibref_prev = NaN;
    hist = C.hist0(:); hpos = 1; hsum = sum(hist);
    ints_prev = ints(:); t_viol = -1; sat_t = 0; synth_t = 0; headroom = Inf;
    Vavg_min = Inf; Vavg_max = -Inf; Iavg_pk = 0; viol_any = 0;
end
i1 = meas(1); i2 = meas(2); ib = meas(3); vb = meas(4); V = meas(5);
t = C.t0 + tau;                         % absolute time at the sample
t_step = C.t0 + C.Tpre;
t_app = t + 1.5 * Tsw;                  % middle of the period the new duties act in
after = t >= t_step;

% ---- references
[u1a, u2a, i1a, i2a, paca] = ac_refs(t_app, t_step, C);
[~, ~, i1s, i2s, ~] = ac_refs(t, t_step, C);

% ---- 1. energies
Eb = 0.5 * C.Cb * vb^2;
Ed = 0.5 * C.Cdc * V^2;
Eb_avg = hsum / C.M;

% ---- 2. supervisor
if C.policy == 1 && after
    trel_k = t - t_step; trel_a = t_app - t_step;
    j = min(floor(trel_k / C.pTs), numel(C.ppd) - 1) + 1;
    Ed_ref = interp_clamp(C.pt, C.pEd, trel_a);
    dEd_ff = interp_clamp(C.pt, C.pdEd, trel_a);
    Etot_plan = interp_clamp(C.pt, C.pEb, trel_k) + interp_clamp(C.pt, C.pEd, trel_k);
    pd_cmd = C.ppd(j) + C.k_tot * (Etot_plan - (Eb + Ed));
else
    Ed_ref = 0.5 * C.Cdc * C.Vnom^2; dEd_ff = 0;
    x = min(max((t - t_step) / C.Tr2, 0), 1);
    P_cmd = C.P_pre + (C.P_post - C.P_pre) * 0.5 * (1 - cos(pi * x));
    if after, Eb_ref = C.Ebref_post; else, Eb_ref = C.Ebref_pre; end
    pd_cmd = P_cmd + C.k_rec * (Eb_ref - Eb_avg);
end

% ---- 3. DC1 power with slew limit
stp = C.slew * Tsw;
pd = min(max(pd + min(max(pd_cmd - pd, -stp), stp), C.Pmin), C.Pmax);

% ---- 4-5. buffer leg
pb_ref = pd - C.Pdc2 - paca - dEd_ff + C.k_bus * (Ed - Ed_ref);
ib_ref = pb_ref / max(vb, 1);
if isnan(ibref_prev), dib_ff = 0; else, dib_ff = (ib_ref - ibref_prev) / Tsw; end
ibref_prev = ib_ref;
vN_cmd = vb + C.Rb * ib + C.Lb * (dib_ff + C.wc * (ib_ref - ib));

% ---- 6. AC legs
u1c = u1a + C.L1 * C.wc_ac * (i1s - i1);
u2c = u2a + C.L2 * C.wc_ac * (i2s - i2);

% ---- 7. pole-N band
lo = max([0, -u1c, -u2c]);
hi = V - max([0, u1c, u2c]);
synth = max(lo - hi, 0);
if synth > 0
    vN = 0.5 * (lo + hi);
else
    vN = min(max(vN_cmd, lo), hi);
end

% ---- 8. duties with dead-time compensation
comp = C.td / Tsw;
iref = [i1a; i2a; ib_ref - i1a - i2a];
vref = [u1c + vN; u2c + vN; vN];
duty = zeros(3, 1);
for k = 1:3
    duty(k) = min(max(vref(k) / V + comp * min(max(iref(k) / 0.3, -1), 1), 0), 1);
end
pd_out = pd;

% ---- statistics for periods after the step (CONTROLLER_SPEC.md section 5)
if after
    sat_t = sat_t + Tsw * double(vN ~= vN_cmd);
    synth_t = synth_t + Tsw * double(synth > 0);
    headroom = min([headroom, vN_cmd - lo, hi - vN_cmd]);
end
avg = (ints(:) - ints_prev) / Tsw;      % averages over the period that just ended
ints_prev = ints(:);
t_prev = t - Tsw;                       % start of that period
stop = 0;
if tau > 0.5 * Tsw && t_prev >= t_step
    Vavg_min = min(Vavg_min, avg(1)); Vavg_max = max(Vavg_max, avg(1));
    Iavg_pk = max([Iavg_pk, abs(avg(3)), abs(avg(4)), abs(avg(5))]);
    bad = avg(1) < C.Vmin || avg(1) > C.Vmax || Iavg_pk > C.Imax || synth > 0;
    if bad && t_viol < 0, t_viol = t_prev - t_step; viol_any = 1; end
    if t_viol >= 0 && (t_prev - t_step) > t_viol + C.trip_after, stop = 1; end
end

% ---- moving-average history (sampled buffer energy)
hsum = hsum + Eb - hist(hpos);
hist(hpos) = Eb;
hpos = mod(hpos, C.M) + 1;

lg = [tau, i1, i2, ib, vb, V, pd, vN_cmd, lo, hi, duty(1), duty(2), duty(3), ...
      avg(1), avg(2), avg(6), double(vN ~= vN_cmd), synth, viol_any, t_viol, ...
      sat_t, synth_t, headroom, Vavg_min, Vavg_max, Iavg_pk, stop];
end

function [u1, u2, i1, i2, pac] = ac_refs(t, t_step, C)
% AC references of CONTROLLER_SPEC.md section 3 at absolute time t.
w = C.w;
I1 = 2 * C.P1 / C.E1;
I2a = 2 * C.P2 / C.E2;
I2b = 2 * (C.P2 + C.dP) / C.E2;
x = min(max((t - t_step) / C.Tr2, 0), 1);
A2 = I2a + (I2b - I2a) * 0.5 * (1 - cos(pi * x));
dA2 = (I2b - I2a) * 0.5 * pi / C.Tr2 * sin(pi * x);
s1 = sin(w * t); c1 = cos(w * t);
s2 = sin(w * t + C.phi2); c2 = cos(w * t + C.phi2);
i1 = I1 * s1;
i2 = A2 * s2;
u1 = C.E1 * s1 + C.L1 * I1 * w * c1 + C.R1 * i1;
u2 = C.E2 * s2 + C.L2 * (dA2 * s2 + A2 * w * c2) + C.R2 * i2;
pac = u1 * i1 + u2 * i2;
end

function y = interp_clamp(xs, ys, xq)
% Linear interpolation that holds the end values outside the table (like numpy.interp).
xq = min(max(xq, xs(1)), xs(end));
y = interp1(xs, ys, xq, 'linear');
end
