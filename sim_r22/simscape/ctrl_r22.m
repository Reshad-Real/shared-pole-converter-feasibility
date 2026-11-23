function [duty, lg] = ctrl_r22(meas, intU, tau, C)
%CTRL_R22 R22 control law (sim_r22/MODEL_SPEC.md), run once per switching period.
%   meas = [i_a i_c u_cac U] sampled at the period start; intU = running integral of U;
%   tau = simulation time (absolute: u_ac = U_ac sin(w tau)). duty = [dA dB dC] for the
%   next period; lg = [tau i_a i_c u_cac U U_avg_prev_period dA dB dC], where the duties are
%   those acting in the current period (computed one period earlier).
%#codegen
persistent init maf mpos msum integ intU_prev dprev
Tsw = C.Tsw; w = C.w;
if isempty(init) || tau < 0.5 * C.Ts
    init = true; maf = C.U_ref * ones(C.n_maf, 1); mpos = 1; msum = sum(maf);
    integ = C.I_d0; intU_prev = intU; dprev = C.duty0(:);
end
i_a = meas(1); i_c = meas(2); u_cac = meas(3); U = meas(4);

% 1. DC voltage: 10 ms moving average, PI -> grid-current amplitude
msum = msum + U - maf(mpos); maf(mpos) = U; mpos = mod(mpos, C.n_maf) + 1;
e = C.U_ref - msum / C.n_maf;
integ = integ + C.Ki * e * Tsw;
I_d = C.Kp * e + integ;
if C.statcom
    I = hypot(I_d, C.I_q); beta = atan2(C.I_q, I_d);
else
    I = I_d; beta = C.phi;
end

% 2-3. grid-current reference and loop
t_app = tau + 1.5 * Tsw;
ia_s = I * sin(w * tau + beta);
ia_r = I * sin(w * t_app + beta);
dia_r = I * w * cos(w * t_app + beta);
u_ab = C.U_ac * sin(w * t_app) - C.R_ac * ia_r - C.L_ac * dia_r - C.L_ac * C.w_ci * (ia_s - i_a);

if C.apd
    % 4. capacitor-voltage reference that absorbs the double-line grid-side power
    X = C.U_ac - (C.R_ac + 1i * w * C.L_ac) * I * exp(1i * beta);
    Y = I * exp(1i * beta);
    A = abs(X * Y) / 2;
    alpha = angle(X * Y) - pi / 2;
    V = sqrt(2 * A / (w * C.C_ac * (1 - w^2 * C.L_c * C.C_ac)));
    th = alpha / 2;
    % 5-6. capacitor-voltage and capacitor-current loops
    uc_s = V * sin(w * tau + th);
    ic_ref = -C.C_ac * w * V * cos(w * tau + th) - C.C_ac * C.w_cv * (uc_s - u_cac);
    dic_ref = C.C_ac * w^2 * V * sin(w * t_app + th);
    u_cb = u_cac - C.R_c * i_c - C.L_c * (dic_ref + C.w_cc * (ic_ref - i_c));
    % 7. phase references (R22 eq. 5)
    u = [(2 * u_ab - u_cb) / 3; (-u_ab - u_cb) / 3; (-u_ab + 2 * u_cb) / 3];
    % 8. zero sequence
    im = [i_a; -i_a - i_c; i_c];
    if C.scheme == 1                                   % SVPWM
        u0 = -(max(u) + min(u)) / 2;
    else
        if C.scheme == 2                               % DPWM1
            [~, k] = max(abs(u));
        else                                           % CB-GDPWM
            [~, m] = min(abs(u)); [~, n] = max(abs(im));
            if m ~= n
                k = n;
            else
                if n == 1, p = 2; q = 3; elseif n == 2, p = 1; q = 3; else, p = 1; q = 2; end
                if abs(im(p)) > abs(im(q)), k = p; else, k = q; end
            end
        end
        u0 = sign(u(k)) * U / 2 - u(k);
    end
else
    u = [u_ab / 2; -u_ab / 2; 0];
    u0 = 0;
end
duty = min(max(0.5 + (u + u0) / U, 0), 1);

U_avg = (intU - intU_prev) / Tsw;      % average over the period that just ended
intU_prev = intU;
lg = [tau, i_a, i_c, u_cac, U, U_avg, dprev(1), dprev(2), dprev(3)];
dprev = duty;
end
