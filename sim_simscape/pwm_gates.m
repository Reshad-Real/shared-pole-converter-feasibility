function g = pwm_gates(duty, tau, C)
%PWM_GATES Symmetric-carrier PWM with dead time (CONTROLLER_SPEC.md section 2).
%   Called every solver step Ts. Returns gate signals [A_hi A_lo B_hi B_lo N_hi N_lo].
%   A gate turns on only after its command has held for C.Nd steps (dead time td = Nd*Ts).
%#codegen
NBUF = 64;
persistent buf pos
if isempty(buf) || tau < 0.5 * C.Ts
    buf = false(3, NBUF); pos = 1;      % legs start low, as in the Python model
end
n = round(tau / C.Ts);
ph = mod(n, C.Nper) / C.Nper;           % position within the carrier period, 0..1
s = false(3, 1);
for k = 1:3
    d = duty(k);
    if d <= 1e-6
        s(k) = false;
    elseif d >= 1 - 1e-6
        s(k) = true;
    else
        s(k) = ph >= 0.5 * (1 - d) && ph < 0.5 * (1 + d);
    end
end
buf(:, pos) = s;
sd = buf(:, mod(pos - 1 - C.Nd, NBUF) + 1);   % command Nd steps ago
pos = mod(pos, NBUF) + 1;
g = zeros(6, 1);
for k = 1:3
    g(2 * k - 1) = double(s(k) && sd(k));        % upper gate
    g(2 * k) = double(~s(k) && ~sd(k));          % lower gate
end
end
