function [g, intU] = fine_r22(duty, U, tau, C)
%FINE_R22 Every solver step: PWM with dead time (leg c off when decoupling is off) and the
%   running integral of the DC voltage used for period averages.
%#codegen
persistent acc
if isempty(acc) || tau < 0.5 * C.Ts
    acc = 0;
end
g = pwm_gates(duty, tau, C);
if ~C.apd
    g(5:6) = 0;                       % leg c held off
end
intU = acc;
acc = acc + C.Ts * U;
end
