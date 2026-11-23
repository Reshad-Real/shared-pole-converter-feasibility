function [ints, inst] = fine_monitor(meas, vN, tau, C)
%FINE_MONITOR Running integrals and instantaneous extremes, every solver step Ts.
%   meas = [i1 i2 ib vb Vdc]; vN = pole-N voltage.
%   ints = integrals of [Vdc vb iA iB iN vN] (the controller differences them into period
%   averages); inst = [min Vdc, max Vdc, max |pole current|] after the step.
%#codegen
persistent acc vmin vmax ipk
if isempty(acc) || tau < 0.5 * C.Ts
    acc = zeros(6, 1); vmin = Inf; vmax = -Inf; ipk = 0;
end
i1 = meas(1); i2 = meas(2); ib = meas(3); vb = meas(4); V = meas(5);
iN = ib - i1 - i2;
ints = acc;                              % integral up to the start of this step
acc = acc + C.Ts * [V; vb; i1; i2; iN; vN];
if tau >= C.Tpre
    vmin = min(vmin, V); vmax = max(vmax, V);
    ipk = max([ipk, abs(i1), abs(i2), abs(iN)]);
end
inst = [vmin; vmax; ipk];
end
