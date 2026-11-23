# Controller and test specification (three-leg converter)

This is the written specification that the Simulink controller in `ctrl_step.m` implements. It describes behaviour, not code. The Python switched model (`sim/switched.py`) was written earlier to the same behaviour; the two implementations share this document and the scenario data, nothing else.

## 1 Circuit

- Three half-bridge legs on one DC bus: A, B and N. Each switch has on-resistance Ron = 50 mΩ and an antiparallel diode (forward drop 0.8 V, 50 mΩ). The bus capacitor is Cdc = 20 µF.
- AC1 is the loop pole A → R1, L1 → source e1 → pole N; AC2 is pole B → R2, L2 → source e2 → pole N. The sources are e1 = E1 sin(ωt) and e2 = E2 sin(ωt + φ2), floating with respect to the bus.
- The buffer is pole N → Rb, Lb → Cb → negative rail.
- DC1 and DC2 act together as one controlled current into the positive rail: (pd − Pdc2) / Vdc. Here pd is the DC1 power chosen by the controller, and DC2 is a constant load Pdc2.
- Parameters are those of `Params` in `sim/model.py`: f = 50 Hz, E1 = E2 = 80 V, L1 = L2 = 2 mH with 0.1 Ω, Lb = 1 mH with 0.1 Ω, Cb = 250 µF, P1 = 150 W, AC2 100 W before the step, Pdc2 = 50 W, and a DC1 slew limit of 50 W/ms.

## 2 Modulation

- There is one symmetric carrier at fsw = 20 kHz (period Tsw = 50 µs). In period k, starting at tk, the upper switch of a leg with duty d is commanded on during [tk + (1 − d)Tsw/2, tk + (1 + d)Tsw/2). A duty within 1e-6 of 0 or 1 means the switch stays off or on for the whole period.
- Dead time td delays every turn-on: the upper gate is on when the command has been high for at least td, and the lower gate when it has been low for at least td. Default td = 1 µs; variants use 0 and 2 µs.
- States are sampled at tk, the middle of the low-side interval. The duties computed from that sample apply during the next period, a one-period delay. The DC1 power applies at once.

## 3 References and scenario

- The test is an AC2 load step of size ΔP. The absolute time t is chosen so that ωt is the AC1 phase. The simulation starts at t0 = θ0/ω and the step occurs one line cycle later, at t_step = t0 + T, so it arrives at AC1 phase θ0.
- AC current references at the port terminals, with unity power factor:
  - AC1: i1* = I1 sin(ωt), with I1 = 2P1/E1.
  - AC2: i2* = A2(t) sin(ωt + φ2). Its amplitude rises from 2P2/E2 to 2(P2 + ΔP)/E2 along a raised cosine over Tr2 = 1 ms starting at t_step.
- Converter-side reference voltages: u1* = e1 + L1 di1*/dt + R1 i1* and u2* = e2 + L2 di2*/dt + R2 i2*, using the exact derivative of the amplitude ramp. The reference AC pole power is pac* = u1* i1* + u2* i2*.
- From the Python orbit calculation, the scenario supplies:
  - the initial states;
  - P_pre and P_post, the average power DC1 must deliver before and after the step;
  - the mean buffer energy of the pre- and post-step orbits;
  - the pre-step buffer-energy history used to start the moving average.

## 4 Control law (executed once per period at tk)

Let x(t) = clip((t − t_step)/Tr2, 0, 1) and let t_app = tk + 1.5 Tsw be the middle of the period in which the new duties act.

1. **Energies.** Eb = Cb vb²/2 and Ed = Cdc Vdc²/2. Eb_avg is the mean of the last M = 200 sampled values of Eb (half a line cycle).
2. **Supervisor.** It sets the DC1 command pd_cmd, the bus-energy reference Ed_ref and its feedforward dEd_ff.
   - *Conventional:* Ed_ref = Cdc · 200² / 2 and dEd_ff = 0. P_cmd = P_pre + (P_post − P_pre)(1 − cos πx(tk))/2. Eb_ref is the pre-step mean before t_step and the post-step mean after it. Then pd_cmd = P_cmd + k_rec (Eb_ref − Eb_avg), with k_rec = 2π · 40 s⁻¹.
   - *LP plan, after the step:* Ed_ref and dEd_ff are the plan's bus energy and its slope, interpolated at t_app. pd_cmd is the plan's DC1 power for the current 100 µs interval plus k_tot times the gap between the plan's total energy and Eb + Ed at tk, with k_tot = 2π · 10 s⁻¹.
3. **DC1.** pd moves toward pd_cmd by at most 50 W/ms · Tsw per period and stays within ±1000 W.
4. **Buffer current reference.** pb_ref = pd − Pdc2 − pac*(t_app) − dEd_ff + k_bus (Ed − Ed_ref), with k_bus = 2π · 200 s⁻¹. Then ib_ref = pb_ref / max(vb, 1 V).
5. **Buffer leg voltage.** vN_cmd = vb + Rb ib + Lb (dib_ff + ωc (ib_ref − ib)), where dib_ff is the change of ib_ref since the previous period divided by Tsw, and ωc = 2π · 2000 rad/s.
6. **AC legs.** u1c = u1*(t_app) + L1 ωac (i1*(tk) − i1), and likewise u2c, with ωac = 2π · 1500 rad/s.
7. **Pole-N band.** lo = max(0, −u1c, −u2c) and hi = Vdc − max(0, u1c, u2c). If lo > hi the band is empty: vN = (lo + hi)/2 and the period counts as a synthesis failure. Otherwise vN = clip(vN_cmd, lo, hi).
8. **Duties.** dA = (u1c + vN)/Vdc, dB = (u2c + vN)/Vdc and dN = vN/Vdc. Each then receives dead-time compensation (td/Tsw) · clip(i_ref/0.3 A, −1, 1), where i_ref is i1*, i2* or ib_ref − i1* − i2* respectively, and is clipped to [0, 1].

## 5 Violation, trip and metrics

- After the step, a period violates the limits if any of these holds:
  - its average bus voltage is outside 190–210 V;
  - the largest period-average pole current so far exceeds 20 A;
  - the pole-N band is empty.
- The pole currents are iA = i1, iB = i2 and iN = ib − i1 − i2.
- A violation is time-stamped at the start of its period. The run stops 1 ms after the first violation (protective trip) or at the end of the horizon, eight line cycles after the step.
- Reported metrics:
  - period-average bus voltage extremes and instantaneous bus voltage extremes;
  - peak period-average and peak instantaneous pole currents;
  - total time pole N was clipped and the time its band was empty;
  - the smallest distance between the commanded pole-N voltage and the band edge;
  - whether any limit was violated, and the time of the first violation.
