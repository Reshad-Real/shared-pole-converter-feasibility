# R22 reproduction: model and control specification

Both implementations follow this document: `r22_switched.py` (Python, event-driven) and `sim_simscape/ctrl_r22.m` with `build_r22.m` (Simscape Electrical). Values marked *assumed* are not given in R22.

## Circuit

- **Grid:** u_ac = 311.1 sin(ωt), 50 Hz, ideal.
- **Grid branch:** node a ← L_ac = 1.44 mH with 0.1 Ω (*assumed*) ← u_ac(+); u_ac(−) connects to node b.
- **Decoupling branch:** node c ← L_c = 1.15 mH with 0.1 Ω (*assumed*) ← C_ac(+), with C_ac = 130 µF (R22, eq. 9; first runs used 135 µF in error); C_ac(−) connects to node b.
- **DC side:** C_dc = 135 µF in parallel with the load resistance (80, 83 or 113 Ω; open in STATCOM mode).
- **Switches:** six switches, on-resistance 0.15 Ω (*assumed*), each with an antiparallel diode of 3 V forward drop (*assumed*, SiC body diode) and 0.15 Ω.
- **Sign conventions:** i_a flows from the grid into node a, i_c from the capacitor branch into node c, and i_b = −i_a − i_c into node b. The capacitor obeys C_ac du_cac/dt = −i_c.

## Modulation and timing

- **Carrier and sampling:** one symmetric carrier at 40 kHz. States are sampled at the start of each period; the duties take effect one period later.
- **Dead time:** 200 ns (*assumed*), with variants of 100 and 300 ns. No dead-time compensation.
- **Leg c with decoupling off:** both gates are held off and the branch carries no current.

## Control law, once per period (t_app = middle of the period the duties act in)

1. **DC voltage:** U_f is U_dc averaged over the last 10 ms. A PI controller with crossover 10 Hz (*assumed*) sets the grid-current amplitude I_d.
   - In STATCOM mode it adds a quadrature amplitude I_q of 5 √2 or 10 √2 A.
2. **Grid-current reference:** i_a* = I_d sin(ωt + φ) in rectifier mode, or I_d sin ωt + I_q sin(ωt + 90°) in STATCOM mode. φ is set as in each experiment.
3. **Grid-current loop:** u_ab* = u_ac − R_ac i_a* − L_ac di_a*/dt − L_ac ω_ci (i_a* − i_a), with ω_ci = 2π · 2 kHz (*assumed*).
4. **Capacitor-voltage reference:**
   - The double-line component of the grid-side power, p_ab = u_ab* i_a*, is computed in closed form as A sin(2ωt + α).
   - The capacitor branch must absorb it, which gives u_cac* = V sin(ωt + θ) with V = √(2A / (ω C_ac (1 − ω² L_c C_ac))) and θ = α/2.
   - This matches R22's eqs. 1–4, which give θ_cac = φ_ab / 2.
5. **Capacitor-voltage loop:** i_c* = −C_ac du_cac*/dt − C_ac ω_cv (u_cac* − u_cac), with ω_cv = 2π · 300 Hz (*assumed*).
6. **Capacitor-current loop:** u_cb* = u_cac − R_c i_c − L_c (di_c*/dt + ω_cc (i_c* − i_c)), with ω_cc = 2π · 2 kHz (*assumed*).
7. **Phase references (R22 eq. 5):**
   - u_a* = (2u_ab* − u_cb*)/3
   - u_b* = (−u_ab* − u_cb*)/3
   - u_c* = (−u_ab* + 2u_cb*)/3
8. **Zero sequence:** u0 is set by the modulation scheme, and duty d_x = clip(0.5 + (u_x* + u0)/U_dc, 0, 1).
   - *SVPWM:* u0 = −(max u* + min u*)/2.
   - *DPWM1:* clamp the phase with the largest |u*|: u0 = sign(u_k*) U_dc/2 − u_k*.
   - *CB-GDPWM (R22 Fig. 2):*
     - Let m be the phase with the smallest |u*| and n the phase with the largest measured |i|.
     - If m ≠ n, clamp n. Otherwise, clamp whichever of the other two phases carries the larger |i|.
     - u0 = sign(u_i*) U_dc/2 − u_i*.
9. **Decoupling off (target T3):** legs a and b only, with u_a* = u_ab*/2, u_b* = −u_ab*/2 and u0 = 0.

## Runs and measurements

- **Start:** each run starts from the computed steady state (U_dc = 400 V, references and branch states at t = 0). It runs 0.25 s; the last 0.1 s (five cycles) is analysed.
- **THD:** of i_a, harmonics 2–50 over the fundamental, as a power analyser would measure it.
- **Amplitudes and phases:** least-squares sine fits at 50 Hz.
- **DC ripple:** peak to peak of the DC voltage averaged over each switching period.
- **Clamp duration:** of S_a, the time per half cycle during which leg a's duty is 0 or 1.
