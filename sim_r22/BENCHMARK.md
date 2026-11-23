# WP3 benchmark: R22 three-leg decoupling converter

Fixed on 22 October 2026, before digitising the measurements or running either model.

**Source.** J. Xu, T. B. Soeiro, F. Gao, H. Tang and P. Bauer, "Carrier-Based Generalized Discontinuous PWM Strategy for Single-Phase Three-Legs Active Power Decoupling Converters," IEEE TIE 68(11):11609–11613, 2021. Local copy in `literature/`; it is read as data only.

## Circuit and operating point (from the paper)

- **Topology:**
  - Three legs a, b and c on a 400 V bus, C_dc = 135 µF.
  - The grid, 220 V rms and 50 Hz, connects through L_ac = 1.44 mH between legs a and b.
  - The decoupling capacitor C_ac = 130 µF connects in series with L_c = 1.15 mH between legs c and b.
  - *Corrected 22 October 2026, after the first runs:* this line first read 135 µF, the DC-link value. R22 designs C_ac = 130 µF with its eq. (9) (p. 11611), as WP0 found. Targets and pass marks are unchanged; the first results are kept in `results/archive_cac135/`.
- **Switching:** SiC MOSFETs C3M0120090J at f_s = 40 kHz.
- **Load and phase angle:** rectifier mode into a resistive DC load. Measured cases:
  - φ = 0°: 80 Ω (about 2 kW);
  - φ = 15°: 83 Ω;
  - φ = 45°: 113 Ω;
  - STATCOM mode at 5 A and 10 A rms.
- **Control (Fig. 1):**
  - A DC-voltage loop sets the grid-current amplitude, applied at the set angle φ.
  - A grid-current loop produces u_ab*.
  - A capacitor-voltage loop with an inner capacitor-current loop produces u_cb*. The voltage reference is a line-frequency sine sized to absorb the double-line ripple (eqs. 1–4).
  - The modulator adds a zero-sequence signal: SVPWM, DPWM1, or CB-GDPWM (Fig. 2).

## Not given in the paper, so assumed (with sensitivity runs)

- Controller gains and structure details. The loops are tuned for clearly separated bandwidths: about 10 Hz for the DC voltage, 2 kHz for the currents, 300 Hz for the capacitor voltage.
- Dead time, set to 200 ns with variants of 100 and 300 ns. Sampling once per carrier period with a one-period delay.
- Device on-resistance of 0.15 Ω and diode drop of 3 V (SiC body diode). Ideal grid and ideal PLL.

## Measured targets

From Table I:
- **T1:** grid-current THD at φ = 0°: SVPWM 2.19 %, DPWM1 2.92 %, GDPWM 2.74 %.
- **T2:** grid-current THD at φ = 15° and 45° for the three schemes, as an extended check.

Digitised from the Fig. 7 oscilloscope traces:
- **T3:** DC-link ripple, peak to peak, before the decoupling is enabled (Fig. 7(a) and (b) main windows).
- **T4:** DC-link ripple, peak to peak, with decoupling enabled (zoom windows). The trace width limits this to an upper bound.
- **T5:** capacitor-voltage amplitude u_cac relative to the grid-voltage amplitude u_ac (same 500 V/div scale), and its phase relative to u_ac.
- **T6:** current amplitudes of i_a, i_b and i_c (20 A/div).
- **T7:** clamped duration of switch S_a per half cycle, CB-GDPWM (Fig. 7(a)) against DPWM1 (Fig. 7(b)).

Losses and efficiency (Table I ΔP) are excluded. The paper's own PLECS switching-loss differences between schemes (about 3 W, Fig. 5(a)) do not account for the measured differences (about 13 W), so they cannot be reproduced from the published information.

## How the traces are digitised

- The vertical scale is calibrated on the grid voltage: 311 V peak is 0.622 div at 500 V/div.
- Time comes from the 20 ms/div and 5 ms/div settings.
- Amplitudes come from least-squares sine fits at 50 Hz (100 Hz for the DC ripple).
- The digitisation uncertainty is the larger of ±2 pixels and the fit residual, converted to physical units.

## Pass criteria for each model (Python switched and Simscape)

| Target | Criterion |
|---|---|
| T1 | THD within ±1.5 percentage points of the measurement at each scheme, and the measured ordering (SVPWM lowest, DPWM1 highest) reproduced |
| T3, T5, T6 | Within ±10 % of the digitised value plus the digitisation uncertainty |
| T4 | Below the digitised upper bound |
| T7 | Clamped duration within ±1 ms, and GDPWM clamping longer than DPWM1 |

T2 is reported without a pass mark.
