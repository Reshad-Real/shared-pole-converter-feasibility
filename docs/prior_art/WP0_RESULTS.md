# WP0: prior-art check — results

22 October 2026. Scored against the claims and decision rule fixed in `WP0_PROTOCOL.md` before any paper was read. Seven sources were read in full; three that could not be obtained were scored through legal stand-ins, and R13 could not be obtained at all.

## Decision

**Proceed.** No source gives C1 (admissibility depending on the line phase at which an event starts) with C2 (a trajectory criterion with a slew-limited DC source), so the pivot rule does not apply. No source gives C3 or C4 for a related architecture either, so nothing has to be narrowed to an extension.

Two claims are anticipated in part and must be worded accordingly:
- **C5 is weakly anticipated:**
  - R18 (Fig. 18) shows a slow PFC making a load step drain the buffer below its steady-state range.
  - R16 (§8.2) notes in prose that a slow buffer-voltage loop can push the buffer out of its operating region during sudden load changes.
  - Neither gives a criterion, a phase dependence or a bandwidth condition. Our paper should present the mechanism as known and claim its quantification.
- **"Phase dependence" in a steady-state or instantaneous sense is well established** (R2, R7, R16, R22 and the nine-switch literature). Every use of "phase" in our paper must mean the line phase at which the event starts, with equal stored energy, buffer current, port powers and power factor.

**Open:**
- R13 (ECCE 2026), the closest topical match, is not yet obtainable. G1's prior-art half is met for every obtainable source, on condition that R13 is checked when it reaches IEEE Xplore.
- R17, R20a and R20b were scored through stand-ins, not their own full texts.

## Scores

Y gives it, P partly, N does not. "Full text" means every page was read: by a subagent per paper, from extracted text and rendered pages, with key passages spot-checked by the main session.

| Source | Evidence basis | C1 | C2 | C3 | C4 | C5 |
|---|---|---|---|---|---|---|
| R2 Heller, Krismer, Kolar, OJPEL 2022 | Full text | N | N | N | N | N |
| R4 Liu, Sun, Su, IET PE 2016 | Review's earlier full-text notes; not re-read (servers unreachable) | N | N | N | N | N |
| R7 Liang et al., TIE 2022 | Full text | N | N | N | N | N |
| R13 ECCE 2026, paper 2424 | Not obtainable | — | — | — | — | — |
| R16 Hwang, thesis 2019 | Full text (103 pages) | N | N | N | N | P (weak) |
| R17 Krein, Balog, Mirjafari, TPEL 2012 | Stand-in: patents US 7,755,916 and US 8,325,499 (same derivation), plus the abstract | N | N | N | N | N |
| R18 Chen, Afridi, Perreault, TPEL 2013 | Full text (authors' manuscript) | N | N | N | N | P (weak) |
| R20a Kominami, Fujimoto, IAS 2007 | Stand-in: R20c | N | N | N | N | N |
| R20b Liu et al., TPEL 2009 | Stand-in: abstract (modulation, loss comparison) | N | N | N | N | N |
| R20c Kominami, Fujimoto, IEEJ 2008 | Full text (Japanese, read from rendered pages) | N | N | N | N | N |
| R21 da Costa et al., TIE 2021 | Full text | N | N | N | N | N |
| R22 Xu et al., TIE 2021 | Full text | N | N | N | N | N |
| Nine-switch phase-shift limits: Bradaschia et al., Eletrônica de Potência 2014; Wang et al., TPEL 2026 | Abstracts | N | N | N | N | N |

## What each source establishes, and how we cite it

| Source | Closest result (page or equation) | Difference to state |
|---|---|---|
| R2 | Sharing switching nodes between PFC and DAB makes transferable power depend on the instantaneous, line-determined duty cycles (eq. 11, Fig. 9a); a modulation makes the total limit phase-independent (eq. 27). Leg sharing costs degrees of freedom (p. 706) | Instantaneous, steady-state limits, worst-cased over phase; no event, buffer state or horizon |
| R4 | Capacitor design with modulation and switch-voltage constraints gives a feasible region (Section 4) | Steady-state design region; no transient |
| R7 | Submodule power limits with capacitor ripple, from a one-cycle periodic orbit (eqs. 8–20); error up to about 9 % depending on the power-factor angle (eq. 21, Fig. 5); insertion-voltage bounds of the same min/max form as our band (eqs. 4–7); a controller-independent steady-state "if and only if" viability condition; limits change only with the operating point (p. 5957) | "Phase" is the power-factor angle; constant-power DC source; no transient. Our C3 novelty rests on transient, multi-cycle outer bounds, not on controller independence alone |
| R16 | Instantaneous DCM/sliding-existence bounds on buffer voltage evaluated at worst-case line instants (eqs. 3.33, 3.35; Figs. 3.7–3.8); energy sizing (eqs. 2.4, 2.6); shared-switch APD losing authority at low buffer voltage (Ch. 7); a slow loop letting the buffer leave its region during sudden load changes (§8.2) | One AC and one DC port, buffer sharing a coupled inductor or switches, not a pole; load steps only simulated or measured; no criterion, no phase of event, no bound for all controllers |
| R17 | Minimum stored energy P₀/ω per half cycle, with a phase condition for the ripple port (patents, eqs. 3, 12–14) | Steady-state double-frequency ripple; ideal DC side. It is the energy-sizing baseline that C4 challenges |
| R18 | Half-cycle energy swing E_b = P_dc/ω_line (eq. 3) and energy buffering ratio (footnote 1); load step with a slow PFC drawing the buffer below its steady-state range (Fig. 18); stored energy fed back to the PFC loop (eqs. 5–10) | Buffer across one DC link, no shared pole; no admissibility condition, no phase dependence, no bound. The transient deficit itself is not new; our energy coordinates are an LP formulation, not a new state choice |
| R20a–c | Carrier-band split with an amplitude-proportional offset between two outputs sharing switches (R20c eqs. 7–13), equivalent to our pole band | Steady-state modulation range with a stiff DC source; no buffer, transient or event |
| Nine-switch line | Attainable amplitudes depend on the phase shift between the two outputs (Bradaschia 2014; linear range "phase-dependent" in Wang 2026) | Steady-state phase between ports, not the phase at which an event starts |
| R21 | Shared-leg pole-voltage condition: synchronised ports, V_l < V_g < v_C (eq. 8); load steps shown experimentally with a 30 Hz bus loop | No DC port, no buffer; steady-state condition only |
| R22 | Buffer amplitude and phase from the double-frequency power, depending on the power-factor angle (eqs. 1–4); zero-sequence clamping rule (eqs. 5–8, Fig. 2); capacitor sizing from steady-state ripple power (eq. 9); enable transient shown (Figs. 4, 7) | Steady-state sizing and instantaneous clamping; no admissibility analysis; "phase" is the power-factor angle |

## Wording rules for the paper

1. **The pole-voltage band and its free offset are background.** Cite R20c, R21, R22, and R7 for the same min/max form; claim only the transient, phase-dependent feasibility built on the band.
2. **Define "phase" once and keep it:** the AC1 line phase at which the event starts. Separately name and set aside the power-factor angle (R7, R22), the phase shift between ports (R21, nine-switch line) and instantaneous line-phase dependence (R2, R16).
3. **The transient energy deficit with a slow source is known qualitatively** (R18 Fig. 18, R16 §8.2). Claim the criterion, its certificates and the condition on the DC ramp, not the deficit.
4. **Do not claim energy coordinates as a new state;** R18 already feeds back stored energy.
5. **Controller-independent limits exist in steady state** (R7). Our proof claim is a transient, multi-cycle outer bound, including the buffer inductor (WP4 follow-up).
6. **Give one sentence on the "global stability" proofs of R16:** they assume unsaturated duties and do not contradict a bound for all controllers.
7. **Present energy sizing as stated in R17, R18 (eq. 3), R16 (eq. 2.6) and R22 (eq. 9)** as the conventional view that C4 tests.

## Still open

- **R13**, "Optimized Power Decoupling for a Three-Leg Single-Phase AC/AC Converter: Reduced DC Bus Ripple with Minimal Capacitance" (ECCE 2026, Vancouver, 4–8 October 2026). Check IEEE Xplore when the proceedings appear, or obtain it from the authors. Its title suggests steady-state capacitance minimisation, but that is not evidence.
- **R17, R20a and R20b full texts:** not accessible to us. The stand-ins cover their substance; obtain them through interlibrary loan or the authors before submission if possible.
- **Wang et al., TPEL 2026** (nine-switch maximum modulation index): its open full text is in the University of Pretoria repository; read it before submission, because it is recent and uses the word "phase-dependent".
- **R4:** re-read when its server or Wiley is reachable.

## Side finding for WP3

R22 states that its buffer capacitor C_ac is 130 µF, designed with its eq. (9) (p. 11611), and that 135 µF is the DC-link capacitor. The WP3 reproduction (`sim_r22/BENCHMARK.md`, `MODEL_SPEC.md`, both models) used 135 µF for C_ac as well. With the paper's formula for the capacitor-voltage amplitude, about √(1/C), the modelled amplitude should rise by about 2 %.

**Resolved the same day.** WP3 was rerun with 130 µF. The amplitude rose 1.9 %, and all 30 waveform checks still pass. The capacitor voltage, capacitor current and GDPWM clamp moved closer to the measurements; see `sim_r22/results/wp3_r22_check.md`.
