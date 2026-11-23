# WP5a: buffer sizing rule and design charts

Fixed on 23 October 2026, before the WP5a campaign. The question is how much buffer capacitance a design needs so that an AC2 load step is admissible at every AC phase, and how that compares with conventional, phase-blind sizing.

## Revision before the campaign

A first version of this protocol let the circuit outer bound start from any periodic pre-step state (a "free start"). A method test at the 120 V base design, rated step-up, showed two things:

1. **The free start is vacuous.** The bound fell to about 95 µF against 723 µF for conventional sizing. At the step, the solution claimed 9.9 J in the buffer capacitor where its voltage supports 4.4 J, plus 1.9 J in an inductor carrying almost no current. The chord relaxations tying stored energy to voltage and current let a free periodic start carry up to 7.8 J of energy that does not exist. Such a bound is valid but useless, so it is dropped.
2. **The operating point matters.** Every start must therefore be an exact state on a real steady-state orbit, and the operating point becomes a design variable (α below).

The same test computed three values at the 120 V base design with the centred operating point: Cb_conv = 723 µF, Cb_nec = 1255 µF (+74 %) and the reduced outer bound 1496 µF. The S1 pass mark below (15 %) was fixed in the plan before that test; the 230 V value was not seen.

## Design space

Both realistic ratings from `ratings.py`. The grid voltage is fixed; the designer chooses the bus voltage.

| Variable | Values | Base |
|---|---|---|
| Modulation index m = E / Vdc | 0.30, 0.35, 0.40, 0.45, 0.50 (Vdc = E / m) | 0.40 (425 V at 120 V; 812.5 V at 230 V) |
| DC1 ramp time T_r for the rated step | 0.5, 2, 5, 10 ms, the bus-loop value (13.3 ms at 120 V, 15.9 ms at 230 V), 20, 30 ms | bus-loop value |
| Operating point α (below) | 0.5, 0.75, 1.0 | 0.5 |
| AC2 power after the step, as a fraction of its rating | 50 %, 75 %, 100 % (from one third before the step) | 100 % |
| AC2 power factor | 1, 0.9 lagging, 0.8 lagging, 0.8 leading | 1 |

- **Main chart:** m × T_r at α = 0.5, rated step, unity power factor (35 points per rating).
- **Operating-point chart:** α = 0.75 and 1.0 along T_r at base m (14 more per rating).
- **One-at-a-time:** step size and power factor at the base point (5 more per rating).
- **Scaling with the bus voltage:** bus limits ±5 %, bus-capacitor stored energy kept (Cdc ∝ 1 / Vdc²), inner-bound headroom 1.5 % of Vdc. Inductors, pole-current limit and port powers are fixed.
- **DC1 slew:** the rated step divided by T_r.

## Operating point

Before the step the buffer runs on a steady-state orbit at the pre-step load. Its mean energy can sit anywhere in the interval that keeps the orbit inside the pole band with the inner headroom.

- **α:** the position in that interval. α = 0.5 is the centred orbit that the WP1–WP4 controller holds; α = 1 is the highest orbit, which favours step-ups.
- **Every start is exact:** the state on that orbit at the step phase.

## Capacitances per design point and α

1. **Conventional, Cb_conv (phase-blind):** the larger of two values.
   - The smallest Cb with a steady-state orbit at the post-step power inside the pole band, at nominal bus and with no headroom.
   - The largest over 24 phases (15° grid) of the smallest Cb for which the energy LP admits the step from the orbit at α. Its energy window is [0, Vdc_max].
2. **Necessary, Cb_nec (circuit outer bound, `outer_circuit.py`):** the largest over the checked phases of the smallest Cb for which the circuit outer LP is feasible from the orbit at α.
   - **Horizon:** three line cycles after the step, with no terminal condition.
   - **Phases checked:** 7.5° steps within ±22.5° of the worst phase of the reduced outer LP on the 24-phase grid. A subset of phases still gives a valid necessary value.
   - **What it means:** with that operating point, no controller can admit the step at every phase with less capacitance.
3. **Sufficient, Cb_suff (reduced model):** the largest over the 24 phases, plus 7.5° and 22.5° on each side of its worst phase, of the smallest Cb for which the inner LP is feasible from the orbit at α. It uses headroom, eight cycles and the terminal cycle. It is certified for the reduced model and checked in the switched model (S4).

Smallest-Cb searches are geometric bisections between 20 µF and 20 mF to 1 %, assuming feasibility rises with Cb (checked under S5).

## Quantities

- **Penalty** P(α) = Cb_nec(α) / Cb_conv(α) − 1: the extra capacitance, at the same operating point, that phase-blind sizing misses.
- **Operating-point lever:** Cb_nec(α) / Cb_nec(0.5), how much moving the operating point recovers.
- **Bracket:** Cb_suff / Cb_nec.

## Criteria

| ID | Criterion | Pass mark |
|---|---|---|
| S1 (H2, IJEPES) | Penalty at the base design, α = 0.5 | ≥ 15 % at both ratings (120 V seen in the method test) |
| S2 | Where in the m × T_r grid the penalty falls below 15 % | reported |
| S3 | **Chart transfer:** penalty from bilinear interpolation of the 120 V grid in (m, T_r · f), with T_r in line cycles, predicts the 230 V penalty at α = 0.5 | within ±10 points at every 230 V point inside the 120 V range |
| S4 | **Switched-model verification** at the base design, α = 0.5, worst phase, both ratings: (a) at Cb_suff the tracked inner-LP plan holds every limit; (b) at Cb_conv the conventional controller trips; (c) at 0.95 Cb_nec the conventional controller trips, since a success would contradict the bound | all, at both ratings |
| S5 | Monotonicity in Cb on 13-point checks at the base points; and three versus five post-step cycles for the circuit bound at the base points | reported |
| S6 | **Trade-off at the base design:** for α = 0.5, 0.75 and 1.0, the smallest Cb that also admits the drop of AC2 to zero load at every phase. Measured by the energy LP and the circuit outer LP; the latter is valid but loose for surplus, so this is indicative | reported |

- **A closed-form approximation** may be fitted after the runs to the 120 V results and tested on 230 V. It is labelled post hoc and reported with its error and range; it has no pass mark.
- **Not covered:** tolerances (WP4 covers them; add the tolerance margin on top), other topologies, and controller design.

## Addendum after the campaign: a current-aware sufficient bound

Written on 23 October 2026, after S4(a) failed: the reduced inner bound's plans exceeded the pole-current limit. This bound is post hoc and labelled so in the report.

- **What it is** (`inner_circuit.py`):
  - **Coordinates:** buffer voltage, buffer current and pole-N voltage, with the exact sampled LC dynamics.
  - **Pole-N band:** with the inner headroom, enforced at both ends of every sample.
  - **Pole-N current:** within κ = 0.97 of its limit, about 2 A of headroom for tracking.
  - **Buffer current:** changes by at most 36 A/ms, a rate a 2 kHz current loop can follow. The transient needs about 7 A/ms and the double-line ripple about 14 A/ms.
  - **Bus energy:** bracketed between worst-case low and high values, from chords and tangents on stored energy and on conduction losses. Both are kept 2 % of nominal voltage inside the bus limits.
  - **Why 2 %:** with 1 % the certified plan at 120 V, m = 0.35 dipped 4 V below the bus limit in the switched model. With 2 % the energy tracker still dipped 0.8 V below it. So the tracker changed as well (below).
- **How it is solved:**
  - Voltage and current are confined to tubes around a reference trajectory. The tubes start wide (100 V, 40 A) and narrow to 8 V and 2.5 A.
  - The reference is the outer bound's trajectory, smoothed, then refined by sequential LP. If the first solve is infeasible, the search restarts from the pre-step orbit continued over the horizon, which is smooth by construction.
  - The bound certifies a step only when a strict LP without slack is feasible. Among certified trajectories, it returns the one closest to nominal bus energy with the smoothest DC1 power and buffer current.
  - **Horizon:** four line cycles, ending in a cycle that repeats in buffer voltage and current at constant DC1 power. For simulation the plan is continued by repeating that cycle.
- **History:** a first chart run without the bus headroom and the current-rate limit gave plans that rode both bus limits and changed the buffer current by about 10 A per sample. The switched model could not track them, so that definition was superseded. Its results are kept in `results/wp5a/*_v1_superseded.json`.
  - A third run with 1 % bus headroom is kept in `results/wp5a/*_v3_bus1pct.json`.
  - A second run with the final definition but without the fallback reference was non-monotone in the ramp time. At large capacitance the first solve was infeasible and the search gave up. It is kept in `results/wp5a/*_v2_no_fallback.json`.
- **"Not certified" is not a proof of infeasibility.** The method is local, so it can miss a certificate. At the base design with the centred operating point, every attempt leaves a pole-N current overshoot of 1–4 A near the line peak, consistent with a real current limit.
- **Cb_suff,c on the main chart** (α = 0.5, m ≤ 0.45):
  - A geometric search at the necessary bound's worst phase, between max(Cb_nec, steady minimum) and three times that, to 4 %.
  - Certification at eight more phases 45° apart, raising Cb by 10 % while any phase fails.
- **S4(a) rerun:** the certified plan at the worst phase is tracked in the switched model.
- **Tracker.** All trackers follow the plan's DC1 power and stored energy. S4(a) is judged with the **smoothed-bus tracker**, policy `plan_smooth` with gains `sync_ib` and `bus_smooth_n = 15`:
  - **Smoothed bus:** it follows a 1.25 ms moving average of the plan's bus. That keeps the plan's 120 Hz use of the bus capacitor and drops its sample-to-sample steps. Fed forward unsmoothed (the energy tracker, `plan`), those steps saturate pole N. Holding the bus at nominal instead (the flat-bus tracker, `plan_flat`) fails where the plan needs the bus capacitor's share of the ripple.
  - **Timing (`sync_ib`):** it compares the buffer current with its reference at the measurement instant, as the AC loops do. The earlier loop compared the measured current with the reference for the period the duty acts in, 1.5 switching periods later. The buffer current then ran about 110 µs ahead of the ripple.
  - **Reported beside it:** the flat-bus and energy trackers.
  - **Variants tried and dropped as worse:**
    - following the plan's buffer current or its pole-N voltage;
    - feeding forward the buffer inductor's power;
    - extrapolating DC1 power to the duty's period;
    - a faster bus loop;
    - an integral bus term;
    - 5–10 % DC1 headroom in the bound.
- **Search floor.** The sufficient search starts at max(Cb_nec, the steady-state minimum with the inner headroom), not without headroom.
  - **Why:** at the steady-state minimum without headroom, the post-step orbit sits on the band edge, and fast-ramp plans could end only by swinging the bus.
  - **A final-cycle bus band was tried first.** It was infeasible, because by the final cycle the bound's bus-energy bracket spans 4–5 % of nominal at 120 V.
  - **Effect:** this changes only the 0.5 ms points of the chart.
- **Checks added** (reported, no pass mark):
  - every certified main-chart point at its worst phase, with every tracker (`--part verify_chart`);
  - S4(b) and (c) again with the timing fix (`--part verify_conv`): the conventional controller must still trip.
- **Archive:**
  - the 2 % runs with the energy tracker only: `results/wp5a/*_v4_energy_tracker.json`;
  - the runs with the flat-bus tracker and the old search floor: `results/wp5a/*_v5_flat_tracker.json`.
