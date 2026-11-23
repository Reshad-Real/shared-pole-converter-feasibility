# Discriminating experiment: equal stored energy, different AC phase

This folder implements the first discriminating experiment of the study's research
plan (the planning documents are not part of this release). It looks for two operating states with
the same buffer energy and the same steady-state port powers whose transient admissibility
differs because of the AC phase. It then checks each state with energy sizing,
instantaneous constraints and a full trajectory test.

## Run

```bash
py sim/run_pair_experiment.py
```

```bash
py sim/run_sensitivity.py
```

```bash
py sim/run_switched_check.py
```

**Realistic ratings (WP1).** Rerun the pair experiment and the sensitivity sweep for a rating defined in `ratings.py`, then write the WP1 report:

```bash
py sim/run_pair_experiment.py --rating r120
```

```bash
py sim/run_sensitivity.py --rating r120
```

```bash
py sim/wp1_report.py
```

- `--rating` accepts `r120` (120 V, 60 Hz on 425 V) or `r230` (230 V, 50 Hz on 800 V). Results go to `results/r120/` and `results/r230/`.
- `--regime fast` on the pair experiment uses the fast DC1 response instead of the bus-loop one.
- The WP1 report needs both ratings and writes `results/wp1_ratings.md` and figure 9.

**Robustness campaign (WP4).** The criteria are in `WP4_CRITERIA.md`, fixed before any run. Run the four parts for each rating after its pair experiment (`bounds` first, because `nonideal` uses the tolerance corners it picks), then write the report:

```bash
py -u sim/run_wp4.py --rating r120 --part bounds events mc nonideal
```

```bash
py -u sim/run_wp4.py --rating r230 --part bounds events mc nonideal
```

```bash
py sim/wp4_report.py
```

- `bounds`: robust phase sweep over the tolerance set (R2, R3) and the θb certificate (R1).
- `events`: step-down, regenerative and ramped step-up sweeps.
- `mc`: 300 random cases per rating scored against the truth (R5).
- `nonideal`: switched-model pair runs under 24 grid, switching, sensor and tolerance conditions (R6).
- Post hoc, run after the verdicts to explain them: `diagnostics` (phase-blind limits per corner plant, one-at-a-time tolerance sensitivity, the R4 cases), `conditional` (R1–R3 with the grid amplitude measured) and `deadtime40` (the 40 kHz pair with a 0.25 µs dead time). Run them after `nonideal`.
- The circuit outer bound (follow-up to R4), run after `conditional` and `deadtime40`:
  - `r4circuit`: rechecks every run that held its limits;
  - `circuit`: phase sweep, θb margin and measured-grid analysis with the circuit bound (about 1.5 hours per rating);
  - `mccircuit`: rescores R5;
  - `circuitcheck`: consistency with a small inductor and with more cuts.
- Each part writes `results/<rating>/wp4_<part>.json`. The report writes `results/wp4_robustness.md`, `wp4_robustness.json` and figures 11–13. The whole campaign takes about an hour on 8 worker processes.

The first two take about two minutes each on a multi-core laptop; the switched check takes
about 15 seconds and reads the pair from `results.json`, so run it after the pair experiment.
Outputs go to `sim/results/`:

- The pair experiment writes `summary.md` (tables), `results.json` (every number) and
  figures 1–3.
- The sensitivity sweep writes `sensitivity.md`, `sensitivity.json` and figure 4. It shows
  how the effect changes with DC1 slew rate, buffer oversizing and modulation index.
- The switched check writes `switched_check.md`, `switched_check.json` and figures 5–6.

## Results at a glance

Full tables are in `results/summary.md` and `results/sensitivity.md`.

- **The pair exists.** AC2 load step of 210 W arriving at θa = 81° or θb = 261°. Both states
  have identical buffer energy (1.406 J), buffer current, bus voltage, DC1 power and AC power.
  - Energy sizing and the instantaneous/steady-state checks admit the step at both phases.
  - The trajectory test certifies θa feasible and proves θb infeasible in the reduced model.
    θb would need the pole-N band relaxed by 12.4 V.
- **Design consequence.** Energy sizing says 221 µF is enough at both phases. θb needs at
  least 312 µF, which is 1.41× more.
- **Independent check agrees.** At 24 phases, a conventional controller (tuned in the
  richer time-domain model) never succeeded where the LP proved infeasibility. It failed at
  3 phases where the LP found a trajectory; those failures are controller limits.
  - Tracking the LP plan at θa held every limit, with 3.2 V of pole-N headroom.
  - At θb the time-domain run tripped at 0.7 ms, slightly earlier than the LP's unavoidable
    violation near 1.5 ms. The buffer inductor needs voltage headroom to change its current,
    which the reduced model ignores.
- **Where it matters.** At the worst phase the phase-blind checks over-admit the step by at
  least:
  - 36–43 % when DC1 slews at 25–100 W/ms, falling to 16 % at 200 W/ms and 6 % at
    500 W/ms. At 1000 W/ms they become conservative at every phase.
  - 17 % to 58 % as the modulation index rises from 0.30 to 0.45.

  So the effect needs a DC source that is slow relative to the line cycle.
- **It survives switching.** The switched model uses a 20 kHz carrier, dead time 0, 1 and
  2 µs, diode freewheeling, real AC current loops and a one-period control delay.
  - At every dead time θa holds all limits (conventional controller and LP plan) and θb
    trips at 0.8–1.0 ms, against 0.7 ms in the averaged model.
  - At 24 phases the switched outcome matches the averaged one at 22. It never succeeds where
    the LP proved infeasibility. The two differences are a controller limit (135°) and an
    LP-unresolved phase (285°).
  - Two implementation costs appear that the averaged model hides:
    - The sampled, delayed control lets the bus swing about ±6 V instead of ±1.5 V.
    - Switching ripple lifts the θa instantaneous pole current to 20.3 A, just over the
      20 A limit (period average 18.1 A).

**Robustness (WP4, 120 V and 230 V ratings).** Full report: `results/wp4_robustness.md`.

- The predeclared criteria R1–R4 and R6 fail. Over the full tolerance box (grid amplitude ±10 %, Cb and Cdc ±20 %, DC1 response ±30 %), the most favourable plant admits nearly as much as the phase-blind checks: 8–9 % over-admission against the 15 % mark.
- The effect is still present in every plant (26–65 % computed plant by plant). Under random tolerances the nominal trajectory bounds admit 2–4 % of infeasible cases, against 20–21 % for energy sizing.
- With the grid amplitude measured (post hoc), the effect survives component and DC1 tolerances at nominal and high grid voltage (27–52 %) but not in a 10 % sag.
- In 3 of 242 successful switched runs, the buffer inductor's volt-seconds let the circuit succeed where the reduced-model outer LP says no trajectory exists.
- **Follow-up:** `outer_circuit.py` adds the buffer inductor to the outer bound. None of the 246 runs that held their limits contradicts it. θb stays provably infeasible with 17.7 V and 34.2 V of margin (32.5 V and 58.5 V before). The nominal G1 over-admission is 37 % and 38 % (45 % and 46 % before).

**Sizing (WP5a).** The protocol is `WP5A_PROTOCOL.md`; the report is `results/wp5a/wp5a_sizing.md`, with figures 14–15.

```bash
py -u sim/run_wp5a.py --part grid extras verify
```

```bash
py sim/wp5a_report.py
```

- The grid takes about 55 minutes on 8 worker processes.
- Post-hoc parts (`suff_scan`, `drop_reduced`) take about a minute together.

Results:
- **Base design, centred operating point:** the necessary capacitance exceeds phase-blind sizing by 73 % (120 V) and 74 % (230 V).
- **Operating point:** a higher pre-step operating point cuts the requirement by about a third, and the penalty remains.
- **The sufficient side** comes from a current-aware bound (`inner_circuit.py`, post hoc):
  - It certifies 36 of the 56 main-chart points, at 1.0–1.3 times the necessary value for m ≤ 0.35.
  - At m ≥ 0.40 with slow DC1 ramps the centred operating point is pole-current limited.
  - **S4(a) passes:** the certified plans at the top operating point and at m = 0.35 hold in the switched model with the smoothed-bus tracker (policy `plan_smooth`, gains `sync_ib`, `bus_smooth_n = 15`), with 8–17 V of bus margin.
  - **Chart-wide check:** 34 of the 36 certified plans hold, the thinnest by 2.5 V. The two failures are m = 0.45 with a 0.5 ms ramp, the feasibility edge.
  - The parts are `suff_circuit`, `verify2`, `alpha1`, `verify_chart` and `verify_conv`. The chart takes about 100 minutes, the checks about 15 minutes.

**System case study (WP6).** The protocol is `WP6_PROTOCOL.md`; the report is `results/wp6/wp6_case_study.md`, with figure 16.

```bash
py -u sim/run_wp6.py --part d4 events
```

```bash
py sim/wp6_report.py
```

- The D4 search takes about 11 minutes and the 2304 events about 110 minutes on 8 workers.
- **Results:**
  - At conventional sizing, an AC2 event that no controller can handle occurs on 98.5 % of days (stylised model, 60 events a day).
  - The WP5 step-up sizing leaves only the full step from zero load, which is beyond its design step.
  - Full load drops need more than twice the conventional capacitance with a bus-loop DC1.
  - The certified controller fails on drops to zero load at 120 V: the plans are certified but trip when tracked (bus undervoltage near the negative line peak), and the cause is unresolved.

**Digital-twin validation (WP9).** The protocol is `WP9_TWIN_PROTOCOL.md` (criteria fixed before any run; revision 1 after the method test). The physical converter is emulated by the switched model with the parameters of each WP4 misclassification case hidden from the twin. The twin identifies Cb, Lb, the port voltage amplitudes and the DC1 slew from two recordings of normal operation, measured with 12-bit, noisy sensors with gain errors (`twin.py`), and the bounds are evaluated on it. Run after the WP4 parts `mc` and `mccircuit`:

```bash
py -u sim/run_wp9.py --rating r120
```

```bash
py -u sim/run_wp9.py --rating r230
```

```bash
py sim/wp9_report.py
```

- Each rating takes about 40 minutes on 12 workers and writes `results/<rating>/wp9_twin.json`; the report writes `results/wp9/wp9_twin.md` and `wp9_summary.json`.

**Paper figures.** `run_paper_sweep.py` recomputes the phase boundary with the circuit-level bounds (`results/paper/phase_sweep_circuit.json`). `paper/make_figures.py` draws every data figure of the paper at print size from the result files. The sweep and its descending-scan pass (`scan`) take about an hour each on 8 workers.

```bash
py -u sim/run_paper_sweep.py
```

```bash
py -u sim/run_paper_sweep.py scan
```

## Circuit and scenario

This is the illustrative candidate from plan section 5. It is not a verified build-ready
circuit.

- Poles A, B and N share one DC bus. AC1 sits between A and N, AC2 between B and N. The
  ripple buffer (Lb, Rb, Cb) runs from pole N to the negative rail.
- DC1 feeds the bus through a slew-limited interface (a battery converter, say). DC2 is a
  constant bus load.
- The averaged pole limits give the pole-N band
  `max(0, -u1, -u2) <= vN <= Vdc - max(0, u1, u2)`.
- With the buffer inductor neglected, the buffer voltage must stay inside that band, which
  moves with the AC phase.
- The test event is an AC2 load step that arrives at AC1 phase θ. The AC current loops are
  taken as ideal, so the requested AC power is always delivered unless the band is empty.

All parameters live in `Params` in `model.py`.

## Why θ and θ + 180° make a clean pair

With both AC ports at the same frequency, every power and energy quantity repeats every
half cycle. A state at θ and a state at θ + 180° on the nominal orbit therefore have the
same buffer energy, buffer current, bus voltage, DC1 power, instantaneous AC power and
command. Only the sign of the AC voltages and currents differs.

Any criterion built from energy and power alone has to give both states the same verdict.
So if the trajectory test separates them, every such criterion misclassifies one of the two.

## Criteria compared

| Criterion | What it checks | Where |
|---|---|---|
| Energy sizing | Same dynamics, buffer energy anywhere in the capacitor rating window `[0, Vdc_max]` | `feasibility.py`, mode `energy` |
| Instantaneous and steady-state checks | AC commands synthesisable along the horizon; buffer voltage inside the band at the step instant; pre- and post-step periodic orbits exist | `run_pair_experiment.py`, `instantaneous_check` |
| Static derating | Buffer energy inside the worst-case intersection of the band over the cycle | mode `static` |
| Trajectory test, inner | Linear sufficient set (tangent/secant bounds on the concave `V(E)` maps, 3 V control headroom). A feasible LP gives a certified trajectory of the reduced model | mode `inner` |
| Trajectory test, outer | Linear necessary set. An infeasible LP proves that no controller can meet the limits in the reduced model | mode `outer` |

Every LP ends with one line cycle at constant DC1 power and periodic buffer energy. Repeating
that cycle keeps the trajectory admissible forever, so "feasible" means the converter
reaches a sustainable post-step state, not just survives the horizon.

## What a result means, and what it does not

- **Proof scope.** An infeasible outer LP is a proof for the reduced averaged model: buffer
  inductor neglected, lossless buffer branch, constraints enforced at 100 µs samples. It is
  not yet a proof for the switched circuit.
- **Independent check.** `timedomain.py` re-runs the pair in a richer averaged model, with the
  buffer inductor and resistance, the nonlinear bus, a finite-bandwidth buffer current loop
  and continuous-time pole limits. It is evidence, not proof.
- **Controller vs physics.** A conventional supervisor, tuned over 24 phases, is compared
  against tracking the LP plan. A conventional failure where the LP finds a trajectory is a
  controller limit; the LP's infeasibility at θb is a physical limit of the model.
- **Switched check.** `switched.py` adds PWM, dead time, diode freewheeling, device drops,
  real AC current loops and a one-period control delay. It confirms the verdicts but proves
  nothing either.
  - Both models are judged on period-averaged bus voltage and pole currents.
  - Instantaneous extremes, including switching ripple, are reported separately.
- **Not yet covered.** (WP4 added minimum pulse width, sensor gain and offset, and ADC quantisation.)
  - Switching transients and parasitic capacitance.
  - The common-mode path of real AC sources, and hardware results.
  - In the reduced LP, the buffer inductor; `outer_circuit.py` covers it. Both LPs still take the AC currents as delivered exactly.
  - DC1 and DC2 are still ideal power ports.
  - Pole currents are checked after the fact, not constrained in the LP.

## Files

- `model.py`: parameters, AC waveforms, pole-N band, steady-orbit interval test.
- `feasibility.py`: the four LP criteria in energy coordinates (scipy HiGHS).
- `timedomain.py`: averaged RK4 simulation with conventional and plan-tracking supervisors.
- `run_pair_experiment.py`: sweep, pair selection, certificates, buffer sizing, figures and
  reports.
- `run_sensitivity.py`: one-at-a-time sweeps of DC1 slew, Cb oversizing and modulation index.
- `ratings.py`: the baseline and the two realistic ratings, buffer sizing rule and DC1 regimes.
- `wp1_report.py`: G1 evaluation across ratings, report and figure 9.
- `switched.py`: event-driven switched model of the three legs with dead time.
- `run_switched_check.py`: the pair and the 24-phase sweep in the switched model, compared
  with the averaged model and the LP verdicts.
- `WP4_CRITERIA.md`: the WP4 tolerance set, conditions and pass marks.
- `robust.py`: tolerance corners and draws, perturbed plants, and the 24 non-ideal conditions.
- `run_wp4.py`: the four WP4 campaign parts.
- `wp4_report.py`: R1–R6 evaluation, report and figures 11–13.
- `outer_circuit.py`: outer bound for the circuit with its buffer inductor (pole-voltage band, exact sampled LC dynamics, energy cuts, bounded losses). `solve_joint` shares one pre-step orbit between several events; its free start is vacuous (see `WP5A_PROTOCOL.md`).
- `WP5A_PROTOCOL.md`, `sizing.py`, `run_wp5a.py`, `wp5a_report.py`: the WP5a sizing study.
- `inner_circuit.py`: current-aware sufficient bound for the circuit. It solves by sequential LP in tubes, and `extend_plan` continues a certified plan for simulation. `switched.py` gained:
  - `pre_orbit` (hold a given pre-step orbit);
  - the `circuit_plan` and `state_plan` policies (both tracked worse);
  - `plan_flat` (the energy tracker with the bus held at nominal) and `plan_smooth` (on a smoothed copy of the plan's bus);
  - the `sync_ib` gain, which compares the buffer current with its reference at the measurement instant. `plan_smooth` with `sync_ib` tracks the certified plans best;
  - `post_orbit`, the orbit the conventional controller recovers to after the step (WP6 operating-point schedule).
- `WP6_PROTOCOL.md`, `wp6.py`, `run_wp6.py`, `wp6_report.py`: the WP6 case study.
- `run_paper_sweep.py`: the paper's phase sweep with the circuit-level bounds.
- `WP9_TWIN_PROTOCOL.md`, `twin.py`, `run_wp9.py`, `wp9_report.py`: the WP9 digital-twin validation. `twin.py` holds the measurement model, the identification (least squares on windowed integral forms of the circuit equations) and the drift monitor.

  `feasibility.solve` also gained an optional pole-current limit for the inner bound (`i_limit`), used only in a post-hoc test.

`model.py` also carries the WP4 port departures (voltage harmonics, phase jump, AC2 power
factor and current harmonics), and `switched.py` can run a plant that differs from the
controller's model, with sensor gain, offset, ADC quantisation and a minimum pulse width.
With the defaults, both reproduce the earlier results exactly.
