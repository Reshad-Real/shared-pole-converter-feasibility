# WP4: robustness criteria

Fixed on 22 October 2026, before any WP4 run. The ratings are `r120` and `r230` from `ratings.py`, with DC1 in the bus-loop regime. The test pairs and worst phases are those found in WP1 (`results/r120/`, `results/r230/`).

## Tolerance set Θ

| Parameter | Range | In the corner set |
|---|---|---|
| Buffer capacitance Cb | ±20 % (film tolerance plus ageing) | yes |
| Bus capacitance Cdc | ±20 % | yes |
| Grid amplitude, both AC ports | ±10 % (sag and swell) | yes |
| DC1 response, through its slew rate | ±30 % | yes |
| Filter inductances L1, L2 and buffer inductance Lb | ±10 % each | random draws only |
| Buffer-branch series resistance Rb (stands for ESR drift) | 0.5–2 × | random draws only |

- **Corners:** the 16 corners of the first four parameters.
- **Draws:** uniform random samples over all eight parameters, log-uniform for Rb (fixed seed). The phase sweeps use the first 8 of the 64 R1 draws.
- **How a perturbed plant is analysed:** its LP is posed with its own parameters, starting from its own centred steady-state orbit (the state a controller tuned to that plant would hold). The port powers are commanded, so a sag raises the currents.
- **Designer's view:** the phase-blind checks (energy sizing, instantaneous checks, static derating) and the nominal LP bounds use nominal values, as a designer would.

## Definitions

- **Robustly infeasible:** the outer LP is infeasible at every corner and every draw.
  - This is a corner-and-sample check, not a proof over the continuous set; feasibility is not known to be monotone in grid amplitude.
- **Robust outer bound** at a phase: the largest outer-bound step over the corners and draws (16 corners and 8 draws).
- **Robust inner bound** at a phase: the smallest inner-bound step over the 16 corners, each with its own trajectory, so the bound is plant by plant.
  - This replaces the plan's stricter "one trajectory for every parameter value". That would need a shared terminal condition the plants cannot all meet, because their losses differ.
  - The practical one-trajectory question is tested instead in simulation: the nominal LP plan is tracked in perturbed plants.

## Criteria

| ID | Criterion | Pass mark |
|---|---|---|
| R1 | The infeasible member θb of the WP1 pair stays outer-infeasible at all 16 corners and all 64 random draws | all infeasible, at both ratings |
| R2 (H3) | Worst-phase over-admission, with the conventional limit at nominal parameters against the **robust** outer bound (9° grid) | ≥ 15 %, at both ratings |
| R3 (H3) | At that worst phase, robust outer bound minus robust inner bound | ≤ 10 % of the conventional limit |
| R4 | In every switched-model run of WP4, no run succeeds where the outer LP for that plant and condition, started from the run's own state at the step, proves infeasibility | zero contradictions |
| R5 | False-feasible and false-infeasible rates of each criterion over 300 random cases per rating | reported, no pass mark |
| R6 | Under every grid, switching and tolerance condition of the campaign below, θb still trips in the switched model, or its outer LP under that condition is infeasible | all conditions, at both ratings |

How R5 decides the truth:
- **Infeasible:** the outer LP of the actual plant is infeasible.
- **Feasible:** the inner-LP plan of the actual plant, tracked in the switched model of that plant, holds every limit.
- **Unresolved:** anything else.

## Campaigns

- **R5, misclassification:** 300 cases per rating, each drawing an AC1 phase, a step size up to the search limit, and a plant from Θ.
  - Every check is scored against the truth.
  - The switched model is the Python one, validated against Simscape in WP2 and against R22 in WP3. Simscape is about 50 times slower; a Simscape subset is a follow-up, not part of this pass mark.
- **R6, non-ideality:** for each condition, three runs on the pair: conventional at θa, plan at θa, conventional at θb. The reference is the nominal run with 1 µs dead time.
  - Grid and loads: the controller measures these, so its feedforward and its LP plan include them.
    - Voltage distortion of about 5 % (3rd 3 %, 5th 3 %, 7th 2 %) on both ports, phased in two ways: all harmonics peaking with the fundamental (raises the crest factor), or all opposing it (flat-topped).
    - Sag −10 % and swell +10 %; frequency −1 % and +1 %.
    - A ±20° phase jump of both port voltages at the step. The current references follow with a 5 ms first-order lag, and the plan is the one computed without the jump.
    - AC2 power factor 0.8, lagging and leading.
    - A nonlinear AC2 load: 3rd-harmonic current 30 % and 5th 15 % of the fundamental, peaking with it as from a capacitor-input rectifier.
  - Switching and measurement:
    - Dead time 0.5 and 2 µs; switching frequency 10 and 40 kHz.
    - Minimum pulse 1 µs (shorter pulses dropped or stretched, whichever is nearer).
    - 12-bit ADC quantisation on every measurement; two-period control delay.
    - Bus-voltage sensor gain −1 % and +1 %; buffer-current sensor offset 0.5 A.
    - When the switching frequency or the delay changes, the current-loop bandwidths are rescaled to keep the nominal delay phase lag at crossover (never above nominal). This is a design-time choice; the supervisory gains are unchanged.
  - Tolerances: the two corners of Θ with the smallest and the largest outer bound at θa. The controller measures the grid amplitude, so its feedforward and plan follow it, but keeps the nominal component values and DC1 model. The runs start with three line cycles of settling.
  - A run whose bus already leaves its limits in the last line cycle before the step is flagged as outside the design envelope; it says nothing about the step.
  - The LP is evaluated under every grid and load condition. Switching and measurement effects exist only in simulation; for those, θb must trip.
- **Events (LP phase sweeps, reported):**
  - step-downs and regenerative steps: AC2 power falling by up to twice its pre-step value, so beyond P2_pre it reverses;
  - step-ups ramped over 10 ms instead of 1 ms.
  - Worst-phase over-admission is reported for each, with bisection monotonicity spot-checked.
- **Also reported, no pass mark:** at the R2 phase, the plant-by-plant gap (largest inner–outer difference of a single corner), which separates the conservatism of the bounds from the spread caused by the tolerances.

## Not covered

- Device switching transients, parasitics and the common-mode path. These need a detailed device model at nanosecond steps.
- AC2 at a different frequency from AC1. The two ports then share no period, so the terminal repeating cycle of the LP is undefined.
- Controller retuning under perturbation. The controller keeps its nominal design.
