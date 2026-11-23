# WP9 protocol: digital-twin validation against an emulated plant

**Fixed 24 October 2026, before any WP9 run.** Revisions after the method test, and before the campaign, are recorded at the end of this file. Analyses defined after the verdicts are known will be labelled *post hoc*.

## Why this work package exists

On 22 October 2026 the twin contribution was withdrawn (`simulation_route_plan.md`, Decision), because a model without a physical counterpart and a live data link would read as overclaiming. On 24 October 2026 the authors asked for the digital-twin validation to be included. No hardware exists, so the twin is validated here against an **emulated plant**: the switched model with plant parameters hidden from the twin. This tests the twin's pipeline (measurement, identification, model update and prediction) and its effect on admissibility predictions. It does not test fidelity to hardware, and the paper must say so.

The validation follows contribution 3 and hypothesis H3 of `converter_literature_review_and_refined_plan.md` (section 7), and the metrics of its section 8:

> **H3:** Calibration reduces false-feasible predictions or unnecessary rejection at a matched confidence target relative to a fixed nominal model. Report both error directions; rejecting everything is not an improvement.

## Set-up

**Plant.** For each case, the plant is the rating's nominal parameter set scaled by the case's factors from the WP4 tolerance set: Cb and Cdc ±20 %, grid amplitude ±10 %, DC1 slew ±30 %, L1, L2 and Lb ±10 %, Rb 0.5–2 times. The plant is simulated by the switched model (`switched.py`) with 20 kHz PWM, 1 µs dead time, 50 mΩ on-resistance and 0.8 V diode drop. Its controller is the conventional one. It uses the measured grid amplitude, as a converter's PLL provides, and nominal values for everything else. Its own sensors have 12-bit ADCs.

**What the twin sees.** The twin uses only what the converter's controller has. It sees the samples taken at the start of each switching period of:

- buffer voltage $v_b$, bus voltage $V_{dc}$, buffer current $i_b$ and port currents $i_1$, $i_2$;
- port voltages $e_1$, $e_2$;
- DC1 power $p_d$ (reported by the battery converter) and DC2 power;
- the applied duties of the three legs.

Every measured channel has a 12-bit ADC, white noise of 1 LSB (standard deviation), and a gain error drawn once per converter, uniform in ±0.5 %. The channel ranges are:

| Channel | Range |
|---|---|
| Unipolar voltage ($v_b$, $V_{dc}$) | 0 to 1.25 $V_{max}$ |
| Port voltage ($e_1$, $e_2$) | ±0.625 $V_{max}$ |
| Current ($i_b$, $i_1$, $i_2$) | ±1.5 $I_{max}$ |
| Power ($p_d$, DC2) | ±$P_{dc1,max}$ |

The twin knows the nominal design values, the device data (on-resistance, diode drop, dead time) and the switching frequency.

**Calibration data.** These come from two recordings of the plant in normal operation. Neither is used for evaluation.

- **C1, steady operation.** Two line cycles at the pre-step load, after three cycles of settling.
- **C2, an earlier load event.** An AC2 step-up of one eighth of the port's rating (300 W at 120 V, 600 W at 230 V), at a phase drawn independently of the evaluation event, recorded for two cycles after the step.

**Identified parameters.** The twin identifies Cb, Lb and Rb from the buffer branch, Cdc from the bus balance, the grid amplitude from the port voltages, and the DC1 slew from the C2 response. It uses least squares on integral (windowed) forms of the circuit equations, with the pole voltages estimated from duties, the measured bus voltage, dead time and on-resistance. L1, L2, R1 and R2 keep their nominal values: their effect on the band is a few tenths of a volt (checked in the method test). Estimates outside a plausible range (0.5 to 2 times nominal; 0.25 to 4 times for Rb) are clamped, and the clamp is reported.

**Twin.** The twin is the averaged circuit model with the identified parameters. Its predictions are the bounds of the paper evaluated on it, from its own centred orbit, exactly as the nominal model is evaluated.

**Evaluation events.** These are the 300 random cases per rating of the WP4 misclassification campaign (`results/<rating>/wp4_mc.json`): random phase, step and plant. Their truth labels are those of Table 3 of the paper:

- *infeasible*: the plant's own reduced and circuit outer bounds are both infeasible;
- *feasible*: the plant's own inner-bound plan held every limit when tracked in the switched model by a controller that knew the plant;
- otherwise *unresolved*, and excluded.

The twin never sees an evaluation event before it predicts it.

**Predictors compared:**

- **N**, the nominal model: the Table 3 rows, all nominal values;
- **NE**, the nominal model with the measured grid amplitude;
- **T**, the twin.

## Metrics

- **M1 Identification.** The relative error of each identified parameter, over the 300 plants per rating.
- **M2 Boundary.** At each case's phase, the largest admissible AC2 step by the circuit outer bound, for N, NE, T and the true plant, by bisection to 0.5 % of the AC2 rating. Three post-step cycles are used, as in the sizing campaign. The error is the predictor's value minus the true plant's, normalised by the AC2 rating.
- **M3 Classification.** False-feasible and false-infeasible rates against the truth labels, for the circuit outer bound and the reduced inner bound (eight post-step cycles, as in Table 3), for N, NE and T.
- **M4 Closed loop.** For every resolved case that T admits (its inner plan exists), the plan is tracked on the plant by the energy tracker of the WP4 truth runs, with the twin's parameters as the controller's model, 12-bit controller sensors and three settling cycles. The same is done for NE (plan and controller model from NE). Outcomes are:
  - *held*;
  - *violated*: a limit is violated on the plant;
  - *rejected*: the predictor's inner plan does not exist.
- **M5 Drift.** Each plant's buffer capacitor is aged after calibration: −15 % capacitance and doubled Rb. The twin's one-line-cycle prediction residual of $v_b$ is evaluated on two cycles of steady operation of the aged plant, and of the unaged plant (C1-like data, recorded separately). The alarm threshold is three times the 95th percentile of the residual on the calibration data of the same plant.
  - **M5a:** the alarm rate on aged plants and the false-alarm rate on unaged ones.
  - **M5b:** M3 for the aged plants with the stale twin and with a recalibrated twin, on the same events. Truth for aged plants is re-established with the aged plant's own bounds and an oracle-tracked plan, as in WP4.
- **M6 Computation.** The time to identify and the time for one twin outer-bound verdict.

## Predeclared criteria

| | Criterion | Pass if, at both ratings |
|---|---|---|
| TW1 | Identification | The 95th percentile of the absolute relative error is at most a quarter of the parameter's tolerance half-width: 5 % for Cb and Cdc, 2.5 % for the grid amplitude and Lb, 7.5 % for the DC1 slew. Rb is reported only. |
| TW2 | Boundary | The twin's mean absolute boundary error (M2) is at most half that of NE. |
| TW3 | Classification (H3) | For both the circuit outer and the reduced inner bound, the twin's false-feasible rate and false-infeasible rate are each no higher than N's, and its total number of misclassified cases is lower than N's. |
| TW4 | Closed loop | At least 90 % of truth-feasible cases are admitted by the twin and held. At most 5 % of the twin's admitted, resolved cases violate a limit. |
| TW5 | Drift | The residual monitor alarms on at least 95 % of aged plants, with at most 5 % false alarms on unaged plants. |

M5b and M6 are reported without a criterion.

## Interpretation rules

- **Emulated plant.** All results hold for an emulated plant. The twin and the plant are different models (averaged and switched), but they were written by the same authors and share the circuit topology and parameterisation. No result here is evidence about hardware.
- **Twin as infrastructure.** If TW3 or TW4 fails while TW1 passes, the twin is reported as identification infrastructure, not as an improvement of admissibility predictions (H3, last sentence).
- **Truth by construction.** On infeasible cases, an oracle twin (the true parameters) agrees with the truth by construction. Its rows are therefore not reported as a predictor.

## Revisions

**Revision 1 (24 October 2026), after the method test and before the campaign.** The method test used four plants per rating drawn from seeds 9000 and 9001, which no evaluation case uses.

It found three problems with the identification as first written:

1. **Pole-voltage estimate near zero current.** The duty-based estimate of the period-average pole-N voltage is accurate (0.4 V RMS) where |i_N| ≥ 8 A. Near the current's zero crossings it is off by 4–9 V: the current changes sign within a period, so the dead-time correction has the wrong sign. With every window used, Lb was biased by +2 to +8 %.
2. **Cdc and Rb.** These were not practically identifiable from normal operation, with errors of +28 % to +711 % (Cdc) and −86 % to +270 % (Rb). The bus capacitor's current is well under 1 A, against tens of amperes of leg current, so the bus balance cannot reveal it; Rb's drop is comparable to the pole-voltage estimation error.
3. **One-line-cycle residual.** The residual of $v_b$ predicted from the sampled buffer current has a baseline of 2.7 V even on exact data, because the small difference between a sample and the period average accumulates over a cycle. A 15 % capacitance loss would not clear the 3× threshold reliably.

Changes:

- **Window selection.** The Lb and bus regressions use only windows in which the relevant pole currents stay above 2.5 times the nominal switching-ripple amplitude. On the method-test plants Lb then came within ±1.4 % (1.5 times gave −1.6 to +2.6 %; 3.5 times left too few windows).
- **Cdc, Rb, L1 and L2 at nominal.** The twin keeps the nominal Cdc and Rb, as already planned for L1 and L2. Rb and Cdc are still estimated and reported. On the method-test check, at three phases per rating, ±20 % Cdc, Rb from 0.5 to 2 times, or +10 % L1 and L2 moved the circuit-outer-bound boundary by at most one bisection step (0.5 % of the AC2 rating).
- **TW1 scope.** TW1 scores Cb, Lb, the grid amplitude and the DC1 slew; Cdc is reported only, like Rb.
- **Drift monitor (M5, TW5).** The statistic is now the RMS, over 8-period windows, of the buffer-capacitor equation error $(v_b[k+8]-v_b[k]) - \int i_b\,dt/\hat C_b$, on a two-cycle steady record. The alarm threshold is three times this statistic on the twin's own steady calibration record (C1). On the method-test plants the statistic was 0.21–0.44 V on calibration and on fresh unaged data, and 1.06–1.94 V on aged plants.

The criteria TW1–TW5 and their thresholds are otherwise unchanged.
