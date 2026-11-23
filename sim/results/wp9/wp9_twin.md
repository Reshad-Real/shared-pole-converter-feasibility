# WP9: digital-twin validation against an emulated plant

Protocol: `sim/WP9_TWIN_PROTOCOL.md` (criteria fixed before any run; revision 1 after the method test, before the campaign). Plant: the switched model with hidden parameters from the WP4 tolerance set; the twin sees 12-bit, noisy measurements with gain errors. Evaluation events: the 300 WP4 misclassification cases per rating, unseen by the twin. **No result here is evidence about hardware.**

## Verdicts

| Criterion | 120 V | 230 V |
|---|---|---|
| TW1 identification | pass | pass |
| TW2 boundary | pass | pass |
| TW3 classification (H3) | pass | pass |
| TW4 closed loop | pass | pass |
| TW5 drift monitor | pass | pass |

## r120

Cases 300, errors 0, scored 300; truth {'feasible': 85, 'infeasible': 198, 'unresolved': 17}. Clamped identifications: 0; calibration runs that violated a limit: 0.

**M1 identification** (relative error)

| Parameter | median abs | 95th pct abs | max abs | mean | TW1 limit |
|---|---|---|---|---|---|
| Cb | 0.34 % | 0.83 % | 1.11 % | -0.15 % | 5.0 % |
| Lb | 0.40 % | 1.12 % | 2.00 % | -0.30 % | 2.5 % |
| E1 (grid amplitude) | 0.27 % | 0.47 % | 0.50 % | -0.00 % | 2.5 % |
| E2 (backup amplitude) | 0.25 % | 0.48 % | 0.50 % | -0.01 % | 2.5 % |
| DC1 slew | 0.76 % | 2.06 % | 4.33 % | -0.10 % | 7.5 % |
| Rb (reported) | 111.65 % | 419.21 % | 674.03 % | 38.92 % | reported |
| Cdc (reported) | 448.24 % | 2590.54 % | 2846.56 % | 768.78 % | reported |

**M2 boundary error** (predicted minus true largest step, / AC2 rating)

| Predictor | mean abs | 95th pct abs | bias | optimistic share | max abs |
|---|---|---|---|---|---|
| N | 5.36 % | 11.52 % | 0.40 % | 49 % | 16.84 % |
| NE | 4.54 % | 10.19 % | 0.05 % | 46 % | 14.62 % |
| T | 0.31 % | 0.89 % | -0.10 % | 2 % | 1.33 % |

**M3 classification** (false feasible / false infeasible)

| Bound | N | NE | T |
|---|---|---|---|
| circuit_outer | 4 (2.0 %) / 1 (1.2 %) | 3 (1.5 %) / 1 (1.2 %) | 0 (0.0 %) / 0 (0.0 %) |
| inner | 1 (0.5 %) / 4 (4.7 %) | 2 (1.0 %) / 5 (5.9 %) | 0 (0.0 %) / 0 (0.0 %) |

**M4 closed loop** (truth-feasible cases; admitted resolved cases)

| Predictor | held | violated | rejected | held share | admitted violated |
|---|---|---|---|---|---|
| NE | 80 | 0 | 5 | 94.1 % | 2 of 82 (2.4 %) |
| T | 84 | 1 | 0 | 98.8 % | 1 of 85 (1.2 %) |

**M5 drift** (aged capacitor: −15 % Cb, 2× Rb): alarms on aged plants 100.0 %, on unaged 0.0 %; median residual 0.23 V unaged, 1.14 V aged, threshold 0.69 V. Aged truth {'feasible': 80, 'infeasible': 210, 'unresolved': 10}.

| Bound | N | stale twin | recalibrated twin |
|---|---|---|---|
| circuit_outer | 11 / 0 | 11 / 0 | 0 / 0 |
| inner | 5 / 1 | 0 / 0 | 0 / 0 |

**M6 computation** (one core): identification 0.00 s, twin verdicts (outer + inner) 11.1 s, whole case 84 s (medians).


## r230

Cases 300, errors 0, scored 300; truth {'feasible': 89, 'infeasible': 186, 'unresolved': 25}. Clamped identifications: 0; calibration runs that violated a limit: 0.

**M1 identification** (relative error)

| Parameter | median abs | 95th pct abs | max abs | mean | TW1 limit |
|---|---|---|---|---|---|
| Cb | 0.28 % | 0.79 % | 1.13 % | -0.09 % | 5.0 % |
| Lb | 0.75 % | 2.19 % | 3.61 % | -0.11 % | 2.5 % |
| E1 (grid amplitude) | 0.23 % | 0.47 % | 0.50 % | 0.01 % | 2.5 % |
| E2 (backup amplitude) | 0.24 % | 0.47 % | 0.50 % | -0.01 % | 2.5 % |
| DC1 slew | 0.67 % | 2.07 % | 3.19 % | -0.07 % | 7.5 % |
| Rb (reported) | 102.35 % | 352.39 % | 554.31 % | 33.00 % | reported |
| Cdc (reported) | 212.98 % | 567.86 % | 101403339158715360.00 % | 431503570888422.00 % | reported |

**M2 boundary error** (predicted minus true largest step, / AC2 rating)

| Predictor | mean abs | 95th pct abs | bias | optimistic share | max abs |
|---|---|---|---|---|---|
| N | 5.25 % | 12.36 % | 0.31 % | 45 % | 17.66 % |
| NE | 4.28 % | 9.71 % | 0.12 % | 45 % | 12.80 % |
| T | 0.33 % | 0.88 % | -0.10 % | 3 % | 1.77 % |

**M3 classification** (false feasible / false infeasible)

| Bound | N | NE | T |
|---|---|---|---|
| circuit_outer | 9 (4.8 %) / 2 (2.2 %) | 4 (2.2 %) / 1 (1.1 %) | 1 (0.5 %) / 0 (0.0 %) |
| inner | 3 (1.6 %) / 5 (5.6 %) | 2 (1.1 %) / 6 (6.7 %) | 0 (0.0 %) / 0 (0.0 %) |

**M4 closed loop** (truth-feasible cases; admitted resolved cases)

| Predictor | held | violated | rejected | held share | admitted violated |
|---|---|---|---|---|---|
| NE | 83 | 0 | 6 | 93.3 % | 2 of 85 (2.4 %) |
| T | 89 | 0 | 0 | 100.0 % | 0 of 89 (0.0 %) |

**M5 drift** (aged capacitor: −15 % Cb, 2× Rb): alarms on aged plants 100.0 %, on unaged 0.0 %; median residual 0.41 V unaged, 1.81 V aged, threshold 1.22 V. Aged truth {'feasible': 79, 'infeasible': 201, 'unresolved': 20}.

| Bound | N | stale twin | recalibrated twin |
|---|---|---|---|
| circuit_outer | 18 / 0 | 17 / 0 | 1 / 2 |
| inner | 7 / 2 | 2 / 0 | 0 / 2 |

**M6 computation** (one core): identification 0.00 s, twin verdicts (outer + inner) 13.9 s, whole case 96 s (medians).

