# WP2: independent Simscape Electrical model

This folder holds the commercial-tool check of the simulation-only plan (WP2 of the study plan). The same three-leg circuit is rebuilt in Simscape Electrical, with the controller re-implemented in Simulink MATLAB Function blocks from `CONTROLLER_SPEC.md`. The 81°/261° pair and the 24-phase sweep are then rerun and compared with the Python switched model (`sim/switched.py`) against predeclared criteria.

## Run

From the repository root:

```bash
py sim_simscape/export_scenarios.py
```

```bash
matlab -batch "cd sim_simscape; run_check"
```

```bash
py sim_simscape/compare_simscape.py
```

1. The export writes `scenarios.json`: initial states, orbit-derived references and the LP plan, shared by both models.
2. `run_check` builds `three_leg_switched.slx` if it is missing and runs 57 cases, about 20 s each.
3. The comparison writes `results/simscape_check.md`, `simscape_check.json` and figures 7–8.

## Results (22 October 2026)

Full report: `results/simscape_check.md`. Figure 7 overlays both models; figure 8 shows the switching-level mechanism at 261°.

| Predeclared criterion (dead time 1 µs) | Result |
|---|---|
| No Simscape run succeeds where the LP proves infeasibility | Met (also at 2 µs) |
| Verdicts match the Python switched model at ≥ 22 of 24 phases | Met: 23 of 24; the exception, 285°, is LP-unresolved |
| Trip times within 0.3 ms, buffer voltage within 2 V RMS | Met: within 0.10 ms and 1.65 V |

- At 81° and 261° the pair behaves the same in both models. 81° holds every limit under both the conventional controller and the LP plan. 261° trips at 0.8 ms in Simscape, against 0.9 ms in Python.
- At 2 µs dead time, five runs at LP-feasible phases trip in Simscape but not in Python. Both models put them within half a volt of the 190 V limit, and Simscape sits 0.2–0.6 V lower. Halving the Simscape step changes its result by 0.27 V and keeps the trip, so this is a genuine small model difference on a knife edge. The practical conclusion holds in both models: 2 µs of dead time leaves the tuned controller no bus margin at those phases.

## Model

| Part | Implementation |
|---|---|
| Switches | Six Ideal Semiconductor Switch blocks (50 mΩ), integral protection diodes (0.8 V, 50 mΩ) |
| Passive parts | Foundation inductors with series resistance, capacitors, floating AC voltage sources |
| DC1 and DC2 | One controlled current source, (pd − Pdc2) / Vdc, into the positive rail |
| Solver | Partitioning local solver, fixed step 0.1 µs (500 steps per 50 µs switching period) |
| Controller | `ctrl_step.m`, once per period, one-period delay to the PWM |
| PWM and dead time | `pwm_gates.m`, every step: symmetric carrier, turn-on delayed by the dead time |
| Monitor | `fine_monitor.m`: running integrals for period averages, instantaneous extremes |

Port order and polarity were established by test circuits (`test_ports.m`), not assumed:
- On the switch block, the first right-hand port is the drain (positive) terminal.
- A positive command into the controlled current source drives current out of its left port.

## What the comparison can and cannot show

- The circuit solver, device models, gate timing and controller code are independent of the Python model. The scenario data and the specification are shared, and the same author wrote both implementations.
- Agreement shows that the conclusions do not depend on one implementation. It is not experimental validation.

## Files

- `CONTROLLER_SPEC.md`: the written specification both implementations follow.
- `export_scenarios.py`: shared scenario data, written to `scenarios.json`.
- `build_model.m`, `lib_paths.m`: build the Simulink/Simscape model from script.
- `ctrl_step.m`, `pwm_gates.m`, `fine_monitor.m`: controller, modulator and monitor.
- `run_case.m`, `run_check.m`: single runs and the full campaign.
- `test_ports.m`: port-order and polarity tests.
- `compare_simscape.py`: criteria, report and figures.
