# WP3: anchoring the tools to a published experiment (R22)

This folder reproduces measurements from R22 (Xu, Soeiro, Gao, Tang and Bauer, IEEE TIE 68(11):11609–11613, 2021), a 2 kVA single-phase three-leg active-power-decoupling converter. It uses both of the project's simulation tools:
- the event-driven Python switched model;
- Simscape Electrical.

The aim is to show that the simulation chain reproduces measured behaviour of a closely related circuit. R22 shares a switching leg between the grid and the decoupling capacitor, as our circuit shares pole N.

## Results (22 October 2026)

Full table: `results/wp3_r22_check.md`. Chart: `results/fig10_r22_benchmark.png`.

| | Python switched | Simscape |
|---|---|---|
| Waveform and power-flow targets (DC ripple with and without decoupling, capacitor-voltage amplitude and phase, three branch currents, S_a clamp durations): 15 checks each | 15 of 15 pass | 15 of 15 pass |
| THD of the grid current at φ = 0° (Table I): 3 checks each | 1 of 3 | 1 of 3 |

- **Close match on waveforms:**
  - Capacitor-voltage phase: −46.1° modelled against −46.5° measured.
  - Capacitor-voltage amplitude, GDPWM: 320 V against 325 V.
  - Capacitor current, GDPWM: 13.1 A against 13.1 A.
  - S_a clamp, GDPWM: 3.66 ms against 3.66 ms.
  - S_a clamp, DPWM1: 2.75 ms against 2.78 ms.
  - DC ripple without decoupling: 116 V against 123 ± 16 V.
- **THD is not reproduced:** both models give 0.8–1.3 % against the measured 2.2–2.9 %, and the measured ordering of the schemes is not recovered. R22's own PLECS simulation also shows 2–4 % with an ideal grid, so the measured THD depends mainly on unpublished controller details.
- **Losses and efficiency:** excluded in advance (see `BENCHMARK.md`).
- **Correction (22 October 2026):** the decoupling capacitor is now R22's 130 µF. The first runs used the DC-link value, 135 µF, in error; that run is kept in `results/archive_cac135/`.
  - All 30 waveform checks passed both times.
  - Moved closer to the measurements: capacitor-voltage amplitude 314 → 320 V (325 V measured), capacitor current 13.3 → 13.1 A (13.1 A), GDPWM clamp 3.75 → 3.66 ms (3.66 ms).
  - Moved slightly away, by 1 % or less: i_b amplitude and the DPWM1 clamp.
  - One THD check changed: Simscape's GDPWM THD fell from 1.34 % to 1.12 %, just outside its band. The total is 32 of 36 checks, against 33 before.

## Run

From the repository root:

```bash
py sim_r22/digitize_fig7.py
```

```bash
py sim_r22/run_python.py
```

```bash
matlab -batch "cd sim_r22/simscape; run_r22"
```

```bash
py sim_r22/compare_r22.py
```

1. The digitiser needs `pdftoppm` (Poppler, for example from MiKTeX) and the published article, which is not included in this release: obtain it from IEEE Xplore (DOI 10.1109/TIE.2020.3038099) and save it as `literature/Carrier-Based_Generalized_Discontinuous_PWM_Strategy_for_Single-Phase_Three-Legs_Active_Power_Decoupling_Converters.pdf`. The digitised traces it produces are included (`results/fig7_digitised.json`), so this step can be skipped.
2. The Python campaign (22 cases) takes under a minute.
3. The Simscape runs take about 10 minutes for 10 cases; `build_r22` runs automatically if the model is missing.

## Files

- `BENCHMARK.md`: targets and pass criteria, fixed before any measurement was digitised or any model run.
- `MODEL_SPEC.md`: circuit and control law that both implementations follow; assumed values are marked.
- `digitize_fig7.py`: Fig. 7 oscilloscope traces to numbers, calibrated on the known grid-voltage amplitude. `work/` holds the page render and an overlay used to check the extraction by eye.
- `r22_switched.py`, `run_python.py`: Python model and its campaign, including STATCOM and dead-time variants.
- `simscape/build_r22.m`, `ctrl_r22.m`, `fine_r22.m`, `run_r22.m`: the Simscape model, built by script. It reuses `sim_simscape/lib_paths.m` and `pwm_gates.m`.
- `r22_metrics.py`: identical metric code for both models.
- `compare_r22.py`: scoring, report and chart.

## Copyright

The paper in `literature/` is an IEEE copy licensed to the downloading institution. Do not include it, or `work/r22p-4.png` rendered from it, in a public code release; cite it instead.
