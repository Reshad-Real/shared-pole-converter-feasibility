# Phase-dependent transient feasibility and buffer sizing of a shared-pole four-port AC/DC converter with digital-twin validation

Code, models and results for the article of the same title, submitted to the *International Journal of Electrical Power & Energy Systems* \[citation and DOI to be added on publication\].

In a converter whose energy buffer shares a switching pole with two AC ports, whether a load step can be survived depends on the AC phase at which it arrives, not only on the stored energy. This repository contains everything needed to reproduce the study:

- the averaged and switched converter models;
- the linear-programming feasibility bounds;
- the campaigns behind every figure and table, with their predeclared protocols and results.

Everything is simulation; no hardware was built. The digital twin of the study is validated against an emulated plant (the switched model with parameters hidden from the twin), not against hardware.

## Contents

| Folder | What it holds |
|---|---|
| `sim/` | Python models and campaigns: averaged model (`model.py`, `timedomain.py`), switched model (`switched.py`), feasibility bounds (`feasibility.py`, `outer_circuit.py`, `inner_circuit.py`), ratings and sizing (`ratings.py`, `sizing.py`), the digital twin (`twin.py`), campaign scripts (`run_*.py`), report generators (`*_report.py`), and the predeclared protocols (`WP4_CRITERIA.md`, `WP5A_PROTOCOL.md`, `WP6_PROTOCOL.md`, `WP9_TWIN_PROTOCOL.md`). `sim/README.md` documents every script and result in detail. |
| `sim/results/` | All results: JSON data, Markdown reports and figures, per study part (see below). |
| `sim_simscape/` | The independent Simscape Electrical model (MATLAB scripts, the `.slx` model, the controller specification) and its comparison with the Python switched model. |
| `sim_r22/` | Reproduction of a published three-leg decoupling experiment \[Xu et al., IEEE Trans. Ind. Electron. 68(11):11609–11613, 2021\]: the Python and Simscape models, the predeclared benchmark (`BENCHMARK.md`), digitised measurement traces and the comparison. |
| `paper/` | `make_figures.py`, which draws every data figure of the paper at print size from the result files. |
| `docs/prior_art/` | Protocol and results of the full-text prior-art check. |
| `MANIFEST.md` | Size and SHA-256 checksum of every file. |

The planning documents and the manuscript sources are not part of this release.

### Not included, and why

- **Published articles.** No article PDF is included. Several were read under institutional licences (IEEE Xplore), and none may be redistributed.
- **Renders of a published figure.** The page renders and extraction aids used to digitise the oscilloscope traces of Xu et al. (2021) are not included either.
- **What is included instead:** the digitised numerical traces (`sim_r22/results/fig7_digitised.json`), so the benchmark comparison runs without the article. To repeat the digitisation, obtain the article (DOI 10.1109/TIE.2020.3038099) and follow `sim_r22/README.md`.
- **Superseded raw output:** the Simscape output of a superseded R22 run that used a wrong capacitor value (35 MB). Its summary is kept in `sim_r22/results/archive_cac135/`.
- **Run logs and build caches:** these contained machine-specific paths and nothing needed to reproduce the results.

## Requirements

- Python 3 with the packages in `requirements.txt`. Tested with Python 3.14.5, NumPy 2.4.6, SciPy 1.17.1 (whose `linprog` provides the HiGHS solver), Matplotlib 3.10.9 and Pillow 12.2.0. Pillow is needed only for the optional digitisation step.
- MATLAB R2025b with Simulink and Simscape Electrical, for `sim_simscape/` and `sim_r22/simscape/` only.
- `pdftoppm` (Poppler), for the optional digitisation step only.

```bash
python -m pip install -r requirements.txt
```

All commands below are run from the repository root (`py` is the Windows launcher; use `python` elsewhere). The campaigns use up to 8 worker processes; the times are for a 20-thread laptop.

## Reproducing the paper

| Paper item | Commands | Output |
|---|---|---|
| Fig. 1 (circuit schematic) and Fig. 6 (digital-twin pipeline) | drawn in the manuscript (TikZ) | — |
| Table 1 (closest prior work) | prior-art check, documented in `docs/prior_art/` | — |
| Table 2 (parameter sets) | `sim/ratings.py`, `sim/sizing.py` (values in code) | — |
| Algorithms 1–2 | `sim/inner_circuit.py` (`certify`), `sim/run_wp5a.py` (`_suff_point`), `sim/sizing.py` | — |
| Fig. 2 (equal-energy pair, 120 V) and the pair numbers of Section 4.1 | `py sim/run_pair_experiment.py --rating r120` and `--rating r230` (about 2 min each); `py sim/run_sensitivity.py --rating r120` (and `r230`); `py sim/wp1_report.py` | `sim/results/r120/`, `sim/results/r230/`, `sim/results/wp1_ratings.md` |
| Fig. 3 (phase boundary with the circuit-level bounds) | `py -u sim/run_paper_sweep.py`, then `py -u sim/run_paper_sweep.py scan` (about 1 h each) | `sim/results/paper/phase_sweep_circuit.json` |
| Lab-scale baseline pair and switched check (Section 4.4.1) | `py sim/run_pair_experiment.py`; `py sim/run_switched_check.py` | `sim/results/summary.md`, `sim/results/switched_check.md` |
| Fig. 4 (Simscape against the Python switched model) | `py sim_simscape/export_scenarios.py`; `matlab -batch "cd sim_simscape; run_check"` (57 cases, about 20 min); `py sim_simscape/compare_simscape.py` | `sim_simscape/results/` |
| Fig. 5 (published-experiment benchmark) | `py sim_r22/run_python.py`; `matlab -batch "cd sim_r22/simscape; run_r22"` (about 10 min); `py sim_r22/compare_r22.py` | `sim_r22/results/wp3_r22_check.md` |
| Sections 5.1–5.4 and Table 3 (robustness) | `py -u sim/run_wp4.py --rating r120 --part bounds events mc nonideal` (and `r230`), then the post hoc and circuit-bound parts listed in `sim/README.md`; `py sim/wp4_report.py` | `sim/results/wp4_robustness.md` |
| Section 5.5, Table 4 and Fig. 7 (digital-twin validation) | after the WP4 parts `mc` and `mccircuit`: `py -u sim/run_wp9.py --rating r120` and `--rating r230` (about 40 min each on 12 workers); `py sim/wp9_report.py` | `sim/results/<rating>/wp9_twin.json`, `sim/results/wp9/wp9_twin.md` |
| Sections 6.1–6.6, Figs. 8–9 and Tables 5–6 (sizing) | `py -u sim/run_wp5a.py --part grid extras verify` (about 55 min), then `--part suff_circuit verify2 alpha1 verify_chart verify_conv` (about 2 h); `py sim/wp5a_report.py` | `sim/results/wp5a/wp5a_sizing.md` |
| Section 6.7, Fig. 10 and Tables 7–8 (case study) | `py -u sim/run_wp6.py --part d4 events` (about 2 h); `py sim/wp6_report.py` | `sim/results/wp6/wp6_case_study.md` |
| The 255° consistency check (Section 4.4.4) | `py sim/run_paper_sweep.py check255` (under a minute) | `sim/results/paper/check_255.json` |
| Figs. 2–5 and 7–10 as printed | `py paper/make_figures.py` (after the campaigns above; a few minutes, as it recomputes the time-domain traces of Figs. 2 and 4) | `paper/figs/*.pdf` |

Every report is regenerated from its JSON files by the `*_report.py` script, so the figures and tables can be rebuilt without rerunning the campaigns. The bounds are solved with HiGHS through SciPy, so results can differ in the last digits between solver versions. The certificates and verdicts should not change.

## How to read the results

- **Predeclared criteria** are fixed in `WP4_CRITERIA.md`, `WP5A_PROTOCOL.md`, `WP6_PROTOCOL.md`, `WP9_TWIN_PROTOCOL.md`, `sim_r22/BENCHMARK.md` and `docs/prior_art/WP0_PROTOCOL.md`, each dated before its campaign.
- **Changes made after a method test, before the campaign,** are recorded in each protocol as a revision.
- **Analyses defined after the verdicts were known** are labelled *post hoc* in the protocols and reports.
- **Superseded result files** are kept with a version suffix (for example `wp5a_suff_circuit_v5_flat_tracker.json`), and the protocols say why each was superseded.

## Licence

- **Code** (`.py`, `.m`, `.slx`): MIT licence, see `LICENSE`.
- **Data, figures and documents:** Creative Commons Attribution 4.0, see `DATA_LICENSE.md`.

## Citation

See `CITATION.cff`. Please cite the article \[reference to be added\] and this repository.
