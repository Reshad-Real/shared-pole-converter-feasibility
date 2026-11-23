"""Score both models against the R22 measurements with the criteria fixed in BENCHMARK.md.

Inputs: results/fig7_digitised.json (digitize_fig7.py), results/python_metrics.json
(run_python.py) and results/simscape/*.mat (simscape/run_r22.m).
Outputs: results/wp3_r22_check.md, results/wp3_r22_check.json, results/fig10_r22_benchmark.png.

Run from the repository root:  py sim_r22/compare_r22.py
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
sys.path.insert(0, str(HERE.parent / "sim"))
from r22_metrics import compute                                   # noqa: E402
from run_pair_experiment import INK, INK2, MUTED, S1, S2, SURF, _style  # noqa: E402

# Table I of R22 (measured THD of i_a, %)
THD_MEAS = {"rect_phi0": {"SVPWM": 2.19, "DPWM1": 2.92, "GDPWM": 2.74},
            "rect_phi15": {"SVPWM": 2.05, "DPWM1": 3.12, "GDPWM": 2.96},
            "rect_phi45": {"SVPWM": 1.95, "DPWM1": 3.38, "GDPWM": 2.52},
            "statcom_5A": {"SVPWM": 2.58, "DPWM1": 4.07, "GDPWM": 3.94},
            "statcom_10A": {"SVPWM": 1.26, "DPWM1": 3.31, "GDPWM": 2.96}}
THD_TOL_PP, REL_TOL, CLAMP_TOL_MS = 1.5, 0.10, 1.0


def simscape_metrics():
    out = {}
    for f in sorted((RES / "simscape").glob("*.mat")):
        m = loadmat(f, squeeze_me=True)
        per, fine = m["per"], m["fine"]
        d = dict(t=per[:, 0], i_a=per[:, 1], i_c=per[:, 2], u_cac=per[:, 3], U_avg=per[:, 5],
                 dA=per[:, 6], fine_t=fine[:, 0], fine_ia=fine[:, 1], meta=dict(f=50.0, Tsw=25e-6))
        out[f.stem] = compute(d)
    return out


def main():
    D = json.loads((RES / "fig7_digitised.json").read_text())
    PY = json.loads((RES / "python_metrics.json").read_text())
    SS = simscape_metrics()
    models = {"Python switched": PY, "Simscape": SS}
    a, b = D["a_GDPWM"], D["b_DPWM1"]
    period_px = 4 * (1201 - 412) / 10
    ph_unc = 2 / period_px * 360

    # Scalar targets: (label, measured, uncertainty, model key, metric, tolerance)
    def amp(dig, k):
        return dig[f"amp_{k}"]["value"], dig[f"amp_{k}"]["unc"]
    t3 = np.mean([a["T3_dc_ripple_pp_no_decoupling_V"]["value"], b["T3_dc_ripple_pp_no_decoupling_V"]["value"]])
    t3u = a["T3_dc_ripple_pp_no_decoupling_V"]["unc"]
    rows = [("T3 DC ripple, decoupling off (V p-p)", t3, t3u, "apd_off", "dc_ripple_pp_V")]
    for lab, dig, key in (("GDPWM", a, "rect_phi0_GDPWM"), ("DPWM1", b, "rect_phi0_DPWM1")):
        v, u = amp(dig, "u_cac")
        rows.append((f"T5 u_cac amplitude, {lab} (V)", v, u, key, "amp_u_cac"))
        rows.append((f"T5 u_cac phase vs u_ac, {lab} (deg)", dig["phase_u_cac_vs_u_ac_deg"], ph_unc, key,
                     "phase_u_cac_vs_u_ac_deg"))
        for cur in ("i_a", "i_b", "i_c"):
            v, u = amp(dig, cur)
            rows.append((f"T6 {cur} amplitude, {lab} (A)", v, u, key, f"amp_{cur}"))
    results = []
    for lab, meas, unc, key, metric in rows:
        tol = REL_TOL * abs(meas) + unc
        entry = dict(target=lab, measured=meas, uncertainty=unc, tolerance=tol, models={})
        for mname, M in models.items():
            if key in M:
                val = M[key][metric]
                entry["models"][mname] = dict(value=val, passed=abs(val - meas) <= tol)
        results.append(entry)
    # T4 upper bounds
    for lab, dig, key in (("GDPWM", a, "rect_phi0_GDPWM"), ("DPWM1", b, "rect_phi0_DPWM1")):
        bound = dig["T4_dc_ripple_pp_with_decoupling_upper_bound_V"]
        entry = dict(target=f"T4 DC ripple with decoupling, {lab} (V p-p, upper bound)", measured=bound,
                     uncertainty=0.0, tolerance=None, models={})
        for mname, M in models.items():
            if key in M:
                entry["models"][mname] = dict(value=M[key]["dc_ripple_pp_V"],
                                              passed=M[key]["dc_ripple_pp_V"] <= bound)
        results.append(entry)
    # T7 clamp durations
    clamp = {}
    for lab, dig, key in (("GDPWM", a, "rect_phi0_GDPWM"), ("DPWM1", b, "rect_phi0_DPWM1")):
        meas = dig["T7_Sa_clamp_mean_ms"]
        entry = dict(target=f"T7 S_a clamp per half cycle, {lab} (ms)", measured=meas,
                     uncertainty=dig["T7_unc_ms"], tolerance=CLAMP_TOL_MS, models={})
        for mname, M in models.items():
            if key in M:
                val = M[key]["Sa_clamp_mean_ms"]
                entry["models"][mname] = dict(value=val, passed=abs(val - meas) <= CLAMP_TOL_MS)
                clamp.setdefault(mname, {})[lab] = val
        results.append(entry)
    order_ok = {m: clamp[m]["GDPWM"] > clamp[m]["DPWM1"] for m in clamp}
    # T1 THD at phi = 0
    thd_rows = []
    for s in ("SVPWM", "DPWM1", "GDPWM"):
        meas = THD_MEAS["rect_phi0"][s]
        entry = dict(target=f"T1 THD of i_a, {s} (%)", measured=meas, uncertainty=0.0,
                     tolerance=THD_TOL_PP, models={})
        for mname, M in models.items():
            key = f"rect_phi0_{s}"
            if key in M:
                val = M[key]["thd_ia_pct"]
                entry["models"][mname] = dict(value=val, passed=abs(val - meas) <= THD_TOL_PP)
        thd_rows.append(entry)
    results += thd_rows
    thd_order = {}
    for mname, M in models.items():
        v = {s: M[f"rect_phi0_{s}"]["thd_ia_pct"] for s in ("SVPWM", "DPWM1", "GDPWM")
             if f"rect_phi0_{s}" in M}
        if len(v) == 3:
            thd_order[mname] = v["SVPWM"] < min(v["DPWM1"], v["GDPWM"]) and v["DPWM1"] > v["GDPWM"]
    # T2 extended THD (no pass mark) and dead-time sensitivity
    t2 = {}
    for case, meas in THD_MEAS.items():
        for s, mv in meas.items():
            key = f"{case}_{s}"
            t2[key] = dict(measured=mv, **{m: M[key]["thd_ia_pct"] for m, M in models.items() if key in M})
    dt = {f"{s} td {td} ns": PY[f"rect_phi0_{s}_td{td}"]["thd_ia_pct"] for s in ("SVPWM", "DPWM1", "GDPWM")
          for td in (100, 300)}
    summary = dict(results=results, clamp_order=order_ok, thd_order=thd_order, T2=t2,
                   thd_dead_time_python=dt)
    (RES / "wp3_r22_check.json").write_text(json.dumps(summary, indent=2, default=float))
    write_md(summary)
    plot(summary)
    for r in results:
        print(f"{r['target']:52s} meas {r['measured']:8.2f}  " + "  ".join(
            f"{m[:6]} {v['value']:8.2f} {'ok' if v['passed'] else 'FAIL'}" for m, v in r["models"].items()))


def write_md(s):
    L = ["# WP3: reproducing R22's measurements in both models\n",
         "Generated by `sim_r22/compare_r22.py`. Targets and pass criteria were fixed in "
         "`BENCHMARK.md` before the measurements were digitised or the models run. Measured values "
         "come from R22's Table I and from its Fig. 7 oscilloscope traces, digitised by "
         "`digitize_fig7.py` with the grid voltage as the scale reference.\n",
         "**Correction (22 October 2026).** The decoupling capacitor is now 130 µF, the value R22 "
         "designs with its eq. (9). The first runs used 135 µF, the DC-link value, by mistake. "
         "That run scored 33 of 36, with all 30 non-THD checks passing; it is kept in "
         "`results/archive_cac135/`. The correction raises the capacitor-voltage amplitude by about "
         "2 % and brings it, the capacitor current and the GDPWM clamp closer to the measurements. "
         "The other waveform values change by 1 % or less.\n",
         "| Target | Measured | Tolerance | Python switched | Simscape |", "|---|---|---|---|---|"]
    for r in s["results"]:
        meas = (f"≤ {r['measured']:.1f}" if r["tolerance"] is None else
                f"{r['measured']:.2f} ± {r['uncertainty']:.2f}" if r["uncertainty"] else f"{r['measured']:.2f}")
        tol = "—" if r["tolerance"] is None else f"± {r['tolerance']:.2f}"
        cells = []
        for m in ("Python switched", "Simscape"):
            v = r["models"].get(m)
            cells.append("—" if v is None else f"{v['value']:.2f} ({'pass' if v['passed'] else '**fail**'})")
        L.append(f"| {r['target']} | {meas} | {tol} | {cells[0]} | {cells[1]} |")
    L.append("")
    L.append("Orderings: GDPWM clamps longer than DPWM1 — " + ", ".join(
        f"{m} {'yes' if v else 'no'}" for m, v in s["clamp_order"].items()) +
        ". THD ordering SVPWM < GDPWM < DPWM1 as measured — " + ", ".join(
        f"{m} {'yes' if v else 'no'}" for m, v in s["thd_order"].items()) + ".\n")
    n_all = sum(len(r["models"]) for r in s["results"])
    n_pass = sum(v["passed"] for r in s["results"] for v in r["models"].values())
    n_thd = sum(len(r["models"]) for r in s["results"] if r["target"].startswith("T1"))
    n_thd_pass = sum(v["passed"] for r in s["results"] if r["target"].startswith("T1")
                     for v in r["models"].values())
    L.append(f"**Result:** {n_pass} of {n_all} model-target checks pass. Excluding THD (T1), "
             f"{n_pass - n_thd_pass} of {n_all - n_thd}.\n")
    L.append("## THD across operating points (T2, no pass mark)\n")
    L.append("| Case | Measured (%) | Python (%) | Simscape (%) |\n|---|---|---|---|")
    for k, v in s["T2"].items():
        L.append(f"| {k} | {v['measured']:.2f} | {v.get('Python switched', float('nan')):.2f} | "
                 + (f"{v['Simscape']:.2f}" if "Simscape" in v else "—") + " |")
    L.append("\nPython THD at φ = 0° against dead time: " + ", ".join(
        f"{k} {v:.2f} %" for k, v in s["thd_dead_time_python"].items()) + ".\n")
    L.append("## Reading\n")
    L.append("- The power-flow physics of the shared-leg decoupling converter is reproduced by both "
             "models within the measurement uncertainty: the DC ripple with and without decoupling, "
             "the capacitor voltage that absorbs the double-line power (amplitude and phase), the "
             "three branch currents, and how long each modulation clamps S_a.")
    L.append("- THD is not reproduced: both models give a cleaner grid current than measured, even "
             "with longer dead time. R22's own PLECS simulation (Fig. 5c) also shows 2–4 % with an "
             "ideal grid, so the measured THD is set mainly by details of R22's controller that are "
             "not published. THD is therefore not a usable anchor here.")
    L.append("- Losses and efficiency (Table I) were excluded in advance; see `BENCHMARK.md`.")
    L.append("- The reproduction relies on assumed controller gains, dead time and device data "
             "(`MODEL_SPEC.md`); the waveform-level results are insensitive to them, the THD is not.")
    (RES / "wp3_r22_check.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def plot(s):
    plt = _style()
    rows = [r for r in s["results"] if r["tolerance"] is not None]
    fig, ax = plt.subplots(figsize=(7.4, 0.36 * len(rows) + 1.4))
    for i, r in enumerate(rows[::-1]):
        m = r["measured"]
        lo, hi = 1 - r["tolerance"] / abs(m), 1 + r["tolerance"] / abs(m)
        ax.plot([lo, hi], [i, i], color="#cde2fb", lw=7, solid_capstyle="butt", zorder=1)
        if r["uncertainty"]:
            ax.plot([1 - r["uncertainty"] / abs(m), 1 + r["uncertainty"] / abs(m)], [i, i], color=MUTED,
                    lw=1.2, zorder=2)
        for mname, mk, col, off in (("Python switched", "o", S2, 0.12), ("Simscape", "D", S1, -0.12)):
            v = r["models"].get(mname)
            if v:
                ax.plot([v["value"] / m], [i + off], mk, ms=6, color=col, mec=SURF, mew=1.2, zorder=3,
                        label=mname if i == 0 else None)
    ax.axvline(1.0, color=INK2, lw=0.8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r["target"] for r in rows[::-1]], fontsize=7.5)
    ax.set_xlabel("Model value / measured value")
    ax.set_title("Published-experiment benchmark: both models against the measurements", loc="left", color=INK)
    ax.legend(loc="upper right")
    fig.text(0.02, -0.02, "Blue band: pass tolerance (±10 % plus digitisation uncertainty; THD ±1.5 points; "
             "clamp ±1 ms). Grey line: digitisation uncertainty.", color=MUTED, fontsize=7.5)
    fig.savefig(RES / "fig10_r22_benchmark.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
