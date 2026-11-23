"""WP5a report: sizing rule, design charts and the S1-S6 checks of WP5A_PROTOCOL.md.

Reads results/wp5a/wp5a_{grid,extras,verify}.json (written by run_wp5a.py) and writes
results/wp5a/wp5a_sizing.md, results/wp5a/wp5a_summary.json and figures 14-15.

Run:  py sim/wp5a_report.py
"""
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from run_pair_experiment import INK, INK2, MUTED, SURF, _style

RES = Path(__file__).resolve().parent / "results" / "wp5a"
RATINGS = ("r120", "r230")
NAMES = {"r120": "120 V rating", "r230": "230 V rating"}
M_GRID = (0.30, 0.35, 0.40, 0.45, 0.50)
BLUES = ("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281")   # ordinal ramp, validated
ORANGES = ("#f7a072", "#eb6834", "#b4461b")                       # alpha 0.5, 0.75, 1.0
S1_PCT, S3_PTS = 15.0, 10.0
TRACKER = "smooth"             # the switched-model tracker S4(a) is judged with
TRACKER_NAMES = {"smooth": "smoothed-bus tracker", "flat": "flat-bus tracker",
                 "energy": "energy tracker (follows the plan's bus)"}


def _load(name):
    f = RES / f"wp5a_{name}.json"
    return json.loads(f.read_text()) if f.exists() else None


TRACKERS_ORDER = ("smooth", "flat", "energy")


def _track_res(v):
    """One switched-model result as text."""
    if v is None:
        return "—"
    if "error" in v:
        return v["error"]
    res = "held" if not v["violation"] else (
        f"trip at {v['t_violation_ms']:.1f} ms ("
        + ", ".join(w for w, b in (("bus", v["bus_violation"]), ("pole current", v["current_violation"])) if b) + ")")
    return res + f"; bus {v['Vdc_min']:.0f}–{v['Vdc_max']:.0f} V, current {v['I_pole_peak']:.1f} A"


def _margin(x):
    """Smallest distance of the tracked bus from its limits (negative: violated)."""
    return min(x["Vdc_min"] - x["Vdc_lim"][0], x["Vdc_lim"][1] - x["Vdc_max"])


def _ramp_txt(x):
    return "bus-loop" if x["ramp_ms"] is None else f"{x['ramp_ms']:g} ms"


def _pol(x):
    """Tracker of a switched-model result (older results name the policy instead)."""
    if "tracker" in x:
        return x["tracker"]
    return {"plan_flat": "flat", "plan": "energy"}.get(x.get("policy", "plan"), x.get("policy"))


def _tau(pt):
    return pt["T_r_s"] * pt["f"]


def _sel(pts, **kw):
    return [p for p in pts if all(p[k] == v for k, v in kw.items())]


def chart_table(pts, rating):
    """Penalty over m x tau at alpha = 0.5 (rated step, unity power factor)."""
    rows = _sel(pts, rating=rating, group="chart")
    taus = sorted({round(_tau(p), 4) for p in rows})
    tab = np.full((len(M_GRID), len(taus)), np.nan)
    for p in rows:
        tab[M_GRID.index(p["m"]), taus.index(round(_tau(p), 4))] = p["penalty_pct"]
    return np.array(taus), tab


def transfer(pts):
    """S3: bilinear interpolation of the 120 V chart predicts the 230 V penalty. Points where
    either rating has no design (no steady orbit at any Cb) are listed but not scored."""
    t1, tab1 = chart_table(pts, "r120")
    rows_ok = [i for i in range(len(M_GRID)) if np.all(np.isfinite(tab1[i]))]
    ms = np.array([M_GRID[i] for i in rows_ok])
    f = RegularGridInterpolator((ms, t1), tab1[rows_ok], method="linear")
    t2, tab2 = chart_table(pts, "r230")
    out = []
    for i, m in enumerate(M_GRID):
        if m > ms[-1]:
            continue
        for j, t in enumerate(t2):
            if t1[0] <= t <= t1[-1]:
                pred = float(f([[m, t]])[0])
                meas = float(tab2[i, j])
                ok = np.isfinite(pred) and np.isfinite(meas)
                out.append(dict(m=m, tau=float(t), measured=meas, predicted=pred,
                                error_pts=meas - pred if ok else None))
    scored = [o for o in out if o["error_pts"] is not None]
    worst = max(abs(o["error_pts"]) for o in scored)
    sub = [o for o in scored if o["m"] <= 0.40]
    worst_sub = max(abs(o["error_pts"]) for o in sub)
    return dict(points=out, n=len(scored), n_no_design=len(out) - len(scored),
                max_abs_error_pts=worst, passed=worst <= S3_PTS,
                max_abs_error_pts_m_le_040=worst_sub, n_m_le_040=len(sub))


def fit_closed_form(pts):
    """Post hoc: penalty = a + b ln(tau) + c (m - 0.4), fitted on 120 V, tested on 230 V,
    over the points with tau >= 0.1 line cycles (below that the effect is small)."""
    def rows(r):
        return [p for p in _sel(pts, rating=r, group="chart")
                if _tau(p) >= 0.1 and p["m"] <= 0.40 and np.isfinite(p["penalty_pct"])]
    A = np.array([[1.0, np.log(_tau(p)), p["m"] - 0.4] for p in rows("r120")])
    y = np.array([p["penalty_pct"] for p in rows("r120")])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)

    def err(r):
        X = np.array([[1.0, np.log(_tau(p)), p["m"] - 0.4] for p in rows(r)])
        z = np.array([p["penalty_pct"] for p in rows(r)])
        return float(np.max(np.abs(X @ coef - z))), float(np.sqrt(np.mean((X @ coef - z) ** 2)))
    return dict(form="penalty % = a + b ln(T_r f) + c (m - 0.4), for T_r f >= 0.1, m <= 0.40",
                a=float(coef[0]), b=float(coef[1]), c=float(coef[2]),
                fit_max_abs_pts=err("r120")[0], fit_rms_pts=err("r120")[1],
                test_max_abs_pts=err("r230")[0], test_rms_pts=err("r230")[1])


# ---------------------------------------------------------------- figures
def fig_chart(pts, path):
    plt = _style()
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.9), sharey=True, gridspec_kw=dict(wspace=0.08))
    for ax, r in zip(axs, RATINGS):
        taus, tab = chart_table(pts, r)
        for i, m in enumerate(M_GRID[:-1]):                 # m = 0.50 has no design
            y = np.where(tab[i] > 175, np.nan, tab[i])     # off-scale values are not drawn
            ax.plot(taus, y, color=BLUES[i], marker="o", ms=4, label=f"m = {m:.2f}")
        ax.set_ylim(-60, 180)
        off = tab[3][np.isfinite(tab[3]) & (tab[3] > 175)]
        if off.size:
            ax.text(taus[0], 165, f"m = 0.45 above {taus[0]:.2f} cycles: {off.min():.0f}–"
                    f"{off.max():.0f} %, off scale", color=MUTED, fontsize=7.5)
        ax.axhline(S1_PCT, color=MUTED, lw=1.0, ls=(0, (4, 3)))
        ax.axhline(0, color=INK2, lw=0.8)
        ax.axvline(1 / (np.pi * 0.4), color=MUTED, lw=0.8, ls=":")
        ax.text(1 / (np.pi * 0.4) * 1.05, 2, "bus loop", color=MUTED, fontsize=7.5, va="bottom")
        ax.set_xscale("log")
        ax.set_title(NAMES[r], loc="left")
        ax.set_xlabel("DC1 ramp time for the rated step (line cycles)")
    axs[0].set_ylabel("Capacitance phase-blind sizing misses (%)")
    axs[0].text(axs[0].get_xlim()[0] * 1.1, S1_PCT + 2, "15 % mark", color=MUTED, fontsize=7.5)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower left", bbox_to_anchor=(0.08, 0.98), ncol=5, title=None)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_alpha(pts, path):
    plt = _style()
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.9), sharey=True, gridspec_kw=dict(wspace=0.08))
    for ax, r in zip(axs, RATINGS):
        for k, a in enumerate((0.5, 0.75, 1.0)):
            grp = "chart" if a == 0.5 else "alpha"
            rows = sorted([p for p in _sel(pts, rating=r, group=grp, alpha=a) if p["m"] == 0.40],
                          key=_tau)
            x = [_tau(p) for p in rows]
            ref = rows[0]["cb_steady"]
            ax.plot(x, [p["cb_nec"] / ref for p in rows], color=ORANGES[k], marker="o", ms=4,
                    label=f"Necessary, α = {a:g}")
            ax.plot(x, [p["cb_conv"] / ref for p in rows], color=ORANGES[k], lw=1.0,
                    ls=(0, (4, 2)), label=f"Phase-blind, α = {a:g}")
        ax.axhline(1.0, color=INK2, lw=0.8)
        ax.set_xscale("log")
        ax.set_title(NAMES[r], loc="left")
        ax.set_xlabel("DC1 ramp time for the rated step (line cycles)")
    axs[0].set_ylabel("Buffer capacitance /\nsteady-state minimum at rated power")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower left", bbox_to_anchor=(0.08, 0.98), ncol=3)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- report
def _uf(x):
    return "> 20 mF" if not np.isfinite(x) else f"{x * 1e6:.0f} µF"


def main():
    pts, ext, ver = _load("grid"), _load("extras"), _load("verify")
    base = {p["rating"]: p for p in pts if p["is_base"]}
    S = {}
    for r in RATINGS:
        S[r] = dict(S1=base[r]["penalty_pct"] >= S1_PCT, penalty_pct=base[r]["penalty_pct"])
    tr = transfer(pts)
    fit = fit_closed_form(pts)
    if ver:
        for r in RATINGS:
            v = {x["case"]: x for x in ver if x.get("rating") == r}
            S[r]["S4"] = (not v["plan_at_suff"].get("violation", True)
                          and v["conventional_at_conv"]["violation"]
                          and v["conventional_at_0.95nec"]["violation"])
    v2 = _load("verify2")
    if v2 is not None:
        for r in RATINGS:
            v = [x for x in v2 if x["rating"] == r and "error" not in x and _pol(x) == TRACKER]
            S[r]["S4a_current_aware"] = (not v[0]["violation"]) if v else None
    a1 = _load("alpha1")
    if a1:
        for r in RATINGS:
            v = [x for x in a1["verify"] if x["rating"] == r and "error" not in x
                 and _pol(x) == TRACKER]
            S[r]["S4a_current_aware_alpha1"] = (not v[0]["violation"]) if v else None
    vcv = _load("verify_conv")
    if vcv:
        for r in RATINGS:
            S[r]["S4bc_with_timing_fix"] = all(x["violation"] for x in vcv if x["rating"] == r)
    vc = _load("verify_chart")
    if vc:
        for r in RATINGS:
            for pol in sorted({_pol(x) for x in vc}):
                v = [x for x in vc if x["rating"] == r and "error" not in x and _pol(x) == pol]
                S[r][f"S4a_chart_{pol}"] = dict(held=sum(not x["violation"] for x in v),
                                                 tested=len(v))
    (RES / "wp5a_summary.json").write_text(json.dumps(dict(summary=S, transfer=tr, fit=fit),
                                                      indent=1, default=float))
    fig_chart(pts, RES / "fig14_wp5a_design_chart.png")
    fig_alpha(pts, RES / "fig15_wp5a_operating_point.png")
    write_md(pts, base, ext, ver, S, tr, fit)
    print(json.dumps(S, indent=1, default=float))


def write_md(pts, base, ext, ver, S, tr, fit):
    L = ["# WP5a: buffer sizing against phase-blind sizing\n",
         "Generated by `sim/wp5a_report.py` from `results/wp5a/wp5a_*.json` (written by "
         "`sim/run_wp5a.py`). The protocol is `sim/WP5A_PROTOCOL.md`, including its revision "
         "before the campaign. Everything is simulation of the illustrative circuit.\n",
         "## Verdicts\n",
         "| Criterion | 120 V rating | 230 V rating |\n|---|---|---|"]

    def pf(b):
        return "**pass**" if b else "**fail**"
    L.append("| S1: penalty at the base design (α = 0.5) ≥ 15 % | " + " | ".join(
        f"{pf(S[r]['S1'])} ({S[r]['penalty_pct']:.0f} %)" for r in RATINGS) + " |")
    L.append(f"| S3: the 120 V chart predicts the 230 V penalty within ±{S3_PTS:.0f} points | "
             f"{pf(tr['passed'])} (largest error {tr['max_abs_error_pts']:.0f} points, at m = 0.45; "
             f"{tr['max_abs_error_pts_m_le_040']:.0f} points at m ≤ 0.40; m = 0.50 has no design) | |")
    if ver:
        L.append("| S4: switched-model verification | " + " | ".join(pf(S[r]["S4"]) for r in RATINGS) + " |")
    if "S4a_current_aware" in S[RATINGS[0]]:
        L.append("| Post hoc: S4(a) with the current-aware sufficient bound (centred; m = 0.35, bus-loop ramp; no base design certified) | " + " | ".join(
            "no certified design (pole-current limited)" if S[r]["S4a_current_aware"] is None
            else pf(S[r]["S4a_current_aware"]) for r in RATINGS) + " |")
    if "S4a_current_aware_alpha1" in S[RATINGS[0]]:
        L.append("| Post hoc: the same at the top operating point (α = 1) | " + " | ".join(
            "no certified design" if S[r]["S4a_current_aware_alpha1"] is None
            else pf(S[r]["S4a_current_aware_alpha1"]) for r in RATINGS) + " |")
    nv = sum(p.get("n_no_verdict", 0) for p in pts)
    L.append(f"\nPhase searches without a solver verdict, even after retrying two other HiGHS "
             f"methods: {nv} of {sum(len(p['energy']) + len(p['outer']) + len(p['inner']) + len(p['circuit']) + len(p['inner_extra']) for p in pts)}"
             f"{' (left out of the maxima)' if nv else ''}.\n")

    L.append("## Base design (m = 0.40, bus-loop DC1 ramp, rated step, unity power factor, α = 0.5)\n")
    L.append("| | 120 V rating | 230 V rating |\n|---|---|---|")
    rows = (("Bus voltage", lambda p: f"{p['Vdc']:.1f} V"),
            ("AC2 step", lambda p: f"{p['dP_W']:.0f} W"),
            ("DC1 ramp time for the rated step", lambda p: f"{p['T_r_s'] * 1e3:.1f} ms ({_tau(p):.2f} cycles)"),
            ("Steady-state minimum at the post-step power", lambda p: _uf(p["cb_steady"])),
            ("Energy LP, worst phase", lambda p: _uf(max(p["energy"]))),
            ("**Phase-blind sizing, Cb_conv**", lambda p: f"**{_uf(p['cb_conv'])}**"),
            ("**Necessary, Cb_nec (circuit outer bound)**",
             lambda p: f"**{_uf(p['cb_nec'])}** at {p['phase_nec_deg']:.1f}°"),
            ("Reduced outer bound (no inductor), for reference", lambda p: _uf(p["cb_outer_reduced"])),
            ("Sufficient, Cb_suff (inner bound)", lambda p: f"{_uf(p['cb_suff'])} at {p['phase_suff_deg']:.1f}°"),
            ("**Penalty, Cb_nec / Cb_conv − 1**", lambda p: f"**{p['penalty_pct']:.0f} %**"),
            ("Bracket, Cb_suff / Cb_nec", lambda p: f"{p['cb_suff'] / p['cb_nec']:.2f}"))
    for lab, fn in rows:
        L.append(f"| {lab} | " + " | ".join(fn(base[r]) for r in RATINGS) + " |")
    L.append("")

    L.append("## Design chart (S2)\n")
    L.append("Penalty in % at α = 0.5, rated step, unity power factor. Columns: DC1 ramp time "
             "for the rated step, in line cycles.\n")
    for r in RATINGS:
        taus, tab = chart_table(pts, r)
        L.append(f"**{NAMES[r]}**\n")
        L.append("| m | " + " | ".join(f"{t:.2f}" for t in taus) + " |")
        L.append("|---|" + "---|" * len(taus))
        for i, m in enumerate(M_GRID):
            L.append(f"| {m:.2f} | " + " | ".join(
                "no design" if not np.isfinite(v) else
                ("**" if v >= S1_PCT else "") + f"{v:.0f}".replace("-", "−") + ("**" if v >= S1_PCT else "")
                for v in tab[i]) + " |")
        fin = np.isfinite(tab)
        below = int(np.sum(tab[fin] < S1_PCT))
        L.append(f"\nBold: at least 15 %. Below 15 % at {below} of {int(fin.sum())} points with a design.\n")
    L.append("- **m = 0.50 has no design:** the pole band closes, so no steady-state orbit fits at any "
             "capacitance up to 20 mF.")
    L.append("- **m = 0.45 is at the feasibility edge:** no sufficient design exists at either rating, and "
             "at 230 V even the reduced outer bound finds none. The penalty there grows steeply and "
             "differs between ratings.")
    L.append("- **Negative values at the fastest ramp (0.5 ms):** the transient bound falls below the "
             "steady-state minimum that conventional sizing already includes. With a fast DC1, phase-aware "
             "sizing adds nothing.")
    L.append("- **Penalty grows with modulation index and with ramp time.** It already passes 15 % at "
             "0.1 line cycle for m ≥ 0.30 at both ratings.\n")
    L.append("![Design chart](fig14_wp5a_design_chart.png)\n")

    L.append("## Chart transfer between ratings (S3)\n")
    big = sorted((o for o in tr["points"] if o["error_pts"] is not None),
                 key=lambda o: -abs(o["error_pts"]))[:3]
    L.append(f"The 120 V chart, interpolated bilinearly in (m, ramp time in line cycles), predicts "
             f"the 230 V penalty with a largest error of {tr['max_abs_error_pts']:.0f} points over "
             f"{tr['n']} scored points, so S3 fails. Largest errors: " + "; ".join(
                 f"m = {o['m']:.2f}, {o['tau']:.2f} cycles: {o['measured']:.0f} % against "
                 f"{o['predicted']:.0f} % predicted" for o in big) + ".\n")
    L.append(f"- **Restricted to m ≤ 0.40 (post hoc),** the largest error is "
             f"{tr['max_abs_error_pts_m_le_040']:.1f} points over {tr['n_m_le_040']} points. It sits "
             f"where the penalty rises steeply at short ramps (0.03 to 0.30 line cycles), which "
             f"the coarse ramp grid there cannot follow. Elsewhere the ratings agree closely "
             f"(base design 73 % against 74 %).")
    L.append("- **At m = 0.45 the ratings diverge strongly.** The effective modulation index, peak "
             "converter voltage over bus voltage, equals the nominal one to within 0.004 at both "
             "ratings, so the filter inductors do not explain this.\n")

    L.append("## Operating point\n")
    L.append("Before the step the buffer runs on a steady orbit whose mean energy sits at position "
             "α of its admissible interval (0.5 centred, 1 highest). m = 0.40, rated step.\n")
    L.append("| Rating | DC1 ramp (cycles) | α | Cb_conv | Cb_nec | Penalty | Cb_nec / Cb_nec(α = 0.5) |\n"
             "|---|---|---|---|---|---|---|")
    for r in RATINGS:
        ref = {round(_tau(p), 3): p["cb_nec"] for p in _sel(pts, rating=r, group="chart") if p["m"] == 0.40}
        for a in (0.5, 0.75, 1.0):
            grp = "chart" if a == 0.5 else "alpha"
            for p in sorted([q for q in _sel(pts, rating=r, group=grp, alpha=a) if q["m"] == 0.40], key=_tau):
                L.append(f"| {NAMES[r]} | {_tau(p):.2f} | {a:g} | {_uf(p['cb_conv'])} | {_uf(p['cb_nec'])} | "
                         f"{p['penalty_pct']:.0f} % | {p['cb_nec'] / ref[round(_tau(p), 3)]:.2f} |".replace("| -", "| −"))
    b0 = {r: base[r]["cb_conv"] for r in RATINGS}
    top = {r: next(p for p in _sel(pts, rating=r, group="alpha", alpha=1.0) if p["ramp_ms"] is None)
           for r in RATINGS}
    L.append("\n**Reading.**\n")
    L.append("- **Raising the operating point to the top of its interval (α = 1) lowers the necessary "
             "capacitance by about a third** at ramps of half a line cycle or more.")
    L.append(f"- **The penalty remains.** At the base ramp it is "
             f"{top['r120']['penalty_pct']:.0f} % and {top['r230']['penalty_pct']:.0f} % against "
             f"phase-blind sizing at the same operating point. Against phase-blind sizing as usually "
             f"computed from the centred orbit ({_uf(b0['r120'])}, {_uf(b0['r230'])}), it is "
             f"{100 * (top['r120']['cb_nec'] / b0['r120'] - 1):.0f} % and "
             f"{100 * (top['r230']['cb_nec'] / b0['r230'] - 1):.0f} %.")
    L.append("- **At ramps of about 0.1 line cycle,** the top operating point removes the penalty "
             "altogether.\n")
    L.append("![Operating point](fig15_wp5a_operating_point.png)\n")

    L.append("## Step size and AC2 power factor (base design, α = 0.5)\n")
    L.append("| Rating | Case | Cb_conv | Cb_nec | Cb_suff | Penalty |\n|---|---|---|---|---|---|")
    for r in RATINGS:
        cases = [(f"step to {int(100 * p['p2_frac'])} % of rating", p) for p in _sel(pts, rating=r, group="step")]
        cases += [("step to 100 % (base)", base[r])]
        cases += [(f"power factor {p['pf_label']}", p) for p in _sel(pts, rating=r, group="pf")]
        for lab, p in cases:
            L.append(f"| {NAMES[r]} | {lab} | {_uf(p['cb_conv'])} | {_uf(p['cb_nec'])} | {_uf(p['cb_suff'])} | "
                     f"{p['penalty_pct']:.0f} % |")
    L.append("")

    if ext:
        L.append("## Load drop against operating point (S6, indicative)\n")
        L.append("A higher operating point helps step-ups and hurts load drops. Smallest Cb that "
                 "admits the drop of AC2 to zero load at every phase, from the same orbit. The circuit "
                 "outer LP is valid but loose for surplus, so its load-drop value is a weak lower bound.\n")
        dr = _load("drop_reduced") or []
        drd = {(d["rating"], d["alpha"]): d["down_reduced_outer"] for d in dr}
        L.append("| Rating | α | Step-up Cb_nec | Load drop, energy LP | Load drop, circuit outer | "
                 "Load drop, reduced outer (post hoc) |\n|---|---|---|---|---|---|")
        for e in ext:
            for t in e["tradeoff"]:
                rd = drd.get((e["rating"], t["alpha"]))
                L.append(f"| {NAMES[e['rating']]} | {t['alpha']:g} | {_uf(t['up_cb_nec'])} | "
                         f"{_uf(t['down_energy'])} | {_uf(t['down_circuit'])} | "
                         f"{_uf(rd) if rd is not None else '—'} |")
        L.append("\n- **The energy LP and the circuit outer LP** both stop at the smallest capacitance for "
                 "which the pre-step orbit exists. They cannot see the cost of a high operating point "
                 "on load drops.")
        L.append("- **The reduced outer LP** keeps the buffer voltage inside the band, so it is tight for "
                 "surplus, though it ignores the inductor. Added post hoc, it shows the load drop needs "
                 "only a little more capacitance as α rises, far less than the step-up needs.")
        L.append("- **So for this event set, a high operating point at light load is nearly free.** A "
                 "load drop from rated load starts from a different operating point (the orbit at rated "
                 "load), which a controller can place low.")
        L.append("\n## Checks (S5)\n")
        for e in ext:
            mono = ", ".join(f"{m['criterion']} {'yes' if m['monotone'] else '**no**'}" for m in e["monotonicity"])
            L.append(f"- {NAMES[e['rating']]}: feasibility monotone in Cb at the worst phase: {mono}. "
                     f"Circuit bound with three versus five post-step cycles: "
                     f"{_uf(e['circuit_3_cycles'])} versus {_uf(e['circuit_5_cycles'])}.")
        L.append("")

    if ver:
        L.append("## Switched-model verification (S4)\n")
        L.append("| Rating | Case | Cb | Phase | Result |\n|---|---|---|---|---|")
        for v in ver:
            if "error" in v:
                L.append(f"| — | {v['case']} | — | — | {v['error']} |")
                continue
            res = "held" if not v["violation"] else f"tripped at {v['t_violation_ms']:.1f} ms"
            L.append(f"| {NAMES[v['rating']]} | {v['case']} | {_uf(v['cb'])} | {v['phase_deg']:.1f}° | {res} |")
        L.append("\n(b) and (c) hold: conventional sizing and anything below the necessary value trip. "
                 "(a) fails: the inner bound leaves out the pole-current limit and the voltage the buffer "
                 "inductor needs, and its plans exceed the current limit by 1–2 A.\n")
        scan = _load("suff_scan")
        if scan:
            L.append("**Post hoc: the same plan with more capacitance.** Inner-LP plan tracked at the "
                     "worst phase with Cb = k × Cb_suff.\n")
            L.append("| Rating | k | Cb | Plan's own peak pole-N current | Switched result |\n|---|---|---|---|---|")
            for s in scan:
                res = "held" if not s["violation"] else (
                    f"trip at {s['t_violation_ms']:.1f} ms ("
                    + ", ".join(w for w, b in (("bus", s["bus_violation"]), ("pole current", s["current_violation"])) if b) + ")")
                L.append(f"| {NAMES[s['rating']]} | {s['k']:g} | {_uf(s['cb'])} | {s['plan_peak_iN']:.1f} A "
                         f"(limit {s['I_pole_max']:.0f} A) | {res} |")
            L.append("\n- **120 V:** a plan holds from 1.25 × Cb_suff.")
            L.append("- **230 V:** no capacitance helps. Near the line peak the buffer must carry the DC1 "
                     "deficit on top of the ripple current, which exceeds the pole-current limit, so the "
                     "rated step at this rating is limited by current, not by energy.")
            L.append("- **Consequence:** a certified sufficient value needs an inner bound with the pole-current "
                     "limit and the inductor voltage. The linear version tried here (current within 90 % of "
                     "the limit, buffer voltage bounded by its secant) is infeasible at every capacitance, "
                     "because at rated power the pole-N current is already within about 10 % of its limit in "
                     "steady state.\n")

    sc, v2 = _load("suff_circuit"), _load("verify2")
    if sc:
        L.append("## Current-aware sufficient bound (post hoc)\n")
        L.append("Defined in the protocol's addendum (`inner_circuit.py`). It certifies a trajectory "
                 "of the averaged circuit with:\n")
        import inner_circuit as ic
        L.append("- the pole-N current within 97 % of its limit;")
        L.append(f"- the bus {100 * ic.BUS_HEADROOM:.0f} % inside its limits (raised from 1 % during the "
                 f"S4(a) checks; see the protocol addendum);")
        L.append("- the buffer current changing by at most 36 A/ms;")
        L.append("- a search that starts at the steady-state minimum with the inner headroom, so the "
                 "plan ends in a steady state a controller can hold;")
        L.append("- the pole-N band with the inner headroom.\n")
        L.append("Main chart, centred operating point. Each cell gives Cb_suff,c in µF and, in brackets, "
                 "its ratio to the larger of Cb_nec and the steady-state minimum. Each value is certified "
                 "at the worst phase and eight more phases 45° apart.\n")
        byid = {s["id"]: s for s in sc}
        for r in RATINGS:
            rows = [p for p in _sel(pts, rating=r, group="chart") if p["m"] <= 0.45]
            taus = sorted({round(_tau(p), 4) for p in rows})
            L.append(f"**{NAMES[r]}**\n")
            L.append("| m | " + " | ".join(f"{t:.2f}" for t in taus) + " |")
            L.append("|---|" + "---|" * len(taus))
            for m in M_GRID[:-1]:
                cells = []
                for t in taus:
                    p = next(q for q in rows if q["m"] == m and round(_tau(q), 4) == t)
                    s = byid.get(p["id"])
                    if s is None:
                        cells.append("—")
                    elif not np.isfinite(s["cb_suff_c"]):
                        cells.append("not certified" if s.get("note") != "no design" else "no design")
                    else:
                        cells.append(f"{s['cb_suff_c'] * 1e6:.0f} ({s['cb_suff_c'] / max(p['cb_nec'], p['cb_steady']):.2f})")
                L.append(f"| {m:.2f} | " + " | ".join(cells) + " |")
            L.append("")
        n_ok = sum(np.isfinite(s["cb_suff_c"]) for s in sc)
        L.append(f"Certified at {n_ok} of {len(sc)} points.\n")
        L.append("- **\"Not certified\":** the search found no certificate up to three times the base "
                 "value. The method is local, so this is not a proof that none exists.")
        L.append("- **At the base design with the centred operating point,** every attempt at both ratings "
                 "leaves a pole-N current overshoot of 1–4 A near the line peak, 2.5–3 ms after the "
                 "step. That is consistent with the step being limited by pole current there, not by "
                 "capacitance.\n")
    cols = " | ".join(TRACKER_NAMES[t].split(" (")[0].capitalize() for t in TRACKERS_ORDER)
    dashes = "|---" * len(TRACKERS_ORDER)
    a1 = _load("alpha1")
    if a1:
        L.append("**Top operating point (α = 1), base design.** The same bound, with the buffer held "
                 "at the top of its band before the step. The last columns track the certified plan "
                 "in the switched model:\n")
        L.append(f"| Rating | Cb_nec (α = 1) | Cb_suff,c (α = 1) | Ratio | {cols} |\n"
                 f"|---|---|---|---{dashes}|")
        for s in a1["suff"]:
            p = next(q for q in pts if q["id"] == s["id"])
            cs = s["cb_suff_c"]
            if not np.isfinite(cs):
                L.append(f"| {NAMES[s['rating']]} | {_uf(p['cb_nec'])} | not certified | — |"
                         + " — |" * len(TRACKERS_ORDER))
                continue
            vs = {_pol(v): v for v in a1["verify"] if v["rating"] == s["rating"]}
            L.append(f"| {NAMES[s['rating']]} | {_uf(p['cb_nec'])} | {_uf(cs)} | "
                     f"{cs / max(p['cb_nec'], p['cb_steady']):.2f} | "
                     + " | ".join(_track_res(vs.get(t)) for t in TRACKERS_ORDER) + " |")
        L.append("")
    if v2:
        L.append("**S4(a) rerun with the certified plans, centred operating point.** No base design is "
                 "certified there, so each rating's certified point at the bus-loop ramp with the largest m "
                 "is tested (worst phase, switched model, plan continued by repeating its last cycle):\n")
        L.append(f"| Rating | m | Cb | Plan's peak pole-N current | {cols} |\n|---|---|---|---{dashes}|")
        for r in RATINGS:
            vs = {_pol(v): v for v in v2 if v["rating"] == r}
            v = vs.get(TRACKER) or next(iter(vs.values()), None)
            if v is None:
                continue
            if "error" in v:
                L.append(f"| {NAMES[r]} | — | — | — | {v['error']} |" + " — |" * (len(TRACKERS_ORDER) - 1))
                continue
            L.append(f"| {NAMES[r]} | {v.get('m', 0.40):.2f} | {_uf(v['cb'])} | {v['plan_peak_iN']:.1f} A | "
                     + " | ".join(_track_res(vs.get(t)) for t in TRACKERS_ORDER) + " |")
        L.append("")
    vcv = _load("verify_conv")
    if vcv:
        ok = all(x["violation"] for x in vcv)
        L.append("**S4(b) and (c) with the timing fix.** The conventional controller with the final "
                 "tracker's buffer-loop timing " + ("still trips" if ok else "does not always trip")
                 + ": " + "; ".join(f"{NAMES[x['rating']]}, {_uf(x['cb'])}: "
                                    + (f"trip at {x['t_violation_ms']:.1f} ms" if x["violation"] else "held")
                                    for x in vcv) + ".\n")
    vc = _load("verify_chart")
    if vc:
        L.append("**S4(a) across the main chart.** Every certified point, centred operating point, worst "
                 "phase, tracked in the switched model:\n")
        L.append("| Rating | Tracker | Held | Tested | Smallest bus margin | Within 1 % of a bus limit | "
                 "Largest pole-N current |\n|---|---|---|---|---|---|---|")
        for r in RATINGS:
            for pol in TRACKERS_ORDER:
                v = [x for x in vc if x["rating"] == r and _pol(x) == pol and "error" not in x]
                if not v:
                    continue
                mg = [_margin(x) for x in v]
                thin = sum(0 <= g_ < 0.01 * np.mean(x["Vdc_lim"]) for g_, x in zip(mg, v))
                L.append(f"| {NAMES[r]} | {TRACKER_NAMES[pol]} | {sum(not x['violation'] for x in v)} | "
                         f"{len(v)} | {min(mg):+.1f} V | {thin} | {max(x['I_pole_peak'] for x in v):.1f} A "
                         f"(limit {v[0]['I_pole_max']:.0f} A) |")
        bad = [x for x in vc if _pol(x) == TRACKER and "error" not in x and x["violation"]]
        err = [x for x in vc if _pol(x) == TRACKER and "error" in x]
        if bad:
            L.append("\nFailures with the " + TRACKER_NAMES[TRACKER] + ": " + "; ".join(
                f"{NAMES[x['rating']]}, m = {x['m']:.2f}, ramp {_ramp_txt(x)}: {_track_res(x)}"
                for x in bad) + ".")
        if err:
            L.append(f"\n{len(err)} points were not certified again on the rerun.")
        L.append("")
    if v2 or vc:
        L.append("- **Trackers.** All three follow the plan's DC1 power and stored energy. They differ in "
                 "the bus reference:")
        L.append("  - the **energy tracker** follows the plan's bus sample by sample and feeds its rate "
                 "forward. The plans move the bus by about ±5 V at 120 Hz in steps on their 83 µs grid, "
                 "and those steps ask for buffer-power pulses that saturate pole N;")
        L.append("  - the **flat-bus tracker** holds the bus at nominal. It fails where the plan needs the "
                 "bus capacitor's share of the ripple: fast DC1 ramps at the steady-state minimum;")
        L.append("  - the **smoothed-bus tracker** (policy `plan_smooth`) follows a 15-sample (1.25 ms) "
                 "moving average of the plan's bus. It keeps the plan's 120 Hz use of the bus capacitor "
                 "and drops the steps.")
        L.append("- **Buffer-loop timing.** The smoothed-bus tracker also compares the buffer current "
                 "with its reference at the measurement instant (gain `sync_ib`), as the AC loops do. The "
                 "earlier loop compared the measured current with the reference for the period the duty "
                 "acts in, 1.5 switching periods later. The buffer current then ran about 110 µs ahead of "
                 "the ripple it carries, which a 23 µF bus feels at once.")
        if v2:
            vs = {_pol(v): v for v in v2 if v["rating"] == "r120" and "error" not in v}
            if "smooth" in vs and "flat" in vs:
                nom = np.mean(vs["smooth"]["Vdc_lim"])
                L.append(f"  At the 120 V S4 point the bus swing went from "
                         f"{vs['flat']['Vdc_min'] - nom:+.0f}/{vs['flat']['Vdc_max'] - nom:+.0f} V "
                         f"(flat-bus tracker) to {vs['smooth']['Vdc_min'] - nom:+.0f}/"
                         f"{vs['smooth']['Vdc_max'] - nom:+.0f} V.")
        if vc:
            f_ = [x for x in vc if _pol(x) == TRACKER and "error" not in x]
            held = [x for x in f_ if not x["violation"]]
            if held:
                exc = {r: [100 * max(x["Vdc_max"] - np.mean(x["Vdc_lim"]), np.mean(x["Vdc_lim"]) - x["Vdc_min"])
                           / np.mean(x["Vdc_lim"]) for x in held if x["rating"] == r] for r in RATINGS}
                thin = [x for x in held if _margin(x) < 0.01 * np.mean(x["Vdc_lim"])]
                L.append(f"- **Margins.** On the plans that held, the bus stays within "
                         f"{max(exc['r120']):.1f} % of nominal at 120 V and {max(exc['r230']):.1f} % at 230 V, "
                         f"against limits of ±5 %. {len(thin)} of the {len(held)} came within 1 % of a "
                         f"limit; the closest margin is {min(_margin(x) for x in held):.2f} V.")
            bad = [x for x in f_ if x["violation"]]
            if bad and all(x["m"] == 0.45 for x in bad):
                L.append("- **The failures are at m = 0.45 with a 0.5 ms ramp,** the edge of feasibility. "
                         "There the plans swing DC1 by several kilowatts within a few milliseconds while "
                         "the buffer runs low, and the tracked buffer cannot follow. In a side test no "
                         "tracker held at 120 V up to 1507 µF, twice the chart value. At 230 V the "
                         "smoothed-bus tracker held at 1163 µF (1.9 times), with a 1.3 V margin. These "
                         "points need a tracker that follows the plan's full state.")
        L.append("- **Variants tried and dropped, all worse:**")
        L.append("  - following the plan's buffer current (`circuit_plan`) or its pole-N voltage (`state_plan`);")
        L.append("  - feeding forward the buffer inductor's power;")
        L.append("  - extrapolating DC1 power to the duty's period;")
        L.append("  - a faster bus loop: 400 Hz helped, 600 Hz was unstable;")
        L.append("  - an integral term on the bus loop;")
        L.append("  - in the bound, 5–10 % DC1 headroom;")
        L.append("  - in the bound, a final-cycle bus band. That was infeasible, because by the final cycle "
                 "the bound's bus-energy bracket spans 4–5 % of nominal at 120 V.\n")

    L.append("## Closed-form approximation (post hoc)\n")
    L.append(f"Fitted to the 120 V chart and tested on the 230 V chart, for ramp times of at least "
             f"0.1 line cycle: penalty (%) ≈ {fit['a']:.1f} + {fit['b']:.1f} ln(T_r f) + "
             f"{fit['c']:.0f} (m − 0.4). Fit error: largest {fit['fit_max_abs_pts']:.1f} points, "
             f"RMS {fit['fit_rms_pts']:.1f}. Test on 230 V: largest {fit['test_max_abs_pts']:.1f} "
             f"points, RMS {fit['test_rms_pts']:.1f}. Range of validity: m from 0.30 to 0.40, ramp "
             f"times from 0.1 to 1.8 line cycles, the rated step, unity power factor and the centred "
             f"operating point.\n")
    (RES / "wp5a_sizing.md").write_text("\n".join(L), encoding="utf-8")


if __name__ == "__main__":
    main()
