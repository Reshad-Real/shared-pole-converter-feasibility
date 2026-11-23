"""Paper figures at print size (Results in Engineering, Elsevier CAS single-column template).

Every figure is drawn at the template's text width (468 pt) and included at 100 %, so the font
sizes set here are the printed sizes: 9 pt for ticks, legends and annotations, 9.5 pt for axis
labels and panel titles. All text is black. Output: vector PDF in paper/figs/.

    py paper/make_figures.py                 # every figure
    py paper/make_figures.py boundary case   # selected figures

Figures and the results they read:
  pair       sim/results/r120/results.json (time-domain runs are recomputed)
  boundary   sim/results/paper/phase_sweep_circuit.json, sim/results/<rating>/results.json
  simscape   sim_simscape/results/simscape_results.json, sim/results/results.json
  r22        sim_r22/results/wp3_r22_check.json
  twin       sim/results/<rating>/wp9_twin.json
  chart      sim/results/wp5a/wp5a_grid.json
  alpha      sim/results/wp5a/wp5a_grid.json
  case       sim/results/wp6/wp6_events.json
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SIM = ROOT / "sim"
sys.path.insert(0, str(SIM))
OUT = Path(__file__).resolve().parent / "figs"

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

W = 468.33 / 72.27                      # text width of the template (in)
# Categorical slots 1-3 of the validated reference palette, plus black for references
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
BLACK = "#000000"
BLUE_LIGHT, BAND, GRID = "#cde2fb", "#e4e4e4", "#d0d0d0"
GOOD, CRIT = "#008300", "#e34948"
SEQ_BLUE = ("#ffffff", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
RAMP_M = ("#6da7ec", "#2a78d6", "#1c5cab", "#0d366b")      # ordinal: modulation index
RAMP_A = ("#f7a072", "#eb6834", "#b4461b")                 # ordinal: operating point
MARKERS = ("o", "s", "^", "D")
NAMES = {"r120": "120 V rating", "r230": "230 V rating"}
RATINGS = ("r120", "r230")


def style():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["STIXGeneral", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 9, "axes.labelsize": 9.5,
        "axes.titlesize": 9.5, "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
        "text.color": BLACK, "axes.labelcolor": BLACK, "axes.titlecolor": BLACK,
        "xtick.color": BLACK, "ytick.color": BLACK, "axes.edgecolor": BLACK,
        "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
        "xtick.major.size": 3, "ytick.major.size": 3,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
        "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True,
        "lines.linewidth": 1.3, "legend.frameon": False, "legend.handlelength": 2.2,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "pdf.fonttype": 42, "figure.constrained_layout.h_pad": 0.03,
        "figure.constrained_layout.w_pad": 0.03,
    })


def save(fig, name):
    path = OUT / f"{name}.pdf"
    fig.savefig(path)
    fig.savefig(OUT / f"{name}.png", dpi=200)
    plt.close(fig)
    print("wrote", path)


def _deg(x):
    return f"{x:.0f}°"


# ---------------------------------------------------------------- Fig. 2: equal-energy pair
def fig_pair():
    import ratings
    from feasibility import nominal_orbit, solve
    from timedomain import simulate
    p = ratings.get("r120")
    R = json.loads((SIM / "results" / "r120" / "results.json").read_text())
    pr = R["pair"]
    th_a, th_b, dP = pr["theta_a_deg"], pr["theta_b_deg"], pr["dP_test_W"]
    gains = dict(k_rec=2 * np.pi * R["timedomain_tuning"]["k_rec_hz_chosen"])
    plan_a = solve(p, np.deg2rad(th_a), dP, "inner", orbit=nominal_orbit(p), smooth=True)
    runs = [(rf"Conventional, $\theta_a$ = {_deg(th_a)}", simulate(p, np.deg2rad(th_a), dP, gains=gains)),
            (rf"Conventional, $\theta_b$ = {_deg(th_b)}", simulate(p, np.deg2rad(th_b), dP, gains=gains)),
            (rf"LP plan tracked, $\theta_a$ = {_deg(th_a)}",
             simulate(p, np.deg2rad(th_a), dP, policy="plan", plan=plan_a, gains=gains))]
    fig, axs = plt.subplots(3, 3, figsize=(W, 5.0), sharex=True, sharey="row", layout="constrained")
    t_lo, t_hi = -5.0, 40.0
    for c, (title, s) in enumerate(runs):
        t = s["t"] * 1e3
        k = (t >= t_lo) & (t <= t_hi)
        a0, a1, a2 = axs[:, c]
        a0.fill_between(t[k], s["lo"][k], s["hi"][k], color=BAND, lw=0, label="Pole-N band")
        a0.plot(t[k], s["vb"][k], color=BLUE, label="Buffer voltage $v_b$")
        a1.plot(t[k], s["Vdc"][k], color=BLUE)
        for v in (p.Vdc_min, p.Vdc_max):
            a1.axhline(v, color=BLACK, lw=0.9, ls=(0, (4, 2)),
                       label="Bus limits (±5 %)" if v == p.Vdc_min and c == 0 else None)
        a2.plot(t[k], s["pd"][k] / 1e3, color=BLUE)
        m = s["metrics"]
        verdict = "violation, trip" if m["violation"] else "all limits held"
        a0.set_title(f"{title}\n{verdict}", loc="left")
        if m["t_violation_ms"] is not None:
            for a in axs[:, c]:
                a.axvline(m["t_violation_ms"], color=CRIT, lw=1.1)
            a1.annotate("first\nviolation", (m["t_violation_ms"], p.Vdc_min - 8),
                        xytext=(4, 0), textcoords="offset points", va="top", fontsize=9)
            print(f"pair {title}: violation at {m['t_violation_ms']:.1f} ms")
        a2.set_xlabel("Time after step (ms)")
        a2.set_xlim(t_lo, t_hi)
        a2.set_xticks([0, 10, 20, 30, 40])
    axs[0, 0].set_ylabel("Voltage (V)")
    axs[1, 0].set_ylabel("Bus voltage (V)")
    axs[2, 0].set_ylabel("DC1 power (kW)")
    axs[1, 0].set_ylim(p.Vdc_min - 30, p.Vdc_max + 10)
    h, l = axs[0, 0].get_legend_handles_labels()
    h2, l2 = axs[1, 0].get_legend_handles_labels()
    fig.legend(h + h2, l + l2, loc="outside upper center", ncol=3)
    fig.align_ylabels(axs[:, 0])
    save(fig, "fig_pair_timedomain")


# ---------------------------------------------------------------- Fig. 3: boundary against phase
def fig_boundary():
    import ratings
    from feasibility import nominal_orbit
    from model import ac_state, band
    sweep = json.loads((SIM / "results" / "paper" / "phase_sweep_circuit.json").read_text())
    fig, axs = plt.subplots(2, 2, figsize=(W, 5.3), sharex=True, layout="constrained",
                            gridspec_kw=dict(height_ratios=[2.6, 1.2]))
    summary = {}
    for c, rating in enumerate(RATINGS):
        R = json.loads((SIM / "results" / rating / "results.json").read_text())
        sw = R["sweep"]
        ph = np.array(sw["phase_deg"])

        def get(kind):
            return np.array([next(x["step_W"] for x in sweep if x["rating"] == rating
                                  and x["kind"] == kind and x["phase_deg"] == q) for q in ph])
        outer, suff = get("outer"), get("suff_scan")
        conv = np.minimum(np.array(sw["energy_W"]), R["steady"]["C2_limit_W"])
        ax = axs[0, c]
        ax.fill_between(ph, np.nan_to_num(suff) / 1e3, outer / 1e3, color=BLUE_LIGHT, lw=0,
                        label="Between the two bounds")
        ax.plot(ph, outer / 1e3, color=BLUE, lw=1.6,
                label="Necessary bound (nothing admissible above)")
        ax.plot(ph, suff / 1e3, color=BLUE, lw=1.1, ls=(0, (4, 2)),
                label="Sufficient bound (certified below)")
        ax.plot(ph, np.array(sw["energy_W"]) / 1e3, color=ORANGE, label="Energy sizing")
        ax.plot(ph, np.array(sw["static_W"]) / 1e3, color=AQUA, ls=(0, (6, 2, 1, 2)),
                label="Static derating")
        ax.axhline(R["steady"]["C2_limit_W"] / 1e3, color=BLACK, lw=1.0, ls=(0, (1, 1.5)),
                   label="Steady-state and instantaneous checks")
        pa, pb, dPt = R["pair"]["theta_a_deg"], R["pair"]["theta_b_deg"], R["pair"]["dP_test_W"]
        ax.plot([pa], [dPt / 1e3], "o", ms=7, color=GOOD, mec="white", mew=1.2, zorder=5)
        ax.plot([pb], [dPt / 1e3], "X", ms=8, color=CRIT, mec="white", mew=1.0, zorder=5)
        ax.annotate(r"$\theta_a$", (pa, dPt / 1e3), xytext=(6, -12), textcoords="offset points")
        ax.annotate(r"$\theta_b$", (pb, dPt / 1e3), xytext=(6, -12), textcoords="offset points")
        ax.set_ylim(0, sw["search_cap_W"] * 1.04 / 1e3)
        ax.set_title(f"{NAMES[rating]} (AC2 step of {dPt:.0f} W tested)", loc="left")
        if c == 0:
            ax.set_ylabel("Largest admissible AC2 step (kW)")
        over = conv - outer
        k = int(np.nanargmax(over))
        summary[rating] = dict(worst_phase_deg=float(ph[k]), conventional_W=float(conv[k]),
                               outer_W=float(outer[k]), over_pct=float(100 * over[k] / conv[k]),
                               suff_at_worst_W=float(suff[k]),
                               outer_at_pair=dict(a=float(outer[list(ph).index(pa)]),
                                                  b=float(outer[list(ph).index(pb)])),
                               suff_at_pair=dict(a=float(suff[list(ph).index(pa)]),
                                                 b=float(suff[list(ph).index(pb)])),
                               max_gap_W=float(np.nanmax(outer - suff)))
        p = ratings.get(rating)
        orb = nominal_orbit(p)
        th_o = np.rad2deg(orb["t"] * p.w)
        u1, u2, _, _, _ = ac_state(p, orb["t"])
        lo, m = band(u1, u2)
        a2 = axs[1, c]
        a2.fill_between(th_o, lo, p.Vdc_nom - m, color=BAND, lw=0, label="Pole-N band")
        a2.plot(th_o, np.sqrt(2 * orb["Eb"] / p.Cb), color=BLUE, label="Buffer voltage, pre-step orbit")
        for th in (pa, pb):
            for a in (ax, a2):
                a.axvline(th, color="#9a9a9a", lw=0.7)
        if c == 0:
            a2.set_ylabel("Voltage (V)")
        a2.set_xlabel("AC1 phase at which the step arrives (deg)")
        a2.set_xticks(np.arange(0, 361, 90))
        a2.set_xlim(0, 360)
    h, l = axs[0, 0].get_legend_handles_labels()
    h2, l2 = axs[1, 0].get_legend_handles_labels()
    fig.legend(h + h2, l + l2, loc="outside upper center", ncol=2)
    fig.align_ylabels(axs[:, 0])
    save(fig, "fig_boundary")
    (Path(__file__).parent / "boundary_summary.json").write_text(
        json.dumps(summary, indent=1))


# ---------------------------------------------------------------- Fig. 4: Simscape against Python
def fig_simscape():
    from feasibility import nominal_orbit, solve
    from model import Params
    from switched import SwitchParams, simulate_switched
    R = json.loads((SIM / "results" / "results.json").read_text())
    M = json.loads((ROOT / "sim_simscape" / "results" / "simscape_results.json").read_text())
    prm = {k: (tuple(v) if isinstance(v, list) else v) for k, v in R["params"].items()}
    p = Params(**prm)
    pr = R["pair"]
    gains = dict(k_rec=2 * np.pi * R["timedomain_tuning"]["k_rec_hz_chosen"])
    plan = solve(p, np.deg2rad(pr["theta_a_deg"]), pr["dP_test_W"], "inner",
                 orbit=nominal_orbit(p), smooth=True)
    cases = {"conventional_theta_a": (pr["theta_a_deg"], "conventional",
                                      rf"Conventional, $\theta_a$ = {_deg(pr['theta_a_deg'])}"),
             "conventional_theta_b": (pr["theta_b_deg"], "conventional",
                                      rf"Conventional, $\theta_b$ = {_deg(pr['theta_b_deg'])}"),
             "plan_theta_a": (pr["theta_a_deg"], "plan",
                              rf"LP plan tracked, $\theta_a$ = {_deg(pr['theta_a_deg'])}")}
    fig, axs = plt.subplots(2, 3, figsize=(W, 3.9), sharex=True, sharey="row", layout="constrained")
    Tsw = 50e-6
    for c, (key, (th, pol, title)) in enumerate(cases.items()):
        py = simulate_switched(p, SwitchParams(), np.deg2rad(th), pr["dP_test_W"], policy=pol,
                               plan=plan if pol == "plan" else None, gains=gains)
        s = M["pair"][key]
        ts = 1e3 * (np.asarray(s["ts"]["t"]) - 0.5 * Tsw)
        tp = 1e3 * (np.asarray(py["t"]) + 0.5 * Tsw)
        ks = (ts >= -5) & (ts <= 40) & (np.asarray(s["ts"]["t"]) > -0.019)
        kp = (tp >= -5) & (tp <= 40)
        axs[0, c].plot(tp[kp], py["vb_avg"][kp], color=ORANGE, lw=1.5, label="Python switched model")
        axs[0, c].plot(ts[ks], np.asarray(s["ts"]["vb_avg"])[ks], color=BLUE, lw=1.0,
                       ls=(0, (3, 1.5)), label="Simscape Electrical model")
        for v in (p.Vdc_min, p.Vdc_max):
            axs[1, c].axhline(v, color=BLACK, lw=0.9, ls=(0, (4, 2)))
        axs[1, c].plot(tp[kp], py["Vdc_avg"][kp], color=ORANGE, lw=1.5)
        axs[1, c].plot(ts[ks], np.asarray(s["ts"]["Vdc_avg"])[ks], color=BLUE, lw=1.0, ls=(0, (3, 1.5)))
        for m, col in ((py["metrics"], ORANGE), (s["metrics"], BLUE)):
            if m["violation"]:
                axs[1, c].axvline(m["t_violation_ms"], color=col, lw=0.9)
        axs[0, c].set_title(title, loc="left")
        axs[1, c].set_xlabel("Time after step (ms)")
        axs[1, c].set_xlim(-5, 40)
        axs[1, c].set_xticks([0, 10, 20, 30, 40])
    axs[0, 0].set_ylabel("Buffer voltage (V)")
    axs[1, 0].set_ylabel("Bus voltage (V)")
    axs[1, 0].set_ylim(p.Vdc_min - 25, p.Vdc_max + 10)
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside upper center", ncol=2)
    fig.align_ylabels(axs[:, 0])
    save(fig, "fig_simscape_vs_python")


# ---------------------------------------------------------------- Fig. 5: published experiment
def _r22_label(t):
    t = t.replace("u_cac", "$u_{cac}$").replace("u_ac", "$u_{ac}$").replace("S_a", "$S_a$")
    t = re.sub(r"\bi_([abc])\b", r"$i_\1$", t)
    return t.replace(" vs ", " vs. ")


def fig_r22():
    s = json.loads((ROOT / "sim_r22" / "results" / "wp3_r22_check.json").read_text())
    rows = [r for r in s["results"] if r["tolerance"] is not None]
    fig, ax = plt.subplots(figsize=(W, 0.24 * len(rows) + 1.2), layout="constrained")
    for i, r in enumerate(rows[::-1]):
        m = r["measured"]
        lo, hi = 1 - r["tolerance"] / abs(m), 1 + r["tolerance"] / abs(m)
        ax.plot([lo, hi], [i, i], color=BLUE_LIGHT, lw=8, solid_capstyle="butt", zorder=1,
                label="Predeclared pass band" if i == 0 else None)
        if r["uncertainty"]:
            ax.plot([1 - r["uncertainty"] / abs(m), 1 + r["uncertainty"] / abs(m)], [i, i],
                    color=BLACK, lw=1.0, zorder=2,
                    label="Digitisation uncertainty" if i == len(rows) - 1 else None)
        for mname, mk, col, off in (("Python switched", "o", ORANGE, 0.15),
                                    ("Simscape", "D", BLUE, -0.15)):
            v = r["models"].get(mname)
            if v:
                ax.plot([v["value"] / m], [i + off], mk, ms=5.5, color=col, mec="white", mew=0.8,
                        zorder=3, label=(mname + " model") if i == 0 else None)
    ax.axvline(1.0, color=BLACK, lw=0.8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([_r22_label(r["target"]) for r in rows[::-1]])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Model value / measured value")
    ax.set_ylim(-0.7, len(rows) - 0.3)
    fig.legend(loc="outside upper center", ncol=4, handlelength=1.6, columnspacing=1.2)
    save(fig, "fig_r22_benchmark")


# ---------------------------------------------------------------- Fig. 7: digital twin
def fig_twin():
    import wp9_report as w9
    rated = w9.RATED
    data = {r: json.loads((SIM / "results" / r / "wp9_twin.json").read_text()) for r in RATINGS}
    truth = {r: w9.truth_labels(r) for r in RATINGS}
    cases = {r: [c for c in data[r]["cases"] if "error" not in c and c["index"] in truth[r]]
             for r in RATINGS}
    fig, axs = plt.subplots(2, 2, figsize=(W, 5.6), layout="constrained")
    cols = {"r120": BLUE, "r230": ORANGE}
    off = {"r120": -0.17, "r230": 0.17}

    def box(ax, pos, vals, col):
        b = ax.boxplot([vals], positions=[pos], widths=0.28, whis=(5, 95), showfliers=False,
                       patch_artist=True, medianprops=dict(color=BLACK, lw=1.0),
                       whiskerprops=dict(color=col, lw=1.0), capprops=dict(color=col, lw=1.0),
                       boxprops=dict(facecolor="white", edgecolor=col, lw=1.2))
        return b

    # (a) identification
    ax = axs[0, 0]
    params = [("Cb", "$C_b$"), ("Lb", "$L_b$"), ("E1", "$E_1$"), ("slew_dc1", "DC1 slew")]
    for r in RATINGS:
        for i, (k, _) in enumerate(params):
            box(ax, i + off[r], [100 * c["param_error"][k] for c in cases[r]], cols[r])
    for i, (k, _) in enumerate(params):
        lim = 100 * w9.TW1_LIMIT[k]
        ax.plot([i - 0.4, i + 0.4], [lim, lim], color=BLACK, lw=0.9, ls=(0, (3, 1.5)))
        ax.plot([i - 0.4, i + 0.4], [-lim, -lim], color=BLACK, lw=0.9, ls=(0, (3, 1.5)))
    ax.set_xticks(range(len(params)))
    ax.set_xticklabels([n for _, n in params])
    ax.axhline(0, color=BLACK, lw=0.6)
    ax.set_ylabel("Identification error (%)")
    ax.set_title("(a) Identified parameters", loc="left")
    ax.grid(axis="x", visible=False)

    # (b) boundary error
    ax = axs[0, 1]
    preds = [("N", "Nominal"), ("NE", "Nominal +\nmeas. $E$"), ("T", "Twin")]
    for r in RATINGS:
        for i, (k, _) in enumerate(preds):
            box(ax, i + off[r], [100 * (c["max_step_W"][k] - c["max_step_W"]["true"]) / rated[r]
                                 for c in cases[r]], cols[r])
    ax.axhline(0, color=BLACK, lw=0.6)
    ax.set_xticks(range(len(preds)))
    ax.set_xticklabels([n for _, n in preds])
    ax.set_ylabel("Boundary error (% of AC2 rating)")
    ax.set_title("(b) Largest admissible step", loc="left")
    ax.grid(axis="x", visible=False)

    # (c) drift monitor
    ax = axs[1, 0]
    rng = np.random.default_rng(0)
    for j, (key, lab) in enumerate((("residual_unaged_V", "unaged"), ("residual_aged_V", "aged"))):
        for r in RATINGS:
            y = np.array([c["drift"][key] / c["drift"]["threshold_V"] for c in cases[r]])
            x = j + off[r] + rng.uniform(-0.08, 0.08, len(y))
            ax.plot(x, y, "o", ms=2.6, color=cols[r], alpha=0.55, mec="none")
    ax.axhline(1.0, color=BLACK, lw=1.0, ls=(0, (4, 2)))
    ax.text(1.5, 1.04, "alarm\nthreshold", ha="right", va="bottom", fontsize=9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Unaged plant", "Aged capacitor\n(−15 % $C_b$, 2× $R_b$)"])
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(0, None)
    ax.set_ylabel("Monitor residual / threshold")
    ax.set_title("(c) Drift monitor", loc="left")
    ax.grid(axis="x", visible=False)

    # (d) predicted against true largest step
    ax = axs[1, 1]
    for r in RATINGS:
        tr = np.array([c["max_step_W"]["true"] / rated[r] for c in cases[r]])
        for k, mk, ms in (("N", "x", 3.5), ("T", "o", 2.6)):
            pr = np.array([c["max_step_W"][k] / rated[r] for c in cases[r]])
            ax.plot(tr, pr, mk, ms=ms, color=cols[r], alpha=0.6 if k == "T" else 0.45,
                    mec=cols[r] if k == "N" else "none", mew=0.7)
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color=BLACK, lw=0.8)
    ax.set_xlabel("True largest step / AC2 rating")
    ax.set_ylabel("Predicted / AC2 rating")
    ax.set_title("(d) Nominal (×) and twin (●)", loc="left")
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=BLUE, marker="s", ls="", ms=7, mfc="white", mew=1.2,
                      label="120 V rating"),
               Line2D([], [], color=ORANGE, marker="s", ls="", ms=7, mfc="white", mew=1.2,
                      label="230 V rating"),
               Line2D([], [], color=BLACK, ls=(0, (3, 1.5)), lw=0.9,
                      label="TW1 limit (a)")]
    fig.legend(handles=handles, loc="outside upper center", ncol=3)
    save(fig, "fig_twin_results")


# ---------------------------------------------------------------- Figs. 8-9: sizing
def fig_chart():
    import wp5a_report as w5
    pts = w5._load("grid")
    fig, axs = plt.subplots(1, 2, figsize=(W, 3.1), sharey=True, layout="constrained")
    for ax, r in zip(axs, RATINGS):
        taus, tab = w5.chart_table(pts, r)
        for i, m in enumerate(w5.M_GRID[:-1]):
            y = np.where(tab[i] > 175, np.nan, tab[i])
            ax.plot(taus, y, color=RAMP_M[i], marker=MARKERS[i], ms=4.5, label=f"$m$ = {m:.2f}")
        off = tab[3][np.isfinite(tab[3]) & (tab[3] > 175)]
        if off.size:
            ax.text(0.03, 0.97, f"$m$ = 0.45: {off.min():.0f}–{off.max():.0f} % off scale",
                    transform=ax.transAxes, va="top", fontsize=9)
        ax.set_ylim(-60, 180)
        ax.axhline(w5.S1_PCT, color=BLACK, lw=0.9, ls=(0, (4, 2)))
        ax.axhline(0, color=BLACK, lw=0.7)
        ax.axvline(1 / (np.pi * 0.4), color=BLACK, lw=0.8, ls=(0, (1, 1.5)))
        ax.text(1 / (np.pi * 0.4) * 1.08, -52, "bus\nloop", fontsize=9, va="bottom")
        ax.set_xscale("log")
        ax.set_title(NAMES[r], loc="left")
        ax.set_xlabel("DC1 ramp time for the rated step (line cycles)")
    axs[0].text(0.035, w5.S1_PCT + 3, "15 % mark", fontsize=9, va="bottom")
    axs[0].set_ylabel("Capacitance penalty (%)")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside upper center", ncol=4)
    save(fig, "fig_design_chart")


def fig_alpha():
    import wp5a_report as w5
    pts = w5._load("grid")
    fig, axs = plt.subplots(1, 2, figsize=(W, 3.2), sharey=True, layout="constrained")
    for ax, r in zip(axs, RATINGS):
        for k, a in enumerate((0.5, 0.75, 1.0)):
            grp = "chart" if a == 0.5 else "alpha"
            rows = sorted([q for q in w5._sel(pts, rating=r, group=grp, alpha=a) if q["m"] == 0.40],
                          key=w5._tau)
            x = [w5._tau(q) for q in rows]
            ref = rows[0]["cb_steady"]
            ax.plot(x, [q["cb_nec"] / ref for q in rows], color=RAMP_A[k], marker=MARKERS[k],
                    ms=4.5, label=rf"Necessary, $\alpha$ = {a:g}")
            ax.plot(x, [q["cb_conv"] / ref for q in rows], color=RAMP_A[k], lw=1.0,
                    ls=(0, (4, 2)), label=rf"Phase-blind, $\alpha$ = {a:g}")
        ax.axhline(1.0, color=BLACK, lw=0.7)
        ax.set_xscale("log")
        ax.set_title(NAMES[r], loc="left")
        ax.set_xlabel("DC1 ramp time for the rated step (line cycles)")
    axs[0].set_ylabel("$C_b$ / steady-state minimum")
    h, l = axs[0].get_legend_handles_labels()     # necessary and phase-blind alternate: two rows
    fig.legend(h, l, loc="outside upper center", ncol=3)
    save(fig, "fig_operating_point")


# ---------------------------------------------------------------- Fig. 10: case study
def fig_case():
    import wp6
    import wp6_report as w6
    events = json.loads((w6.RES / "wp6_events.json").read_text())
    D = wp6.designs()
    cols = [("D1", "unavoidable"), ("D1", "conv"), ("D3", "unavoidable"), ("D3", "conv"),
            ("D2", "unavoidable"), ("D2", "conv"), ("D2", "cert"),
            ("D4", "unavoidable"), ("D4", "conv"), ("D4", "cert")]
    cols = [c for c in cols if c[0] in D]
    head = {"unavoidable": "bound", "conv": "conv.", "cert": "cert."}
    types = w6._types()
    cmap = LinearSegmentedColormap.from_list("seq_blue", SEQ_BLUE)
    fig, axes = plt.subplots(2, 1, figsize=(W, 7.4), layout="constrained")
    for ax, rating in zip(axes, RATINGS):
        A = np.full((len(types), len(cols)), np.nan)
        for c, (d, kind) in enumerate(cols):
            Mx = w6.p_matrix(events, d, rating, kind)
            for r, (i, j) in enumerate(types):
                A[r, c] = 100 * Mx[i][j]
        im = ax.imshow(A, cmap=cmap, vmin=0, vmax=100, aspect="auto")
        for r in range(A.shape[0]):
            for c in range(A.shape[1]):
                v = A[r, c]
                if np.isfinite(v):
                    ax.text(c, r, f"{v:.0f}", ha="center", va="center", fontsize=9,
                            color="white" if v > 55 else BLACK)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([f"{d}\n{head[k]}" for d, k in cols])
        ax.xaxis.tick_top()
        ax.set_yticks(range(len(types)))
        ax.set_yticklabels([w6._evtxt(i, j) for i, j in types])
        ax.set_ylabel("AC2 load before → after\n(fraction of rating)")
        ax.set_title(NAMES[rating], loc="left", pad=4)
        ax.tick_params(length=0)
        ax.grid(False)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xticks(np.arange(-0.5, len(cols)), minor=True)
        ax.set_yticks(np.arange(-0.5, len(types)), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.5)
        ax.tick_params(which="minor", length=0)
        for x in (1.5, 3.5, 6.5):
            ax.axvline(x, color=BLACK, lw=1.0)
    cb = fig.colorbar(im, ax=axes, shrink=0.5, pad=0.02)
    cb.set_label("Events that fail at a uniformly random phase (%)")
    cb.outline.set_visible(False)
    save(fig, "fig_case_events")


FIGS = dict(pair=fig_pair, boundary=fig_boundary, simscape=fig_simscape, r22=fig_r22,
            twin=fig_twin, chart=fig_chart, alpha=fig_alpha, case=fig_case)


def main():
    style()
    OUT.mkdir(exist_ok=True)
    names = sys.argv[1:] or list(FIGS)
    for n in names:
        FIGS[n]()


if __name__ == "__main__":
    main()
