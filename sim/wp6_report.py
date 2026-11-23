"""WP6 report: per-event failure probabilities, daily probabilities, criteria W1-W6 and fig16.

    py sim/wp6_report.py

Reads results/wp6/wp6_events.json and wp6_d4.json; writes results/wp6/wp6_case_study.md,
wp6_summary.json and fig16_wp6_events.png.
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

from sizing import design
import wp6

RES = wp6.RES6
RATINGS = ("r120", "r230")
NAMES = {"r120": "120 V rating", "r230": "230 V rating"}
LEVEL_TXT = ("0", "1/3", "2/3", "1")
SEQ = ("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")  # blue 100-700
INK, INK2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"
W1_PCT, W3_PCT = 10.0, 95.0


def _types():
    n = len(wp6.LEVELS)
    return [(i, j) for i in range(n) for j in range(n) if i != j]


def _evtxt(i, j):
    return f"{LEVEL_TXT[i]} → {LEVEL_TXT[j]}"


def outcome_fail(e, kind):
    """Did event e fail under outcome `kind`? None when not evaluated."""
    if "error" in e:
        return None
    if kind == "unavoidable":
        return not e["outer_feasible"]
    if kind == "conv":
        return bool(e["conv"]["violation"])
    if kind == "cert":
        c = e.get("cert")
        if c is None:
            return None
        return not (c["certified"] and not c["violation"])
    raise ValueError(kind)


def p_matrix(events, design, rating, kind, phases=None):
    """Failure probability per event type over the phases (all grid phases by default)."""
    n = len(wp6.LEVELS)
    M = [[0.0] * n for _ in range(n)]
    for i, j in _types():
        v = [outcome_fail(e, kind) for e in events
             if e["design"] == design and e["rating"] == rating and e["i"] == i and e["j"] == j
             and (phases is None or e["phase_deg"] in phases)]
        v = [x for x in v if x is not None]
        M[i][j] = float(np.mean(v)) if v else float("nan")
    return M


def _nan_free(M):
    return all(np.isfinite(M[i][j]) for i, j in _types())


def daily_table(events, design, rating, kind):
    """Daily probability under the base model and each sensitivity case."""
    M = p_matrix(events, design, rating, kind)
    if not _nan_free(M):
        return None
    zc = p_matrix(events, design, rating, kind, phases=(0, 180))
    w0 = {1: 1.0, 2: 0.0, 3: 0.0}
    w2 = {1: 1.0, 2: 2 * wp6.WEIGHTS[2], 3: 2 * wp6.WEIGHTS[3]}
    out = {}
    for name, (lam, w, mat) in dict(base=(wp6.LAMBDA, wp6.WEIGHTS, M),
                                    lam20=(20.0, wp6.WEIGHTS, M), lam150=(150.0, wp6.WEIGHTS, M),
                                    small_only=(wp6.LAMBDA, w0, M),
                                    large_x2=(wp6.LAMBDA, w2, M),
                                    zero_crossing=(wp6.LAMBDA, wp6.WEIGHTS, zc)).items():
        rate, P = wp6.daily(mat, lam=lam, weights=w)
        out[name] = dict(per_day=rate, P_day=P)
    return out


def contradictions(events):
    """Held by a controller but infeasible by the necessary bound: would contradict the bound."""
    bad = []
    for e in events:
        if "error" in e or e["outer_feasible"]:
            continue
        held_conv = not e["conv"]["violation"]
        c = e.get("cert")
        held_cert = c is not None and c["certified"] and not c["violation"]
        if held_conv or held_cert:
            bad.append(e)
    return bad


def fig_events(events, D, path):
    cols = [("D1", "unavoidable"), ("D1", "conv"), ("D3", "unavoidable"), ("D3", "conv"),
            ("D2", "unavoidable"), ("D2", "conv"), ("D2", "cert"),
            ("D4", "unavoidable"), ("D4", "conv"), ("D4", "cert")]
    cols = [c for c in cols if c[0] in D]
    head = {"unavoidable": "bound", "conv": "conv.", "cert": "cert."}
    types = _types()
    cmap = LinearSegmentedColormap.from_list("seq_blue", (SURFACE,) + SEQ)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2), constrained_layout=True)
    for ax, rating in zip(axes, RATINGS):
        A = np.full((len(types), len(cols)), np.nan)
        for c, (d, kind) in enumerate(cols):
            M = p_matrix(events, d, rating, kind)
            for r, (i, j) in enumerate(types):
                A[r, c] = 100 * M[i][j]
        ax.imshow(A, cmap=cmap, vmin=0, vmax=100, aspect="auto")
        for r in range(A.shape[0]):
            for c in range(A.shape[1]):
                v = A[r, c]
                if np.isfinite(v):
                    ax.text(c, r, f"{v:.0f}", ha="center", va="center", fontsize=7.5,
                            color="white" if v > 55 else INK)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([f"{d}\n{head[k]}" for d, k in cols], fontsize=8, color=INK)
        ax.set_yticks(range(len(types)))
        ax.set_yticklabels([_evtxt(i, j) for i, j in types], fontsize=8, color=INK)
        ax.set_ylabel("AC2 load before → after the event (fraction of rating)", fontsize=9, color=INK2)
        ax.set_title(NAMES[rating], fontsize=11, color=INK)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xticks(np.arange(-0.5, len(cols)), minor=True)
        ax.set_yticks(np.arange(-0.5, len(types)), minor=True)
        ax.grid(which="minor", color=SURFACE, linewidth=2)
        ax.tick_params(which="minor", length=0)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 100))
    cb = fig.colorbar(sm, ax=axes, shrink=0.6, pad=0.01)
    cb.set_label("Events that fail at a uniformly random phase (%)", fontsize=9, color=INK2)
    cb.outline.set_visible(False)
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)


def main():
    events = json.loads((RES / "wp6_events.json").read_text())
    d4 = json.loads((RES / "wp6_d4.json").read_text()) if (RES / "wp6_d4.json").exists() else {}
    D = wp6.designs()
    errs = [e for e in events if "error" in e]
    S = dict(designs={k: dict(name=v["name"], cb=v["cb"]) for k, v in D.items()}, daily={},
             criteria={}, errors=len(errs))
    for d in D:
        for r in RATINGS:
            if r not in D[d]["cb"]:
                continue
            for kind in ("unavoidable", "conv", "cert"):
                t = daily_table(events, d, r, kind)
                if t is not None:
                    S["daily"][f"{d}/{r}/{kind}"] = t
    def pday(d, r, kind, case="base"):
        t = S["daily"].get(f"{d}/{r}/{kind}")
        return None if t is None else 100 * t[case]["P_day"]

    def share_ok(d, r):
        v = [outcome_fail(e, "cert") for e in events if e["design"] == d and e["rating"] == r]
        v = [x for x in v if x is not None]
        return None if not v else 100 * (1 - np.mean(v)), len(v)

    crit = S["criteria"]
    crit["W1"] = {r: dict(P_day_pct=pday("D1", r, "unavoidable"),
                          pass_=(pday("D1", r, "unavoidable") or 0) >= W1_PCT) for r in RATINGS}
    crit["W2"] = {r: dict(infeasible_events=sum(1 for e in events if e["design"] == "D2"
                                                and e["rating"] == r and outcome_fail(e, "unavoidable")),
                          ) for r in RATINGS}
    for r in RATINGS:
        crit["W2"][r]["pass_"] = crit["W2"][r]["infeasible_events"] == 0
    for key, d in (("W3", "D2"), ("W3b", "D4")):
        if d in D:
            crit[key] = {}
            for r in RATINGS:
                pct, n = share_ok(d, r)
                crit[key][r] = dict(share_pct=pct, n=n, pass_=pct is not None and pct >= W3_PCT)
    S["contradictions"] = len(contradictions(events))
    (RES / "wp6_summary.json").write_text(json.dumps(S, indent=1, default=float))
    fig_events(events, D, RES / "fig16_wp6_events.png")
    write_md(events, D, d4, S, pday, share_ok, errs)
    print(json.dumps(crit, indent=1, default=float))


def _pct(x):
    return "—" if x is None or not np.isfinite(x) else f"{x:.1f} %"


def write_md(events, D, d4, S, pday, share_ok, errs):
    L = ["# WP6: system case study\n",
         "Generated by `sim/wp6_report.py` from `results/wp6/wp6_*.json` (written by "
         "`sim/run_wp6.py`). The protocol, with its revision before the campaign, is "
         "`sim/WP6_PROTOCOL.md`.\n"]
    L.append("## Designs\n")
    L.append("| Design | | 120 V | 230 V | Operating point |\n|---|---|---|---|---|")
    for k, v in D.items():
        L.append(f"| {k} | {v['name']} | " + " | ".join(
            f"{v['cb'][r] * 1e6:.0f} µF" if r in v["cb"] else "—" for r in RATINGS)
            + f" | {'centred' if v['schedule'] is None else 'load schedule'} |")
    L.append("")
    if d4:
        L.append("D4 search (full load drop certified at the 24 grid phases):\n")
        for r in RATINGS:
            if r in d4:
                L.append(f"- {NAMES[r]}: " + "; ".join(
                    f"×{t['factor']:g} ({t['cb'] * 1e6:.0f} µF): {sum(t['certified'])}/24"
                    for t in d4[r]["tried"]) + ".")
        L.append("")
    L.append("## Failure probability per event type\n")
    L.append("Share of the 24 grid phases at which the event fails, in %. *No controller*: the "
             "circuit outer bound is infeasible. *Conventional*: the conventional controller trips "
             "in the switched model. *Certified*: not certified, or the tracked certified plan "
             "trips (D2 and D4). See `fig16_wp6_events.png`.\n")
    cols = [(d, k) for d in D for k in (("unavoidable", "conv", "cert") if d in ("D2", "D4")
                                       else ("unavoidable", "conv"))]
    head = {"unavoidable": "no controller", "conv": "conventional", "cert": "certified"}
    for r in RATINGS:
        L.append(f"**{NAMES[r]}**\n")
        L.append("| Event | " + " | ".join(f"{d} {head[k]}" for d, k in cols) + " |")
        L.append("|---|" + "---|" * len(cols))
        mats = {c: p_matrix(events, c[0], r, c[1]) for c in cols}
        for i, j in _types():
            L.append(f"| {_evtxt(i, j)} | " + " | ".join(
                _pct(100 * mats[c][i][j]) for c in cols) + " |")
        L.append("")
    L.append("## Daily probability (stylised load model)\n")
    pi, P = wp6.transition_model()
    L.append(f"λ = {wp6.LAMBDA:.0f} events per day; stationary load distribution "
             + ", ".join(f"{LEVEL_TXT[i]}: {pi[i]:.2f}" for i in range(len(pi)))
             + ". Probability of at least one failing event in a day:\n")
    L.append("| Design | Outcome | 120 V | 230 V |\n|---|---|---|---|")
    for d in D:
        for kind in ("unavoidable", "conv", "cert"):
            if any(f"{d}/{r}/{kind}" in S["daily"] for r in RATINGS):
                L.append(f"| {d} | {head[kind]} | " + " | ".join(
                    _pct(pday(d, r, kind)) for r in RATINGS) + " |")
    L.append("")
    L.append("Sensitivity of the daily probability (120 V / 230 V):\n")
    cases = dict(lam20="λ = 20 per day", lam150="λ = 150 per day",
                 small_only="no jumps larger than 1/3", large_x2="large-jump weights doubled",
                 zero_crossing="zero-crossing switching (0° and 180°)")
    L.append("| Case | " + " | ".join(f"{d} {head[k]}" for d, k in
                                       (("D1", "unavoidable"), ("D1", "conv"), ("D2", "cert"),
                                        ("D4", "cert")) if f"{d}/r120/{k}" in S["daily"]) + " |")
    sel = [(d, k) for d, k in (("D1", "unavoidable"), ("D1", "conv"), ("D2", "cert"), ("D4", "cert"))
           if f"{d}/r120/{k}" in S["daily"]]
    L.append("|---|" + "---|" * len(sel))
    for c, txt in cases.items():
        L.append(f"| {txt} | " + " | ".join(
            f"{_pct(pday(d, 'r120', k, c))} / {_pct(pday(d, 'r230', k, c))}" for d, k in sel) + " |")
    L.append("")
    L.append("## Criteria\n")
    cr = S["criteria"]
    L.append("| ID | Criterion | 120 V | 230 V | Verdict |\n|---|---|---|---|---|")
    L.append(f"| W1 | Unavoidable violations at conventional sizing, P_day ≥ {W1_PCT:.0f} % | "
             + " | ".join(_pct(cr["W1"][r]["P_day_pct"]) for r in RATINGS)
             + f" | {'pass' if all(cr['W1'][r]['pass_'] for r in RATINGS) else 'fail'} |")
    L.append("| W2 | Proposed sizing (D2): no event infeasible by the necessary bound | "
             + " | ".join(f"{cr['W2'][r]['infeasible_events']} events" for r in RATINGS)
             + f" | {'pass' if all(cr['W2'][r]['pass_'] for r in RATINGS) else 'fail'} |")
    for key, d in (("W3", "D2"), ("W3b", "D4")):
        if key in cr:
            L.append(f"| {key} | {d} with the certified controller: certified and held ≥ {W3_PCT:.0f} % | "
                     + " | ".join(f"{_pct(cr[key][r]['share_pct'])} of {cr[key][r]['n']}" for r in RATINGS)
                     + f" | {'pass' if all(cr[key][r]['pass_'] for r in RATINGS) else 'fail'} |")
    L.append(f"| W4 | Conventional-controller trips per day, D1 against D2 | "
             + " | ".join(f"{_pct(pday('D1', r, 'conv'))} / {_pct(pday('D2', r, 'conv'))}" for r in RATINGS)
             + " | reported |")
    L.append(f"| W5 | Operating point alone: D3 against D1, no-controller / conventional | "
             + " | ".join(f"{_pct(pday('D3', r, 'unavoidable'))} vs {_pct(pday('D1', r, 'unavoidable'))}; "
                          f"{_pct(pday('D3', r, 'conv'))} vs {_pct(pday('D1', r, 'conv'))}" for r in RATINGS)
             + " | reported |")
    L.append("")
    L.append("**W6, cost of the proposed sizing.** Extra capacitance against D1, and the extra "
             "energy it stores at the bus maximum voltage (the buffer's voltage rating):\n")
    L.append("| Design | 120 V | 230 V |\n|---|---|---|")
    for d in ("D2", "D4"):
        if d in D:
            cells = []
            for r in RATINGS:
                if r not in D[d]["cb"]:
                    cells.append("—")
                    continue
                dc = D[d]["cb"][r] - D["D1"]["cb"][r]
                Vmax = design(r)[0].Vdc_max
                cells.append(f"+{dc * 1e6:.0f} µF (+{100 * dc / D['D1']['cb'][r]:.0f} %), "
                             f"+{0.5 * dc * Vmax ** 2:.0f} J at {Vmax:.0f} V")
            L.append(f"| {d} | " + " | ".join(cells) + " |")
    L.append("")
    L.append(f"**Consistency check.** Events held by a controller but infeasible by the necessary "
             f"bound: {S['contradictions']}. Events with errors: {len(errs)}.\n")
    (RES / "wp6_case_study.md").write_text("\n".join(L), encoding="utf-8")


if __name__ == "__main__":
    main()
