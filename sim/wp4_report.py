"""WP4 report: robustness criteria R1-R6 at both realistic ratings.

Reads results/<rating>/wp4_{bounds,events,mc,nonideal}.json written by run_wp4.py and writes
results/wp4_robustness.md, results/wp4_robustness.json and figures 11-13.

Run:  py sim/wp4_report.py
"""
import json
from pathlib import Path

import numpy as np

from run_pair_experiment import INK, INK2, MUTED, S1, S1_LIGHT, S2, SURF, _style

RES = Path(__file__).resolve().parent / "results"
RATINGS = ("r120", "r230")
NAMES = {"r120": "120 V rating", "r230": "230 V rating"}
CRITERIA = (("energy_sizing", "Energy sizing + steady-state check"),
            ("instantaneous", "Instantaneous checks"),
            ("static_derating", "Static derating"),
            ("nominal_inner", "Inner bound, nominal values"),
            ("nominal_outer", "Outer bound, nominal values"),
            ("robust_inner", "Inner bound, robust (all 16 corners)"),
            ("robust_outer", "Outer bound, robust (any of 24 plants)"))
RUNS = (("conventional_theta_a", "Conv. θa"), ("plan_theta_a", "Plan θa"),
        ("conventional_theta_b", "Conv. θb"))


def _load(rating, part):
    f = RES / rating / f"wp4_{part}.json"
    return json.loads(f.read_text()) if f.exists() else None


def _pf(b):
    return "**pass**" if b else "**fail**"


# ---------------------------------------------------------------- analysis
def misclassification(M, criteria=CRITERIA):
    cases = [c for c in M["cases"] if "error" not in c]
    truth = [c["truth"] for c in cases]
    n = {k: truth.count(k) for k in ("feasible", "infeasible", "unresolved")}
    rows = {}
    for key, _ in criteria:
        adm = [c["admits"][key] for c in cases]
        ff = sum(a and t == "infeasible" for a, t in zip(adm, truth))
        fi = sum((not a) and t == "feasible" for a, t in zip(adm, truth))
        adm_res = sum(a and t != "unresolved" for a, t in zip(adm, truth))
        rows[key] = dict(false_feasible=ff, false_infeasible=fi,
                         ff_rate=ff / n["infeasible"] if n["infeasible"] else np.nan,
                         fi_rate=fi / n["feasible"] if n["feasible"] else np.nan,
                         admitted_resolved=adm_res,
                         infeasible_share_of_admitted=ff / adm_res if adm_res else np.nan)
    errors = len(M["cases"]) - len(cases)
    return dict(counts=n, errors=errors, rows=rows)


def misclassification_circuit(M, MCC):
    """R5 rescored: infeasible only when the plant's circuit outer LP is infeasible too; cases
    the reduced LP excluded but the circuit LP does not become unresolved."""
    cc = {r["index"]: r for r in MCC["cases"]}
    cases = []
    for c in M["cases"]:
        if "error" in c:
            continue
        x = dict(c)
        x["admits"] = dict(c["admits"], nominal_circuit_outer=cc[c["index"]]["nominal_circuit_outer"])
        if c["truth"] == "infeasible" and cc[c["index"]]["plant_circuit_outer"]:
            x["truth"] = "unresolved"
        cases.append(x)
    out = misclassification(dict(cases=cases), CRITERIA + (("nominal_circuit_outer", ""),))
    out["errors"] = len(M["cases"]) - len(cases)
    return out


def contradictions(M, N):
    """R4: runs that held every limit although the outer LP from their own state says no."""
    bad, total = [], 0
    if M:
        for c in M["cases"]:
            if "outer_from_run_state" in c:
                total += 1
                if not c["outer_from_run_state"]:
                    bad.append(f"MC case {c['index']}")
    if N:
        for r in N["runs"]:
            m = r.get("metrics")
            if m is None or m["violation"] or "outer_from_run_state" not in r:
                continue
            total += 1
            if not r["outer_from_run_state"]:
                bad.append(f"{r['condition']} / {r['run']}")
    return dict(successful_runs_checked=total, contradictions=bad)


def nonideal_table(N):
    by = {}
    for r in N["runs"]:
        by.setdefault(r["condition"], {})[r["run"]] = r
    rows = []
    for c in N["conditions"]:
        rr = by.get(c["key"], {})
        pre = any(r.get("metrics", {}).get("pre_step_violation") for r in rr.values())
        b = rr.get("conventional_theta_b", {})
        b_trips = bool(b.get("metrics", {}).get("violation"))
        b_lp_inf = b.get("outer_from_orbit") is False
        a_ok = any(rr.get(k, {}).get("metrics", {}).get("violation") is False
                   for k in ("conventional_theta_a", "plan_theta_a"))
        a_lp = rr.get("conventional_theta_a", {}).get("outer_from_orbit")
        rows.append(dict(key=c["key"], label=c["label"], group=c["group"], runs=rr,
                         pre_step_violation=pre, theta_b_trips=b_trips,
                         theta_b_outer_infeasible=b_lp_inf,
                         R6=None if pre else (b_trips or b_lp_inf),
                         theta_a_held=a_ok, theta_a_outer_feasible=a_lp,
                         theta_a_failure=None if (a_ok or pre) else
                         ("physical" if a_lp is False else "controller"),
                         effect_visible=a_ok and b_trips and not pre))
    return rows


def _short(row):
    """Condition label without the rating-specific corner factors."""
    return row["label"].split(" (")[0] if row["group"] == "tolerance" else row["label"]


def _run_word(r):
    if r is None:
        return "—"
    if r.get("no_plan"):
        return "no plan"
    m = r["metrics"]
    if not m["violation"]:
        return "held"
    return f"trip {m['t_violation_ms']:.1f} ms"


# ---------------------------------------------------------------- figures
def fig_bounds(B, path):
    plt = _style()
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.9), sharex=True,
                            gridspec_kw=dict(wspace=0.18))
    for ax, rating in zip(axs, RATINGS):
        b = B[rating]
        ph = np.array(b["phases_deg"])
        outer = np.array(b["outer_W"])
        ax.fill_between(ph, outer.min(axis=0), outer.max(axis=0), color=S1_LIGHT, lw=0,
                        label="Outer bound across the 24 plants")
        ax.plot(ph, b["conventional_W"], color=S2, label="Phase-blind checks (nominal values)")
        ax.plot(ph, b["nominal_outer_W"], color=S1, lw=1.0, ls=(0, (4, 2)),
                label="Outer bound, nominal plant")
        ax.plot(ph, b["robust_outer_W"], color=S1, label="Robust outer bound (largest)")
        ax.plot(ph, b["robust_inner_W"], color=INK2, lw=1.2,
                label="Robust inner bound (smallest of 16 corners)")
        r2 = b["R2"]
        ax.plot([r2["phase_deg"]], [r2["robust_outer_W"]], "o", ms=6, color=S1, mec=SURF,
                mew=1.5, zorder=5)
        ax.annotate(f"worst phase {r2['phase_deg']:.0f}°: phase-blind\nlimit {r2['over_pct']:.0f} % "
                    "above the robust\nouter bound", (r2["phase_deg"], r2["robust_outer_W"]),
                    xytext=(r2["phase_deg"] - 20, 0.78 * b["cap_W"]), ha="right", va="center",
                    color=INK, fontsize=7.5,
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=2, shrinkB=4))
        ax.set_title(NAMES[rating], loc="left")
        ax.set_xlim(0, 351)
        ax.set_xticks(np.arange(0, 361, 90))
        ax.set_ylim(0, b["cap_W"] * 1.02)
        ax.set_xlabel("AC1 phase at which the step arrives (deg)")
    axs[0].set_ylabel("Largest admissible AC2 step (W)")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower left", bbox_to_anchor=(0.07, 0.98), ncol=3, columnspacing=1.2)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_misclassification(MC, path):
    plt = _style()
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.8), sharey=True, gridspec_kw=dict(wspace=0.08))
    y = np.arange(len(CRITERIA))
    h = 0.36
    for ax, (key, title) in zip(axs, (("ff_rate", "False feasible: admitted, but provably infeasible"),
                                      ("fi_rate", "False infeasible: rejected, but shown feasible"))):
        for j, (rating, col) in enumerate(zip(RATINGS, (S1, S2))):
            v = np.array([100 * MC[rating]["rows"][c][key] for c, _ in CRITERIA])
            ax.barh(y + (j - 0.5) * h, v, height=h - 0.04, color=col, label=NAMES[rating])
            for yi, vi in zip(y + (j - 0.5) * h, v):
                ax.text(vi + 1.0, yi, f"{vi:.0f}", va="center", fontsize=7.5, color=INK2)
        ax.set_title(title, loc="left", fontsize=9.5)
        ax.set_axisbelow(True)
        ax.set_xlim(0, 105)
        ax.set_xlabel("Share of cases (%)")
        ax.grid(axis="y", visible=False)
    axs[0].set_yticks(y)
    axs[0].set_yticklabels([lab for _, lab in CRITERIA])
    axs[0].invert_yaxis()
    h_, l_ = axs[0].get_legend_handles_labels()
    fig.legend(h_, l_, loc="lower left", bbox_to_anchor=(0.30, 0.98), ncol=2)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_nonideal(NT, lim, path):
    plt = _style()
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 7.6), sharey=True, gridspec_kw=dict(wspace=0.06))
    marks = {"conventional_theta_a": ("o", S1, -0.2), "plan_theta_a": ("s", INK2, 0.0),
             "conventional_theta_b": ("D", S2, 0.2)}
    for ax, rating in zip(axs, RATINGS):
        Vn = 0.5 * (lim[rating]["Vdc_min"] + lim[rating]["Vdc_max"])
        for i, row in enumerate(NT[rating]):
            if row["pre_step_violation"]:
                ax.text(-11.8, i, "bus out of limits before the step", fontsize=7, color=MUTED,
                        va="center")
                continue
            for run, lab in RUNS:
                r = row["runs"].get(run)
                if r is None or "metrics" not in r:
                    continue
                m = r["metrics"]
                margin = min(m["Vdc_min"] - lim[rating]["Vdc_min"],
                             lim[rating]["Vdc_max"] - m["Vdc_max"])
                x = max(100 * margin / Vn, -12.0)
                mk, col, dy = marks[run]
                ax.plot([x], [i + dy], mk, ms=6, color=col if not m["violation"] else SURF,
                        mec=col, mew=1.5, label=lab if i == 0 else None, zorder=3)
        ax.axvline(0, color=INK2, lw=1.0)
        ax.set_xlim(-12.5, 4)
        ax.set_title(NAMES[rating], loc="left")
        ax.set_xlabel("Smallest bus-voltage margin after the step (% of nominal)")
    axs[0].set_yticks(range(len(NT[RATINGS[0]])))
    axs[0].set_yticklabels([_short(r) for r in NT[RATINGS[0]]], fontsize=7.5)
    axs[0].invert_yaxis()
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower left", bbox_to_anchor=(0.30, 0.985), ncol=3)
    fig.text(0.30, 0.975, "Filled: every limit held.  Hollow: a limit was violated and the run "
             "tripped.  Margins below −12 % are drawn at −12 %.", fontsize=7.5, color=MUTED,
             va="top")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- report
def main():
    B = {r: _load(r, "bounds") for r in RATINGS}
    E = {r: _load(r, "events") for r in RATINGS}
    M = {r: _load(r, "mc") for r in RATINGS}
    N = {r: _load(r, "nonideal") for r in RATINGS}
    lim = {}
    for r in RATINGS:
        W = json.loads((RES / r / "results.json").read_text())
        lim[r] = dict(Vdc_min=W["params"]["Vdc_min"], Vdc_max=W["params"]["Vdc_max"])
    MC = {r: misclassification(M[r]) for r in RATINGS if M[r]}
    NT = {r: nonideal_table(N[r]) for r in RATINGS if N[r]}
    R4 = {r: contradictions(M[r], N[r]) for r in RATINGS}
    D = {r: _load(r, "diagnostics") for r in RATINGS}
    F40 = {r: _load(r, "deadtime40") for r in RATINGS}
    F40 = F40 if all(F40.values()) else None
    CD = {r: _load(r, "conditional") for r in RATINGS}
    C4 = {r: _load(r, "r4circuit") for r in RATINGS}
    CC = {r: _load(r, "circuit") for r in RATINGS}
    CK = {r: _load(r, "circuitcheck") for r in RATINGS}
    MCC = {r: _load(r, "mccircuit") for r in RATINGS}
    MCc = ({r: misclassification_circuit(M[r], MCC[r]) for r in RATINGS}
           if all(MCC.values()) else {})

    S = {}
    for r in RATINGS:
        s = {}
        if B[r]:
            s.update(R1=B[r]["R1"]["pass"], R2=B[r]["R2"]["pass"], R3=B[r]["R3"]["pass"])
        s["R4"] = not R4[r]["contradictions"]
        if r in NT:
            vals = [row["R6"] for row in NT[r] if row["R6"] is not None]
            s["R6"] = all(vals)
            s["R6_failed"] = [row["key"] for row in NT[r] if row["R6"] is False]
            s["R6_outside_envelope"] = [row["key"] for row in NT[r] if row["R6"] is None]
        S[r] = s
    (RES / "wp4_robustness.json").write_text(json.dumps(
        dict(summary=S, misclassification=MC, R4=R4,
             nonideal={r: [{k: v for k, v in row.items() if k != "runs"} for row in NT[r]]
                       for r in NT}), indent=1, default=float))

    if all(B.values()):
        fig_bounds(B, RES / "fig11_wp4_robust_bounds.png")
    if len(MC) == 2:
        fig_misclassification(MC, RES / "fig12_wp4_misclassification.png")
    if len(NT) == 2:
        fig_nonideal(NT, lim, RES / "fig13_wp4_nonideal.png")
    if all(C4.values()):
        for r in RATINGS:
            S[r]["R4_circuit_bound"] = not C4[r]["contradictions_circuit"]
    write_md(B, E, M, N, MC, NT, R4, S, D, F40, CD, lim, (C4, CC, CK, MCc))
    print(json.dumps(S, indent=1))


def diagnostics_md(B, D, F40):
    """Post-hoc diagnostics: they explain the verdicts and change none of them."""
    L = ["## Diagnostics after the verdicts\n",
         "These were run after the verdicts were known, to explain them. They change no verdict.\n"]
    if all(D.values()):
        L.append("**Does every plant show the phase effect?** Here the phase-blind limits (energy "
                 "sizing plus the steady-state check) are recomputed for each corner plant itself "
                 "and compared with that plant's own outer bound, at that plant's worst phase.\n")
        L.append("| | 120 V rating | 230 V rating |\n|---|---|---|")
        pp = {r: sorted(D[r]["per_plant"], key=lambda x: x["over_pct"]) for r in RATINGS}
        L.append("| Worst-phase over-admission, smallest over the 16 corners | " + " | ".join(
            f"{pp[r][0]['over_pct']:.0f} % ({pp[r][0]['plant']})" for r in RATINGS) + " |")
        L.append("| Median over the 16 corners | " + " | ".join(
            f"{np.median([x['over_pct'] for x in pp[r]]):.0f} %" for r in RATINGS) + " |")
        L.append("| Largest over the 16 corners | " + " | ".join(
            f"{pp[r][-1]['over_pct']:.0f} % ({pp[r][-1]['plant']})" for r in RATINGS) + " |\n")
        L.append("**Which tolerance moves the bound?** Outer bound at the R2 phase with one "
                 "parameter at a time at the edge of its range.\n")
        L.append("| Change | 120 V rating | 230 V rating |\n|---|---|---|")
        L.append("| None (nominal plant) | " + " | ".join(
            f"{D[r]['oat_nominal_outer_W']:.0f} W at {D[r]['oat_phase_deg']:.0f}°" for r in RATINGS) + " |")
        for i, row in enumerate(D["r120"]["oat"]):
            L.append(f"| {row['change']} | " + " | ".join(
                f"{D[r]['oat'][i]['delta_W']:+.0f} W" for r in RATINGS) + " |")
        L.append("\n**Why R4 failed.** For each run that held every limit although the outer LP "
                 "says no trajectory exists: the band relaxation the LP needs before it finds one, "
                 "and how far the run's buffer-capacitor voltage actually left the band.\n")
        L.append("| Rating | Run | Relaxation the LP needs | Buffer voltage outside the band | "
                 "Pole-N command clipped |\n|---|---|---|---|---|")
        for r in RATINGS:
            for x in D[r]["r4"]:
                L.append(f"| {NAMES[r]} | {x['condition']}, {x['run']} | {x['relaxation_needed_V']:.1f} V | "
                         f"{x['vb_outside_band_max_V']:.1f} V | {x['sat_time_ms']:.2f} ms |")
        L.append("")
    if F40:
        L.append("**40 kHz with a dead time sized for it.** With the nominal 1 µs dead time (4 % of "
                 "the 25 µs period), the dead-time compensation pushes near-full duties into the "
                 "clamp, the AC currents overshoot by 3–4 A and drain the bus. Rerun with 0.25 µs:\n")
        L.append("| Rating | Conv. θa | Plan θa | Conv. θb |\n|---|---|---|---|")
        for r in RATINGS:
            L.append(f"| {NAMES[r]} | " + " | ".join(_run_word(x) for x in F40[r]["runs"]) + " |")
        L.append("")
    return L


def conditional_md(CD):
    """Post hoc: R1-R3 with the grid amplitude measured."""
    L = ["**† With the grid amplitude measured.** A converter measures its grid voltage, so here "
         "the amplitude is fixed at each level and known to every check; only Cb (±20 %), Cdc "
         "(±20 %) and the DC1 response (±30 %) are unknown. Phase-blind checks use the measured "
         "amplitude and nominal components. Robust bounds range over the 8 corners of those "
         "three; the θb check adds 64 random draws. Pass marks as for R1–R3.\n",
         "| Rating | Grid amplitude | θb outer-feasible plants (R1) | Over-admission against the "
         "robust outer bound (R2) | Robust gap (R3) |\n|---|---|---|---|---|"]
    for r in RATINGS:
        if not CD.get(r):
            continue
        for key, lv in CD[r]["levels"].items():
            level = ("nominal" if abs(lv["E"] - 1) < 1e-9 else
                     f"{100 * (lv['E'] - 1):+.0f} %".replace("-", "−"))
            L.append(f"| {NAMES[r]} | {level} | {lv['R1']['n_outer_feasible']} of "
                     f"{lv['R1']['n_plants']} | {lv['R2']['over_pct']:.0f} % at {lv['R2']['phase_deg']:.0f}° "
                     f"({lv['R2']['over_W']:.0f} W) | {lv['R3']['gap_W']:.0f} W ({lv['R3']['gap_pct']:.0f} %) |")
    L.append("")
    return L


def circuit_md(C4, CC, CK, MCc, CD):
    """Follow-up to R4: the outer bound with the buffer inductor (sim/outer_circuit.py)."""
    L = ["## Follow-up: an outer bound that includes the buffer inductor\n",
         "Added after the verdicts, because of R4. `sim/outer_circuit.py` keeps the buffer branch "
         "in its own linear coordinates (capacitor voltage, inductor current, pole-N voltage) with "
         "the exact sampled LC dynamics, puts the band on the pole voltage instead of the "
         "capacitor voltage, and keeps the bus in energy. The stored energies are tied to the "
         "voltage and current by tangent cuts and chords, losses are bounded, and the pole-N "
         "current limit bounds the inductor current. It drops the terminal cycle, which a run "
         "that holds its limits over the horizon need not reach exactly.\n"]
    if all(CK.values()):
        L.append("**Consistency.** θb relaxation the circuit bound needs, against the reduced LP's:\n")
        L.append("| Rating | Reduced LP | Lb = 20 µH | Lb = 100 µH | Nominal Lb | Nominal Lb, twice the cuts |\n"
                 "|---|---|---|---|---|---|")
        for r in RATINGS:
            k = CK[r]
            L.append(f"| {NAMES[r]} | {k['reduced_relaxation_V']:.1f} V | " + " | ".join(
                f"{x['relaxation_V']:.1f} V" for x in k["by_Lb"]) + f" | {k['nominal_Lb_24_16_cuts_V']:.1f} V |")
        L.append("\nWith a small inductor it reproduces the reduced LP; the real inductor is what "
                 "lowers the margin.\n")
    if all(C4.values()):
        L.append("**R4 again.** Every WP4 run that held its limits was re-simulated (the runs "
                 "reproduce bit for bit) and checked against the circuit bound from its own state "
                 "at the step, including the inductor current.\n")
        L.append("| | 120 V rating | 230 V rating |\n|---|---|---|")
        L.append("| Runs that held every limit | " + " | ".join(str(C4[r]["n_held"]) for r in RATINGS) + " |")
        L.append("| Re-simulations identical to the stored runs | " + " | ".join(
            f"{C4[r]['n_reproduced']} of {C4[r]['n_checked']}" for r in RATINGS) + " |")
        L.append("| Contradictions with the reduced outer LP | " + " | ".join(
            str(len(C4[r]["contradictions_reduced"])) for r in RATINGS) + " |")
        L.append("| **Contradictions with the circuit outer LP** | " + " | ".join(
            f"**{len(C4[r]['contradictions_circuit'])}**" for r in RATINGS) + " |\n")
    if all(CC.values()):
        L.append("**Headline numbers with the circuit bound.**\n")
        L.append("| | 120 V rating | 230 V rating |\n|---|---|---|")
        L.append("| θb relaxation needed: reduced → circuit | " + " | ".join(
            f"{CC[r]['theta_b']['reduced_relaxation_V']:.1f} → {CC[r]['theta_b']['relaxation_V']:.1f} V"
            for r in RATINGS) + " |")
        L.append("| Nominal worst-phase over-admission (G1 metric): reduced → circuit | " + " | ".join(
            f"{_nominal_over(CC[r]):.0f} % → {CC[r]['worst']['over_pct']:.0f} % at "
            f"{CC[r]['worst']['phase_deg']:.0f}°" for r in RATINGS) + " |")
        for key, word in (("0.9", "−10 %"), ("1.0", "nominal"), ("1.1", "+10 %")):
            L.append(f"| Measured grid {word}: over-admission against the robust outer bound (R2), "
                     f"reduced → circuit | " + " | ".join(
                         f"{_signed(CD[r]['levels'][key]['R2']['over_pct'])} → "
                         f"{_signed(CC[r]['measured_grid'][key]['R2']['over_pct'])}" for r in RATINGS) + " |")
            L.append(f"| Measured grid {word}: θb outer-feasible plants (R1), reduced → circuit | " + " | ".join(
                f"{CD[r]['levels'][key]['R1']['n_outer_feasible']} → "
                f"{CC[r]['measured_grid'][key]['R1']['n_outer_feasible']} of "
                f"{CC[r]['measured_grid'][key]['R1']['n_plants']}" for r in RATINGS) + " |")
        L.append("\nThe measured-grid robust bounds with the circuit LP were searched in an 11-phase "
                 "window (±45°) around each level's worst phase from the reduced analysis.\n")
    if MCc:
        L.append("**R5 rescored with the circuit bound.** A case is infeasible only if its own "
                 "plant's circuit outer LP is infeasible too; cases only the reduced LP excluded "
                 "become unresolved. The nominal circuit outer bound is added as a criterion.\n")
        L.append("| | 120 V rating | 230 V rating |\n|---|---|---|")
        for k in ("feasible", "infeasible", "unresolved"):
            L.append(f"| Truth: {k} | " + " | ".join(str(MCc[r]["counts"][k]) for r in RATINGS) + " |")
        L.append("\n| Criterion | False feasible, 120 V | False infeasible, 120 V | False feasible, 230 V | "
                 "False infeasible, 230 V |\n|---|---|---|---|---|")
        for key, name in CRITERIA + (("nominal_circuit_outer", "Circuit outer bound, nominal values"),):
            a, b = MCc["r120"]["rows"][key], MCc["r230"]["rows"][key]

            def pct(x):
                return "—" if np.isnan(x) else f"{100 * x:.0f} %"
            L.append(f"| {name} | {pct(a['ff_rate'])} ({a['false_feasible']}) | {pct(a['fi_rate'])} "
                     f"({a['false_infeasible']}) | {pct(b['ff_rate'])} ({b['false_feasible']}) | "
                     f"{pct(b['fi_rate'])} ({b['false_infeasible']}) |")
        L.append("")
    return L


def _signed(x):
    """Percentage with a typographic minus, and no '-0'."""
    s = f"{x:.0f}"
    return ("0" if s in ("-0", "0") else s.replace("-", "−")) + " %"


def _nominal_over(cc):
    conv = np.array(cc["conventional_W"])
    red = np.array(cc["reduced_outer_W"])
    k = int(np.argmax(conv - red))
    return 100 * (conv[k] - red[k]) / conv[k]


def findings_md(B, D, MC, NT, R4, CD, lim, CIRC=None):
    """The reading of the verdicts, with every number taken from the result files."""
    L = ["## What the campaign shows\n",
         "Items marked † use diagnostics run after the verdicts were known (see the end of this "
         "report). They explain the verdicts and change none of them.\n"]
    pp = {r: [x["over_pct"] for x in D[r]["per_plant"]] for r in RATINGS} if all(D.values()) else None
    sp = {r: B[r]["R3"]["spread_outer_corners_W"] for r in RATINGS}
    L.append("1. **The effect is present in every plant, but the plant uncertainty is larger than "
             "the effect.**")
    if pp:
        L.append(f"   - † Computed plant by plant, worst-phase over-admission stays at "
                 f"{min(pp['r120']):.0f}–{max(pp['r120']):.0f} % (120 V) and "
                 f"{min(pp['r230']):.0f}–{max(pp['r230']):.0f} % (230 V) across all 16 corners.")
    L.append(f"   - A designer who knows only nominal values faces an outer bound of "
             f"{sp['r120'][0]:.0f}–{sp['r120'][1]:.0f} W (120 V) and {sp['r230'][0]:.0f}–"
             f"{sp['r230'][1]:.0f} W (230 V) at the worst phase. Against the most favourable plant "
             f"the phase-blind limit over-admits by only {B['r120']['R2']['over_pct']:.0f} % and "
             f"{B['r230']['R2']['over_pct']:.0f} % (R2), and θb becomes feasible in "
             f"{B['r120']['R1']['n_outer_feasible']} of {B['r120']['R1']['n_plants']} plants (R1).")
    if pp:
        L.append("   - † The grid amplitude moves the bound most, then Cb, then the DC1 response; Cdc, "
                 "the inductors and Rb hardly matter. Every corner plant that makes θb feasible "
                 "combines a 10 % low grid voltage with a 30 % faster DC1.")
    L.append(f"2. **Bounds made robust over the whole tolerance box are too conservative to use** "
             f"(R3). The robust inner bound rejects {100 * MC['r120']['rows']['robust_inner']['fi_rate']:.0f}–"
             f"{100 * MC['r230']['rows']['robust_inner']['fi_rate']:.0f} % of feasible cases, and the "
             f"robust outer bound admits {100 * MC['r120']['rows']['robust_outer']['ff_rate']:.0f}–"
             f"{100 * MC['r230']['rows']['robust_outer']['ff_rate']:.0f} % of infeasible ones, no better "
             f"than energy sizing.")
    rr = {r: MC[r]["rows"] for r in RATINGS}

    def rng(key, rate):
        a, b = sorted((100 * rr["r120"][key][rate], 100 * rr["r230"][key][rate]))
        return f"{a:.0f}–{b:.0f} %" if round(a) != round(b) else f"{a:.0f} %"
    L.append(f"3. **Under random tolerances the nominal trajectory bounds misclassify far less than "
             f"the phase-blind checks.** Nominal inner and outer bounds: {rng('nominal_inner', 'ff_rate')} "
             f"and {rng('nominal_outer', 'ff_rate')} false feasible, {rng('nominal_inner', 'fi_rate')} and "
             f"{rng('nominal_outer', 'fi_rate')} false infeasible. Energy sizing: {rng('energy_sizing', 'ff_rate')} "
             f"false feasible; instantaneous checks: {rng('instantaneous', 'ff_rate')}; static derating: "
             f"{rng('static_derating', 'fi_rate')} false infeasible.")
    n_ok = sum(R4[r]["successful_runs_checked"] for r in RATINGS)
    n_bad = sum(len(R4[r]["contradictions"]) for r in RATINGS)
    if all(D.values()):
        x = [y for r in RATINGS for y in D[r]["r4"]]
        L.append(f"4. **The infeasibility proof does not carry over to the circuit when its margin is "
                 f"small** (R4). † In {n_bad} of {n_ok} successful runs the circuit did what the "
                 f"reduced model forbids: its buffer-capacitor voltage left the pole-N band by "
                 f"{min(y['vb_outside_band_max_V'] for y in x):.0f}–{max(y['vb_outside_band_max_V'] for y in x):.0f} V, "
                 f"using the volt-seconds of the buffer inductor, while the LP needed only "
                 f"{min(y['relaxation_needed_V'] for y in x):.0f}–{max(y['relaxation_needed_V'] for y in x):.0f} V "
                 f"of relaxation. The nominal θb certificates need {D['r120']['theta_b_relaxation_nominal_V']:.1f} V "
                 f"and {D['r230']['theta_b_relaxation_nominal_V']:.1f} V.")
        C4, CC = (CIRC[0], CIRC[1]) if CIRC else ({}, {})
        if all(C4.values()) and all(CC.values()):
            nh = sum(C4[r]["n_held"] for r in RATINGS)
            nc = sum(len(C4[r]["contradictions_circuit"]) for r in RATINGS)
            L.append(f"   - **Follow-up: fixed.** An outer LP that includes the buffer inductor "
                     f"(last section) is contradicted by {nc} of the {nh} runs that held their limits. "
                     f"θb stays provably infeasible, with {CC['r120']['theta_b']['relaxation_V']:.1f} V "
                     f"and {CC['r230']['theta_b']['relaxation_V']:.1f} V of margin instead of "
                     f"{CC['r120']['theta_b']['reduced_relaxation_V']:.1f} V and "
                     f"{CC['r230']['theta_b']['reduced_relaxation_V']:.1f} V. The nominal worst-phase "
                     f"over-admission (the G1 metric) becomes {CC['r120']['worst']['over_pct']:.0f} % "
                     f"and {CC['r230']['worst']['over_pct']:.0f} % instead of "
                     f"{_nominal_over(CC['r120']):.0f} % and {_nominal_over(CC['r230']):.0f} %.")
        else:
            L.append("   - A proof for the circuit needs the inductor's volt-second budget in the "
                     "outer LP, or a stated minimum margin.")
    if NT:
        valid = {r: [row for row in NT[r] if not row["pre_step_violation"]] for r in RATINGS}
        trips = {r: sum(row["theta_b_trips"] for row in valid[r]) for r in RATINGS}
        vis = {r: sum(row["effect_visible"] for row in valid[r]) for r in RATINGS}
        phys = {r: [_short(row) for row in valid[r] if row["theta_a_failure"] == "physical"] for r in RATINGS}
        ctrl = {r: [_short(row) for row in valid[r] if row["theta_a_failure"] == "controller"] for r in RATINGS}
        ref_max = next(row for row in NT["r120"] if row["key"] == "reference")["runs"][
            "conventional_theta_a"]["metrics"]["Vdc_max"]
        L.append(f"5. **Grid, switching and sensor conditions.** θb tripped under {trips['r120']} of "
                 f"{len(valid['r120'])} valid conditions (120 V) and {trips['r230']} of "
                 f"{len(valid['r230'])} (230 V); it held under the 10 % sag, the +20° phase jump and the "
                 f"most favourable corner. The full split (θa held, θb tripped) appeared under "
                 f"{vis['r120']} and {vis['r230']} conditions.")
        L.append(f"   - θa failed for physical reasons (its own outer LP infeasible under the condition): "
                 f"{'; '.join(phys['r120']) or 'none'} (120 V); {'; '.join(phys['r230']) or 'none'} (230 V).")
        L.append(f"   - θa failed although the LP leaves room, so the controllers are the limit: "
                 f"{'; '.join(ctrl['r120']) or 'none'} (120 V); {'; '.join(ctrl['r230']) or 'none'} (230 V). "
                 f"In the nominal run the 120 V conventional controller keeps only "
                 f"{lim['r120']['Vdc_max'] - ref_max:.1f} V below the overvoltage limit, and the "
                 f"inner LP keeps no headroom on the bus limits, so plan tracking fails on small errors.")
    L.append("")
    if CD and all(CD.values()):
        L += conditional_md(CD)
        nom = {r: CD[r]["levels"]["1.0"] for r in RATINGS}
        sag = {r: CD[r]["levels"]["0.9"] for r in RATINGS}
        L.append("**What is and is not supported.**\n")
        L.append("- **Not supported:** H3 as predeclared. Over the full declared tolerance set the "
                 "phase effect is smaller than the designer's uncertainty about the plant. Bounds "
                 "made robust over that set are too wide to use: the inner one rejects about half "
                 "of the feasible steps, and the outer one admits as many infeasible steps as "
                 "energy sizing does.")
        CC = CIRC[1] if CIRC else {}
        MCc = CIRC[3] if CIRC else {}
        if all(CC.values()) and MCc:
            g = {r: {k: CC[r]["measured_grid"][k]["R2"]["over_pct"] for k in ("0.9", "1.0", "1.1")}
                 for r in RATINGS}

            def rc(key, rate):
                a, b = sorted((100 * MCc["r120"]["rows"][key][rate], 100 * MCc["r230"]["rows"][key][rate]))
                return f"{a:.0f}–{b:.0f} %" if round(a) != round(b) else f"{a:.0f} %"
            L.append(f"- **Supported, with the circuit bound:** the effect is a property of every plant "
                     f"in the set. With the grid amplitude measured it survives component and DC1 "
                     f"tolerances clearly at +10 % grid voltage ({g['r120']['1.1']:.0f} % and "
                     f"{g['r230']['1.1']:.0f} %), only narrowly at nominal voltage "
                     f"({g['r120']['1.0']:.0f} % and {g['r230']['1.0']:.0f} % against the 15 % mark), and "
                     f"not in a 10 % sag ({_signed(g['r120']['0.9'])} and {_signed(g['r230']['0.9'])}). Under "
                     f"random tolerances the nominal circuit outer bound admits "
                     f"{rc('nominal_circuit_outer', 'ff_rate')} of the infeasible cases and rejects "
                     f"{rc('nominal_circuit_outer', 'fi_rate')} of the feasible ones, against "
                     f"{rc('energy_sizing', 'ff_rate')} admitted by energy sizing.")
            L.append("- **Proof scope:** the outer bound now covers the buffer inductor and no run "
                     "contradicts it. It remains a proof for the averaged circuit sampled at Ts, with "
                     "the AC currents delivered as commanded.\n")
        else:
            L.append(f"- **Supported:** the effect is a property of every plant in the set; with the grid "
                     f"amplitude measured it survives component and DC1 tolerances at nominal and high "
                     f"grid voltage (over-admission {nom['r120']['R2']['over_pct']:.0f} % and "
                     f"{nom['r230']['R2']['over_pct']:.0f} % at nominal), but not in a 10 % sag "
                     f"({sag['r120']['R2']['over_pct']:.0f} % and {sag['r230']['R2']['over_pct']:.0f} %); "
                     f"and under random tolerances the nominal trajectory bounds admit "
                     f"{rng('nominal_inner', 'ff_rate')} (inner) and {rng('nominal_outer', 'ff_rate')} "
                     f"(outer) of the infeasible cases, against {rng('energy_sizing', 'ff_rate')} for "
                     f"energy sizing.")
            L.append("- **To fix before any proof claim:** the outer LP must allow for the buffer "
                     "inductor's volt-seconds (R4).\n")
    return L


def write_md(B, E, M, N, MC, NT, R4, S, D=None, F40=None, CD=None, lim=None, CIRC=None):
    L = ["# WP4: robustness of the phase effect\n",
         "Generated by `sim/wp4_report.py` from `results/<rating>/wp4_*.json` (written by "
         "`sim/run_wp4.py`). Criteria and pass marks were fixed in `sim/WP4_CRITERIA.md` before "
         "any campaign run. Everything is simulation; nothing is measured.\n",
         "## Verdicts\n",
         "| Criterion | 120 V rating | 230 V rating |\n|---|---|---|"]
    lab = dict(R1="R1: θb infeasible across the tolerance set",
               R2="R2: worst-phase over-admission against the robust outer bound ≥ 15 %",
               R3="R3: robust inner–outer gap ≤ 10 % of the conventional limit",
               R4="R4: no run succeeds where the outer LP proves infeasibility",
               R6="R6: θb still trips (or is LP-infeasible) under every condition")
    for k in ("R1", "R2", "R3", "R4", "R6"):
        L.append(f"| {lab[k]} | " + " | ".join(_pf(S[r][k]) if k in S[r] else "not run"
                                                for r in RATINGS) + " |")
    L.append("| R5: misclassification rates | reported | reported |")
    if "R4_circuit_bound" in S[RATINGS[0]]:
        L.append("| Follow-up, not predeclared: R4 with the circuit outer bound (buffer inductor "
                 "included) | " + " | ".join(_pf(S[r]["R4_circuit_bound"]) for r in RATINGS) + " |")
    L.append("")
    if all(B.values()) and len(MC) == 2 and D:
        L += findings_md(B, D, MC, NT, R4, CD, lim, CIRC)

    if all(B.values()):
        L.append("## Robust bounds (R1–R3)\n")
        L.append("| | 120 V rating | 230 V rating |\n|---|---|---|")

        def row(t, fn):
            L.append(f"| {t} | " + " | ".join(fn(B[r]) for r in RATINGS) + " |")
        row("WP1 pair and step", lambda b: f"θb = {b['R1']['theta_b_deg']:.0f}°, "
            f"{b['R1']['dP_W']:.0f} W")
        row("R1: plants where θb is outer-feasible", lambda b: f"{b['R1']['n_outer_feasible']} of "
            f"{b['R1']['n_plants']}")
        row("Worst phase against the robust outer bound", lambda b: f"{b['R2']['phase_deg']:.0f}°")
        row("Phase-blind limit (nominal)", lambda b: f"{b['R2']['conventional_W']:.0f} W")
        row("Outer bound, nominal plant", lambda b: f"{b['R2']['nominal_outer_W']:.0f} W")
        row("Robust outer bound", lambda b: f"{b['R2']['robust_outer_W']:.0f} W")
        row("**R2: over-admission against the robust outer bound**",
            lambda b: f"**{b['R2']['over_pct']:.0f} %** ({b['R2']['over_W']:.0f} W)")
        row("Over-admission against the nominal outer bound (WP1)",
            lambda b: f"{b['nominal_worst']['over_pct']:.0f} % at {b['nominal_worst']['phase_deg']:.0f}°")
        row("Robust inner bound at that phase", lambda b: f"{b['R3']['robust_inner_W']:.0f} W")
        row("**R3: robust gap**", lambda b: f"**{b['R3']['gap_W']:.0f} W "
            f"({b['R3']['gap_pct_of_conventional']:.0f} %)**")
        row("Outer bound over the 16 corners", lambda b: f"{b['R3']['spread_outer_corners_W'][0]:.0f}–"
            f"{b['R3']['spread_outer_corners_W'][1]:.0f} W")
        row("Plant-by-plant gap (largest / mean over corners)",
            lambda b: f"{b['R3']['plant_by_plant_gap_max_W']:.0f} / "
                      f"{b['R3']['plant_by_plant_gap_mean_W']:.0f} W")
        L.append("\n![Robust bounds against phase](fig11_wp4_robust_bounds.png)\n")

    if all(E.values()):
        L.append("## Other events (reported)\n")
        L.append("Worst-phase over-admission: the phase-blind limit minus the outer bound, at the "
                 "phase where that difference is largest.\n")
        L.append("| Event | Rating | Nominal plant | Robust outer bound | Bisection monotone |\n"
                 "|---|---|---|---|---|")
        for ev, name in (("step_down", "Step-down and regenerative (AC2 falls by up to 2 × P2_pre)"),
                         ("ramp10", "Step-up ramped over 10 ms")):
            for r in RATINGS:
                e = E[r][ev]
                wn, wr = e["worst_nominal"], e["worst_robust"]
                L.append(f"| {name} | {NAMES[r]} | {wn['over_pct']:.0f} % ({wn['over_W']:.0f} W at "
                         f"{wn['phase_deg']:.0f}°) | {wr['over_pct']:.0f} % ({wr['over_W']:.0f} W at "
                         f"{wr['phase_deg']:.0f}°) | {'yes' if e['monotone_spot_check'] else 'no'} |")
        L.append("")

    if len(MC) == 2:
        L.append("## Misclassification (R5)\n")
        L.append("Each case draws an AC1 phase, a step up to the search limit and a plant from the "
                 "tolerance set. Truth: *infeasible* if the plant's own outer LP is infeasible; "
                 "*feasible* if its inner-LP plan, tracked by a controller that knows the plant, "
                 "holds every limit in the switched model; otherwise *unresolved*.\n")
        L.append("| | 120 V rating | 230 V rating |\n|---|---|---|")
        for k in ("feasible", "infeasible", "unresolved"):
            L.append(f"| Truth: {k} | " + " | ".join(str(MC[r]["counts"][k]) for r in RATINGS) + " |")
        L.append("| Solver errors (excluded) | " + " | ".join(str(MC[r]["errors"]) for r in RATINGS) + " |\n")
        L.append("| Criterion | False feasible, 120 V | False infeasible, 120 V | False feasible, 230 V | "
                 "False infeasible, 230 V | Admitted cases that are infeasible, 120 V / 230 V |\n"
                 "|---|---|---|---|---|---|")
        for key, name in CRITERIA:
            a, b = MC["r120"]["rows"][key], MC["r230"]["rows"][key]

            def pct(x):
                return "—" if np.isnan(x) else f"{100 * x:.0f} %"
            L.append(f"| {name} | {pct(a['ff_rate'])} ({a['false_feasible']}) | {pct(a['fi_rate'])} "
                     f"({a['false_infeasible']}) | {pct(b['ff_rate'])} ({b['false_feasible']}) | "
                     f"{pct(b['fi_rate'])} ({b['false_infeasible']}) | "
                     f"{pct(a['infeasible_share_of_admitted'])} / {pct(b['infeasible_share_of_admitted'])} |")
        L.append("\nFalse-feasible rate: share of the truly infeasible cases the criterion admits. "
                 "False-infeasible rate: share of the truly feasible cases it rejects. Counts in "
                 "brackets.\n")
        L.append("![Misclassification rates](fig12_wp4_misclassification.png)\n")

    if len(NT) == 2:
        L.append("## Non-idealities (R6)\n")
        L.append("Switched-model runs of the WP1 pair (θa = 270°, θb = 90°; 1000 W at the 120 V "
                 "rating, 1975 W at the 230 V rating). *LP at θb* is the outer LP of the actual "
                 "plant under the condition.\n")
        L.append("| Condition | 120 V: conv. θa / plan θa / conv. θb | LP at θb | 230 V: conv. θa / "
                 "plan θa / conv. θb | LP at θb |\n|---|---|---|---|---|")
        by230 = {row["key"]: row for row in NT["r230"]}
        for row in NT["r120"]:
            cells = []
            for rr in (row, by230[row["key"]]):
                if rr["pre_step_violation"]:
                    cells += ["outside the design envelope (bus out of limits before the step)", "—"]
                    continue
                cells.append(" / ".join(_run_word(rr["runs"].get(k)) for k, _ in RUNS))
                b = rr["runs"].get("conventional_theta_b", {})
                cells.append("infeasible" if b.get("outer_from_orbit") is False else
                             ("not excluded" if b.get("outer_from_orbit") else "—"))
            L.append(f"| {_short(row)} | " + " | ".join(cells) + " |")
        for r in RATINGS:
            corners = [row["label"].split(" (", 1)[1].rstrip(")") for row in NT[r]
                       if row["group"] == "tolerance"]
            L.append(f"\nTolerance corners at the {NAMES[r]} (least / most favourable outer bound "
                     f"at θa): {' / '.join(corners)}.")
        L.append("\n![Non-ideality runs](fig13_wp4_nonideal.png)\n")

    L.append("## No contradictions (R4)\n")
    for r in RATINGS:
        c = R4[r]
        L.append(f"- {NAMES[r]}: {c['successful_runs_checked']} runs held every limit; in "
                 f"{len(c['contradictions'])} of them the outer LP from the run's own state at "
                 f"the step says no trajectory exists"
                 + (f" ({', '.join(c['contradictions'])})." if c["contradictions"] else "."))
    L.append("")
    if D:
        L += diagnostics_md(B, D, F40)
    if CIRC and any(any(x.values()) for x in CIRC):
        L += circuit_md(*CIRC, CD=CD)
    (RES / "wp4_robustness.md").write_text("\n".join(L), encoding="utf-8")


if __name__ == "__main__":
    main()
