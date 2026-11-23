"""WP9 report: digital-twin validation against an emulated plant (WP9_TWIN_PROTOCOL.md).

    py sim/wp9_report.py

Reads results/<rating>/wp9_twin.json with the WP4 truth (wp4_mc.json, wp4_mccircuit.json) and
writes results/wp9/wp9_summary.json and results/wp9/wp9_twin.md.
"""
import json
from pathlib import Path

import numpy as np

RES = Path(__file__).resolve().parent / "results"
OUT = RES / "wp9"
RATINGS = ("r120", "r230")
RATED = {"r120": 2400.0, "r230": 4800.0}
TW1_LIMIT = dict(Cb=0.05, Lb=0.025, E1=0.025, E2=0.025, slew_dc1=0.075)
PARAM_NAMES = dict(Cb="Cb", Lb="Lb", E1="E1 (grid amplitude)", E2="E2 (backup amplitude)",
                   slew_dc1="DC1 slew", Rb="Rb (reported)", Cdc="Cdc (reported)")


def _load(rating, name):
    return json.loads((RES / rating / name).read_text())


def truth_labels(rating):
    """Table 3 truth: infeasible only if the plant's reduced and circuit outer bounds both are."""
    M, MCC = _load(rating, "wp4_mc.json"), _load(rating, "wp4_mccircuit.json")
    cc = {r["index"]: r for r in MCC["cases"]}
    out = {}
    for c in M["cases"]:
        if "error" in c:
            continue
        t = c["truth"]
        if t == "infeasible" and cc[c["index"]]["plant_circuit_outer"]:
            t = "unresolved"
        out[c["index"]] = dict(truth=t, N=dict(circuit_outer=cc[c["index"]]["nominal_circuit_outer"],
                                               inner=c["admits"]["nominal_inner"]))
    return out


def _q(x, q):
    return float(np.quantile(np.asarray(x, float), q)) if len(x) else float("nan")


def classification(cases, truth, get):
    """False-feasible and false-infeasible counts and rates of a verdict function."""
    ff = fi = n_inf = n_feas = 0
    for c in cases:
        t = truth[c["index"]]["truth"]
        v = get(c)
        if t == "infeasible":
            n_inf += 1
            ff += bool(v)
        elif t == "feasible":
            n_feas += 1
            fi += (not v)
    return dict(false_feasible=ff, false_infeasible=fi, n_infeasible=n_inf, n_feasible=n_feas,
                ff_rate=ff / n_inf if n_inf else np.nan, fi_rate=fi / n_feas if n_feas else np.nan,
                total=ff + fi)


def analyse(rating):
    R = _load(rating, "wp9_twin.json")
    truth = truth_labels(rating)
    cases = [c for c in R["cases"] if "error" not in c and c["index"] in truth]
    S = dict(n_cases=len(R["cases"]), n_errors=sum("error" in c for c in R["cases"]),
             n_scored=len(cases),
             counts={k: sum(truth[c["index"]]["truth"] == k for c in cases)
                     for k in ("feasible", "infeasible", "unresolved")})
    # M1 identification
    m1 = {}
    for k in list(TW1_LIMIT) + ["Rb", "Cdc"]:
        src = "param_error" if k in TW1_LIMIT else "reported_error"
        e = np.array([c[src][k] for c in cases if np.isfinite(c[src].get(k, np.nan))])
        a = np.abs(e)
        m1[k] = dict(n=len(e), median_abs=_q(a, 0.5), p95_abs=_q(a, 0.95), max_abs=float(a.max()),
                     mean=float(e.mean()), limit=TW1_LIMIT.get(k),
                     passed=bool(_q(a, 0.95) <= TW1_LIMIT[k]) if k in TW1_LIMIT else None)
    S["M1"] = m1
    S["clamped"] = sum(bool(c["calibration"]["clamped"]) for c in cases)
    S["calibration_runs_failed"] = sum(not c["calibration"]["calibration_runs_held"] for c in cases)
    S["TW1"] = all(v["passed"] for k, v in m1.items() if k in TW1_LIMIT)
    # M2 boundary
    m2 = {}
    for name in ("N", "NE", "T"):
        e = np.array([(c["max_step_W"][name] - c["max_step_W"]["true"]) / RATED[rating] for c in cases])
        m2[name] = dict(mae=float(np.mean(np.abs(e))), p95_abs=_q(np.abs(e), 0.95),
                        bias=float(np.mean(e)), share_optimistic=float(np.mean(e > 0.0075)),
                        max_abs=float(np.max(np.abs(e))))
    S["M2"] = m2
    S["TW2"] = bool(m2["T"]["mae"] <= 0.5 * m2["NE"]["mae"])
    # M3 classification
    m3 = {}
    for bound in ("circuit_outer", "inner"):
        m3[bound] = dict(
            N=classification(cases, truth, lambda c: truth[c["index"]]["N"][bound]),
            NE=classification(cases, truth, lambda c: c["verdicts"]["NE"][bound]),
            T=classification(cases, truth, lambda c: c["verdicts"]["T"][bound]))
    S["M3"] = m3
    S["TW3"] = all(m3[b]["T"]["false_feasible"] <= m3[b]["N"]["false_feasible"]
                   and m3[b]["T"]["false_infeasible"] <= m3[b]["N"]["false_infeasible"]
                   and m3[b]["T"]["total"] < m3[b]["N"]["total"] for b in m3)
    # M4 closed loop
    m4 = {}
    for name in ("NE", "T"):
        held = viol = rej = 0
        adm_res = viol_res = 0
        for c in cases:
            t = truth[c["index"]]["truth"]
            r = c["closed_loop"][name]
            if t == "feasible":
                if r is None:
                    rej += 1
                elif r["held"]:
                    held += 1
                else:
                    viol += 1
            if r is not None and t != "unresolved":
                adm_res += 1
                viol_res += not r["held"]
        nf = S["counts"]["feasible"]
        m4[name] = dict(feasible_held=held, feasible_violated=viol, feasible_rejected=rej,
                        held_share_of_feasible=held / nf if nf else np.nan,
                        admitted_resolved=adm_res, admitted_violated=viol_res,
                        violated_share_of_admitted=viol_res / adm_res if adm_res else np.nan)
    S["M4"] = m4
    S["TW4"] = bool(m4["T"]["held_share_of_feasible"] >= 0.90
                    and m4["T"]["violated_share_of_admitted"] <= 0.05)
    # M5 drift
    d = [c["drift"] for c in cases]
    S["M5a"] = dict(alarm_aged=float(np.mean([x["alarm_aged"] for x in d])),
                    alarm_unaged=float(np.mean([x["alarm_unaged"] for x in d])),
                    residual_unaged_median_V=_q([x["residual_unaged_V"] for x in d], 0.5),
                    residual_aged_median_V=_q([x["residual_aged_V"] for x in d], 0.5),
                    threshold_median_V=_q([x["threshold_V"] for x in d], 0.5),
                    aged_runs_failed=sum(not x["aged_runs_held"] for x in d))
    S["TW5"] = bool(S["M5a"]["alarm_aged"] >= 0.95 and S["M5a"]["alarm_unaged"] <= 0.05)
    aged_truth = {c["index"]: dict(truth=c["drift"]["truth"]) for c in cases}
    m5b = dict(counts={k: sum(v["truth"] == k for v in aged_truth.values())
                       for k in ("feasible", "infeasible", "unresolved")})
    for bound in ("circuit_outer", "inner"):
        m5b[bound] = {name: classification(cases, aged_truth,
                                           lambda c, n=name: c["drift"]["verdicts"][n][bound])
                      for name in ("N", "stale", "recal")}
    S["M5b"] = m5b
    S["M5b_recal_param_p95"] = {k: _q(np.abs([c["drift"]["recal_param_error"][k] for c in cases]), 0.95)
                                for k in TW1_LIMIT}
    # M6 computation
    S["M6"] = dict(t_identify_median_s=_q([c["t_identify_s"] for c in cases], 0.5),
                   t_verdicts_median_s=_q([c["t_twin_verdicts_s"] for c in cases], 0.5),
                   t_case_median_s=_q([c["t_case_s"] for c in cases], 0.5))
    return S


def _pct(x, d=1):
    return "—" if x is None or not np.isfinite(x) else f"{100 * x:.{d}f} %"


def write_md(A):
    L = ["# WP9: digital-twin validation against an emulated plant\n",
         "Protocol: `sim/WP9_TWIN_PROTOCOL.md` (criteria fixed before any run; revision 1 after the "
         "method test, before the campaign). Plant: the switched model with hidden parameters from "
         "the WP4 tolerance set; the twin sees 12-bit, noisy measurements with gain errors. "
         "Evaluation events: the 300 WP4 misclassification cases per rating, unseen by the twin. "
         "**No result here is evidence about hardware.**\n"]
    L.append("## Verdicts\n")
    L.append("| Criterion | 120 V | 230 V |\n|---|---|---|")
    for k, name in (("TW1", "TW1 identification"), ("TW2", "TW2 boundary"),
                    ("TW3", "TW3 classification (H3)"), ("TW4", "TW4 closed loop"),
                    ("TW5", "TW5 drift monitor")):
        L.append(f"| {name} | {'pass' if A['r120'][k] else '**fail**'} | "
                 f"{'pass' if A['r230'][k] else '**fail**'} |")
    for r in RATINGS:
        S = A[r]
        L.append(f"\n## {r}\n")
        L.append(f"Cases {S['n_cases']}, errors {S['n_errors']}, scored {S['n_scored']}; truth "
                 f"{S['counts']}. Clamped identifications: {S['clamped']}; calibration runs that "
                 f"violated a limit: {S['calibration_runs_failed']}.\n")
        L.append("**M1 identification** (relative error)\n")
        L.append("| Parameter | median abs | 95th pct abs | max abs | mean | TW1 limit |\n|---|---|---|---|---|---|")
        for k, v in S["M1"].items():
            L.append(f"| {PARAM_NAMES[k]} | {_pct(v['median_abs'], 2)} | {_pct(v['p95_abs'], 2)} | "
                     f"{_pct(v['max_abs'], 2)} | {_pct(v['mean'], 2)} | "
                     f"{_pct(v['limit'], 1) if v['limit'] else 'reported'} |")
        L.append("\n**M2 boundary error** (predicted minus true largest step, / AC2 rating)\n")
        L.append("| Predictor | mean abs | 95th pct abs | bias | optimistic share | max abs |\n|---|---|---|---|---|---|")
        for n, v in S["M2"].items():
            L.append(f"| {n} | {_pct(v['mae'], 2)} | {_pct(v['p95_abs'], 2)} | {_pct(v['bias'], 2)} | "
                     f"{_pct(v['share_optimistic'], 0)} | {_pct(v['max_abs'], 2)} |")
        L.append("\n**M3 classification** (false feasible / false infeasible)\n")
        L.append("| Bound | N | NE | T |\n|---|---|---|---|")
        for b, row in S["M3"].items():
            cells = [f"{x['false_feasible']} ({_pct(x['ff_rate'])}) / {x['false_infeasible']} "
                     f"({_pct(x['fi_rate'])})" for x in (row["N"], row["NE"], row["T"])]
            L.append(f"| {b} | " + " | ".join(cells) + " |")
        L.append("\n**M4 closed loop** (truth-feasible cases; admitted resolved cases)\n")
        L.append("| Predictor | held | violated | rejected | held share | admitted violated |\n|---|---|---|---|---|---|")
        for n, v in S["M4"].items():
            L.append(f"| {n} | {v['feasible_held']} | {v['feasible_violated']} | {v['feasible_rejected']} | "
                     f"{_pct(v['held_share_of_feasible'])} | {v['admitted_violated']} of "
                     f"{v['admitted_resolved']} ({_pct(v['violated_share_of_admitted'])}) |")
        a = S["M5a"]
        L.append(f"\n**M5 drift** (aged capacitor: −15 % Cb, 2× Rb): alarms on aged plants "
                 f"{_pct(a['alarm_aged'])}, on unaged {_pct(a['alarm_unaged'])}; median residual "
                 f"{a['residual_unaged_median_V']:.2f} V unaged, {a['residual_aged_median_V']:.2f} V "
                 f"aged, threshold {a['threshold_median_V']:.2f} V. Aged truth {S['M5b']['counts']}.\n")
        L.append("| Bound | N | stale twin | recalibrated twin |\n|---|---|---|---|")
        for b in ("circuit_outer", "inner"):
            row = S["M5b"][b]
            cells = [f"{row[n]['false_feasible']} / {row[n]['false_infeasible']}"
                     for n in ("N", "stale", "recal")]
            L.append(f"| {b} | " + " | ".join(cells) + " |")
        m = S["M6"]
        L.append(f"\n**M6 computation** (one core): identification {m['t_identify_median_s']:.2f} s, "
                 f"twin verdicts (outer + inner) {m['t_verdicts_median_s']:.1f} s, whole case "
                 f"{m['t_case_median_s']:.0f} s (medians).\n")
    (OUT / "wp9_twin.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main():
    OUT.mkdir(exist_ok=True)
    A = {r: analyse(r) for r in RATINGS if (RES / r / "wp9_twin.json").exists()}
    (OUT / "wp9_summary.json").write_text(json.dumps(A, indent=1, default=float))
    if len(A) == len(RATINGS):
        write_md(A)
    print(json.dumps({r: {k: A[r][k] for k in ("TW1", "TW2", "TW3", "TW4", "TW5")} for r in A},
                     indent=1))


if __name__ == "__main__":
    main()
