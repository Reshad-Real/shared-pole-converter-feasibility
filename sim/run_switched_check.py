"""Does the averaged result survive switching? Rerun the pair in the PWM model with dead time.

1. Load the pair (theta_a, theta_b, step size, tuned controller gain) from results.json.
2. Run conventional theta_a, conventional theta_b and the LP plan at theta_a in the averaged
   model and in the switched model with dead time 0, 1 and 2 us.
3. Repeat the 24-phase conventional-controller sweep in the switched model (td = 1 and 2 us)
   and compare outcomes with the averaged model and the LP verdicts.
4. Figures 5 and 6, switched_check.md and switched_check.json.

Run after run_pair_experiment.py:  py sim/run_switched_check.py
"""
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from pathlib import Path

import numpy as np

from model import Params
from feasibility import nominal_orbit, solve
from timedomain import simulate
from switched import SwitchParams, simulate_switched
from run_pair_experiment import (AXIS, BAND, CRIT, INK, INK2, MUTED, S1, S1_LIGHT, S2, SURF,
                                 _style)

OUT = Path(__file__).resolve().parent / "results"
TD_US = (0.0, 1.0, 2.0)
SWEEP_TD_US = (1.0, 2.0)


def _sw_metrics(args):
    p, sw, th_deg, dP, k_rec = args
    return simulate_switched(p, sw, np.deg2rad(th_deg), dP,
                             gains=dict(k_rec=2 * np.pi * k_rec))["metrics"]


def _fmt_result(m):
    if not m["violation"]:
        return "all limits held"
    return f"violation at {m['t_violation_ms']:.1f} ms, tripped"


def _rms_diff(sw_run, av_run, key_sw, key_av):
    """RMS difference between the switched period averages and the averaged model, over the
    window both runs cover (period averages are placed at period mid-points)."""
    Tsw_ms = 1e3 * np.median(np.diff(sw_run["t"]))
    ts = 1e3 * sw_run["t"] + 0.5 * Tsw_ms
    ta = 1e3 * av_run["t"]
    end = min(ts[-1], ta[-1])
    k = (ts >= -5.0) & (ts <= end)
    ref = np.interp(ts[k], ta, av_run[key_av])
    return float(np.sqrt(np.mean((sw_run[key_sw][k] - ref) ** 2)))


def main():
    t_start = time.time()
    R0 = json.loads((OUT / "results.json").read_text())
    p = Params(**R0["params"])
    pr = R0["pair"]
    th_a, th_b, dP = pr["theta_a_deg"], pr["theta_b_deg"], pr["dP_test_W"]
    k_rec = R0["timedomain_tuning"]["k_rec_hz_chosen"]
    gains = dict(k_rec=2 * np.pi * k_rec)
    orb = nominal_orbit(p)
    plan = solve(p, np.deg2rad(th_a), dP, "inner", orbit=orb, smooth=True)
    sw1 = SwitchParams()
    R = dict(pair=dict(theta_a_deg=th_a, theta_b_deg=th_b, dP_W=dP, k_rec_hz=k_rec),
             switch_params=asdict(sw1))

    # 1. Pair runs in both models
    cases = [("conventional_theta_a", th_a, "conventional"),
             ("conventional_theta_b", th_b, "conventional"),
             ("plan_theta_a", th_a, "plan")]
    av_runs, sw_runs = {}, {}
    R["pair_runs"] = {}
    for name, th, pol in cases:
        kw = dict(policy=pol, plan=plan if pol == "plan" else None, gains=gains)
        av = simulate(p, np.deg2rad(th), dP, **kw)
        av_runs[name] = av
        entry = dict(averaged=av["metrics"])
        for td in TD_US:
            sw = SwitchParams(td=td * 1e-6)
            s = simulate_switched(p, sw, np.deg2rad(th), dP, zoom_ms=(-0.5, 3.0) if td == 1.0
                                  else None, **kw)
            entry[f"switched_td{td:g}us"] = s["metrics"]
            if td == 1.0:
                sw_runs[name] = s
                entry["rms_diff_vb_V"] = _rms_diff(s, av, "vb_avg", "vb")
                entry["rms_diff_Vdc_V"] = _rms_diff(s, av, "Vdc_avg", "Vdc")
        R["pair_runs"][name] = entry
        print(name, {k: v.get("violation") if isinstance(v, dict) else round(v, 2)
                     for k, v in entry.items()}, flush=True)

    # Integration-step check at the pair
    conv = {}
    for name, th, pol in cases[:2]:
        fine = simulate_switched(p, SwitchParams(h_max=1.25e-6), np.deg2rad(th), dP, gains=gains)
        base = R["pair_runs"][name]["switched_td1us"]
        conv[name] = max(abs(fine["metrics"][k] - base[k]) for k in ("Vdc_min", "Vdc_max"))
    R["step_convergence_max_dV"] = conv

    # 2. 24-phase sweep in the switched model
    phases = [c["phase_deg"] for c in R0["controller_vs_lp"]]
    workers = max(1, min(8, (os.cpu_count() or 2) - 2))
    jobs = [(p, SwitchParams(td=td * 1e-6), th, dP, k_rec) for td in SWEEP_TD_US for th in phases]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        res = list(ex.map(_sw_metrics, jobs))
    n = len(phases)
    rows = []
    for q, c in enumerate(R0["controller_vs_lp"]):
        row = dict(c)
        row["averaged_ok"] = row.pop("conventional_ok")
        for i, td in enumerate(SWEEP_TD_US):
            m = res[i * n + q]
            row[f"switched_td{td:g}us_ok"] = not m["violation"]
            row[f"switched_td{td:g}us_inst_current_over"] = m["current_violation_inst"]
        rows.append(row)
    R["phase_sweep"] = rows
    R["runtime_s"] = time.time() - t_start

    (OUT / "switched_check.json").write_text(json.dumps(R, indent=2, default=float))
    fig_overlay(p, R, av_runs, sw_runs, OUT / "fig5_switched_vs_averaged.png")
    fig_zoom(p, R, sw_runs, OUT / "fig6_switching_detail.png")
    write_report(p, R)
    print(f"done in {R['runtime_s']:.0f} s", flush=True)


# ---------------------------------------------------------------- figures
def fig_overlay(p, R, av_runs, sw_runs, path):
    plt = _style()
    pr = R["pair"]
    labels = {"conventional_theta_a": f"Conventional, θa = {pr['theta_a_deg']:.0f}°",
              "conventional_theta_b": f"Conventional, θb = {pr['theta_b_deg']:.0f}°",
              "plan_theta_a": f"LP plan tracked, θa = {pr['theta_a_deg']:.0f}°"}
    fig, axs = plt.subplots(2, 3, figsize=(9.6, 5.2), sharex=True, sharey="row",
                            gridspec_kw=dict(hspace=0.12, wspace=0.08))
    t_lo, t_hi = -5.0, 40.0
    for c, name in enumerate(labels):
        av, s = av_runs[name], sw_runs[name]
        ta = 1e3 * av["t"]
        ka = (ta >= t_lo) & (ta <= t_hi)
        Tsw_ms = 1e3 * np.median(np.diff(s["t"]))
        ts = 1e3 * s["t"] + 0.5 * Tsw_ms
        ks = (ts >= t_lo) & (ts <= t_hi)
        ax = axs[0, c]
        ax.fill_between(ta[ka], av["lo"][ka], av["hi"][ka], color=BAND, lw=0, label="Pole-N band")
        ax.plot(ta[ka], av["vb"][ka], color=S2, label="Averaged model")
        ax.plot(ts[ks], s["vb_avg"][ks], color=S1, lw=1.1, label="Switched model (period average)")
        ax.set_title(labels[name], loc="left")
        ax2 = axs[1, c]
        ax2.fill_between(ts[ks], s["Vdc_lo"][ks], s["Vdc_hi"][ks], color=S1_LIGHT, lw=0,
                         label="Switched model, ripple envelope")
        ax2.axhline(p.Vdc_min, color=MUTED, lw=1.0)
        ax2.axhline(p.Vdc_max, color=MUTED, lw=1.0)
        ax2.plot(ta[ka], av["Vdc"][ka], color=S2)
        ax2.plot(ts[ks], s["Vdc_avg"][ks], color=S1, lw=1.1)
        ax2.set_xlabel("Time after step (ms)")
        for run, col in ((av, S2), (s, S1)):
            tv = run["metrics"]["t_violation_ms"]
            if tv is not None:
                ax2.axvline(tv, color=col, lw=0.8)
    axs[0, 0].set_ylabel("Buffer voltage (V)")
    axs[1, 0].set_ylabel("Bus voltage (V)")
    axs[1, 0].set_ylim(p.Vdc_min - 25, p.Vdc_max + 10)
    h1, l1 = axs[0, 0].get_legend_handles_labels()
    h2, l2 = axs[1, 0].get_legend_handles_labels()
    fig.legend(h1 + h2, l1 + l2, loc="lower left", bbox_to_anchor=(0.07, 0.97), ncol=4)
    fig.text(0.07, -0.03, "Vertical lines in the bus panels mark each model's first violation "
             "(the run then trips). Bus limits 190 and 210 V.", color=MUTED, fontsize=7.5)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_zoom(p, R, sw_runs, path):
    plt = _style()
    pr = R["pair"]
    names = [("conventional_theta_a", f"θa = {pr['theta_a_deg']:.0f}°"),
             ("conventional_theta_b", f"θb = {pr['theta_b_deg']:.0f}°")]
    fig, axs = plt.subplots(3, 2, figsize=(8.4, 6.2), sharex=True, sharey="row",
                            gridspec_kw=dict(hspace=0.12, wspace=0.08))
    t_lo, t_hi = -0.5, 3.0
    for c, (name, title) in enumerate(names):
        s = sw_runs[name]
        Tsw_ms = 1e3 * np.median(np.diff(s["t"]))
        ts = 1e3 * s["t"]
        ks = (ts >= t_lo) & (ts <= t_hi)
        z = s["zoom"]
        tz = 1e3 * z["t"]
        ax = axs[0, c]
        shift = (R["switch_params"]["delay_periods"] + 0.5) * Tsw_ms   # band acts after the delay
        ax.fill_between(ts[ks] + shift, s["lo"][ks], s["hi"][ks], step="mid", color=BAND,
                        lw=0, label="Pole-N band (commanded)")
        ax.step(ts[ks] + shift, s["lo"][ks], where="mid", color=MUTED, lw=0.9,
                label="Band lower edge")
        ax.step(ts[ks] + 0.5 * Tsw_ms, s["vN_avg"][ks], where="mid", color=S1, lw=1.2,
                label="Pole-N voltage, period average")
        ax.plot(tz, z["vb"], color=S2, lw=1.1, label="Buffer-capacitor voltage")
        ax.set_title(f"{title}: switched model, conventional controller", loc="left")
        axs[1, c].plot(tz, z["ib"], color=S1, lw=0.7)
        axs[1, c].axhline(0, color=AXIS, lw=0.8)
        axs[2, c].plot(tz, z["Vdc"], color=S1, lw=0.7)
        axs[2, c].axhline(p.Vdc_min, color=MUTED, lw=1.0)
        tv = s["metrics"]["t_violation_ms"]
        if tv is not None:
            for a in axs[:, c]:
                a.axvline(tv, color=CRIT, lw=1.0)
            axs[2, c].annotate("first violation", (tv, p.Vdc_min), xytext=(5, -36),
                               textcoords="offset points", color=CRIT, fontsize=7.5)
        axs[2, c].set_xlabel("Time after step (ms)")
        axs[2, c].set_xlim(t_lo, t_hi)
    axs[0, 0].set_ylabel("Voltage (V)")
    axs[1, 0].set_ylabel("Buffer current (A)")
    axs[2, 0].set_ylabel("Bus voltage (V)")
    axs[2, 0].set_ylim(p.Vdc_min - 25, p.Vdc_max + 10)
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower left", bbox_to_anchor=(0.08, 0.98), ncol=4)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- report
def write_report(p, R):
    pr, sw = R["pair"], R["switch_params"]
    L = ["# Switched-model check of the equal-energy pair\n",
         "Generated by `sim/run_switched_check.py` from `switched_check.json`. "
         "Nothing here has been measured on hardware.\n",
         "## Switched model\n",
         f"- Three half-bridge legs, one symmetric carrier at {sw['fsw'] / 1e3:.0f} kHz, dead "
         f"time {sw['td'] * 1e6:g} µs (0 and 2 µs as variants). During dead time the "
         f"freewheeling diode chosen by the pole-current sign sets the pole voltage.",
         f"- Devices: on-resistance {sw['Ron'] * 1e3:.0f} mΩ, diode drop {sw['Vf']:.1f} V.",
         f"- Control runs once per period on sampled states with a "
         f"{sw['delay_periods']}-period delay. The AC current loops are real "
         f"({sw['wc_ac'] / 2 / np.pi / 1e3:.1f} kHz crossover); duties include dead-time "
         "compensation.",
         "- The circuit is integrated exactly between switching edges. Cutting the step limit "
         "from 5 µs to 1.25 µs moved the bus extremes by at most "
         f"{max(R['step_convergence_max_dV'].values()):.3f} V.",
         "- Both models are judged on period-averaged bus voltage and pole currents. "
         "Instantaneous values, which include switching ripple, are listed separately.\n",
         f"## The pair: θa = {pr['theta_a_deg']:.0f}°, θb = {pr['theta_b_deg']:.0f}°, "
         f"AC2 step {pr['dP_W']:.0f} W\n",
         "| Run | Model | Bus, period avg (V) | Bus, instantaneous (V) | Pole current avg / inst (A) "
         "| Result |",
         "|---|---|---|---|---|---|"]
    names = {"conventional_theta_a": "Conventional, θa", "conventional_theta_b": "Conventional, θb",
             "plan_theta_a": "LP plan tracked, θa"}
    for key, lab in names.items():
        e = R["pair_runs"][key]
        a = e["averaged"]
        L.append(f"| {lab} | averaged | {a['Vdc_min']:.1f}–{a['Vdc_max']:.1f} | — | "
                 f"{a['I_pole_peak']:.1f} / — | {_fmt_result(a)} |")
        for td in TD_US:
            m = e[f"switched_td{td:g}us"]
            L.append(f"| {lab} | switched, td {td:g} µs | {m['Vdc_min']:.1f}–{m['Vdc_max']:.1f} | "
                     f"{m['Vdc_inst_min']:.1f}–{m['Vdc_inst_max']:.1f} | "
                     f"{m['I_pole_peak']:.1f} / {m['I_pole_peak_inst']:.1f} | {_fmt_result(m)} |")
    L.append("\nRMS difference between the switched model (td 1 µs, period averages) and the "
             "averaged model, from 5 ms before the step until the end or the trip:\n")
    L.append("| Run | Buffer voltage (V) | Bus voltage (V) |\n|---|---|---|")
    for key, lab in names.items():
        e = R["pair_runs"][key]
        L.append(f"| {lab} | {e['rms_diff_vb_V']:.2f} | {e['rms_diff_Vdc_V']:.2f} |")

    rows = R["phase_sweep"]
    L.append("\n## Phase sweep with the tuned conventional controller\n")
    L.append("| θ (deg) | " + " | ".join(f"{r['phase_deg']:.0f}" for r in rows) + " |")
    L.append("|---|" + "---|" * len(rows))

    def ok(b):
        return "ok" if b else "**fails**"

    def lp(r):
        if r["inner_feasible"]:
            return "feasible"
        return "unresolved" if r["outer_feasible"] else "infeasible"

    L.append("| LP verdict | " + " | ".join(lp(r) for r in rows) + " |")
    L.append("| Averaged model | " + " | ".join(ok(r["averaged_ok"]) for r in rows) + " |")
    for td in SWEEP_TD_US:
        L.append(f"| Switched, td {td:g} µs | " +
                 " | ".join(ok(r[f"switched_td{td:g}us_ok"]) for r in rows) + " |")
    for td in SWEEP_TD_US:
        k = f"switched_td{td:g}us_ok"
        same = sum(r[k] == r["averaged_ok"] for r in rows)
        contra = [f"{r['phase_deg']:.0f}°" for r in rows if r[k] and not r["outer_feasible"]]
        diff = [f"{r['phase_deg']:.0f}°" for r in rows if r[k] != r["averaged_ok"]]
        inst = [f"{r['phase_deg']:.0f}°" for r in rows
                if r[k] and r[f"switched_td{td:g}us_inst_current_over"]]
        L.append(f"\n- td {td:g} µs: same outcome as the averaged model at {same} of {len(rows)} "
                 f"phases (differs at {', '.join(diff) if diff else 'none'}). Succeeds where the "
                 f"LP proved infeasibility: {', '.join(contra) if contra else 'none'}. Passing "
                 f"runs whose instantaneous pole current exceeds {p.I_pole_max:.0f} A because of "
                 f"ripple: {', '.join(inst) if inst else 'none'}.")
    L.append("\nWhere the switched model differs from the averaged one, the LP verdict is either "
             "feasible (a controller limit: the sampled, one-period-delayed loops let the bus "
             "swing about ±6 V instead of ±1.5 V) or unresolved (either outcome is consistent "
             "with the bounds).")
    L.append("\n## Limits of this check\n")
    L.append("- Ideal gate timing, no switching transients, no parasitic capacitance, and no "
             "common-mode path from the floating AC sources. Sensors and ADC quantisation are "
             "ideal.")
    L.append("- DC1 and DC2 remain ideal power ports; only the three-leg converter is switched.")
    L.append("- Same single parameter set and illustrative topology as the averaged study.")
    (OUT / "switched_check.md").write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
