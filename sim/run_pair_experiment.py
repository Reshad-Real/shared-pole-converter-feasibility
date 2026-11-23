"""First discriminating experiment (plan section 11): equal stored energy, different AC phase.

1. Steady-state design facts for the chosen parameters.
2. Sweep the phase at which an AC2 load step arrives; find the largest admissible step under
   four criteria (trajectory test inner/outer, energy sizing, static derating).
3. Pick the pair (theta_a, theta_a + 180 deg) that best separates the trajectory test from
   the phase-blind checks. Both members have the same buffer energy, buffer current, bus
   state, DC1 power and port powers; only the AC voltage and current polarity differ.
4. At the pair: verdicts of every criterion, the band relaxation theta_b would need, and the
   buffer capacitance each phase would need.
5. Independent time-domain averaged runs (buffer inductor, losses, nonlinear bus).
6. Figures, results.json and summary.md in sim/results/.

Run:  py sim/run_pair_experiment.py
"""
import argparse
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np

from model import Params, ac_state, band, bisect, orbit_energy_at, steady_orbit
from feasibility import feasible, nominal_orbit, pole_currents, post_step_limit, solve
from timedomain import GAINS, simulate
import ratings

OUT = Path(__file__).resolve().parent / "results"
PHASE_STEP_DEG = 9
TOL_W = 1.0
K_REC_GRID_HZ = (5, 10, 20, 40, 80)


# ---------------------------------------------------------------- parallel workers
def _max_step(args):
    p, th_deg, mode, cap = args[:4]
    tol = args[4] if len(args) > 4 else TOL_W
    orb = nominal_orbit(p)
    return bisect(lambda d: feasible(p, np.deg2rad(th_deg), d, mode, orbit=orb), 0.0, cap, tol)


def _feasible(args):
    p, th_deg, dP, mode, delta = args
    return feasible(p, np.deg2rad(th_deg), dP, mode, delta=delta)


def _min_cb(args):
    """Smallest Cb (F) for which `mode` is feasible at (theta, dP); C2 adds the steady check."""
    p, th_deg, dP, mode = args

    def infeasible(cb):
        q = replace(p, Cb=cb)
        if mode == "conventional":
            ok = (feasible(q, np.deg2rad(th_deg), dP, "energy")
                  and steady_orbit(q, q.P2_pre + dP, q.Vdc_nom, 0.0) is not None)
        else:
            ok = feasible(q, np.deg2rad(th_deg), dP, mode)
        return not ok

    return bisect(infeasible, 20e-6, 3e-3, 2e-6)


def _sim_metrics(args):
    p, th_deg, dP, k_rec_hz = args
    s = simulate(p, np.deg2rad(th_deg), dP, gains=dict(k_rec=2 * np.pi * k_rec_hz))
    return s["metrics"]


# ---------------------------------------------------------------- checks
def initial_state(p, theta0, orbit):
    t0 = theta0 / p.w
    Eb0 = orbit_energy_at(orbit, p, t0)
    vb0 = np.sqrt(2 * Eb0 / p.Cb)
    pb0 = float(np.interp(np.mod(t0, p.T), orbit["t"], orbit["pb"], period=p.T))
    u1, u2, i1, i2, pac = ac_state(p, t0)
    return dict(Eb0_J=Eb0, vb0_V=vb0, ib0_A=pb0 / vb0, pb0_W=pb0, Vdc0_V=p.Vdc_nom,
                Pdc1_0_W=orbit["Pd"], pac0_W=float(pac), u1_V=float(u1), u2_V=float(u2),
                i1_A=float(i1), i2_A=float(i2))


def instantaneous_check(p, theta0, dP, orbit):
    """Conventional pointwise checks: AC commands synthesisable along the horizon at nominal
    bus, initial buffer voltage inside the band now, and pre/post steady orbits exist."""
    t_step = theta0 / p.w
    t = t_step + p.Ts * np.arange(p.n_cycles * p.Np + 1)
    u1, u2, _, _, _ = ac_state(p, t, t_step, dP)
    lo, m = band(u1, u2)
    st = initial_state(p, theta0, orbit)
    checks = dict(
        modulation_ok=bool(np.all(lo <= p.Vdc_nom - m)),
        vb0_in_band=bool(lo[0] <= st["vb0_V"] <= p.Vdc_nom - m[0]),
        pre_orbit_ok=steady_orbit(p, p.P2_pre, p.Vdc_nom, 0.0) is not None,
        post_orbit_ok=steady_orbit(p, p.P2_pre + dP, p.Vdc_nom, 0.0) is not None,
        pole_current_ok=bool(max(abs(st["i1_A"]), abs(st["i2_A"]),
                                 abs(st["ib0_A"] - st["i1_A"] - st["i2_A"])) <= p.I_pole_max),
        band_headroom_now_V=float(min(st["vb0_V"] - lo[0], p.Vdc_nom - m[0] - st["vb0_V"])),
    )
    checks["pass"] = all(v for k, v in checks.items() if k.endswith("_ok") or k == "vb0_in_band")
    return checks


# ---------------------------------------------------------------- figures
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
S1_LIGHT, BAND = "#cde2fb", "#ecebe6"
GOOD, CRIT = "#0ca30c", "#d03b3b"


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"], "font.size": 9,
        "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.labelcolor": INK2,
        "axes.titlecolor": INK, "axes.titlesize": 10, "axes.titleweight": "semibold",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2,
        "ytick.labelcolor": INK2, "lines.linewidth": 1.6, "legend.frameon": False,
        "legend.fontsize": 8, "text.color": INK,
    })
    return plt


def fig_boundary(p, R, path):
    plt = _style()
    ph = np.array(R["sweep"]["phase_deg"])
    sw = R["sweep"]
    cap = R["sweep"]["search_cap_W"]
    pa, pb_, dPt = R["pair"]["theta_a_deg"], R["pair"]["theta_b_deg"], R["pair"]["dP_test_W"]
    orb = nominal_orbit(p)

    fig, axs = plt.subplots(3, 1, figsize=(7.2, 7.6), sharex=True,
                            gridspec_kw=dict(height_ratios=[3.2, 1.1, 1.1], hspace=0.12))
    ax = axs[0]
    inner, outer = np.array(sw["inner_W"]), np.array(sw["outer_W"])
    ax.fill_between(ph, inner, outer, color=S1_LIGHT, lw=0, label="Gap between the two bounds")
    ax.plot(ph, outer, color=S1, label="Trajectory test, outer bound (nothing feasible above)")
    ax.plot(ph, inner, color=S1, lw=0.9, label="Trajectory test, inner bound (certified below)")
    ax.plot(ph, sw["energy_W"], color=S2, label="Energy sizing (blind to AC polarity)")
    ax.plot(ph, sw["static_W"], color=S3, label="Static derating (worst-case band)")
    ax.axhline(R["steady"]["C2_limit_W"], color=MUTED, lw=1.0,
               label="Steady-state and instantaneous checks")
    for th in (pa, pb_):
        ax.axvline(th, color=AXIS, lw=0.8)
    ax.plot([pa], [dPt], "o", ms=8, color=GOOD, mec=SURF, mew=2, zorder=5)
    ax.plot([pb_], [dPt], "X", ms=9, color=CRIT, mec=SURF, mew=1.5, zorder=5)
    ax.annotate(f"θa = {pa:.0f}°\nadmissible", (pa, dPt), xytext=(9, -4),
                textcoords="offset points", color=INK, fontsize=8, va="center")
    ax.annotate(f"θb = {pb_:.0f}°\nno trajectory exists", (pb_, dPt), xytext=(-9, -36),
                textcoords="offset points", color=INK, fontsize=8, va="center", ha="right",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=0, shrinkB=5))
    ax.set_ylim(0, cap * 1.04)
    ax.set_ylabel("Largest admissible AC2 step (W)")
    fig.suptitle(f"Same stored energy, different phase: AC2 step of {dPt:.0f} W",
                 x=0.12, ha="left", y=0.995, fontsize=10.5, fontweight="semibold", color=INK)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, columnspacing=1.2,
              handlelength=1.8)
    ax.text(0.995, 0.02, f"{cap:.0f} W is the search limit", transform=ax.transAxes,
            ha="right", va="bottom", color=MUTED, fontsize=7.5)

    th_o = np.rad2deg(orb["t"] * p.w)
    axs[1].plot(th_o, orb["Eb"], color=S1)
    e_pair = [orbit_energy_at(orb, p, np.deg2rad(x) / p.w) for x in (pa, pb_)]
    axs[1].plot([pa, pb_], e_pair, "o", ms=7, color=INK2, mec=SURF, mew=2)
    axs[1].set_ylabel("Buffer energy\nbefore step (J)")
    axs[1].annotate("equal energy", (pb_, e_pair[1]), xytext=(8, 6),
                    textcoords="offset points", color=INK2, fontsize=8)

    u1, u2, _, _, _ = ac_state(p, orb["t"])
    lo, m = band(u1, u2)
    axs[2].fill_between(th_o, lo, p.Vdc_nom - m, color=BAND, lw=0)
    axs[2].plot(th_o, lo, color=MUTED, lw=1.0)
    axs[2].plot(th_o, p.Vdc_nom - m, color=MUTED, lw=1.0)
    axs[2].plot(th_o, np.sqrt(2 * orb["Eb"] / p.Cb), color=S1)
    axs[2].set_ylabel("Buffer voltage\nand pole-N band (V)")
    axs[2].set_xlabel("AC1 phase at which the step arrives (deg)")
    axs[2].set_xticks(np.arange(0, 361, 45))
    axs[2].set_xlim(0, 360)
    for a in axs[1:]:
        for th in (pa, pb_):
            a.axvline(th, color=AXIS, lw=0.8)
    fig.align_ylabels(axs)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_pair_lp(p, R, ra, rb, path):
    plt = _style()
    fig, axs = plt.subplots(3, 2, figsize=(7.6, 6.4), sharex=True, sharey="row",
                            gridspec_kw=dict(hspace=0.12, wspace=0.08))
    t_max = 60.0
    titles = [f"θa = {R['pair']['theta_a_deg']:.0f}°: certified trajectory",
              f"θb = {R['pair']['theta_b_deg']:.0f}°: least-violation trajectory\n"
              f"(needs the band relaxed by {R['pair']['theta_b_relaxation_V']:.1f} V)"]
    for c, (r, title) in enumerate(zip((ra, rb), titles)):
        t = r["t"] * 1e3
        k = t <= t_max
        vb = np.sqrt(2 * r["Eb"] / p.Cb)
        Vdc = np.sqrt(2 * r["Ed"] / p.Cdc)
        hi = Vdc - r["m"]
        ax = axs[0, c]
        ax.fill_between(t[k], r["lo"][k], hi[k], color=BAND, lw=0, label="Pole-N band")
        ax.plot(t[k], vb[k], color=S1, label="Buffer voltage")
        bad = (vb < r["lo"] - 1e-6) | (vb > hi + 1e-6)
        if bad[k].any():
            ax.plot(t[k][bad[k]], vb[k][bad[k]], ".", ms=3, color=CRIT, label="Outside true band")
        ax.set_title(title, loc="left")
        tm = t[:-1]
        km = tm <= t_max
        dem = np.full_like(tm, r["P_post"])
        axs[1, c].plot(tm[km], dem[km], color=S2, label="Average demand after step")
        axs[1, c].plot(tm[km], r["pd"][km], color=S1, label="DC1 power (slew limited)")
        axs[2, c].axhline(p.Vdc_min, color=MUTED, lw=1.0)
        axs[2, c].axhline(p.Vdc_max, color=MUTED, lw=1.0)
        axs[2, c].plot(t[k], Vdc[k], color=S1)
        axs[2, c].set_xlabel("Time after step (ms)")
    axs[0, 0].set_ylabel("Voltage (V)")
    axs[1, 0].set_ylabel("Power (W)")
    axs[2, 0].set_ylabel("Bus voltage (V)")
    axs[0, 1].legend(loc="upper right")
    axs[1, 1].legend(loc="lower right")
    axs[2, 0].text(t_max, p.Vdc_min, "bus limits ", ha="right", va="bottom", color=MUTED, fontsize=7.5)
    fig.align_ylabels(axs[:, 0])
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_timedomain(p, R, runs, path):
    plt = _style()
    fig, axs = plt.subplots(3, len(runs), figsize=(9.6, 6.4), sharex=True, sharey="row",
                            gridspec_kw=dict(hspace=0.12, wspace=0.08))
    t_lo, t_hi = -5.0, 40.0
    for c, (label, s) in enumerate(runs):
        t = s["t"] * 1e3
        k = (t >= t_lo) & (t <= t_hi)
        ax = axs[0, c]
        ax.fill_between(t[k], s["lo"][k], s["hi"][k], color=BAND, lw=0, label="Pole-N band")
        ax.plot(t[k], s["vb"][k], color=S1, label="Buffer voltage")
        axs[1, c].axhline(p.Vdc_min, color=MUTED, lw=1.0)
        axs[1, c].axhline(p.Vdc_max, color=MUTED, lw=1.0)
        axs[1, c].plot(t[k], s["Vdc"][k], color=S1)
        axs[2, c].plot(t[k], s["pd"][k], color=S1)
        m = s["metrics"]
        verdict = "violation, trip" if m["violation"] else "all limits held"
        ax.set_title(f"{label}\n{verdict}", loc="left", color=INK)
        if m["t_violation_ms"] is not None:
            for a in axs[:, c]:
                a.axvline(m["t_violation_ms"], color=CRIT, lw=1.0)
            axs[1, c].annotate("first violation", (m["t_violation_ms"], p.Vdc_min),
                               xytext=(4, -14), textcoords="offset points", color=CRIT,
                               fontsize=7.5)
        axs[2, c].set_xlabel("Time after step (ms)")
    axs[0, 0].set_ylabel("Voltage (V)")
    axs[1, 0].set_ylabel("Bus voltage (V)")
    axs[2, 0].set_ylabel("DC1 power (W)")
    axs[1, 0].set_ylim(p.Vdc_min - 25, p.Vdc_max + 10)
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower right", bbox_to_anchor=(0.9, 0.995), ncol=2)
    fig.align_ylabels(axs[:, 0])
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    global OUT
    ap = argparse.ArgumentParser(description="Equal-energy pair experiment")
    ap.add_argument("--rating", choices=ratings.RATINGS, default="baseline")
    ap.add_argument("--regime", choices=("bus-loop", "fast"), default="bus-loop",
                    help="DC1 response regime for realistic ratings")
    args = ap.parse_args()
    if args.rating != "baseline":
        OUT = OUT / (args.rating + ("" if args.regime == "bus-loop" else "_fast"))
    t_start = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    p = ratings.get(args.rating, args.regime)
    orb = nominal_orbit(p)
    R = dict(params=asdict(p), rating=args.rating, regime=args.regime)

    # 1. Steady-state design facts
    cb_min = bisect(lambda c: steady_orbit(replace(p, Cb=c), p.P2_pre, p.Vdc_nom, 0.0) is None,
                    1e-6, 3e-3, 1e-7)
    C2 = post_step_limit(p, hi=max(2000.0, 4 * (p.P1 + p.P2_pre)))
    vb_orb = np.sqrt(2 * orb["Eb"] / p.Cb)
    R["steady"] = dict(Cb_min_steady_uF=cb_min * 1e6, Cb_over_min=p.Cb / cb_min,
                       orbit_vb_min_V=float(vb_orb.min()), orbit_vb_max_V=float(vb_orb.max()),
                       orbit_margin_J=float(orb["margin"]), Pdc1_pre_W=float(orb["Pd"]),
                       C2_limit_W=float(C2))
    print("steady:", R["steady"], flush=True)

    workers = max(1, min(16, (os.cpu_count() or 2) - 2))
    phases = np.arange(0, 360, PHASE_STEP_DEG)
    cap = 1.2 * C2
    modes = ("inner", "outer", "energy", "static")
    with ProcessPoolExecutor(max_workers=workers) as ex:
        # 2. Phase sweep
        tol = max(TOL_W, 0.002 * cap)
        jobs = [(p, float(th), mode, cap, tol) for mode in modes for th in phases]
        res = list(ex.map(_max_step, jobs))
        sweep = {f"{mode}_W": res[i * len(phases):(i + 1) * len(phases)]
                 for i, mode in enumerate(modes)}
        sweep.update(phase_deg=phases.tolist(), search_cap_W=cap)
        R["sweep"] = sweep
        print(f"sweep done ({time.time() - t_start:.0f} s)", flush=True)

        # 3. Pair selection
        n = len(phases)
        half = n // 2
        best = None
        for i in range(n):
            j = (i + half) % n
            lo_ = sweep["outer_W"][j]
            hi_ = min(sweep["inner_W"][i], sweep["energy_W"][i], sweep["energy_W"][j], C2)
            if np.isfinite(lo_) and np.isfinite(hi_) and hi_ > lo_:
                if best is None or hi_ - lo_ > best[2] - best[1]:
                    best = (i, lo_, hi_)
        if best is None:
            raise SystemExit("No separating equal-energy pair: the hypothesis fails here.")
        i, lo_, hi_ = best
        th_a, th_b = float(phases[i]), float(phases[(i + half) % n])
        dP = float(5 * round(0.5 * (lo_ + hi_) / 5))
        print(f"pair: {th_a} / {th_b} deg, dP = {dP} W (window {lo_:.1f}-{hi_:.1f})", flush=True)

        # 4. Verdicts and certificates at the pair
        verdicts = {}
        for th in (th_a, th_b):
            v = {mode: feasible(p, np.deg2rad(th), dP, mode, orbit=orb) for mode in modes}
            v["instantaneous"] = instantaneous_check(p, np.deg2rad(th), dP, orb)
            v["state"] = initial_state(p, np.deg2rad(th), orb)
            verdicts[f"{th:.0f}"] = v

        def needs_relax(d):
            return not feasible(p, np.deg2rad(th_b), dP, "outer", delta=d, orbit=orb)
        relax = bisect(needs_relax, 0.0, max(60.0, 0.3 * p.Vdc_nom), 0.05)

        cb_jobs = [(p, th_a, dP, "inner"), (p, th_a, dP, "outer"),
                   (p, th_b, dP, "inner"), (p, th_b, dP, "outer"), (p, th_a, dP, "conventional")]
        cb = list(ex.map(_min_cb, cb_jobs))
        R["pair"] = dict(theta_a_deg=th_a, theta_b_deg=th_b, dP_test_W=dP,
                         window_W=[float(lo_), float(hi_)], theta_b_relaxation_V=float(relax),
                         verdicts=verdicts,
                         Cb_required_uF=dict(
                             theta_a_sufficient=cb[0] * 1e6, theta_a_necessary=cb[1] * 1e6,
                             theta_b_sufficient=cb[2] * 1e6, theta_b_necessary=cb[3] * 1e6,
                             energy_sizing_plus_steady=cb[4] * 1e6, installed=p.Cb * 1e6))
        print("relaxation", relax, "Cb", R["pair"]["Cb_required_uF"], flush=True)

        # Monotonicity spot check behind the bisection
        grid = np.linspace(0, cap, 13)
        mono_jobs = [(p, th, float(d), mode, None) for th in (th_a, th_b)
                     for mode in ("inner", "outer") for d in grid]
        mono = list(ex.map(_feasible, mono_jobs))
        mono_ok = all(
            all(seq[k] >= seq[k + 1] for k in range(len(seq) - 1))
            for seq in (mono[q * 13:(q + 1) * 13] for q in range(4)))
        R["pair"]["monotone_in_dP"] = bool(mono_ok)

        # 5. Time-domain runs: tune the conventional baseline over phases first
        tune_ph = np.arange(0, 360, 15)
        nt = len(tune_ph)
        tune_jobs = [(p, float(th), dP, k) for k in K_REC_GRID_HZ for th in tune_ph]
        tune = list(ex.map(_sim_metrics, tune_jobs))
        score = {}
        for q, k in enumerate(K_REC_GRID_HZ):
            score[k] = sum(not m["violation"] for m in tune[q * nt:(q + 1) * nt])
        k_best = max(K_REC_GRID_HZ, key=lambda k: (score[k], -k))
        best_runs = tune[K_REC_GRID_HZ.index(k_best) * nt:(K_REC_GRID_HZ.index(k_best) + 1) * nt]
        R["timedomain_tuning"] = dict(k_rec_hz_grid=list(K_REC_GRID_HZ), n_phases=nt,
                                      successes=[score[k] for k in K_REC_GRID_HZ],
                                      k_rec_hz_chosen=k_best, dP_W=dP)
        print("tuning", score, "->", k_best, flush=True)

        # Does the tuned controller fail where the LP says no trajectory exists?
        lp_jobs = [(p, float(th), dP, mode, None) for mode in ("inner", "outer") for th in tune_ph]
        lp = list(ex.map(_feasible, lp_jobs))
        R["controller_vs_lp"] = [
            dict(phase_deg=float(th), conventional_ok=not best_runs[q]["violation"],
                 inner_feasible=bool(lp[q]), outer_feasible=bool(lp[nt + q]))
            for q, th in enumerate(tune_ph)]

    gains = dict(k_rec=2 * np.pi * k_best)
    plan_a = solve(p, np.deg2rad(th_a), dP, "inner", orbit=orb, smooth=True)
    pc = pole_currents(p, plan_a)
    sim_ca = simulate(p, np.deg2rad(th_a), dP, gains=gains)
    sim_cb = simulate(p, np.deg2rad(th_b), dP, gains=gains)
    sim_pa = simulate(p, np.deg2rad(th_a), dP, policy="plan", plan=plan_a, gains=gains)
    R["timedomain"] = {"conventional_theta_a": sim_ca["metrics"],
                       "conventional_theta_b": sim_cb["metrics"],
                       "plan_theta_a": sim_pa["metrics"],
                       "plan_theta_a_lp_pole_current_peak_A": float(np.max(np.abs(pc["iN"])))}

    # 6. Figures and reports
    rb = solve(p, np.deg2rad(th_b), dP, "outer", delta=relax + 0.05, orbit=orb,
               min_violation=True, smooth=True)
    fig_boundary(p, R, OUT / "fig1_boundary_vs_phase.png")
    fig_pair_lp(p, R, plan_a, rb, OUT / "fig2_pair_lp_trajectories.png")
    fig_timedomain(p, R, [(f"Conventional, θa = {th_a:.0f}°", sim_ca),
                          (f"Conventional, θb = {th_b:.0f}°", sim_cb),
                          (f"LP plan tracked, θa = {th_a:.0f}°", sim_pa)],
                   OUT / "fig3_pair_timedomain.png")
    R["runtime_s"] = time.time() - t_start
    (OUT / "results.json").write_text(json.dumps(R, indent=2, default=float))
    write_summary(p, R)
    print(f"done in {R['runtime_s']:.0f} s -> {OUT}", flush=True)


def _yn(b):
    return "admissible" if b else "not admissible"


def write_summary(p, R):
    pr, st = R["pair"], R["steady"]
    va, vb = pr["verdicts"][f"{pr['theta_a_deg']:.0f}"], pr["verdicts"][f"{pr['theta_b_deg']:.0f}"]
    sa, sb = va["state"], vb["state"]
    cb = pr["Cb_required_uF"]
    td = R["timedomain"]
    tu = R["timedomain_tuning"]
    L = []
    L.append("# Discriminating experiment: equal stored energy, different AC phase\n")
    L.append("Generated by `sim/run_pair_experiment.py`. Every number below comes from "
             "`results.json`. Nothing here has been measured on hardware.\n")
    L.append("## Scenario\n")
    L.append(f"- Bus {p.Vdc_nom:.0f} V (limits {p.Vdc_min:.0f}–{p.Vdc_max:.0f} V), "
             f"Cdc = {p.Cdc * 1e6:.0f} µF; AC1 and AC2 at {p.f:.0f} Hz in phase, "
             f"{p.E1:.0f} V and {p.E2:.0f} V peak, L = {p.L1 * 1e3:.0f} mH.")
    L.append(f"- Before the step: AC1 {p.P1:.0f} W, AC2 {p.P2_pre:.0f} W, DC2 {p.Pdc2:.0f} W, "
             f"DC1 supplies {st['Pdc1_pre_W']:.1f} W. DC1 slew limit "
             f"{p.slew_dc1 / 1e3:.0f} W/ms.")
    L.append(f"- Buffer Cb = {p.Cb * 1e6:.0f} µF, which is {st['Cb_over_min']:.2f}× the "
             f"smallest capacitance with a steady-state orbit ({st['Cb_min_steady_uF']:.0f} µF). "
             f"Nominal orbit: vb {st['orbit_vb_min_V']:.1f}–{st['orbit_vb_max_V']:.1f} V.")
    L.append(f"- The AC2 load steps up at an AC1 phase θ; its current amplitude ramps over "
             f"{p.Tr2 * 1e3:.0f} ms.\n")
    L.append("## The pair\n")
    L.append(f"θa = {pr['theta_a_deg']:.0f}° and θb = {pr['theta_b_deg']:.0f}°, "
             f"AC2 step of **{pr['dP_test_W']:.0f} W**.\n")
    L.append("| State at the step | θa | θb |\n|---|---|---|")
    for key, lab in (("Eb0_J", "Buffer energy (J)"), ("vb0_V", "Buffer voltage (V)"),
                     ("ib0_A", "Buffer current (A)"), ("Pdc1_0_W", "DC1 power (W)"),
                     ("pac0_W", "Instantaneous AC power (W)"), ("u1_V", "AC1 converter voltage (V)"),
                     ("i1_A", "AC1 current (A)")):
        L.append(f"| {lab} | {sa[key]:.4f} | {sb[key]:.4f} |")
    L.append("\nEnergy, current and power match; only the AC polarity differs.\n")
    L.append("| Criterion | θa | θb |\n|---|---|---|")
    L.append(f"| Energy sizing (rating window 0–{p.Vdc_max:.0f} V) | {_yn(va['energy'])} | "
             f"{_yn(vb['energy'])} |")
    L.append(f"| Instantaneous and steady-state checks | {_yn(va['instantaneous']['pass'])} | "
             f"{_yn(vb['instantaneous']['pass'])} |")
    L.append(f"| Static derating (worst-case band) | {_yn(va['static'])} | {_yn(vb['static'])} |")
    L.append(f"| Trajectory test, inner bound | {'certified feasible' if va['inner'] else 'not certified'}"
             f" | {'certified feasible' if vb['inner'] else 'not certified'} |")
    L.append(f"| Trajectory test, outer bound | {'not excluded' if va['outer'] else 'proved infeasible'}"
             f" | {'not excluded' if vb['outer'] else 'proved infeasible'} |")
    L.append(f"\nAt θb the band must be relaxed by **{pr['theta_b_relaxation_V']:.1f} V** before "
             "any reduced-model trajectory exists.\n")
    L.append("## Buffer needed for this step\n")
    L.append("| Requirement | Cb (µF) |\n|---|---|")
    L.append(f"| Installed | {cb['installed']:.0f} |")
    L.append(f"| Energy sizing plus steady-state check (same for both phases) | "
             f"{cb['energy_sizing_plus_steady']:.0f} |")
    L.append(f"| θa, necessary (outer) / sufficient (inner) | {cb['theta_a_necessary']:.0f} / "
             f"{cb['theta_a_sufficient']:.0f} |")
    L.append(f"| θb, necessary (outer) / sufficient (inner) | {cb['theta_b_necessary']:.0f} / "
             f"{cb['theta_b_sufficient']:.0f} |")
    L.append("\n## Phase sweep\n")
    ph = R["sweep"]["phase_deg"]
    L.append("| θ (deg) | inner (W) | outer (W) | energy sizing (W) | static (W) |\n|---|---|---|---|---|")
    for q, th in enumerate(ph):
        L.append(f"| {th} | {R['sweep']['inner_W'][q]:.0f} | {R['sweep']['outer_W'][q]:.0f} | "
                 f"{R['sweep']['energy_W'][q]:.0f} | {R['sweep']['static_W'][q]:.0f} |")
    L.append(f"\nValues of {R['sweep']['search_cap_W']:.0f} W hit the search limit "
             f"(1.2× the {st['C2_limit_W']:.0f} W steady-state limit). "
             f"Feasibility was monotone in the step size on a 13-point check at the pair: "
             f"{pr['monotone_in_dP']}.\n")
    L.append("## Time-domain check (buffer inductor, losses, nonlinear bus)\n")
    L.append(f"The conventional supervisor's recovery gain was tuned over {tu['n_phases']} phases "
             f"at {tu['dP_W']:.0f} W. Successes per gain (Hz: count): "
             f"{dict(zip(tu['k_rec_hz_grid'], tu['successes']))}; chosen {tu['k_rec_hz_chosen']} Hz.\n")
    cv = R["controller_vs_lp"]

    def lp_word(c):
        if c["inner_feasible"]:
            return "feasible"
        return "unresolved" if c["outer_feasible"] else "infeasible"

    L.append("Tuned conventional controller compared with the LP verdict at each phase:\n")
    L.append("| θ (deg) | " + " | ".join(f"{c['phase_deg']:.0f}" for c in cv) + " |")
    L.append("|---|" + "---|" * len(cv))
    L.append("| Conventional (time domain) | " +
             " | ".join("ok" if c["conventional_ok"] else "**fails**" for c in cv) + " |")
    L.append("| LP verdict (reduced model) | " + " | ".join(lp_word(c) for c in cv) + " |")
    agree = sum(c["conventional_ok"] == c["outer_feasible"] for c in cv if
                c["inner_feasible"] or not c["outer_feasible"])
    resolved = sum(1 for c in cv if c["inner_feasible"] or not c["outer_feasible"])
    ctrl_lim = [f"{c['phase_deg']:.0f}°" for c in cv if c["inner_feasible"] and not c["conventional_ok"]]
    contra = [f"{c['phase_deg']:.0f}°" for c in cv if c["conventional_ok"] and not c["outer_feasible"]]
    L.append(f"\nAt phases the LP resolves, the controller outcome matches the LP verdict at "
             f"{agree} of {resolved}. Controller-limited phases (LP feasible, controller fails): "
             f"{', '.join(ctrl_lim) if ctrl_lim else 'none'}. Contradictions (controller succeeds "
             f"where the LP proved infeasibility): {', '.join(contra) if contra else 'none'}.\n")
    L.append("Pair runs:\n")
    L.append("| Run | Bus min–max (V) | Peak pole current (A) | Smallest vN headroom (V) | "
             "Pole-N saturation (ms) | Result |\n|---|---|---|---|---|---|")
    for key, lab in (("conventional_theta_a", "Conventional, θa"),
                     ("conventional_theta_b", "Conventional, θb"),
                     ("plan_theta_a", "LP plan tracked, θa")):
        m = td[key]
        res = "all limits held" if not m["violation"] else (
            f"violation at {m['t_violation_ms']:.1f} ms, tripped")
        L.append(f"| {lab} | {m['Vdc_min']:.1f}–{m['Vdc_max']:.1f} | {m['I_pole_peak']:.1f} | "
                 f"{m['min_vN_headroom_V']:.1f} | {m['sat_time_ms']:.2f} | {res} |")
    L.append("\nvN headroom is the distance of the commanded pole-N voltage from the band edge; "
             "a negative value means the command was clipped. The buffer-capacitor voltage itself "
             "may leave the band briefly, because the inductor drop separates it from vN.\n")
    L.append("## Limits of this result\n")
    L.append("- The infeasibility proof covers the reduced averaged model (buffer inductor "
             "neglected, lossless buffer branch, 100 µs sampling). The time-domain runs test the "
             "richer averaged model but do not prove anything about it.")
    L.append("- No switched-model, dead-time or hardware evidence yet. AC current loops are ideal.")
    L.append("- One topology (the illustrative candidate of plan section 5) and one parameter set. "
             "Whether the gap survives at practical ratings across the design space is still open.")
    (OUT / "summary.md").write_text("\n".join(L), encoding="utf-8")


if __name__ == "__main__":
    main()
