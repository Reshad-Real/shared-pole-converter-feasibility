"""WP4 robustness campaign. The criteria were fixed in WP4_CRITERIA.md before any run.

Parts, each writing results/<rating>/wp4_<part>.json:
  bounds    robust phase sweep (R2, R3), the theta_b certificate over the tolerance set (R1),
            and the two tolerance corners used by `nonideal`
  events    step-down/regenerative and ramped step-up phase sweeps (reported)
  mc        misclassification campaign over random cases (R5, part of R4)
  nonideal  switched-model pair runs under grid, switching, sensor and tolerance conditions
            (R6, part of R4)

Run from the repository root after run_pair_experiment.py for the rating, bounds first:
  py -u sim/run_wp4.py --rating r120 --part bounds
"""
import argparse
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np

from model import bisect, steady_orbit
from feasibility import feasible, nominal_orbit, solve
from switched import SwitchParams, simulate_switched
from timedomain import GAINS
from run_pair_experiment import TOL_W, instantaneous_check
import ratings
import robust

RES = Path(__file__).resolve().parent / "results"
PHASES = np.arange(0, 360, 9)
N_MC = 300
MC_SEED = {"r120": 2026, "r230": 2027}
R2_PCT, R3_PCT = 15.0, 10.0
WORKERS = 8


def _json_default(o):
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def _save(rating, part, R):
    path = RES / rating / f"wp4_{part}.json"
    path.write_text(json.dumps(R, indent=1, default=_json_default))
    print("wrote", path, flush=True)


def _load(rating, part):
    return json.loads((RES / rating / f"wp4_{part}.json").read_text())


# ---------------------------------------------------------------- parallel workers
def _max_step(args):
    """Largest step magnitude (W) in direction `sign` that `mode` admits; 0 if none."""
    p, th_deg, mode, cap, tol, sign = args
    orb = nominal_orbit(p)
    if orb is None:
        return 0.0
    r = bisect(lambda d: feasible(p, np.deg2rad(th_deg), sign * d, mode, orbit=orb),
               0.0, cap, tol)
    return 0.0 if np.isnan(r) else float(r)


def _steady_limit(args):
    """Largest step magnitude in direction `sign` with a post-step steady orbit (C2)."""
    p, cap, tol, sign = args
    r = bisect(lambda d: steady_orbit(p, p.P2_pre + sign * d, p.Vdc_nom, 0.0) is not None,
               0.0, cap, tol)
    return 0.0 if np.isnan(r) else float(r)


def _feasible(args):
    p, th_deg, dP, mode = args
    return bool(feasible(p, np.deg2rad(th_deg), dP, mode))


def _outer_from(p, th, dP, x0):
    """Outer LP of plant p started from a given state (Eb0, Ed0, pd0)."""
    r = solve(p, th, dP, "outer", x0=dict(Eb0=float(x0["Eb0"]), Ed0=float(x0["Ed0"]),
                                          pd0=float(x0["pd0"])))
    if "warning" in r:
        raise RuntimeError(r["warning"])
    return bool(r["feasible"])


def _metrics(s):
    m = dict(s["metrics"])
    m["state_at_step"] = s["state_at_step"]
    return m


def _mc_case(args):
    """Truth and every criterion's verdict for one random case."""
    p, c, C, D8, C2, k_rec = args
    th, dP = np.deg2rad(c["theta_deg"]), c["dP_W"]
    out = dict(c)
    try:
        pp = robust.perturb(p, **c["factors"])
        orb_p = nominal_orbit(pp)
        if orb_p is None:
            out["truth"] = "unresolved"
            out["truth_note"] = "no centred pre-step orbit"
        elif not feasible(pp, th, dP, "outer", orbit=orb_p):
            out["truth"] = "infeasible"
        else:
            plan = solve(pp, th, dP, "inner", orbit=orb_p, smooth=True)
            if not plan["feasible"]:
                out["truth"] = "unresolved"
                out["truth_note"] = "inner LP infeasible"
            else:
                s = simulate_switched(pp, SwitchParams(), th, dP, policy="plan", plan=plan,
                                      gains=dict(k_rec=2 * np.pi * k_rec))
                m = _metrics(s)
                out["plan_run"] = m
                out["truth"] = "unresolved" if m["violation"] else "feasible"
                if not m["violation"]:
                    out["outer_from_run_state"] = _outer_from(pp, th, dP, m["state_at_step"])

        orb = nominal_orbit(p)
        v = {}
        v["energy_sizing"] = bool(feasible(p, th, dP, "energy", orbit=orb) and dP <= C2)
        v["instantaneous"] = bool(instantaneous_check(p, th, dP, orb)["pass"])
        v["static_derating"] = bool(feasible(p, th, dP, "static", orbit=orb))
        v["nominal_outer"] = bool(feasible(p, th, dP, "outer", orbit=orb))
        v["nominal_inner"] = bool(feasible(p, th, dP, "inner", orbit=orb))
        v["robust_outer"] = any(feasible(robust.perturb(p, **f), th, dP, "outer")
                                for f in C + D8)
        v["robust_inner"] = all(feasible(robust.perturb(p, **f), th, dP, "inner") for f in C)
        out["admits"] = v
    except RuntimeError as e:
        out["error"] = str(e)
    return out


def _ni_run(args):
    """One switched-model run of the non-ideality campaign, with its LP checks."""
    cond, run, th_deg, dP, k_rec = args
    th = np.deg2rad(th_deg)
    out = dict(condition=cond["key"], run=run, theta_deg=th_deg)
    policy = "plan" if run == "plan_theta_a" else "conventional"
    plan = None
    if policy == "plan":
        plan = solve(cond["plan_p"], th, dP, "inner", smooth=True)
        if not plan["feasible"]:
            out["no_plan"] = True
            return out
    gains = dict(k_rec=2 * np.pi * k_rec, wc=GAINS["wc"] * cond["wc_scale"])
    s = simulate_switched(cond["pc"], cond["sw"], th, dP, policy=policy, plan=plan,
                          pre_cycles=cond["pre_cycles"], gains=gains, plant=cond["pp"])
    out["metrics"] = _metrics(s)
    actual = cond["pp"] if cond["pp"] is not None else cond["pc"]
    try:
        out["outer_from_run_state"] = _outer_from(actual, th, dP, s["state_at_step"])
        out["outer_from_orbit"] = bool(feasible(actual, th, dP, "outer"))
    except RuntimeError as e:
        out["lp_error"] = str(e)
    return out


# ---------------------------------------------------------------- parts
def _wp1(rating):
    W = json.loads((RES / rating / "results.json").read_text())
    pr = W["pair"]
    return W, pr["theta_a_deg"], pr["theta_b_deg"], pr["dP_test_W"]


def part_bounds(rating, p, ex):
    W, th_a, th_b, dP = _wp1(rating)
    cap = W["sweep"]["search_cap_W"]
    tol = max(TOL_W, 0.002 * cap)
    C, D = robust.corners(), robust.draws(robust.N_DRAWS_R1)
    plants_out = C + D[:robust.N_DRAWS_SWEEP]
    jobs = [(robust.perturb(p, **f), float(th), "outer", cap, tol, 1.0)
            for f in plants_out for th in PHASES]
    jobs += [(robust.perturb(p, **f), float(th), "inner", cap, tol, 1.0)
             for f in C for th in PHASES]
    t0 = time.time()
    res = list(ex.map(_max_step, jobs, chunksize=2))
    n_o = len(plants_out) * len(PHASES)
    outer = np.array(res[:n_o]).reshape(len(plants_out), len(PHASES))
    inner = np.array(res[n_o:]).reshape(len(C), len(PHASES))
    print(f"robust sweep {time.time() - t0:.0f} s", flush=True)

    sw = W["sweep"]
    C2 = W["steady"]["C2_limit_W"]
    conv = np.minimum(np.array(sw["energy_W"], float), C2)
    rob_out, rob_in = outer.max(axis=0), inner.min(axis=0)
    over = conv - rob_out
    k = int(np.argmax(over))
    gap_corner = outer[:len(C), k] - inner[:, k]
    R2 = dict(phase_deg=float(PHASES[k]), conventional_W=float(conv[k]),
              robust_outer_W=float(rob_out[k]), nominal_outer_W=float(sw["outer_W"][k]),
              over_W=float(over[k]), over_pct=float(100 * over[k] / conv[k]))
    R2["pass"] = R2["over_pct"] >= R2_PCT
    R3 = dict(phase_deg=float(PHASES[k]), robust_outer_W=float(rob_out[k]),
              robust_inner_W=float(rob_in[k]), gap_W=float(rob_out[k] - rob_in[k]),
              gap_pct_of_conventional=float(100 * (rob_out[k] - rob_in[k]) / conv[k]),
              plant_by_plant_gap_max_W=float(gap_corner.max()),
              plant_by_plant_gap_mean_W=float(gap_corner.mean()),
              spread_outer_corners_W=[float(outer[:len(C), k].min()),
                                      float(outer[:len(C), k].max())])
    R3["pass"] = R3["gap_pct_of_conventional"] <= R3_PCT
    # Nominal over-admission on the same grid, for comparison
    k_n = int(np.argmax(conv - np.array(sw["outer_W"], float)))

    r1 = list(ex.map(_feasible, [(robust.perturb(p, **f), th_b, dP, "outer") for f in C + D]))
    R1 = dict(theta_b_deg=th_b, dP_W=dP, n_plants=len(r1), n_outer_feasible=int(sum(r1)),
              feasible_plants=[robust.label(f) for f, ok in zip(C + D, r1) if ok])
    R1["pass"] = R1["n_outer_feasible"] == 0

    i_a = int(np.where(PHASES == th_a)[0][0])
    oa = outer[:len(C), i_a]
    corner_lo, corner_hi = C[int(np.argmin(oa))], C[int(np.argmax(oa))]
    R = dict(rating=rating, phases_deg=PHASES.tolist(), cap_W=cap, tol_W=tol, C2_W=C2,
             corners=C, draws=D, plants_outer=[robust.label(f) for f in plants_out],
             outer_W=outer, inner_W=inner, robust_outer_W=rob_out, robust_inner_W=rob_in,
             conventional_W=conv, nominal_outer_W=sw["outer_W"], nominal_inner_W=sw["inner_W"],
             nominal_worst=dict(phase_deg=float(PHASES[k_n]),
                                over_pct=float(100 * (conv[k_n] - sw["outer_W"][k_n]) / conv[k_n])),
             R1=R1, R2=R2, R3=R3,
             corner_lo=corner_lo, corner_hi=corner_hi,
             corner_outer_at_theta_a_W=[float(oa.min()), float(oa.max())])
    print("R1", R1["pass"], R1["n_outer_feasible"], "| R2", R2, "| R3", R3, flush=True)
    _save(rating, "bounds", R)


def part_events(rating, p, ex):
    W, th_a, th_b, dP = _wp1(rating)
    cap_up = W["sweep"]["search_cap_W"]
    events = {"step_down": dict(sign=-1.0, Tr2=p.Tr2, cap=2.0 * p.P2_pre),
              "ramp10": dict(sign=1.0, Tr2=10e-3, cap=cap_up)}
    C, D = robust.corners(), robust.draws(robust.N_DRAWS_R1)
    plants = C + D[:robust.N_DRAWS_SWEEP]
    R = dict(rating=rating, phases_deg=PHASES.tolist())
    for name, e in events.items():
        t0 = time.time()
        q = replace(p, Tr2=e["Tr2"])
        tol = max(TOL_W, 0.002 * e["cap"])
        modes = ("inner", "outer", "energy", "static")
        jobs = [(q, float(th), m, e["cap"], tol, e["sign"]) for m in modes for th in PHASES]
        jobs += [(robust.perturb(q, **f), float(th), "outer", e["cap"], tol, e["sign"])
                 for f in plants for th in PHASES]
        res = list(ex.map(_max_step, jobs, chunksize=2))
        n = len(PHASES)
        nom = {m: res[i * n:(i + 1) * n] for i, m in enumerate(modes)}
        rob = np.array(res[4 * n:]).reshape(len(plants), n).max(axis=0)
        C2 = _steady_limit((q, e["cap"], tol, e["sign"]))
        conv = np.minimum(np.array(nom["energy"]), C2)
        over_n = conv - np.array(nom["outer"])
        over_r = conv - rob
        kn, kr = int(np.argmax(over_n)), int(np.argmax(over_r))
        # Monotonicity spot check at the nominal worst phase
        grid = np.linspace(0, e["cap"], 13)
        mono_jobs = [(q, float(PHASES[kn]), float(e["sign"] * d), m)
                     for m in ("inner", "outer") for d in grid]
        mono = list(ex.map(_feasible, mono_jobs))
        mono_ok = all(all(s[i] >= s[i + 1] for i in range(12)) for s in (mono[:13], mono[13:]))
        R[name] = dict(sign=e["sign"], Tr2_s=e["Tr2"], cap_W=e["cap"], C2_W=C2,
                       nominal=nom, robust_outer_W=rob, conventional_W=conv,
                       worst_nominal=dict(phase_deg=float(PHASES[kn]), conventional_W=float(conv[kn]),
                                          outer_W=float(nom["outer"][kn]), over_W=float(over_n[kn]),
                                          over_pct=float(100 * over_n[kn] / conv[kn]) if conv[kn] > 0 else 0.0),
                       worst_robust=dict(phase_deg=float(PHASES[kr]), conventional_W=float(conv[kr]),
                                         robust_outer_W=float(rob[kr]), over_W=float(over_r[kr]),
                                         over_pct=float(100 * over_r[kr] / conv[kr]) if conv[kr] > 0 else 0.0),
                       monotone_spot_check=bool(mono_ok))
        print(name, f"{time.time() - t0:.0f} s", R[name]["worst_nominal"], R[name]["worst_robust"],
              "monotone", mono_ok, flush=True)
    _save(rating, "events", R)


def part_mc(rating, p, ex):
    W, th_a, th_b, dP = _wp1(rating)
    cap = W["sweep"]["search_cap_W"]
    C2 = W["steady"]["C2_limit_W"]
    k_rec = W["timedomain_tuning"]["k_rec_hz_chosen"]
    rng = np.random.default_rng(MC_SEED[rating])
    th = rng.uniform(0, 360, N_MC)
    d = rng.uniform(0, cap, N_MC)
    plants = robust.draws(N_MC, seed=MC_SEED[rating] + 100)
    C, D8 = robust.corners(), robust.draws(robust.N_DRAWS_R1)[:robust.N_DRAWS_SWEEP]
    cases = [dict(index=i, theta_deg=float(th[i]), dP_W=float(d[i]), factors=plants[i])
             for i in range(N_MC)]
    t0 = time.time()
    res = []
    for r in ex.map(_mc_case, [(p, c, C, D8, C2, k_rec) for c in cases], chunksize=1):
        res.append(r)
        if len(res) % 50 == 0:
            print(f"  {len(res)} cases, {time.time() - t0:.0f} s", flush=True)
    R = dict(rating=rating, n=N_MC, seed=MC_SEED[rating], cap_W=cap, C2_W=C2, k_rec_hz=k_rec,
             cases=res)
    _save(rating, "mc", R)


def part_nonideal(rating, p, ex):
    W, th_a, th_b, dP = _wp1(rating)
    k_rec = W["timedomain_tuning"]["k_rec_hz_chosen"]
    B = _load(rating, "bounds")
    conds = robust.conditions(p, B["corner_lo"], B["corner_hi"])
    runs = (("conventional_theta_a", th_a), ("plan_theta_a", th_a), ("conventional_theta_b", th_b))
    jobs = [(c, r, th, dP, k_rec) for c in conds for r, th in runs]
    t0 = time.time()
    res = list(ex.map(_ni_run, jobs, chunksize=1))
    print(f"non-ideality runs {time.time() - t0:.0f} s", flush=True)
    R = dict(rating=rating, theta_a_deg=th_a, theta_b_deg=th_b, dP_W=dP, k_rec_hz=k_rec,
             conditions=[dict(key=c["key"], label=c["label"], group=c["group"],
                              wc_scale=c["wc_scale"], pre_cycles=c["pre_cycles"],
                              switch=asdict(c["sw"])) for c in conds],
             runs=res)
    _save(rating, "nonideal", R)


def _r4_diagnosis(args):
    """For a run that held every limit although the outer LP says no: the band relaxation the
    LP needs, and how far the run's buffer voltage actually left the band."""
    cond, r, dP, k_rec = args
    m = r["metrics"]
    act = cond["pp"] if cond["pp"] is not None else cond["pc"]
    th = np.deg2rad(r["theta_deg"])
    x0 = {k: float(m["state_at_step"][k]) for k in ("Eb0", "Ed0", "pd0")}
    dstar = bisect(lambda d: not solve(act, th, dP, "outer", delta=d, x0=x0)["feasible"],
                   0.0, 80.0, 0.05)
    s = simulate_switched(cond["pc"], cond["sw"], th, dP, pre_cycles=cond["pre_cycles"],
                          gains=dict(k_rec=2 * np.pi * k_rec, wc=GAINS["wc"] * cond["wc_scale"]),
                          plant=cond["pp"])
    from model import ac_state, band
    k = s["t"] >= 0
    t_step = (th % (2 * np.pi)) / cond["pc"].w + cond["pre_cycles"] * cond["pc"].T
    u1, u2, _, _, _ = ac_state(act, s["t"][k] + t_step, t_step, dP)
    lo, mm = band(u1, u2)
    hi = s["Vdc_avg"][k] - mm
    exc_vb = np.maximum(lo - s["vb_avg"][k], s["vb_avg"][k] - hi)
    exc_vN = np.maximum(lo - s["vN_avg"][k], s["vN_avg"][k] - hi)
    return dict(condition=r["condition"], run=r["run"], relaxation_needed_V=float(dstar),
                vb_outside_band_max_V=float(exc_vb.max()),
                vN_outside_reference_band_max_V=float(exc_vN.max()),
                sat_time_ms=m["sat_time_ms"], ac_current_rms_error_A=m["ac_current_rms_error_A"])


def part_diagnostics(rating, p, ex):
    """Post-hoc diagnostics (not part of the pass marks): phase-blind limits per corner plant,
    one-at-a-time tolerance sensitivity of the outer bound, and the R4 cases."""
    W, th_a, th_b, dP = _wp1(rating)
    B = _load(rating, "bounds")
    cap, tol = B["cap_W"], B["tol_W"]
    C = robust.corners()
    t0 = time.time()
    jobs = [(robust.perturb(p, **f), float(th), "energy", cap, tol, 1.0) for f in C for th in PHASES]
    en = np.array(list(ex.map(_max_step, jobs, chunksize=2))).reshape(len(C), len(PHASES))
    c2 = np.array(list(ex.map(_steady_limit, [(robust.perturb(p, **f), cap, tol, 1.0) for f in C])))
    outer = np.array(B["outer_W"])[:len(C)]
    conv = np.minimum(en, c2[:, None])
    over = conv - outer
    per_plant = []
    for i, f in enumerate(C):
        k = int(np.argmax(over[i]))
        per_plant.append(dict(plant=robust.label(f), phase_deg=float(PHASES[k]),
                              conventional_W=float(conv[i, k]), outer_W=float(outer[i, k]),
                              over_W=float(over[i, k]), over_pct=float(100 * over[i, k] / conv[i, k])))
    print(f"per-plant sweep {time.time() - t0:.0f} s", flush=True)

    ph = B["R2"]["phase_deg"]
    keys = ("Cb", "Cdc", "E", "slew", "L1", "L2", "Lb", "Rb")
    oat = []
    for key in keys:
        for s in (-1, 1):
            if key == "Rb":
                fac = robust.RB_RANGE[0] if s < 0 else robust.RB_RANGE[1]
            else:
                fac = 1.0 + s * robust.TOL["L" if key in ("L1", "L2", "Lb") else key]
            oat.append({key: fac})
    res = list(ex.map(_max_step, [(robust.perturb(p, **f), ph, "outer", cap, tol, 1.0) for f in oat]))
    nominal = float(np.array(B["nominal_outer_W"])[list(PHASES).index(ph)])
    sens = [dict(change=robust.label(f), outer_W=float(r), delta_W=float(r - nominal))
            for f, r in zip(oat, res)]

    conds = {c["key"]: c for c in robust.conditions(p, B["corner_lo"], B["corner_hi"])}
    N = _load(rating, "nonideal")
    bad = [r for r in N["runs"] if r.get("metrics") and not r["metrics"]["violation"]
           and r.get("outer_from_run_state") is False]
    k_rec = W["timedomain_tuning"]["k_rec_hz_chosen"]
    r4 = list(ex.map(_r4_diagnosis, [(conds[r["condition"]], r, dP, k_rec) for r in bad]))

    R = dict(rating=rating, per_plant=per_plant, per_plant_min_over_pct=min(x["over_pct"] for x in per_plant),
             oat_phase_deg=ph, oat_nominal_outer_W=nominal, oat=sens, r4=r4,
             theta_b_relaxation_nominal_V=W["pair"]["theta_b_relaxation_V"])
    print("per-plant min over %", R["per_plant_min_over_pct"], "| oat", sens, "| r4", r4, flush=True)
    _save(rating, "diagnostics", R)


E_LEVELS = (0.9, 1.0, 1.1)


def part_conditional(rating, p, ex):
    """Post hoc: R1-R3 with the grid amplitude measured (fixed at each level) and only the
    component and DC1 tolerances unknown. Phase-blind checks use the measured amplitude and
    nominal components; robust bounds range over the 8 corners of (Cb, Cdc, slew)."""
    import itertools
    W, th_a, th_b, dP = _wp1(rating)
    B = _load(rating, "bounds")
    cap, tol = B["cap_W"], B["tol_W"]
    keys = ("Cb", "Cdc", "slew")
    C8 = [{k: 1.0 + s * robust.TOL[k] for k, s in zip(keys, signs)}
          for signs in itertools.product((-1, 1), repeat=3)]
    R = dict(rating=rating, levels={})
    for e in E_LEVELS:
        t0 = time.time()
        q = robust.perturb(p, E=e)
        plants = [robust.perturb(q, **f) for f in C8]
        jobs = [(q, float(th), "energy", cap, tol, 1.0) for th in PHASES]
        jobs += [(pl, float(th), "outer", cap, tol, 1.0) for pl in plants for th in PHASES]
        jobs += [(pl, float(th), "inner", cap, tol, 1.0) for pl in plants for th in PHASES]
        res = np.array(list(ex.map(_max_step, jobs, chunksize=2)))
        n = len(PHASES)
        en = res[:n]
        outer = res[n:n + 8 * n].reshape(8, n)
        inner = res[n + 8 * n:].reshape(8, n)
        C2 = _steady_limit((q, cap, tol, 1.0))
        conv = np.minimum(en, C2)
        rob_out, rob_in = outer.max(axis=0), inner.min(axis=0)
        over = conv - rob_out
        k = int(np.argmax(over))
        draws = [dict(f, E=e) for f in robust.draws(robust.N_DRAWS_R1)]
        r1_plants = C8 + draws
        r1 = list(ex.map(_feasible, [(robust.perturb(p, **dict(f, E=e)), th_b, dP, "outer")
                                     for f in r1_plants]))
        R["levels"][f"{e:.1f}"] = dict(
            E=e, conventional_W=conv, robust_outer_W=rob_out, robust_inner_W=rob_in,
            R1=dict(n_plants=len(r1), n_outer_feasible=int(sum(r1))),
            R2=dict(phase_deg=float(PHASES[k]), conventional_W=float(conv[k]),
                    robust_outer_W=float(rob_out[k]), over_W=float(over[k]),
                    over_pct=float(100 * over[k] / conv[k])),
            R3=dict(robust_inner_W=float(rob_in[k]), gap_W=float(rob_out[k] - rob_in[k]),
                    gap_pct=float(100 * (rob_out[k] - rob_in[k]) / conv[k])))
        print(f"E x{e}: {time.time() - t0:.0f} s", R["levels"][f"{e:.1f}"]["R1"],
              R["levels"][f"{e:.1f}"]["R2"], R["levels"][f"{e:.1f}"]["R3"], flush=True)
    _save(rating, "conditional", R)


def part_deadtime40(rating, p, ex):
    """Post hoc: the 40 kHz pair again with a dead time sized for it (0.25 us instead of 1 us)."""
    W, th_a, th_b, dP = _wp1(rating)
    B = _load(rating, "bounds")
    c = dict(next(c for c in robust.conditions(p, B["corner_lo"], B["corner_hi"])
                  if c["key"] == "fsw40"))
    c["sw"] = replace(c["sw"], td=0.25e-6)
    c["key"] = "fsw40_td025"
    k_rec = W["timedomain_tuning"]["k_rec_hz_chosen"]
    runs = (("conventional_theta_a", th_a), ("plan_theta_a", th_a), ("conventional_theta_b", th_b))
    res = list(ex.map(_ni_run, [(c, r, th, dP, k_rec) for r, th in runs]))
    _save(rating, "deadtime40", dict(rating=rating, td_us=0.25, fsw_hz=c["sw"].fsw, runs=res))


# ---------------------------------------------------------------- circuit outer bound (R4 follow-up)
def _circuit_check(act, th, dP, st):
    """Circuit outer LP of plant `act` from a run's state at the step."""
    import outer_circuit as oc
    x0 = dict(vb=float(st["vb"]), ib=float(st["ib"]), Ed0=float(st["Ed0"]), pd0=float(st["pd0"]))
    return bool(oc.feasible(act, th, dP, x0=x0))


def _r4c_nonideal(args):
    """Re-run one non-ideality run and, if it held, check it against the circuit bound."""
    cond, run, th_deg, dP, k_rec, stored = args
    r = _ni_run((cond, run, th_deg, dP, k_rec))
    out = dict(source="nonideal", condition=cond["key"], run=run)
    if "metrics" not in r:
        out["skipped"] = "no plan"
        return out
    m = r["metrics"]
    out.update(held=not m["violation"], reproduced=abs(m["Vdc_min"] - stored["Vdc_min"]) < 1e-6,
               reduced_outer=r.get("outer_from_run_state"))
    if out["held"]:
        act = cond["pp"] if cond["pp"] is not None else cond["pc"]
        out["circuit_outer"] = _circuit_check(act, np.deg2rad(th_deg), dP, m["state_at_step"])
    return out


def _r4c_mc(args):
    """Re-run the truth plan run of one MC case and check it against the circuit bound."""
    p, c, k_rec = args
    th, dP = np.deg2rad(c["theta_deg"]), c["dP_W"]
    pp = robust.perturb(p, **c["factors"])
    plan = solve(pp, th, dP, "inner", orbit=nominal_orbit(pp), smooth=True)
    s = simulate_switched(pp, SwitchParams(), th, dP, policy="plan", plan=plan,
                          gains=dict(k_rec=2 * np.pi * k_rec))
    m = s["metrics"]
    out = dict(source="mc", index=c["index"], held=not m["violation"],
               reproduced=abs(m["Vdc_min"] - c["plan_run"]["Vdc_min"]) < 1e-6,
               reduced_outer=c.get("outer_from_run_state"))
    if out["held"]:
        out["circuit_outer"] = _circuit_check(pp, th, dP, s["state_at_step"])
    return out


def part_r4circuit(rating, p, ex):
    """R4 again, with the circuit outer bound: every WP4 run that held its limits."""
    W, th_a, th_b, dP = _wp1(rating)
    B, N, M = _load(rating, "bounds"), _load(rating, "nonideal"), _load(rating, "mc")
    k_rec = W["timedomain_tuning"]["k_rec_hz_chosen"]
    conds = {c["key"]: c for c in robust.conditions(p, B["corner_lo"], B["corner_hi"])}
    jobs = [(conds[r["condition"]], r["run"], r["theta_deg"], dP, k_rec, r["metrics"])
            for r in N["runs"] if r.get("metrics") and not r["metrics"]["violation"]]
    T = _load(rating, "deadtime40")
    c40 = dict(conds["fsw40"])
    c40["sw"], c40["key"] = replace(c40["sw"], td=0.25e-6), "fsw40_td025"
    jobs += [(c40, r["run"], r["theta_deg"], dP, k_rec, r["metrics"])
             for r in T["runs"] if r.get("metrics") and not r["metrics"]["violation"]]
    mc_jobs = [(p, c, k_rec) for c in M["cases"] if c.get("truth") == "feasible"]
    t0 = time.time()
    res = list(ex.map(_r4c_nonideal, jobs, chunksize=1))
    res += list(ex.map(_r4c_mc, mc_jobs, chunksize=1))
    held = [r for r in res if r.get("held")]
    R = dict(rating=rating, runs=res, n_held=len(held),
             n_reproduced=sum(r["reproduced"] for r in res if "reproduced" in r),
             n_checked=len([r for r in res if "reproduced" in r]),
             contradictions_reduced=[r for r in held if r.get("reduced_outer") is False],
             contradictions_circuit=[r for r in held if r.get("circuit_outer") is False])
    print(f"r4circuit {time.time() - t0:.0f} s: held {len(held)}, reproduced {R['n_reproduced']}/"
          f"{R['n_checked']}, contradictions reduced {len(R['contradictions_reduced'])}, "
          f"circuit {len(R['contradictions_circuit'])}", flush=True)
    _save(rating, "r4circuit", R)


def _max_step_circuit(args):
    import outer_circuit as oc
    p, th_deg, cap, tol = args
    r = bisect(lambda d: oc.feasible(p, np.deg2rad(th_deg), d), 0.0, cap, tol)
    return 0.0 if np.isnan(r) else float(r)


def _feasible_circuit(args):
    import outer_circuit as oc
    p, th_deg, dP = args
    return bool(oc.feasible(p, np.deg2rad(th_deg), dP))


def _relax_circuit(args):
    import outer_circuit as oc
    p, th_deg, dP = args
    return float(bisect(lambda d: not oc.feasible(p, np.deg2rad(th_deg), dP, delta=d),
                        0.0, 120.0, 0.1))


def part_circuit(rating, p, ex):
    """Headline numbers again with the circuit outer bound: the nominal phase sweep and G1
    over-admission, the theta_b margin, and the measured-grid tolerance analysis (R1, R2)."""
    import itertools
    W, th_a, th_b, dP = _wp1(rating)
    B, CD = _load(rating, "bounds"), _load(rating, "conditional")
    cap, tol = B["cap_W"], B["tol_W"]
    t0 = time.time()
    outer = np.array(list(ex.map(_max_step_circuit, [(p, float(th), cap, tol) for th in PHASES],
                                 chunksize=1)))
    conv = np.array(B["conventional_W"])
    over = conv - outer
    k = int(np.argmax(over))
    relax = _relax_circuit((p, th_b, dP))
    print(f"circuit sweep {time.time() - t0:.0f} s: worst {PHASES[k]} deg, over "
          f"{100 * over[k] / conv[k]:.1f} %, theta_b relaxation {relax:.1f} V", flush=True)
    keys = ("Cb", "Cdc", "slew")
    C8 = [{kk: 1.0 + s * robust.TOL[kk] for kk, s in zip(keys, signs)}
          for signs in itertools.product((-1, 1), repeat=3)]
    levels = {}
    for key, lv in CD["levels"].items():
        e = lv["E"]
        i0 = list(PHASES).index(lv["R2"]["phase_deg"])
        win = [float(PHASES[(i0 + j) % len(PHASES)]) for j in range(-5, 6)]
        plants = [robust.perturb(p, **dict(f, E=e)) for f in C8]
        res = np.array(list(ex.map(_max_step_circuit, [(pl, th, cap, tol) for pl in plants
                                                       for th in win], chunksize=1)))
        rob = res.reshape(len(plants), len(win)).max(axis=0)
        cv = np.array([lv["conventional_W"][list(PHASES).index(th)] for th in win])
        ov = cv - rob
        j = int(np.argmax(ov))
        draws = [dict(f, E=e) for f in robust.draws(robust.N_DRAWS_R1)]
        r1 = list(ex.map(_feasible_circuit, [(robust.perturb(p, **dict(f, E=e)), th_b, dP)
                                             for f in C8 + draws], chunksize=1))
        levels[key] = dict(E=e, window_deg=win, robust_outer_W=rob, conventional_W=cv,
                           R2=dict(phase_deg=win[j], conventional_W=float(cv[j]),
                                   robust_outer_W=float(rob[j]), over_W=float(ov[j]),
                                   over_pct=float(100 * ov[j] / cv[j])),
                           R1=dict(n_plants=len(r1), n_outer_feasible=int(sum(r1))))
        print(f"  E x{e}: {time.time() - t0:.0f} s", levels[key]["R2"], levels[key]["R1"], flush=True)
    R = dict(rating=rating, phases_deg=PHASES.tolist(), outer_W=outer, conventional_W=conv,
             reduced_outer_W=B["nominal_outer_W"],
             worst=dict(phase_deg=float(PHASES[k]), conventional_W=float(conv[k]),
                        outer_W=float(outer[k]), over_W=float(over[k]),
                        over_pct=float(100 * over[k] / conv[k])),
             theta_b=dict(theta_deg=th_b, dP_W=dP, relaxation_V=relax,
                          reduced_relaxation_V=W["pair"]["theta_b_relaxation_V"]),
             measured_grid=levels)
    _save(rating, "circuit", R)


def _mcc_case(args):
    """Circuit outer LP of the case's own plant (from its orbit) and of the nominal plant."""
    import outer_circuit as oc
    p, c = args
    th, dP = np.deg2rad(c["theta_deg"]), c["dP_W"]
    pp = robust.perturb(p, **c["factors"])
    own = bool(oc.feasible(pp, th, dP)) if c.get("truth") == "infeasible" else None
    return dict(index=c["index"], plant_circuit_outer=own,
                nominal_circuit_outer=bool(oc.feasible(p, th, dP)))


def part_mccircuit(rating, p, ex):
    """R5 again with the circuit bound: a case is infeasible only if its own plant's circuit
    outer LP is infeasible; the nominal circuit outer bound is scored as a criterion."""
    M = _load(rating, "mc")
    t0 = time.time()
    res = list(ex.map(_mcc_case, [(p, c) for c in M["cases"] if "error" not in c], chunksize=1))
    print(f"mccircuit {time.time() - t0:.0f} s: reduced-infeasible cases still infeasible: "
          f"{sum(r['plant_circuit_outer'] is False for r in res)} of "
          f"{sum(r['plant_circuit_outer'] is not None for r in res)}", flush=True)
    _save(rating, "mccircuit", dict(rating=rating, cases=res))


def part_circuitcheck(rating, p, ex):
    """Consistency of the circuit bound: with a small buffer inductor, the theta_b relaxation
    it needs should approach the reduced LP's; with more tangent cuts it should barely move."""
    import outer_circuit as oc
    W, th_a, th_b, dP = _wp1(rating)
    jobs = [(replace(p, Lb=L), th_b, dP) for L in (20e-6, 100e-6, p.Lb)]
    res = list(ex.map(_relax_circuit, jobs))
    oc.N_CUTS_V, oc.N_CUTS_I = 24, 16
    fine = _relax_circuit((p, th_b, dP))
    R = dict(rating=rating, reduced_relaxation_V=W["pair"]["theta_b_relaxation_V"],
             by_Lb=[dict(Lb_H=L, relaxation_V=r) for (L, r) in
                    zip((20e-6, 100e-6, p.Lb), res)],
             nominal_Lb_24_16_cuts_V=fine)
    print("circuitcheck", R, flush=True)
    _save(rating, "circuitcheck", R)


PARTS = dict(bounds=part_bounds, events=part_events, mc=part_mc, nonideal=part_nonideal,
             diagnostics=part_diagnostics, conditional=part_conditional,
             deadtime40=part_deadtime40, r4circuit=part_r4circuit, circuit=part_circuit,
             circuitcheck=part_circuitcheck, mccircuit=part_mccircuit)


def main():
    ap = argparse.ArgumentParser(description="WP4 robustness campaign")
    ap.add_argument("--rating", choices=("r120", "r230"), required=True)
    ap.add_argument("--part", choices=tuple(PARTS), nargs="+", required=True)
    args = ap.parse_args()
    p = ratings.get(args.rating)
    workers = max(1, min(WORKERS, (os.cpu_count() or 2) - 2))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for part in args.part:
            t0 = time.time()
            print(f"== {args.rating} {part}", flush=True)
            PARTS[part](args.rating, p, ex)
            print(f"== {part} done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
