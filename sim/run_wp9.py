"""WP9 digital-twin validation against an emulated plant (protocol: WP9_TWIN_PROTOCOL.md).

For every case of the WP4 misclassification campaign (results/<rating>/wp4_mc.json) the twin
is calibrated on the case's plant from two recordings of normal operation, then predicts the
case it has not seen. Writes results/<rating>/wp9_twin.json.

Run from the repository root after the WP4 parts mc and mccircuit:
  py -u sim/run_wp9.py --rating r120
"""
import argparse
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path

import numpy as np

import outer_circuit as oc
import ratings
import robust
import twin as tw
from feasibility import feasible, nominal_orbit, solve
from model import bisect
from run_wp4 import _json_default, _wp1
from switched import SwitchParams, simulate_switched

RES = Path(__file__).resolve().parent / "results"
SEED = {"r120": 3026, "r230": 3027}
RATED = {"r120": 2400.0, "r230": 4800.0}
CAL_FRACTION = 1.0 / 8.0          # earlier load event C2: one eighth of the AC2 rating
AGE = dict(Cb=0.85, Rb=2.0)       # aged buffer capacitor (M5)
ALARM_FACTOR = 3.0
WORKERS = 12


def _outer(q, th, dP, **kw):
    try:
        return bool(oc.feasible(q, th, dP, **kw))
    except (RuntimeError, TypeError, ValueError):
        return None


def _inner_plan(q, th, dP):
    orb = nominal_orbit(q)
    if orb is None:
        return None
    r = solve(q, th, dP, "inner", orbit=orb, smooth=True)
    return r if r.get("feasible") else None


def _max_step(q, th, cap, tol):
    """Largest admissible step at phase th by the circuit outer bound (3 post-step cycles)."""
    r = bisect(lambda d: bool(_outer(q, th, d, n_post=3)), 0.0, cap, tol)
    return 0.0 if np.isnan(r) else float(r)


def _track(q, plant, th, dP, plan, k_rec):
    """Track a plan on the plant, with q as the controller's model and 12-bit sensors."""
    s = simulate_switched(q, SwitchParams(adc_bits=tw.ADC_BITS), th, dP, policy="plan",
                          plan=plan, pre_cycles=3, gains=dict(k_rec=2 * np.pi * k_rec),
                          plant=plant)
    m = s["metrics"]
    return dict(held=not m["violation"], t_violation_ms=m["t_violation_ms"],
                Vdc_min=m["Vdc_min"], Vdc_max=m["Vdc_max"], I_pole_peak=m["I_pole_peak"])


def _truth(pa, th, dP, k_rec):
    """Truth label of the WP4 campaign (Table 3) for a plant: own outer bounds and an oracle-
    tracked inner plan."""
    orb = nominal_orbit(pa)
    if orb is None:
        return "unresolved"
    if not feasible(pa, th, dP, "outer", orbit=orb):
        return "infeasible" if _outer(pa, th, dP) is False else "unresolved"
    plan = solve(pa, th, dP, "inner", orbit=orb, smooth=True)
    if not plan["feasible"]:
        return "unresolved"
    s = simulate_switched(pa, SwitchParams(), th, dP, policy="plan", plan=plan,
                          gains=dict(k_rec=2 * np.pi * k_rec))
    return "unresolved" if s["metrics"]["violation"] else "feasible"


def _verdicts(q, th, dP):
    return dict(circuit_outer=_outer(q, th, dP), inner=_inner_plan(q, th, dP) is not None)


def _case(args):
    rating, p, c, cap, k_rec = args
    t_case = time.time()
    rated = RATED[rating]
    th, dP = np.deg2rad(c["theta_deg"]), c["dP_W"]
    rng = np.random.default_rng(SEED[rating] * 1000 + c["index"])
    pp = robust.perturb(p, **c["factors"])
    out = dict(index=c["index"])
    try:
        # ---- calibration on the plant (C1 steady, C2 an earlier load event)
        gains = tw.sensor_gains(rng)
        th_cal = rng.uniform(0.0, 2 * np.pi)
        dP_cal = CAL_FRACTION * rated
        T, info, (m1, m2), _ = tw.calibrate(p, pp, th_cal, dP_cal, gains, rng)
        t0 = time.time()
        tw.identify(p, [m1, m2], dP_cal)
        out["t_identify_s"] = time.time() - t0
        out["calibration"] = info
        out["param_error"] = {k: float(getattr(T, k) / getattr(pp, k) - 1.0) for k in tw.USED}
        out["reported_error"] = {k: float(v / getattr(pp, k) - 1.0)
                                 for k, v in info["reported"].items()}
        NE = replace(p, E1=T.E1, E2=T.E2)
        preds = dict(NE=NE, T=T)

        # ---- M2: boundary at the case's phase
        tol = 0.005 * rated
        out["max_step_W"] = {name: _max_step(q, th, cap, tol)
                             for name, q in (("true", pp), ("N", p), ("NE", NE), ("T", T))}

        # ---- M3: verdicts on the unseen case (8 post-step cycles, as in Table 3)
        t0 = time.time()
        out["verdicts"] = {name: _verdicts(q, th, dP) for name, q in preds.items()}
        out["t_twin_verdicts_s"] = time.time() - t0

        # ---- M4: closed loop on the plant, for each predictor that admits the case
        out["closed_loop"] = {}
        for name, q in preds.items():
            plan = _inner_plan(q, th, dP)
            out["closed_loop"][name] = (None if plan is None
                                        else _track(q, pp, th, dP, plan, k_rec))

        # ---- M5: drift (aged buffer capacitor after calibration)
        pa = replace(pp, Cb=AGE["Cb"] * pp.Cb, Rb=AGE["Rb"] * pp.Rb)
        pc = replace(p, E1=pp.E1, E2=pp.E2)
        sw = SwitchParams(adc_bits=tw.ADC_BITS)
        th_mon = rng.uniform(0.0, 2 * np.pi)
        ru = tw.record(pc, pp, sw, th_mon, 0.0, pre_cycles=3, n_cycles=2)
        ra = tw.record(pc, pa, sw, th_mon, 0.0, pre_cycles=3, n_cycles=2)
        thr = ALARM_FACTOR * info["residual_cal_V"]
        s_u = tw.parity_residual(T, tw.measure(ru, p, pp, gains, rng))
        s_a = tw.parity_residual(T, tw.measure(ra, p, pa, gains, rng))
        Ta, info_a, _, _ = tw.calibrate(p, pa, rng.uniform(0.0, 2 * np.pi), dP_cal, gains, rng)
        out["drift"] = dict(threshold_V=thr, residual_unaged_V=s_u, residual_aged_V=s_a,
                            alarm_unaged=bool(s_u > thr), alarm_aged=bool(s_a > thr),
                            aged_runs_held=bool(not ra["metrics"]["violation"]),
                            truth=_truth(pa, th, dP, k_rec),
                            verdicts=dict(N=_verdicts(p, th, dP), stale=_verdicts(T, th, dP),
                                          recal=_verdicts(Ta, th, dP)),
                            recal_param_error={k: float(getattr(Ta, k) / getattr(pa, k) - 1.0)
                                               for k in tw.USED})
    except RuntimeError as e:
        out["error"] = str(e)
    out["t_case_s"] = time.time() - t_case
    return out


def main():
    ap = argparse.ArgumentParser(description="WP9 digital-twin validation")
    ap.add_argument("--rating", choices=("r120", "r230"), required=True)
    ap.add_argument("--n", type=int, default=None, help="first n cases only (testing)")
    args = ap.parse_args()
    rating = args.rating
    p = ratings.get(rating)
    W, _, _, _ = _wp1(rating)
    cap = W["sweep"]["search_cap_W"]
    k_rec = W["timedomain_tuning"]["k_rec_hz_chosen"]
    M = json.loads((RES / rating / "wp4_mc.json").read_text())
    cases = [c for c in M["cases"] if "error" not in c][:args.n]
    workers = max(1, min(WORKERS, (os.cpu_count() or 2) - 2))
    t0 = time.time()
    res = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for r in ex.map(_case, [(rating, p, c, cap, k_rec) for c in cases], chunksize=1):
            res.append(r)
            if len(res) % 25 == 0:
                print(f"  {rating}: {len(res)} cases, {time.time() - t0:.0f} s", flush=True)
    R = dict(rating=rating, n=len(res), seed=SEED[rating], cal_fraction=CAL_FRACTION, age=AGE,
             alarm_factor=ALARM_FACTOR, adc_bits=tw.ADC_BITS, noise_lsb=tw.NOISE_LSB,
             gain_tol=tw.GAIN_TOL, window=tw.WINDOW, zc_factor=tw.ZC_FACTOR, used=tw.USED,
             cap_W=cap, k_rec_hz=k_rec, cases=res)
    name = "wp9_twin.json" if args.n is None else "wp9_twin_test.json"
    path = RES / rating / name
    path.write_text(json.dumps(R, indent=1, default=_json_default))
    print(f"wrote {path} in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
