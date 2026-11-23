"""Phase sweep for the paper with the circuit-level bounds (WP1 settings, 9 degree grid).

For each rating and phase: the largest AC2 step the circuit outer bound admits (necessary,
eight post-step cycles as in WP4) and the largest the current-aware sufficient bound certifies
(kappa = 0.97, four cycles, as in WP5), both by bisection. The phase-blind criteria come
from the WP1 sweep (results/<rating>/results.json).

    py -u sim/run_paper_sweep.py
"""
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import time

import numpy as np

OUT = Path(__file__).parent / "results" / "paper"
PHASES = tuple(range(0, 360, 9))
KAPPA, N_CYC = 0.97, 4


def _bisect(fn, lo, hi, tol):
    """Largest x in [lo, hi] with fn(x) True, assuming fn is True below and False above."""
    if not fn(lo):
        return float("nan")
    if fn(hi):
        return hi
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if fn(mid) else (lo, mid)
    return lo


def _job(args):
    rating, kind, th_deg, cap = args
    import inner_circuit as ic
    import outer_circuit as oc
    import ratings
    from feasibility import nominal_orbit
    p = ratings.get(rating)
    th = np.deg2rad(th_deg)
    if kind == "outer":
        f = lambda d: bool(oc.feasible(p, th, d))
        tol = 0.002 * cap
    else:
        orb = nominal_orbit(p)
        f = lambda d: bool(ic.certify(p, th, d, orb, kappa=KAPPA, n_cycles=N_CYC)["feasible"])
        tol = 0.01 * cap
    t0 = time.time()
    return dict(rating=rating, kind=kind, phase_deg=th_deg, step_W=_bisect(f, 0.0, cap, tol),
                seconds=time.time() - t0)


def _scan_job(args):
    """Largest certified step on a descending grid below the necessary bound (2.5 % steps).
    Each candidate is certified on its own, so one local miss does not end the search, as it
    does in a bisection."""
    rating, th_deg, outer_W = args
    import inner_circuit as ic
    import ratings
    from feasibility import nominal_orbit
    p = ratings.get(rating)
    orb = nominal_orbit(p)
    th = np.deg2rad(th_deg)
    t0, tried = time.time(), 0
    for k in range(41):
        d = outer_W * (1 - 0.025 * k)
        tried += 1
        if ic.certify(p, th, d, orb, kappa=KAPPA, n_cycles=N_CYC)["feasible"]:
            return dict(rating=rating, kind="suff_scan", phase_deg=th_deg, step_W=d, tried=tried,
                        seconds=time.time() - t0)
    return dict(rating=rating, kind="suff_scan", phase_deg=th_deg, step_W=0.0, tried=tried,
                seconds=time.time() - t0)


def scan():
    f = OUT / "phase_sweep_circuit.json"
    res = json.loads(f.read_text())
    res = [r for r in res if r["kind"] != "suff_scan"]
    jobs = [(r["rating"], r["phase_deg"], r["step_W"]) for r in res if r["kind"] == "outer"]
    with ProcessPoolExecutor(8) as ex:
        for r in ex.map(_scan_job, jobs, chunksize=1):
            res.append(r)
            print(f"  {r['rating']} scan {r['phase_deg']:3d} deg: {r['step_W']:.0f} W "
                  f"({r['tried']} tried, {r['seconds']:.0f} s)", flush=True)
    f.write_text(json.dumps(res, indent=1, default=float))
    print("wrote", f)


def check255():
    """The WP1 pair settings at 240, 255 and 270 degrees with the reduced and the circuit outer
    bound: the averaged runs at 255 degrees hold although the reduced bound calls the step
    infeasible, and the circuit bound (eight cycles) resolves that."""
    import outer_circuit as oc
    import ratings
    from feasibility import feasible
    out = []
    for rating, dP in (("r120", 1000.0), ("r230", 1975.0)):
        p = ratings.get(rating)
        for th in (240, 255, 270):
            t = np.deg2rad(th)
            out.append(dict(rating=rating, phase_deg=th, step_W=dP,
                            reduced_outer=bool(feasible(p, t, dP, "outer")),
                            circuit_outer_8cyc=bool(oc.feasible(p, t, dP))))
            print(out[-1], flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "check_255.json").write_text(json.dumps(out, indent=1))


def main():
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "scan":
        return scan()
    if len(sys.argv) > 1 and sys.argv[1] == "check255":
        return check255()
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = []
    for rating in ("r120", "r230"):
        cap = json.loads((Path(__file__).parent / "results" / rating / "results.json").read_text())[
            "sweep"]["search_cap_W"]
        jobs += [(rating, kind, th, cap) for kind in ("suff", "outer") for th in PHASES]
    res = []
    with ProcessPoolExecutor(8) as ex:
        for r in ex.map(_job, jobs, chunksize=1):
            res.append(r)
            print(f"  {r['rating']} {r['kind']} {r['phase_deg']:3d} deg: {r['step_W']:.0f} W "
                  f"({r['seconds']:.0f} s)", flush=True)
    (OUT / "phase_sweep_circuit.json").write_text(json.dumps(res, indent=1, default=float))
    print("wrote", OUT / "phase_sweep_circuit.json")


if __name__ == "__main__":
    main()
