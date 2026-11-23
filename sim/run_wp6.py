"""WP6 campaign (see WP6_PROTOCOL.md).

    py -u sim/run_wp6.py --part d4 events

d4:     the D4 capacitance: the smallest listed multiple of D2's capacitance at which the full
        load drop is certified at all 24 grid phases (about 15 minutes on 8 workers).
events: every design, rating, event type and phase (about 2 hours).
Results go to results/wp6/.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import os
import time

import numpy as np

import wp6

WORKERS = 8


def _save(name, obj):
    wp6.RES6.mkdir(parents=True, exist_ok=True)
    f = wp6.RES6 / f"wp6_{name}.json"
    f.write_text(json.dumps(obj, indent=1, default=float))
    print(f"wrote {f}", flush=True)


def _cert_job(args):
    rating, cb, ph = args
    import inner_circuit as ic
    from sizing import orbit_at
    i, j = wp6.FULL_DROP
    p, dP = wp6.event_params(rating, cb, i, j)
    orb = orbit_at(p, wp6.SCHEDULE[i])
    if orb is None:
        return False
    return bool(ic.certify(p, np.deg2rad(ph), dP, orb, kappa=wp6.KAPPA, n_cycles=wp6.N_CYC)["feasible"])


def part_d4(ex):
    D = wp6.designs()
    out = {}
    for rating in ("r120", "r230"):
        c2 = D["D2"]["cb"][rating]
        out[rating] = dict(cb=None, tried=[])
        for f in wp6.D4_FACTORS:
            ok = list(ex.map(_cert_job, [(rating, f * c2, ph) for ph in wp6.PHASES_DEG]))
            out[rating]["tried"].append(dict(factor=f, cb=f * c2, certified=ok))
            print(f"  {rating} x{f}: {sum(ok)}/{len(ok)} phases certified", flush=True)
            if all(ok):
                out[rating]["cb"] = f * c2
                break
    _save("d4", out)


def _event_job(args):
    rating, dname, cb, schedule, i, j, ph, cert = args
    try:
        return wp6.evaluate(rating, dname, cb, schedule, i, j, ph, certified=cert)
    except Exception as e:                       # keep the campaign going; report it
        return dict(rating=rating, design=dname, cb=cb, i=i, j=j, phase_deg=ph, error=repr(e))


def part_events(ex):
    D = wp6.designs()
    n = len(wp6.LEVELS)
    jobs = []
    for dname, d in D.items():
        cert = dname in ("D2", "D4")
        for rating, cb in d["cb"].items():
            for i in range(n):
                for j in range(n):
                    if i != j:
                        jobs += [(rating, dname, cb, d["schedule"], i, j, ph, cert)
                                 for ph in wp6.PHASES_DEG]
    jobs.sort(key=lambda a: not a[-1])           # certified (slow) jobs first
    t0, res = time.time(), []
    for r in ex.map(_event_job, jobs, chunksize=1):
        res.append(r)
        if len(res) % 100 == 0:
            print(f"  {len(res)}/{len(jobs)} ({time.time() - t0:.0f} s)", flush=True)
    _save("events", res)


PARTS = dict(d4=part_d4, events=part_events)


def main():
    ap = argparse.ArgumentParser(description="WP6 system case study")
    ap.add_argument("--part", choices=tuple(PARTS), nargs="+", required=True)
    args = ap.parse_args()
    with ProcessPoolExecutor(max_workers=max(1, min(WORKERS, (os.cpu_count() or 2) - 2))) as ex:
        for part in args.part:
            t0 = time.time()
            print(f"== {part}", flush=True)
            PARTS[part](ex)
            print(f"== {part} done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
