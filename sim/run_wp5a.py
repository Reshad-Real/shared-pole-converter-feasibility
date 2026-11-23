"""WP5a sizing campaign. The protocol was fixed in WP5A_PROTOCOL.md before the campaign.

Parts, each writing results/wp5a/wp5a_<part>.json:
  grid     every design point: Cb_conv, Cb_nec (circuit outer), Cb_suff (inner), reduced outer
  extras   base points: monotonicity in Cb, three versus five post-step cycles, and the
           step-down trade-off against the operating point (S5, S6)
  verify   base points in the switched model (S4)

Run from the repository root:  py -u sim/run_wp5a.py --part grid extras verify
"""
import argparse
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path

import numpy as np

from feasibility import solve
from switched import SwitchParams, simulate_switched
from sizing import admits, cb_phase, cb_steady, design
from run_wp4 import _json_default

RES = Path(__file__).resolve().parent / "results" / "wp5a"
RATINGS = ("r120", "r230")
PHASES = np.arange(0, 360, 15)
M_GRID = (0.30, 0.35, 0.40, 0.45, 0.50)
TR_GRID = (0.5, 2.0, 5.0, 10.0, None, 20.0, 30.0)       # None: the bus-loop value
PF_LAG = {"0.9 lagging": float(np.degrees(np.arccos(0.9))),
          "0.8 lagging": float(np.degrees(np.arccos(0.8))),
          "0.8 leading": -float(np.degrees(np.arccos(0.8)))}
WINDOW = (-22.5, -15.0, -7.5, 0.0, 7.5, 15.0, 22.5)
K_REC_HZ = 20
WORKERS = 8


def points():
    base = dict(m=0.40, ramp_ms=None, p2_frac=1.0, pf2_deg=0.0, pf_label="1", alpha=0.5)
    pts = []
    for r in RATINGS:
        for m in M_GRID:
            for tr in TR_GRID:
                pts.append(dict(base, rating=r, m=m, ramp_ms=tr, group="chart"))
        for a in (0.75, 1.0):
            for tr in TR_GRID:
                pts.append(dict(base, rating=r, ramp_ms=tr, alpha=a, group="alpha"))
        for f in (0.5, 0.75):
            pts.append(dict(base, rating=r, p2_frac=f, group="step"))
        for lab, deg in PF_LAG.items():
            pts.append(dict(base, rating=r, pf2_deg=deg, pf_label=lab, group="pf"))
    for i, pt in enumerate(pts):
        pt["id"] = i
        pt["is_base"] = (pt["group"] == "chart" and pt["m"] == 0.40 and pt["ramp_ms"] is None)
    return pts


def _p(pt):
    return design(pt["rating"], pt["m"], pt["ramp_ms"], pt["p2_frac"], pt["pf2_deg"])


def _steady_job(pt):
    p, dP, T_r = _p(pt)
    return cb_steady(p, dP), T_r, p.f, p.Vdc_nom, dP


def _phase_job(args):
    """Smallest Cb for one criterion and phase; NaN if a solver gives no verdict."""
    pt, crit, th = args
    p, dP, _ = _p(pt)
    try:
        return cb_phase(crit, p, float(th), dP, pt["alpha"])
    except RuntimeError as e:
        print(f"no verdict: point {pt['id']} {crit} {th} deg: {e}", flush=True)
        return float("nan")


def _save(name, obj):
    RES.mkdir(parents=True, exist_ok=True)
    (RES / f"wp5a_{name}.json").write_text(json.dumps(obj, indent=1, default=_json_default))
    print("wrote", RES / f"wp5a_{name}.json", flush=True)


def _load(name):
    return json.loads((RES / f"wp5a_{name}.json").read_text())


def part_grid(ex):
    pts = points()
    t0 = time.time()
    st = list(ex.map(_steady_job, pts))
    jobs = [(pt, c, float(th)) for pt in pts for c in ("energy", "outer", "inner") for th in PHASES]
    res = list(ex.map(_phase_job, jobs, chunksize=2))
    print(f"stage 1: {len(jobs)} searches, {time.time() - t0:.0f} s", flush=True)
    n = len(PHASES)
    for i, pt in enumerate(pts):
        cb_st, T_r, f, Vdc, dP = st[i]
        blk = res[i * 3 * n:(i + 1) * 3 * n]
        pt.update(T_r_s=T_r, f=f, Vdc=Vdc, dP_W=dP, cb_steady=cb_st,
                  energy=blk[:n], outer=blk[n:2 * n], inner=blk[2 * n:])
    _save("grid_stage1", pts)

    jobs2 = []
    for pt in pts:
        wo = float(PHASES[int(np.argmax(pt["outer"]))])
        wi = float(PHASES[int(np.argmax(pt["inner"]))])
        pt["circuit_phases"] = [(wo + d) % 360 for d in WINDOW]
        pt["inner_extra_phases"] = [(wi + d) % 360 for d in (-22.5, -7.5, 7.5, 22.5)]
        jobs2 += [(pt, "circuit", th) for th in pt["circuit_phases"]]
        jobs2 += [(pt, "inner", th) for th in pt["inner_extra_phases"]]
    t1 = time.time()
    res2 = iter(ex.map(_phase_job, jobs2, chunksize=1))
    for pt in pts:
        pt["circuit"] = [next(res2) for _ in pt["circuit_phases"]]
        pt["inner_extra"] = [next(res2) for _ in pt["inner_extra_phases"]]
    print(f"stage 2: {len(jobs2)} searches, {time.time() - t1:.0f} s", flush=True)

    for pt in pts:
        pt["n_no_verdict"] = int(sum(np.isnan(v) for k in ("energy", "outer", "inner", "circuit",
                                                           "inner_extra") for v in pt[k]))
        pt["cb_conv"] = max(pt["cb_steady"], float(np.nanmax(pt["energy"])))
        pt["cb_nec"] = float(np.nanmax(pt["circuit"]))
        pt["phase_nec_deg"] = pt["circuit_phases"][int(np.nanargmax(pt["circuit"]))]
        allin = pt["inner"] + pt["inner_extra"]
        allph = list(map(float, PHASES)) + pt["inner_extra_phases"]
        pt["cb_suff"] = float(np.nanmax(allin))
        pt["phase_suff_deg"] = allph[int(np.nanargmax(allin))]
        pt["cb_outer_reduced"] = float(np.nanmax(pt["outer"]))
        pt["penalty_pct"] = 100 * (pt["cb_nec"] / pt["cb_conv"] - 1)
    _save("grid", pts)


def _mono_job(args):
    pt, crit, th, cb_found = args
    p, dP, _ = _p(pt)
    grid = cb_found * np.logspace(np.log10(0.5), np.log10(2.0), 13)
    seq = [admits(crit, replace(p, Cb=float(c)), np.deg2rad(th), dP, pt["alpha"]) for c in grid]
    return dict(criterion=crit, phase_deg=th, cb=grid.tolist(), feasible=seq,
                monotone=all(seq[i] <= seq[i + 1] for i in range(12)))


def _down_job(args):
    pt, alpha, crit, th = args
    p, dP, _ = _p(pt)
    return cb_phase(crit, p, float(th), -p.P2_pre, alpha)


def part_extras(ex):
    pts = _load("grid")
    base = [pt for pt in pts if pt["is_base"]]
    mono = list(ex.map(_mono_job, [
        (pt, c, ph, cb) for pt in base for c, ph, cb in (
            ("energy", float(PHASES[int(np.argmax(pt["energy"]))]), max(pt["energy"])),
            ("inner", pt["phase_suff_deg"], pt["cb_suff"]),
            ("circuit", pt["phase_nec_deg"], pt["cb_nec"]))]))
    five = list(ex.map(_phase_job, [(pt, "circuit5", pt["phase_nec_deg"]) for pt in base]))
    down_jobs = [(pt, a, c, float(th)) for pt in base for a in (0.5, 0.75, 1.0)
                 for c in ("energy", "circuit") for th in PHASES]
    down = list(ex.map(_down_job, down_jobs, chunksize=1))
    out, n = [], len(PHASES)
    for i, pt in enumerate(base):
        rec = dict(rating=pt["rating"], monotonicity=mono[3 * i:3 * i + 3],
                   circuit_3_cycles=pt["cb_nec"], circuit_5_cycles=five[i], tradeoff=[])
        for j, a in enumerate((0.5, 0.75, 1.0)):
            k = (i * 3 + j) * 2 * n
            up = next(q for q in pts if q["rating"] == pt["rating"] and q["m"] == 0.40
                      and q["ramp_ms"] is None and q["alpha"] == a and q["p2_frac"] == 1.0
                      and q["pf2_deg"] == 0.0)
            rec["tradeoff"].append(dict(
                alpha=a, up_cb_nec=up["cb_nec"], up_cb_conv=up["cb_conv"],
                down_energy=max(down[k:k + n]), down_circuit=max(down[k + n:k + 2 * n])))
        out.append(rec)
    _save("extras", out)


def _verify_job(args):
    pt, case, cb, th_deg = args
    p, dP, _ = _p(pt)
    q = replace(p, Cb=cb)
    th = np.deg2rad(th_deg)
    gains = dict(k_rec=2 * np.pi * K_REC_HZ)
    if case == "plan_at_suff":
        plan = solve(q, th, dP, "inner", smooth=True)
        if not plan["feasible"]:
            return dict(case=case, error="inner LP infeasible at Cb_suff")
        s = simulate_switched(q, SwitchParams(), th, dP, policy="plan", plan=plan, gains=gains)
    else:
        s = simulate_switched(q, SwitchParams(), th, dP, gains=gains)
    m = s["metrics"]
    return dict(case=case, rating=pt["rating"], cb=cb, phase_deg=th_deg,
                violation=m["violation"], t_violation_ms=m["t_violation_ms"],
                Vdc_min=m["Vdc_min"], Vdc_max=m["Vdc_max"], I_pole_peak=m["I_pole_peak"])


def part_verify(ex):
    pts = _load("grid")
    jobs = []
    for pt in (q for q in pts if q["is_base"]):
        jobs += [(pt, "plan_at_suff", pt["cb_suff"], pt["phase_suff_deg"]),
                 (pt, "conventional_at_conv", pt["cb_conv"], pt["phase_nec_deg"]),
                 (pt, "conventional_at_0.95nec", 0.95 * pt["cb_nec"], pt["phase_nec_deg"])]
    res = list(ex.map(_verify_job, jobs))
    for r in res:
        print(r, flush=True)
    _save("verify", res)


def _scan_job(args):
    """Post hoc: track the plain inner-LP plan at k x Cb_suff in the switched model."""
    pt, k = args
    p, dP, _ = _p(pt)
    q = replace(p, Cb=k * pt["cb_suff"])
    th = np.deg2rad(pt["phase_suff_deg"])
    from sizing import orbit_at
    from feasibility import pole_currents
    plan = solve(q, th, dP, "inner", orbit=orbit_at(q, 0.5), smooth=True)
    pc = pole_currents(q, plan)
    s = simulate_switched(q, SwitchParams(), th, dP, policy="plan", plan=plan,
                          gains=dict(k_rec=2 * np.pi * K_REC_HZ))
    m = s["metrics"]
    return dict(rating=pt["rating"], k=k, cb=q.Cb, plan_peak_iN=float(np.max(np.abs(pc["iN"]))),
                violation=m["violation"], t_violation_ms=m["t_violation_ms"],
                bus_violation=m["bus_violation"], current_violation=m["current_violation"],
                I_pole_peak=m["I_pole_peak"], I_pole_max=p.I_pole_max)


def part_suff_scan(ex):
    pts = _load("grid")
    jobs = [(pt, k) for pt in pts if pt["is_base"] for k in (1.0, 1.25, 1.5, 2.0, 3.0)]
    res = list(ex.map(_scan_job, jobs))
    for r in res:
        print(r, flush=True)
    _save("suff_scan", res)


def part_drop_reduced(ex):
    """Post hoc: the load drop to zero with the reduced outer LP, which keeps the buffer voltage
    inside the band and is therefore tight for surplus (it ignores the inductor)."""
    pts = _load("grid")
    base = [pt for pt in pts if pt["is_base"]]
    jobs = [(pt, a, "outer", float(th)) for pt in base for a in (0.5, 0.75, 1.0) for th in PHASES]
    res = list(ex.map(_down_job, jobs, chunksize=1))
    n, out = len(PHASES), []
    for i, pt in enumerate(base):
        for j, a in enumerate((0.5, 0.75, 1.0)):
            blk = res[(i * 3 + j) * n:(i * 3 + j + 1) * n]
            out.append(dict(rating=pt["rating"], alpha=a, down_reduced_outer=float(np.nanmax(blk)),
                            phase_deg=float(PHASES[int(np.nanargmax(blk))])))
    _save("drop_reduced", out)


KAPPA = 0.97          # current-aware sufficient bound: pole-N current within 97 % of its limit
N_CYC_SUFF = 4        # line cycles after the step for the current-aware sufficient bound


def _certify(pt, cb, th_deg, alpha=0.5):
    import inner_circuit as ic
    from sizing import orbit_at
    p, dP, _ = _p(pt)
    q = replace(p, Cb=cb)
    orb = orbit_at(q, alpha)
    if orb is None:
        return dict(feasible=False, status="no pre-step orbit")
    try:
        return ic.certify(q, np.deg2rad(th_deg), dP, orb, kappa=KAPPA, n_cycles=N_CYC_SUFF)
    except RuntimeError as e:
        return dict(feasible=False, status=f"solver: {e}")


def _suff_point(pt, alpha=0.5):
    """Post hoc: smallest Cb the current-aware sufficient bound certifies at the worst phase
    (geometric search between max(Cb_nec, steady minimum with the inner headroom) and three
    times that), then certified at eight more phases 45 degrees apart, raising Cb by 10 %
    while any fails. The floor makes the plan end in a steady state a controller holding the
    bus at nominal can keep; without the headroom the post-step orbit sits on the band edge."""
    t0 = time.time()
    th_w = pt["phase_nec_deg"]
    p, dP, _ = _p(pt)
    cb_st_h = cb_steady(p, dP, p.delta_in) if np.isfinite(pt["cb_steady"]) else float("inf")
    base = max(pt["cb_nec"], cb_st_h)
    calls = 0
    out = dict(id=pt["id"], rating=pt["rating"], m=pt["m"], ramp_ms=pt["ramp_ms"], worst_deg=th_w,
                alpha=alpha, cb_steady_h=cb_st_h)
    if not np.isfinite(base):
        return dict(out, cb_suff_c=float("inf"), note="no design")
    hi = 3.0 * base
    calls += 1
    if not _certify(pt, hi, th_w, alpha)["feasible"]:
        return dict(out, cb_suff_c=float("inf"), note="not certified at 3 x base", calls=calls,
                    seconds=time.time() - t0)
    lo = base
    calls += 1
    if _certify(pt, lo, th_w, alpha)["feasible"]:
        hi = lo
    else:
        while hi / lo > 1.04:
            mid = np.sqrt(lo * hi)
            calls += 1
            lo, hi = (lo, mid) if _certify(pt, mid, th_w, alpha)["feasible"] else (mid, hi)
    cb = hi
    others = [(th_w + d) % 360 for d in range(45, 360, 45)]
    ok, raises = False, 0
    for raises in range(9):
        failed = []
        for th in others:
            calls += 1
            if not _certify(pt, cb, th, alpha)["feasible"]:
                failed.append(th)
        if not failed:
            ok = True
            if raises:                                   # the worst phase at the final Cb too
                calls += 1
                ok = bool(_certify(pt, cb, th_w, alpha)["feasible"])
            break
        cb *= 1.10
    return dict(out, cb_suff_c=cb if ok else float("inf"), raises=raises, calls=calls,
                phases_certified=[th_w] + others if ok else [],
                note="" if ok else "not certified at every phase", seconds=time.time() - t0)


def part_suff_circuit(ex):
    pts = [pt for pt in _load("grid") if pt["group"] == "chart" and pt["m"] <= 0.45]
    t0 = time.time()
    res = []
    for r in ex.map(_suff_point, pts, chunksize=1):
        res.append(r)
        print(f"  {len(res)}/{len(pts)} {r['rating']} m {r['m']} ramp {r['ramp_ms']}: "
              f"{r['cb_suff_c'] * 1e6 if np.isfinite(r['cb_suff_c']) else r['cb_suff_c']} "
              f"({r.get('calls')} calls, {time.time() - t0:.0f} s)", flush=True)
    _save("suff_circuit", res)


# Trackers for the certified plans in the switched model, as (policy, extra gains). All follow
# the plan's DC1 power and stored energy. "smooth" (the one S4(a) is judged with) follows a
# smoothed copy of the plan's bus and compares the buffer current with its reference at the
# measurement instant; "flat" holds the bus at nominal; "energy" follows the plan's bus sample
# by sample. The last two are the earlier trackers, kept for comparison.
TRACKERS = dict(smooth=("plan_smooth", dict(sync_ib=True, bus_smooth_n=15)),
                flat=("plan_flat", {}),
                energy=("plan", {}))


def _track_job(args):
    """Certify the step at the worst phase and track the plan with each tracker. Returns one
    dict per tracker."""
    pt, s, alpha = args
    import inner_circuit as ic
    from sizing import orbit_at
    p, dP, _ = _p(pt)
    q = replace(p, Cb=s["cb_suff_c"])
    th = np.deg2rad(s["worst_deg"])
    orb = orbit_at(q, alpha)
    head = dict(id=pt["id"], rating=pt["rating"], m=pt["m"], ramp_ms=pt["ramp_ms"], alpha=alpha,
                cb=q.Cb, phase_deg=s["worst_deg"])
    cert = ic.certify(q, th, dP, orb, kappa=KAPPA, n_cycles=N_CYC_SUFF)
    if not cert["feasible"]:
        return [dict(head, tracker=name, error="not certified on rerun") for name in TRACKERS]
    plan = ic.extend_plan(cert["plan"], q, q.n_cycles)
    peak = float(np.max(np.abs(plan["ib"] - _s_of(q, th, dP, plan["t"]))))
    raw = cert["plan"]
    head.update(plan_vb_min=float(raw["vb"].min()),
                plan_dc1_at_limit=float(np.mean(np.abs(raw["pd"]) > 0.99 * q.Pdc1_max)))
    out = []
    for name, (pol, extra) in TRACKERS.items():
        sim = simulate_switched(q, SwitchParams(), th, dP, policy=pol, plan=plan,
                                gains=dict(k_rec=2 * np.pi * K_REC_HZ, **extra),
                                pre_orbit=None if alpha == 0.5 else orb)
        m = sim["metrics"]
        out.append(dict(head, tracker=name, violation=m["violation"],
                        t_violation_ms=m["t_violation_ms"], bus_violation=m["bus_violation"],
                        current_violation=m["current_violation"], Vdc_min=m["Vdc_min"],
                        Vdc_max=m["Vdc_max"], Vdc_lim=(q.Vdc_min, q.Vdc_max),
                        I_pole_peak=m["I_pole_peak"], I_pole_max=q.I_pole_max,
                        sat_time_ms=m["sat_time_ms"], plan_peak_iN=peak))
    return out


def _s_of(p, th, dP, t_rel):
    from model import ac_state
    t_step = (th % (2 * np.pi)) / p.w
    _, _, i1, i2, _ = ac_state(p, t_rel + t_step, t_step, dP)
    return i1 + i2


def part_verify2(ex):
    """S4(a) with the current-aware bound at the centred operating point: per rating, the
    certified point at the bus-loop ramp with the largest m (the base design if certified)."""
    pts = {pt["id"]: pt for pt in _load("grid")}
    sc = _load("suff_circuit")
    jobs = []
    for r in RATINGS:
        cand = [s for s in sc if s["rating"] == r and s["ramp_ms"] is None
                and np.isfinite(s["cb_suff_c"])]
        if cand:
            s = max(cand, key=lambda s: s["m"])
            jobs.append((pts[s["id"]], s, 0.5))
    res = [r for rs in ex.map(_track_job, jobs) for r in rs]
    for r in res:
        print(r, flush=True)
    _save("verify2", res)


def _suff_alpha1(pt):
    return _suff_point(pt, alpha=1.0)


def part_alpha1(ex):
    """Post hoc: the current-aware sufficient value at the base design with alpha = 1, and its
    switched-model check (S4(a) at the top operating point)."""
    pts = [pt for pt in _load("grid") if pt["group"] == "alpha" and pt["alpha"] == 1.0
           and pt["ramp_ms"] is None]
    res = list(ex.map(_suff_alpha1, pts, chunksize=1))
    for r in res:
        print(r, flush=True)
    byid = {pt["id"]: pt for pt in pts}
    ver = [v for vs in ex.map(_track_job, [(byid[r["id"]], r, 1.0) for r in res
                                            if np.isfinite(r["cb_suff_c"])]) for v in vs]
    for v in ver:
        print(v, flush=True)
    _save("alpha1", dict(suff=res, verify=ver))


def _conv_job(args):
    pt, case, cb = args
    p, dP, _ = _p(pt)
    m = simulate_switched(replace(p, Cb=cb), SwitchParams(), np.deg2rad(pt["phase_nec_deg"]), dP,
                          gains=dict(k_rec=2 * np.pi * K_REC_HZ, sync_ib=True))["metrics"]
    return dict(case=case, rating=pt["rating"], cb=cb, phase_deg=pt["phase_nec_deg"],
                violation=m["violation"], t_violation_ms=m["t_violation_ms"])


def part_verify_conv(ex):
    """S4(b) and (c) again with the buffer-loop timing fix of the final tracker: the
    conventional controller must still trip at Cb_conv and at 0.95 Cb_nec."""
    jobs = []
    for pt in (q for q in _load("grid") if q["is_base"]):
        jobs += [(pt, "conventional_at_conv", pt["cb_conv"]),
                 (pt, "conventional_at_0.95nec", 0.95 * pt["cb_nec"])]
    res = list(ex.map(_conv_job, jobs))
    for r in res:
        print(r, flush=True)
    _save("verify_conv", res)


def part_verify_chart(ex):
    """S4(a) across the main chart: every certified point, centred operating point, worst
    phase, each tracker."""
    pts = {pt["id"]: pt for pt in _load("grid")}
    jobs = [(pts[s["id"]], s, 0.5) for s in _load("suff_circuit") if np.isfinite(s["cb_suff_c"])]
    res = []
    for rs in ex.map(_track_job, jobs, chunksize=1):
        res += rs
        for r in rs:
            print(f"  {r['rating']} m {r['m']} ramp {r['ramp_ms']} {r['tracker']}: "
                  + (r["error"] if "error" in r else
                     f"{'held' if not r['violation'] else 'TRIP %.1f ms' % r['t_violation_ms']} "
                     f"(bus {r['Vdc_min']:.1f}-{r['Vdc_max']:.1f}, I {r['I_pole_peak']:.1f})"),
                  flush=True)
    _save("verify_chart", res)


PARTS = dict(grid=part_grid, extras=part_extras, verify=part_verify, suff_scan=part_suff_scan,
             drop_reduced=part_drop_reduced, suff_circuit=part_suff_circuit, verify2=part_verify2,
             alpha1=part_alpha1, verify_chart=part_verify_chart, verify_conv=part_verify_conv)


def main():
    ap = argparse.ArgumentParser(description="WP5a sizing campaign")
    ap.add_argument("--part", choices=tuple(PARTS), nargs="+", required=True)
    args = ap.parse_args()
    workers = max(1, min(WORKERS, (os.cpu_count() or 2) - 2))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for part in args.part:
            t0 = time.time()
            print(f"== {part}", flush=True)
            PARTS[part](ex)
            print(f"== {part} done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
