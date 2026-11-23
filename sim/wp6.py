"""WP6 system case study (see WP6_PROTOCOL.md).

Event outcomes per design, event type and AC phase, and the daily probabilities of the
stylised backup-load model. An event is an AC2 power step between two load levels, starting
from the design's pre-step steady orbit; the outcomes are the circuit outer bound (can any
controller handle it?), the conventional controller in the switched model, and, for the
proposed design, the certified plan tracked in the switched model.
"""
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

import inner_circuit as ic
import outer_circuit as oc
import ratings
from sizing import design, orbit_at
from switched import SwitchParams, simulate_switched

LEVELS = (0.0, 1 / 3, 2 / 3, 1.0)               # AC2 load, fraction of its rating
PHASES_DEG = tuple(range(0, 360, 15))
SCHEDULE = (1.0, 1.0, 0.5, 0.0)                 # operating point alpha at each load level
K_REC_HZ = 20                                   # conventional mean-energy recovery (as WP5)
KAPPA, N_CYC = 0.97, 4                          # certified controller (as WP5)
TRACKER = ("plan_smooth", dict(sync_ib=True, bus_smooth_n=15))
WEIGHTS = {1: 1.0, 2: 0.25, 3: 0.05}            # next-level weight by jump size (thirds)
LAMBDA = 60.0                                   # events per day
RES5 = Path(__file__).parent / "results" / "wp5a"
RES6 = Path(__file__).parent / "results" / "wp6"
D4_FACTORS = (1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0)   # D4 search, times D2's capacitance
FULL_DROP = (3, 0)


def designs():
    """D1-D4 of the protocol: capacitances from the WP5 results and, for D4, the WP6 search."""
    grid = json.loads((RES5 / "wp5a_grid.json").read_text())
    a1 = json.loads((RES5 / "wp5a_alpha1.json").read_text())
    base = {pt["rating"]: pt for pt in grid if pt["is_base"]}
    suff = {s["rating"]: s["cb_suff_c"] for s in a1["suff"]}
    out = dict(
        D1=dict(name="conventional sizing", cb={r: base[r]["cb_conv"] for r in base}, schedule=None),
        D2=dict(name="proposed sizing", cb=suff, schedule=SCHEDULE),
        D3=dict(name="conventional capacitance, proposed operating point",
                cb={r: base[r]["cb_conv"] for r in base}, schedule=SCHEDULE))
    f = RES6 / "wp6_d4.json"
    if f.exists():
        d4 = json.loads(f.read_text())
        out["D4"] = dict(name="proposed sizing with load drops",
                         cb={r: d4[r]["cb"] for r in d4 if d4[r]["cb"] is not None},
                         schedule=SCHEDULE)
    return out


def event_params(rating, cb, i, j):
    """Params with the pre-step load at level i, and the step to level j."""
    p, _, _ = design(rating)
    P2r = ratings._raw(rating)["P2_rated"]
    return replace(p, Cb=cb, P2_pre=LEVELS[i] * P2r), (LEVELS[j] - LEVELS[i]) * P2r


def alpha(schedule, i):
    return 0.5 if schedule is None else schedule[i]


def evaluate(rating, dname, cb, schedule, i, j, phase_deg, certified=False):
    """All outcomes of one event (see the protocol)."""
    p, dP = event_params(rating, cb, i, j)
    a_pre, a_post = alpha(schedule, i), alpha(schedule, j)
    orb = orbit_at(p, a_pre)
    orb_post = orbit_at(replace(p, P2_pre=p.P2_pre + dP), a_post)
    th = np.deg2rad(phase_deg)
    out = dict(rating=rating, design=dname, cb=cb, i=i, j=j, phase_deg=phase_deg,
               alpha_pre=a_pre, alpha_post=a_post)
    if orb is None or orb_post is None:
        return dict(out, error="no steady orbit at this load")
    out["outer_feasible"] = bool(oc.feasible(p, th, dP, orbit=orb, n_post=3))
    keep = ("violation", "t_violation_ms", "bus_violation", "current_violation",
            "synth_violation", "Vdc_min", "Vdc_max", "I_pole_peak")
    m = simulate_switched(p, SwitchParams(), th, dP, gains=dict(k_rec=2 * np.pi * K_REC_HZ),
                          pre_orbit=None if a_pre == 0.5 else orb,
                          post_orbit=None if a_post == 0.5 else orb_post)["metrics"]
    out["conv"] = {k: m[k] for k in keep}
    if certified:
        cert = ic.certify(p, th, dP, orb, kappa=KAPPA, n_cycles=N_CYC)
        if not cert["feasible"]:
            out["cert"] = dict(certified=False)
        else:
            plan = ic.extend_plan(cert["plan"], p, p.n_cycles)
            pol, extra = TRACKER
            m = simulate_switched(p, SwitchParams(), th, dP, policy=pol, plan=plan,
                                  gains=dict(k_rec=2 * np.pi * K_REC_HZ, **extra),
                                  pre_orbit=None if a_pre == 0.5 else orb)["metrics"]
            out["cert"] = dict(certified=True, **{k: m[k] for k in keep})
    return out


def transition_model(weights=WEIGHTS):
    """Stationary level distribution pi and next-level probabilities P of the stylised chain."""
    n = len(LEVELS)
    W = np.array([[0.0 if a == b else weights.get(abs(a - b), 0.0) for b in range(n)]
                  for a in range(n)])
    P = W / W.sum(axis=1, keepdims=True)
    vals, vecs = np.linalg.eig(P.T)
    pi = np.real(vecs[:, np.argmin(abs(vals - 1))])
    return pi / pi.sum(), P


def daily(p_fail, lam=LAMBDA, weights=WEIGHTS):
    """Expected failing events per day and P_day = 1 - exp(-expected), from p_fail[i][j]."""
    pi, P = transition_model(weights)
    n = len(LEVELS)
    rate = lam * sum(pi[i] * P[i, j] * p_fail[i][j] for i in range(n) for j in range(n) if i != j)
    return rate, 1.0 - np.exp(-rate)
