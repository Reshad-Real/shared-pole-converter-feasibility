"""Buffer sizing for WP5a (see WP5A_PROTOCOL.md).

A design point fixes the rating, modulation index, DC1 ramp time, step size and AC2 power
factor. For each criterion and AC phase, the smallest buffer capacitance that admits the step is
found by geometric bisection, starting from an exact state on the pre-step steady-state orbit at
operating point alpha (0.5 = centred, 1 = highest).
"""
from dataclasses import replace

import numpy as np

from model import steady_orbit
from feasibility import feasible, solve
import outer_circuit as oc
import ratings

CB_LO, CB_HI, RTOL = 20e-6, 20e-3, 0.01
N_POST = 3                      # line cycles after the step for the circuit outer bound


def design(rating, m=0.40, ramp_ms=None, p2_frac=1.0, pf2_deg=0.0):
    """Params (Cb unset) and AC2 step for a design point; also returns the ramp time in s."""
    raw = ratings._raw(rating)
    p, P2r = raw["p"], raw["P2_rated"]
    Vdc = p.E1 / m
    p = replace(p, Vdc_nom=Vdc, Vdc_min=0.95 * Vdc, Vdc_max=1.05 * Vdc,
                Cdc=p.Cdc * (p.Vdc_nom / Vdc) ** 2, delta_in=0.015 * Vdc, pf2_deg=pf2_deg)
    rated_step = P2r - p.P2_pre
    T_r = 1.0 / (np.pi * 0.4 * p.f) if ramp_ms is None else ramp_ms * 1e-3   # bus loop: 2 tau
    p = replace(p, slew_dc1=rated_step / T_r)
    return p, p2_frac * P2r - p.P2_pre, T_r


def orbit_at(p, alpha):
    """Pre-step steady orbit with its mean energy at position alpha of the admissible interval
    (the interval used for the centred orbit; alpha is clipped to [0.001, 0.999] so the start
    is not exactly on the edge)."""
    Ed_lo, Ed_hi = p.Ed(p.Vdc_min), p.Ed(p.Vdc_max)
    V_sec = p.Vdc_min + (p.Vdc_max - p.Vdc_min) * (p.Ed(p.Vdc_nom) - Ed_lo) / (Ed_hi - Ed_lo)
    o = steady_orbit(p, p.P2_pre, V_sec, p.delta_in)
    if o is None:
        return None
    a = min(max(alpha, 0.001), 0.999)
    shift = (a - 0.5) * (o["upper"] - o["lower"])
    return dict(o, Eb=o["Eb"] + shift)


I_LIMIT = 0.9                   # corrected inner bound: pole-N current within 90 % of its limit


def inner_i_delta(p):
    """Band headroom of the corrected inner bound: the inner headroom plus the voltage the
    buffer inductor needs to swing a current of I_LIMIT * I_pole_max at twice line frequency."""
    return p.delta_in + p.Lb * 2 * p.w * I_LIMIT * p.I_pole_max


def admits(criterion, p, theta, dP, alpha):
    """Does `criterion` admit the step at phase theta (rad) with capacitance p.Cb?"""
    orb = orbit_at(p, alpha)
    if orb is None:
        return False
    if criterion == "inner_i":                  # post hoc: current-limited, inductor headroom
        r = solve(p, theta, dP, "inner", delta=inner_i_delta(p), orbit=orb, i_limit=I_LIMIT)
        if "warning" in r:
            raise RuntimeError(r["warning"])
        return bool(r["feasible"])
    if criterion == "circuit":
        return bool(oc.feasible(p, theta, dP, orbit=orb, n_post=N_POST))
    if criterion == "circuit5":
        return bool(oc.feasible(p, theta, dP, orbit=orb, n_post=5))
    return bool(feasible(p, theta, dP, criterion, orbit=orb))


def min_cb(fn, lo=CB_LO, hi=CB_HI, rtol=RTOL):
    """Smallest capacitance for which fn(Cb) holds, assuming it holds above some value."""
    if not fn(hi):
        return float("inf")
    if fn(lo):
        return lo
    while hi / lo > 1 + rtol:
        mid = np.sqrt(lo * hi)
        lo, hi = (lo, mid) if fn(mid) else (mid, hi)
    return hi


def cb_steady(p, dP, delta=0.0):
    """Smallest Cb with a post-step steady orbit inside the band at nominal bus, `delta` volts
    inside each band edge (default: no headroom)."""
    P2 = p.P2_pre + dP
    return min_cb(lambda c: steady_orbit(replace(p, Cb=c), P2, p.Vdc_nom, delta) is not None)


def cb_phase(criterion, p, theta_deg, dP, alpha):
    """Smallest Cb that `criterion` needs at one phase."""
    th = np.deg2rad(theta_deg)
    return min_cb(lambda c: admits(criterion, replace(p, Cb=c), th, dP, alpha))
