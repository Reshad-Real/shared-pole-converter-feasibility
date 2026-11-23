"""Tolerance set and non-ideal operating conditions of WP4 (see WP4_CRITERIA.md).

A plant is the nominal rating with each uncertain parameter scaled by a factor:
  corners: the 16 corners of (Cb, Cdc, grid amplitude, DC1 slew)
  draws:   random factors over all eight parameters, fixed seed
Each perturbed plant is analysed with its own LP from its own centred orbit; the designer's
checks keep the nominal values.
"""
import itertools
from dataclasses import replace

import numpy as np

from switched import SwitchParams

CORNER_KEYS = ("Cb", "Cdc", "E", "slew")
TOL = dict(Cb=0.20, Cdc=0.20, E=0.10, slew=0.30, L=0.10)
RB_RANGE = (0.5, 2.0)
SEED = 4
N_DRAWS_R1 = 64
N_DRAWS_SWEEP = 8

# Harmonic sets: (order, relative amplitude, phase). "peak" adds at the fundamental's peak,
# "flat" opposes it.
HARM_PEAK = ((3, 0.03, -np.pi), (5, 0.03, 0.0), (7, 0.02, -np.pi))
HARM_FLAT = ((3, 0.03, 0.0), (5, 0.03, -np.pi), (7, 0.02, 0.0))
HARM_LOAD = ((3, 0.30, -np.pi), (5, 0.15, 0.0))  # peaky AC2 load current (rectifier)


def perturb(p, Cb=1.0, Cdc=1.0, E=1.0, slew=1.0, L1=1.0, L2=1.0, Lb=1.0, Rb=1.0):
    """Plant with each parameter scaled by its factor (E scales both port amplitudes)."""
    return replace(p, Cb=p.Cb * Cb, Cdc=p.Cdc * Cdc, E1=p.E1 * E, E2=p.E2 * E,
                   slew_dc1=p.slew_dc1 * slew, L1=p.L1 * L1, L2=p.L2 * L2, Lb=p.Lb * Lb,
                   Rb=p.Rb * Rb)


def corners():
    """Factor sets of the 16 corners, in a fixed order."""
    return [{k: 1.0 + s * TOL[k] for k, s in zip(CORNER_KEYS, signs)}
            for signs in itertools.product((-1, 1), repeat=len(CORNER_KEYS))]


def draws(n, seed=SEED):
    """n random factor sets over all eight parameters."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        f = {k: 1.0 + TOL[k] * rng.uniform(-1, 1) for k in CORNER_KEYS}
        f.update({k: 1.0 + TOL["L"] * rng.uniform(-1, 1) for k in ("L1", "L2", "Lb")})
        f["Rb"] = float(np.exp(rng.uniform(*np.log(RB_RANGE))))
        out.append({k: float(v) for k, v in f.items()})
    return out


def label(f):
    """Short text for a factor set, for example 'Cb-20 Cdc+20 E-10 slew+30'."""
    return " ".join(f"{k}{100 * (v - 1):+.0f}" for k, v in f.items() if abs(v - 1) > 1e-9) \
        or "nominal"


# ---------------------------------------------------------------- non-ideal conditions
REF_DELAY_S = 1.5 / 20e3       # nominal sample-to-mid-period delay (one period at 20 kHz)


def loop_scale(sw):
    """Current-loop bandwidth factor that keeps the nominal delay phase lag at crossover."""
    return min(1.0, REF_DELAY_S / ((sw.delay_periods + 0.5) / sw.fsw))


def _switch(**kw):
    sw = SwitchParams(**kw)
    s = loop_scale(sw)
    return replace(sw, wc_ac=sw.wc_ac * s), s


def conditions(p, corner_lo, corner_hi):
    """The R6 conditions. Each entry gives the controller's model `pc`, the plant `pp`
    (None: the same as pc), the switching and sensor model `sw`, a factor on the buffer
    current-loop bandwidth, the parameters the plan is solved with, and the settling time."""
    def freq(x):
        f = p.f * x
        return replace(p, f=f, Ts=1 / (200 * f))

    out = []

    def add(key, label_, group, pc=p, pp=None, plan_p=None, pre_cycles=1, **swkw):
        sw, s = _switch(**swkw)
        out.append(dict(key=key, label=label_, group=group, pc=pc, pp=pp,
                        plan_p=pc if plan_p is None else plan_p, sw=sw, wc_scale=s,
                        pre_cycles=pre_cycles))

    add("reference", "Nominal, dead time 1 µs", "reference")
    add("harm_peak", "Voltage harmonics 4.7 %, peaking", "grid", pc=replace(p, harm_v=HARM_PEAK))
    add("harm_flat", "Voltage harmonics 4.7 %, flat-topped", "grid",
        pc=replace(p, harm_v=HARM_FLAT))
    add("sag", "Sag −10 %", "grid", pc=replace(p, E1=0.9 * p.E1, E2=0.9 * p.E2))
    add("swell", "Swell +10 %", "grid", pc=replace(p, E1=1.1 * p.E1, E2=1.1 * p.E2))
    add("f_low", "Frequency −1 %", "grid", pc=freq(0.99))
    add("f_high", "Frequency +1 %", "grid", pc=freq(1.01))
    add("jump_pos", "Phase jump +20° at the step", "grid", pc=replace(p, jump_deg=20.0),
        plan_p=p)
    add("jump_neg", "Phase jump −20° at the step", "grid", pc=replace(p, jump_deg=-20.0),
        plan_p=p)
    add("pf_lag", "AC2 power factor 0.8 lagging", "grid",
        pc=replace(p, pf2_deg=float(np.degrees(np.arccos(0.8)))))
    add("pf_lead", "AC2 power factor 0.8 leading", "grid",
        pc=replace(p, pf2_deg=-float(np.degrees(np.arccos(0.8)))))
    add("nonlinear", "Nonlinear AC2 load (3rd 30 %, 5th 15 %)", "grid",
        pc=replace(p, harm_i2=HARM_LOAD))
    add("td05", "Dead time 0.5 µs", "switching", td=0.5e-6)
    add("td2", "Dead time 2 µs", "switching", td=2e-6)
    add("fsw10", "Switching frequency 10 kHz", "switching", fsw=10e3)
    add("fsw40", "Switching frequency 40 kHz", "switching", fsw=40e3)
    add("tmin", "Minimum pulse 1 µs", "switching", t_min=1e-6)
    add("adc12", "12-bit ADC on every measurement", "measurement", adc_bits=12)
    add("delay2", "Two-period control delay", "measurement", delay_periods=2)
    add("gainV_lo", "Bus sensor gain −1 %", "measurement", gain_V=0.99)
    add("gainV_hi", "Bus sensor gain +1 %", "measurement", gain_V=1.01)
    add("offset_ib", "Buffer-current sensor offset 0.5 A", "measurement", offset_ib=0.5)
    for key, word, f in (("corner_lo", "least", corner_lo), ("corner_hi", "most", corner_hi)):
        # The controller measures the grid amplitude but keeps the nominal components and DC1
        pp = perturb(p, **f)
        add(key, f"Tolerance corner, {word} favourable ({label(f)})", "tolerance",
            pc=replace(p, E1=pp.E1, E2=pp.E2), pp=pp, pre_cycles=3)
    return out
