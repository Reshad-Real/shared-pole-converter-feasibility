"""Parameter sets for the study: the lab-scale baseline and two realistic residential ratings.

Realistic sets (WP1 of simulation_route_plan.md) keep the baseline's modulation index (about
0.4), bus band (+/-5 %), bus-capacitor energy ratio and port-power proportions, and choose:

- r120: 120 V, 60 Hz residential system on a 425 V bus, 2.0 kW AC before the step,
  AC2 (backup output) rated 2.4 kW.
- r230: 230 V, 50 Hz system on an 800 V bus (1200 V SiC), 4.0 kW AC before the step,
  AC2 rated 4.8 kW.

The buffer is sized as a designer would: 1.2 x the smallest capacitance that admits a
steady-state orbit with AC2 at its rating. DC1 is slew-limited with the rate that leaves
the same energy deficit as a first-order response of time constant tau to the rated AC2
step (slew = step / (2 tau)):

- "bus-loop" (default): DC1 follows a DC-bus voltage loop kept below the double-line
  frequency, crossover 24 Hz at 60 Hz and 20 Hz at 50 Hz (fc = f2 / 5);
- "fast": DC1 receives a power feedforward through a 600 Hz current loop.
"""
from dataclasses import replace

import numpy as np

from model import Params, bisect, steady_orbit

RATINGS = ("baseline", "r120", "r230")
F_CURRENT_LOOP = 600.0          # Hz, DC1 current-loop crossover for the "fast" regime


def _raw(name):
    if name == "r120":
        return dict(p=Params(f=60.0, Vdc_nom=425.0, Vdc_min=403.75, Vdc_max=446.25, Cdc=30e-6,
                             E1=170.0, E2=170.0, L1=0.5e-3, R1=0.05, L2=0.5e-3, R2=0.05,
                             P1=1200.0, P2_pre=800.0, Pdc2=400.0, Lb=1e-3, Rb=0.05,
                             Pdc1_min=-6000.0, Pdc1_max=6000.0, I_pole_max=66.0,
                             delta_in=6.4, Ts=1 / (60 * 200)),
                    P2_rated=2400.0)
    if name == "r230":
        return dict(p=Params(f=50.0, Vdc_nom=800.0, Vdc_min=760.0, Vdc_max=840.0, Cdc=20e-6,
                             E1=325.0, E2=325.0, L1=1.5e-3, R1=0.08, L2=1.5e-3, R2=0.08,
                             P1=2400.0, P2_pre=1600.0, Pdc2=800.0, Lb=2e-3, Rb=0.08,
                             Pdc1_min=-12000.0, Pdc1_max=12000.0, I_pole_max=69.0,
                             delta_in=12.0, Ts=1 / (50 * 200)),
                    P2_rated=4800.0)
    raise ValueError(name)


def cb_min_steady(p, P2, V=None, delta=0.0, hi=20e-3):
    """Smallest buffer capacitance with a steady-state orbit at AC2 power P2."""
    V = p.Vdc_nom if V is None else V
    return bisect(lambda c: steady_orbit(replace(p, Cb=c), P2, V, delta) is None,
                  1e-6, hi, 1e-7) + 1e-7


def dc1_slew(p, step_W, regime="bus-loop"):
    """Slew limit (W/s) equivalent to a first-order DC1 response to `step_W`."""
    fc = p.f * 2 / 5 if regime == "bus-loop" else F_CURRENT_LOOP
    tau = 1 / (2 * np.pi * fc)
    return step_W / (2 * tau)


def get(name, regime="bus-loop"):
    """Params for a rating. The baseline is the lab-scale set used so far."""
    if name == "baseline":
        return Params()
    raw = _raw(name)
    p = raw["p"]
    step = raw["P2_rated"] - p.P2_pre
    cb = 1.2 * cb_min_steady(p, raw["P2_rated"])
    return replace(p, Cb=cb, slew_dc1=dc1_slew(p, step, regime))


def describe(name):
    """Summary numbers for reports."""
    if name == "baseline":
        p = Params()
        return dict(name=name, p=p, P2_rated=None)
    raw = _raw(name)
    p = get(name)
    step = raw["P2_rated"] - raw["p"].P2_pre
    return dict(name=name, p=p, P2_rated=raw["P2_rated"], rated_step_W=step,
                slew_bus_W_per_ms=dc1_slew(p, step) / 1e3,
                slew_fast_W_per_ms=dc1_slew(p, step, "fast") / 1e3,
                ramp_bus_ms=1e3 * step / dc1_slew(p, step),
                ramp_fast_ms=1e3 * step / dc1_slew(p, step, "fast"),
                modulation_index=p.E1 / p.Vdc_nom)
