"""Averaged model of the illustrative shared-pole four-port converter (plan section 5).

Poles A, B and N share one DC bus. AC1 sits between poles A and N, AC2 between B and N,
and the ripple buffer (Lb, Rb, Cb) runs from pole N to the negative rail. DC1 feeds the
bus through a slew-limited interface; DC2 is a constant bus load.

Sign convention: positive AC power flows from the converter into the AC port.
"""
from dataclasses import dataclass, replace

import numpy as np


@dataclass(frozen=True)
class Params:
    # Line and bus
    f: float = 50.0
    Vdc_nom: float = 200.0
    Vdc_min: float = 190.0          # predeclared bus limits (+/-5 %)
    Vdc_max: float = 210.0
    Cdc: float = 20e-6              # small bus capacitor (the point of active buffering)
    # AC ports: e = E sin(wt + phi) at the port terminal, unity power factor there
    E1: float = 80.0
    E2: float = 80.0
    phi2: float = 0.0
    L1: float = 2e-3
    R1: float = 0.1
    L2: float = 2e-3
    R2: float = 0.1
    P1: float = 150.0               # AC1 power (grid export)
    P2_pre: float = 100.0           # AC2 power before the step (backup load)
    Pdc2: float = 50.0              # DC2 constant load
    Tr2: float = 1e-3               # AC2 current-amplitude ramp time at the step
    # Buffer
    Cb: float = 250e-6
    Lb: float = 1e-3
    Rb: float = 0.1
    # DC1 interface
    slew_dc1: float = 50e3          # W/s  (50 W/ms)
    Pdc1_min: float = -1000.0
    Pdc1_max: float = 1000.0
    # Limits and discretisation
    I_pole_max: float = 20.0
    delta_in: float = 3.0           # control headroom kept inside the band by the inner bound
    Ts: float = 100e-6
    n_cycles: int = 8
    # Departures from ideal sinusoidal ports (WP4); the defaults give the ideal case
    harm_v: tuple = ()              # ((order, rel. amplitude, phase rad), ...) in both port voltages
    harm_i2: tuple = ()             # ((order, rel. amplitude, phase rad), ...) in the AC2 current
    pf2_deg: float = 0.0            # AC2 current lag behind its voltage (deg); negative = leading
    jump_deg: float = 0.0           # phase jump of both port voltages at the step
    pll_tau: float = 5e-3           # time constant with which the current references follow it

    @property
    def w(self):
        return 2 * np.pi * self.f

    @property
    def T(self):
        return 1.0 / self.f

    @property
    def Np(self):
        return int(round(self.T / self.Ts))

    def Ed(self, V):
        return 0.5 * self.Cdc * np.asarray(V) ** 2

    def Eb(self, v):
        return 0.5 * self.Cb * np.asarray(v) ** 2


def ideal_ports(p):
    return not (p.harm_v or p.harm_i2 or p.pf2_deg or p.jump_deg)


def _jump(p, t, t_step):
    """Voltage phase offset, current-reference phase offset and its time derivative."""
    J = np.deg2rad(p.jump_deg)
    after = t >= t_step
    decay = np.exp(-np.where(after, t - t_step, 0.0) / p.pll_tau)
    return J * after, J * after * (1.0 - decay), J * after * decay / p.pll_tau


def port_voltages(p, t, t_step=np.inf):
    """Port-terminal voltages e1, e2, with harmonics and a phase jump at t_step."""
    t = np.asarray(t, dtype=float)
    th = p.w * t + _jump(p, t, t_step)[0]
    e1, e2 = np.sin(th), np.sin(th + p.phi2)
    for h, a, ph in p.harm_v:
        e1 = e1 + a * np.sin(h * th + ph)
        e2 = e2 + a * np.sin(h * (th + p.phi2) + ph)
    return p.E1 * e1, p.E2 * e2


def port_voltage_fn(p, t_step):
    """Scalar e1(t), e2(t) for the switched plant."""
    import math
    E1, E2, w, ph2 = p.E1, p.E2, p.w, p.phi2
    if ideal_ports(p):
        return lambda t: (E1 * math.sin(w * t), E2 * math.sin(w * t + ph2))
    J, H = math.radians(p.jump_deg), tuple(p.harm_v)

    def ev(t):
        th = w * t + (J if t >= t_step else 0.0)
        e1, e2 = math.sin(th), math.sin(th + ph2)
        for h, a, ph in H:
            e1 += a * math.sin(h * th + ph)
            e2 += a * math.sin(h * (th + ph2) + ph)
        return E1 * e1, E2 * e2
    return ev


def ac_state(p, t, t_step=np.inf, dP=0.0):
    """Converter-side AC voltages u1, u2, port currents i1, i2 and AC pole power u1*i1 + u2*i2.

    AC2 power steps from P2_pre to P2_pre + dP at t_step; its current amplitude ramps smoothly
    over Tr2.
    Currents are prescribed (ideal current loops); u includes the filter L and R drops.
    With the WP4 port departures the voltages carry harmonics and may jump in phase at the
    step; the current references stay sinusoidal on AC1, follow the jump with lag pll_tau, and
    on AC2 may lag their voltage and carry harmonics. P1 and P2 set the fundamental power.
    """
    if not ideal_ports(p):
        return _ac_state_general(p, t, t_step, dP)
    t = np.asarray(t, dtype=float)
    th = p.w * t
    I1 = 2 * p.P1 / p.E1
    I2a = 2 * p.P2_pre / p.E2
    I2b = 2 * (p.P2_pre + dP) / p.E2
    # Raised-cosine amplitude ramp: dA2/dt is continuous, so u2 has no step
    x = np.clip((t - t_step) / p.Tr2, 0.0, 1.0)
    A2 = I2a + (I2b - I2a) * 0.5 * (1.0 - np.cos(np.pi * x))
    dA2 = (I2b - I2a) * 0.5 * np.pi / p.Tr2 * np.sin(np.pi * x)
    s1, c1 = np.sin(th), np.cos(th)
    s2, c2 = np.sin(th + p.phi2), np.cos(th + p.phi2)
    i1 = I1 * s1
    i2 = A2 * s2
    u1 = p.E1 * s1 + p.L1 * I1 * p.w * c1 + p.R1 * i1
    u2 = p.E2 * s2 + p.L2 * (dA2 * s2 + A2 * p.w * c2) + p.R2 * i2
    return u1, u2, i1, i2, u1 * i1 + u2 * i2


def _ac_state_general(p, t, t_step, dP):
    t = np.asarray(t, dtype=float)
    _, Ji, dJi = _jump(p, t, t_step)
    e1, e2 = port_voltages(p, t, t_step)
    thi = p.w * t + Ji                           # phase of the current references
    dth = p.w + dJi
    cpf = np.cos(np.deg2rad(p.pf2_deg))
    I1 = 2 * p.P1 / p.E1
    I2a = 2 * p.P2_pre / (p.E2 * cpf)
    I2b = 2 * (p.P2_pre + dP) / (p.E2 * cpf)
    x = np.clip((t - t_step) / p.Tr2, 0.0, 1.0)
    A2 = I2a + (I2b - I2a) * 0.5 * (1.0 - np.cos(np.pi * x))
    dA2 = (I2b - I2a) * 0.5 * np.pi / p.Tr2 * np.sin(np.pi * x)
    psi = thi + p.phi2 - np.deg2rad(p.pf2_deg)
    g, dg = np.sin(psi), np.cos(psi)
    for h, b, ph in p.harm_i2:
        g = g + b * np.sin(h * psi + ph)
        dg = dg + b * h * np.cos(h * psi + ph)
    i1 = I1 * np.sin(thi)
    i2 = A2 * g
    u1 = e1 + p.L1 * I1 * np.cos(thi) * dth + p.R1 * i1
    u2 = e2 + p.L2 * (dA2 * g + A2 * dg * dth) + p.R2 * i2
    return u1, u2, i1, i2, u1 * i1 + u2 * i2


def band(u1, u2):
    """Pole-N limits implied by 0 <= vA, vB, vN <= Vdc:  lo <= vN <= Vdc - m."""
    lo = np.maximum(0.0, np.maximum(-u1, -u2))
    m = np.maximum(0.0, np.maximum(u1, u2))
    return lo, m


def steady_orbit(p, P2, V_bus, delta, n=4000):
    """Centred periodic buffer-energy orbit for constant port powers and a constant bus.

    Applies the scalar interval test of plan section 6 over one line cycle:
    an initial energy exists iff max(L - W) <= min(U - W). The returned orbit sits at the
    midpoint of that interval. Returns None when no periodic orbit exists.
    """
    q = replace(p, P2_pre=P2)
    t = np.arange(n) * (p.T / n)
    u1, u2, _, _, pac = ac_state(q, t)
    Pd = pac.mean() + p.Pdc2
    pb = Pd - p.Pdc2 - pac                      # bus held constant: buffer takes all ripple
    W = np.concatenate(([0.0], np.cumsum(0.5 * (pb[1:] + pb[:-1]) * (p.T / n))))
    lo, m = band(u1, u2)
    v_lo, v_hi = lo + delta, V_bus - m - delta
    if np.any(v_hi < v_lo):
        return None
    lower = np.max(p.Eb(v_lo) - W)
    upper = np.min(p.Eb(v_hi) - W)
    if lower > upper:
        return None
    Ebar = 0.5 * (lower + upper)
    return dict(t=t, Eb=Ebar + W, Pd=Pd, pb=pb, margin=upper - lower, lower=lower, upper=upper)


def orbit_energy_at(orbit, p, t):
    """Buffer energy of a periodic orbit at absolute time t."""
    return float(np.interp(np.mod(t, p.T), orbit["t"], orbit["Eb"], period=p.T))


def bisect(fn, lo, hi, tol):
    """Largest x in [lo, hi] with fn(x) True, assuming fn is monotone (True below, False above)."""
    if not fn(lo):
        return np.nan
    if fn(hi):
        return hi
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if fn(mid):
            lo = mid
        else:
            hi = mid
    return lo
