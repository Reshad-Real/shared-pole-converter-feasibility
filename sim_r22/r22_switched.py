"""Event-driven switched model of R22's three-leg decoupling converter (MODEL_SPEC.md).

Legs a, b, c on a DC bus with C_dc and a resistive load. The grid (with L_ac) connects
between legs a and b; the decoupling capacitor C_ac (with L_c) between legs c and b.
Switching edges, dead time and diode freewheeling follow sim/switched.py; the circuit
is integrated exactly between edges with Heun steps.
"""
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "sim"))
from switched import DEAD, HIGH, LOW, PWM, _leg_state, _spec   # noqa: E402,F401


@dataclass(frozen=True)
class R22Params:
    f: float = 50.0
    U_ac: float = 220 * math.sqrt(2)
    L_ac: float = 1.44e-3
    R_ac: float = 0.1
    L_c: float = 1.15e-3
    R_c: float = 0.1
    C_ac: float = 130e-6        # R22 eq. (9); 135 uF is the DC-link value
    C_dc: float = 135e-6
    U_ref: float = 400.0
    fsw: float = 40e3
    td: float = 200e-9
    Ron: float = 0.15
    Vf: float = 3.0
    w_ci: float = 2 * math.pi * 2000.0
    w_cc: float = 2 * math.pi * 2000.0
    w_cv: float = 2 * math.pi * 300.0
    w_dc: float = 2 * math.pi * 10.0
    h_max: float = 2.5e-6

    @property
    def w(self):
        return 2 * math.pi * self.f


def cap_reference(p, I, beta):
    """Capacitor-voltage amplitude and phase that absorb the double-line power of the grid side.

    Grid current i* = I sin(wt + beta); u_ab* = u_ac - (R_ac + jwL_ac) i* (phasors, sine basis).
    The 2w part of x*y with phasors X, Y is (|XY|/2) sin(2wt + angle(XY) - 90 deg).
    """
    X = p.U_ac - (p.R_ac + 1j * p.w * p.L_ac) * I * np.exp(1j * beta)
    Y = I * np.exp(1j * beta)
    A = abs(X * Y) / 2
    alpha = np.angle(X * Y) - np.pi / 2
    V = math.sqrt(2 * A / (p.w * p.C_ac * (1 - p.w ** 2 * p.L_c * p.C_ac)))
    return V, alpha / 2


def zero_sequence(scheme, u, i_meas, U):
    """Zero-sequence signal for SVPWM, DPWM1 or CB-GDPWM (R22 Fig. 2)."""
    if scheme == "SVPWM":
        return -(max(u) + min(u)) / 2
    if scheme == "DPWM1":
        k = int(np.argmax(np.abs(u)))
    elif scheme == "GDPWM":
        m = int(np.argmin(np.abs(u)))
        n = int(np.argmax(np.abs(i_meas)))
        if m != n:
            k = n
        else:
            p_, q_ = [x for x in range(3) if x != n]
            k = p_ if abs(i_meas[p_]) > abs(i_meas[q_]) else q_
    else:
        raise ValueError(scheme)
    return math.copysign(U / 2, u[k]) - u[k]


def initial_state(p, phi, R_load, I_q, apd):
    """Steady-state start: grid current amplitude from the load power (2 % loss allowance)."""
    if R_load is None:
        I_d = 0.3
    else:
        I_d = 2 * (p.U_ref ** 2 / R_load) * 1.02 / (p.U_ac * math.cos(phi))
    I, beta = (math.hypot(I_d, I_q), math.atan2(I_q, I_d)) if R_load is None else (I_d, phi)
    i_a = I * math.sin(beta)
    if apd:
        V, th = cap_reference(p, I, beta)
        u_cac, i_c = V * math.sin(th), -p.C_ac * p.w * V * math.cos(th)
    else:
        u_cac, i_c = 0.0, 0.0
    return I_d, [i_a, i_c, u_cac, p.U_ref]


def simulate(p, scheme="GDPWM", phi_deg=0.0, R_load=80.0, I_q_rms=0.0, apd=True, T_end=0.25,
             T_rec=0.1):
    """Run one case. R_load=None selects STATCOM mode with quadrature current I_q_rms."""
    Tsw = 1.0 / p.fsw
    td, Ron, Vf = p.td, p.Ron, p.Vf
    w = p.w
    phi = math.radians(phi_deg)
    I_q = math.sqrt(2) * I_q_rms
    statcom = R_load is None
    I_d0, x = initial_state(p, phi, R_load, I_q, apd)
    n_per = int(round(T_end / Tsw))
    k_rec0 = n_per - int(round(T_rec / Tsw))
    # DC-voltage PI (plant gain dU/dt per A of grid amplitude, about U_ac/(2 C_dc U))
    g = p.U_ac / (2 * p.C_dc * p.U_ref)
    Kp = p.w_dc / g
    Ki = Kp * p.w_dc / 4
    integ = I_d0
    n_maf = int(round(0.01 / Tsw))
    maf = np.full(n_maf, p.U_ref)
    msum, mpos = float(maf.sum()), 0
    queue = []
    prev = [(LOW, 0.0, 0.0)] * 3
    enabled = [True, True, apd]

    def pole(state, i_out, U):
        if state == HIGH:
            return U - Ron * i_out, 1.0
        if state == LOW:
            return -Ron * i_out, 0.0
        if i_out > 0.0:
            return -Vf - Ron * i_out, 0.0
        return U + Vf - Ron * i_out, 1.0

    def deriv(xx, t, st):
        i_a, i_c, u_cac, U = xx
        io = (-i_a, i_a + i_c, -i_c)
        va, sa = pole(st[0], io[0], U)
        vb, sb = pole(st[1], io[1], U)
        if apd:
            vc, scc = pole(st[2], io[2], U)
        else:
            vc, scc = vb, 0.0
        u_ac = p.U_ac * math.sin(w * t)
        di_a = (u_ac - p.R_ac * i_a - (va - vb)) / p.L_ac
        di_c = (u_cac - p.R_c * i_c - (vc - vb)) / p.L_c if apd else 0.0
        du = -i_c / p.C_ac
        i_dc = sa * io[0] + sb * io[1] + scc * io[2]
        load = 0.0 if statcom else U / R_load
        dU = (-i_dc - load) / p.C_dc
        return (di_a, di_c, du, dU)

    rec = {k: [] for k in ("t", "i_a", "i_c", "u_cac", "U", "U_avg", "dA", "dB", "dC")}
    fine_t, fine_ia = [], []
    for k in range(n_per):
        tk = k * Tsw
        i_a, i_c, u_cac, U = x
        # ---- control (sampled at tk, applied next period)
        msum += U - maf[mpos]; maf[mpos] = U; mpos = (mpos + 1) % n_maf
        U_f = msum / n_maf
        e = p.U_ref - U_f
        integ += Ki * e * Tsw
        I_d = Kp * e + integ
        if statcom:
            I, beta = math.hypot(I_d, I_q), math.atan2(I_q, I_d)
        else:
            I, beta = I_d, phi
        t_app = tk + 1.5 * Tsw
        ia_s = I * math.sin(w * tk + beta)
        ia_r = I * math.sin(w * t_app + beta)
        dia_r = I * w * math.cos(w * t_app + beta)
        u_ab = p.U_ac * math.sin(w * t_app) - p.R_ac * ia_r - p.L_ac * dia_r - p.L_ac * p.w_ci * (ia_s - i_a)
        if apd:
            V, th = cap_reference(p, I, beta)
            uc_s = V * math.sin(w * tk + th)
            ic_ref = -p.C_ac * w * V * math.cos(w * tk + th) - p.C_ac * p.w_cv * (uc_s - u_cac)
            dic_ref = p.C_ac * w * w * V * math.sin(w * t_app + th)
            u_cb = u_cac - p.R_c * i_c - p.L_c * (dic_ref + p.w_cc * (ic_ref - i_c))
            u = [(2 * u_ab - u_cb) / 3, (-u_ab - u_cb) / 3, (-u_ab + 2 * u_cb) / 3]
            u0 = zero_sequence(scheme, u, (i_a, -i_a - i_c, i_c), U)
        else:
            u = [u_ab / 2, -u_ab / 2, 0.0]
            u0 = 0.0
        duties = [min(max(0.5 + (ux + u0) / U, 0.0), 1.0) for ux in u]
        queue.append(duties)
        applied = queue.pop(0) if len(queue) > 1 else duties

        # ---- integrate the period edge to edge
        specs = [_spec(d, tk, Tsw) if en else (LOW, 0.0, 0.0) for d, en in zip(applied, enabled)]
        t_end = tk + Tsw
        cand = {tk + td}
        for sp_ in specs + prev:
            if sp_[0] == PWM:
                cand.update((sp_[1], sp_[1] + td, sp_[2], sp_[2] + td))
        edges = [tk] + sorted(c for c in cand if tk < c < t_end) + [t_end]
        intU = 0.0
        recording = k >= k_rec0
        for a, b in zip(edges[:-1], edges[1:]):
            n_sub = max(1, math.ceil((b - a) / p.h_max - 1e-9))
            h = (b - a) / n_sub
            tm = 0.5 * (a + b)
            st = [_leg_state(prev[j], specs[j], tk, td, tm) for j in range(3)]
            ta = a
            for _ in range(n_sub):
                k1 = deriv(x, ta, st)
                xp = [x[q] + h * k1[q] for q in range(4)]
                k2 = deriv(xp, ta + h, st)
                xn = [x[q] + 0.5 * h * (k1[q] + k2[q]) for q in range(4)]
                intU += 0.5 * h * (x[3] + xn[3])
                x = xn
                ta += h
                if recording:
                    fine_t.append(ta); fine_ia.append(x[0])
        prev = specs
        if recording:
            for key, val in zip(rec, (tk, i_a, i_c, u_cac, U, intU / Tsw, *applied)):
                rec[key].append(val)
    out = {k: np.asarray(v) for k, v in rec.items()}
    out["fine_t"] = np.asarray(fine_t)
    out["fine_ia"] = np.asarray(fine_ia)
    out["meta"] = dict(scheme=scheme, phi_deg=phi_deg, R_load=R_load, I_q_rms=I_q_rms, apd=apd,
                       td_ns=p.td * 1e9, f=p.f, Tsw=Tsw)
    return out
