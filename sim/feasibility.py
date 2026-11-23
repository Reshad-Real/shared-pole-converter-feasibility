"""Transient feasibility tests posed as linear programs in energy coordinates.

Reduced averaged model (Lb -> 0, lossless buffer branch), sampled every Ts from the step:

    Eb[k+1] = Eb[k] + Ts*pb[k]
    Ed[k+1] = Ed[k] + Ts*(pd[k] - Pdc2 - pac[k] - pb[k])

Eb is the buffer-capacitor energy, Ed the bus-capacitor energy, pd the DC1 power (slew
limited) and pac the prescribed AC pole power. The buffer voltage vb = sqrt(2 Eb / Cb) must
stay inside the shared-pole band lo(t) <= vb <= Vdc(t) - m(t).

Because Vdc and vb are concave in energy, the band is bracketed by two linear sets:

  inner  (sufficient):  tangent_b(Eb) <= secant_d(Ed) - m - delta,  Eb >= Eb(lo + delta)
  outer  (necessary):   secant_b(Eb) <= tangent_d(Ed) - m + delta,  Eb >= Eb(lo - delta)

An inner-feasible trajectory satisfies the band (with headroom delta) at every sample.
An infeasible outer LP proves that no trajectory of the reduced model exists, for any
controller, because every true trajectory satisfies the outer constraints. Pole-current
limits are left out of the LP (a relaxation for the outer bound) and checked afterwards.

Two phase-blind baselines use the same dynamics with a fixed energy window:
  energy  : capacitor rating window [0, Vdc_max]  (nominal energy sizing)
  static  : worst-case intersection of the band over the cycle (static derating)

The horizon ends with one line cycle at constant DC1 power and a periodic buffer energy, so
repeating that cycle continues the trajectory forever (a terminal invariant set).
"""
from dataclasses import replace

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

from model import ac_state, band, orbit_energy_at, steady_orbit

MODES = ("inner", "outer", "energy", "static")


def nominal_orbit(p):
    """Pre-step orbit, centred inside the inner set so the start is inner-admissible."""
    Ed_lo, Ed_hi = p.Ed(p.Vdc_min), p.Ed(p.Vdc_max)
    V_sec = p.Vdc_min + (p.Vdc_max - p.Vdc_min) * (p.Ed(p.Vdc_nom) - Ed_lo) / (Ed_hi - Ed_lo)
    return steady_orbit(p, p.P2_pre, V_sec, p.delta_in)


def solve(p, theta0, dP, mode, delta=None, orbit=None, smooth=False, min_violation=False,
          x0=None, i_limit=None):
    """Pose and solve the trajectory LP for an AC2 step dP applied at AC1 phase theta0 (rad).

    smooth: among feasible trajectories, keep the bus near nominal and DC1 smooth.
    min_violation (outer mode only): keep the band relaxed by `delta` as a hard limit and
    minimise the total violation of the unrelaxed band, so violations appear only where needed.
    With smooth as well, the smoothing terms only break ties behind the violation penalty.
    x0: dict overriding the initial state taken from the orbit (Eb0, Ed0, pd0), for example
    the state a simulation had at the step.
    i_limit (inner mode only): keep the pole-N current within i_limit * I_pole_max. With
    ib = pb / vb and vb bounded below by its secant over the inner interval, the rows
    (s - k Imax) vb_low <= pb <= (s + k Imax) vb_low are sufficient for |ib - s| <= k Imax.
    """
    if mode not in MODES:
        raise ValueError(mode)
    if min_violation and mode != "outer":
        raise ValueError("min_violation needs mode='outer'")
    if delta is None:
        delta = p.delta_in if mode == "inner" else 0.0
    x0 = x0 or {}
    if orbit is None:
        orbit = nominal_orbit(p)
    if orbit is None and not {"Eb0", "pd0"} <= set(x0):
        return dict(feasible=False, status="no pre-step orbit")

    N, Np, Ts = p.n_cycles * p.Np, p.Np, p.Ts
    t_step = (theta0 % (2 * np.pi)) / p.w
    tn = t_step + Ts * np.arange(N + 1)
    tm = tn[:-1] + 0.5 * Ts
    u1, u2, _, _, _ = ac_state(p, tn, t_step, dP)
    pac = ac_state(p, tm, t_step, dP)[4]
    lo, m = band(u1, u2)

    Eb0 = x0["Eb0"] if "Eb0" in x0 else orbit_energy_at(orbit, p, t_step)
    Ed_nom = float(p.Ed(p.Vdc_nom))
    Ed0 = float(x0.get("Ed0", Ed_nom))
    pd0 = x0["pd0"] if "pd0" in x0 else orbit["Pd"]
    P_post = pac[-Np:].mean() + p.Pdc2
    Ed_lo, Ed_hi = float(p.Ed(p.Vdc_min)), float(p.Ed(p.Vdc_max))

    iEb = np.arange(N + 1)
    iEd = N + 1 + np.arange(N + 1)
    iPd = 2 * (N + 1) + np.arange(N)
    nv = 3 * N + 2

    # Box bounds
    Eb_lb = np.zeros(N + 1)
    Eb_ub = np.full(N + 1, float(p.Eb(p.Vdc_max + delta)))
    rows, cols, vals, rhs = [], [], [], []      # coupled band rows (A_ub x <= b_ub)

    def add_row(r, cs, vs, b):
        rows.extend([r] * len(cs))
        cols.extend(cs)
        vals.extend(vs)
        rhs.append(b)

    if mode == "inner":
        v_lo = lo + delta
        vs_ = p.Vdc_min - m - delta              # tangent point: the band edge at the lowest bus
        v_tan = p.Vdc_nom - m - delta
        if np.any(vs_ <= v_lo):
            return dict(feasible=False, status="band empty (inner)")
        Eb_lb = p.Eb(v_lo)
        sd1 = (p.Vdc_max - p.Vdc_min) / (Ed_hi - Ed_lo)
        sd0 = p.Vdc_min - sd1 * Ed_lo
        for k in range(N + 1):
            # tangent_b(Eb) = v_tan/2 + Eb/(Cb v_tan)  <=  sd0 + sd1 Ed - m - delta
            add_row(len(rhs), [iEb[k], iEd[k]], [1.0 / (p.Cb * v_tan[k]), -sd1],
                    sd0 - 0.5 * v_tan[k] - m[k] - delta)
        if i_limit is not None:
            # vb >= secant over [v_lo, v_top], with v_top = Vdc_max - m - delta
            v_top = p.Vdc_max - m - delta
            E_lo, E_top = p.Eb(v_lo), p.Eb(v_top)
            B = (v_top - v_lo) / (E_top - E_lo)
            A = v_lo - B * E_lo
            _, _, i1m, i2m, _ = ac_state(p, tm, t_step, dP)
            s = i1m + i2m
            Ik = i_limit * p.I_pole_max
            for k in range(N):
                hi_c, lo_c = Ts * (s[k] + Ik), Ts * (s[k] - Ik)
                for j in (k, k + 1):
                    # Eb[k+1] - Eb[k] <= hi_c (A_j + B_j Eb[j]);  >= lo_c (A_j + B_j Eb[j])
                    cs, vs = [iEb[k + 1], iEb[k]], [1.0, -1.0]
                    vs = [v - (hi_c * B[j] if c == iEb[j] else 0.0) for c, v in zip(cs, vs)]
                    add_row(len(rhs), cs, vs, hi_c * A[j])
                    cs, vs = [iEb[k + 1], iEb[k]], [-1.0, 1.0]
                    vs = [v + (lo_c * B[j] if c == iEb[j] else 0.0) for c, v in zip(cs, vs)]
                    add_row(len(rhs), cs, vs, -lo_c * A[j])
    elif mode == "outer":
        va = np.maximum(lo - delta, 0.0)
        vc = p.Vdc_max - m + delta
        if np.any(vc <= va):
            return dict(feasible=False, status="band empty (outer)")
        Eb_lb = p.Eb(va)
        Eb_ub = p.Eb(vc)
        td1 = 1.0 / (p.Cdc * p.Vdc_nom)
        td0 = 0.5 * p.Vdc_nom
        sb1 = 2.0 / (p.Cb * (vc + va))
        sb0 = va - sb1 * p.Eb(va)
        for k in range(N + 1):
            # secant_b(Eb) = sb0 + sb1 Eb  <=  td0 + td1 Ed - m + delta
            add_row(len(rhs), [iEb[k], iEd[k]], [sb1[k], -td1], td0 - sb0[k] - m[k] + delta)
    elif mode == "energy":
        Eb_lb = np.zeros(N + 1)
        Eb_ub = np.full(N + 1, float(p.Eb(p.Vdc_max)))
    elif mode == "static":
        u1c, u2c, _, _, _ = ac_state(replace(p, P2_pre=p.P2_pre + dP), orbit["t"])
        lo_c, m_c = band(u1c, u2c)
        lo_s = max(lo.max(), lo_c.max())
        hi_s = p.Vdc_nom - max(m.max(), m_c.max())
        if hi_s <= lo_s:
            return dict(feasible=False, status="static window empty")
        Eb_lb = np.full(N + 1, float(p.Eb(lo_s)))
        Eb_ub = np.full(N + 1, float(p.Eb(hi_s)))

    if Eb0 < Eb_lb[0] - 1e-12 or Eb0 > Eb_ub[0] + 1e-12:
        return dict(feasible=False, status="initial state outside set", Eb0=Eb0)

    # Slew rows on DC1 power
    rTs = p.slew_dc1 * Ts
    for k in range(N):
        if k == 0:
            add_row(len(rhs), [iPd[0]], [1.0], pd0 + rTs)
            add_row(len(rhs), [iPd[0]], [-1.0], -pd0 + rTs)
        else:
            add_row(len(rhs), [iPd[k], iPd[k - 1]], [1.0, -1.0], rTs)
            add_row(len(rhs), [iPd[k], iPd[k - 1]], [-1.0, 1.0], rTs)

    # Equalities: initial state, energy balance, terminal cycle
    er, ec, ev, eb = [], [], [], []

    def add_eq(cs, vs, b):
        r = len(eb)
        er.extend([r] * len(cs))
        ec.extend(cs)
        ev.extend(vs)
        eb.append(b)

    add_eq([iEb[0]], [1.0], Eb0)
    add_eq([iEd[0]], [1.0], Ed0)
    for k in range(N):
        add_eq([iEd[k + 1], iEd[k], iEb[k + 1], iEb[k], iPd[k]],
               [1.0, -1.0, 1.0, -1.0, -Ts], -Ts * (pac[k] + p.Pdc2))
    for k in range(N - Np, N):
        add_eq([iPd[k]], [1.0], P_post)
    add_eq([iEb[N], iEb[N - Np]], [1.0, -1.0], 0.0)

    lb = np.concatenate([Eb_lb, np.full(N + 1, Ed_lo), np.full(N, p.Pdc1_min)])
    ub = np.concatenate([Eb_ub, np.full(N + 1, Ed_hi), np.full(N, p.Pdc1_max)])
    c = np.zeros(nv)
    A_ub = sp.csr_matrix((vals, (rows, cols)), shape=(len(rhs), nv))
    A_eq = sp.csr_matrix((ev, (er, ec)), shape=(len(eb), nv))
    b_ub, b_eq = np.array(rhs), np.array(eb)

    if smooth:
        # Keep the bus near nominal and the DC1 command smooth: |Ed - Ed_nom| and |dpd| slacks
        na, nb = N + 1, N
        A_ub = sp.hstack([A_ub, sp.csr_matrix((A_ub.shape[0], na + nb))]).tocsr()
        A_eq = sp.hstack([A_eq, sp.csr_matrix((A_eq.shape[0], na + nb))]).tocsr()
        extra_r, extra_c, extra_v, extra_b = [], [], [], []
        r0 = 0
        for k in range(N + 1):
            a = nv + k
            for s in (1.0, -1.0):
                extra_r += [r0, r0]
                extra_c += [iEd[k], a]
                extra_v += [s, -1.0]
                extra_b.append(s * Ed_nom)
                r0 += 1
        for k in range(1, N):
            b_ = nv + na + k
            for s in (1.0, -1.0):
                extra_r += [r0, r0, r0]
                extra_c += [iPd[k], iPd[k - 1], b_]
                extra_v += [s, -s, -1.0]
                extra_b.append(0.0)
                r0 += 1
        extra = sp.csr_matrix((extra_v, (extra_r, extra_c)), shape=(r0, nv + na + nb))
        A_ub = sp.vstack([A_ub, extra]).tocsr()
        b_ub = np.concatenate([b_ub, extra_b])
        c = np.concatenate([c, np.full(na, 1.0 / Ed_nom), np.full(nb, 1e-3)])
        lb = np.concatenate([lb, np.zeros(na + nb)])
        ub = np.concatenate([ub, np.full(na + nb, np.inf)])

    if min_violation:
        # Slacks in volts below the true lower edge and above the true upper edge
        base = A_ub.shape[1]
        na = 2 * (N + 1)
        A_ub = sp.hstack([A_ub, sp.csr_matrix((A_ub.shape[0], na))]).tocsr()
        A_eq = sp.hstack([A_eq, sp.csr_matrix((A_eq.shape[0], na))]).tocsr()
        er_, ec_, ev_, eb_ = [], [], [], []
        r0 = 0
        for k in range(N + 1):
            s_lo, s_hi = base + k, base + N + 1 + k
            if lo[k] > 0:
                # Eb(lo) - Eb <= Cb * lo * s_lo, so s_lo ~ lo - vb for small violations
                er_ += [r0, r0]
                ec_ += [iEb[k], s_lo]
                ev_ += [-1.0, -p.Cb * lo[k]]
                eb_.append(-float(p.Eb(lo[k])))
                r0 += 1
            er_ += [r0, r0, r0]
            ec_ += [iEb[k], iEd[k], s_hi]
            ev_ += [sb1[k], -td1, -1.0]
            eb_.append(td0 - sb0[k] - m[k])
            r0 += 1
        extra = sp.csr_matrix((ev_, (er_, ec_)), shape=(r0, base + na))
        A_ub = sp.vstack([A_ub, extra]).tocsr()
        b_ub = np.concatenate([b_ub, eb_])
        c = np.concatenate([c, np.full(na, 1e3)])
        lb = np.concatenate([lb, np.zeros(na)])
        ub = np.concatenate([ub, np.full(na, np.inf)])

    for method in ("highs", "highs-ipm", "highs-ds"):     # retry when HiGHS gives no verdict
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                      bounds=np.column_stack([lb, ub]), method=method,
                      options=dict(time_limit=60.0))
        if res.status in (0, 2):
            break
    out = dict(feasible=res.status == 0, status=res.message, code=res.status,
               t=tn - t_step, t_abs=tn, lo=lo, m=m, pac=pac, P_post=P_post, pd0=pd0,
               Eb0=Eb0, mode=mode, delta=delta, theta0=theta0, dP=dP)
    if res.status == 0:
        x = res.x
        out.update(Eb=x[iEb], Ed=x[iEd], pd=x[iPd])
    elif res.status not in (2,):
        out["warning"] = f"solver status {res.status}: {res.message}"
    return out


def feasible(p, theta0, dP, mode, delta=None, orbit=None):
    r = solve(p, theta0, dP, mode, delta=delta, orbit=orbit)
    if "warning" in r:
        raise RuntimeError(r["warning"])
    return r["feasible"]


def post_step_limit(p, V_bus=None, delta=0.0, hi=2000.0, tol=0.5):
    """Largest AC2 step for which a post-step periodic orbit exists (steady-state check)."""
    from model import bisect
    V = p.Vdc_nom if V_bus is None else V_bus
    return bisect(lambda d: steady_orbit(p, p.P2_pre + d, V, delta) is not None, 0.0, hi, tol)


def pole_currents(p, r):
    """Pole currents implied by an LP trajectory (reduced model, ib = pb / vb)."""
    vb = np.sqrt(2 * r["Eb"] / p.Cb)
    pb = np.diff(r["Eb"]) / p.Ts
    vb_mid = 0.5 * (vb[1:] + vb[:-1])
    ib = pb / np.maximum(vb_mid, 1e-3)
    tm = r["t_abs"][:-1] + 0.5 * p.Ts
    t_step = r["t_abs"][0]
    _, _, i1, i2, _ = ac_state(p, tm, t_step, r["dP"])
    return dict(ib=ib, iA=i1, iB=i2, iN=ib - i1 - i2)
