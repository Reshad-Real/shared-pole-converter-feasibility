"""Outer (necessary) feasibility bound for the circuit with its buffer inductor (WP4 follow-up).

The reduced LP in feasibility.py drops the buffer inductor, so it requires the buffer-capacitor
voltage vb itself to stay in the pole-N band. In the circuit only the pole voltage vN must stay
there; the inductor carries vN - vb with its volt-seconds, and stores energy. WP4 found three
switched runs that used this to succeed where the reduced outer LP says no trajectory exists.

This LP keeps the buffer branch in its own linear coordinates and the bus in energy:

    [vb, ib](k+1) = Phi [vb, ib](k) + Gamma vN(k)        exact for vN constant over a sample
    lo(k) - delta <= vN(k) <= tangent_d(Ed(k)) - m(k) + delta
    ib(k) - i1(k) - i2(k) in [-I_pole_max, I_pole_max]   pole-N current limit
    S = Ed + Eb + EL,  S(k+1) = S(k) + Ts (pd(k) - Pdc2 - pac(k) - loss(k))
    Eb >= Cb vb^2 / 2,  EL >= Lb ib^2 / 2                 tangent cuts (necessary)
    Eb <= chord over [0, Vcap],  EL <= chord over the ib range
    0 <= loss(k) <= a bound on resistive, switch and diode losses

The cuts make the stored energies at least what vb and ib require, so the bound is tight when
the limit is an energy deficit (every case met so far) and only valid, not tight, for surplus.
The AC currents are prescribed, as in the reduced LP. Every true trajectory of the averaged
circuit that holds the bus, pole-current and pole-band limits at the samples satisfies these
constraints, apart from the within-sample shape of vN (Ts * w0 is about 0.14 rad).

By default there is no terminal cycle: a run that holds its limits over the horizon need not
be exactly periodic at its end, so the cycle is not a necessary condition.
"""
import numpy as np
import scipy.sparse as sp
from scipy.linalg import expm
from scipy.optimize import linprog

from model import ac_state, band, orbit_energy_at
from feasibility import nominal_orbit

N_CUTS_V = 12          # tangent cuts on Eb >= Cb vb^2 / 2 (24 tightens by under 1 V)
N_CUTS_I = 8           # tangent cuts on EL >= Lb ib^2 / 2
V_CAP_FACTOR = 1.5     # vb range for the chord: [0, 1.5 Vdc_max]


def _zoh(p):
    """Exact sampled dynamics of the Lb-Rb-Cb branch driven by vN."""
    A = np.array([[0.0, 1.0 / p.Cb], [-1.0 / p.Lb, -p.Rb / p.Lb]])
    B = np.array([[0.0], [1.0 / p.Lb]])
    M = np.zeros((3, 3))
    M[:2, :2], M[:2, 2:] = A * p.Ts, B * p.Ts
    E = expm(M)
    return E[:2, :2], E[:2, 2]


def initial_from_orbit(p, theta0, orbit=None):
    """Circuit state at the step on the reduced model's pre-step orbit (ib = pb / vb)."""
    orbit = nominal_orbit(p) if orbit is None else orbit
    t = (theta0 % (2 * np.pi)) / p.w
    Eb0 = orbit_energy_at(orbit, p, t)
    vb0 = np.sqrt(2 * Eb0 / p.Cb)
    pb0 = float(np.interp(np.mod(t, p.T), orbit["t"], orbit["pb"], period=p.T))
    return dict(vb=float(vb0), ib=float(pb0 / vb0), Ed0=float(p.Ed(p.Vdc_nom)), pd0=float(orbit["Pd"]))


def solve(p, theta0, dP, delta=0.0, x0=None, orbit=None, terminal=False, Ron=0.05, Vf=0.8,
          free_start=False, n_post=None):
    """Circuit outer LP for an AC2 step dP at AC1 phase theta0 (rad).

    x0: dict with vb, ib, Ed0, pd0 (for example a run's state at the step); default: the
    pre-step orbit. Ron, Vf: switch on-resistance and diode drop for the loss bound.
    free_start: instead of a given state, start one line cycle before the step from any state
    that repeats over that cycle at constant DC1 power, so the bound covers every pre-step
    operating point. n_post: line cycles after the step (default p.n_cycles).
    """
    free = bool(free_start)
    if not free and x0 is None:
        orbit = nominal_orbit(p) if orbit is None else orbit
        if orbit is None:
            return dict(feasible=False, status="no pre-step orbit")
        x0 = initial_from_orbit(p, theta0, orbit)
    Np, Ts = p.Np, p.Ts
    N0 = Np if free else 0                       # samples before the step
    N = (p.n_cycles if n_post is None else n_post) * Np + N0
    t_step = (theta0 % (2 * np.pi)) / p.w
    tn = t_step - N0 * Ts + Ts * np.arange(N + 1)
    tm = tn[:-1] + 0.5 * Ts
    u1, u2, i1, i2, _ = ac_state(p, tn, t_step, dP)
    pac = ac_state(p, tm, t_step, dP)[4]
    lo, m = band(u1, u2)
    s = i1 + i2
    Imax = p.I_pole_max
    ib_lo, ib_hi = s - Imax, s + Imax
    Vcap = V_CAP_FACTOR * p.Vdc_max
    Ed_lo, Ed_hi = float(p.Ed(p.Vdc_min)), float(p.Ed(p.Vdc_max))
    td1, td0 = 1.0 / (p.Cdc * p.Vdc_nom), 0.5 * p.Vdc_nom
    Phi, Gam = _zoh(p)

    n1 = N + 1
    ivb, iib, iEd, iEb, iEL = (np.arange(n1) + j * n1 for j in range(5))
    ivN, ipd, iloss = (5 * n1 + np.arange(N) + j * N for j in range(3))
    nv = 5 * n1 + 3 * N

    lb = np.full(nv, -np.inf)
    ub = np.full(nv, np.inf)
    lb[ivb], ub[ivb] = 0.0, Vcap
    lb[iib], ub[iib] = ib_lo, ib_hi
    lb[iEd], ub[iEd] = Ed_lo, Ed_hi
    lb[iEb], lb[iEL] = 0.0, 0.0
    lb[ivN], ub[ivN] = np.maximum(lo[:-1] - delta, 0.0), p.Vdc_max + delta
    lb[ipd], ub[ipd] = p.Pdc1_min, p.Pdc1_max
    loss_max = (p.Rb * np.maximum(ib_lo[:-1] ** 2, ib_hi[:-1] ** 2)
                + Ron * (i1[:-1] ** 2 + i2[:-1] ** 2 + Imax ** 2)
                + Vf * (np.abs(i1[:-1]) + np.abs(i2[:-1]) + Imax))
    lb[iloss], ub[iloss] = 0.0, loss_max

    if not free:
        vb0, ib0 = float(x0["vb"]), float(x0["ib"])
        if not (ib_lo[0] - 1e-9 <= ib0 <= ib_hi[0] + 1e-9) or not (Ed_lo <= x0["Ed0"] <= Ed_hi):
            return dict(feasible=False, status="initial state outside limits")
        fix = {ivb[0]: vb0, iib[0]: ib0, iEd[0]: float(x0["Ed0"]),
               iEb[0]: 0.5 * p.Cb * vb0 ** 2, iEL[0]: 0.5 * p.Lb * ib0 ** 2}
        for i, v in fix.items():
            lb[i] = ub[i] = v

    er, ec, ev, eb = [], [], [], []

    def eq(cs, vs, b):
        er.extend([len(eb)] * len(cs))
        ec.extend(cs)
        ev.extend(vs)
        eb.append(b)

    for k in range(N):
        eq([ivb[k + 1], ivb[k], iib[k], ivN[k]], [1.0, -Phi[0, 0], -Phi[0, 1], -Gam[0]], 0.0)
        eq([iib[k + 1], ivb[k], iib[k], ivN[k]], [1.0, -Phi[1, 0], -Phi[1, 1], -Gam[1]], 0.0)
        eq([iEd[k + 1], iEb[k + 1], iEL[k + 1], iEd[k], iEb[k], iEL[k], ipd[k], iloss[k]],
           [1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -Ts, Ts], -Ts * (pac[k] + p.Pdc2))
    if free:                                     # the pre-step cycle repeats at constant DC1 power
        for idx in (ivb, iib, iEd, iEb, iEL):
            eq([idx[N0], idx[0]], [1.0, -1.0], 0.0)
        for k in range(1, N0):
            eq([ipd[k], ipd[k - 1]], [1.0, -1.0], 0.0)
    if terminal:
        P_post = pac[-Np:].mean() + p.Pdc2
        for idx in (ivb, iib, iEd, iEb, iEL):
            eq([idx[N], idx[N - Np]], [1.0, -1.0], 0.0)
        for k in range(N - Np + 1, N):
            eq([ipd[k], ipd[k - 1]], [1.0, -1.0], 0.0)
        lb[ipd[N - 1]] = max(lb[ipd[N - 1]], P_post)

    rows, cols, vals, rhs = [], [], [], []

    def le(cs, vs, b):
        rows.extend([len(rhs)] * len(cs))
        cols.extend(cs)
        vals.extend(vs)
        rhs.append(b)

    a_v = np.linspace(0.0, Vcap, N_CUTS_V)
    for k in range(0 if free else 1, n1):
        for a in a_v:                       # Eb >= Cb (2 a vb - a^2) / 2
            le([iEb[k], ivb[k]], [-1.0, p.Cb * a], 0.5 * p.Cb * a * a)
        le([iEb[k], ivb[k]], [1.0, -0.5 * p.Cb * Vcap], 0.0)          # chord on [0, Vcap]
        for c in np.linspace(ib_lo[k], ib_hi[k], N_CUTS_I):            # EL >= Lb (2 c ib - c^2) / 2
            le([iEL[k], iib[k]], [-1.0, p.Lb * c], 0.5 * p.Lb * c * c)
        le([iEL[k], iib[k]], [1.0, -0.5 * p.Lb * (ib_lo[k] + ib_hi[k])],
           -0.5 * p.Lb * ib_lo[k] * ib_hi[k])                          # chord on the ib range
    for k in range(N):                      # vN <= Vdc - m + delta, with Vdc <= tangent_d(Ed)
        le([ivN[k], iEd[k]], [1.0, -td1], td0 - m[k] + delta)
    rTs = p.slew_dc1 * Ts
    if not free:
        pd0 = float(x0["pd0"])
        le([ipd[0]], [1.0], pd0 + rTs)
        le([ipd[0]], [-1.0], -pd0 + rTs)
    for k in range(1, N):
        le([ipd[k], ipd[k - 1]], [1.0, -1.0], rTs)
        le([ipd[k], ipd[k - 1]], [-1.0, 1.0], rTs)

    A_ub = sp.csr_matrix((vals, (rows, cols)), shape=(len(rhs), nv))
    A_eq = sp.csr_matrix((ev, (er, ec)), shape=(len(eb), nv))
    for method in ("highs", "highs-ipm", "highs-ds"):     # retry when HiGHS gives no verdict
        res = linprog(np.zeros(nv), A_ub=A_ub, b_ub=np.array(rhs), A_eq=A_eq, b_eq=np.array(eb),
                      bounds=np.column_stack([lb, ub]), method=method,
                      options=dict(time_limit=120.0))
        if res.status in (0, 2):
            break
    out = dict(feasible=res.status == 0, status=res.message, code=res.status, t=tn - t_step,
               lo=lo, m=m, delta=delta)
    if res.status == 0:
        x = res.x
        out.update(vb=x[ivb], ib=x[iib], Ed=x[iEd], Eb=x[iEb], EL=x[iEL], vN=x[ivN], pd=x[ipd],
                   loss=x[iloss])
    elif res.status != 2:
        out["warning"] = f"solver status {res.status}: {res.message}"
    return out


def feasible(p, theta0, dP, **kw):
    r = solve(p, theta0, dP, **kw)
    if "warning" in r:
        raise RuntimeError(r["warning"])
    return r["feasible"]


# ---------------------------------------------------------------- one orbit, several events
class _LP:
    """Sparse feasibility LP assembled piece by piece."""

    def __init__(self):
        self.nv, self.lb, self.ub = 0, [], []
        self.E = ([], [], [], [])        # rows, cols, vals, rhs
        self.U = ([], [], [], [])

    def var(self, n, lo, hi):
        idx = self.nv + np.arange(n)
        self.nv += n
        self.lb.append(np.broadcast_to(np.asarray(lo, float), (n,)).copy())
        self.ub.append(np.broadcast_to(np.asarray(hi, float), (n,)).copy())
        return idx

    @staticmethod
    def _add(store, cs, vs, b):
        r, c, v, rhs = store
        r.extend([len(rhs)] * len(cs))
        c.extend(int(x) for x in cs)
        v.extend(vs)
        rhs.append(b)

    def eq(self, cs, vs, b):
        self._add(self.E, cs, vs, b)

    def le(self, cs, vs, b):
        self._add(self.U, cs, vs, b)

    def solve(self, time_limit=300.0):
        lb, ub = np.concatenate(self.lb), np.concatenate(self.ub)
        A_eq = sp.csr_matrix((self.E[2], (self.E[0], self.E[1])), shape=(len(self.E[3]), self.nv))
        A_ub = sp.csr_matrix((self.U[2], (self.U[0], self.U[1])), shape=(len(self.U[3]), self.nv))
        for method in ("highs", "highs-ipm", "highs-ds"):
            res = linprog(np.zeros(self.nv), A_ub=A_ub, b_ub=np.array(self.U[3]), A_eq=A_eq,
                          b_eq=np.array(self.E[3]), bounds=np.column_stack([lb, ub]),
                          method=method, options=dict(time_limit=time_limit))
            if res.status in (0, 2):
                break
        return res


def _segment(lp, p, tn, t_step, dP, delta, Ron, Vf, Phi, Gam):
    """Variables and constraints of one sampled trajectory (states at tn, inputs between)."""
    N = len(tn) - 1
    Ts = p.Ts
    u1, u2, i1, i2, _ = ac_state(p, tn, t_step, dP)
    pac = ac_state(p, tn[:-1] + 0.5 * Ts, t_step, dP)[4]
    lo, m = band(u1, u2)
    s = i1 + i2
    Imax = p.I_pole_max
    ib_lo, ib_hi = s - Imax, s + Imax
    Vcap = V_CAP_FACTOR * p.Vdc_max
    td1, td0 = 1.0 / (p.Cdc * p.Vdc_nom), 0.5 * p.Vdc_nom
    loss_max = (p.Rb * np.maximum(ib_lo[:-1] ** 2, ib_hi[:-1] ** 2)
                + Ron * (i1[:-1] ** 2 + i2[:-1] ** 2 + Imax ** 2)
                + Vf * (np.abs(i1[:-1]) + np.abs(i2[:-1]) + Imax))
    x = dict(vb=lp.var(N + 1, 0.0, Vcap), ib=lp.var(N + 1, ib_lo, ib_hi),
             Ed=lp.var(N + 1, float(p.Ed(p.Vdc_min)), float(p.Ed(p.Vdc_max))),
             Eb=lp.var(N + 1, 0.0, np.inf), EL=lp.var(N + 1, 0.0, np.inf),
             vN=lp.var(N, np.maximum(lo[:-1] - delta, 0.0), p.Vdc_max + delta),
             pd=lp.var(N, p.Pdc1_min, p.Pdc1_max), loss=lp.var(N, 0.0, loss_max))
    vb, ib, Ed, Eb, EL, vN, pd, loss = (x[k] for k in ("vb", "ib", "Ed", "Eb", "EL", "vN", "pd",
                                                        "loss"))
    for k in range(N):
        lp.eq([vb[k + 1], vb[k], ib[k], vN[k]], [1.0, -Phi[0, 0], -Phi[0, 1], -Gam[0]], 0.0)
        lp.eq([ib[k + 1], vb[k], ib[k], vN[k]], [1.0, -Phi[1, 0], -Phi[1, 1], -Gam[1]], 0.0)
        lp.eq([Ed[k + 1], Eb[k + 1], EL[k + 1], Ed[k], Eb[k], EL[k], pd[k], loss[k]],
              [1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -Ts, Ts], -Ts * (pac[k] + p.Pdc2))
        lp.le([vN[k], Ed[k]], [1.0, -td1], td0 - m[k] + delta)
    a_v = np.linspace(0.0, Vcap, N_CUTS_V)
    for k in range(N + 1):
        for a in a_v:
            lp.le([Eb[k], vb[k]], [-1.0, p.Cb * a], 0.5 * p.Cb * a * a)
        lp.le([Eb[k], vb[k]], [1.0, -0.5 * p.Cb * Vcap], 0.0)
        for c in np.linspace(ib_lo[k], ib_hi[k], N_CUTS_I):
            lp.le([EL[k], ib[k]], [-1.0, p.Lb * c], 0.5 * p.Lb * c * c)
        lp.le([EL[k], ib[k]], [1.0, -0.5 * p.Lb * (ib_lo[k] + ib_hi[k])],
              -0.5 * p.Lb * ib_lo[k] * ib_hi[k])
    rTs = p.slew_dc1 * Ts
    for k in range(1, N):
        lp.le([pd[k], pd[k - 1]], [1.0, -1.0], rTs)
        lp.le([pd[k], pd[k - 1]], [-1.0, 1.0], rTs)
    return x


def solve_joint(p, events, n_post=3, delta=0.0, Ron=0.05, Vf=0.8, time_limit=300.0):
    """One pre-step operating orbit shared by several events: is there a periodic orbit from
    which every event (theta_j, dP_j) can be followed for n_post cycles within the limits?

    The pre-step cycle starts at AC1 phase 0, repeats exactly, and keeps DC1 power constant;
    event j starts from the orbit's state at its phase (rounded to the sample grid). Infeasible
    means that no operating point and no controller handles all the events, so the capacitance
    is too small; with one event this is the free-start bound of solve().
    """
    Np, Ts = p.Np, p.Ts
    Phi, Gam = _zoh(p)
    lp = _LP()
    t_pre = Ts * np.arange(Np + 1)
    pre = _segment(lp, p, t_pre, np.inf, 0.0, delta, Ron, Vf, Phi, Gam)
    for key in ("vb", "ib", "Ed", "Eb", "EL"):
        lp.eq([pre[key][Np], pre[key][0]], [1.0, -1.0], 0.0)
    for k in range(1, Np):
        lp.eq([pre["pd"][k], pre["pd"][k - 1]], [1.0, -1.0], 0.0)
    rTs = p.slew_dc1 * Ts
    for theta, dP in events:
        kj = int(round((theta % (2 * np.pi)) / (2 * np.pi) * Np)) % Np
        t0 = kj * Ts
        tn = t0 + Ts * np.arange(n_post * Np + 1)
        br = _segment(lp, p, tn, t0, dP, delta, Ron, Vf, Phi, Gam)
        for key in ("vb", "ib", "Ed", "Eb", "EL"):
            lp.eq([br[key][0], pre[key][kj]], [1.0, -1.0], 0.0)
        lp.le([br["pd"][0], pre["pd"][0]], [1.0, -1.0], rTs)
        lp.le([br["pd"][0], pre["pd"][0]], [-1.0, 1.0], rTs)
    res = lp.solve(time_limit)
    out = dict(feasible=res.status == 0, code=res.status, status=res.message, n_vars=lp.nv)
    if res.status not in (0, 2):
        out["warning"] = f"solver status {res.status}: {res.message}"
    return out


def feasible_joint(p, events, **kw):
    r = solve_joint(p, events, **kw)
    if "warning" in r:
        raise RuntimeError(r["warning"])
    return r["feasible"]
