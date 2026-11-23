"""Sufficient (inner) feasibility bound for the circuit, with the pole-current limit and the
buffer-inductor voltage (WP5a follow-up).

The reduced inner LP (feasibility.py) has neither, and WP5a's switched runs showed its plans
exceeding the pole-current limit. Here the buffer branch is in its own linear coordinates, as in
outer_circuit.py, and every constraint is made conservative:

    [vb, ib](k+1) = Phi [vb, ib](k) + Gamma vN(k)            exact for vN constant over a sample
    lo + delta <= vN(k) <= Vdc_low - m - delta                at both ends of each sample
    |ib - i1 - i2| <= kappa * I_pole_max                      pole-N current, with headroom
    bus energy bracketed:  Ed_low <= Ed_true <= Ed_high, both inside the bus limits

with S = Ed + Eb + EL advanced twice, once with the largest and once with the smallest possible
losses, and Ed_low = S_low - Eb_upper - EL_upper, Ed_high = S_high - Eb_lower - EL_lower.

Stored energies and losses are convex in vb and ib. Their lower bounds are tangents; their upper
bounds are chords, which are valid only inside an interval, so vb and ib are confined to tubes
of half-width w_v, w_i around a reference trajectory. That restriction is allowed in a
sufficient bound, and narrow tubes keep the chords tight. The reference is refined by sequential
linear programming: an elastic LP minimises the total violation, its solution becomes the next
reference, and the bound certifies the step only when a strict LP (no slack) is feasible.

A certified trajectory holds every limit of the averaged circuit at the samples, with the given
headroom, over the horizon, and ends in a cycle that repeats in vb and ib at constant DC1 power.
The bus energy is held inside its limits for the worst-case losses but is not forced to repeat:
the slow bus loop absorbs the residual loss mismatch.
"""
import numpy as np

from model import ac_state, band, orbit_energy_at
from outer_circuit import _LP, _zoh
import outer_circuit as oc

DT_FRAC = 0.04          # diode conduction fraction: 2 * dead time * fsw (1 us at 20 kHz)
RON, VF = 0.05, 0.8     # switch and diode drops, as in the switched model


def _loss(p, ib, s, i1, i2):
    """Conduction loss (W) of the three poles and the buffer branch, convex in ib."""
    iN = ib - s
    return (p.Rb * ib ** 2 + RON * (i1 ** 2 + i2 ** 2 + iN ** 2)
            + VF * DT_FRAC * (np.abs(i1) + np.abs(i2) + np.abs(iN)))


def _dloss(p, ib, s):
    iN = ib - s
    return 2 * p.Rb * ib + 2 * RON * iN + VF * DT_FRAC * np.sign(iN)


BUS_HEADROOM = 0.02     # bus kept 2 % of nominal inside each limit (1 % proved too thin in WP5a S4)
DIB_MAX = 36e3          # A/s: buffer-current rate a 2 kHz current loop can follow (3 A per 83 us);
                        # the transient needs about 7 A/ms and the 2f ripple about 14 A/ms


def _build(p, ref, x0, pre, t_step, dP, kappa, delta, wv, wi, elastic, smooth=False):
    N = len(ref["vb"]) - 1
    Np, Ts = p.Np, p.Ts
    tn = t_step + Ts * np.arange(N + 1)
    u1, u2, i1, i2, _ = ac_state(p, tn, t_step, dP)
    pac = ac_state(p, tn[:-1] + 0.5 * Ts, t_step, dP)[4]
    lo, m = band(u1, u2)
    s = i1 + i2
    Ik = kappa * p.I_pole_max
    Ed_lo, Ed_hi = float(p.Ed(p.Vdc_min)), float(p.Ed(p.Vdc_max))
    sd1 = (p.Vdc_max - p.Vdc_min) / (Ed_hi - Ed_lo)             # Vdc >= Vdc_min + sd1 (Ed - Ed_lo)
    hb = BUS_HEADROOM * p.Vdc_nom
    Eb_lim_lo, Eb_lim_hi = float(p.Ed(p.Vdc_min + hb)), float(p.Ed(p.Vdc_max - hb))
    Vcap = oc.V_CAP_FACTOR * p.Vdc_max
    Phi, Gam = _zoh(p)

    cv, ci = ref["vb"], ref["ib"]
    va, vbnd = np.maximum(cv - wv, 0.0), np.minimum(cv + wv, Vcap)
    ia, ibnd = ci - wi, ci + wi
    lp = _LP()
    vb = lp.var(N + 1, va, vbnd)
    ib = lp.var(N + 1, ia, ibnd)
    Slo = lp.var(N + 1, -np.inf, np.inf)
    Shi = lp.var(N + 1, -np.inf, np.inf)
    vN = lp.var(N, 0.0, p.Vdc_max)
    pd = lp.var(N, p.Pdc1_min, p.Pdc1_max)
    sig = lp.var(6 * (N + 1), 0.0, np.inf if elastic else 0.0)  # band lo/hi, current hi/lo, bus lo/hi
    sg = sig.reshape(6, N + 1)

    vb0, ib0 = float(x0["vb"]), float(x0["ib"])
    S0 = float(x0["Ed0"]) + 0.5 * p.Cb * vb0 ** 2 + 0.5 * p.Lb * ib0 ** 2
    for idx, val in ((vb[0], vb0), (ib[0], ib0), (Slo[0], S0), (Shi[0], S0)):
        lp.eq([idx], [1.0], val)

    # Linear bounds on Eb(vb), EL(ib) and the loss, per sample: tangent (lower), chord (upper)
    def chord(f, a, b):
        sl = (f(b) - f(a)) / np.maximum(b - a, 1e-9)
        return sl, f(a) - sl * a                                 # f <= sl x + c on [a, b]

    fEb = lambda v: 0.5 * p.Cb * v * v
    fEL = lambda i: 0.5 * p.Lb * i * i
    eb_s, eb_c = chord(fEb, va, vbnd)
    el_s, el_c = chord(fEL, ia, ibnd)
    ebt_s, ebt_c = p.Cb * cv, -0.5 * p.Cb * cv * cv              # tangent at cv
    elt_s, elt_c = p.Lb * ci, -0.5 * p.Lb * ci * ci
    fl = [lambda x, k=k: _loss(p, x, s[k], i1[k], i2[k]) for k in range(N)]
    for k in range(N):
        a, b = ia[k], ibnd[k]
        ls = (fl[k](b) - fl[k](a)) / (b - a)
        lc = fl[k](a) - ls * a
        lt_s = _dloss(p, ci[k], s[k])
        lt_c = fl[k](ci[k]) - lt_s * ci[k]
        lp.eq([vb[k + 1], vb[k], ib[k], vN[k]], [1.0, -Phi[0, 0], -Phi[0, 1], -Gam[0]], 0.0)
        lp.eq([ib[k + 1], vb[k], ib[k], vN[k]], [1.0, -Phi[1, 0], -Phi[1, 1], -Gam[1]], 0.0)
        # S_low uses the largest loss (chord), S_high the smallest (tangent, floored at zero
        # by not letting it exceed... the tangent is a valid lower bound of the convex loss)
        lp.eq([Slo[k + 1], Slo[k], pd[k], ib[k]], [1.0, -1.0, -Ts, Ts * ls], -Ts * (pac[k] + p.Pdc2 + lc))
        lp.eq([Shi[k + 1], Shi[k], pd[k], ib[k]], [1.0, -1.0, -Ts, Ts * lt_s], -Ts * (pac[k] + p.Pdc2 + lt_c))
        for j in (k, k + 1):
            # vN >= lo + delta
            lp.le([vN[k], sg[0, j]], [-1.0, -1.0], -(lo[j] + delta))
            # vN <= Vdc_low - m - delta, Vdc_low = Vdc_min + sd1 (Ed_low - Ed_lo),
            # Ed_low = Slo - (eb_s vb + eb_c) - (el_s ib + el_c)
            lp.le([vN[k], Slo[j], vb[j], ib[j], sg[1, j]],
                  [1.0, -sd1, sd1 * eb_s[j], sd1 * el_s[j], -1.0],
                  p.Vdc_min - sd1 * Ed_lo - sd1 * (eb_c[j] + el_c[j]) - m[j] - delta)
    for k in range(N + 1):
        lp.le([ib[k], sg[2, k]], [1.0, -1.0], s[k] + Ik)          # ib - s <= Ik
        lp.le([ib[k], sg[3, k]], [-1.0, -1.0], -(s[k] - Ik))      # ib - s >= -Ik
        # Ed_low >= Ed(Vdc_min + headroom)
        lp.le([Slo[k], vb[k], ib[k], sg[4, k]], [-1.0, eb_s[k], el_s[k], -1.0],
              -Eb_lim_lo - eb_c[k] - el_c[k])
        # Ed_high = Shi - (ebt_s vb + ebt_c) - (elt_s ib + elt_c) <= Ed(Vdc_max - headroom)
        lp.le([Shi[k], vb[k], ib[k], sg[5, k]], [1.0, -ebt_s[k], -elt_s[k], -1.0],
              Eb_lim_hi + ebt_c[k] + elt_c[k])
    dI = DIB_MAX * Ts
    for k in range(N):                                            # trackable buffer current
        lp.le([ib[k + 1], ib[k]], [1.0, -1.0], dI)
        lp.le([ib[k + 1], ib[k]], [-1.0, 1.0], dI)
    rTs = p.slew_dc1 * Ts
    lp.le([pd[0]], [1.0], pre["pd0"] + rTs)
    lp.le([pd[0]], [-1.0], -pre["pd0"] + rTs)
    for k in range(1, N):
        lp.le([pd[k], pd[k - 1]], [1.0, -1.0], rTs)
        lp.le([pd[k], pd[k - 1]], [-1.0, 1.0], rTs)
    term = lp.var(4, 0.0, np.inf if elastic else 0.0)             # terminal mismatch (elastic)
    for j, idx in enumerate((vb, ib)):                            # terminal cycle repeats
        lp.eq([idx[N], idx[N - Np], term[2 * j], term[2 * j + 1]], [1.0, -1.0, 1.0, -1.0], 0.0)
    for k in range(N - Np + 1, N):
        lp.eq([pd[k], pd[k - 1]], [1.0, -1.0], 0.0)

    # Objective: total violation (normalised), plus a small push towards nominal bus energy
    w = np.concatenate([np.full(N + 1, 1.0 / max(p.Vdc_nom * 0.01, 1.0))] * 2
                       + [np.full(N + 1, 1.0 / max(p.I_pole_max * 0.01, 1e-3))] * 2
                       + [np.full(N + 1, 1.0 / (0.01 * Ed_hi))] * 2)
    if smooth and not elastic:
        # Pick a trackable certified plan: bus near nominal and DC1 smooth. Deviation variables
        # d >= |Ed_mid - Ed_nom| and e >= |pd(k) - pd(k-1)|; feasibility is unchanged.
        Ed_nom = float(p.Ed(p.Vdc_nom))
        dvar = lp.var(N + 1, 0.0, np.inf)
        evar = lp.var(N, 0.0, np.inf)
        for k in range(N + 1):
            # Ed_mid = 0.5 (Slo - eb_s vb - eb_c - el_s ib - el_c) + 0.5 (Shi - ebt_s vb - ebt_c - elt_s ib - elt_c)
            cs = [Slo[k], Shi[k], vb[k], ib[k]]
            vs = [0.5, 0.5, -0.5 * (eb_s[k] + ebt_s[k]), -0.5 * (el_s[k] + elt_s[k])]
            const = -0.5 * (eb_c[k] + el_c[k] + ebt_c[k] + elt_c[k])
            lp.le(cs + [dvar[k]], vs + [-1.0], Ed_nom - const)
            lp.le(cs + [dvar[k]], [-v for v in vs] + [-1.0], -(Ed_nom - const))
        for k in range(1, N):
            lp.le([pd[k], pd[k - 1], evar[k]], [1.0, -1.0, -1.0], 0.0)
            lp.le([pd[k], pd[k - 1], evar[k]], [-1.0, 1.0, -1.0], 0.0)
        fvar = lp.var(N, 0.0, np.inf)              # |ib(k+1) - ib(k)|: a current a loop can follow
        for k in range(N):
            lp.le([ib[k + 1], ib[k], fvar[k]], [1.0, -1.0, -1.0], 0.0)
            lp.le([ib[k + 1], ib[k], fvar[k]], [-1.0, 1.0, -1.0], 0.0)
    c = np.zeros(lp.nv)
    c[sig] = w if elastic else 0.0
    if elastic:
        c[term] = [1.0 / max(p.Vdc_nom * 0.01, 1.0)] * 2 + [1.0 / max(p.I_pole_max * 0.01, 1e-3)] * 2
    if smooth and not elastic:
        c[dvar] = 1.0 / Ed_nom
        c[evar] = 1e-3 / max(p.P2_pre, 1.0)
        c[fvar] = 1.0 / p.I_pole_max
    idx = dict(vb=vb, ib=ib, Slo=Slo, Shi=Shi, vN=vN, pd=pd, sig=sig,
               eb=(eb_s, eb_c), el=(el_s, el_c), ebt=(ebt_s, ebt_c), elt=(elt_s, elt_c))
    return lp, c, idx, tn


def _solve_lp(lp, c, time_limit=120.0):
    import scipy.sparse as sp
    from scipy.optimize import linprog
    lb, ub = np.concatenate(lp.lb), np.concatenate(lp.ub)
    A_eq = sp.csr_matrix((lp.E[2], (lp.E[0], lp.E[1])), shape=(len(lp.E[3]), lp.nv))
    A_ub = sp.csr_matrix((lp.U[2], (lp.U[0], lp.U[1])), shape=(len(lp.U[3]), lp.nv))
    for method in ("highs", "highs-ipm", "highs-ds"):
        res = linprog(c, A_ub=A_ub, b_ub=np.array(lp.U[3]), A_eq=A_eq, b_eq=np.array(lp.E[3]),
                      bounds=np.column_stack([lb, ub]), method=method,
                      options=dict(time_limit=time_limit))
        if res.status in (0, 2):
            break
    return res


TUBES = ((100.0, 40.0), (60.0, 20.0), (30.0, 10.0), (15.0, 5.0), (8.0, 2.5))   # (w_v V, w_i A)


def _smooth(x, n=9):
    """Moving average with the first value kept: a starting reference the rate limit can follow."""
    x = np.asarray(x, float)
    y = np.convolve(np.pad(x, (n // 2, n // 2), mode="edge"), np.ones(n) / n, mode="valid")
    return y - y[0] + x[0]


def extend_plan(plan, p, n_cycles):
    """Continue a certified plan to n_cycles by repeating its last line cycle, which repeats in
    vb and ib at constant DC1 power (the bus energy repeats only approximately)."""
    Np = p.Np
    n_have = len(plan["vb"]) - 1
    n_need = n_cycles * Np
    if n_need <= n_have:
        return plan
    reps = int(np.ceil((n_need - n_have) / Np))
    out = dict(plan)
    for key in ("vb", "ib", "Eb", "Ed", "Ed_low", "Ed_high"):
        x = np.asarray(plan[key])
        out[key] = np.concatenate([x] + [x[-Np:]] * reps)[:n_need + 1]
    for key in ("vN", "pd"):
        x = np.asarray(plan[key])
        out[key] = np.concatenate([x] + [x[-Np:]] * reps)[:n_need]
    out["t"] = plan["t"][0] + p.Ts * np.arange(n_need + 1)
    return out


def certify(p, theta0, dP, orbit, kappa=0.95, delta=None, tubes=TUBES, max_iter=40,
            ref=None, n_cycles=None):
    """Try to certify the step at phase theta0 (rad) from the orbit's state.

    Tubes start wide, for large moves of the reference, and narrow when the violation stops
    falling fast; a certificate at any width is valid (wider is only more conservative).
    Returns dict(feasible, iterations, violation history, plan). The plan has t, Eb, Ed, pd, vb,
    ib for simulate_switched's 'plan' policy (Ed is the middle of the certified bracket).
    """
    delta = p.delta_in if delta is None else delta
    if n_cycles is not None and n_cycles != p.n_cycles:
        from dataclasses import replace
        p = replace(p, n_cycles=n_cycles)
    N = p.n_cycles * p.Np
    t_step = (theta0 % (2 * np.pi)) / p.w
    Eb0 = orbit_energy_at(orbit, p, t_step)
    vb0 = float(np.sqrt(2 * Eb0 / p.Cb))
    pb0 = float(np.interp(np.mod(t_step, p.T), orbit["t"], orbit["pb"], period=p.T))
    x0 = dict(vb=vb0, ib=pb0 / vb0, Ed0=float(p.Ed(p.Vdc_nom)))
    pre = dict(pd0=float(orbit["Pd"]))
    if ref is not None and (len(ref["vb"]) != N + 1 or abs(ref["vb"][0] - vb0) > 1e-6
                            or abs(ref["ib"][0] - x0["ib"]) > 1e-6):
        ref = dict(vb=np.asarray(ref["vb"]) - ref["vb"][0] + vb0,       # warm start: shift
                   ib=np.asarray(ref["ib"]) - ref["ib"][0] + x0["ib"]) if len(ref["vb"]) == N + 1 else None
    if ref is None:                         # start from the circuit outer bound's trajectory
        r = oc.solve(p, theta0, dP, orbit=orbit, n_post=p.n_cycles)
        if r.get("feasible"):
            ref = dict(vb=_smooth(r["vb"]), ib=_smooth(r["ib"]))
        else:
            ref = dict(vb=np.full(N + 1, vb0), ib=np.full(N + 1, x0["ib"]))
    # Fallback reference, smooth by construction: the pre-step orbit continued over the horizon
    tt = t_step + p.Ts * np.arange(N + 1)
    Eo = np.array([orbit_energy_at(orbit, p, t) for t in tt])
    vbo = np.sqrt(2 * Eo / p.Cb)
    ibo = np.interp(np.mod(tt, p.T), orbit["t"], orbit["pb"], period=p.T) / vbo
    fallback = dict(vb=vbo - vbo[0] + vb0, ib=ibo - ibo[0] + x0["ib"])
    hist, level, used_fallback = [], 0, False
    for it in range(max_iter):
        wv, wi = tubes[level]
        lp, c, idx, tn = _build(p, ref, x0, pre, t_step, dP, kappa, delta, wv, wi, elastic=True)
        res = _solve_lp(lp, c)
        if res.status != 0 and not used_fallback:
            used_fallback, ref, level = True, fallback, 0     # restart from the orbit reference
            lp, c, idx, tn = _build(p, ref, x0, pre, t_step, dP, kappa, delta, *tubes[0],
                                    elastic=True)
            res = _solve_lp(lp, c)
        if res.status != 0:
            return dict(feasible=False, iterations=it + 1, history=hist,
                        status=f"elastic LP status {res.status}")
        v = float(res.fun)
        hist.append(v)
        ref = dict(vb=res.x[idx["vb"]], ib=res.x[idx["ib"]])
        if v >= 1e-7 and len(hist) >= 2 and v > 0.6 * hist[-2] and level < len(tubes) - 1:
            level += 1                       # slow progress: narrow the tubes
        if v < 1e-7:
            lp, c, idx, tn = _build(p, ref, x0, pre, t_step, dP, kappa, delta, wv, wi,
                                    elastic=False, smooth=True)
            res = _solve_lp(lp, c)
            if res.status == 0:
                x = res.x
                vbx, ibx = x[idx["vb"]], x[idx["ib"]]
                eb_s, eb_c = idx["eb"]
                el_s, el_c = idx["el"]
                ebt_s, ebt_c = idx["ebt"]
                elt_s, elt_c = idx["elt"]
                Ed_low = x[idx["Slo"]] - (eb_s * vbx + eb_c) - (el_s * ibx + el_c)
                Ed_high = x[idx["Shi"]] - (ebt_s * vbx + ebt_c) - (elt_s * ibx + elt_c)
                plan = dict(t=tn - t_step, vb=vbx, ib=ibx, vN=x[idx["vN"]], pd=x[idx["pd"]],
                            Eb=0.5 * p.Cb * vbx ** 2, Ed=0.5 * (Ed_low + Ed_high),
                            Ed_low=Ed_low, Ed_high=Ed_high)
                return dict(feasible=True, iterations=it + 1, history=hist, plan=plan,
                            tube=(wv, wi))
        if level == len(tubes) - 1 and len(hist) >= 4 and v > 0.8 * hist[-4]:
            break                            # narrowest tubes and no longer improving
    sg = res.x[idx["sig"]].reshape(6, -1)
    names = ("band low (V)", "band high (V)", "current high (A)", "current low (A)",
             "bus low (J)", "bus high (J)")
    where = {n: dict(max=float(sg[i].max()), at_ms=float(1e3 * (tn[int(np.argmax(sg[i]))] - t_step)),
                     samples=int(np.sum(sg[i] > 1e-6))) for i, n in enumerate(names)}
    return dict(feasible=False, iterations=len(hist), history=hist, status="violation remains",
                residual=where)
