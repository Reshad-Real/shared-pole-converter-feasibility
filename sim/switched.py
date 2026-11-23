"""Switched (PWM) model of the three-leg circuit with dead time.

Circuit (the same as the averaged model, now with real switches):
  AC1 loop: pole A -> L1, R1 -> e1 -> pole N      AC2 loop: pole B -> L2, R2 -> e2 -> pole N
  buffer:   pole N -> Lb, Rb -> Cb -> negative rail
  bus:      Cdc, fed by DC1 (power source, slew limited) and loaded by DC2 (constant power)

Each leg compares its duty with one shared symmetric carrier at fsw. Every turn-on is delayed
by the dead time td. While both gates of a leg are off, the freewheeling diode selected by the
sign of the pole current sets the pole voltage. Switches have on-resistance Ron; diodes add Vf.

Controller: the averaged model's supervisor and buffer loops (timedomain.py), run once per
switching period on states sampled at the period start. That is the middle of the low-side
interval, where the sample equals the period average. The AC current loops are now real
(reference-voltage feedforward plus a proportional term). Duties include dead-time compensation
and take effect one period after sampling (computation delay).

Integration: within a period the switch states change only at known edges, so the circuit is
integrated edge to edge with Heun steps no longer than h_max.
"""
import math
from dataclasses import dataclass

import numpy as np

from model import ac_state, orbit_energy_at, port_voltage_fn
from timedomain import GAINS, _post_orbit


@dataclass(frozen=True)
class SwitchParams:
    fsw: float = 20e3
    td: float = 1e-6                    # dead time
    Ron: float = 0.05                   # switch and diode on-resistance
    Vf: float = 0.8                     # diode forward drop (dead-time freewheeling)
    delay_periods: int = 1              # sample-to-duty delay
    wc_ac: float = 2 * np.pi * 1500.0   # AC current loop crossover (rad/s)
    dt_comp: bool = True                # dead-time compensation on the duties
    h_max: float = 5e-6
    # Measurement and modulator non-idealities (WP4); the defaults are ideal
    t_min: float = 0.0                  # minimum gate pulse; shorter ones dropped or stretched
    adc_bits: int = 0                   # ADC resolution on every measurement (0 = exact)
    gain_V: float = 1.0                 # bus-voltage sensor gain
    offset_ib: float = 0.0              # buffer-current sensor offset (A)


_EPS = 1e-6
LOW, HIGH, DEAD = 0, 1, 2      # conduction states of a leg
PWM = 3                        # command kind: high in [t_on, t_off) within the period


def _spec(d, tk, Tsw):
    """PWM command of one leg in one period: always low, always high, or high in [t_on, t_off)."""
    if d <= _EPS:
        return (LOW, 0.0, 0.0)
    if d >= 1.0 - _EPS:
        return (HIGH, 0.0, 0.0)
    return (PWM, tk + 0.5 * (1.0 - d) * Tsw, tk + 0.5 * (1.0 + d) * Tsw)


def _cmd(spec, tau):
    kind, a, b = spec
    if kind == LOW:
        return False
    if kind == HIGH:
        return True
    return a <= tau < b


def _leg_state(prev, cur, tk, td, tau):
    """Gate state at tau: a leg conducts high or low only once its command has held for td."""
    sa = _cmd(cur if tau >= tk else prev, tau)
    tb = tau - td
    sb = _cmd(cur if tb >= tk else prev, tb)
    if sa and sb:
        return HIGH
    if not sa and not sb:
        return LOW
    return DEAD


def _sensors(sw, p):
    """Sensor model: bus-voltage gain, buffer-current offset and ADC quantisation.

    Voltage channels span 0 to 1.25 Vdc_max and current channels +/-1.5 I_pole_max.
    Returns None when the measurements are exact.
    """
    if sw.adc_bits == 0 and sw.gain_V == 1.0 and sw.offset_ib == 0.0:
        return None
    n = 2 ** sw.adc_bits
    qv, qi = 1.25 * p.Vdc_max / n, 3.0 * p.I_pole_max / n

    def meas(vb, V, ib, i1, i2):
        V, ib = V * sw.gain_V, ib + sw.offset_ib
        if sw.adc_bits:
            vb, V = round(vb / qv) * qv, round(V / qv) * qv
            ib, i1, i2 = (round(x / qi) * qi for x in (ib, i1, i2))
        return vb, V, ib, i1, i2
    return meas


def _min_pulse(d, Tsw, t_min):
    """Drop or stretch a gate pulse (high or low) shorter than t_min, whichever is nearer."""
    on, off = d * Tsw, (1.0 - d) * Tsw
    if 0.0 < on < t_min:
        return 0.0 if on < 0.5 * t_min else t_min / Tsw
    if 0.0 < off < t_min:
        return 1.0 if off < 0.5 * t_min else 1.0 - t_min / Tsw
    return d


def simulate_switched(p, sw, theta0, dP, policy="conventional", plan=None, pre_cycles=1,
                      gains=None, trip_after=1e-3, zoom_ms=None, plant=None, pre_orbit=None,
                      post_orbit=None):
    """Run one step test in the switched model. Mirrors timedomain.simulate.

    Violations (and the trip) use period-averaged bus voltage and pole currents, the quantities
    the averaged model and the LP describe, so the two models are judged alike. Instantaneous
    extremes, which add switching ripple, are reported separately (`*_inst` metrics).
    zoom_ms=(a, b) records every integration point between a and b ms after the step.
    plant: Params of the actual circuit when it differs from the controller's model `p`
    (component tolerances, DC1 response). The controller's gains, references, energy
    estimates and AC feedforward use `p`; the circuit, port voltages and DC1 slew use `plant`.
    pre_orbit: the pre-step orbit to start on and hold (default: the centred one).
    post_orbit: the post-step orbit the conventional controller recovers to (default: the
    centred one).
    """
    g = dict(GAINS, **(gains or {}))
    q = p if plant is None else plant
    pre, post = _post_orbit(p, dP)
    if pre_orbit is not None:
        pre = pre_orbit
    if post_orbit is not None:
        post = post_orbit
    if policy == "plan_flat":
        # The energy tracker with the bus held at nominal: it follows the plan's DC1 power and
        # stored energy but not the plan's bus excursions. Fed forward as dEd/dt, those
        # excursions (about +-5 V at 120 Hz, on the plan's 83 us grid) ask for buffer-power
        # pulses that saturate pole N (WP5a S4)
        plan = dict(plan, Ed=np.full_like(plan["Ed"], float(p.Ed(p.Vdc_nom))))
        policy = "plan"
    elif policy == "plan_smooth":
        # The energy tracker on a smoothed copy of the plan's bus: a moving average over
        # gains["bus_smooth_n"] plan samples keeps the plan's 120 Hz use of the bus capacitor
        # and drops its sample-to-sample steps
        n = int(g.get("bus_smooth_n", 9))
        Ed = np.pad(plan["Ed"], n // 2, mode="edge")
        plan = dict(plan, Ed=np.convolve(Ed, np.ones(n) / n, mode="valid")[:len(plan["Ed"])])
        policy = "plan"
    Tsw = 1.0 / sw.fsw
    td = sw.td
    t_step = (theta0 % (2 * np.pi)) / p.w + pre_cycles * p.T
    t0 = t_step - pre_cycles * p.T
    n_per = int(round((pre_cycles + p.n_cycles) * p.T / Tsw))
    D = sw.delay_periods

    tk_all = t0 + Tsw * np.arange(n_per)
    t_app = tk_all + (D + 0.5) * Tsw              # middle of the period a new duty acts in
    u1f, u2f, i1r, i2r, pacr = ac_state(p, t_app, t_step, dP)
    _, _, i1s, i2s, _ = ac_state(p, tk_all, t_step, dP)
    P_pre, P_post = pre["Pd"], post["Pd"]
    ramp = np.clip((tk_all - t_step) / p.Tr2, 0.0, 1.0)
    P_cmd = P_pre + (P_post - P_pre) * 0.5 * (1.0 - np.cos(np.pi * ramp))
    Eb_ref = np.where(tk_all < t_step, pre["Eb"].mean(), post["Eb"].mean())
    if policy in ("plan", "circuit_plan", "state_plan"):
        tp = plan["t"] + t_step
        dEd_plan = np.gradient(plan["Ed"], p.Ts)
    # Optional (off by default; WP5a tracker): sync_ib compares the buffer current with its
    # reference at the measurement instant, as the AC loops do. The reference is computed for
    # the period the duty acts in, (D + 0.5) Tsw later, so without it the buffer current runs
    # about 110 us ahead of the ripple it should carry
    sync_ib = bool(g.get("sync_ib", False))

    Kp1, Kp2 = p.L1 * sw.wc_ac, p.L2 * sw.wc_ac
    ev = port_voltage_fn(q, t_step)
    L1, L2, Lb, R1, R2, Rb = q.L1, q.L2, q.Lb, q.R1, q.R2, q.Rb       # circuit
    Cb, Cdc, Ron, Vf = q.Cb, q.Cdc, sw.Ron, sw.Vf
    Cb_c, Cdc_c, Lb_c, Rb_c = p.Cb, p.Cdc, p.Lb, p.Rb                 # controller's model
    Ed_nom = float(p.Ed(p.Vdc_nom))
    comp = td / Tsw if sw.dt_comp else 0.0
    meas = _sensors(sw, p)

    # Initial state on the pre-step orbit (the buffer voltage the controller holds)
    Eb_i = orbit_energy_at(pre, p, t0)
    vb = math.sqrt(2 * Eb_i / Cb_c)
    ib = float(np.interp(np.mod(t0, p.T), pre["t"], pre["pb"], period=p.T)) / vb
    i1, i2 = float(i1s[0]), float(i2s[0])
    V = p.Vdc_nom
    pd = P_pre
    M = int(round(0.5 * p.T / Tsw))
    hist = np.array([orbit_energy_at(pre, p, t0 - (M - j) * Tsw) for j in range(M)])
    hsum, hpos = float(hist.sum()), 0

    def pole(state, i, Vb):
        """Pole voltage and whether the pole is tied to the positive rail."""
        if state == HIGH:
            return Vb - Ron * i, 1.0
        if state == LOW:
            return -Ron * i, 0.0
        if i > 0.0:                               # lower diode freewheels
            return -Vf - Ron * i, 0.0
        return Vb + Vf - Ron * i, 1.0             # upper diode freewheels

    def deriv(i1_, i2_, ib_, vb_, V_, t, sA, sB, sN):
        iN = ib_ - i1_ - i2_
        vA, cA = pole(sA, i1_, V_)
        vB, cB = pole(sB, i2_, V_)
        vN, cN = pole(sN, iN, V_)
        e1, e2 = ev(t)
        idc = cA * i1_ + cB * i2_ + cN * iN
        return ((vA - vN - e1 - R1 * i1_) / L1, (vB - vN - e2 - R2 * i2_) / L2,
                (vN - vb_ - Rb * ib_) / Lb, ib_ / Cb,
                ((pd - q.Pdc2) / V_ - idc) / Cdc, vN)

    keys = ("t", "vb", "ib", "Vdc", "i1", "i2", "pd", "vN_cmd", "lo", "hi", "dA", "dB", "dN",
            "Vdc_avg", "vb_avg", "vN_avg", "Vdc_lo", "Vdc_hi", "ib_ref", "ib_avg", "pN_avg")
    rec = {k: [] for k in keys}
    zoom = {k: [] for k in ("t", "vN", "ib", "iN", "Vdc", "vb", "i1", "i2")}
    queue = []
    prev_specs = [(LOW, 0.0, 0.0)] * 3
    ib_ref_prev = None
    headroom, sat_time, synth_time, clamps = math.inf, 0.0, 0.0, 0
    t_viol = None
    V_inst_lo, V_inst_hi, I_pk = math.inf, -math.inf, 0.0
    V_avg_lo, V_avg_hi, I_avg_pk = math.inf, -math.inf, 0.0
    V_pre_lo, V_pre_hi = math.inf, -math.inf     # last line cycle before the step
    ac_err2, ac_n = 0.0, 0
    state0 = None

    for k in range(n_per):
        tk = float(tk_all[k])
        after = tk >= t_step
        if after and state0 is None:
            state0 = dict(Eb0=0.5 * Cb * vb * vb, Ed0=0.5 * Cdc * V * V, pd0=pd, vb=vb, V=V,
                          ib=ib, i1=i1, i2=i2, t_ms=1e3 * (tk - t_step))

        # ---- sample and control (the duty takes effect D periods later)
        if meas is None:
            vb_m, V_m, ib_m, i1_m, i2_m = vb, V, ib, i1, i2
        else:
            vb_m, V_m, ib_m, i1_m, i2_m = meas(vb, V, ib, i1, i2)
        Eb = 0.5 * Cb_c * vb_m * vb_m
        Ed = 0.5 * Cdc_c * V_m * V_m
        Eb_avg = hsum / M
        in_plan = policy in ("circuit_plan", "state_plan") and after and t_app[k] <= tp[-1]
        if policy in ("plan", "circuit_plan", "state_plan") and after:
            j = min(int((tk - t_step) / p.Ts), len(plan["pd"]) - 1)
            Ed_ref = float(np.interp(t_app[k], tp, plan["Ed"]))
            Eb_pl = float(np.interp(tk, tp, plan["Eb"]))
            dEd_ff = float(np.interp(t_app[k], tp, dEd_plan)) if policy == "plan" else 0.0
            pd_cmd = plan["pd"][j] + g["k_tot"] * ((Eb_pl + float(np.interp(tk, tp, plan["Ed"])))
                                                  - (Eb + Ed))
        else:
            Ed_ref, dEd_ff = Ed_nom, 0.0
            pd_cmd = P_cmd[k] + g["k_rec"] * (Eb_ref[k] - Eb_avg)
        stp = q.slew_dc1 * Tsw                      # DC1 is part of the plant
        pd = min(max(pd + min(max(pd_cmd - pd, -stp), stp), q.Pdc1_min), q.Pdc1_max)

        if in_plan:
            # Circuit plan: follow its buffer current, with bus-energy feedback
            ib_ref = (float(np.interp(t_app[k], tp, plan["ib"]))
                      + g["k_bus"] * (Ed - Ed_ref) / max(vb_m, 1.0))
        else:
            pb_ref = pd - p.Pdc2 - pacr[k] - dEd_ff + g["k_bus"] * (Ed - Ed_ref)
            ib_ref = pb_ref / max(vb_m, 1.0)
        dib_ff = 0.0 if ib_ref_prev is None else (ib_ref - ib_ref_prev) / Tsw
        ib_ref_prev = ib_ref
        if in_plan and policy == "state_plan":
            # Follow the plan's own state: its pole-N voltage as feedforward, corrected for the
            # buffer-voltage and buffer-current errors (the plan is a trajectory of the
            # averaged circuit, so it nearly reproduces itself open loop)
            jp = min(int((t_app[k] - t_step) / p.Ts), len(plan["vN"]) - 1)
            vb_pl = float(np.interp(tk, tp, plan["vb"]))
            ib_pl = float(np.interp(tk, tp, plan["ib"]))
            vN_cmd = (plan["vN"][jp] + (vb_m - vb_pl) + Rb_c * (ib_m - ib_pl)
                      + Lb_c * g["wc"] * (ib_ref - ib_m))
        else:
            e_ib = ib_ref - (D + 0.5) * Tsw * dib_ff - ib_m if sync_ib else ib_ref - ib_m
            vN_cmd = vb_m + Rb_c * ib_m + Lb_c * (dib_ff + g["wc"] * e_ib)
        u1c = u1f[k] + Kp1 * (i1s[k] - i1_m)
        u2c = u2f[k] + Kp2 * (i2s[k] - i2_m)
        lo = max(0.0, -u1c, -u2c)
        hi = V_m - max(0.0, u1c, u2c)
        synth = max(lo - hi, 0.0)
        vN_c = 0.5 * (lo + hi) if synth > 0 else min(max(vN_cmd, lo), hi)
        duties = []
        for v_ref, i_ref in ((u1c + vN_c, i1r[k]), (u2c + vN_c, i2r[k]),
                             (vN_c, ib_ref - i1r[k] - i2r[k])):
            d = v_ref / V_m + comp * max(-1.0, min(1.0, i_ref / 0.3))
            dc = min(max(d, 0.0), 1.0)
            clamps += after and dc != d
            if sw.t_min > 0.0:
                dc = _min_pulse(dc, Tsw, sw.t_min)
            duties.append(dc)
        queue.append(duties)
        applied = queue.pop(0) if len(queue) > D else duties

        if after:
            sat_time += Tsw * (vN_c != vN_cmd)
            synth_time += Tsw * (synth > 0)
            headroom = min(headroom, vN_cmd - lo, hi - vN_cmd)
            ac_err2 += (i1s[k] - i1) ** 2 + (i2s[k] - i2) ** 2
            ac_n += 2

        # ---- integrate this period edge to edge
        specs = [_spec(d, tk, Tsw) for d in applied]
        t_end = tk + Tsw
        cand = {tk + td}
        for sp_ in specs + prev_specs:
            if sp_[0] == PWM:
                cand.update((sp_[1], sp_[1] + td, sp_[2], sp_[2] + td))
        edges = [tk] + sorted(c for c in cand if tk < c < t_end) + [t_end]
        x = [i1, i2, ib, vb, V]
        intV = intvb = intvN = intiA = intiB = intiN = intib = intpN = 0.0
        pV_lo, pV_hi = math.inf, -math.inf
        z_on = zoom_ms is not None and zoom_ms[0] <= 1e3 * (tk - t_step) <= zoom_ms[1]
        for a, b in zip(edges[:-1], edges[1:]):
            n_sub = max(1, math.ceil((b - a) / sw.h_max - 1e-9))
            h = (b - a) / n_sub
            tm = 0.5 * (a + b)
            sA = _leg_state(prev_specs[0], specs[0], tk, td, tm)
            sB = _leg_state(prev_specs[1], specs[1], tk, td, tm)
            sN = _leg_state(prev_specs[2], specs[2], tk, td, tm)
            ta = a
            for _ in range(n_sub):
                k1 = deriv(*x, ta, sA, sB, sN)
                xp = [x[q] + h * k1[q] for q in range(5)]
                k2 = deriv(*xp, ta + h, sA, sB, sN)
                xn = [x[q] + 0.5 * h * (k1[q] + k2[q]) for q in range(5)]
                intV += 0.5 * h * (x[4] + xn[4])
                intvb += 0.5 * h * (x[3] + xn[3])
                intvN += 0.5 * h * (k1[5] + k2[5])
                intiA += 0.5 * h * (x[0] + xn[0])
                intiB += 0.5 * h * (x[1] + xn[1])
                intiN += 0.5 * h * (x[2] - x[0] - x[1] + xn[2] - xn[0] - xn[1])
                intib += 0.5 * h * (x[2] + xn[2])
                intpN += 0.5 * h * (k1[5] * x[2] + k2[5] * xn[2])
                if z_on:
                    for key, val in zip(("t", "vN", "ib", "iN", "Vdc", "vb", "i1", "i2"),
                                        (ta - t_step, k1[5], x[2], x[2] - x[0] - x[1], x[4],
                                         x[3], x[0], x[1])):
                        zoom[key].append(val)
                x = xn
                ta += h
                pV_lo, pV_hi = min(pV_lo, x[4]), max(pV_hi, x[4])
                if after:
                    I_pk = max(I_pk, abs(x[0]), abs(x[1]), abs(x[2] - x[0] - x[1]))
            x[4] = max(x[4], 1.0)
        prev_specs = specs

        V_avg = intV / Tsw
        for key, val in zip(keys, (tk - t_step, vb, ib, V, i1, i2, pd, vN_cmd, lo, hi,
                                   applied[0], applied[1], applied[2], V_avg, intvb / Tsw,
                                   intvN / Tsw, pV_lo, pV_hi, ib_ref, intib / Tsw, intpN / Tsw)):
            rec[key].append(val)
        i1, i2, ib, vb, V = x
        hsum += Eb - hist[hpos]
        hist[hpos] = Eb
        hpos = (hpos + 1) % M

        if not after and tk >= t_step - p.T:
            V_pre_lo, V_pre_hi = min(V_pre_lo, V_avg), max(V_pre_hi, V_avg)
        if after:
            V_inst_lo, V_inst_hi = min(V_inst_lo, pV_lo), max(V_inst_hi, pV_hi)
            V_avg_lo, V_avg_hi = min(V_avg_lo, V_avg), max(V_avg_hi, V_avg)
            I_avg_pk = max(I_avg_pk, abs(intiA) / Tsw, abs(intiB) / Tsw, abs(intiN) / Tsw)
            bad = (V_avg < p.Vdc_min or V_avg > p.Vdc_max or I_avg_pk > p.I_pole_max
                   or synth > 0)
            if bad and t_viol is None:
                t_viol = tk
            if t_viol is not None and tk > t_viol + trip_after:
                break

    out = {k: np.asarray(v) for k, v in rec.items()}
    out["zoom"] = {k: np.asarray(v) for k, v in zoom.items()}
    out["state_at_step"] = state0
    m = dict(Vdc_min=V_avg_lo, Vdc_max=V_avg_hi, Vdc_inst_min=V_inst_lo, Vdc_inst_max=V_inst_hi,
             I_pole_peak=I_avg_pk, I_pole_peak_inst=I_pk, sat_time_ms=1e3 * sat_time,
             synth_time_ms=1e3 * synth_time, duty_clamps=int(clamps), min_vN_headroom_V=headroom,
             ac_current_rms_error_A=math.sqrt(ac_err2 / max(ac_n, 1)),
             bus_violation=bool(V_avg_lo < p.Vdc_min or V_avg_hi > p.Vdc_max),
             current_violation=bool(I_avg_pk > p.I_pole_max),
             current_violation_inst=bool(I_pk > p.I_pole_max),
             synth_violation=bool(synth_time > 0),
             t_violation_ms=None if t_viol is None else 1e3 * (t_viol - t_step),
             pre_Vdc_min=V_pre_lo, pre_Vdc_max=V_pre_hi,
             pre_step_violation=bool(V_pre_lo < p.Vdc_min or V_pre_hi > p.Vdc_max))
    m["violation"] = m["bus_violation"] or m["current_violation"] or m["synth_violation"]
    out["metrics"] = {k: (float(v) if isinstance(v, (float, np.floating)) else v)
                      for k, v in m.items()}
    return out
