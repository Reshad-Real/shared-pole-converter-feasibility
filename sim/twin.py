"""Digital twin of the shared-pole converter (WP9; protocol in WP9_TWIN_PROTOCOL.md).

There is no hardware, so the physical converter is emulated: the switched model (switched.py)
runs a plant whose parameters are hidden from the twin. The twin sees only what the converter's
controller has: samples at the start of every switching period of the buffer, bus and port
voltages and currents, the DC1 and DC2 powers, and the duties it applied. Every measured channel
has a 12-bit ADC, white noise and a gain error fixed per converter.

The twin identifies the buffer capacitance and inductance, the port voltage amplitudes and the
DC1 slew from two recordings of normal operation (steady operation and an earlier load event),
by least squares on windowed integral forms of the circuit equations. It also estimates the
buffer resistance and the bus capacitance but keeps their nominal values (see USED). The twin
is the averaged circuit model with the identified parameters; the paper's bounds are evaluated
on it exactly as on the nominal model.
"""
from dataclasses import replace

import numpy as np

from model import port_voltages
from switched import SwitchParams, simulate_switched

ADC_BITS = 12
NOISE_LSB = 1.0             # white sensor noise, standard deviation in LSB
GAIN_TOL = 0.005            # sensor gain error, uniform within +/- this, fixed per converter
WINDOW = 8                  # switching periods per integration window
# Windows in which a pole current comes closer to zero than ZC_FACTOR times the nominal
# ripple amplitude are left out of the Lb and bus regressions: there the current changes sign
# within a period, and the dead-time correction of the pole-voltage estimate has the wrong sign
ZC_FACTOR = 2.5
CLAMP = dict(Cb=(0.5, 2.0), E1=(0.5, 2.0), E2=(0.5, 2.0), slew_dc1=(0.5, 2.0), Lb=(0.5, 2.0))
# Identified and used by the twin. Cdc and Rb are estimated and reported, but the twin keeps
# their nominal values: the method test found neither practically identifiable from normal
# operation, and neither moves the admissibility boundary by more than 0.5 % of the AC2
# rating within its tolerance (WP9_TWIN_PROTOCOL.md, revision 1)
USED = ("Cb", "Lb", "E1", "E2", "slew_dc1")
CHANNELS = ("vb", "Vdc", "ib", "i1", "i2", "e1", "e2", "pd", "pdc2")


def _span(p, ch):
    """ADC input range of a channel."""
    if ch in ("vb", "Vdc"):
        return 0.0, 1.25 * p.Vdc_max
    if ch in ("e1", "e2"):
        return -0.625 * p.Vdc_max, 0.625 * p.Vdc_max
    if ch in ("ib", "i1", "i2"):
        return -1.5 * p.I_pole_max, 1.5 * p.I_pole_max
    return p.Pdc1_min, p.Pdc1_max


def sensor_gains(rng):
    """Gain errors of one converter's sensors."""
    return {ch: 1.0 + GAIN_TOL * rng.uniform(-1.0, 1.0) for ch in CHANNELS}


def record(pc, plant, sw, theta0, dP, pre_cycles, n_cycles, gains=None):
    """Run the plant in the switched model and return the run with absolute sample times.

    pc: the controller's model (the converter's own parameter values); plant: the actual circuit.
    """
    q = replace(pc, n_cycles=n_cycles)
    run = simulate_switched(q, sw, theta0, dP, pre_cycles=pre_cycles, plant=plant, gains=gains)
    t_step = (theta0 % (2 * np.pi)) / pc.w + pre_cycles * pc.T
    run["t_abs"] = run["t"] + t_step
    run["t_step"] = t_step
    run["Tsw"] = 1.0 / sw.fsw
    run["td"] = sw.td
    return run


def measure(run, p, plant, gains, rng):
    """The samples the twin receives from one run: quantised, noisy, with gain errors.

    p: nominal design values (channel ranges); plant: the actual circuit (port voltages, DC2).
    """
    e1, e2 = port_voltages(plant, run["t_abs"], run["t_step"])
    true = dict(vb=run["vb"], Vdc=run["Vdc"], ib=run["ib"], i1=run["i1"], i2=run["i2"], e1=e1,
                e2=e2, pd=run["pd"], pdc2=np.full_like(run["vb"], plant.Pdc2))
    out = {}
    for ch, x in true.items():
        lo, hi = _span(p, ch)
        q = (hi - lo) / 2 ** ADC_BITS
        y = gains[ch] * np.asarray(x, float) + NOISE_LSB * q * rng.standard_normal(len(x))
        out[ch] = np.clip(lo + np.round((y - lo) / q) * q, lo, hi)
    out.update(dA=run["dA"], dB=run["dB"], dN=run["dN"], t=run["t"], t_abs=run["t_abs"],
               Tsw=run["Tsw"], td=run["td"])
    return out


def pole_voltages(m, Ron=0.05):
    """Period-average pole voltages estimated from the applied duties and the measured bus.

    Dead time costs (adds) td of high time per period when the pole current is positive
    (negative); the conducting device drops Ron i.
    """
    k = m["td"] / m["Tsw"]
    iN = m["ib"] - m["i1"] - m["i2"]
    V = m["Vdc"]
    vA = (m["dA"] - k * np.sign(m["i1"])) * V - Ron * m["i1"]
    vB = (m["dB"] - k * np.sign(m["i2"])) * V - Ron * m["i2"]
    vN = (m["dN"] - k * np.sign(iN)) * V - Ron * iN
    return vA, vB, vN


def _win(x, w):
    """Window increments x[k+w] - x[k] and trapezoid integrals (in sample units) per window."""
    n = (len(x) - 1) // w
    idx = np.arange(n) * w
    inc = x[idx + w] - x[idx]
    trap = 0.5 * (x[:-1] + x[1:])
    integ = np.add.reduceat(trap[:n * w], idx)
    return inc, integ


def _period_sum(x, w):
    """Sum over each window of a per-period quantity x[k] (k = start of the period)."""
    n = (len(x) - 1) // w
    idx = np.arange(n) * w
    return np.add.reduceat(x[:n * w], idx)


def _amplitude(t, e, w):
    """Amplitude of the fundamental of e(t) at angular frequency w (least squares)."""
    A = np.column_stack([np.sin(w * t), np.cos(w * t), np.ones_like(t)])
    c, *_ = np.linalg.lstsq(A, e, rcond=None)
    return float(np.hypot(c[0], c[1]))


def _slew(m, dP_cal):
    """DC1 slew from the response to an earlier load step: slope of the measured DC1 power
    between 20 % and 80 % of its rise."""
    t, pd = m["t"], m["pd"]
    after = t >= 0.0
    p0 = float(np.median(pd[~after][-50:]))
    rise = pd[after] - p0
    ta = t[after]
    i20 = np.argmax(rise >= 0.2 * dP_cal)
    i80 = np.argmax(rise >= 0.8 * dP_cal)
    if i80 <= i20 + 3:
        return np.nan
    c = np.polyfit(ta[i20:i80 + 1], rise[i20:i80 + 1], 1)
    return float(c[0])


def identify(p, records, dP_cal, Ron=0.05, w=WINDOW):
    """Identified parameters from measured records (list of dicts from measure()).

    records[0] is steady operation, records[1] the earlier load event of size dP_cal.
    Least squares on windowed integral forms of the circuit equations:
      buffer capacitor   Cb (vb[k+w] - vb[k])          = int ib dt
      buffer inductor    Lb (ib[k+w] - ib[k]) + Rb int ib dt = int (vN - vb) dt
      bus                Cdc (V[k+w] - V[k]) - b w Tsw  = int ((pd - pdc2)/V - idc) dt
    with the pole voltages and the bus current idc estimated from the duties.
    Returns (twin params, details).
    """
    a_C, b_C, y_L, X_L, y_D, X_D, E1, E2 = [], [], [], [], [], [], [], []
    for m in records:
        Tsw = m["Tsw"]
        thr = ZC_FACTOR * 0.25 * p.Vdc_nom * Tsw / p.Lb     # nominal ripple amplitude x factor
        n = (len(m["vb"]) - 1) // w
        idx = np.arange(n) * w
        dvb, _ = _win(m["vb"], w)
        dib, ib_int = _win(m["ib"], w)
        a_C.append(dvb)
        b_C.append(Tsw * ib_int)
        vA, vB, vN = pole_voltages(m, Ron)
        iN = m["ib"] - m["i1"] - m["i2"]
        vb_mid = 0.5 * (m["vb"][:-1] + m["vb"][1:])
        ok_L = np.minimum.reduceat(np.abs(iN[:n * w]), idx) >= thr
        y_L.append((Tsw * _period_sum(vN[:-1] - vb_mid, w))[ok_L])
        X_L.append(np.column_stack([dib, Tsw * ib_int])[ok_L])
        k = m["td"] / Tsw
        idc = ((m["dA"] - k * np.sign(m["i1"])) * m["i1"] + (m["dB"] - k * np.sign(m["i2"])) * m["i2"]
               + (m["dN"] - k * np.sign(iN)) * iN)
        src = (m["pd"] - m["pdc2"]) / m["Vdc"] - idc
        dV, _ = _win(m["Vdc"], w)
        i_min = np.minimum(np.minimum(np.abs(m["i1"]), np.abs(m["i2"])), np.abs(iN))
        ok_D = np.minimum.reduceat(i_min[:n * w], idx) >= thr
        y_D.append(dV[ok_D])
        X_D.append(np.column_stack([Tsw * _period_sum(src[:-1], w), np.full(n, w * Tsw)])[ok_D])
        E1.append(_amplitude(m["t_abs"], m["e1"], p.w))
        E2.append(_amplitude(m["t_abs"], m["e2"], p.w))
    a_C, b_C = np.concatenate(a_C), np.concatenate(b_C)
    Cb = float(b_C @ b_C / (b_C @ a_C))
    yL, XL = np.concatenate(y_L), np.vstack(X_L)
    (Lb, Rb), *_ = np.linalg.lstsq(XL, yL, rcond=None)
    yD, XD = np.concatenate(y_D), np.vstack(X_D)
    Cdc = np.nan
    if len(yD) >= 3:
        (inv_Cdc, _), *_ = np.linalg.lstsq(XD, yD, rcond=None)
        Cdc = 1.0 / inv_Cdc if inv_Cdc > 0 else np.nan
    est = dict(Cb=Cb, Lb=float(Lb), E1=float(np.mean(E1)), E2=float(np.mean(E2)),
               slew_dc1=_slew(records[1], dP_cal))
    clamped = {}
    for key, (lo, hi) in CLAMP.items():
        nom = getattr(p, key)
        x = est[key]
        if not np.isfinite(x) or x < lo * nom or x > hi * nom:
            clamped[key] = x
            est[key] = float(np.clip(x if np.isfinite(x) else nom, lo * nom, hi * nom))
    twin = replace(p, **{k: est[k] for k in USED})
    return twin, dict(estimates=est, reported=dict(Rb=float(Rb), Cdc=float(Cdc)),
                      clamped=clamped, windows=dict(C=len(a_C), L=len(yL), D=len(yD)))


def parity_residual(twin, m, w=WINDOW):
    """Drift monitor: RMS over windows of the buffer-capacitor equation error
    (vb[k+w] - vb[k]) - int ib dt / Cb, with the twin's Cb (V)."""
    dvb, _ = _win(m["vb"], w)
    _, ib_int = _win(m["ib"], w)
    return float(np.sqrt(np.mean((dvb - m["Tsw"] * ib_int / twin.Cb) ** 2)))


def calibrate(p, plant, theta_cal, dP_cal, gains, rng, sw=None):
    """Record the two calibration runs on the plant and identify the twin.

    The converter's controller uses the nominal values with the measured port voltage
    amplitudes (here the plant's, as its PLL measures them).
    """
    sw = SwitchParams(adc_bits=ADC_BITS) if sw is None else sw
    pc = replace(p, E1=plant.E1, E2=plant.E2)
    r1 = record(pc, plant, sw, theta_cal, 0.0, pre_cycles=3, n_cycles=2)
    r2 = record(pc, plant, sw, theta_cal, dP_cal, pre_cycles=3, n_cycles=2)
    held = not (r1["metrics"]["violation"] or r2["metrics"]["violation"])
    m1, m2 = measure(r1, p, plant, gains, rng), measure(r2, p, plant, gains, rng)
    twin, info = identify(p, [m1, m2], dP_cal)
    info["calibration_runs_held"] = held
    info["residual_cal_V"] = parity_residual(twin, m1)
    return twin, info, (m1, m2), (r1, r2)
