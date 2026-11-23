"""Time-domain averaged simulation used to check LP predictions independently.

Unlike the LP, this model keeps the buffer inductor and its resistance, the nonlinear bus
(Cdc * Vdc * dVdc/dt = power balance), a cascaded buffer controller with a finite current-loop
bandwidth, and continuous-time evaluation of the pole limits. AC currents are prescribed
(ideal AC current loops), so the AC ports always receive the requested power while pole N
stays inside its band; when the band is empty the AC voltages cannot be synthesised.

Two supervisory policies are available:
  conventional : bus-energy loop on the buffer, DC1 = commanded power + buffer-energy recovery
  plan         : same inner loops, but DC1 and the bus energy follow an LP trajectory
"""
import numpy as np

from model import ac_state, band, orbit_energy_at, steady_orbit
from feasibility import nominal_orbit

GAINS = dict(wc=2 * np.pi * 2000.0,      # buffer current loop (rad/s)
             k_bus=2 * np.pi * 200.0,    # bus energy loop (1/s)
             k_rec=2 * np.pi * 10.0,     # buffer mean-energy recovery through DC1 (1/s)
             k_tot=2 * np.pi * 10.0)     # plan: total-energy correction on DC1 (1/s)


def _post_orbit(p, dP):
    pre = nominal_orbit(p)
    Ed_lo, Ed_hi = p.Ed(p.Vdc_min), p.Ed(p.Vdc_max)
    V_sec = p.Vdc_min + (p.Vdc_max - p.Vdc_min) * (p.Ed(p.Vdc_nom) - Ed_lo) / (Ed_hi - Ed_lo)
    post = steady_orbit(p, p.P2_pre + dP, V_sec, p.delta_in)
    return pre, post if post is not None else pre


def simulate(p, theta0, dP, policy="conventional", plan=None, pre_cycles=1, dt=5e-6,
             gains=None, record_every=4, trip_after=1e-3):
    """Run one step test. A protective trip stops the run trip_after seconds after the first
    predeclared violation (bus outside its band, pole over-current, or an empty pole band)."""
    g = dict(GAINS, **(gains or {}))
    pre, post = _post_orbit(p, dP)
    t_step = (theta0 % (2 * np.pi)) / p.w + pre_cycles * p.T
    t0 = t_step - pre_cycles * p.T
    n = int(round((pre_cycles + p.n_cycles) * p.T / dt))
    tt = t0 + 0.5 * dt * np.arange(2 * n + 1)
    u1, u2, i1, i2, pac = ac_state(p, tt, t_step, dP)
    lo_all, m_all = band(u1, u2)

    # Commanded average power: what a conventional supervisor feeds forward to DC1
    P_pre, P_post = pre["Pd"], post["Pd"]
    ramp = np.clip((tt - t_step) / p.Tr2, 0.0, 1.0)
    P_cmd = P_pre + (P_post - P_pre) * 0.5 * (1.0 - np.cos(np.pi * ramp))
    Eb_ref = np.where(tt < t_step, pre["Eb"].mean(), post["Eb"].mean())

    if policy == "plan":
        tp = plan["t"] + t_step                  # node times of the LP trajectory
        dEd_plan = np.gradient(plan["Ed"], p.Ts)  # node slopes, interpolated for a smooth feedforward

    # Initial state on the pre-step orbit
    Eb_i = orbit_energy_at(pre, p, t0)
    vb = np.sqrt(2 * Eb_i / p.Cb)
    pb_i = float(np.interp(np.mod(t0, p.T), pre["t"], pre["pb"], period=p.T))
    ib = pb_i / vb
    Vdc = p.Vdc_nom
    pd = P_pre
    M = int(round(0.5 * p.T / dt))               # half-cycle moving average of Eb
    hist = np.array([orbit_energy_at(pre, p, t0 - (M - j) * dt) for j in range(M)])
    hsum, hpos = hist.sum(), 0
    Ed_nom = float(p.Ed(p.Vdc_nom))

    keys = ("t", "vb", "ib", "Vdc", "pd", "vN", "lo", "hi", "sat", "synth", "iN", "iA", "iB")
    rec = {k: [] for k in keys}
    sat_time = 0.0
    synth_time = 0.0
    t_viol = None
    ib_ref_prev = None
    headroom = np.inf

    def f(x, vN, k2):
        ib_, vb_, V_ = x
        return np.array([(vN - vb_ - p.Rb * ib_) / p.Lb,
                         ib_ / p.Cb,
                         (pd - p.Pdc2 - pac[k2] - vN * ib_) / (p.Cdc * V_)])

    for k in range(n):
        k2 = 2 * k
        t = tt[k2]
        Eb = 0.5 * p.Cb * vb * vb
        Ed = 0.5 * p.Cdc * Vdc * Vdc
        Eb_avg = hsum / M
        after = t >= t_step

        if policy == "plan" and after:
            j = min(int((t - t_step) / p.Ts), len(plan["pd"]) - 1)
            Ed_ref = float(np.interp(t, tp, plan["Ed"]))
            Eb_pl = float(np.interp(t, tp, plan["Eb"]))
            dEd_ff = float(np.interp(t, tp, dEd_plan))
            pd_cmd = plan["pd"][j] + g["k_tot"] * ((Eb_pl + Ed_ref) - (Eb + Ed))
        else:
            Ed_ref, dEd_ff = Ed_nom, 0.0
            pd_cmd = P_cmd[k2] + g["k_rec"] * (Eb_ref[k2] - Eb_avg)

        step = p.slew_dc1 * dt
        pd = float(np.clip(pd + np.clip(pd_cmd - pd, -step, step), p.Pdc1_min, p.Pdc1_max))

        pb_ref = pd - p.Pdc2 - pac[k2] - dEd_ff + g["k_bus"] * (Ed - Ed_ref)
        ib_ref = pb_ref / max(vb, 1.0)
        dib_ff = 0.0 if ib_ref_prev is None else (ib_ref - ib_ref_prev) / dt
        ib_ref_prev = ib_ref
        vN_cmd = vb + p.Rb * ib + p.Lb * (dib_ff + g["wc"] * (ib_ref - ib))
        lo_t, hi_t = lo_all[k2], Vdc - m_all[k2]
        synth = max(lo_t - hi_t, 0.0)
        if synth > 0:
            vN = 0.5 * (lo_t + hi_t)
            sat = True
        else:
            vN = min(max(vN_cmd, lo_t), hi_t)
            sat = vN != vN_cmd
        if after:
            sat_time += dt * sat
            synth_time += dt * (synth > 0)
            headroom = min(headroom, vN_cmd - lo_t, hi_t - vN_cmd)
            i_pk = max(abs(i1[k2]), abs(i2[k2]), abs(ib - i1[k2] - i2[k2]))
            bad = (Vdc < p.Vdc_min or Vdc > p.Vdc_max or i_pk > p.I_pole_max or synth > 0)
            if bad and t_viol is None:
                t_viol = t
            if t_viol is not None and t > t_viol + trip_after:
                break

        if k % record_every == 0:
            for key, val in zip(keys, (t - t_step, vb, ib, Vdc, pd, vN, lo_t, hi_t, sat, synth,
                                       ib - i1[k2] - i2[k2], i1[k2], i2[k2])):
                rec[key].append(val)

        x = np.array([ib, vb, Vdc])
        a = f(x, vN, k2)
        b = f(x + 0.5 * dt * a, vN, k2 + 1)
        c = f(x + 0.5 * dt * b, vN, k2 + 1)
        d = f(x + dt * c, vN, k2 + 2)
        ib, vb, Vdc = x + dt / 6 * (a + 2 * b + 2 * c + d)
        Vdc = max(Vdc, 1.0)

        Eb_new = 0.5 * p.Cb * vb * vb
        hsum += Eb_new - hist[hpos]
        hist[hpos] = Eb_new
        hpos = (hpos + 1) % M

    out = {k: np.asarray(v) for k, v in rec.items()}
    post_mask = out["t"] >= 0
    Ipk = np.max(np.abs(np.stack([out["iN"], out["iA"], out["iB"]])[:, post_mask]))
    Vmin_, Vmax_ = out["Vdc"][post_mask].min(), out["Vdc"][post_mask].max()
    out["metrics"] = dict(
        Vdc_min=float(Vmin_), Vdc_max=float(Vmax_), I_pole_peak=float(Ipk),
        sat_time_ms=1e3 * sat_time, synth_time_ms=1e3 * synth_time,
        bus_violation=bool(Vmin_ < p.Vdc_min or Vmax_ > p.Vdc_max),
        current_violation=bool(Ipk > p.I_pole_max),
        synth_violation=bool(synth_time > 0),
        t_violation_ms=None if t_viol is None else 1e3 * (t_viol - t_step),
        min_vN_headroom_V=float(headroom),   # negative: the commanded vN left the band
    )
    m_ = out["metrics"]
    m_["violation"] = m_["bus_violation"] or m_["current_violation"] or m_["synth_violation"]
    return out
