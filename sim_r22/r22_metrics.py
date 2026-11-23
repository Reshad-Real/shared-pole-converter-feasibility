"""Benchmark metrics computed the same way for both models (BENCHMARK.md, MODEL_SPEC.md).

Input: a dict with per-period arrays t (period start, s), i_a, i_c, u_cac (sampled at the
period start, which with a symmetric carrier is the period average), U_avg (DC voltage
averaged over each period), dA (leg-a duty applied in that period), plus fine_t/fine_ia for
the THD, and meta with f and Tsw.
"""
import numpy as np


def sine_fit(t, y, f):
    w = 2 * np.pi * f
    A = np.column_stack([np.ones_like(t), np.sin(w * t), np.cos(w * t)])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(np.hypot(c[1], c[2])), float(np.degrees(np.arctan2(c[2], c[1])))


def thd(t, i, f, cycles=5, h_max=50):
    """THD of harmonics 2..h_max over exactly `cycles` periods, resampled to 1 us."""
    T = cycles / f
    tu = t[0] + np.arange(0, T, 1e-6)
    iu = np.interp(tu, t, i)
    X = np.abs(np.fft.rfft(iu)) / len(iu)
    fund = X[cycles]
    harm = X[[h * cycles for h in range(2, h_max + 1)]]
    return float(100 * np.sqrt(np.sum(harm ** 2)) / fund)


def clamp_runs(d, Tsw, eps=1e-6):
    """Complete runs of a clamped duty (0 or 1): list of (level, duration_ms)."""
    lvl = np.where(d >= 1 - eps, 1, np.where(d <= eps, -1, 0))
    runs, cur, start = [], lvl[0], 0
    for j in range(1, len(lvl)):
        if lvl[j] != cur:
            runs.append((cur, start, j)); cur, start = lvl[j], j
    runs.append((cur, start, len(lvl)))
    return [("high" if c > 0 else "low", (b - a) * Tsw * 1e3) for c, a, b in runs
            if c != 0 and a > 0 and b < len(lvl)]


def compute(d):
    f, Tsw = d["meta"]["f"], d["meta"]["Tsw"]
    t = np.asarray(d["t"])
    out = {}
    i_a, i_c = np.asarray(d["i_a"]), np.asarray(d["i_c"])
    for name, y in (("i_a", i_a), ("i_b", -i_a - i_c), ("i_c", i_c), ("u_cac", np.asarray(d["u_cac"]))):
        amp, ph = sine_fit(t, y, f)
        out[f"amp_{name}"] = amp
        out[f"phase_{name}_vs_u_ac_deg"] = float((ph + 180) % 360 - 180)
    out["ucac_over_uac"] = out["amp_u_cac"] / (220 * np.sqrt(2))
    U = np.asarray(d["U_avg"])
    out["dc_ripple_pp_V"] = float(U.max() - U.min())
    out["dc_mean_V"] = float(U.mean())
    out["thd_ia_pct"] = thd(np.asarray(d["fine_t"]), np.asarray(d["fine_ia"]), f)
    runs = clamp_runs(np.asarray(d["dA"]), Tsw)
    out["Sa_clamps"] = runs
    out["Sa_clamp_mean_ms"] = float(np.mean([r[1] for r in runs])) if runs else None
    return out
