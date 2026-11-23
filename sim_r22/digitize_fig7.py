"""Digitise the R22 Fig. 7(a) and 7(b) oscilloscope traces (targets T3-T7 of BENCHMARK.md).

The page is rendered at 300 dpi with pdftoppm. Each trace is separated by hue; for every pixel
column the median row of that hue is taken as the trace. Scales:
- horizontal: 10 divisions across the plot area, 20 ms/div (main) and 5 ms/div (zoom);
- vertical: calibrated in each window on the grid voltage u_ac, whose amplitude is known
  (311 V peak = 0.622 div at 500 V/div), so the result does not rely on the division count.
Amplitudes and phases come from least-squares sine fits; uncertainties combine the fit's
standard error with +/-2 pixels.

Run from the repository root:  py sim_r22/digitize_fig7.py
"""
import colorsys
import json
import subprocess
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
PDF = HERE.parent / "literature" / ("Carrier-Based_Generalized_Discontinuous_PWM_Strategy_for_"
                                    "Single-Phase_Three-Legs_Active_Power_Decoupling_Converters.pdf")
U_AC_PEAK = 220 * np.sqrt(2)
PANELS = {"a_GDPWM": (412, 1201), "b_DPWM1": (1307, 2095)}      # plot-area x limits (px)
MAIN_Y, ZOOM_Y = (312, 505), (514, 702)
HUES = {   # degrees; S_a has a dark fill, so its value threshold is lower
    "i_b": [(350, 360), (0, 12)], "i_c": [(28, 42)], "U_dc": [(48, 72)],
    "u_cac": [(105, 135)], "i_a": [(168, 192)], "u_ac": [(288, 312)], "S_a": [(250, 285)],
}
SCALE = {"U_dc": 100.0, "u_ac": 500.0, "u_cac": 500.0, "i_a": 20.0, "i_b": 20.0, "i_c": 20.0}


def render():
    WORK.mkdir(exist_ok=True)
    png = WORK / "r22p-4.png"
    if not png.exists():
        subprocess.run(["pdftoppm", "-r", "300", "-f", "4", "-l", "4", "-png", str(PDF),
                        str(WORK / "r22p")], check=True, capture_output=True)
    return np.asarray(Image.open(png).convert("RGB"))


def hue_mask(rgb, name):
    f = rgb.reshape(-1, 3) / 255.0
    hsv = np.array([colorsys.rgb_to_hsv(*p) for p in f]).reshape(rgb.shape)
    h, s, v = hsv[..., 0] * 360, hsv[..., 1], hsv[..., 2]
    m = np.zeros(h.shape, bool)
    for lo, hi in HUES[name]:
        m |= (h >= lo) & (h < hi)
    vmin = 0.25 if name == "S_a" else 0.6
    return m & (s > 0.35) & (v > vmin)


def trace(mask, x0, y0, y_lim=None):
    """Median row of the mask in each column (NaN where absent), as -y so up is positive."""
    ys = np.full(mask.shape[1], np.nan)
    for j in range(mask.shape[1]):
        r = np.where(mask[:, j])[0]
        if y_lim is not None:
            r = r[(r + y0 >= y_lim[0]) & (r + y0 <= y_lim[1])]
        if len(r):
            ys[j] = -(np.median(r) + y0)
    return np.arange(mask.shape[1]) + x0, ys


def sine_fit(x, y, period):
    """Fit y = c + a sin(2 pi x / P) + b cos(2 pi x / P); returns amplitude, phase, se, rms resid."""
    k = ~np.isnan(y)
    x, y = x[k], y[k]
    w = 2 * np.pi / period
    A = np.column_stack([np.ones_like(x), np.sin(w * x), np.cos(w * x)])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ coef
    s2 = res @ res / max(len(y) - 3, 1)
    cov = s2 * np.linalg.inv(A.T @ A)
    amp = np.hypot(coef[1], coef[2])
    se = np.sqrt((coef[1] ** 2 * cov[1, 1] + coef[2] ** 2 * cov[2, 2]) / amp ** 2)
    return dict(amp=amp, phase=np.arctan2(coef[2], coef[1]), se=se, resid=np.sqrt(s2),
                offset=coef[0], n=len(y))


def window(img, xl, xr, yw):
    return img[yw[0]:yw[1], xl + 2:xr - 1], xl + 2, yw[0]


def analyse(img, name, xl, xr):
    px_div = (xr - xl) / 10.0
    out, overlay = {}, []
    # ---- main window: calibration, enable instant, DC ripple before decoupling
    sub, x0, y0 = window(img, xl, xr, MAIN_Y)
    xs, uac = trace(hue_mask(sub, "u_ac"), x0, y0)
    fit = sine_fit(xs, uac, 4 * px_div / 4 * 1.0 * 1)   # 20 ms = 1 div at 20 ms/div
    fit = sine_fit(xs, uac, px_div)
    vdiv_main = fit["amp"] / (U_AC_PEAK / 500.0)
    _, ucac = trace(hue_mask(sub, "u_cac"), x0, y0)
    base = np.nanmedian(ucac[:40])
    moved = np.where(np.abs(ucac - base) > 4)[0]
    x_enable = int(xs[moved[0]]) if len(moved) else int(xs[len(xs) // 3])
    _, udc = trace(hue_mask(sub, "U_dc"), x0, y0, y_lim=(MAIN_Y[0], MAIN_Y[0] + 70))
    pre = (xs > xs[0] + 3) & (xs < x_enable - 10)
    rfit = sine_fit(xs[pre], udc[pre], px_div / 2)               # 100 Hz
    unc_px = max(2.0, rfit["se"] * 2)
    out["T3_dc_ripple_pp_no_decoupling_V"] = dict(
        value=2 * rfit["amp"] / vdiv_main * 100.0, unc=2 * unc_px / vdiv_main * 100.0)
    out["calib_main_px_per_div"] = vdiv_main
    out["enable_time_ms"] = (x_enable - xl) / px_div * 20.0
    overlay.append(("main", xs, udc, pre))
    # ---- zoom window (decoupling on, steady state, 5 ms/div)
    sub, x0, y0 = window(img, xl, xr, ZOOM_Y)
    period = 4 * px_div                                          # 20 ms = 4 div
    xs, uac = trace(hue_mask(sub, "u_ac"), x0, y0)
    f_ac = sine_fit(xs, uac, period)
    vdiv = f_ac["amp"] / (U_AC_PEAK / 500.0)
    out["calib_zoom_px_per_div"] = vdiv
    fits = {"u_ac": f_ac}
    for tr in ("u_cac", "i_a", "i_b", "i_c"):
        _, y = trace(hue_mask(sub, tr), x0, y0)
        fits[tr] = sine_fit(xs, y, period)
    for tr in ("u_cac", "i_a", "i_b", "i_c"):
        f = fits[tr]
        unc_px = max(2.0, 2 * f["se"]) + f["amp"] * 2.0 / f_ac["amp"]   # includes calibration
        out[f"amp_{tr}"] = dict(value=f["amp"] / vdiv * SCALE[tr], unc=unc_px / vdiv * SCALE[tr],
                                fit_resid_px=f["resid"])
        out[f"phase_{tr}_vs_u_ac_deg"] = float(np.rad2deg(
            np.angle(np.exp(1j * (f["phase"] - f_ac["phase"])))))
    out["T5_ucac_over_uac"] = dict(value=fits["u_cac"]["amp"] / f_ac["amp"],
                                   unc=2.0 / f_ac["amp"] * (1 + fits["u_cac"]["amp"] / f_ac["amp"]))
    _, udc = trace(hue_mask(sub, "U_dc"), x0, y0, y_lim=(ZOOM_Y[0], ZOOM_Y[0] + 40))
    k = ~np.isnan(udc)
    pp_px = np.percentile(udc[k], 99) - np.percentile(udc[k], 1)
    out["T4_dc_ripple_pp_with_decoupling_upper_bound_V"] = (pp_px + 2) / vdiv * 100.0
    # ---- S_a clamping in the zoom window (gate trace occupies the bottom band)
    band = 60
    m = hue_mask(sub, "S_a")[-band:, :]
    counts = m.sum(axis=1)
    lines = np.where(counts > 0.3 * m.shape[1])[0]      # the two gate levels are long lines
    y_hi, y_lo = lines.min(), lines.max()
    mid = slice(y_hi + 4, y_lo - 3)
    switching = m[mid, :].mean(axis=0) > 0.5
    min_gap = int(1.0 / 5.0 * px_div)                  # gaps under 1 ms are artefact streaks
    runs, cur, start = [], switching[0], 0
    for j in range(1, len(switching)):
        if switching[j] != cur:
            runs.append([cur, start, j]); cur, start = switching[j], j
    runs.append([cur, start, len(switching)])
    for r in runs:
        if not r[0] and (r[2] - r[1]) < min_gap:
            r[0] = True
    clamps = []
    for sw, a, b in runs:
        if not sw and 0 < a and b < len(switching):          # complete clamped intervals only
            level = "high" if m[y_hi:y_hi + 3, a:b].mean() > m[y_lo - 2:y_lo + 1, a:b].mean() else "low"
            clamps.append(dict(level=level, ms=(b - a) / px_div * 5.0))
    out["T7_Sa_clamps"] = clamps
    out["T7_Sa_clamp_mean_ms"] = float(np.mean([c["ms"] for c in clamps])) if clamps else None
    out["T7_unc_ms"] = 2 * 2.0 / px_div * 5.0
    overlay.append(("zoom", xs, uac, None))
    return out


def main():
    img = render()
    res = {}
    for name, (xl, xr) in PANELS.items():
        res[name] = analyse(img, name, xl, xr)
    (HERE / "results").mkdir(exist_ok=True)
    (HERE / "results" / "fig7_digitised.json").write_text(json.dumps(res, indent=2, default=float))
    for k, v in res.items():
        print(k)
        for kk, vv in v.items():
            print("  ", kk, vv)
    draw_check(img)


def draw_check(img):
    """Overlay the extracted traces on the image for a visual check."""
    im = Image.fromarray(img.copy())
    d = ImageDraw.Draw(im)
    col = {"U_dc": (0, 0, 0), "u_ac": (0, 120, 0), "u_cac": (200, 0, 0), "i_a": (0, 0, 200),
           "i_b": (120, 0, 120), "i_c": (0, 120, 120)}
    for xl, xr in PANELS.values():
        for yw in (MAIN_Y, ZOOM_Y):
            sub, x0, y0 = window(img, xl, xr, yw)
            for tr, c in col.items():
                lim = (yw[0], yw[0] + (70 if yw == MAIN_Y else 40)) if tr == "U_dc" else None
                xs, ys = trace(hue_mask(sub, tr), x0, y0, y_lim=lim)
                pts = [(int(x), int(-y)) for x, y in zip(xs, ys) if not np.isnan(y)]
                d.point(pts, fill=c)
    im.crop((380, 290, 2110, 720)).save(WORK / "fig7_digitised_check.png")


if __name__ == "__main__":
    main()
