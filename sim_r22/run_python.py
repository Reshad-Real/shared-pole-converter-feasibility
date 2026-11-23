"""Run the R22 benchmark cases in the Python switched model and save the metrics.

Cases: decoupling off; rectifier at phi = 0, 15, 45 deg (80, 83, 113 ohm) with SVPWM, DPWM1
and CB-GDPWM; STATCOM at 5 and 10 A rms; dead-time variants 100 and 300 ns at phi = 0.
Output: results/python_metrics.json and results/python_waveforms.json (phi = 0 cases).

Run from the repository root:  py sim_r22/run_python.py
"""
import json
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path

import numpy as np

from r22_metrics import compute
from r22_switched import R22Params, simulate

OUT = Path(__file__).resolve().parent / "results"
SCHEMES = ("SVPWM", "DPWM1", "GDPWM")
LOAD = {0: 80.0, 15: 83.0, 45: 113.0}


def cases():
    c = [dict(key="apd_off", td_ns=200, kw=dict(scheme="GDPWM", apd=False, phi_deg=0, R_load=80.0))]
    for phi, R in LOAD.items():
        for s in SCHEMES:
            c.append(dict(key=f"rect_phi{phi}_{s}", td_ns=200, kw=dict(scheme=s, phi_deg=phi, R_load=R)))
    for iq in (5, 10):
        for s in SCHEMES:
            c.append(dict(key=f"statcom_{iq}A_{s}", td_ns=200,
                          kw=dict(scheme=s, phi_deg=90, R_load=None, I_q_rms=iq)))
    for td in (100, 300):
        for s in SCHEMES:
            c.append(dict(key=f"rect_phi0_{s}_td{td}", td_ns=td, kw=dict(scheme=s, phi_deg=0, R_load=80.0)))
    return c


def run(case):
    p = replace(R22Params(), td=case["td_ns"] * 1e-9)
    o = simulate(p, **case["kw"])
    m = compute(o)
    wave = None
    if case["key"] in ("rect_phi0_GDPWM", "rect_phi0_DPWM1", "apd_off"):
        keep = slice(-int(round(0.05 / o["meta"]["Tsw"])), None)    # last 50 ms
        wave = {k: o[k][keep].tolist() for k in ("t", "i_a", "i_c", "u_cac", "U_avg", "dA")}
    return case["key"], m, wave


def main():
    OUT.mkdir(exist_ok=True)
    with ProcessPoolExecutor(max_workers=max(1, min(8, (os.cpu_count() or 2) - 2))) as ex:
        res = list(ex.map(run, cases()))
    metrics = {k: m for k, m, _ in res}
    waves = {k: w for k, _, w in res if w is not None}
    (OUT / "python_metrics.json").write_text(json.dumps(metrics, indent=2, default=float))
    (OUT / "python_waveforms.json").write_text(json.dumps(waves))
    for k, m in metrics.items():
        print(f"{k:26s} THD {m['thd_ia_pct']:5.2f} %  ripple {m['dc_ripple_pp_V']:6.1f} V  "
              f"u_cac {m['amp_u_cac']:6.1f} V  clamp {m['Sa_clamp_mean_ms']}")


if __name__ == "__main__":
    main()
