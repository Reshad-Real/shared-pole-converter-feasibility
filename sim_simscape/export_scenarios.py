"""Export the shared scenario data that the Simscape model starts from.

Both implementations must begin from the same state and use the same orbit-derived
references, so these come from the Python orbit calculation; the circuit and the controller
are implemented separately (see CONTROLLER_SPEC.md).

Run from the repository root:  py sim_simscape/export_scenarios.py
"""
import json
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "sim"))

from model import Params, ac_state, orbit_energy_at        # noqa: E402
from feasibility import nominal_orbit, solve                # noqa: E402
from switched import SwitchParams                           # noqa: E402
from timedomain import GAINS, _post_orbit                   # noqa: E402

OUT = Path(__file__).resolve().parent / "scenarios.json"
PRE_CYCLES = 1


def main():
    R = json.loads((ROOT / "sim" / "results" / "results.json").read_text())
    p = Params(**R["params"])
    sw = SwitchParams()
    dP = R["pair"]["dP_test_W"]
    th_a, th_b = R["pair"]["theta_a_deg"], R["pair"]["theta_b_deg"]
    k_rec_hz = R["timedomain_tuning"]["k_rec_hz_chosen"]
    sweep = [c["phase_deg"] for c in R["controller_vs_lp"]]
    Tsw = 1.0 / sw.fsw
    M = int(round(0.5 * p.T / Tsw))
    pre, post = _post_orbit(p, dP)

    phases = sorted(set(sweep) | {th_a, th_b})
    scen = {}
    for th in phases:
        t0 = np.deg2rad(th) / p.w                     # absolute time of the simulation start
        Eb0 = orbit_energy_at(pre, p, t0)
        vb0 = float(np.sqrt(2 * Eb0 / p.Cb))
        pb0 = float(np.interp(np.mod(t0, p.T), pre["t"], pre["pb"], period=p.T))
        _, _, i1, i2, _ = ac_state(p, t0)
        hist = [orbit_energy_at(pre, p, t0 - (M - j) * Tsw) for j in range(M)]
        scen[f"{th:g}"] = dict(theta0_deg=th, t0=t0, i1=float(i1), i2=float(i2), ib=pb0 / vb0,
                               vb=vb0, Vdc=p.Vdc_nom, pd=pre["Pd"], hist=hist)

    orb = nominal_orbit(p)
    plan = solve(p, np.deg2rad(th_a), dP, "inner", orbit=orb, smooth=True)
    data = dict(
        params=asdict(p), switch=asdict(sw), dP=dP, theta_a_deg=th_a, theta_b_deg=th_b,
        sweep_deg=sweep, pre_cycles=PRE_CYCLES, M=M,
        gains=dict(wc=GAINS["wc"], k_bus=GAINS["k_bus"], k_rec=2 * np.pi * k_rec_hz,
                   k_tot=GAINS["k_tot"]),
        P_pre=pre["Pd"], P_post=post["Pd"], Eb_ref_pre=float(pre["Eb"].mean()),
        Eb_ref_post=float(post["Eb"].mean()),
        plan=dict(theta_deg=th_a, Ts=p.Ts, t=plan["t"].tolist(), Eb=plan["Eb"].tolist(),
                  Ed=plan["Ed"].tolist(), pd=plan["pd"].tolist()),
        scenarios=scen)
    OUT.write_text(json.dumps(data))
    print(f"wrote {OUT}: {len(scen)} phases, plan nodes {len(plan['t'])}, M={M}")


if __name__ == "__main__":
    main()
