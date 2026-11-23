# WP6: system case study

Fixed on 24 October 2026, before any WP6 run. The question is how often, over a day of operation of a residential hybrid inverter, a backup-load event arrives at a phase that conventional buffer sizing cannot handle, and what the proposed sizing changes.

## Revision before the campaign

A method test on 24 October 2026 timed the outcomes at off-grid phases (7°, 97°, 187° and 277°). It showed that WP5's sizing does not cover load drops:

- **The full drop is not certified.** The drop from full AC2 load to zero is not certified at D2's capacitance at any of the four phases, at either rating, even from the lowest operating point. The conventional controller trips on it with a bus overvoltage. At D2's capacitance the 120 V buffer band holds about 19 J at full load. DC1, slew-limited at the bus-loop rate, leaves about 24 J of surplus to absorb.
- **It certifies with more capacitance:** from about 1.5 times D2's capacitance at 120 V and 2 times at 230 V.
- **The other drops are not binding.** They are certified at D2's capacitance, apart from isolated phases that look like misses of the local search.
- **The necessary bound misses this.** The circuit outer bound calls the full drop feasible, because it is loose for surplus (WP5 S6).

The revision therefore adds:

- **Design D4, proposed sizing with load drops:** the smallest of 1.0, 1.25, 1.5, 1.75, 2.0, 2.5 and 3.0 times D2's capacitance at which the full drop is certified at all 24 grid phases, with the operating-point schedule. It is evaluated on every event like the other designs.
- **Criterion W3b:** the same as W3, for D4.
- **D2 is kept** as declared. W3 for D2 is expected to fail on the full drop (24 of 288 events).
- **For drops, the unavoidable-violation outcome is a weak lower bound.** It is reported, but the controllers' results carry the drop events.

## System and designs

- **Converter and ratings:** the two realistic ratings of `ratings.py` at the WP5 base design: m = 0.40, the bus-loop DC1 ramp and unity AC2 power factor.
  - AC2 is the backup (critical-load) port.
  - The other ports keep their rating values (AC1 1.2 kW or 2.4 kW, DC2 0.4 kW or 0.8 kW), and DC1 (the battery converter) supplies the balance at its WP5 slew rate.
- **Designs** (capacitances from WP5, `results/wp5a/`):

| ID | Design | Capacitance (120 V / 230 V) | Pre-step operating point |
|---|---|---|---|
| D1 | Conventional sizing | Cb_conv: 724 / 480 µF | centred (α = 0.5) at every load |
| D2 | Proposed sizing | Cb_suff,c at α = 1: 1076 / 601 µF | the load schedule below |
| D3 | Conventional capacitance, proposed operating point | Cb_conv: 724 / 480 µF | the load schedule below |

- **Operating-point schedule** (D2 and D3), by AC2 load before the event as a fraction of its rating:

| Load | 0 | 1/3 | 2/3 | 1 |
|---|---|---|---|---|
| α | 1.0 | 1.0 | 0.5 | 0.0 |

  - **Why this schedule:** a high orbit at light load admits step-ups, the large events there (WP5). At full load the only events are drops, which a low orbit admits.
  - After an event the controller recovers to the schedule's orbit at the new load.

## Events

- **Load levels:** 0, 1/3, 2/3 and 1 of the AC2 rating (0, 800, 1600, 2400 W at 120 V; twice that at 230 V).
- **Event:** an AC2 power step between two different levels (12 types: 6 step-ups, 6 drops), from the pre-step steady orbit of the design.
- **Phase:** the AC phase at which the event starts, on a 15° grid (24 phases). Appliance switching is not synchronised to the grid, so the base case takes every phase as equally likely.
- **Independence:** events are taken to be far enough apart (minutes) for the converter to return to steady state between them.

## Outcomes per design, event type and phase

1. **Unavoidable violation:** the circuit outer bound (`outer_circuit.py`, three post-step cycles) is infeasible from the design's pre-step orbit, so no controller can handle the event within the limits.
2. **Conventional controller:** the switched model (`switched.py`, default switching parameters) with the conventional controller holding the design's operating point. Recorded as held or tripped, with the first violation's time and type.
3. **Certified controller** (D2 and D4):
   - the current-aware sufficient bound (`inner_circuit.py`, κ = 0.97, four cycles) certifies the event, and the smoothed-bus tracker of WP5 (`plan_smooth`, `sync_ib`, 15 samples) follows the plan in the switched model;
   - recorded as held, tripped, or not certified.

The fraction of failing phases gives each event type's failure probability at a uniformly random phase.

## Daily model (stylised)

- **Arrivals:** events arrive as a Poisson process at λ events per day. The base case is λ = 60 (backup-panel appliances switching on or off).
- **Next level:** from each level, the next is drawn with weight w(|Δ|) for a jump of |Δ| thirds of the rating, with w(1) = 1, w(2) = 0.25 and w(3) = 0.05. The chain is symmetric, so its stationary distribution follows.
- **Daily probability:** P_day = 1 − exp(−λ Σ π_i P_ij p_ij), with
  - π the stationary level distribution,
  - P_ij the next-level probability,
  - p_ij the failure probability of event type i → j.
- **Sensitivity:**
  - λ = 20 and 150;
  - large-jump weights w(2) and w(3) set to zero, and doubled;
  - zero-crossing switching: all events at 0° or 180°, half each.
- **The model is stylised.** The per-event-type failure probabilities are the primary result, so that any other load model can be applied to them.

## Criteria

| ID | Criterion | Pass mark |
|---|---|---|
| W1 | Unavoidable violations at conventional sizing (D1): P_day under the base model | ≥ 10 % at both ratings |
| W2 | Proposed sizing (D2): no event type and phase is infeasible by the necessary bound | at both ratings |
| W3 | Proposed sizing (D2) with the certified controller: share of the 288 events (12 types × 24 phases) that are certified and held | ≥ 95 % at both ratings |
| W3b | The same for D4 (proposed sizing with load drops) | ≥ 95 % at both ratings |
| W4 | Conventional controller trip probability per day, D1 against D2 | reported |
| W5 | Operating point alone (D3 against D1): unavoidable and conventional-trip P_day | reported |
| W6 | Cost of the proposed sizing: extra capacitance, and extra stored energy at the bus maximum voltage | reported |

- **Consistency check:** every event that the conventional or certified controller holds must be feasible by the outer bound. A held event that the bound calls infeasible would contradict the bound, and is counted and reported.
- **Not covered:** tolerances (WP4), other ratings and modulation indices, appliance inrush currents, non-unity power factor, and events closer together than the recovery time.
