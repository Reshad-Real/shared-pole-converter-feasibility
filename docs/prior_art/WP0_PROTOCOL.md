# WP0: prior-art full-text check

Fixed on 22 October 2026, before any of these papers was read in full. Task (plan section 4.1): read the full texts of R13, R17, R18 and R20–R22, compare equations with R2, R4, R7 and R16, and confirm that no prior work already gives the phase-dependent transient criterion. If one does, the pivot rule of review section 11 applies.

## The claims under test

The paper would claim the following. Each is scored separately for every source.

| ID | Claim | What would count as already published |
|---|---|---|
| C1 | In a converter whose buffer shares a pole with AC ports, whether a transient is admissible depends on the AC phase at which it starts, not only on stored energy and port powers | A statement or result that admissibility, a capability limit or a feasible region changes with the line phase at which an event occurs |
| C2 | The test is a trajectory criterion over a horizon, with a slew-limited DC source, not an instantaneous or steady-state band check | A feasibility condition that follows the buffer and bus states through a transient, with a DC-side rate or bandwidth limit |
| C3 | Inner (sufficient) and outer (necessary) certificates, so infeasibility can be proved for any controller | Any necessary condition proved for all controllers, or a reachable or viable set, for a decoupling buffer |
| C4 | Two states with equal stored energy and equal port powers but different admissibility; energy sizing over-admits steps and under-sizes the buffer | An example or theorem showing that the energy-based minimum (for example R17's) or the energy buffering ratio (R18) is not sufficient during transients, depending on phase |
| C5 | The effect needs a DC source that is slow relative to the line cycle, and vanishes with a fast one | A statement of when transient buffer limits matter, in terms of DC-side bandwidth |

**Background that does not count against the claims:** the pole-voltage interval and its free offset (R20, R22), steady-state minimum energy (R17), usable-energy ratios for fixed voltage windows (R18), steady-state feasible regions for capacitor design (R4), ripple-dependent power limits (R7), feasible sliding-motion regions (R16). The paper must cite these and state its difference from each.

## Decision rule

- **Pivot** (review section 11) if any source gives C1 together with C2 for a shared-pole or shared-leg buffer. Then the criterion is not new, and the paper needs a stronger result for a different architecture or uncertainty setting.
- **Narrow the claim** if a source gives C3 or C4 for a related architecture (for example a dedicated buffer leg). Then that part becomes an extension, cited as such.
- **Proceed** if no source gives C1 with C2; record for each source the closest result and the sentence that distinguishes ours.

## What to extract from each source

For each paper, record in `WP0_RESULTS.md` (written after reading):

1. Bibliographic record, confirmed from the PDF (title, authors, venue, volume, pages, DOI).
2. Architecture: which legs are shared, where the buffer sits, which ports are AC or DC.
3. The constraint or limit it derives, with equation numbers: pole-voltage interval, energy minimum, buffering ratio, feasible region, power limit.
4. Whether the result is instantaneous, steady-state (one line cycle) or transient (a horizon after an event), and whether it depends on the line phase.
5. Whether the DC side is ideal, slew-limited or bandwidth-limited.
6. Any infeasibility statement for all controllers, or any reachability or viability argument.
7. A score of Y (gives it), P (partly) or N (does not) for C1–C5, each with a page or equation reference.
8. The sentence our paper will use to cite it and state the difference.

## Reading list

Suggested filenames are for this folder. IEEE, IET and Wiley copies obtained through an institution are licensed to it: do not include them in the public code release (WP7).

### Needs institutional access (IEEE Xplore or Wiley)

| Ref | Paper | DOI | Suggested filename |
|---|---|---|---|
| R17 | P. T. Krein, R. S. Balog, M. Mirjafari, "Minimum Energy and Capacitance Requirements for Single-Phase Inverters and Rectifiers Using a Ripple Port," IEEE TPEL 27(11):4690–4698, 2012 | [10.1109/TPEL.2012.2186640](https://doi.org/10.1109/TPEL.2012.2186640) | `R17_Krein_2012_TPEL.pdf` |
| R20a | T. Kominami, Y. Fujimoto, "A Novel Nine-Switch Inverter for Independent Control of Two Three-Phase Loads," IEEE IAS Annual Meeting 2007, pp. 2346–2350 | [10.1109/07IAS.2007.354](https://doi.org/10.1109/07IAS.2007.354) | `R20a_Kominami_2007_IAS.pdf` |
| R20b | C. Liu, B. Wu, N. R. Zargari, D. Xu, J. Wang, "A Novel Three-Phase Three-Leg AC/AC Converter Using Nine IGBTs," IEEE TPEL 24(5):1151–1160, 2009 | [10.1109/TPEL.2008.2004038](https://doi.org/10.1109/TPEL.2008.2004038) | `R20b_Liu_2009_TPEL.pdf` |
| R21 | A. E. L. da Costa, C. B. Jacobina, N. Rocha, E. R. C. da Silva, A. V. M. Lacerda Filho, "A Single-Phase ac–dc–ac Unidirectional Three-Leg Converter," IEEE TIE 68(5):3876–3886, 2021 | [10.1109/TIE.2020.2987261](https://doi.org/10.1109/TIE.2020.2987261) | `R21_daCosta_2021_TIE.pdf` |
| R7 | G. Liang, H. D. Tafti, G. G. Farivar, J. Pou, C. D. Townsend, G. Konstantinou, S. Ceballos, "Effect of Capacitor Voltage Ripples on Submodule Active Power Control Limits of Cascaded Multilevel Converters," IEEE TIE 69(6):5952–5961, 2022 | [10.1109/TIE.2021.3091934](https://doi.org/10.1109/TIE.2021.3091934) | `R7_Liang_2022_TIE.pdf` |

No open copy was found for any of these. The R7 repository records (UWA, Tecnalia) give the DOI only or are marked closed access.

### Open copies (downloaded 22 October 2026)

| Ref | Paper | Source | File in this folder |
|---|---|---|---|
| R18 | M. Chen, K. K. Afridi, D. J. Perreault, "Stacked Switched Capacitor Energy Buffer Architecture," IEEE TPEL 28(11):5183–5195, 2013, DOI 10.1109/TPEL.2013.2245682 | [author copy (Princeton)](https://swh.princeton.edu/~minjie/files/chen-sscbuffer-2013.pdf); the authors' manuscript, 14 pages, without the journal's page numbers | `R18_Chen_2013_TPEL.pdf` (2.2 MB) |
| R2 | M. J. Heller, F. Krismer, J. W. Kolar, "Duty-Cycle Dependent Phase Shift Modulation of Dual Three-Phase Active Bridge Four-Port AC–DC/DC–AC Converter …," IEEE OJPEL 3:705–722, 2022 (open access) | [ETH copy](https://www.ams-publications.ee.ethz.ch/uploads/tx_ethpublications/10_Duty-Cycle_Dependent_Phase_FINAL_Heller.pdf); the published version, 18 pages | `R2_Heller_2022_OJPEL.pdf` (9.4 MB) |
| R16 | S. Hwang, "Study on optimization and robust control of active power decoupling circuit based converters," doctoral thesis, Hokkaido University, 2019, DOI 10.14943/doctoral.k13531 | [HUSCAP](https://eprints.lib.hokudai.ac.jp/repo/huscap/all/74332/Hwang_Soonhwan.pdf); 103 pages | `R16_Hwang_2019_thesis.pdf` (18.3 MB) |
| R20c | T. Kominami, Y. Fujimoto, "Proposal of a Nine-Switch Inverter that can Independently Control Amplitude and Frequency for Two Loads," IEEJ Trans. IA 128(5):561–568, 2008, DOI 10.1541/ieejias.128.561 | [J-STAGE, free access](https://www.jstage.jst.go.jp/article/ieejias/128/5/128_5_561/_article/-char/en); 9 pages. Only the one-page extended summary is in English; the main text is Japanese, though its equations and figures are usable. R20a remains the primary English source | `R20c_Kominami_2008_IEEJ.pdf` (1.4 MB) |

### Open copies whose servers did not respond from here

| Ref | Paper | Source | Suggested filename |
|---|---|---|---|
| R4 | Y. Liu, Y. Sun, M. Su, "Active power compensation method for single-phase current source rectifier without extra active switches," IET Power Electronics 9(8):1719–1726, 2016, DOI 10.1049/iet-pel.2015.0899 | [authors' lab page (CSU)](https://faculty.csu.edu.cn/sunyao/zh_CN/lwcg/3808/content/25949.htm), or through Wiley | `R4_Liu_2016_IETPE.pdf` |
| R21 (related) | B. S. Gehrke, master's dissertation on shared-leg converters, UFCG, 2019 (Portuguese) | [UFCG repository](https://bdtd.ufcg.edu.br/bitstream/riufcg/29361/1/BRUNA%20SEIBEL%20GEHRKE%20-%20DISSERTA%c3%87%c3%83O%20PPGEE%202019.pdf) | `R21x_Gehrke_2019_UFCG.pdf` |

### Not yet obtainable

- **R13**, "Optimized Power Decoupling for a Three-Leg Single-Phase AC/AC Converter: Reduced DC Bus Ripple with Minimal Capacitance," ECCE 2026, paper 2424.
  - ECCE 2026 was held in Vancouver on 4–8 October 2026. The paper is not yet in Crossref or IEEE Xplore.
  - The program page now needs an ECCE login, so the authors could not be identified.
  - Options: the attendee proceedings, if you have conference access; the Xplore listing when the proceedings appear; or the authors, once identified.
  - This is the closest topical match. WP0 cannot close without it, so it is re-checked before submission.

### Already here

- **R22**: `Carrier-Based_Generalized_Discontinuous_PWM_Strategy_for_Single-Phase_Three-Legs_Active_Power_Decoupling_Converters.pdf` (read for WP3; still to be scored against C1–C5).

### Optional companion

- J. Yang, X. Li, M. Su, Y. Liu, Y. Sun, "Active Power Decoupling Method for Single-Phase Current-Source Rectifier With No Additional Active Switches," IEEE TPEL, 2016. A sibling of R4 by the same group, also linked from the [CSU lab page](https://faculty.csu.edu.cn/sunyao/zh_CN/lwcg/3808/content/25946.htm).
