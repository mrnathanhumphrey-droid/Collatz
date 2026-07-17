# Result 64.B — The 1:4 class-mass ratio (provenance reconstruction)

**Date:** 2026-07-14. **Status: CONFIRMED — elementary and exact.** Written to close a provenance gap: R64.B is cited as load-bearing for c=7/45 in `THEOREM_C_745.md` (H64.B), `D3_DERIVATION_AUDIT.md` (Probe 4.5, step 4), and `result_77_T_lead_spectrum.md` §1, but had **no result file**. This is that file.

Probe: `result_64B_verify.py`; log `result_64B_log.txt`.

## Statement

The R77 T_diag leading eigenvector is **(1, 4)**, and the asymptotic class-resolved Plancherel masses satisfy `P_+ : P_− = 1 : 4`. This ratio is the **squared class-mass ratio** `(1/3)² : (2/3)² = 1 : 4`, where 1/3 and 2/3 are the class probabilities. It is an elementary parity identity of the Geom(½) halving law — not a deep or separately-conjectural input.

## Class definition (from R77 §, R66 chain rule)

`result_77_T_lead_spectrum.md` line 14: the R66 chain rule assigns **class + = v even, class − = v odd**, where `v` is the 2-adic valuation of the halving step (`v ~ Geom(½)`, `P(v=k)=2^{−k}`, k≥1). The class-resolved characteristic functions are `μ̂^+ = Σ_{v even}…`, `μ̂^− = Σ_{v odd}…` (R77 lines 15–16).

## Derivation (two equivalent elementary routes, both exact)

**(a) Class probabilities → squared ratio.** Parity split of `Σ_{v≥1} 2^{−v}`:
- `P(v even) = Σ_{v even≥2} 2^{−v} = (1/4)/(1−1/4) = 1/3`
- `P(v odd)  = Σ_{v odd≥1} 2^{−v} = (1/2)/(1−1/4) = 2/3`

The Plancherel mass is `|μ̂|²`, quadratic in the class amplitude, so `P_± ∝ (class probability)²`, giving `P_+ : P_− = (1/3)² : (2/3)² = **1 : 4**`. This is the (1,4) eigenvector direction (R77 line 119).

**(b) Plancherel weight split → T_diag prefactors.** Each step contributes `2^{−v}` to `μ̂`, hence `4^{−v}` to `|μ̂|²`. Parity split of `Σ_{v≥1} 4^{−v}`:
- `Σ_{v even} 4^{−v} = (1/16)/(1−1/16) = **1/15**`
- `Σ_{v odd}  4^{−v} = (1/4)/(1−1/16) = **4/15**`

Ratio `1/15 : 4/15 = **1 : 4**`. Times the Plancherel factor 3 these give the T_diag prefactors `P^{++}` coeff `= 3·(1/15) = 1/5` and `P^{−−}` coeff `= 3·(4/15) = 4/5`, i.e. `T_diag = (1/5)·[[1,1],[4,4]]` (R77 Thm 77.1) — eigenvalues {1,0}, λ=1 eigenvector (1,4). ✓ verified.

## Numerical confirmation (real Collatz values)

`v = v₂(3n+1)` over 40M uniform odd, coprime-to-3 starts (invariant single-step marginal): `P(v even) = 0.33333`, `P(v odd) = 0.66667`; empirical squared-class-mass ratio `(p_odd)²:(p_even)² = 4.0000 : 1`. Matches exactly.

## Disposition

R64.B is **elementary and exact** — the parity split of Geom(½). c=7/45's dependence on it (via R77's (1,4) eigenvector) rests on this one-line identity, not on any unverified input. The gap flagged by the 2026-07-14 thread audit was **write-up only** (no result file existed); the value was always correct. `D3_DERIVATION_AUDIT.md`'s treatment of R64.B as a rigorous pillar stands.

_Reporting discipline: written to supply missing provenance for a load-bearing constant; the value is confirmed, not merely asserted. A self-check bug (double-counting the halving weight in an empirical 4^{−v} sub-test) was caught and corrected before the verdict — the surviving checks are all exact-rational or exact-eigenvector._
