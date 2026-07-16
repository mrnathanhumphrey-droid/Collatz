# Result 43 (qx+1 paper) — the corrected Chang↔Siegel edge FITS at single-step resolution: Chang's mod-8 return-class law = the generator's depth-3 2-adic truncation. The platform is a PIPELINE (Chang→Siegel→Nathan), not a tiling.

**Date:** 2026-07-16. **Verdicts: [pre-check] AXIS CORRECTION (Chang index K = one-step run-length; Nathan index i = q-adic level — orthogonal, opposite (p,q) sides) / ★ H_EDGE CONFIRMED (generator's single-step law reproduces Chang's mod-8 invariants exactly) / H_SIEGEL: α_H carries the branch structure / H_MAPBAL: finer (Appendix B), flagged not faked.**

**Headline: the platform brief posited Chang (finite 2-adic counts) + Siegel (2-adic Fourier) + Nathan (3-adic 2nd moment) as three MARGINALS of one generator matched index-to-index. The read-only pre-check found that mis-coordinatizes two orthogonal axes (Chang's run-length K vs Nathan's q-adic level i) on opposite (p,q) sides — the exact §6 false-pass trap. Corrected shape: a PIPELINE — Chang (2-adic INPUT: one-step law) → Siegel (TRANSFORM: numen χ_H: Z_2→Z_q) → Nathan (q-adic OUTPUT: level-i 2nd moment, r_q). Edge 1 (Chang↔Siegel), fired at the corrected single-step coordinate, FITS: Chang's mod-8 return-class law is reproduced EXACTLY from the generator's 2^{-k} halving weights × the ord_8(3)=2 mod-8 persistence selector (Pr[persistent]=1/4). So Chang's finite 2-adic object IS the generator's depth-3 2-adic truncation, and Siegel's α_H carries the branch structure. Pipeline first leg holds.**

Probe: `probe_43_chang_siegel_edge.py`. Log: `result_43_chang_siegel_edge_log.txt`. Runtime: instant.

## The pre-check (read-only, per brief §6) — axis correction

The brief's H_TILE matched Chang's `gap-K` (K=3,4,5) to Nathan's `S_0(i)` (i=3,4,5) index-to-index, and the "K=5 Map-Balance defect" to "i=5". But:
- **Chang's index K = `v_2(n+1)` = odd-run length of ONE step** (a single 2-adic valuation; law `2^{−k}`, verified = our `p_v`), refined by `μ = m mod 8`. It is a **single-coordinate VALUE** on the **2-adic domain**, first-moment (counts).
- **Nathan's index i = q-adic LEVEL = number of address coordinates** (each an independent halving depth). It is a **count of coordinates** on the **q-adic codomain**, second moment.

These are **orthogonal axes on opposite (p,q) sides** — matching them manufactures the false pass §6 warns of. The corrected structure is a **pipeline**, not a tiling: Chang (2-adic input statistics) → Siegel (numen transform Z_2→Z_q) → Nathan (q-adic output second moment). And this is already concrete: `r_q` **is** the pair-correlation of the `2^{−v}`-weighted address measure (build_M) — Chang's single-step input, carried through the transform.

## Edge 1 (Chang↔Siegel), corrected coordinate — the test

Both are single-step, first-moment, 2-adic-domain. Test: does the generator's single-step law (`2^{−k}` weights, multiplier 3, tracked mod 8) reproduce Chang's mod-8 invariants (Def 2.4: state `(k,μ)`, persistent iff `3^k·μ ≡ 7 mod 8`)?

| check | result |
|---|---|
| (i) `3^k mod 8` cycle | `3,1,3,1,…` ⇒ **`ord_8(3)=2`** (Chang's engine); our side `ord_3(2)=2` |
| (ii) persistence `3^k·μ≡7 mod 8` | selects **exactly one** `μ∈{1,3,5,7}` per k (k even→μ=7, k odd→μ=5) |
| (iii) `Pr[persistent | k]` | `1/4` for every k |
| (iv) `Pr[persistent]` | `Σ_k 2^{−k}·(1/4) = 1/4` — **exactly Chang's value**, derived his way |

**H_EDGE CONFIRMED.** Chang's mod-8 return-class law is the generator's depth-3 (`mod 2^3`) 2-adic truncation: the `2^{−k}` halving law is the single-coordinate weight, and the `ord_8(3)=2` period-2 cycling is the mod-8 branch selector. Nothing beyond the generator's single step is needed.

## H_SIEGEL — α_H carries the branch structure

`α_H(t) = ½(½ + (3/2)e^{−2πit})` (q=3, R41). Evaluated at 8th-roots `t=n/8`, it is the **mod-2 (2-branch) symbol**; the **mod-8 branch structure is its depth-3 iterate** (Siegel's finite-product Fourier shells — the numen at 2-adic depth 3). So Siegel's analytic symbol is the carrier of the mod-8 combinatorics Chang counts. Chang (counts) and Siegel (symbol) are the finite and analytic versions of the same single-step 2-adic object.

## The (2,3) order reciprocity — the shared arithmetic

Chang's mod-8 face is governed by **`ord_8(3)=2`**; our q-adic face by **`ord_3(2)=2`**. Both are order-2 facts of the *reciprocal* prime (`2²−1=3`, `3²−1=8=2³`). This is the same `⟨2⟩`/`⟨3⟩` order structure seen from the two adic sides — the arithmetic spine the pipeline runs on. Reported as the shared coordinate, not over-claimed as a deep reciprocity theorem (two order-2 facts at the specific (2,3)).

## H_MAPBAL — held (guard §6)

Chang's **Map Balance** (`#{gap-start ≡3} − #{≡7 mod 8} = exactly 1` for K≥5) is Appendix-B **gap-word / burst-gap** combinatorics — it counts odd-run *words* and their gap-start residues, a **multi-step refinement** finer than a single-step marginal. Reconstructing his exact word-definition from the summary risks a category slip (§6: the most likely false-pass route). **Not computed.** Flagged as the finer object to verify separately with Chang's precise Appendix-B definitions. The single-step edge is established; the exact-1 imbalance is a distinct, later check.

## What this places (the sudoku move)

- **Edge 1 fits ⇒ the pipeline's first leg (Chang→Siegel) holds.** Chang's finite 2-adic corner and Siegel's 2-adic Fourier body are the finite and analytic versions of ONE single-step generator.
- **Nathan places by composition, not by parallel marginal.** Our `r_q` is the second moment of the q-adic OUTPUT of that generator through Siegel's transform — the pipeline's last leg. So the three frameworks share a generator, but as **input→transform→output**, and the correct claim is "Chang, Siegel, Nathan are three stages of one pipeline on the (2,3)-adic Syracuse object," NOT "three marginals of one sequence."
- **The brief's §0 intuition is vindicated, its §3 mechanics corrected:** Siegel is the center piece (the transform) touching both — which is why he flagged the decay (Nathan's output) open while building the transform (Chang's input carrier).

## Not at stake
R1–R42. This is a framework-placement result; it changes no `r_q` value or L3 statement. `r_q` remains the q-adic output second-moment rate (build_M / R42).

_Reporting discipline: the axis mis-coordinatization was caught by the mandated §6 pre-check (read-only) BEFORE firing, not after a false pass. H_EDGE is exact set-equality (mod-8 invariants reproduced exactly, no fit). H_MAPBAL was explicitly NOT computed rather than faked from an uncertain definition — the guard against manufacturing a fit. The (2,3) reciprocity is flagged as observation, not claimed as theorem. Nathan's "same object" prior remains humble (0-for-9): the fit is a PIPELINE composition, not the clean tiling H_TILE posited._
