# NISOLI_CLOSURE_CORRECTED_DISPOSITION — top-level disposition

**Date:** 2026-05-12. Route A of the T_lead Nisoli-bypass re-evaluation. Wilson reporting to Nathan.

---

## DISPOSITION: **H_A_EXTRACTION_HARD**

> **Headline.** The Nisoli closure inequality at the corrected rate λ = 43/45 fires mathematically at the smallest plausible cell (r=3, A=3, K=6) with product `|K_bil(3)| · M_3'' · K^{−A} = 6 · 24.43 · 1/216 = 0.679 < 1` — a clean, sub-1, single-digit-A, single-digit-K firing cell. **HOWEVER**, the A required is not extractable from any project source: Tao Prop 1.17's effective `C_A` is INFEASIBLE per BOOKKEEPING_PHASE1_DISPOSITION (super-exponential A-dependence under any bookkeeping), and the cross-freq machinery / T_lead structure / V_M closure deliver structurally different decay objects that don't substitute for Tao's A. So the closure inequality is **mathematically parameterised** (the firing cell exists) but **rigorously unevaluable** (no usable A delivery).
>
> **The corrected-rate Nisoli closure is the first POSITIVE STRUCTURAL ANCHOR** for c = 7/45 from the eight-probe trajectory (T_lead's exact rational eigenvalue 43/45 + matching tight M_3'' = 24.43), but it is **NOT yet a first rigorous spectral closure** — the Tao-A obstruction is rate-invariant and persists.

---

## Pre-registered hypotheses, decided

| Hypothesis | Status |
|---|---|
| H_NISOLI_CLOSES_AT_R3 (closure inequality fires at r=3) | **CONDITIONAL YES** — fires at (r=3, A=3, K=6) under `C_A = 1`; product 0.679. NOT a rigorous closure because the A delivery is missing. |
| H_NISOLI_CLOSES_AT_LARGER_R (fires at r ∈ {4, 5}) | **CONDITIONAL YES** — fires at (r=4, A=3, K=8) product 0.859; (r=5, A=3, K=10) product 0.762. Same A-extraction caveat. |
| H_NISOLI_NEEDS_TIGHTENING | **NO** — M_3'' = 24.43 is exact (closed-form 2×2 SVD); not loose. |
| H_NISOLI_DOESNT_CLOSE | **NO** — closure fires at many cells under `C_A = 1`. |
| **H_A_EXTRACTION_HARD** | **CHOSEN** — A not extractable from Tao (BOOKKEEPING_PHASE1 INFEASIBLE), nor from cross-freq (different decay object), nor from T_lead's null eigenvalue (within-level, not asymptotic). |

The pre-registration favoured H_NISOLI_NEEDS_TIGHTENING. The actual outcome is more structural: tightening isn't needed (M_3'' is already exact); the obstacle is the analytic prerequisite (Tao's A), which is rate-invariant.

---

## Phase summary

### Phase 1: closure inequality articulated; M_3'' verified

`(I − T_lead)⁻¹ = (1/2)·[[9, 9], [28, 38]]`. A^T A = (1/4)·[[865, 1145], [1145, 1525]]. trace = 597.5, det = 506.25 (exact). σ²_max = 596.6515, **M_3'' = σ_max ≈ 24.426**. Matches T_LEAD_CORRECTED_CLOSURE.md.

Closure inequality: `|K_bil(r)| · K^{−A} · M_3'' < 1`.

### Phase 2-3: (r, A, K) tabulation; firing cells

Thresholds `|K_bil(r)| · M_3''`:

| r | |K_bil(r)| | threshold |
|---|---|---|
| 3 | 6 (strict 2√N) | 146.6 |
| 4 | 18 (polylog-free 2√p·√N) | 439.7 |
| 5 | 31.18 | 761.7 |

Smallest firing cells under `C_A = 1`:

| Cell | product | fires |
|---|---|---|
| (r=3, A=3, K=6) | 0.679 | YES — smallest K |
| (r=3, A=2, K=13) | 0.867 | YES — smallest A=2 |
| (r=4, A=3, K=8) | 0.859 | YES |
| (r=5, A=3, K=10) | 0.762 | YES |

### Phase 4: honest read on A

Four sources of A:

- (a) Tao Prop 1.17 effective `C_A`: **INFEASIBLE** (BOOKKEEPING_PHASE1; iterated-cubic recursion forces super-exponential A-dependence).
- (b) Cross-freq machinery: delivers geometric decay in `g` (cross-frequency weights), not in `n` (Tao's decay). Structurally different object.
- (c) T_lead's null eigenvalue: within-level rank-1 collapse on (9, −7), not the asymptotic across-n decay Nisoli needs.
- (d) Novel polynomial-in-A bound outside Tao's method: research project; not in scope.

None deliver A. Closure inequality is unevaluable rigorously.

### Phase 5: tightening M_3''

Skipped — M_3'' = 24.43 is already exact (closed-form 2×2 SVD, exact rational A^T A). No tightening available.

### Phase 6: structural meaning for seven-probe trajectory

The seven-probe "Nisoli structurally inapplicable at rate-1/2" framing combined TWO obstructions:

1. **Operator-theoretic** (rate-specific): rate-1/2 has no clean Q-rational eigenvalue (R77.3 falsified). **LIFTS at corrected rate** (T_lead's 43/45 exists).
2. **Analytic / Tao-A** (rate-invariant): no polynomial-in-A `C_A` from Tao's method. **PERSISTS at corrected rate**.

The eighth probe **separates these**. The operator-theoretic obstruction is rate-specific; the Tao-A obstruction is rate-invariant. The corrected rate improves Item 1, doesn't improve Item 2.

The seven-probe bequest "Nisoli inapplicable" reframes to: **"Nisoli's analytic prerequisite (polynomial-in-A Fourier bound) is structurally unavailable via Tao's renewal-process method, at any rate."**

---

## Comparison: rate-1/2 vs corrected rate, side-by-side

| Quantity | rate-1/2 (R77.3, M3_DISPOSITION) | corrected rate 43/45 (this probe) |
|---|---|---|
| Discrete eigenvalue at rate exists? | NO (3-mode falsified) | **YES (43/45 over Q exact)** |
| Operator-theoretic anchor? | NO (T uncharacterised) | **YES (T_lead = (1/45)·[[7,9],[28,36]])** |
| M_3 exact value | N/A (no operator); loose 800-1000 vs falsified T_3 | **24.43 exact (closed-form SVD)** |
| Smallest firing cell under `C_A = 1` | (r=3, A=3, K=6) if M_3=100, product 0.926 | (r=3, A=3, K=6), product **0.679** |
| Tao-A `C_A` extractable rigorously? | NO (BOOKKEEPING_PHASE1) | NO (same obstruction) |
| Rigorous closure achievable? | NO | NO |
| Structural advance? | partial (operator falsified) | **POSITIVE (T_lead + M_3'' both clean)** |

The corrected rate **structurally lifts** the rate-1/2 operator-theoretic block. It does **not** lift the rate-invariant Tao-A block.

---

## Honest read on what's a "first rigorous spectral closure"

A naive reading of the firing cell (r=3, A=3, K=6) would say: "the corrected-rate Nisoli closure of c = 7/45 is achieved." This is **not honest** because:

1. The `C_A = 1` parameterisation is not what Tao Prop 1.17 delivers. Tao's actual `C_A` is super-exponential in A (BOOKKEEPING_PHASE1's optimistic floor `A^{O(A)}`, faithful bookkeeping `exp(exp(A²))`).

2. Substituting any realistic `C_A`: at A=3, optimistic-floor `C_A ≥ 3^3 = 27`, so the inequality becomes `6 · 27 · 24.43 · 6^{−3} = 3958/216 ≈ 18.3 > 1`. **Fails by a factor of 18.**

3. Faithful bookkeeping: `C_3 ≥ exp(exp(9)) ≈ 10^{3500}`. Astronomical failure.

So the **honest claim** is:

> **The corrected-rate Nisoli closure inequality is FIRING CELL READY** (a single-digit (A, K, r) cell exists with sub-1 product under `C_A = 1`). **It is NOT closed** because the `C_A = 1` parameterisation is not delivered by Tao's method. The rigorous closure remains blocked by the same Tao-A obstruction that BOOKKEEPING_PHASE1 identified at rate-1/2 — but now with the operator-theoretic anchor (T_lead, M_3'' = 24.43) cleanly in hand.

---

## What's actually needed for closure now

The corrected-rate probe **reduces the c = 7/45 Nisoli closure roadmap from three independent obstructions down to one**:

**Pre-eighth-probe** (rate-1/2 setting):
1. Tao C_A (INFEASIBLE).
2. Bilinear bound (DELIVERED at family level).
3. Characterised operator T (uncharacterised; R77.3 falsified T_3).

**Post-eighth-probe** (corrected rate setting):
1. Tao C_A (INFEASIBLE; rate-invariant).
2. Bilinear bound (DELIVERED).
3. Characterised operator (DELIVERED; T_lead with M_3'' = 24.43).

**The remaining block is ONLY Item 1.** Resolving it would close c = 7/45. The path is:

- **Novel Fourier-analytic bound on `|μ̂_n(ξ)|` with polynomial-in-A constant**, outside Tao's renewal-process method. (Mentioned as the lone unblocker in BOOKKEEPING_PHASE1 §"Alternative routes.")

This is the load-bearing open question for c = 7/45 closure post-this-probe.

---

## Trajectory placement (eight probes + this Route-A test)

| Probe | Object | Disposition |
|---|---|---|
| T_3 (R77.3) | 3×3 companion at rate-1/2 | FALSIFIED |
| R_k | Inter-level residual | INTRACTABLE |
| Candidate A | W_k φ_n form | FALSIFIES_F2 |
| R76 §11 2D | T_diag + Off | INCONCLUSIVE |
| T_N | T_diag + Off_lin at rate-1/2 | UNDERSPECIFIED |
| Cross-freq closure | V_M closure space | CLOSES_ON_ENLARGED_SPAN |
| T_V | T_V on V_M at rate-1/2 | RECURSION_UNDERSPECIFIED |
| T_lead corrected (eighth) | T_lead at corrected rate | **DIFFERENT_RATE (43/45 over Q)** |
| **Nisoli closure at corrected rate (this)** | Closure inequality at λ = 43/45 | **H_A_EXTRACTION_HARD** |

This is the **ninth landing**. It establishes that the operator-theoretic anchor at the corrected rate (T_lead) supports a Nisoli framework setup with explicit M_3'' = 24.43 and a single-digit-(A, K) firing cell under `C_A = 1`. The blocker is **Tao-A**, which is rate-invariant.

---

## Adversarial checks

**(A1) M_3'' computation:** Verified via closed-form (I − T_lead)⁻¹ = (1/2)·[[9,9],[28,38]], A^T A = (1/4)·[[865,1145],[1145,1525]], trace = 597.5, det = 506.25 exact, σ²_max = 596.6515. σ_max = 24.426. **Matches T_LEAD_CORRECTED_CLOSURE.md §2.** ✓

**(A2) Bilinear bound:** PATH2's 2√N at r ≤ 3 strict (family-level, p-uniform; PATH2_DISPOSITION); HENSEL's 2√p·√N polylog-free at r ≥ 4 (Hensel-digit-extraction; HENSEL_DISPOSITION's H_HENSEL_CLOSES + H_HENSEL_PARTIAL_TIGHTENING). **Cited from memory `project_collatz_r78_bilinear_cracked` and the two disposition docs.** ✓

**(A3) Firing cell mathematics:** (r=3, A=3, K=6): |K_bil(3)| = 6, M_3'' = 24.43, K^{−A} = 1/216 = 0.00463. Product = 6 · 24.43 · 0.00463 = 0.679. < 1. ✓

**(A4) Tao-A extraction:** BOOKKEEPING_PHASE1's INFEASIBLE finding rests on the iterated-cubic recursion in Tao §7.4 Case 3 — a method-internal obstruction, not a bookkeeping-quality issue. Cross-checked against BOOKKEEPING_PHASE1's optimistic-floor `A^{O(A)}` and faithful-bookkeeping `exp(exp(A²))` bounds. **The obstruction is rate-invariant** (same Tao-A would be needed at any target rate). ✓

**(A5) "First rigorous spectral closure" claim:** I do NOT make this claim. The disposition is **H_A_EXTRACTION_HARD**, not H_NISOLI_CLOSES_AT_R3. The firing cell exists mathematically; rigorous closure requires Tao-A, which is unavailable. ✓ honestly flagged.

**(A6) Structural-vs-numerical:** All M_3'' work uses exact rationals (Fractions) for matrix entries; the σ_max computation is the only numerical step (and exact rational `σ²_max = 596.6515` is preserved before sqrt). All thresholds in §"Phase 2-3 firing cells" are derived from exact rationals; only the final products are reported as floats. ✓

---

## Recommendation (surfaced for Nathan, not chosen)

### Route A1: Document for paper (Recommended for low effort)

The structural finding (T_lead's exact eigenvalue 43/45 with closed-form origin, M_3'' = 24.43 exact, smallest firing cell at (r=3, A=3, K=6) under `C_A = 1`, Tao-A as the lone remaining rate-invariant obstruction) is **directly publishable** as the eighth-probe positive anchor. Combined with PATH2/HENSEL bilinear bound, this is a structurally complete spectral framework for c = 7/45 modulo the Tao-A input.

Estimated effort: 0.5 session.

### Route A2: Pursue novel polynomial-in-A bound (substantial effort)

The lone unblocker (per Phase 4 source (d) and BOOKKEEPING_PHASE1 §"Alternative routes") is a non-Tao Fourier bound on `|μ̂_n(ξ)|` with polynomial-in-A constant. This is a research project, not an extraction. Possibilities to explore: Kalafatelis-style cubic-character bounds, Bourgain-Konyagin sum-product on relevant groups, or novel framework integrating cross-freq's V_M closure with per-frequency Fourier bounds.

Estimated effort: 2-5 sessions to scout candidates; many more if any candidate path opens up.

### Route A3: Pivot to non-Nisoli rigor (alternative)

PRECISE_ASK §4 (c) flagged "phase-cancellation tighter than C_emp · √N" as an alternative to Nisoli. This would bypass the closure inequality entirely. Plausibility: empirical β = 0.522 already incorporates all empirical phase cancellation, so this route is unlikely to deliver more.

Estimated effort: 1-2 sessions to confirm unlikely closure.

---

## Synopsis (one paragraph for Nathan)

The ninth landing in the c = 7/45 spectral trajectory: at the corrected rate λ = 43/45 (T_lead's clean Q-rational eigenvalue), the Nisoli closure inequality `|K_bil(r)| · K^{−A} · M_3'' < 1` has resolvent norm `M_3'' = 24.43` exact (closed-form 2×2 SVD on (I − T_lead)⁻¹). Combined with PATH2's strict `|K_bil(r=3)| = 2√N = 6` bilinear bound, the inequality fires at `(r=3, A=3, K=6)` with product `6 · 24.43 / 216 = 0.679` — a clean sub-1 single-digit firing cell. **However**, Tao Prop 1.17's effective `C_A` is INFEASIBLE per BOOKKEEPING_PHASE1's 2026-05-11 disposition (iterated-cubic recursion forces super-exponential A-dependence regardless of bookkeeping quality), and the cross-freq / T_lead-null-eigenvalue / V_M-closure machinery delivers structurally different decay objects (decay in `g` for cross-freq, within-level rank-1 collapse for T_lead) that do NOT substitute for Tao's asymptotic A. So the closure is mathematically parameterised at the firing cell but **rigorously unevaluable**. The disposition lands **H_A_EXTRACTION_HARD**: the eighth probe's operator-theoretic anchor (T_lead) lifts the rate-specific obstruction that blocked the rate-1/2 probes, but the rate-invariant Tao-A obstruction persists. The corrected-rate roadmap now reduces from three obstructions to one — Tao-A as a novel polynomial-in-A Fourier bound outside Tao's renewal-process method. **The corrected rate finding stands as a positive structural anchor (43/45 over Q with M_3'' = 24.43 exact), but is NOT the first rigorous spectral closure of c = 7/45.**

---

## Deliverables (this Route A)

- NISOLI_CLOSURE_CORRECTED_PHASE1.md — closure inequality articulated; M_3'' verified
- NISOLI_CLOSURE_CORRECTED_TABLE.md — (r, A, K) tabulation
- NISOLI_CLOSURE_CORRECTED_FIRING.md — firing cells documented
- NISOLI_CLOSURE_CORRECTED_HONEST_READ.md — Phase 4 on A's extractability
- NISOLI_CLOSURE_CORRECTED_STRUCTURAL.md — Phase 6 implication for seven-probe trajectory
- NISOLI_CLOSURE_CORRECTED_DISPOSITION.md (this file) — top-level
- nisoli_closure_corrected.py — main-thread verification script (Phase 1 M_3'' + Phase 2-3 tabulation)

---

End of Route A Phase 6 disposition.
