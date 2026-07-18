# NISOLI_CLOSURE_CORRECTED_STRUCTURAL — Phase 6: implication for seven-probe trajectory

**Date:** 2026-05-12. Phase 6 of Route A.

---

## 1. The seven-probe trajectory's "structural inapplicability" framing

Per `T_LEAD_CORRECTED_DISPOSITION.md` §"Trajectory placement," the eight-probe trajectory at rate-1/2 lands:

| Probe | Object | Disposition |
|---|---|---|
| T_3 (R77.3) | 3x3 companion at rate-1/2 | FALSIFIED |
| R_k (R77.4 erratum §1) | Inter-level residual operator | H_R_K_INTRACTABLE |
| Candidate A | W_k φ_n bilinear-pair-form | FALSIFIES_F2 |
| R76 §11 2D | T_diag + Off conjectural at rate-1/2 | INCONCLUSIVE |
| T_N construction | T_diag + Off_lin as 2x2 at rate-1/2 | UNDERSPECIFIED |
| Cross-freq closure | Closure space for Off_lin | CLOSES_ON_ENLARGED_SPAN |
| T_V spectrum | T_V on V_M at rate-1/2 | RECURSION_UNDERSPECIFIED |
| **T_lead corrected (eighth)** | T_lead at corrected rate | **H_T_LEAD_CARRIES_DIFFERENT_RATE (43/45 over Q)** |

The collective seven-probe bequest at rate-1/2 reads as: **"no Q-constructable finite-rank operator carries a discrete eigenvalue at rate 1/2 → Nisoli framework structurally inapplicable at the rate-1/2 target."**

The eighth probe (T_lead at the corrected rate) flipped this with a POSITIVE finding: **the corrected rate IS carried by a clean Q-rational eigenvalue 43/45.** So the rate-1/2 "wrong target" reading was correct, and the structural-inapplicability bequest was specific to the rate-1/2 target.

---

## 2. What the corrected-rate Nisoli closure test tells us

The Route A Nisoli closure inequality test at λ = 43/45 (this work) lands:

- **Phase 1-2-3**: closure inequality fires at (r=3, A=3, K=6) under `C_A = 1` parameterisation. Smallest firing cell. M_3'' = 24.43 exact.
- **Phase 4**: A is not extractable from any project source. The (r=3, A=3, K=6) firing is conditional on a `C_A = 1` Tao Prop 1.17 that **is not what Tao delivers**.

So the eighth probe split the structural-inapplicability question into TWO parts:

### Part 1: discrete eigenvalue exists at the rate of interest

- Rate-1/2 (seven prior probes): NO. R77.3 falsified; R_k intractable; K_k has no eigenvalue near 1/2.
- Corrected rate 43/45 (eighth probe): **YES**. T_lead has clean Q-rational eigenvalue 43/45.

### Part 2: Tao Prop 1.17's A is extractable to make Nisoli's closure inequality fire

- Rate-1/2: NO (BOOKKEEPING_PHASE1 INFEASIBLE).
- Corrected rate 43/45: NO (**SAME OBSTRUCTION** — BOOKKEEPING_PHASE1 INFEASIBLE generalises to any rate).

---

## 3. The structural bequest's re-reading

The seven-probe framing "Nisoli structurally inapplicable at rate-1/2" combined TWO obstructions that should be **separated**:

- **Obstruction 1 (operator-theoretic):** The R77 framework conjectured T_3 (3x3 companion) but R77.3 falsified the 3-mode model that T_3 was built on. So there was no operator at rate-1/2 to apply Nisoli to. ← **rate-specific to 1/2.**

- **Obstruction 2 (analytic / Tao-A):** Even with a hypothetical operator, the Nisoli closure inequality requires Tao Prop 1.17's effective `C_A` to satisfy `η < 1`. BOOKKEEPING_PHASE1 shows this is INFEASIBLE structurally — the obstruction is in Tao's renewal-process METHOD, not in the bookkeeping quality. ← **rate-invariant.**

The seven-probe trajectory conflated these. The eighth probe (this work) **separates them cleanly**:

- Obstruction 1 lifts at the corrected rate: T_lead exists with clean eigenvalue 43/45 over Q.
- **Obstruction 2 PERSISTS at the corrected rate.** This is the load-bearing finding.

---

## 4. Reframed trajectory bequest

Pre-eighth-probe: "Nisoli framework structurally inapplicable at rate-1/2."

Post-eighth-probe (post-this-work): **"Nisoli framework's analytic prerequisite — polynomial-in-A bound on `|μ̂_n(ξ)|` from Tao Prop 1.17 — is structurally infeasible via Tao's method. This obstruction is rate-invariant; it persists at any candidate rate including the corrected rate 43/45 where the operator-theoretic anchor (T_lead) IS clean."**

So:

- The seven-probe "rate-1/2 wrong target" framing is **correct** (rate-1/2 is wrong target; corrected rate is 43/45 over Q).
- The seven-probe "Nisoli inapplicable" framing is **partly correct, partly conflated**:
  - **Operator-theoretic anchor** at the corrected rate: AVAILABLE (T_lead).
  - **Tao-A analytic input**: NOT AVAILABLE (BOOKKEEPING_PHASE1's rate-invariant infeasibility).

The actual blocker is the Tao-A side, which generalises across all rates. The operator-theoretic side opens up at the corrected rate.

---

## 5. Routing implications

### For the c = 7/45 closure roadmap

Three independent obstructions (per M3_DISPOSITION.md and this Phase 6 reframe):

1. **Tao Prop 1.17 effective `C_A`**: INFEASIBLE via Tao's method (BOOKKEEPING_PHASE1). Persists at corrected rate.
2. **|K| bilinear bound**: DELIVERED at family level (PATH2 strict r ≤ 3, HENSEL polylog-free r ≥ 4).
3. **M_3 = ‖R(z, T)‖ for a characterised T**: At rate-1/2: T uncharacterised (R77.3 falsified). **At corrected rate: T_lead characterised, M_3'' = 24.43 exact.** ← **CHANGED by this probe.**

So Item 3 is solved at the corrected rate. Items 1 and 2 are unchanged. The block is now **only** at Item 1 (Tao-A).

### For the Tao communication / paper bequest

The (publishable) structural findings, in order of confidence:

1. **Bilinear bound at family level**: PATH2 strict 2√N at r ≤ 3 + HENSEL polylog-free 2√p·√N at r ≥ 4. Already in publication scope per `project_collatz_r78_bilinear_cracked` memory note.

2. **T_lead's exact eigenvalue 43/45 over Q with structural derivation**: 1 − Σ_g W_+(g) = 1 − 2/45. Closed-form within-level cross-frequency derivation. **First positive spectral anchor for c = 7/45 from the eight-probe trajectory.**

3. **Nisoli closure inequality M_3'' = 24.43, smallest firing cell (r=3, A=3, K=6) under `C_A = 1`**: rigorous up to the Tao-A blocker. **Conditional closure** with explicit constants — useful as a target for a future Fourier-analytic bound (outside Tao's method) that would deliver polynomial-in-A constants.

4. **Tao Prop 1.17 effective `C_A` is structurally INFEASIBLE via Tao's renewal-process method**: BOOKKEEPING_PHASE1's negative result is publishable as a no-go (it shows where Tao's specific method fails to deliver Nisoli's analytic prerequisite).

---

## 6. The honest seven-probe trajectory + eighth probe verdict

> **The seven-probe trajectory established that no Q-constructable rate-1/2 operator captures the actual c = 7/45 dynamics — the rate-1/2 target was wrong. The eighth probe (T_lead at corrected rate) establishes that the corrected rate 43/45 IS carried by a clean Q-rational eigenvalue, and that the resolvent norm M_3'' = 24.43 is tight enough that Nisoli's closure inequality fires at the smallest plausible cell (r=3, A=3, K=6) — under the optimistic `C_A = 1` parameterisation. The remaining rigorous blocker is Tao Prop 1.17's effective `C_A`, which BOOKKEEPING_PHASE1 establishes is INFEASIBLE via Tao's renewal-process method. This blocker is rate-invariant: it persists at the corrected rate just as it did at rate-1/2.**

The structural finding is: **the rate-1/2 obstruction had a rate-specific component (no eigenvalue) and a rate-invariant component (no Tao-A). The eighth probe shows the rate-specific component LIFTS at the corrected rate, but the rate-invariant component PERSISTS. Closure remains blocked, but now by a single identifiable obstruction (Tao-A) rather than the compound rate-1/2 obstruction.**

---

## 7. Files

- T_LEAD_CORRECTED_DISPOSITION.md (eighth probe disposition)
- BOOKKEEPING_PHASE1_DISPOSITION.md (Tao-A INFEASIBLE)
- NISOLI_CLOSURE_CORRECTED_PHASE1.md (closure inequality articulated)
- NISOLI_CLOSURE_CORRECTED_TABLE.md (Phase 2 tabulation)
- NISOLI_CLOSURE_CORRECTED_FIRING.md (Phase 3 firing cells)
- NISOLI_CLOSURE_CORRECTED_HONEST_READ.md (Phase 4 A-extractability)
- NISOLI_CLOSURE_CORRECTED_STRUCTURAL.md (this file)
