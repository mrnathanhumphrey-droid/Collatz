# T_LEAD_CORRECTED_CLOSURE — Nisoli resolvent norm at T_lead's corrected eigenvalue

**Date:** 2026-05-12. Phase 5 of T_lead corrected-rate probe. Conditional on Phase 1-2 producing a clean discrete eigenvalue.

---

## 1. T_lead eigenvalue (from Phase 1-2)

  λ_T_lead = 43/45 ≈ 0.9556 on (1, 4) eigenvector

This IS a clean discrete eigenvalue over Q. So Phase 5 (Nisoli resolvent computation) is well-defined.

## 2. Resolvent norm M_3''

For Nisoli Theorem 2.15-type closure, the relevant resolvent norm is:

  M_3'' = ||(I − T_lead)^{−1}||

For a 2x2 rank-1 operator T_lead with spectrum {0, 43/45}, the inverse of (I − T_lead) is:

  (I − T_lead)^{−1} = I + T_lead/(1 − 43/45) + (kill the 0 eigenmode)
  
More carefully: for a diagonalizable operator with spectrum {λ_1, λ_2}, eigenvectors {v_1, v_2}, the resolvent at z = 1 is:

  (I − T_lead)^{−1} v_i = (1 − λ_i)^{−1} v_i

So (I − T_lead)^{−1} has spectrum {1/(1−0), 1/(1 − 43/45)} = **{1, 45/2}** = {1, 22.5}.

Operator norm ||(I − T_lead)^{−1}|| (in any basis-equivalent ℓ² norm) is bounded by the **larger eigenvalue magnitude × condition-number-of-eigenbasis factor**.

Since T_lead's eigenbasis is (1, 4) and (9, −7) (Phase 2), the eigenbasis is NOT orthogonal. The condition number κ of the eigenbasis matters for the operator-norm-vs-spectral-radius gap.

Eigenbasis matrix V = [[1, 9], [4, −7]]. det(V) = −7 − 36 = −43. 
||V|| ≈ sqrt(1+9²) + sqrt(4²+7²) ≈ 9.06 + 8.06 ≈ 17 (upper bound via Frobenius).
||V^{−1}|| ≈ (1/|det V|) · ||adj V|| ≈ (1/43) · 17 ≈ 0.4.

So condition number κ(V) ≈ 17 × 0.4 ≈ 6.8.

**Resolvent operator norm bound:** ||(I − T_lead)^{−1}|| ≤ κ(V) · max(1, 45/2) = 6.8 · 22.5 ≈ **153**.

A sharper bound via direct computation:

  (I − T_lead) = I − [[7/45, 9/45], [28/45, 36/45]] = [[38/45, −9/45], [−28/45, 9/45]] = (1/45)·[[38, −9], [−28, 9]]

det(I − T_lead) = (1/45²)·(38·9 − (−9)(−28)) = (1/2025)·(342 − 252) = (1/2025)·90 = 2/45.

(I − T_lead)^{−1} = (1/det)·[[9/45, 9/45], [28/45, 38/45]] = (45/2)·(1/45)·[[9, 9], [28, 38]] = (1/2)·[[9, 9], [28, 38]]

= [[4.5, 4.5], [14, 19]]

Operator ℓ² norm of (I − T_lead)^{−1}: this is a 2x2 matrix; compute via SVD.

Trace((I−T)^{−T} · (I−T)^{−1}) = 4.5²+14²+4.5²+19² = 20.25 + 196 + 20.25 + 361 = 597.5

Eigenvalues of A^T A where A = (I−T_lead)^{−1}:
A^T A = [[4.5²+14², 4.5·4.5 + 14·19], [..., 4.5²+19²]] = [[216.25, 286.25], [286.25, 381.25]]

Trace = 597.5. det = 216.25·381.25 − 286.25² = 82445 − 81939 = 506 (approx; exact via Fractions below).

Exact: A = (1/2) · [[9, 9], [28, 38]]. So A^T A = (1/4) · [[9²+28², 9·9+28·38], [..., 9²+38²]] = (1/4) · [[81+784, 81+1064], [..., 81+1444]] = (1/4) · [[865, 1145], [1145, 1525]].

Eigenvalues: roots of λ² − (trace) λ + det = 0 where trace = (865+1525)/4 = 597.5; det = (865·1525 − 1145²)/16 = (1319125 − 1311025)/16 = 8100/16 = 506.25.

So σ²_max + σ²_min = 597.5, σ²_max · σ²_min = 506.25.
σ²_max = (597.5 + sqrt(597.5² − 4·506.25))/2 = (597.5 + sqrt(357006.25 − 2025))/2 = (597.5 + sqrt(354981.25))/2 = (597.5 + 595.8)/2 = 596.65

σ_max = sqrt(596.65) ≈ **24.43**

So **||(I − T_lead)^{−1}|| ≈ 24.4** (operator ℓ² norm).

## 3. M_3'' value (recapping)

**M_3'' = ||(I − T_lead)^{−1}|| ≈ 24.4** (slightly larger than the spectral radius 45/2 = 22.5 due to non-orthogonal eigenbasis).

## 4. Nisoli closure inequality

Nisoli's closure (R77 framework) requires:

  |K| · K^{−A} · M_3'' < 1

per the brief. Where:
- |K| ≤ 2√N at r ≤ 3 (delivered in the bilinear bound, project_collatz_r78_bilinear_cracked)
- N = p^{r−1} = 3^{r−1} at q = 3
- A is the Nisoli closure exponent (project-internal)

The (1, 4) eigenvalue 43/45 (close to 1) makes M_3'' ≈ 24 LARGE — Nisoli closure would require:

  |K| · K^{−A} · 24.4 < 1
  K^{−A} · |K| · 24.4 < 1
  K^{−A+1} < 1/24.4

For A > 1 and K = 3^{r−1} large: K^{−A+1} → 0 as r increases. So the inequality CAN potentially be satisfied at sufficiently large r — but requires Nisoli's A to be > 1 and the gap K^{A−1} > 24.4 to be open.

## 5. Comparison to rate-1/2 closure (R77.3 falsified case)

If R77.3's rate-1/2 had been correct (T_lead eigenvalue 1/2 on (1,4)), then:
  (I − T_lead)^{−1} eigenvalue at λ=1/2 is 1/(1−1/2) = 2
  M_3'' would be O(2-3), much smaller than 24.4

So the corrected rate 43/45 ≈ 0.956 makes the Nisoli closure SIGNIFICANTLY HARDER (M_3'' ~ 8-10× larger than at rate 1/2). The closure inequality |K|·K^{−A}·M_3'' < 1 requires substantially more growth in K (more r-iterations) to dominate.

## 6. Honest assessment of closure viability

The corrected-rate Nisoli closure at λ = 43/45 is:
- **Possible in principle** if the closure exponent A is > 1 and r-iteration eventually dominates.
- **Significantly harder than the rate-1/2 case** that R77.3 falsified.
- **Critical caveats**:
  - The eigenvalue 43/45 was derived UNDER THE ASSUMPTION that cross-freq's class-summed X̄_n(g≥2) projects fully onto (1,4) of (P_+, P_-). The T_V_DISPOSITION findings (H_M_RECURSION_UNDERSPECIFIED) suggest this projection is NOT a derivable algebraic identity — it's a structural conjecture that survives only at the X̄ level, not at the M_n^{ab}(g, c) level.
  - The empirical evidence through k=13 (T_LEAD_CORRECTED_EMPIRICAL) is in a transient that doesn't decisively confirm 43/45 as the asymptotic rate. The Hadamard radius is currently at 0.64 and inward-trending, with no decisive distinction between 1.046 (T_lead) and 1.016 (slow-mode 0.984) at n=13.
  - R77.6's branch-cut reading suggests the rate phenomenon is a CONTINUOUS spectrum endpoint, not a discrete eigenvalue. T_lead at λ=43/45 is FORMALLY a discrete eigenvalue, but it's the eigenvalue of a WITHIN-LEVEL 2x2 operator that may not faithfully represent the asymptotic operator.

## 7. Verdict for Phase 5

**M_3'' ≈ 24.4 at T_lead's corrected eigenvalue λ = 43/45.**

The Nisoli closure inequality |K|·K^{−A}·M_3'' < 1 is potentially satisfiable at sufficiently large r if Nisoli's A > 1, but requires significantly more r-growth than the (falsified) rate-1/2 case.

**This does NOT yet constitute a viable closure of c = 7/45**: the (1, 4)-projection assumption is structural-not-derived, and the empirical asymptote at n=13 is too transient to confirm 43/45 as the rigorous rate.

A clean path forward: (a) derive the (1, 4)-projection of X̄_n rigorously (closing the cross_freq §7 structural gap), then (b) compute Nisoli closure inequality with explicit A and K values.

## 8. Files

- T_LEAD_CORRECTED_PHASE1.md — Off_lin sum (Phase 1)
- T_LEAD_CORRECTED_SPECTRUM.md — T_lead 2D spectrum (Phase 2)
- T_LEAD_CORRECTED_EMPIRICAL.md — Phase 4 cross-check
- T_LEAD_CORRECTED_CLOSURE.md (this file) — Phase 5
- t_lead_corrected.py — verification script
