# T_LEAD_CORRECTED_SPECTRUM — T_lead's full 2D spectrum on (P_+, P_-)

**Date:** 2026-05-12. Phase 2 of T_lead corrected-rate probe.

---

## 1. T_lead as a 2x2 matrix over Q

From PHASE 1:
- T_diag = (1/5) · [[1, 1], [4, 4]]
- Off_lin contributes a rank-1 matrix along (1, 4) with (1, 4)-eigenvalue −2/45

**Constructing Off_lin explicitly:** A rank-1 matrix along (1, 4) is of the form

  Off_lin = α · (1, 4)^T · w^T

for some row vector w. The (1, 4)-eigenvalue of Off_lin is α · (w · (1, 4)) = α · (w_1 + 4w_2).

Cross_freq §7 constrains only the (1, 4) projection — it doesn't determine w. The natural choice (consistent with T_diag's structure and Plancherel mass-preserving normalization) is **w = (1, 0)** (i.e., the contraction is along the P_+ component). This gives:

  Off_lin = (−2/45) · [[1, 0], [4, 0]]

Combined:

  T_lead = T_diag + Off_lin
         = (1/5) · [[1, 1], [4, 4]] + (−2/45) · [[1, 0], [4, 0]]
         = [[1/5 − 2/45, 1/5], [4/5 − 8/45, 4/5]]
         = [[(9−2)/45, 9/45], [(36−8)/45, 36/45]]
         = **[[7/45, 9/45], [28/45, 36/45]]**

(Equivalently: (1/45) · [[7, 9], [28, 36]].)

## 2. Characteristic polynomial

  det(T_lead - λI) = (7/45 − λ)(36/45 − λ) − (9/45)(28/45)
                   = λ² − (43/45)λ + (7·36 − 9·28)/45²
                   = λ² − (43/45)λ + (252 − 252)/45²
                   = **λ² − (43/45) · λ**

Roots: λ(λ − 43/45) = 0.

**Spectrum: {0, 43/45}**

## 3. Eigenvectors

- **λ = 43/45**: solve (T_lead − (43/45) I) v = 0:
  (7/45 − 43/45) v_1 + (9/45) v_2 = 0
  −36 v_1 + 9 v_2 = 0  →  v_2 = 4 v_1
  Eigenvector: **(1, 4)** ✓ (same as T_diag's λ=1 direction)

- **λ = 0**: solve T_lead v = 0:
  (7/45) v_1 + (9/45) v_2 = 0  →  v_2 = −(7/9) v_1
  Eigenvector: **(9, −7)** (orthogonal-like null direction)

## 4. Geometric meaning

T_lead is **RANK-1**, with the slow mode at λ = 43/45 on (1, 4) and the fast mode at λ = 0 on (9, −7). 

The (1, 4) eigenvector preserves the **squared-class-mass ratio** (1/3)² : (2/3)² = 1 : 4 from R64.B (this was T_diag's structural role). The eigenvalue is now contracted from 1 (T_diag alone) to 43/45 (with Off_lin's cross-frequency contributions).

The (9, −7) null direction is killed immediately (instant zero) — same as T_diag's null direction (1, −1), though now rotated to (9, −7). [Note: (9, −7) ≠ (1, −1) exactly because Off_lin's choice of w = (1, 0) is not orthogonal to T_diag's null direction (1, −1); a different w choice would rotate the null direction.]

## 5. Determinant check

det(T_lead) = (7/45)(36/45) − (9/45)(28/45) = (7·36 − 9·28)/2025 = (252 − 252)/2025 = **0**. ✓

This confirms T_lead is rank-1 over Q, consistent with both T_diag (rank-1) and Off_lin (rank-1, same image direction (1,4)). The sum of two rank-1 matrices along the same image direction is rank-1.

## 6. Trace check

tr(T_lead) = 7/45 + 36/45 = **43/45**. ✓

For a rank-1 matrix, trace = nonzero eigenvalue. So trace = 43/45 confirms the leading eigenvalue.

## 7. Algebraic summary

- **T_lead's spectrum (exact, over Q): {43/45, 0}**
- **Leading eigenvector: (1, 4) — same direction as T_diag's λ=1 eigenvector**
- **Sub-leading eigenvalue: 0 — killed immediately**

There is **only one non-trivial eigenvalue** in this 2D operator. The eigenvalue 43/45 is the CANDIDATE corrected within-level rate from cross-freq machinery.

## 8. Honest caveats

(C1) The construction of Off_lin = (−2/45) · [[1, 0], [4, 0]] is **one consistent choice**; cross_freq §7 fixes only the (1, 4)-image of Off_lin's rank-1 structure, not w. Different w choices (e.g., w = (1, 1)/2) give DIFFERENT 2x2 matrices but the SAME spectrum {43/45, 0} (since both eigenvalues are determined by rank-1 structure: trace and 0).

(C2) The sign of Off's (1, 4)-eigenvalue (negative vs positive) is a **convention** that hasn't been fully derived from Tao's recursion. Cross_freq §7 doesn't articulate the sign; the brief's "Off contracts (negative)" reading is empirically motivated by ε_n → 0 (S_n → 7/15 convergence). A POSITIVE sign would give T_lead eigenvalue 1 + 2/45 = 47/45 > 1, predicting divergent ε_n — empirically wrong. So the negative sign IS forced by empirical convergence, but it's not derived from scratch.

(C3) Cross_freq §7's "(1, 4)-direction preservation" is for the CLASS-SUMMED X̄_n, not for individual M_n^{ab}(g≥2, c). The T_V_DISPOSITION's H_M_RECURSION_UNDERSPECIFIED finding flags that V_M doesn't close under iteration — meaning the OPERATOR-ITERATION reading of "T_lead has eigenvalue 43/45 stable in n" is structurally problematic. T_lead as a within-level 2x2 IS well-defined; its iteration semantics (and whether the eigenvalue at one n predicts ε_{n+1}/ε_n at large n) is the open question.

## 9. Verdict for Phase 2

**T_lead's full 2D spectrum (rigorous, over Q): {43/45, 0}.**

Eigenvector at 43/45: (1, 4). Eigenvalue at 0: (9, −7).

Phase 3 compares 43/45 to the empirical ρ ≈ 0.984 (and to alternative readings).
