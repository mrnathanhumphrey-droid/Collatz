# CANDIDATE_A_LIFT_CONSTRUCTION — construction of lift_n(R_k) for all (n, k)

**Date:** 2026-05-12. Wilson (analyst) reporting to Nathan. Phase 2 of the Reading A scoping probe.

---

## Lift operator T_{k→k+1}

R77.5 Stage 1 anchor: the natural lift T : V_k → V_{k+1} acts via uniform-1/3 split across coprime preimages. For u ∈ V_k and r' coprime in Z/3^{k+1}:

  T(u)(r') := u(r' mod 3^k) / 3.

This is well-defined because every coprime r' in Z/3^{k+1} has r' mod 3^k coprime in Z/3^k (the 3 lifts of each coprime r ∈ Z/3^k are r, r + 3^k, r + 2 · 3^k, all coprime in Z/3^{k+1}).

`result_77_5_compute_R_k.py::lift_pi(pi_k, k)` implements this for π_k. Inner-product property:

  ⟨T(u), T(v)⟩_{V_{k+1}} = ⟨u, v⟩_{V_k} / 3   (isometric up to √3).

---

## Iterated lift T^{n−k−1}

For 0 ≤ k < n, define lift_n : V_{k+1} → V_n as T^{n−k−1} (iterated lift). Algebraically, for w ∈ V_{k+1}:

  **lift_n(w)(r') := w(r' mod 3^{k+1}) / 3^{n−k−1}**

for r' coprime in Z/3^n.

`candidate_a_compute.py::lift_n_of_R_k(R_k, k, n)` implements this.

---

## Residuals R_k

R_k for k = 1, ..., 5 are taken from R77.5 anchor: R_k := π_{k+1} − T(π_k) ∈ W_k ⊂ V_{k+1}. All exact rational.

For k = 0, define R_0 := π_1 − π_∞^{(1)} where π_∞^{(1)} is the uniform distribution 1/2 on V_1 = {1, 2}. R_0 ∈ V_1; sum is 0 (mean-zero on V_1 = W_0 ∪ constants split). At level 1, π_1(1) = 2/3, π_1(2) = 1/3 from the project's R64 anchor, so R_0(1) = 1/6, R_0(2) = −1/6. ‖R_0‖² = 2 · (1/6)² = 1/18 ≈ 5.556e-2.

---

## Verification (a): lift orthogonality across distinct k

Pre-registered adversarial check A2: ⟨lift_n(R_{k1}), lift_n(R_{k2})⟩ = 0 for k1 ≠ k2.

`candidate_a_compute.py` verifies this over Q for all valid pairs at n = 3, 4, 5, 6 — total 34 (k1, k2) pairs with k1 < k2 ≤ n−1:

> **Result: ⟨lift_n(R_{k1}), lift_n(R_{k2})⟩ = 0/1 exactly for all 34 checked pairs.** ✓

CSV record: `candidate_a_lift_orthogonality.csv`.

This is automatic from the structural identity R_k ∈ W_k ⊥ T(V_k) (R77.5 §3.1): T-iterates of R_{k1} for k1 < k2 land in T-iterates of T(V_{k1}) ⊂ T(V_{k2}), which is orthogonal to W_{k2} containing R_{k2}'s lift. The verification confirms the lift apparatus preserves the orthogonality.

---

## Verification (b): lift norm scaling

R77.5 Stage 1: ‖T(u)‖² = ‖u‖² / 3. Iterating gives ‖lift_n(R_k)‖² = ‖R_k‖² / 3^{n−k−1}.

`candidate_a_compute.py` verifies this over Q at all valid (n, k). All passed exactly. Implementation: each ‖lift_n(R_k)‖² is computed directly and compared to ‖R_k‖² · 1/3^{n−k−1}; equality holds in Q.

---

## Verified residual norms (R77.5 Stage 1 reproduction)

| k | dim V_{k+1} | ‖R_k‖² (exact) | float |
|---|-------------|----------------|-------|
| 0 | 2           | 1/18           | 5.556e-2 |
| 1 | 6           | 10/189         | 5.291e-2 |
| 2 | 18          | 31370/1835001  | 1.710e-2 |
| 3 | 54          | (R77.5 exact)  | 5.731e-3 |
| 4 | 162         | (R77.5 exact)  | 1.916e-3 |
| 5 | 486         | (R77.5 exact)  | 6.395e-4 |

Ratios ‖R_k‖² / ‖R_{k−1}‖² ≈ 1/3 for k ≥ 2 (R77.5 finding). Sum check Σ_{r'} R_k(r') = 0 for all k ≥ 1 verified exactly.

---

## Files

- `candidate_a_compute.py` — full computation
- `candidate_a_lift_orthogonality.csv` — exact-rational lift orthogonality table (34 rows, all zero)
