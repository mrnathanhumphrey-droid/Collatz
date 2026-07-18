# CANDIDATE_A_C_NK — table of c_{n,k} := ⟨φ_n, lift_n(R_k)⟩ over Q

**Date:** 2026-05-12. Wilson (analyst) reporting to Nathan. Phase 3 of the Reading A scoping probe.

---

## (a) Full table — exact rationals and float approximation

All computations via `fractions.Fraction` over Q. CSV: `candidate_a_c_nk_table.csv`.

| n | k | c_{n,k} (exact)                  | float            | sign |
|---|---|----------------------------------|------------------|------|
| 1 | 0 | 1/6                              | +1.667e-1        | +    |
| 2 | 0 | 0/1                              | 0                | 0    |
| 2 | 1 | 10/21                            | +4.762e-1        | +    |
| 3 | 0 | 0/1                              | 0                | 0    |
| 3 | 1 | 0/1                              | 0                | 0    |
| 3 | 2 | 31370/67963                      | +4.616e-1        | +    |
| 4 | 0 | 0/1                              | 0                | 0    |
| 4 | 1 | 0/1                              | 0                | 0    |
| 4 | 2 | 0/1                              | 0                | 0    |
| 4 | 3 | 143195649659456490 / 308468774477179141 | +4.642e-1 | + |
| 5 | 0 | 0/1                              | 0                | 0    |
| 5 | 1 | 0/1                              | 0                | 0    |
| 5 | 2 | 0/1                              | 0                | 0    |
| 5 | 3 | 0/1                              | 0                | 0    |
| 5 | 4 | (60-digit num)/(61-digit den)    | +4.655e-1        | +    |
| 6 | 0 | 0/1                              | 0                | 0    |
| 6 | 1 | 0/1                              | 0                | 0    |
| 6 | 2 | 0/1                              | 0                | 0    |
| 6 | 3 | 0/1                              | 0                | 0    |
| 6 | 4 | 0/1                              | 0                | 0    |
| 6 | 5 | (217-digit num)/(218-digit den)  | +4.662e-1        | +    |

**The most important load-bearing fact:** for every n ≥ 2 and every k with 0 ≤ k < n − 1, **c_{n,k} = 0/1 exactly over Q.** Only c_{n, n−1} is nonzero.

For n = 1, c_{1, 0} = 1/6 (the only k = 0 case; the n = 1 transient).

For n ≥ 2, c_{n, n−1} = S_n (since ⟨φ_n, π_∞^{(n)}⟩ = 0/1 — see decomposition section).

---

## (b) Dominant-k k*(n)

| n | k*(n) | |c_{n, k*(n)}| | share of Σ_k |c_{n,k}| |
|---|-------|----------------|------------------------|
| 1 | 0     | 1/6            | 1.0000 (only k = 0)    |
| 2 | 1     | 10/21          | **1.0000** (all other k zero) |
| 3 | 2     | 31370/67963    | **1.0000**             |
| 4 | 3     | (S_4)          | **1.0000**             |
| 5 | 4     | (S_5)          | **1.0000**             |
| 6 | 5     | (S_6)          | **1.0000**             |

k*(n) = n − 1 for all n in the test range. The dominant carries 100% of the total |c|.

---

## (c) Signed sum S_n^signed := Σ_k c_{n,k}

| n | Σ_k c_{n,k}              |
|---|--------------------------|
| 1 | 1/6 ≈ +1.667e-1          |
| 2 | 10/21 ≈ +4.762e-1        |
| 3 | 31370/67963 ≈ +4.616e-1  |
| 4 | S_4 ≈ +4.642e-1          |
| 5 | S_5 ≈ +4.655e-1          |
| 6 | S_6 ≈ +4.662e-1          |

Note: Σ_k c_{n,k} = S_n exactly (because ⟨φ_n, π_∞^{(n)}⟩ = 0 for n ≥ 2). So the signed sum converges to S_∞ = 7/15 ≈ 0.4667, NOT to zero.

### Decomposition sanity check (verified over Q)

The decomposition Σ_k c_{n,k} = ⟨φ_n, π_n − π_∞^{(n)}⟩ requires comparing against S_n − ⟨φ_n, π_∞^{(n)}⟩:

| n | ⟨φ_n, π_∞^{(n)}⟩ |
|---|--------------------|
| 1 | 1/2                |
| 2 | 0/1                |
| 3 | 0/1                |
| 4 | 0/1                |
| 5 | 0/1                |
| 6 | 0/1                |

For n ≥ 2, ⟨φ_n, π_∞^{(n)}⟩ = 0 exactly. This is because φ_n ∈ W_{n−1} (zero 3-fiber-mean at scale 3^{n−1}, see CANDIDATE_A_PHI_CONSTRUCTION.md) and π_∞^{(n)} is the uniform / constant function on V_n, which lies in T^{n−1}(V_1) (the "scale-0" lift), orthogonal to W_{n−1}.

For n = 1: φ_1 ∈ V_1; ⟨φ_1, π_∞^{(1)}⟩ = 1/2 (the n = 1 setting is different — V_1 has only 2 states and the W_0 / constants split is trivial).

The decomposition Σ_k c_{n,k} = ⟨φ_n, π_n − π_∞^{(n)}⟩ holds exactly at every n in the test:
- n = 1: 1/6 = 2/3 − 1/2 = 1/6 ✓
- n = 2: 10/21 = 10/21 − 0 = 10/21 ✓
- (n = 3..6 all verified ✓)

So **the W_k filtration is structurally complete and exact** as a decomposition apparatus. The pre-registered (C1) sanity check **passes**. The failure to localize rate-1/2 is **NOT** a structural bug in the filtration — it's a property of φ_n.

---

## (d) A1 ratios |c_{n, k*(n)}| / |c_{n−1, k*(n−1)}|

| n | ratio    | target |
|---|----------|--------|
| 2 | 2.857143 | 0.5    |
| 3 | 0.969307 | 0.5    |
| 4 | 1.005719 | 0.5    |
| 5 | 1.002802 | 0.5    |
| 6 | 1.001405 | 0.5    |

The ratios **converge to 1.0**, not to 0.5. Because k*(n) = n − 1 always carries the entire ε_n + 7/15 = S_n, and S_n → 7/15 ≈ 0.4667, the ratio c_{n, n−1} / c_{n−1, n−2} → 7/15 / (7/15) = 1.

**A1 (rate-0.5 decay of dominant-k) is decisively rejected.** The dominant doesn't decay; it converges to 7/15.

---

## (e) A2 ratios |S_n^signed| / |S_{n−1}^signed|

| n | ratio    | target |
|---|----------|--------|
| 2 | 2.857143 | 0.5    |
| 3 | 0.969307 | 0.5    |
| 4 | 1.005719 | 0.5    |
| 5 | 1.002802 | 0.5    |
| 6 | 1.001405 | 0.5    |

Identical to (d) because at each n only k = n − 1 contributes. **A2 (signed sum rate-0.5) decisively rejected.**

---

## (f) Where rate-1/2 actually lives (diagnostic supplement)

Rate-1/2 lives in **the deviation of c_{n, n−1} from its asymptote 7/15**, which is by construction ε_n = S_n − 7/15. Empirical ratios ε_n / ε_{n−1}:

| n | ε_n / ε_{n−1} |
|---|---------------|
| 2 | +0.048 (n=1 transient) |
| 3 | −0.535 |
| 4 | +0.482 |
| 5 | +0.470 |
| 6 | +0.432 |

These match R77.4 / R77.6 anchors. The rate-1/2 phenomenon is REAL — it just doesn't decompose into the W_k filtration in any nontrivial way. It lives in the n-evolution of a single coefficient c_{n, n−1} → 7/15, where each c_{n, n−1} is the full bilinear pair-form moment at level n.

---

## (g) Sign pattern of c_{n,k}

| n | signs (k = 0, 1, …, n−1) |
|---|---------------------------|
| 2 | 0 +                       |
| 3 | 0 0 +                     |
| 4 | 0 0 0 +                   |
| 5 | 0 0 0 0 +                 |
| 6 | 0 0 0 0 0 +               |

All nonzero c are positive (single sign). No cross-k cancellation possible (only one term per row).

---

## Files

- `candidate_a_compute.py` — full computation
- `candidate_a_c_nk_table.csv` — full Q + float table
- `candidate_a_phi_n_verify.csv` — S_n / ε_n / decomposition verification
- `candidate_a_lift_orthogonality.csv` — lift orthogonality verification
