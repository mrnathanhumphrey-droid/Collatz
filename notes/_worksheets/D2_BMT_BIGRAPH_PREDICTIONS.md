# D2 — BMT and Bigraph predicted moment values for Syracuse rows (b), (d), (f)

**Date:** 2026-05-15
**Mode:** E — derivations from verbatim Defns; verbatim source see `D2_BMT_BIGRAPH_VERBATIM.md`.

---

## 0. Setup

Syracuse indices j ∈ N have natural linear order j_1 < j_2 < .... We test the BMT/bigraph moment formulas on a 2-vertex induced sub(di)graph with vertices {j_1, j_2}, taking various digraph choices.

We compare to Syracuse's measured moments at level n = 3, (Z/27)*, 18 states, V_TRUNC = 16, Reading B (marginal centering X̃_j = Off_j − E[Off_j | b_prior]), four scalar reductions {tr_pi, vac_pi, delta_1, sum_entries}.

Throughout, **ϕ(X̃_j) = 0** by Reading B marginal centering (this is what "centering" means).

---

## 1. BMT predictions

### 1.1 Row (b): m = 2, i = (j_1, j_2), j_1 ≠ j_2

**ker[i] = {{1}, {2}}** (positions 1, 2 carry distinct labels).

Defn 2.8 kerG[i]: two positions p, p′ go in the same block iff i_p = i_{p′}. Since i_1 = j_1 ≠ j_2 = i_2, they never can be in the same block regardless of digraph. **kerG[i] = {{1}, {2}}** for every digraph G.

**BMT prediction (every digraph):**
> ϕ(a_1 a_2) = ϕ(X̃_{j_1}) · ϕ(X̃_{j_2}) = 0 · 0 = **0**.

### 1.2 Row (d): m = 3, i = (j_1, j_2, j_1) with j_1 < j_2

**ker[i] = {{1, 3}, {2}}.**

kerG[i]: positions 1 ∼ 3 requires (i_2, i_1) = (j_2, j_1) ∈ E (only one ℓ = 2 between k=1, k′=3).

**Case 1: Monotone digraph (Prop 3.5(iii)).** E = {(i, j) : j < i}. Then (j_2, j_1) ∈ E ⇔ j_1 < j_2. **TRUE.** So **kerG[i] = {{1, 3}, {2}}.**

> BMT prediction: ϕ(a_1 a_2 a_3) = ϕ[(a_k)|_{1,3}] · ϕ[(a_k)|_2] = ϕ(X̃_{j_1} · X̃_{j_1}) · ϕ(X̃_{j_2}) = ϕ(X̃_{j_1}²) · 0 = **0**.

**Case 2: Tensor digraph (Prop 3.5(i)).** Complete digraph; kerG[i] = ker[i] = {{1, 3}, {2}}. **Same prediction = 0.**

**Case 3: Boolean digraph (empty).** kerG[i] = singletons = {{1}, {2}, {3}}. ϕ(a_1)·ϕ(a_2)·ϕ(a_3) = 0·0·0 = **0**.

**Case 4: Anti-monotone digraph (E = {(i, j) : i < j}).** (j_2, j_1) ∈ E ⇔ j_2 < j_1 = FALSE. So 1 ≁ 3, kerG[i] = {{1}, {2}, {3}}, ϕ = 0·0·0 = **0**.

**Across ALL four standard digraph choices for BMT, row (d) prediction = 0.**

### 1.3 Row (f): m = 4, i = (j_1, j_2, j_1, j_2), j_1 < j_2

**ker[i] = {{1, 3}, {2, 4}}.**

kerG[i] (Defn 2.8):
- **1 ∼ 3:** ℓ = 2 between k=1, k′=3 with i_ℓ = j_2 ≠ i_1 = j_1. Requires (i_2, i_1) = (j_2, j_1) ∈ E.
- **2 ∼ 4:** ℓ = 3 between k=2, k′=4 with i_ℓ = j_1 ≠ i_2 = j_2. Requires (i_3, i_2) = (j_1, j_2) ∈ E.

**Case 1: Monotone digraph.** (j_2, j_1) ∈ E (since j_1<j_2). (j_1, j_2) ∈ E requires j_2 < j_1, **FALSE**. So **1 ∼ 3 but 2 ≁ 4**:

> **kerG[i] = {{1, 3}, {2}, {4}}.**
>
> BMT prediction: ϕ(X̃_{j_1}²) · ϕ(X̃_{j_2}) · ϕ(X̃_{j_2}) = ϕ(X̃_{j_1}²) · 0 · 0 = **0**.

**This contradicts the auditor's claim** that BMT under monotone digraph gives kerG = {{1,3},{2,4}}. The auditor's example was wrong: under "monotone digraph" with E = {(i,j) : j < i}, only ONE direction (j_2 → j_1) is an edge, not both directions. Position 4 stays a singleton because the "earlier" position 2 (with label j_2) requires the intervening position 3 (with label j_1) to have (j_1, j_2) ∈ E — which is the WRONG direction for monotone.

**Case 2: Tensor digraph (complete).** Both edges present. **kerG[i] = ker[i] = {{1, 3}, {2, 4}}.**

> BMT prediction: ϕ(X̃_{j_1} · X̃_{j_1}) · ϕ(X̃_{j_2} · X̃_{j_2}) = **ϕ(X̃_{j_1}²) · ϕ(X̃_{j_2}²)**.

Qualitatively non-zero. Quantitative magnitude requires numerical compute (see `d2_extract_phi_X_squared.py`).

**Case 3: Anti-monotone digraph (E = {(i,j): i<j}).** (j_2, j_1) ∈ E requires j_2<j_1 FALSE; (j_1, j_2) ∈ E requires j_1<j_2 TRUE. So **1 ≁ 3 but 2 ∼ 4**:

> kerG[i] = {{1}, {2, 4}, {3}}. BMT prediction: ϕ(X̃_{j_1})·ϕ(X̃_{j_2}²)·ϕ(X̃_{j_1}) = 0·ϕ(X̃_{j_2}²)·0 = **0**.

**Case 4: Boolean (empty).** All singletons → **0**.

### 1.4 Summary of BMT predictions

| Digraph | Row (b) | Row (d) | Row (f) |
|---|---|---|---|
| Monotone (Prop 3.5(iii)) | 0 | 0 | 0 |
| Anti-monotone | 0 | 0 | 0 |
| Tensor (complete, Prop 3.5(i)) | 0 | 0 | **ϕ(X̃_1²) · ϕ(X̃_2²)** ≠ 0 |
| Boolean (empty, Prop 3.5(ii)) | 0 | 0 | 0 |

**No BMT digraph predicts row (d) ≠ 0.** The Syracuse row (d) value 0.108 is incompatible with BMT under ANY digraph choice on a 2-vertex graph.

---

## 2. Bigraph predictions (Gilliers–Jekel Defn 1.4)

Defn 1.4: φ(a_1 ··· a_k) = ∑_{π ∈ P(c, G)} K^free_π(a_1, ..., a_k).

K^free_π denotes the **partitioned free cumulant** (NS06). For a partition π = {B_1, ..., B_r}, K^free_π is the product over blocks of free cumulants applied to each block's variables.

For centered variables (κ_1^free = ϕ(X) = 0):
- κ_1^free(X̃_j) = 0.
- κ_2^free(X̃_j, X̃_j) = ϕ(X̃_j²) − ϕ(X̃_j)² = ϕ(X̃_j²).
- κ_2^free(X̃_{j_1}, X̃_{j_2}) for j_1 ≠ j_2 = ϕ(X̃_{j_1} X̃_{j_2}) − ϕ(X̃_{j_1}) ϕ(X̃_{j_2}) = ϕ(X̃_{j_1} X̃_{j_2}) (Syracuse row (b), measured ≈ 0).

### 2.1 Row (b): m = 2, c = (j_1, j_2), j_1 ≠ j_2

P(c, G) candidates: π ∈ P(2) = {{{1},{2}}, {{1,2}}}.
- π = {{1,2}} requires c(1) = c(2) per Defn 1.2 (1): FALSE.
- π = {{1},{2}}: trivially in P(c, G).

> Bigraph prediction: φ(a_1 a_2) = K^free_{ {{1},{2}} }(a_1, a_2) = κ_1^free(a_1) · κ_1^free(a_2) = ϕ(X̃_{j_1}) · ϕ(X̃_{j_2}) = **0**.

Matches Syracuse row (b) ≈ 0 regardless of bigraph choice.

### 2.2 Row (d): m = 3, c = (j_1, j_2, j_1)

P(c, G) candidates:
- π = {{1},{2},{3}}: always admissible.
- π = {{1,3},{2}}: requires (per Defn 1.2 (1)) c(1) = c(3): TRUE. Per (2): i_1 = 1 < j = 2 < i_2 = 3, c(i_1) = j_1, c(j) = j_2; needs (j_1, j_2) ∈ E_1.
- π = {{1,2,3}}: requires c(1) = c(2) = c(3), but c(2) = j_2 ≠ j_1: FALSE.
- π = {{1,2},{3}}: requires c(1) = c(2): FALSE.
- π = {{1},{2,3}}: requires c(2) = c(3): FALSE.

**Case A: (j_1, j_2) ∈ E_1** (any bigraph admitting this edge, e.g. E_1 includes (j_1, j_2)). Then P(c, G) = {{{1},{2},{3}}, {{1,3},{2}}}.

> Bigraph: φ = K^free_{singletons}(...) + K^free_{ {{1,3},{2}} }(...)
>         = κ_1(X̃_{j_1})·κ_1(X̃_{j_2})·κ_1(X̃_{j_1}) + κ_2^free(X̃_{j_1}, X̃_{j_1})·κ_1(X̃_{j_2})
>         = 0 + ϕ(X̃_{j_1}²) · 0 = **0**.

**Case B: (j_1, j_2) ∉ E_1.** P(c, G) = {{{1},{2},{3}}} only.

> φ = 0.

**Across all bigraph choices, row (d) prediction = 0.** **MISMATCH with measured 0.108.**

### 2.3 Row (f): m = 4, c = (j_1, j_2, j_1, j_2)

Candidate partitions in P(4) respecting (1) (same-block ⇒ same color):
- π_0 = {{1},{2},{3},{4}} (singletons): always admissible.
- π_a = {{1,3},{2},{4}}: needs (j_1, j_2) ∈ E_1 (for i_1=1<j=2<i_2=3 with j_1∼j_1).
- π_b = {{1},{2,4},{3}}: needs (j_2, j_1) ∈ E_1 (for i_1=2<j=3<i_2=4 with j_2∼j_2).
- π_c = {{1,3},{2,4}} (crossing): needs BOTH conditions (2) at (i_1=1,j=2,i_2=3) → (j_1,j_2) ∈ E_1; AND (i_1=2,j=3,i_2=4) → (j_2,j_1) ∈ E_1; AND condition (3) at (i_1=1,j_1=2,i_2=3,j_2=4) with i_1 ∼_π i_2 and j_1 ∼_π j_2 and i_1 ≁_π j_1 → (c(i_2), c(j_1)) = (j_1, j_2) ∈ E_2.
- π_d = {{1,2,3,4}}: requires c constant — FALSE.
- π_e = {{1,2,3},{4}}, {{1,2,4},{3}}, etc.: require same-color within block — FALSE since j_1 ≠ j_2.

So at most 4 admissible partitions {π_0, π_a, π_b, π_c} depending on bigraph G.

**Case A: (j_1, j_2), (j_2, j_1) ∈ E_1 AND (j_1, j_2) ∈ E_2.** Then all 4 partitions admissible.

> φ = K^free_{π_0} + K^free_{π_a} + K^free_{π_b} + K^free_{π_c}
>   = 0 (singletons all centered) + κ_2^free(X̃_{j_1}, X̃_{j_1})·κ_1(X̃_{j_2})·κ_1(X̃_{j_2}) + κ_1(X̃_{j_1})·κ_2^free(X̃_{j_2}, X̃_{j_2})·κ_1(X̃_{j_1}) + κ_2^free(X̃_{j_1}, X̃_{j_1})·κ_2^free(X̃_{j_2}, X̃_{j_2})
>   = 0 + ϕ(X̃_{j_1}²)·0·0 + 0·ϕ(X̃_{j_2}²)·0 + ϕ(X̃_{j_1}²)·ϕ(X̃_{j_2}²)
>   = **ϕ(X̃_{j_1}²) · ϕ(X̃_{j_2}²)**.

Per Prop 1.5(5): (v,w) ∈ E_ten := E_1 ∩ E_1^op ∩ E_2 corresponds to **tensor independence**. So Case A = tensor between j_1 and j_2.

**Case B: (j_1, j_2), (j_2, j_1) ∈ E_1 AND (j_1, j_2) ∉ E_2.** Per Prop 1.5(4) = **freely independent**. Admissible partitions: π_0, π_a, π_b (not π_c because (3) requires the E_2 edge).

> φ = 0 + ϕ(X̃_{j_1}²)·0·0 + 0·ϕ(X̃_{j_2}²)·0 = **0**.

**Case C: (j_1, j_2) ∈ E_1, (j_2, j_1) ∉ E_1.** Per Prop 1.5(2): monotone (E_mono = E_1 \ E_1^op). Admissible: π_0, π_a (since π_b requires (j_2, j_1) ∈ E_1 which fails; π_c requires both, also fails).

> φ = 0 + ϕ(X̃_{j_1}²)·0·0 = **0**.

**Case D: (j_2, j_1) ∈ E_1, (j_1, j_2) ∉ E_1.** Per Prop 1.5(3): anti-monotone. Admissible: π_0, π_b.

> φ = 0 + 0·ϕ(X̃_{j_2}²)·0 = **0**.

**Case E: Boolean (neither edge in E_1).** Only π_0. φ = 0.

### 2.4 Summary of bigraph predictions

| Pairwise relation between A_{j_1}, A_{j_2} | Row (b) | Row (d) | Row (f) |
|---|---|---|---|
| Boolean | 0 | 0 | 0 |
| Monotone (j_1 → j_2) | 0 | 0 | 0 |
| Anti-monotone (j_2 → j_1) | 0 | 0 | 0 |
| Free | 0 | 0 | 0 |
| Tensor | 0 | 0 | **ϕ(X̃_{j_1}²) · ϕ(X̃_{j_2}²)** ≠ 0 |

**Bigraph also predicts row (d) = 0 for every pairwise relation.** Same MISMATCH with measured 0.108.

**Bigraph matches BMT tensor case at row (f).** This makes sense because BMT-tensor and bigraph-tensor are the same: tensor independence is the same notion in both frameworks (Prop 1.5(7) shows BMT independence is recovered when E_2 ⊇ (E_1 ∩ E_1^op) \ ∆).

---

## 3. Numerical check needed (row (f) tensor case only)

Both BMT (tensor digraph) and bigraph (tensor pairwise) predict row (f) = ϕ(X̃_{j_1}²) · ϕ(X̃_{j_2}²). Script `C:/Collatz/d2_extract_phi_X_squared.py` computes these under Reading B at the same V_TRUNC = 16, returning results in `C:/Collatz/experiments_output/d2_phi_X_squared.json`.

If the product ϕ(X̃_1²)·ϕ(X̃_2²) ≈ 6.089 × 10^{−1} at the sum_entries reduction, BMT/bigraph tensor matches row (f).

If the product is materially different, the framework's row (f) prediction is the wrong number (a quantitative failure even after the qualitative match noted by the auditor).

---

## 4. What can be concluded WITHOUT the numerical row (f) check

Rows (b) and (d) are settled analytically:

- **Row (b):** every BMT digraph and every bigraph predicts 0. Matches Syracuse ≈ 0. ✓
- **Row (d):** every BMT digraph and every bigraph predicts 0. **Syracuse measures 0.108.** **MISMATCH.**

The row (d) mismatch is decisive against ANY BMT or bigraph fit for Syracuse on a 2-vertex graph. Even if row (f) under tensor matches numerically, the framework still fails at row (d).

The only escape: use a LARGER vertex set V with auxiliary algebras, or a different categorical extension of the formula. But Syracuse's natural data is 2 vertices per row, so this is not a clean fit.

---

## 5. Files

- This file: `C:/Collatz/D2_BMT_BIGRAPH_PREDICTIONS.md`
- Companion verbatim: `C:/Collatz/D2_BMT_BIGRAPH_VERBATIM.md`
- Numerical script: `C:/Collatz/d2_extract_phi_X_squared.py` (sandbox-denied; ran by user externally if needed)
- Numerical output: `C:/Collatz/experiments_output/d2_phi_X_squared.json` (pending user run)
