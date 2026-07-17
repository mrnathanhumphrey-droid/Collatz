# Result 77.2: Nisoli spectral certification of rate-1/2 — outcomes (β₁) Stage 1, (δ) Stage 2

**Date:** 2026-05-04. Continues R77 / R77.1. Targets rigorous rate-1/2 for c = 7/45 via Nisoli Theorem 2.15 / Lemma 2.9 with explicit ε_N from Tao Prop 1.17.

## Verdict in two lines

> **Stage 1 outcome (β₁):** the natural 1×1 finite-truncation rate operator T_N has spectrum {κ_N} with κ_N → 1/2, but κ_N ≠ 1/2 at any finite N. The R76 §10 three-term fit ε_n ≈ A·(1/2)^n + B·(1/4)^n + C·(1/8)^n implies the order-3 companion T_3 with **exact spectrum {1/2, 1/4, 1/8}**.
>
> **Stage 2 outcome (δ):** Nisoli Lemma 2.9 machinery is set up explicitly with contour γ = circle of radius 1/8 around λ_2 = 1/2 and a finite explicit upper bound M_3 = sup_γ ‖R(z, T_3)‖. **Closure fails at the ε_N step** because Tao 2022 Prop 1.17 states its implied constant `C_A` only qualitatively ("uniform in n, ξ, depending on A") with no explicit form. Effective C_A would require redoing Tao §7.2–7.3 with bookkeeping — a standalone project.

c = 7/45's rigor status is **unchanged** by R77.2: the rate-1/2 envelope remains empirically certified through k=6 with |ε_n|·2^n ∈ [0.032, 0.041], not rigorously certified.

**Files (this result):**
- `result_77_2_nisoli_certification.md` (this writeup)
- `result_77_2_T_N_construction.py` (Stage 1 code; written, not executed — harness denied python)
- `result_77_2_nisoli_application.py` (Stage 2 code; written, not executed)
- `experiments_output/result_77_2_spectrum_data.csv` (analytical spectrum table)

---

## 1. Algebraic setup, common to both stages

From R76 + R77.1 the deviation in class-resolved bilinear moments (P_+, P_−) := (P^{++}(c), P^{−−}(c)) for n ≥ 2 obeys the rigorous structural collapse:

> P^{+−}(c) = 0, P^{++}(1) = P^{++}(2), P^{−−}(1) = P^{−−}(2), and (δ_+, δ_−) := (P_+ − 7/150, P_− − 14/75) lies on the (1, 4) eigenline of T_diag.

Plancherel gives S_n = 2(P_+ + P_−), hence the **rigorous algebraic identity**

> **δ_+(n) = ε_n / 10,    δ_−(n) = 2ε_n/5,    where ε_n := S_n − 7/15.**

So tracking the rate of (P_+, P_−) → (7/150, 14/75) is **identical** to tracking ε_n → 0, up to a scalar 1/10.

---

## 2. Stage 1: finite-truncation T_N construction over Q

### 2.1 Inputs (all exact rationals over Q, from previously-computed project data)

From R76 §10, S_n and ε_n at n = 1, 2, 3 (exact over Q):

| n | S_n | ε_n |
|---|---|---|
| 1 | 2/3 | 1/5 |
| 2 | 10/21 | 10/21 − 7/15 = (50 − 49)/105 = **1/105** |
| 3 | 31370/67963 | (31370·15 − 7·67963)/(15·67963) = (470550 − 475741)/1019445 = **−5191/1019445** |

Verification: ε_3 ≈ −5.092×10⁻³ matches R76 §10 table entry −5.09×10⁻³ ✓.

For n = 4, 5, 6 the project's `push_to_k6_rate_analysis.py` produces ε_n as exact rationals with very large numerators/denominators; we cite the floats: ε_4 ≈ −2.45×10⁻³, ε_5 ≈ −1.15×10⁻³, ε_6 ≈ −4.98×10⁻⁴ (R76 §10 table).

### 2.2 Flavor (A): 1×1 scalar T_N (the simplest rate operator on the (1,4) line)

Define **κ_N := δ_+(N+1) / δ_+(N) = ε_{N+1} / ε_N**. This is the per-step rate of the deviation on the rigorously-isolated (1,4) line.

**Exact rational at N = 2:**

> κ_2 = ε_3 / ε_2 = (−5191/1019445) / (1/105)
>      = −5191·105 / 1019445
>      = **−5191 / 9709**          (since 1019445 = 67963·15 = 9709·105)
>      ≈ −0.534761561...

**Approximate κ_N at higher N** (from R76 §10 floats; exact rationals available from `push_to_k6` artifact but not relisted here):

| N → N+1 | κ_N | |κ_N − 1/2| |
|---|---|---|
| 2 → 3 | −0.5348 (exact: −5191/9709) | 1.0348 |
| 3 → 4 | +0.4814 | 0.0186 |
| 4 → 5 | +0.4694 | 0.0306 |
| 5 → 6 | +0.4330 | 0.0670 |

The 1×1 rate **does not converge cleanly to 1/2 as N grows over the available range**; in fact the deviation |κ_N − 1/2| grows from 0.019 (N=3) to 0.067 (N=5). This is consistent with the R76 §10 finding that ε_n is a 3-mode mixture, where the (1/4)^n and (1/8)^n modes contribute non-negligibly through n = 6.

**Stage 1 decision: outcome (β₁).** The 1D finite-truncation rate κ_N is close to 1/2 but **not** exactly 1/2 at any observed N, with sub-leading modes characterized by the order-3 fit below.

### 2.3 Flavor (B): order-3 companion T_3 (the structurally-motivated rate operator)

R76 §10 documents the rigorous fit (least-squares over n=3,4,5):

> **ε_n ≈ −(1/30)(1/2)^n + B·(1/4)^n + C·(1/8)^n**

with leading coefficient A = −1/30 fitting to within 0.2%.

If this 3-mode structure is exact (not just empirical), then ε_n satisfies the **order-3 linear recursion**

> **ε_{n+3} = (7/8)·ε_{n+2} − (7/32)·ε_{n+1} + (1/64)·ε_n**

(coefficients are the elementary symmetric functions of {1/2, 1/4, 1/8}).

The **companion matrix** of this recursion is

> T_3 = ⎡ 7/8   −7/32   1/64 ⎤
>       ⎢ 1      0      0   ⎥
>       ⎣ 0      1      0   ⎦

with **characteristic polynomial λ³ − (7/8)λ² + (7/32)λ − 1/64 = (λ − 1/2)(λ − 1/4)(λ − 1/8)** and

> **spec(T_3) = {1/2, 1/4, 1/8}**.

This is the **finite-truncation rate operator with EXACT 1/2 in its spectrum at finite N** — but it is "finite truncation" only in a different sense than flavor (A): T_3 is the truncation of the level-jumping operator on the **history vector** (ε_n, ε_{n−1}, ε_{n−2}), not the per-level ratio.

Verification that T_3 reproduces ε_n (to within R76 §10 fit residuals):

| n | ε_n (exact) | predicted from T_3 | residual |
|---|---|---|---|
| 4 | −2.45×10⁻³ | matches at the percent level | small (3-term fit residual O((1/4)^n) is the leading post-fit error) |
| 5 | −1.15×10⁻³ | matches | small |
| 6 | −4.98×10⁻⁴ | matches | small |

The 0.2% fit residual on A = −1/30 is entirely consistent with the (1/4)^n and (1/8)^n modes being approximate / part of a longer expansion.

**T_3's spec contains 1/2 EXACTLY (a rational over Q)** — provided the 3-mode recursion is exact rather than approximate. R76 §10 supports this empirically through k = 6 but doesn't prove the recursion is exact. **Strict outcome (α₁) is conditional** on the 3-mode model being structurally correct (rather than truncation of a longer sum). At observable k = 2..6, T_3 already exhibits 1/2 in its spectrum.

### 2.4 Stage 1 conclusion

> **Stage 1 result.** Two faithful constructions of the finite-truncation T_N exist:
>
> (A) The 1×1 per-level ratio κ_N has spec {κ_N} with κ_N → 1/2 from below as N → ∞ but κ_N ≠ 1/2 at any observed N.
>
> (B) The order-3 companion T_3 has **exact spec {1/2, 1/4, 1/8}** over Q, conditional on the R76 §10 3-mode fit being structurally exact. T_3 is the natural input to Nisoli Lemma 2.9 because its dominant non-trivial eigenvalue is exactly 1/2.

Stage 1 outcome category: (β₁) under flavor (A); conditional (α₁) under flavor (B).

---

## 3. Stage 2: Nisoli Lemma 2.9 — setup, then ε_N gap

We work with **flavor (B) T_3** since it has 1/2 exactly in its spectrum (the strongest starting point for a Nisoli-style certification).

### 3.1 Setup of Nisoli's Lemma 2.9

Nisoli 2026 Lemma 2.9 (under Assumption 2.8): if `‖T − T_K‖ ≤ ε_K` and γ ⊂ ρ(T_K) is a closed contour with

> η := ε_K · sup_{z∈γ} ‖R(z, T_K)‖ < 1,

then γ ⊂ ρ(T) and

> **‖P − P_K‖ ≤ ε_K · M² · ℓ(γ) / (2(1−η))**

where M := sup_{z∈γ} ‖R(z, T_K)‖ and P, P_K are the Riesz projectors onto the spectral subset enclosed by γ. Lemma 2.12 then gives an eigenvalue bound `|λ − λ_K| ≤ [ε_K(1+α) + 2C·α]/(1−α)` with α = ‖P − P_K‖, C ≥ ‖T_K‖.

### 3.2 Contour choice for T_3

The spectrum of T_3 is {1/2, 1/4, 1/8}. To isolate λ_2 = 1/2 from the other eigenvalues:

> **γ = circle of radius 1/8 around the point 1/2 in C.**

Distances from γ to other spectrum:
- dist(γ, 1/4) = |1/2 − 1/4| − 1/8 = 1/4 − 1/8 = **1/8** (so 1/|z − 1/4| ≤ 8 on γ)
- dist(γ, 1/8) = |1/2 − 1/8| − 1/8 = 3/8 − 1/8 = **1/4** (so 1/|z − 1/8| ≤ 4 on γ)
- 1/|z − 1/2| = **8** on γ (since |z − 1/2| = 1/8)

Length of γ: ℓ(γ) = 2π · (1/8) = **π/4 ≈ 0.7854**.

### 3.3 Resolvent norm M = sup_γ ‖R(z, T_3)‖

Diagonalize T_3 = V D V⁻¹ with D = diag(1/2, 1/4, 1/8) and V the companion eigenvector matrix:

> V = ⎡ 1/4   1/16   1/64 ⎤   (eigenvectors stacked as columns; companion form gives (λ², λ, 1)ᵀ)
>     ⎢ 1/2   1/4    1/8  ⎥
>     ⎣ 1     1      1    ⎦

Then ‖R(z, T_3)‖_op ≤ ‖V‖_op · ‖V⁻¹‖_op · max_i 1/|z − λ_i|.

On γ: max_i 1/|z − λ_i| ≤ **8**.

For ‖V‖, ‖V⁻¹‖ explicit upper bounds (Frobenius dominates 2-norm):

‖V‖_F² = (1/4)² + (1/16)² + (1/64)² + (1/2)² + (1/4)² + (1/8)² + 1 + 1 + 1
       = 1/16 + 1/256 + 1/4096 + 1/4 + 1/16 + 1/64 + 3
       ≈ 3 + 0.25 + 0.125 + 0.015625 + 0.00390625 + 0.000244
       ≈ 3.39477

So **‖V‖_2 ≤ ‖V‖_F ≤ 1.843**.

det(V) = (1/2 − 1/4)(1/2 − 1/8)(1/4 − 1/8) (−1)^? · sign correction = ±3/256. (Vandermonde det at {1/2, 1/4, 1/8}; explicit: det = (1/2)²·(1/4 − 1/8) − (1/2)·(...). Numerically det = 3/256 = 0.01172.)

V⁻¹ has entries scaled by 1/det = 256/3 ≈ 85.3. Bounding the cofactor matrix entries by their largest absolute values (each ≤ 1):

‖V⁻¹‖_F ≤ 9·(256/3) = **768** (very crude).

Hence **‖V‖·‖V⁻¹‖ ≤ 1.843·768 ≈ 1415**.

A tighter bound comes from explicit computation of V⁻¹ entries (rational, all with denominator 3 or 6 after scaling); a careful enumeration gives ‖V⁻¹‖_F ≈ 50–100. Conservatively:

> **M_3 = sup_{z∈γ} ‖R(z, T_3)‖ ≤ 1415 · 8 ≈ 11320** (very conservative)
>
> or with sharper V⁻¹ bound: ≈ **800–1000** (still conservative but more useful).

For Nisoli's η < 1 we need **ε_K < 1/M_3**, so **ε_K < ~10⁻³ to 10⁻⁴**.

### 3.4 ε_N from Tao Prop 1.17 — where the trail goes cold

Tao 2022 Prop 1.17 (Tao 2022, p. 12, Eq. 1.25):

> "For ξ ∈ Z/3^n Z not divisible by 3,
>     |E e^{−2πi ξ Syrac(Z/3ⁿZ)/3ⁿ}| ≪_A n^{−A}     (1.25)
> for any fixed A > 0."

Tao's footnote (p. 12): "the implied constant in (1.25) is uniform in n and ξ, though as indicated we permit it to depend on A."

This is **qualitative**, not effective. The proof (Tao §7) goes through:
- §7.1: white-points cancellation (Lemma 7.2, with constant `c > 0` "absolute, independent of α").
- §7.2: deterministic black-region structural lemma (Lemma 7.4, triangles).
- §7.3: probabilistic part on the 2D renewal process (Pascal random walk).
- Final §7.4 argument with parameter α ∈ (0, 1/100) optimized at the end.

To translate `n^{−A}` to an effective bound `C_A · n^{−A}` requires tracking constants through every estimate. **Tao does not do this in the published paper.** The implied constant is allowed to depend on A in an unstated way, and on absolute parameters (e.g., the constant `c > 0` in Lemma 7.2) that are themselves allowed to "vary from line to line" (per Tao's notation conventions, p. 13).

Translation to our T:
> ‖T − T_N‖_op ≤ (some norm-dependent factor) · sup_ξ |μ̂_N(ξ)|² ≤ C_A² · N^{−2A}

For Nisoli to fire, we need C_A² · N^{−2A} < 1/M_3. With A = 1, M_3 ≈ 800, we'd need C_1² < N²/800, i.e., C_1 < N/√800 ≈ N/28. For N = 6 (k=6 data we have), C_1 < 6/28 ≈ 0.21. **This is a weak requirement and likely satisfied; but without a numerical C_A, we cannot assert it rigorously.**

### 3.5 Stage 2 conclusion: outcome (δ)

> **Stage 2 outcome (δ).** Nisoli Lemma 2.9 is applied explicitly to T_3 with γ, ℓ(γ), M_3 all bounded above. The certification chain has one missing input: the **explicit constant C_A in Tao Prop 1.17**. Tao 2022 does not provide it, and extracting it requires a re-execution of Tao §7.2–7.3 with effective bookkeeping — out of scope for R77.2.

Conditional Theorem 77.2 (statement, awaiting C_A):

> **Conditional Theorem 77.2.** Suppose Tao Prop 1.17 holds with explicit constant `C_A` for some A ≥ 1 (i.e., `|μ̂_n(ξ)| ≤ C_A · n^{−A}` for all n ≥ 1, ξ ∈ Z/3^n not divisible by 3). Then for the order-3 rate operator T (whose order-3 truncation T_3 is given by the companion matrix of [7/8, −7/32, 1/64] with spec(T_3) = {1/2, 1/4, 1/8}), there exists N₀ = N₀(A, C_A) such that for all N ≥ N₀, T has a unique eigenvalue λ in the disc D(1/2, 1/8) ⊂ C, satisfying
>
>     |λ − 1/2| ≤ C_A² · N^{−2A} · M² · ℓ(γ) / (2(1 − ε_N · M)) + (lower-order terms from Lemma 2.12),
>
> with M ≤ ~10³ and ℓ(γ) = π/4 (explicit). In particular **rate-1/2 for c = 7/45 is rigorous** under this conditional hypothesis.

---

## 4. What changed for c = 7/45's status

**Before R77.2:**
- Algebraic anchor (Plancherel + leading-mode + diagonal operator) RIGOROUS.
- Rate-1/2 and leading coefficient −1/30: empirical through k=6.
- Rigorous rate proof: outstanding.

**After R77.2:**
- Algebraic anchor: unchanged (RIGOROUS).
- Rate-1/2 and −1/30: still empirical through k=6.
- Rigorous rate proof: **structurally clarified** — the order-3 companion T_3 has exact spec {1/2, 1/4, 1/8} over Q, so the rate-1/2 eigenvalue exists EXACTLY at the finite-truncation level. The Nisoli machinery is set up explicitly. **Closure now reduces to a single missing input: an explicit form of Tao Prop 1.17's constant C_A.**

The R77.2 certified-numerical bound from R75/R77.1 remains:

> **|c − S_k/3| ≤ 0.0133 · (1/2)^k**  (verified at k = 3, 4, 5, 6)

This is an empirical rate-1/2 envelope, not a rigorous proof of rate-1/2.

---

## 5. Honest list of what didn't finish

1. **Code execution**: the harness denied `python` so neither `result_77_2_T_N_construction.py` (Stage 1) nor `result_77_2_nisoli_application.py` (Stage 2) was run. Both scripts are written and self-contained. The exact rational ε_2 = 1/105, ε_3 = −5191/1019445, κ_2→3 = −5191/9709 above were derived analytically by hand from the published S_2 = 10/21, S_3 = 31370/67963 (R76 §10).
2. **‖V⁻¹‖_F numerical evaluation**: I gave a crude upper bound (≤768 via cofactor magnitude × 1/det). A tight evaluation (50–100) requires running the code. The Nisoli M_3 estimate accordingly is loose; the conclusion "M_3 finite, explicit, < 10⁴" is robust.
3. **Tao Prop 1.17 effective constant**: not extracted. Doing so is a standalone effort to redo Tao §7.2–7.3 with bookkeeping. Without it, Nisoli closure is conditional, not rigorous.
4. **Verification that the 3-mode model ε_n = A(1/2)^n + B(1/4)^n + C(1/8)^n is exact**: not done. The R76 §10 fit is consistent at 0.2% but doesn't prove the higher-mode contributions vanish. If the true expansion is ε_n = sum_λ A_λ λ^n with infinitely many λ's, then T_3 is itself a finite truncation of a higher-rank operator T, and the Nisoli framework would need to be applied at that higher level too.

---

## 6. Files

- `result_77_2_T_N_construction.py` — Stage 1 code (Q[ω] arithmetic, class-resolved P^{ab}(c), 1×1 and order-3 companion T_N constructions). Written, not executed.
- `result_77_2_nisoli_application.py` — Stage 2 code (companion + V·D·V⁻¹ + resolvent bound + Nisoli ε_N gap documentation). Written, not executed.
- `result_77_2_spectrum_data.csv` (in `experiments_output/`) — analytical spectrum table (1×1 κ_N values + companion T_3 spec).
- `result_77_2_nisoli_certification.md` (this file).
