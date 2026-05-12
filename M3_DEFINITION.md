# M3_DEFINITION — explicit definition of M_3 in project notation

**Date:** 2026-05-11. Phase 1 of `M3_*` probe. Articulates what M_3 *is* in the R77.2 Nisoli framework, what operator it is built from, and what is/isn't known about it numerically — including the load-bearing impact of the R77.3 falsification on the spectrum that M_3 was originally computed against.

## 1. M_3 in project notation

From `result_77_2_nisoli_certification.md` §3.1–3.3 (Nisoli Lemma 2.9 setup):

> **M_3  :=  sup_{z ∈ γ}  ‖R(z, T_3)‖_op**
>
> where:
> - `T_3` is the order-3 companion matrix of the (conjectured) recursion `ε_{n+3} = (7/8)ε_{n+2} − (7/32)ε_{n+1} + (1/64)ε_n` (R77.2 §2.3 flavor B),
> - `R(z, T_3) := (z·I − T_3)^{-1}` is the resolvent of `T_3` at the complex parameter `z`,
> - `γ` is the closed contour `|z − 1/2| = 1/8` (circle of radius 1/8 around the dominant eigenvalue λ₂ = 1/2),
> - `‖·‖_op` is the operator 2-norm on C^3.

Subscript "3" refers both to the dimension of T_3 (3×3) and to the radius selection (1/8 vs 1/4 separation from the next eigenvalue) — they coincide because the conjectured spectrum {1/2, 1/4, 1/8} has dimension 3.

M_3 is the **operator-norm resolvent ceiling** on the Riesz-projection contour. It is the load-bearing constant in Nisoli Lemma 2.9's closure inequality.

## 2. T_3 — domain, action, conjectured spectrum, falsification

### 2.1 Domain and action (R77.2 §2.3)

`T_3 : C^3 → C^3` acts on the "history vector" `v_n := (ε_n, ε_{n−1}, ε_{n−2})ᵀ` by

> T_3 = ⎡ 7/8   −7/32   1/64 ⎤
>       ⎢ 1      0      0   ⎥
>       ⎣ 0      1      0   ⎦

The first row encodes the recursion coefficients (elementary symmetric functions of the conjectured eigenvalues {1/2, 1/4, 1/8}); the bottom two rows shift history. Acting on `v_n` produces `(ε_{n+1}, ε_n, ε_{n−1})ᵀ = v_{n+1}`.

### 2.2 Conjectured spectrum (R77.2)

> spec(T_3) = {1/2, 1/4, 1/8}    [R77.2 §2.3, equivalent to companion-form roots of `λ³ − (7/8)λ² + (7/32)λ − 1/64 = (λ−1/2)(λ−1/4)(λ−1/8)`]

This is the spectrum **assuming the 3-mode model ε_n = A(1/2)^n + B(1/4)^n + C(1/8)^n is exact in Q**. R77.2 §2.3 flagged this as "conditional on the 3-mode fit being structurally exact rather than truncation of a longer sum."

### 2.3 R77.3 falsification of the spectrum

From `result_77_3_nisoli_bypass.md` §3–§5 (and confirmed by `result_77_4_K_spectrum_erratum.md`):

> **The 3-mode model is FALSIFIED in Q.** Solving for `(A, B, C)` exactly over Q from `{ε_1, ε_2, ε_3}` gives `A = −157462/3058335 ≈ −0.05148`, **not** `A = −1/30 ≈ −0.0333`. Predictions from this exact `(A, B, C)` miss the actual `ε_n` at `n = 4, 5, 6` by **28%, 33%, 41%** relative — orders of magnitude above any rounding tolerance.

R77.3 §5 confirms independently via the direct recursion test at `n=1`: predicted `ε_4 = −222733/65244480 ≈ −3.41×10⁻³`; actual `ε_4 ≈ −2.45×10⁻³`. Residual ≈ +0.96×10⁻³, ~28% relative, nonzero in Q.

R77.3 §6 (B.2): the **empirical** order-3 companion fit to actual data has approximate spectrum `{0.534, 0.144, −0.084}`, **not** `{1/2, 1/4, 1/8}`. The clean rational eigenvalues that motivated M_3's contour are an **artifact of the falsified conjecture**, not a property of any underlying operator.

R77.4 erratum §"Empirical evidence" hammers this further: the actual within-level Markov transition `K_k` has no eigenvalue anywhere near 1/2 at any `k ∈ {3,4,5,6}` — closest eigenvalue at all levels is `λ₂ ≈ 0` with `|λ_2 − 1/2| ≈ 0.5`. So the "1/2 in the spectrum" is also absent from the natural operator on the state space.

### 2.4 What is T_3 actually, post-R77.3?

After R77.3 + R77.4 erratum:

- **T_3 as defined in R77.2 (companion of `(7/8, −7/32, 1/64)`)** still has the literal spectrum `{1/2, 1/4, 1/8}` — that's a mathematical fact about that particular 3×3 matrix.
- But T_3 **no longer corresponds to any operator-theoretic statement about ε_n's actual dynamics**. R77.3 shows ε_n does **not** satisfy that recursion; R77.4 shows the natural level-k transition has spectrum nowhere near 1/2.
- M_3 = sup_γ ‖R(z, T_3)‖ is therefore the resolvent norm of a 3×3 matrix that was a **proposed** model, now falsified.
- The "true" T (if any well-defined operator exists whose spectrum captures the (1/2)^n envelope) is **uncharacterized**. R77.4 erratum §"What this DOES change" lists candidates (inter-level residual operator, generating-function singularities) as *parked pending direction*.

## 3. The Nisoli closure inequality, three constants

From `result_77_2_nisoli_certification.md` §3.1, Nisoli Lemma 2.9 closure requires:

> **η  :=  ε_K · M_3  <  1**    (the resolvent–perturbation closure condition)

where `ε_K := ‖T − T_K‖_op` is the **operator-norm distance** from the full T to the finite truncation T_K. From `PRECISE_ASK.md` §4 the project relates `ε_K` to the bilinear sum `|K|`:

> roughly `ε_K  ~  |K(K)| / √q`    (via Plancherel + R76 §10 decomposition; from PRECISE_ASK.md §4)

Together the three constants enter the closure as:

> **|K| · q^{−1/2} · M_3  <  1**

In the user's task statement this is written `|K| · K^{−A} · M_3 < 1`. The mapping is:

- `|K|` = bilinear sum bound at level K (= `K(r=K−1, c=1, m=0)` magnitude),
- `K^{−A}` = `q^{−1/2}` factor → with `q = 3^{r+1} = 3^K`, this is `3^{−K/2}`. Reparameterizing the base, `K^{−A}` is a placeholder for Tao Prop 1.17's `n^{−A}` decay of `|μ̂_n(ξ)|`,
- `M_3` = sup_γ ‖R(z, T_3)‖ as defined above.

**Strict form (R77.2 §3.5 / PRECISE_ASK §4):**

> Need  `ε_K  <  1/M_3`    i.e.  `|K(K)| / √q  <  1/M_3`    i.e.  `|K(K)|  <  √q / M_3`.

With `M_3 ≈ 800–1000` and at K=6 (q=3^7=2187, √q ≈ 46.8): need `|K| < 46.8/1000 ≈ 0.047` to fire Nisoli. Empirically `|K_max| ≈ 16.58` at K=6 — three orders of magnitude too large.

## 4. What's known about M_3 numerically

### 4.1 R77.2's estimate (built on falsified spectrum)

R77.2 §3.3 derives, **assuming spec(T_3) = {1/2, 1/4, 1/8}**:

- On γ = `|z − 1/2| = 1/8`:
  - `1/|z − 1/2| = 8` (exact, on the circle)
  - `1/|z − 1/4| ≤ 8` (closest approach 1/8)
  - `1/|z − 1/8| ≤ 4` (closest approach 1/4)
  - `max_i 1/|z − λ_i| ≤ 8`
- Diagonalization T_3 = V D V⁻¹ with V = Vandermonde at {1/2, 1/4, 1/8}:
  - `‖V‖_F ≤ 1.843` (Frobenius dominates 2-norm)
  - `det V = ±3/256` ≈ 0.0117
  - `‖V⁻¹‖_F ≤ 768` (very crude cofactor bound) or **50–100** (sharper enumeration; R77.2 §3.3)
- M_3 = sup_γ ‖R(z, T_3)‖ ≤ `‖V‖ · ‖V⁻¹‖ · max_i 1/|z − λ_i|`:
  - **Crude:** 1.843 × 768 × 8 ≈ **11320**
  - **Sharper:** 1.843 × (50..100) × 8 ≈ **800–1500**

R77.2 §3.3 quotes "**M_3 ≈ 800–1000**" as the best loose estimate; PRECISE_ASK.md §4 inherits this.

### 4.2 Status under R77.3 falsification

The numerical bound M_3 ≈ 800–1000 is **computed against a falsified spectrum**. If T_3 means "the companion of `(7/8, −7/32, 1/64)`" as a literal 3×3 matrix, then M_3 ≈ 800–1000 stands as an algebraic fact about that matrix. **But that matrix does not describe ε_n's actual dynamics**, so the Nisoli closure built on top is hollow: η < 1 for that T_3 doesn't certify anything about the actual T (which is undefined operator-theoretically).

### 4.3 What would need to be true for M_3 to mean something

For Nisoli closure to actually fire for c=7/45, **all three** of the following must hold:

1. There exists a well-defined finite-rank operator T_K converging to T in operator norm with explicit `ε_K → 0` rate (Tao Prop 1.17, gated on effective C_A — Route 1 INFEASIBLE this session).
2. T_K has a spectral feature near 1/2 that's isolatable on a contour γ with `sup_γ ‖R(z, T_K)‖` explicitly finite (R77.3 falsifies the candidate 3-mode T_3; R77.4 erratum rules out the natural K_k).
3. The closure inequality `ε_K · M_K · ℓ(γ) / (2(1−η)) < (target precision)` holds with explicit numbers.

Items (1) and (2) are both currently open. Item (3) is moot until (1) and (2) are resolved.

### 4.4 Other project documents on M_3

- `result_77_2_nisoli_certification.md` §3.3, §3.5 — only place M_3 is computed.
- `PRECISE_ASK.md` §4 — quotes M_3 ≈ 800–1000 from R77.2, uses it to derive required `|K| < 0.05` at K=6.
- `result_77_3_nisoli_bypass.md` — falsifies the spectrum M_3 was built against, doesn't re-derive M_3.
- `result_77_4_K_spectrum_erratum.md` — rules out K_k as the operator M_3 should have been built from.
- No other project document gives a numerical value for M_3.

## 5. Summary for Phase 2 entry

Going into Phase 2 (extraction of M_3 via three approaches), the entering state is:

- **Definition:** M_3 = sup_{z ∈ γ} ‖R(z, T_3)‖_op with T_3 a specific 3×3 matrix and γ = circle of radius 1/8 around 1/2.
- **Spectrum status:** spec(T_3) = {1/2, 1/4, 1/8} as algebraic fact about the matrix; but this spectrum is **not the spectrum of any operator describing ε_n** (R77.3 falsification). There is **no currently-characterized operator** T whose spectrum captures the (1/2)^n envelope (R77.4 erratum lists candidates as parked).
- **Numerical bound:** R77.2 gives `M_3 ∈ [800, 1500]` (crude 11320 also derivable, sharper 800–1000) **for the falsified T_3**.
- **Normality:** T_3 is companion-form, hence **non-normal** in general. The spectral-radius lower bound (8 on γ) is therefore not tight for the operator norm; the V·V⁻¹ multiplier (50..768) reflects the non-normality.

This sets up Phase 2 Approach A (direct spectral) as a **lower bound only** for the operator norm, and Approach C (numerical) as the natural route — though either approach is **moot for the closure question** because the underlying T_3 doesn't describe ε_n's dynamics.

The load-bearing obstruction is therefore at the **definition layer**, not the computation layer.
