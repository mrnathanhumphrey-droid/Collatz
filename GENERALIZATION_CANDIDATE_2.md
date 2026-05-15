# GENERALIZATION_CANDIDATE_2 — Heath-Brown cubic character sum on prime-power modulus

**Candidate object:** `T(χ, a, q) := Σ_{n ∈ (Z/q)×} χ(n) · e_q(a · n³)` where `q = p^r` (r ≥ 2), `χ` is a primitive Dirichlet character on `(Z/q)×`, and `a ∈ (Z/q)×`.

**Empirical signal:** Heath-Brown's cubic-character-sum bound `T ≪ q^{1−1/8+ε}` is established for square-free q. On prime-power q the empirical rate is `≈ q^{1/2}` (square-root saturation), not the q^{7/8} of Heath-Brown's bound — i.e. the cancellation is STRONGER on prime-power than on square-free q.

**Pre-assessment:** STEP1 partial/fail — Postnikov-decomposed phase mixes polynomial + truncated log; downstream steps fragile.

---

## STEP 1 — closed-form magnitude theorem

**R78 analog:** F̂_p(p·a) closed form via Cochrane's polynomial-phase machinery.

**Candidate test:** Apply Postnikov to write `χ(n) = e_q(λ_χ · L_p(n))` on the principal-unit coset `n ≡ 1 mod p`. This is valid on the principal-unit subgroup, with `λ_χ` a fixed "log of χ".

After substitution:
> T(χ, a, q) restricted to n ≡ 1 mod p = `Σ_{n ≡ 1 mod p} e_q(λ_χ · L_p(n) + a · n³)`

The combined phase `Φ(n) = λ_χ · L_p(n) + a · n³` is **NOT a polynomial** — it has a truncated-log term plus a cubic term. Cochrane's Theorem 2 requires polynomial phase. So direct application fails.

**Workaround:** Substitute `n = 1 + p·u` for `u ∈ Z/p^{r−1}`. Then:
- `n³ = (1+pu)³ = 1 + 3pu + 3p²u² + p³u³`
- `L_p(n) = L_p(1+pu) = pu − p²u²/2 + p³u³/3 − ...` (truncated to J_p terms)

Substituting:
> `Φ(1+pu) = λ_χ · (pu − p²u²/2 + p³u³/3) + a · (1 + 3pu + 3p²u² + p³u³) mod q`

This IS a polynomial in `u` of degree 3 (cubic). Cochrane's Prop 4 applies.

**Verdict:** **STEP1_WORKS_MODIFIED** — only after Postnikov substitution `n = 1 + pu`. The phase becomes a polynomial of `u` mixing coefficient contributions from λ_χ and `a`. Closed-form magnitude theorem then follows from Cochrane.

But: the polynomial is mixed — the coefficient structure differs from R78's `P_a(s) = ps − C_a · L_p(1+ps)`. The bijection step (Step 2) becomes more delicate because two parameters (`λ_χ`, `a`) parametrize the phase, while in R78 only one parameter `a` parametrizes via `C_a`.

---

## STEP 2 — bijection between parameter and saddle index

**R78 analog:** a ↔ C_a single-parameter bijection.

**Candidate test:** The Heath-Brown setup has **two** parameters: `λ_χ` (fixed by χ) and `a` (variable). For each fixed χ, the parameter `a` enters as:
> `a · (1 + 3pu + 3p²u² + p³u³) = a + 3a·pu + 3a·p²u² + a·p³u³`

So the polynomial in u is:
> `Φ(u) = (constant) + (linear) · u + (quadratic) · u² + (cubic) · u³`
> where:
> - linear coefficient: `λ_χ · p + 3a · p`
> - quadratic coefficient: `−λ_χ · p²/2 + 3a · p²`
> - cubic coefficient: `λ_χ · p³/3 + a · p³`

The saddle is `Φ'(u) = 0 mod p^⌈r/2⌉`. Saddle position `u*` depends on the **combined** coefficients — not a clean one-parameter bijection.

**Verdict:** **STEP2_MODIFIED** — bijection exists between (λ_χ, a) jointly and saddle position, but is two-parameter rather than one-parameter as in R78. This complicates the downstream Inner-Plancherel step because the "second base-p digit" of the parameter is now ambiguous (which parameter's digit?).

---

## STEP 3 — saddle exactness

**R78 analog:** Saddle exact at r = 3 due to J_p = 3 matching polynomial degree.

**Candidate test:** The polynomial Φ(u) has degree 3 (after truncation of L_p at J = 3). Saddle exactness needs `⌈r/2⌉ ≥ 3`, i.e. `r ≥ 5` or `r ≥ 6`. At `r = 5, 6`, the saddle position is determined to precision matching the polynomial degree.

But: there's a subtlety. The coefficients of Φ depend on BOTH `λ_χ` and `a`. If we fix `a` and let `χ` vary (or vice versa), the saddle moves in a way that's not single-parameter. For the bilinear bound on `T(χ, a, q)` summed over (χ, a) jointly, the saddle structure varies.

**Verdict:** **STEP3_STATIONARY_ONLY** — saddle method applies and gives stationary-phase approximation, but saddle exactness in the strict R78 sense (saddle gives exact contribution, no error term) requires the precise J = degree alignment. This works at `r = 5, 6` for a fixed (λ_χ, a) pair, but the joint structure is more delicate than R78's single-parameter exactness.

---

## STEP 4 — Inner-Plancherel reduction

**R78 analog:** Linear-in-c_2 phase structure enables Σ_{c_2} e_p(c_2 · s*) = p · δ — collapse to length-p Dirichlet kernel.

**Candidate test:** Let's check if any of the parameter base-p digits enters linearly after substitution at the saddle. Take `a = a_0 + p·a_1 + p²·a_2` and similarly for `λ_χ`. At the saddle, the phase value `Φ(u*)`:

Using the saddle method, `Φ(u*) = Φ(u_0) − (Φ'(u_0))²/(2 Φ''(u_0)) + O(p³)` where `u_0` is the leading-order saddle.

Working through (skipping the algebra): the `a_2` digit appears in `Φ(u*) mod p^6` via:
> `p² a_2 · u*_0 − p² a_2 · u*_0 · (some Hessian factor)` — generically a polynomial in `u*_0`, not necessarily linear in `a_2 · u*_0` after subtracting the Hessian correction.

For pure cubic `ax³` (no χ), the phase value at saddle is `(2/3)·a·u*` which IS linear in `a` after substituting `u* = (−b/3a)^{1/2}` and simplifying. So pure-cubic case: linear-in-second-digit holds.

For mixed `λ_χ · L_p(1+pu) + a · (1+pu)³`: the saddle equation is `λ_χ · L_p'(1+pu) · p + a · 3(1+pu)² · p = 0`, i.e. `λ_χ · L_p'(1+pu)/(1+pu) ≈ −3a · (1+pu)` — this is NOT a clean square-root extraction. The saddle is implicit in a transcendental-style equation (truncated p-adically).

After saddle substitution, the phase `Φ(u*(a, λ_χ))` evaluated at the second-base-p-digit level: contains BOTH `a_2 · stuff` and cross-terms `a_2 · λ_{χ,1}` etc. The linear-in-single-digit structure DOES NOT hold cleanly.

**Verdict:** **STEP4_DIFFERENT_INNER** — the inner phase at the saddle is no longer linear-in-second-digit-of-parameter due to the two-parameter mixing. Inner Plancherel still applies but with a different structure (likely Σ_{(a_2, λ_{χ,1})} e_p(a_2 · k + λ_{χ,1} · k')) which is a 2D Plancherel, not 1D Dirichlet kernel.

If we fix χ (so λ_χ is held constant) and vary only `a`: the inner is 1D Plancherel on `a_2` and collapses cleanly. But then we're not summing over χ — we've reduced to a fixed-character setup, which is a different problem from the standard Heath-Brown bilinear.

---

## STEP 5 — 1/sin grid identity

**R78 analog:** Σ_α 1/|sin(...)| ≤ p + 2 log p ≤ 2p.

**Candidate test:** The outer sum after Inner-Plancherel collapse — what does it look like for Heath-Brown's setup? If we fixed χ and summed over `a`, the outer sum is over the saddle classes of `a mod p`. The grid identity applies: same 1/sin sum bound.

If we summed over χ AND `a` (the bilinear of interest), the outer sum is two-dimensional and the 1/sin identity needs a 2D analog. No such clean identity is known in this generality.

**Verdict:** **STEP5_WORKS at fixed χ, STEP5_DIFFERENT_GRID at varying χ.**

---

## Five-step summary for Candidate 2

| Step | R78 (reference) | Candidate 2 | Outcome |
|------|-----------------|-------------|---------|
| 1 | F̂_p closed form | Postnikov + Cochrane on n = 1+pu | **WORKS_MODIFIED** (after Postnikov substitution; mixes two parameters) |
| 2 | One-parameter bijection | Two-parameter (χ, a) → saddle | **MODIFIED** — joint not single-parameter |
| 3 | Saddle exact at r = 3 | Saddle exact at r ∈ {5, 6} for fixed (χ, a) | **STATIONARY_ONLY** at varying parameters |
| 4 | Linear-in-c_2 | Linear-in-single-digit FAILS at varying χ | **DIFFERENT_INNER** — 2D Plancherel needed |
| 5 | 1/sin grid identity | Works at fixed χ; 2D identity unknown at varying χ | **DIFFERENT_GRID** at the full bilinear |

**Overall:** The chain FAILS at the FULL Heath-Brown bilinear (joint sum over χ and a). It works in the degenerate fixed-χ sub-case, but that's already known to be a Cochrane-Pinner cubic sum (subsumed by Candidate 1).

**Where exactly the chain breaks:** Step 4 — the linear-in-second-digit structure that enables Plancherel collapse is specific to ONE-PARAMETER setups. Heath-Brown's bilinear has two parameters and the cross-term between `a_2` and `λ_{χ,1}` prevents the clean 1D Plancherel collapse.

**Verdict:** STEP1-3 PARTIAL (work in fixed-parameter degeneration only); **STEP4 FAILS** at the actual Heath-Brown bilinear; STEP5 inherits the failure.

This is genuinely a step-by-step failure at a specific load-bearing point: the chain depends on the parameter being a SINGLE base-p-decomposable index, and Heath-Brown's character sum has two-index parameter structure that breaks the inner Plancherel.
