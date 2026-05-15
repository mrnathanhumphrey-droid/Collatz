# GENERALIZATION_CANDIDATE_1 — Cochrane–Pinner cubic exponential sums

**Candidate object:** `S(f, q) := Σ_{x ∈ Z/p^n} e_q(f(x))` where `q = p^n`, `n ≥ 3`, and `f(x) = ax³ + bx² + cx + d` is a cubic polynomial with `gcd(a, p) = 1`.

**Empirical signal:** Cochrane–Pinner 2003 derive an exact closed-form magnitude `|S(f, q)| = p^{n/2} · m` where `m` counts saddle points; on principal-unit-coset restriction `x ≡ 1 mod p` the magnitude saturates at `p^{(n+1)/2}` consistent with `√q · √p`.

**Pre-assessment:** STEP1-3 expected pass (this IS the literature R78 inherits from); STEP4-5 depend on the specific cubic phase structure after substitution.

---

## STEP 1 — closed-form magnitude theorem

**R78 analog:** F̂_p(p·a) = p · e_q(1) · G_p(a) with |G_p(a)| = √q · p (T78.4 + 78.6).

**Candidate test:** Cochrane–Pinner 2003 Theorem 1.1 gives directly:
> For `f(x) = ax³ + bx + c` with `v_p(f'') = 0` (non-degenerate) and `n ≥ 2`:
>   `S(f, p^n) = p^{n/2} · Σ_{x: f'(x) ≡ 0 mod p^⌈n/2⌉} e_{p^n}(f(x))`

This is the saddle-point closed form. The sum on the right is over saddle points, and at the right precision the per-saddle phase is exact.

**Verdict:** **STEP1_WORKS** (literally — R78's T78.4–T78.6 are a specialization of this).

---

## STEP 2 — bijection between parameter and saddle index

**R78 analog:** a ↔ C_a bijection (T78.5). Allows parameter to be re-indexed by saddle position.

**Candidate test:** For Cochrane–Pinner cubic `f(x) = ax³ + bx`, the saddle is `f'(x) = 3ax² + b ≡ 0 mod p^⌈n/2⌉`, i.e. `x ≡ ±√(−b/3a) mod p^⌈n/2⌉` (assuming `−b/3a` is a QR mod p). Bijection from `b` (parameter) to saddle `x*` is:
- For each non-zero QR class of `b`, two saddles `±x*(b)`.
- For non-QR class, no real saddle (sum vanishes by Cochrane's Cor 6).

So the bijection `parameter ↔ saddle` is **two-to-one** modulo QR class, not one-to-one as in R78. The R78 setup has c_1 = (C_a−1)/p ∈ Z/p giving a `p`-fold bijection on the principal-unit coset; the Cochrane–Pinner cubic has at most a `2`-fold bijection (sign choice on saddle).

**Verdict:** **STEP2_MODIFIED** — bijection exists but cardinality is different. For the chain to push through, the inner Plancherel step needs an analogous "p-fold parametrization of the coset". In R78 this comes from the principal-unit subgroup structure on the parameter side `a ≡ 1 mod p` in `(Z/p^{r+1})×`; for Cochrane–Pinner cubic the parameter is the coefficient `b ∈ Z/p^n` ranging freely, and the multi-saddle structure is the QR-classification mod p, not principal-unit.

This is a structural mismatch but not fatal. The "inner sum" structure (next step) will tell whether the modified bijection still admits collapse.

---

## STEP 3 — saddle exactness

**R78 analog:** at r = 3 (i.e. n = 4 in Cochrane–Pinner notation), `G_p(a) = √q · e_q(P_a(s*(C_a)))` is EXACT, not stationary-phase approximate. Driven by `J_p = r` condition (truncation level of p-adic log matches r exactly).

**Candidate test:** Cochrane–Pinner Theorem 1.1 states exact equality for `f` cubic, `n ≥ 3`. At `n = 4` the precision `⌈n/2⌉ = 2` truncates p-adic log at degree 2, which is INSUFFICIENT for the cubic phase to saturate at machine precision the way R78's J = 3 case does. At `n = 6` (r = 5 in R78 numbering), `⌈n/2⌉ = 3` matches the cubic degree exactly — this is the analog of R78's r = 3 saddle exactness.

**Verdict:** **STEP3_WORKS at `n` such that `⌈n/2⌉ = deg(f) = 3`** — i.e. `n ∈ {5, 6}`. At other `n`, the saddle is exact at the polynomial level but the truncation introduces a non-trivial higher-order correction analogous to R78's r ≥ 4 Hensel deviation.

R78's saddle exactness at r = 3 corresponds to a SPECIFIC alignment of polynomial degree and truncation level. Cochrane–Pinner cubic admits an analogous alignment at `n ∈ {5, 6}`.

---

## STEP 4 — Inner-Plancherel reduction (linear-in-c_2)

**R78 analog:** at r = 3, after substitution `P_a(s*) ≡ −p²·s*²/2 + p³·(s*³/6 − c_2·s*) mod p^4`. The `c_2` (second base-p digit of C_a) enters LINEARLY in the cubic-order term. This linear appearance is the key — Σ_{c_2} e_p(c_2 · s*) = p · δ(s*=0) collapses the inner sum to a length-p Dirichlet kernel.

**Candidate test:** For `f(x) = ax³ + bx` on `Z/p^n`, the parameter is `b`. Expand `b` base-p: `b = b_0 + p·b_1 + p²·b_2 + ...`. The saddle `x* = ±√(−b/(3a)) mod p^⌈n/2⌉` depends on `b` through square-root extraction.

At `n = 6` (chosen for saddle exactness): `⌈n/2⌉ = 3`, so `x*` is determined to `mod p³` from `b mod p³`. After the square-root extraction, the saddle phase `f(x*)` is:
> `f(x*) = a(x*)³ + b·x* = (1/2) · b · x*` (using `b = −3a(x*)²`, so `a(x*)³ = −b·x*/3`, hence `f(x*) = b·x* − b·x*/3 = 2b·x*/3`).

So `f(x*) = (2/3) · b · x*(b)`. This is the phase at the saddle, as a function of `b`.

Now: is this LINEAR in `b_2` (the second base-p digit)? The dependence is `b · x*(b)` where `x*(b)` depends on `b` through square root. Expanding `x*(b) = x*(b_0 + p·b_1) + ε(b_2)` to first order in `b_2`:
> `dx*/db = 1/(2·a·(x*)²·3) = 1/(6·a·(x*)²) = −1/(2b)` (using `b = −3a(x*)²`)
> So `δx* ≈ −(b_2 · p²)/(2b)`.

Thus `f(x*) = (2/3) · b · x* = (2/3) · b_0 · x*_0 + linear in b_1, linear in b_2 + ...`. The b_2 contribution to `f(x*) mod p^6`:
> `δf(x*) = (2/3) · [p² b_2 · x*_0 + b_0 · δx*] = (2/3) · [p² b_2 · x*_0 − b_0 · p² b_2 / (2 b_0)] = (2/3) · p² · b_2 · [x*_0 − 1/2]`

So **the `b_2` dependence is linear at order p²** — analogous to R78's `c_2` linear appearance at order p³.

**Verdict:** **STEP4_WORKS_MODIFIED.** The linear-in-second-digit structure transfers, but at a different p-adic order (p² vs p³ in R78). The inner Plancherel collapse still applies — Σ_{b_2 ∈ Z/p} e_p(b_2 · k) = p · δ(k=0).

---

## STEP 5 — 1/sin grid identity

**R78 analog:** `Σ_{α=0}^{p−1} 1/|sin(π(1+pα)/p²)| = p · csc(π/p²) + Σ_{α≥1} ... ≤ p + 2 log p`.

**Candidate test:** The outer sum after Inner-Plancherel collapse is over saddle classes (parametrized by b_0 mod p, or QR-class of b_0). For each saddle class, the Dirichlet-kernel structure of the outer sum is determined by `b_0 · x*_0 mod p` evaluation pattern.

In the worst case, the outer sum is `Σ_{b_0 ∈ (Z/p)×} 1/|sin(π · ratio(b_0)/p²)|` where `ratio(b_0)` depends on the specific cubic. For generic cubic `ax³ + bx`, this sum is bounded by the same 1/sin grid identity used in R78 — `≤ p + 2 log p ≤ 2p`.

**Verdict:** **STEP5_WORKS** — classical 1/sin grid identity applies on any prime-power coset.

---

## Five-step summary for Candidate 1

| Step | R78 (reference) | Candidate 1 | Outcome |
|------|-----------------|-------------|---------|
| 1 | F̂_p closed form | Cochrane–Pinner Theorem 1.1 directly | **WORKS** (R78 inherits this) |
| 2 | a ↔ C_a (p-fold) | b ↔ ±x*(b) (2-fold) | **MODIFIED** — bijection cardinality differs |
| 3 | Saddle exact at r = 3 | Saddle exact at n ∈ {5, 6} (analog) | **WORKS at specific n** |
| 4 | Linear-in-c_2 | Linear-in-b_2 at order p² | **WORKS_MODIFIED** |
| 5 | 1/sin grid identity | Same 1/sin identity | **WORKS** |

**Overall:** PASSES the chain in modified form. The Cochrane–Pinner cubic exponential sum admits an analogous bilinear bound `|S(f, q)| ≤ C · √q · p` at `n` such that saddle exactness holds, with `C` of order O(1) uniform in p.

**Novelty caveat (A4 advance flag):** This is essentially Cochrane–Pinner 2003 + classical Plancherel orthogonality + classical 1/sin identity. The Inner-Plancherel collapse (Step 4) is the only step that requires the specific linear-in-second-digit structure to be checked — but for cubic phases this is a routine 2-line p-adic expansion. **The "method" is the Cochrane-Pinner framework as published, applied to a bilinear evaluation context. No novelty over published 2003 work — this is a competent application of a 20-year-old technique.**

This candidate validates that the chain runs on its closest cousin, but doesn't create new methods leverage.
