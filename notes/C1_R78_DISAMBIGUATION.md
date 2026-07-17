# C1 Phase 4 — R78 D=0 obstruction disambiguation

**Date:** 2026-05-12
**Sources:** `C:/Collatz/result_78.md`, `result_78_FINAL.md`, `cochrane.txt`.

---

## 4(a) What was the specific setting in R78?

R78 attempted to bound the **partial / incomplete sum**:

  S_partial = Σ_{u=0}^{3^{r-1} - 1} e_{3^{r+1}}( c · 4^u − 9 m u )

via Cochrane's **Theorem 2 + Corollary 6** machinery — the *coefficient-extraction / complete-sum-vanishing* theorem, NOT Theorem 1.2 (the mixed-sum upper bound).

From `result_78.md`:
- "Cochrane Theorem 2 with D = 0 gives 'sum vanishes by counting H(a) ≡ 0 condition' — but this is a **complete-sum result**. Our incomplete sum (range 3^{r-1}, modulus 3^{r+1}) does NOT inherit this trivial vanishing."
- "**Cochrane's machinery is designed to extract cancellation from polynomial structure of the character / phase modulo p.** The Syracuse map's specific 3-adic structure of (1+3)^u places ALL cancellation at the constant level after scaling by 3^τ — the higher-derivative polynomial behavior is 'trivial mod 3.' Cochrane's bound on D = 0 sums is the trivial bound."

The polynomial identification was: f(u) = 1 (trivial additive part), g(u) = c · Σ_{k=0}^r C(u, k) 3^k − 9mu (degree r polynomial in u). Theorem 2 / Corollary 6 give vanishing-of-complete-sum when the **degree mod p** (denoted D = deg_p H+ in Cochrane's notation) is non-zero. For R78, **D = 0** because after factoring τ = 1, H(u) mod 3 has only its constant term non-zero.

R78_FINAL extended this to **Theorem 78.1 (Complete-sum vanishing) for all (r, ℓ, ε, m)** — verified to 10^{-14} machine precision. So R78 closed *complete-sum* vanishing affirmatively; the obstruction was at the **incomplete → complete translation**, where Pólya-Vinogradov gave only a weak bound that didn't beat the trivial length-N bound for r ≥ 3.

## 4(b) Does the D=0 obstruction propagate to Theorem 1.2's setting?

**Cochrane Theorem 2** is a coefficient-extraction / vanishing-condition theorem (the "S_α = 0 if α not a critical point" part of Theorem 1.1, and its degree-of-derivative variant Cor 6). It tells you the sum vanishes when D ≠ 0; when D = 0, you fall back to the *trivial* bound.

**Cochrane Theorem 1.2** is the *upper-bound* theorem (the (1.14)–(1.17) bounds). Its hypothesis is structural (m ≥ t+2, χ multiplicative, f, g rational of degrees d_1, d_2). It does NOT directly require D ≠ 0; rather, the bound (1.17) is:

  |S(χ, g, f, p^m)| ≤ 4.41 · p^{m(1 − 1/(d_1+d_2))}.

This is non-trivial *as long as d_1 + d_2 is bounded uniformly in m*.

**Does R78's D=0 phenomenon propagate to Theorem 1.2?** *Indirectly yes, via a different mechanism:*

In the natural R78-style polynomial identification, g(u) = c · Σ_{k=0}^r C(u, k) 3^k − 9mu has degree d = r in u. For r = level = n (or some constant offset), d ≍ n. Then Cor 1.1's bound becomes:

  4.41 · (3^{r+1})^{1 − 1/r} = 4.41 · 3^{(r+1)(r−1)/r} = 4.41 · 3^{r − 1/r} ≍ 3^r

which is **the trivial bound** (the sum has length 3^{r+1}, and Cor 1.1 gives 3^r — only a factor of 3 saved). The saving 4.41 · p^{m(1 − 1/(d_1+d_2))} requires d_1 + d_2 bounded as p^m → ∞; for Syracuse the natural degree d_1 + d_2 grows linearly with m, so Theorem 1.2's bound *degenerates to trivial* in the limit n → ∞.

This is a **structurally different obstruction from R78's D=0**: R78's D=0 is *the polynomial structure becomes trivial mod 3 after factoring τ=1*; Theorem 1.2's failure is *the polynomial degree grows with n, killing the saving exponent*.

Both obstructions trace back to the **same root cause**: the natural polynomial identification of (1+3)^u (or 4^u = (1+3)^u via 4 = 1+3) gives a polynomial whose p-adic structure is "trivial in the limit" — coefficients have v_3 ≥ 1 except the constant (R78's D=0), AND degree grows with level r (Theorem 1.2's degree blowup).

## 4(c) Disposition for Cochrane route

- **R78's D=0 (Theorem 2): real obstruction, closed for the polynomial-identification approach.**
- **Theorem 1.2 / Cor 1.1: separately obstructed by degree blowup d ≍ n, gives trivial bound.**

Neither closes Syracuse. The two failures are **structurally distinct but trace to the same algebraic root** — 4 = 1 + 3 in Z_3 makes the natural polynomial expansion p-adically trivial at the top, p-adically degree-large at the lift.

**Could a different (non-polynomial) identification rescue Cochrane?** Phase 3 (3b.i, 3b.ii, 3b.iii) showed: at the *full μ̂_n* level, the weight is not a multiplicative character; at the *conditional f(x, b)* level, the sum is sub-trivially short (2 terms for b=3); at the *one-step recursion* level, A_v is a single phase, not a sum. **There is no level at which a non-trivial Cochrane S(χ, g, f, 3^n) form emerges naturally.**

## 4(d) Audit against the inherited-claim trap (Mode E)

The brief explicitly asked: "don't assume R78 closed Cochrane entirely without verifying which theorem and which setting was tested."

Verified:
- R78 tested Theorem 2 (coefficient extraction) on a specific Kalafatelis polynomial (1+3)^u family at modulus 3^{r+1}.
- R78 did NOT test Theorem 1.2 directly. The 1.2 inheritance is via the degree-blowup analysis here (Phase 3b.iv, Phase 4b).
- R78's D=0 is **not** inheritable wholesale — Theorem 1.2 fails on a different (degree) mechanism.

**But both failures route Cochrane-as-a-whole to the SAME closed disposition**: Cochrane's prime-power machinery is structurally not the right tool for Syracuse's renewal-walk-product form. The closure is not "R78 closed Cochrane" but rather "**every natural way to force Syracuse into Cochrane's setup hits a structural obstruction** — either D=0 (Theorem 2) or degree blowup (Theorem 1.2) or domain-category mismatch (full μ̂_n vs Z/p^m complete sum)."

## Strategic takeaway

R78's outcome (γ → β-strong) opened **Fourier-sparsity structure of e_q(c · 4^u)** as a new R76/R77-compatible structural lemma (R78.1–R78.3). That structural by-product is genuinely useful for the c=7/45 closure arc. But it does **not** unlock the polynomial-in-A bound on |μ̂_n(ξ)|.

For the polynomial-in-A bound specifically, Cluster 1 (Cochrane/BC/HB/Konyagin) is **closed-negative** by combination of:
- Phase 3a, 3b, 3c, 3d failures (this probe)
- R78 Theorem 2 D=0 obstruction
- Theorem 1.2 degree blowup (this probe, Phase 3b.iv)
