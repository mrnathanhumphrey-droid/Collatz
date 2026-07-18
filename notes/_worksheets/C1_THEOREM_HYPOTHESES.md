# C1 Phase 2 — Cluster 1 theorem hypotheses, verbatim

**Date:** 2026-05-12

---

## 2.1 Bourgain–Chang Theorem 4.7 (composite moduli, multiplicative subgroup)

Source: `C:/tmp/crystal/cochrane_mccarthy.txt` (this is BC composite-moduli, not Cochrane-McCarthy — filename misleading per brief).

**Setup (p. 42, §4 lead-in, lines 1936–1942):**

> "Let q = ∏_{α=1}^β p_α^{ν_α} ∈ ℤ_+. We say q has **few prime factors**, if
> ∑_{α ≤ β} ν_α < C_0 for some constant C_0."

**Theorem 4.7 (p. 42, lines 2173–2184):**

> "Let q ∈ ℤ_+ have few prime factors and let H < Z_q^* satisfy |H| = q^δ. Then
> max_{ξ ∈ Z_q^*} | Σ_{x ∈ H} e_q(ξ x) | < q^{-ε} |H|.    (4.21)
> where ε > 0 depends on the numbers of prime factors of q and δ > 0 only."

**Earlier scaffolding (Introduction, p. 4, lines 90–95):**

> "we start investigating … in the case of composite moduli q = p_1^{α_1} ⋯ p_r^{α_r}. We assume q has only a bounded number of prime factors p_1, …, p_r, which are moreover **'large'**, i.e. p_i > q^ε for some ε > 0 (hence α_1 + ⋯ + α_r < 1/ε)."

**Verdict on "few prime factors":** the constraint is **two-pronged** — bounded number of distinct primes AND the exponents (ν_α / α_i) themselves bounded so that the primes are "large" relative to q (p_i > q^ε ⇒ α_i < 1/ε). For q = 3^n with α_1 = n → ∞: the *number* of distinct primes is 1 (bounded), but the *exponent* α_1 = n violates p_1 > q^ε (since 3 > 3^{nε} fails for any fixed ε > 0 as n grows). **Hence q = 3^n with growing n VIOLATES BC 4.7's "few prime factors" hypothesis.**

**Conclusion structure (the additive-character sum on H):** the bound is on Σ_{x ∈ H} e_q(ξ x) — a **complete sum over a multiplicative subgroup H of Z_q^***, using the *additive* character e_q.

---

## 2.2 Cochrane Theorem 1.2 and Corollary 1.1 (mixed exp sums, prime power)

Source: `C:/tmp/crystal/cochrane.txt`.

**Setup (p. 1, lines 8–17):**

> "S(χ, g, f, p^m) = Σ_{x=1}^{p^m} χ(g(x)) e_{p^m}(f(x)),
> where p^m is a prime power with m ≥ 2, χ is a multiplicative character (mod p^m), e_{p^m}(·) is the additive character, e_{p^m}(x) = e^{2πi x / p^m}, and f, g are rational functions with integer coefficients."

(The sum is over x such that g, f defined on Z/(p^m) and g nonzero mod p; either f or g nonconstant.)

**Theorem 1.2 (p. 3, lines 167–179):**

> "Suppose that p is odd, f, g are rational functions over Z, not both constant, χ is any multiplicative character (mod p^m), and that m ≥ t + 2. Put λ = (5/4)^5 = 3.05…. If α is a critical point of multiplicity ν_α ≥ 1 then
> |S_α(χ, g, f, p^m)| ≤ λ_α p^{t / (ν_α + 1)} p^{m (1 − 1/(ν_α+1))},    (1.14)
> where λ_α = min(ν_α, λ)."

**Corollary 1.1 (p. 3, lines 203–212):**

> "Suppose that f, g are polynomials over Z of degrees d_1, d_2 respectively, p is an odd prime, m ≥ 1, χ is a multiplicative character (mod p^m). If m = 1 suppose that χ(g) e_p(f) is not constant on F_p (wherever it is defined), and if m ≥ 2 suppose that the sum S(χ, g, f, p^m) does not degenerate to one of smaller modulus. Then we have
> |S(χ, g, f, p^m)| ≤ 4.41 p^{m (1 − 1/(d_1+d_2))}.    (1.17)
> If p = 2 the same bound holds with constant 8.82 on the right-hand side."

**Critical structural feature for fit:** the sum is over **x = 1, …, p^m** — a *single complete sum* (running through all residues mod p^m), with **integrand factoring as multiplicative-character-of-g times additive-character-of-f**, with f, g rational functions of **bounded total degree** d_1 + d_2. The bound's quality depends on d_1 + d_2 staying bounded; if d_1 + d_2 → ∞, the bound → p^m (trivial).

---

## 2.3 Heath-Brown Theorem 1 (Heilbronn sum)

Source: `C:/tmp/crystal/heathbrown.txt`.

**Setup (p. 1, lines 8–17):**

> "Heilbronn's exponential sum is defined by
> S(a) = Σ_{n=1}^p e( a n^p / p^2 ),
> for any integer a coprime to p. It is important to note here that if n ≡ n' (mod p), then n^p ≡ n'^p (mod p^2). Thus the summand in S(a) has period p with respect to n, so that S(a) is a 'complete sum' to modulus p."

**Theorem 1 (p. 1, line 25):**

> "If p is a prime and p ∤ a then S(a) ≪ p^{11/12}, uniformly in a."

**Critical structural feature:** the phase function is **n^p mod p^2** — Fermat-quotient-related. This is a *specific algebraic structure* (the p-th power map mod p^2 lifts the identity mod p with non-trivial p-adic behavior given by q(n) = (n^{p-1}−1)/p). Stepanov's method works precisely because n^p − n ≡ 0 mod p but the lift to mod p^2 is non-trivial.

**Theorem 2 (Fermat quotient form, p. 3, lines 113–125):**

> "For any integer a coprime to p we have
> Σ_{M < n ≤ M+N, p ∤ n} e( a q(n) / p ) ≪ N^{1/2} p^{3/8},
> uniformly for M, N ≥ 1. In particular
> Σ_{n=1}^{p-1} e( a q(n) / p ) ≪ p^{7/8}, uniformly for p ∤ a."

(where q(n) = (n^{p-1} − 1)/p is the Fermat quotient).

**Hypotheses load-bearing for fit:** Heath-Brown requires the phase to be of the form n^p / p^2 or a q(n) / p — the **Fermat-quotient lift**. Syracuse's phase function (Tao eq 7.8) is θ(j, l) = {ξ · 3^{2j-2} · (2^{-l+1} mod 3^n) / 3^n} — a **2-adic exponential 2^{-l}** viewed mod 3^n. This is not a power-residue / Fermat quotient.

---

## 2.4 Kowalski 2024 (BGK expository, additive character on multiplicative subgroup of F_p^*)

Source: `C:/tmp/crystal/2401.04756v2.txt`.

**Theorem 1.1 (p. 1, lines 22–34):**

> "Bourgain, Glibichuk and Konyagin. Let γ > 0 be a real number. There exists a real number ν > 0, depending only on γ, such that for any prime number p and any subgroup H ⊂ F_p^* with |H| ≥ p^γ, we have
> | Σ_{x ∈ H} e( a x / p ) | ≪ |H| p^{-ν}
> for any a ∈ F_p^*, where the implied constant depends only on γ."

**Critical hypothesis:** **prime modulus p**, multiplicative subgroup of F_p^*, additive character. **Does not cover p^m directly.** BC 4.7 is the p^m / composite-moduli extension. Kowalski 2024 sketches the BGK proof for the prime case; the composite extension is BC (Bourgain-Chang) per Remark 4.8 of `cochrane_mccarthy.txt` (line 2461).

---

## 2.5 Shparlinski Problem 1 (multiplicative-character BGK, OPEN)

Source: `C:/tmp/crystal/shparlinski.txt`.

**Problem 1 (p. 3, lines 79–89):**

> "Obtain analogues of the results of J. Bourgain, A. A. Glibichuk and S. V. Konyagin [35] for **multiplicative character sums**
> Σ_{x_1,…,x_k ∈ X} χ(x_1 ⋯ x_k + a) and Σ_{x=1}^N χ(g^x + a)
> with very small values of N relative to p, where X ⊂ Z_p, gcd(g, p) = 1, and χ is a nonprincipal multiplicative character modulo p, see also [22, 23, 24, 25, 27, 28, 30]."

**Status:** OPEN per Shparlinski's 2016 list. Syracuse's character χ (Tao eq 7.1) is **additive**, so Problem 1 does not bite directly. But if any reduction route forced introducing a multiplicative twist, this open-problem barrier would activate.

---

## Summary table

| Theorem | Modulus | Domain | Character | Phase function | Critical hypothesis violated for Syracuse? |
|---|---|---|---|---|---|
| BC 4.7 | q (few prime factors, p_i > q^ε) | H < Z_q^* multiplicative subgroup, \|H\|=q^δ | additive e_q | linear ξ x | **YES** — q=3^n violates p_i > q^ε; also Syracuse domain is tuple-space not multiplicative subgroup |
| Cochrane 1.2 / Cor 1.1 | p^m | complete x=1..p^m | mult χ · additive e_{p^m} | f, g rational functions deg d_1, d_2 bounded | **YES** — Syracuse phase is 2-adic exponential, not rational function of bounded degree |
| HB Heilbronn Thm 1 | p | n=1..p | additive e_{p^2} | n^p / p^2 (Fermat quotient) | **YES** — Syracuse phase is 2^{-l} mod 3^n, not p-th power residue |
| Kowalski/BGK 1.1 | p prime | H < F_p^* | additive e_p | linear | **YES** — Syracuse modulus is p^n not p; domain not multiplicative subgroup of (Z/3^n)^* |
| Shparlinski Prob 1 | p | X ⊂ Z_p | **multiplicative** χ | various | irrelevant — Syracuse uses additive, not multiplicative; but route open even if needed |
