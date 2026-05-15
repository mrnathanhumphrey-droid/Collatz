# IGUSA_H — Bernstein-Sato b-function

## Phase 0 — verbatim

**Theorem 2.2 (Veys).** Let f ∈ K[x_1,…,x_n]\K. There exists a nonzero operator P(s) ∈ D_n[s] and a nonzero polynomial b(s) ∈ K[s] satisfying P(s) · f^{s+1} = b(s) · f^s. The unique monic polynomial of smallest degree satisfying this is called the Bernstein-Sato polynomial b_f(s).

**Theorem 2.3.** Z(f, φ; s) has poles of the form λ − k with λ a root of b_f and k ∈ Z_{≥0}.

**Theorem 2.6 (Kashiwara-Malgrange).** If s_0 is a root of b_f, then e^{2πi s_0} is a monodromy eigenvalue at some point of {f=0}, and conversely.

## Hypothesis types

- (i) f ∈ K[x_1,…,x_n]\K, nonzero non-constant polynomial.
- (ii) **Roots of b_f are rational and negative.** (Kashiwara theorem — Bourbaki §2.3 / Veys "the roots of b_f(s) are negative rational numbers".)

## Phase 1 — substrate check

R78 substrate g(u) ∈ Z[u] (or Q[u]). g(0) = c ≠ 0 (g is a unit) — but for the b-function we don't need f(0)=0.

For substrate (1) g(u) is **constant mod 3** but not constant as a polynomial in u over Q. So b_g(s) is well-defined as the b-function of the polynomial g ∈ Q[u].

For univariate g of degree r, b_g(s) is computable. If g has no multiple roots (generic), b_g(s) = s + 1 (the b-function of a smooth-as-a-polynomial univariate is s + 1, since g_u = g'(u) ≠ 0 generically).

## Phase 2 — conclusion shape

Roots of b_g(s) = s + 1 are at s = -1. Igusa poles ⊂ {s = -1 - k : k ∈ Z_{≥0}}.

**Real parts of Igusa poles are negative rational, by Kashiwara theorem.** log_3(2), log_3(45/43), log_3(1/0.984) are all **positive AND irrational** — doubly excluded.

## Phase 3 — substrate match

Bernstein-Sato confirms the structural barrier: **all roots of b_f are negative rational**. log_3(2) ∉ negative rationals. **CATEGORICAL BAR.**

## Disposition: NO_FIT (categorical, irrational positive target inaccessible to b-function roots)

This is the third independent confirmation (after Igusa rationality, after Monodromy Conjecture) that **log_3(2) cannot be an Igusa pole real-part for any polynomial f ∈ Q[x_1,…,x_n]**.
