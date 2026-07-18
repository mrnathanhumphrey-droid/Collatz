# IGUSA_J — Open literature search for c·(1+p)^u-type substrates

## Phase 0 — corpus search

Searched all 10 PDFs at C:/tmp/igusa/ for: "(1+p)^u", "principal unit", "exponential of u", "Postnikov", "log_p", "character sum on Z_p^*".

Results:
- **Denef Bourbaki § 1.4**: connects Igusa local zeta to exponential sums E(z) = ∫ Φ(x) Ψ(zf(x)) |dx| — but f is the POLYNOMIAL, the exponential is on the AMBIENT additive character side, not the integrand. No (1+p)^u inside.
- **Veys §6**: ideals of polynomials; "monomial ideal", not "exponential of u".
- **Nicaise §2-3**: Newton polyhedron + p-adic; standard polynomial f.
- **Potemans-Veys**: lectures on p-adic Igusa — references "p-adic log canonical threshold", but as a property of the polynomial f, not as an integrand component.

No paper in corpus treats f = c·(1+p)^u − target_value (or related) as an Igusa-zeta integrand. The Igusa-school object is uniformly f ∈ Z_p[x_1,…,x_n] (algebraic polynomial), not f(u) involving exponentials in u.

**Negative result: no specific c·(1+p)^u-type Igusa theorem in corpus.**

## Phase 1 — could one EXIST?

Hypothetically, one could define a "p-adic exponential Igusa zeta":
Z_exp(s, c, p) = ∫_{Z_p} |c · (1+p)^u − 1|^s du

with the modification that the integrand is an exponential function of u (not a polynomial). Examining this:

|c·(1+p)^u − 1|_p as a function of u ∈ Z_p:
- u = 0: c·(1+p)^0 − 1 = c − 1. v_3(c−1) ≥ 0 (since c ∈ (Z/p^{r+1})*, c−1 can have any v_3 ≥ 0).
- u = 1: c·(1+p) − 1 = c + cp − 1.

For c=1: at u=0, c−1 = 0 (singular!). At u≠0: (1+p)^u − 1 = pu + p²·u(u-1)/2 + … = pu·(1 + p·(u-1)/2 + …), so v_3((1+p)^u − 1) = v_3(p) + v_3(u) = 1 + v_3(u).

So for c=1: |c(1+p)^u − 1|_3 = 3^{-(1 + v_3(u))}, and
Z_exp(s, 1, 3) = ∫_{Z_3} 3^{-(1 + v_3(u)) s} du = 3^{-s} · ∫_{Z_3} 3^{-v_3(u) s} du = 3^{-s} · (1 − 3^{-1}) / (1 − 3^{-(s+1)}).

This is the **standard** Igusa zeta of f(u) = pu (a linear monomial in u with extra factor p)! No surprise: (1+p)^u − 1 ≈ pu to leading order in 3-adic norm. Pole at s = -1.

For general c with v_3(c−1) = k₀ ≥ 0: the integrand factors as a similar computation, giving pole at s = -1 (or shifted by integer if there are constant factors).

**No new pole locations from exponential-in-u reformulation.** The Igusa-of-exponential collapses (after 3-adic linearization) to Igusa-of-linear-form, which has pole at s=-1.

## Phase 2 — does the F̂_p magnitude formula carry pole information?

F̂_p^full(ξ) has |F̂_p^full(ξ)| = p^{(r+3)/2} on a support of size p^{r-1}. Plancherel: Σ |F̂(ξ)|² = p^{r+1} · (period normalization).

Magnitude p^{(r+3)/2} is **algebraic** in p (rational exponent (r+3)/2). The Mellin transform of |F̂(ξ)|^s would have pole at s such that p^{s·(r+3)/2} balances the support cardinality p^{r-1}: specifically, ∫|F̂|^s dξ = (#support) · p^{s(r+3)/2} = p^{r-1+s(r+3)/2}, no pole as a function of s (it's a single exponential), so no Mellin pole structure here either.

**The F̂_p saturation theorem does not give Mellin pole structure compatible with log_3(2) directly.**

## Disposition: NO_FIT (corpus negative; structural reformulation gives same standard poles)

No paper in corpus treats the (1+p)^u substrate as an Igusa input. Hypothetical reformulation collapses to standard linear-monomial Igusa with pole at s = -1.

## Note on ε_7 / branch-cut consistency

Even if the s = log_3(2) location were a pole, it would be **irrational**, hence not an Igusa pole real-part for any polynomial f ∈ Q[x] (excluded by Igusa rationality and by Monodromy/b-function arguments — see IGUSA_A, IGUSA_D, IGUSA_H).

The fact that R77.6's z=2 branch-cut is **REFUTED** by PADE_NUMERICAL_DISPOSITION (Hadamard radius 1.57 at n=13, NOT 2) further suggests that log_3(2) is **not even the correct numerical anchor** anymore. The actual leading singularity is at z ≈ 1.5..1.7 → s = log_3(z) ≈ 0.37..0.48 (still positive, still irrational, still inaccessible to Igusa).
