# IGUSA_I — Igusa for monomials / linear forms / "easy" classes

## Phase 0 — verbatim

For f = x_1^{a_1} · x_2^{a_2} · … · x_n^{a_n} (monomial), the Igusa zeta is the product of single-variable Mellin integrals:

Z(s, f, p) = ∏_{j=1}^n ∫_{Z_p} |x_j|^{a_j s} dx_j = ∏_j (1 − p^{-1}) / (1 − p^{-(a_j s + 1)}).

Poles at a_j s + 1 = 0, i.e., **s = −1/a_j** for each j.

For f = x_1 + x_2 + … + x_n (linear form): Z(s, f, p) = (1 − p^{-1}) / (1 − p^{-s - n}), single pole at **s = -n**.

(Reference: Denef Bourbaki §1.3, Veys Thm 3.4 example.)

## Hypothesis types

- (i) **Specific f form** (monomial, linear, or after change of variable reduces to one of these).

## Phase 1 — substrate check

R78 substrate g(u) = c + (cu·p) + (c·u(u-1)/2 · p²) + … − p² m u:
- Mod p: g ≡ c (constant).
- Reduces to a monomial / linear after change of variables? Try u → u (no change). Not a monomial: has all degrees from 0 to r.
- Try expansion at u=0: g(0) = c (nonzero), g'(0) = c·p − p²·m (3-adic v_3 = 1 for generic m). Tangent line at u=0 has slope unit·p.

Closest "easy" reduction: write g(u) = c · ∏_{k} (1 − u/α_k) for some 3-adic roots α_k. Since g mod 3 = c (constant), g has **no 3-adic roots** — the roots all live in extensions of Q_3, with valuation v_3(α_k) ≤ 0 (i.e., outside Z_3).

Result: g is a 3-adic unit on Z_3, equivalent (up to 3-adic-units) to the constant 1. Z(s, g, 3) = 1 (no poles).

## Phase 2 — conclusion shape

If we artificially convert to a monomial by lifting to higher dimension (e.g., considering f(u,v) = u^a · v^b for some auxiliary v variable), we manufacture poles at s = -1/a, -1/b. But the substrate has nothing to do with this auxiliary structure.

## Phase 3 — substrate match

Poles −1/a_j are all of the form −1/(positive integer) ∈ {-1, -1/2, -1/3, …, -1/a}. log_3(2) ≈ 0.631 is not of this form (positive and irrational).

## Disposition: NO_FIT (categorical)

Monomial / linear-form Igusa zetas all produce poles at negative rationals of form -k/m for integers k,m. log_3(2) is not of this form.
