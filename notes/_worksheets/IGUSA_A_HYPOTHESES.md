# IGUSA_A — Igusa's rationality theorem

## Phase 0 — verbatim hypotheses (Denef Bourbaki §1.3.2)

**Theorem (Igusa [28],[30]).** Assume the notation of 1.1 and 1.3.1. Let f ∈ R[x_1,…,x_n] non-zero (with K a p-adic field, R its valuation ring), Φ a Schwartz-Bruhat function on K^n, χ a character of R*. Let (Y, h) be an embedded resolution of (f) over K. Let E_i, i ∈ T be the irreducible components of (h^{-1}(D))_red, with (N_i, ν_i) the numerical data (multiplicity of E_i in div(f∘h) and in div(h*(dx_1∧…∧dx_n))+1). Then:

(i) Z_Φ(s, χ, f) := ∫_{K^n} Φ(x) χ(ac f(x)) |f(x)|^s |dx| is a **rational function of q^{-s}**. Its poles are among the values s = -ν_i/N_i + (2π√-1 k)/(N_i log q) with k ∈ Z and i ∈ T such that the order of χ divides N_i.

(ii) If C_f ∩ Supp Φ ⊂ f^{-1}(0), then Z_Φ(s,χ,f) = 0 for almost all χ.

## Hypothesis types

- (i) **Polynomial f**: any non-zero polynomial f ∈ R[x_1,…,x_n] (the valuation ring of a p-adic field K). No degree, smoothness, or non-degeneracy restriction.
- (ii) **Local ring**: R = Z_p for K = Q_p; works for any prime p.
- (iii) **Integration measure**: standard Haar |dx|, multiplied by Schwartz-Bruhat function Φ.
- (iv) **Conclusion**: Z is rational in q^{-s}; pole list **inexplicit without embedded resolution data** (the (N_i, ν_i) come from a resolution which has to be constructed for the specific f).

## Phase 1 — hypothesis × R78 substrate

**Substrate (1): g(u) = c · Σ_{k=0}^r C(u,k) · p^k − p²·m·u mod p^{r+1}, viewed as element of Z_p[u].**

- (i) Polynomial: SATISFIED (g(u) is a univariate polynomial in u over Z_p of degree r).
- (ii) Local ring: SATISFIED (Z_3).
- (iii) Measure: SATISFIED (Haar on Z_3, Φ = char of Z_3 or of u-range).
- (iv) Conclusion: theorem applies, Z(s, g, 3) is rational in 3^{-s}.

**Substrate (2): cubic phase polynomial P_a(s) = p s − C_a · L_p(1+ps) from R78 Feature 7 (Postnikov-style).**

- (i) Polynomial: SATISFIED if we restrict to a fixed a (Postnikov-truncated p-adic log gives a polynomial of fixed degree).
- (ii) Local ring: SATISFIED.
- (iii) Measure: SATISFIED.
- (iv) Conclusion: theorem applies.

## Phase 2 — conclusion-shape check

For substrate (1) — g(u):

The critical question is the candidate pole list −ν_i/N_i.

**KEY 3-ADIC OBSERVATION (R78's D=0 finding, rederived):** g(u) mod 3 = c (a non-zero unit, since c ∈ (Z/3^{r+1})*). Therefore **v_3(g(u)) = 0 for all u ∈ Z_3**.

Consequence:
> Z(s, g, 3) = ∫_{Z_3} |g(u)|_3^s du = ∫_{Z_3} 1 · du = 1.

The Igusa zeta of g is **TRIVIAL** (the constant function 1, no poles anywhere).

Equivalently in the resolution language: the singular locus C_g = {u : g(u) = 0, g'(u) = 0} ∩ Z_3 is empty (g is a unit, has no zeros in Z_3 mod 3). The embedded resolution is the identity. There are no E_i, no (N_i, ν_i), no candidate poles.

For substrate (2) — cubic phase P_a(s):

After Postnikov substitution and saddle-point manipulation, P_a(s) = p s − C_a · L_p(1+ps) where L_p is the truncated p-adic logarithm. Expanding: L_p(1+ps) = ps − (ps)²/2 + (ps)³/3 − … truncated at the relevant level. For p=3 and Cochrane-Pinner-cubic the relevant polynomial is degree 3 in s with leading coefficient 3³/3 = 9 (or similar), with the constant term being a unit.

Computing v_3(P_a(s)): the constant term in s is **0** (P_a(0) = 0). The first-order term is **p·s = 3s** (after subtracting C_a · L_p contribution, which gives at order s: −C_a · p · s, so total coefficient is p · (1 − C_a)). For generic C_a ≠ 1, the linear coefficient is a unit times p, so v_3(P_a(s)) for small s is v_3(s) + 1.

At s = u (3-adic integer), v_3(P_a(u)) = v_3(u) + 1 generically. The Igusa zeta is:
> Z(s, P_a, 3) = ∫_{Z_3} |P_a(u)|_3^s du = Σ_{k≥0} (Vol{u : v_3(u) = k}) · 3^{-(k+1)s}
>             = (1 − 3^{-1}) · 3^{-s} · Σ_{k≥0} 3^{-k} · 3^{-ks}
>             = (2/3) · 3^{-s} / (1 − 3^{-1-s})

This has a **single simple pole at 3^{-1-s} = 1, i.e., s = -1** (with shifts by 2πi k / log 3).

**The pole location is s = -1**, NOT s = log_3(2) ≈ 0.631 or any of the alternative anchors.

## Phase 3 — substrate match: is log_3(2) in the pole list?

**Substrate (1):** pole list is **EMPTY**. Conclusion: **POLYNOMIAL_FORM_WRONG** (or "the polynomial degenerates to a unit, Igusa zeta trivial, no information").

**Substrate (2):** pole list is {s = -1 + 2πi k / log 3 : k ∈ Z}. Real part is -1, **not** log_3(2) ≈ 0.631, log_3(45/43) ≈ 0.041, or log_3(1/0.984) ≈ 0.015. **POLE_LIST_MISMATCH.**

## Disposition: NO_FIT

Igusa rationality applies (vacuously for substrate 1; non-trivially with pole s=-1 for substrate 2), but the pole list does NOT include any of the target locations.

For substrate 1 the issue is **structural**: the natural univariate polynomial g(u) is a 3-adic unit (R78's D=0), so its Igusa zeta is trivial. The structural reason: 4 = 1 + 3 in Z_3 makes (1+3)^u = 4^u 3-adically just oscillate within Z_3*, never approaching 0.

For substrate 2 the issue is **pole location**: the natural pole of the cubic-Postnikov phase integral is at s = -1 (the standard "log canonical threshold" for a tangent-line-touching polynomial), not at log_3(2).

**Disposition: NO_FIT (pole list mismatch, structurally driven by D=0).**
