# IGUSA_B — Denef-Hoornaert Newton-polyhedron explicit formula

## Phase 0 — verbatim hypotheses

**Theorem 0.27 (Denef-Hoornaert) [Bories-Veys §0.4, citing Denef-Hoornaert].** Let f(x) = f(x_1,…,x_n) be a nonzero polynomial in Z_p[x_1,…,x_n] satisfying f(0) = 0. Suppose that f is **non-degenerated over F_p with respect to all the compact faces of its Newton polyhedron Γ_f**. Then the local Igusa p-adic zeta function associated to f is the meromorphic complex function:

Z_f^0 = Σ_{τ compact face of Γ_f} L_τ · S(Δ_τ),

with

L_τ(s) = ((p−1)/p)^n − N_τ/p^{n-1} · (p^s − 1)/(p^{s+1} − 1),
N_τ = #{x ∈ (F_p×)^n : f̄_τ(x) = 0},
S(Δ_τ)(s) = Σ_{k ∈ Z^n ∩ Δ_τ} p^{−σ(k) − m(k) s},

and σ(k) = k_1 + … + k_n.

**Candidate pole theorem (Thm 0.26 Denef-Hoornaert):** If s_0 is a pole of Z_{f,Φ}, then either
  (1) s_0 = −1 + 2kπi / log p for some k ∈ Z, or
  (2) s_0 = −σ(v_j)/m(v_j) + 2kπi/(m(v_j) log p) for some j and some k ∈ Z.

## Hypothesis types

- (i) **Polynomial f**: nonzero, **multivariable** in n ≥ 1 variables over Z_p, with f(0) = 0.
- (ii) **Non-degeneracy condition**: for every compact face τ of Γ_f, the zero locus f̄_τ in (F_p×)^n has no singularities.
- (iii) **Local ring**: Z_p, any prime.
- (iv) **Conclusion**: explicit formula in terms of compact faces of Γ_f and combinatorics of cones Δ_τ; pole list **explicitly readable** as s = −σ(v_j)/m(v_j) where v_j are primitive vectors perpendicular to facets.

## Phase 1 — hypothesis × R78 substrate

**Substrate (1): g(u) ∈ Z_3[u], univariate, degree r.**

- (i) Polynomial: SATISFIED but **n=1 univariate is a degenerate case**. Γ_g ⊂ R≥0 is a half-line. Compact faces: the single vertex at the lowest exponent k_min with non-zero coefficient.
- R78 D=0 finding: g(u) = c · (1 + 3u + 9·u(u-1)/2 + …) − 9mu. Lowest-degree term in u is c (the constant term), coefficient is a 3-adic unit (c ∈ (Z/3^{r+1})*).
- But the theorem requires f(0) = 0; here g(0) = c ≠ 0. **HYPOTHESIS h_FAILS at f(0)=0.**

  *Workaround:* substitute g(u) − c (shift by constant) — but then the Igusa integrand |g(u) − c|^s ≠ |g(u)|^s; we are computing a different integral.

  *Alternative:* take f(u) := u (linear monomial, n=1, f(0)=0, Newton polyhedron is the ray [1,∞)). This isn't R78's substrate; it's a placeholder.

  *True polynomial-substrate for which f(0)=0 and which carries R78 content:* g(u) − c = c·u·(p + p²·(u-1)/2 + p³·(u-1)(u-2)/6 + …) − p²·m·u. Factor out p: g(u) − c = p · [c·u + (cp/2)·u(u-1) + … − p·m·u]. So g − c = p · g̃(u) where g̃ is a polynomial in u of degree r with leading term involving 1/r!. **Newton polyhedron of g̃ in 1D is the ray [1, ∞), single facet at v = 1, m(v=1) = 1, σ(v=1) = 1.**
  
- (ii) Non-degeneracy: in 1D the only compact face is the vertex at u^1. f̄_τ = c̄ · ū = c̄u (after reduction mod p, c̄ ∈ F_p× since c is a unit). Zero locus in (F_p×)^1 is **empty** (c̄u = 0 with u ∈ F_p× has no solution since c̄ ≠ 0). So non-degeneracy is **SATISFIED VACUOUSLY**.
- (iii) Z_p with p=3: SATISFIED.
- (iv) Conclusion: pole list reads off facets.

  **Pole list (univariate):** the only primitive vector perpendicular to the single facet of Γ_{g̃} = [1,∞) is v = 1. σ(v=1) = 1, m(v=1) = 1.
  Candidate poles: s = −σ(v)/m(v) + 2kπi/(m(v) log p) = **−1 + 2kπi / log 3**.
  Plus the universal candidate s = −1 from Thm 0.26 (1).

  Conclusion: **only real-part pole is s = −1.**

**Substrate (2): cubic Postnikov phase P_a(s).**

- (i) Polynomial: SATISFIED (polynomial in s of degree 3).
- (ii) Non-degeneracy: depends on a; for generic a, the cubic P_a(s)/p has non-degenerate Newton polyhedron (single vertex at u^1, since lowest non-zero term is linear). Same as substrate (1).
- (iii) Z_3: SATISFIED.
- (iv) Pole list: same (single facet v=1, m(v)=1, σ(v)=1) → s = -1.

## Phase 2 — conclusion-shape check

For both natural univariate substrates, Denef-Hoornaert gives:

> Z(s, f, 3) has unique real-part pole at s = -1.

The explicit formula:
> Z(s, g̃, 3) = (2/3) · 3^{−s} / (1 − 3^{−1−s})

confirms a simple pole at s = -1, residue = (2/3) · 3 / log 3 (up to normalization).

## Phase 3 — substrate match: is log_3(2) in the pole list?

**NO.** The only pole at real part is s = −1.

- Real part of pole: **−1**.
- Target log_3(2) ≈ **+0.631**.
- Target log_3(45/43) ≈ **+0.041**.
- Target log_3(1/0.984) ≈ **+0.015**.

**POLE_LIST_MISMATCH.** Real parts of poles in Igusa local zeta are always **≤ −1/d** (where d is the max degree appearing in the Newton polyhedron); they cannot be positive. The log_3(2) and other targets are **POSITIVE** real numbers — categorically outside the Igusa pole locus.

This is a **fundamental categorical mismatch**: Igusa local zeta poles have **negative real part** (by the basic estimate |s_0| ≤ 0 from convergence of ∫|f|^s for Re(s) > 0). The targets log_3(2), log_3(45/43), log_3(1/0.984) are all **positive** — they cannot be Igusa poles at all.

## Disposition: NO_FIT (categorical)

The Igusa local zeta of any polynomial f ∈ Z_p[x] has all poles with Re(s) < 0 (in fact, by Igusa's theorem, ≤ −ν_min/N_max < 0). The R77.6 branch-cut at z=2 corresponds to s = log_3(2) > 0; **positive real-part Mellin variable, NOT an Igusa pole location**.

This is the **categorical barrier** for the entire Igusa probe: positive Mellin variable ↔ negative Igusa pole, off by a sign.

**Sign-flip variant:** if instead we identified the Mellin variable as z = 3^{-s} (with z=2 ↔ −s = log_3(2), i.e., s = -log_3(2) ≈ -0.631), then we'd be looking for a pole at NEGATIVE s = -0.631. This is in the Igusa range (negative). Substrate-1 gives s = -1; substrate-2 cubic Postnikov gives s = -1 — neither produces -0.631.

The closest Newton-polyhedron rational that gives Re(s_0) = -log_3(2) would need σ(v)/m(v) = log_3(2). But σ and m are integer-valued; their ratio is rational, **not** log_3(2) (irrational). No Newton polyhedron over Z gives an irrational pole real-part.

**Disposition: NO_FIT — STRUCTURAL CATEGORICAL BARRIER.**
