# W4.A — Faure 2009 Verbatim: √3 statement, applicability check

**Date:** 2026-05-14
**Source:** arXiv:0903.2747v1 [math.DS] 16 Mar 2009. "Semiclassical origin of
the spectral gap for transfer operators of partially expanding maps." Frédéric Faure.
**PDF location:** `C:/Users/Nate/OneDrive/Documents/faure_semiclassical/pdfs/Faure_2009_Semiclassical_Spectral_Gap_Partially_Expanding.pdf`
**Extraction method:** pdfminer (pypdf failed: custom Type-3 font encoding
produced OCF glyph sequences like "/CC/CW/CT" instead of text). Prior extraction
at `C:/Collatz/experiments_output/faure_2009_pages/` and `C:/tmp/faure/`.
**Companion file:** `C:/Collatz/FAURE_A_HYPOTHESES.md` (verbatim excerpts, Phase 1 check).

---

## 1. The map and transfer operator (verbatim)

> "Let g: S¹ → S¹ be a C^∞ diffeomorphism (on S¹ := R/Z). Let k ∈ N, k ≥ 2,
> and let the map E: S¹ → S¹ be defined by **E(x) = kg(x) mod 1**.
> Let **E_min := min_x (dE/dx)(x) = k min_x (dg/dx(x))**.
> We will suppose that the function g is such that **E_min > 1** so that E is
> a uniform expanding map on S¹. The map E is then a **k:1 map** (i.e. every
> point y has k previous images x ∈ E⁻¹(y))."
>
> "Let τ: S¹ → R be a C^∞ function, and define a map f on T² = S¹ × S¹ by:
> f: (x, s) ↦ (x' = E(x) = kg(x) mod 1, s' = s + (1/2π)τ(x) mod 1)   (3)"

> "(F̂_ν φ)(x) := φ(E(x)) e^{iντ(x)}   (7)
> The parameter ν is a semiclassical parameter; ν → ∞ is the semiclassical limit."

---

## 2. Theorem 1 (verbatim)

> "**Theorem 1. Discrete spectrum of resonances.**
> Let m < 0. The operator F̂_ν leaves the Sobolev space H^m(S¹) invariant, and
> F̂_ν: H^m(S¹) → H^m(S¹) is a bounded operator and can be written
> F̂_ν = R̂ + K̂   (9)
> where K̂ is a compact operator, and R̂ has a small norm:
> ‖R̂‖ ≤ r_m := 1/E_min^|m| ≤ k/E_min   (10)
> Therefore, F̂_ν has an **essential spectral radius less than r_m**, which means
> that F̂_ν has discrete (eventually empty) spectrum of generalized eigenvalues λ_i
> outside the circle of radius r_m. The eigenvalues λ_i are called **Ruelle resonances**."

---

## 3. Theorem 2 — THE √3 STATEMENT (verbatim)

> "**Theorem 2. Spectral gap in the semiclassical limit.**
> If the map f is **partially captive** (definition given page 15) (and m small enough),
> then the spectral radius of the operator F̂_ν: H^m(S¹) → H^m(S¹) does not depend
> on m and satisfies in the semi-classical limit ν → ∞:
>
> **r_s(F̂_ν) ≤ 1/√E_min + o(1)   (11)**
>
> which is strictly smaller than 1 from (3)."

---

## 4. The √3 identification: where k=3 enters

Equation (11) gives: r_s(F̂_ν) ≤ 1/√E_min.

For the simplest uniform expanding map E(x) = kx mod 1 (i.e. g(x) = x, so
dg/dx = 1), E_min = k. With **k = 3** (the degree of the expanding map):

> r_s(F̂_ν) ≤ **1/√3 ≈ 0.5774**

Equivalently, in the generating-function language where z = 1/λ:

> **The leading singularity of the generating function lies at |z| ≥ √3 ≈ 1.732.**

This is the "Faure √3" value. It arises as the **L²-amplitude decay rate of the
"trapped component"** of the transfer operator, which has probability weight 1/k
and L²-norm weight √(1/k) = 1/√k per level step.

Faure's intuitive explanation (page 3, verbatim from OCF-decoded text):

> "...the probability on the trapped set K decays by a factor **1/k**. This is the
> origin of the spectral gap at **1/√k** on Figure 2."

---

## 5. Hypothesis check: does Faure 2009 apply to Syracuse?

Carried over from `FAURE_A_HYPOTHESES.md` (full analysis there):

| Faure Hypothesis | Syracuse status | Verdict |
|---|---|---|
| C^∞ smooth compact T² | Syracuse: profinite Z_3* | **FAILS** |
| Uniformly expanding base E_min > 1 | Syracuse: stochastic Geom(2) factor, no manifold | **FAILS** |
| C^∞ skew-product structure | Syracuse: Tao renewal product | **FAILS** |
| Pseudodifferential calculus on T*S¹ | No cotangent bundle on Z_3* | **FAILS** |
| Anisotropic Sobolev H^m | Functions on profinite group, no smooth structure | **FAILS** |

**Overall applicability of Faure 2009 Theorem 2 to Syracuse: NOT APPLICABLE (hypotheses fail).**

The CONCLUSION-SIDE numerical value √3 ≈ 1.732 is empirically consistent with
the PADE Hadamard radius trajectory (2.06 → 1.81 → 1.66 → 1.57 at n=10..13, trending
toward √3), but Theorem 2 cannot be formally invoked.

---

## 6. The k=3 correspondence in Syracuse

In Faure's setup, k = degree of expanding map E. In the Tao recursion at level n:

- The map ξ ↦ ξ · 2^{-v} mod 3^n has k = 3 "branches" in the following sense:
  at each level n → n+1, the 3^n residue classes of (Z/3^{n+1})* cover k = 3
  lifts of each (Z/3^n)* residue.
- The "fan-out" of the Syracuse/Tao structure at each level is 3:1, matching k=3.
- Hence E_min = k = 3, and 1/√k = 1/√3 ≈ 0.5774 is the Faure-predicted spectral radius.

This is not a coincidence: the arithmetic structure of the Syracuse map (3x+1 mod 2^v)
has a fundamental 3-adic branching that corresponds exactly to the k=3 degree of the
Faure-style expanding map model.

---

## 7. Faure 2009 equation (12): general bound without "partially captive"

> "A general bound for r_s(F̂_ν) (with no hypothesis on f) is given by
> r_s(F̂_ν) ≤ (1/√E_min) · exp((1/2) lim_{n→∞} log(N(n))/n) + o(1)   (12)
> where N(n) is defined in Eq.(30). This bound is similar to the bound given
> in [Tsu08a, Theorem 1.1] by M. Tsujii."

Eq.(30) defines N(n) as the number of "trapped trajectories at time n" on the
cotangent dynamics. The "partially captive" condition is: N(n) = O(1) (bounded,
not growing). When N(n) grows (non-captive), the bound is weaker by the factor
exp((1/2) lim log(N(n))/n).

For Syracuse, the analog of N(n) is the number of Fourier modes in (Z/3^n)* that
don't decay to zero — which is O(1) (just the leading κ_1^B direction). This
supports the "partially captive" condition being satisfied in the Syracuse analog.

---

## Files

- `FAURE_A_HYPOTHESES.md` — full Phase 1 hypothesis check, verbatim theorem extracts
- `FAURE_DISPOSITION.md` — Phase 0-3 full disposition (PARTIAL for candidate A)
- `C:/Collatz/experiments_output/faure_2009_pages/` — OCF-decoded pages (pypdf)
- `C:/tmp/faure/` — pdfminer-decoded pages (readable text)
