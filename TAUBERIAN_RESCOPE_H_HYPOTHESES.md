# TAUBERIAN_RESCOPE_H_HYPOTHESES (Secondary candidates)

**Date:** 2026-05-13. Brief verbatim extraction from each H-class candidate, just enough to assess fit / nonfit. Detailed hypothesis-check delivered inline (no separate _CHECK.md files for H-class — they are all transparent NO_FIT or PARTIAL).

---

## H3: Bringmann-Jennings-Shaffer-Mahlburg 1910.03036 — Ingham + Euler-Maclaurin Tauberian

Source: `C:/Users/Nate/OneDrive/Documents/tauberian/pdfs/arxiv_1910.03036_Ingham_Euler_Maclaurin.pdf`.

**Theorem 1.1 (PDF lines 77-110):** verbatim

> "Suppose that B(q) = Σ_{n≥0} b_n q^n is a power series with **non-negative real coefficients** and radius of convergence at least one. If λ, α, β, and γ are real numbers with **γ > 0** such that
>   B(e^{-t}) ~ λ log(1/t)^α t^β e^{γ/t} as t → 0+,
>   B(e^{-z}) ≪ log(1/|z|)^α |z|^β e^{γ/|z|} as z → 0,    (1.2)
> with z = x + iy in each region of the form |y| ≤ Δx for Δ > 0, then
>   Σ_{n=0}^N b_n ~ λ γ^{β/2 − 1/4} log(N)^α / (2^{α+1} √π N^{β/2 + 1/4}) e^{2√(γN)}."

**Hypotheses:**
- b_n ≥ 0 (NON-NEGATIVE coefficients — load-bearing).
- Radius of convergence ≥ 1.
- B(e^{-t}) has explicit exponential-type singularity λ log(1/t)^α t^β e^{γ/t} with γ > 0.
- Bound on B(e^{-z}) in restricted angle.

**Fit against inputs:**
- (1) ε_k: signed (NOT non-negative). FAILS h_1 (non-negativity).
- (1) |ε_k|·2^k: non-negative; could be tried. But the conclusion's e^{2√(γN)} super-polynomial growth predicts a partition-like growing partial sum — not the bounded-or-decaying pattern of Σ |ε_k|·2^k. So even if non-negativity is fixed by |·|, γ > 0 (required for the theorem) gives super-polynomial conclusion — Syracuse ε_k is not exponential-growth.

**Disposition: NO_FIT** — wrong category of asymptotic (Ingham targets partition functions with e^{c√N} growth; Syracuse ε_k has decay).

---

## H4: Granville (?) — Selberg-Delange remarks

Source: `C:/Users/Nate/OneDrive/Documents/tauberian/pdfs/Selberg_Delange_Remarks_On_SD.pdf` (16 pp).

**Relevance:** SD-method targets Σ_{n ≤ x} f(n) for multiplicative f whose Dirichlet series has analytic continuation of the form ζ(s)^z · h(s). Syracuse's ε_k is not a multiplicative-arithmetic-function partial sum. Same h_1-style failure as F (LSD): no Dirichlet-character / multiplicative-arithmetic-function structure.

**Disposition: NO_FIT** (same categorical mismatch as F).

---

## H6: Häggström — Basics of Tauberian (118 pp, textbook-style)

Source: `C:/Users/Nate/OneDrive/Documents/tauberian/pdfs/Haggstrom_Basics_of_Tauberian.pdf`.

**Content:** Pedagogical / textbook treatment of classical Tauberian theory (Hardy-Littlewood, Karamata, Wiener, Ikehara). Contains no novel theorem that isn't already covered by C, D, E above.

**Disposition: NO_FIT (subsumed by C/D/E)** — same theorems, weaker formulations.

---

## H9: 2504.16233 — Guide to Tauberian arithmetic applications

Source: `C:/Users/Nate/OneDrive/Documents/tauberian/pdfs/arxiv_2504.16233_Guide_Tauberian_Arithmetic_Apps.pdf` (83 pp).

**Theorems A, B (PDF lines 22-30, 322-365):**

> "Hypothesis A for the Tauberian theorem requires {a_n} ⊂ ℂ such that the Dirichlet series A(s) = Σ a_n n^{-s} has analytic continuation to a region containing the line Re(s) = σ_a (the abscissa of convergence) with at most poles at s = σ_a; explicit growth conditions on |A(s)| in vertical strips and the residues at s = σ_a are required for Theorem A. Theorem A then gives Σ_{n ≤ x} a_n = (main term from residues) + O(x^{σ_a − δ})."

**Fit against inputs:** Same Mode H circular trap as C, D, E — requires *a priori* analytic continuation of the Dirichlet series, which IS the closure target.

**Disposition: BLOCKER (Mode H circular)** — subsumed by C, D, E.

---

## H10-H12: Lagarias 3x+1 references

Sources:
- `Lagarias_2111.02635_3x_plus_1_Overview.pdf` (27 pp) — survey.
- `Lagarias_3x_plus_1_Generalizations.pdf` (21 pp) — generalizations review.
- `Lagarias_Stochastic_Models_3x_plus_1.pdf` (74 pp) — stochastic models.

**Content:** Surveys / expository accounts of the 3x+1 problem. References to Tauberian-type results are at the literature-pointer level (e.g., Lagarias et al. discuss Wirsing/Friedman density results and the heuristic average ratio 3/4 per iterate, but do not state a Tauberian theorem applicable to Tao's Syracuse Fourier coefficient closure).

**Disposition: BLOCKER (pointer-only, not theorem statements)** — same disposition as G.

(NOTE: The Lagarias surveys do contain bibliographic pointers to *other* papers that might be Tauberian — Davison 1985 "Generalized 3x+1 mappings: Markov chains" cites a measure-rigidity argument, and Tao 2019 cites Crandall, Lagarias, Eliahou for 3-adic measure construction. These pointers are not theorems themselves.)

---

## H1, H2, H7, H8 — Borwein, Holland, Mandrekar, Riemenschneider

- H1 Borwein survey (17 pp): same century-of-Tauberian content as Korevaar C. **Subsumed by C; NO_FIT/BLOCKER.**
- H2 Holland Abel limit theorem (8 pp): Abel summability → ordinary convergence under Tauberian condition n a_n bounded. **Requires (n ε_n) bounded — not verifiable; same Mode H trap.**
- H7 Mandrekar on Wiener (PDF on Wiener's work): expository. **Pointer-only.**
- H8 Riemenschneider PNT simple analytic proofs: subsumed by D (Newman-Zagier). **NO_FIT/BLOCKER.**

**All H1, H2, H7, H8: NO_FIT or BLOCKER (subsumed by primary candidates).**

---

## Aggregate H disposition

No H candidate produces a SELECTED outcome. All are either:
- NO_FIT (wrong category of object, e.g., H3 requires partition-type non-negative growth),
- BLOCKER (Mode H circular, same as primary; or pointer-only).
