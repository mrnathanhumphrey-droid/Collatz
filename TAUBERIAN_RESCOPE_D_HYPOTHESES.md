# TAUBERIAN_RESCOPE_D_HYPOTHESES (Newman 1980 / Zagier 1997 Analytic Theorem)

**Source PDFs:** `C:/Users/Nate/OneDrive/Documents/tauberian/pdfs/Newman_1980_Simple_Analytic_Proof_PNT.pdf` (image-only — body could not be extracted via pypdf), `C:/Users/Nate/OneDrive/Documents/tauberian/pdfs/Zagier_Newmans_Short_Proof_PNT.pdf`
**Extracted text:** `C:/Collatz/tauberian_extract/D_Newman.txt` (text empty — JSTOR cover page only), `C:/Collatz/tauberian_extract/D_Zagier.txt` (full text extracted).
**Mode:** E — verbatim from Zagier's text (Newman's body unavailable; Zagier's statement is identical in content).

---

## Analytic Theorem — verbatim from Zagier (PDF line 98-100)

> **Analytic Theorem.** Let f(t) (t ≥ 0) be a bounded and locally integrable function and suppose that the function g(z) = ∫_0^∞ f(t) e^{−zt} dt (Re(z) > 0) extends holomorphically to Re(z) ≥ 0. Then ∫_0^∞ f(t) dt exists (and equals g(0)).

---

## Hypotheses extracted (load-bearing list)

| # | Hypothesis | Source |
|---|---|---|
| h_1 | **f: [0, ∞) → ℝ (or ℂ) bounded:** ∃ B > 0 with |f(t)| ≤ B for all t ≥ 0. | Zagier line 98 |
| h_2 | **f locally integrable on [0, ∞)**. | Zagier line 98 |
| h_3 | **g(z) := ∫_0^∞ f(t) e^{−zt} dt** exists for Re(z) > 0. | Zagier line 99 |
| h_4 | **g extends holomorphically to Re(z) ≥ 0** (the closed right half-plane). | Zagier line 99 |

Conclusion: **∫_0^∞ f(t) dt converges and equals g(0)**.

(In Zagier's PNT application: f(t) = ϑ(e^t) e^{−t} − 1 where ϑ is Chebyshev's prime-counting function; g(z) = Φ(z + 1)/(z + 1) − 1/z where Φ(s) = Σ_p (log p)/p^s.)

---

## Notational mapping for our use case

To apply: construct f: [0, ∞) → ℝ bounded and locally integrable from inputs (1)-(4), whose exponential Laplace transform g(z) admits holomorphic extension to Re(z) ≥ 0.

Possible target functions:
  T1. f(t) = (something built from |μ̂_n(ξ)|² with n = ⌊e^t⌋).
  T2. f(t) = ε_⌊e^t⌋ · e^t / 5 (some normalized version of the ε_k sequence on a logarithmic time scale).

The Newman-Zagier theorem is the **weakest** Tauberian theorem in the corpus: it concludes only that ∫_0^∞ f(t) dt **converges** to g(0). It does NOT give an asymptotic expansion of finite partial sums or coefficients.

For c=7/45 closure: would need to identify f, g such that g(0) = 7/45 (or a rational related to 7/45 by an explicit normalization).

---

## End of D HYPOTHESES extraction.
