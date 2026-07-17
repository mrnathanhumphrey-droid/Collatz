# TAUBERIAN_RESCOPE_C_HYPOTHESES (Korevaar 2002 — Wiener-Ikehara and Newman-style)

**Source PDF:** `C:/Users/Nate/OneDrive/Documents/tauberian/pdfs/Korevaar_2002_Century_Complex_Tauberian.pdf`
**Extracted text:** `C:/Collatz/tauberian_extract/C.txt`
**Mode:** E — verbatim from PDF, no inheritance.

---

## Theorem 4.2 (Wiener-Ikehara) — verbatim (PDF lines 442-457)

> **Theorem 4.2.** Let S(t) vanish for t < 0, be nonnegative and nondecreasing for t ≥ 0 and such that the Laplace-Stieltjes transform
>   f(z) = LdS(z) = ∫_{0−}^∞ e^{−zt} dS(t) = z ∫_0^∞ S(t) e^{−zt} dt,    z = x + iy    (4.4)
> exists for Re z = x > 1. Suppose that for some constant A, the analytic function
>   g(z) = f(z) − A/(z − 1),   x > 1    (4.5)
> has a boundary function g(1 + iy) in the following sense. For x ↘ 1, the function g_x(iy) = g(x + iy) converges to g(1 + iy) uniformly or in L¹ on every finite interval −λ < y < λ. Then
>   e^{−t} S(t) → A   as t → ∞.    (4.6)

---

## Theorem 4.1 (Landau-style precursor, Dirichlet/Mellin form) — verbatim (PDF lines 424-437)

> **Theorem 4.1.** Let s(v) vanish for v < 1, be nonnegative and nondecreasing on [1, ∞) and such that the Mellin-Stieltjes transform
>   f(z) = ∫_{1−}^∞ v^{−z} ds(v) = z ∫_1^∞ s(v) v^{−z−1} dv,    z = x + iy    (4.1)
> exists for Re z = x > 1. Suppose that for some constant A, the analytic function
>   g(z) = f(z) − A/(z − 1),   x > 1    (4.2)
> has a continuous extension to the closed half-plane x ≥ 1. Then s(u)/u → A as u → ∞.

(Remark: it is sufficient to have s(v) + Cv nonnegative and nondecreasing for some constant C.)

---

## Theorem 6.1 (Newman / Korevaar's "analytic theorem") — verbatim (PDF lines 1099+)

> **Theorem 6.1.** Let f be given for Re z > 1 by a Dirichlet series Σ a_n n^{-z}, with bounded coefficients (|a_n| ≤ C). Suppose that the analytic function g(z) = f(z) − A/(z − 1) (for some constant A) has a holomorphic extension to Re z ≥ 1. Then the series Σ a_n / n^z is convergent at z = 1, with sum f(1) − A.

(Section heading at PDF line 1183-1199 makes clear: Theorem 6.1's hypothesis "|a_n| ≤ C" can be relaxed to "a_n ≥ −C provided s(v) = O(v)" — see Theorem 8.1.)

---

## Theorem 8.1 (one-sided / Laplace form) — verbatim (PDF lines 1438-1453)

> **Theorem 8.1.** Let f(z) be given for Re z > 1 by an absolutely convergent Dirichlet series f(z) = Σ a_n n^{−z} with a_n ≥ −C for all n, where C is a constant. Suppose that for some constant A, the function g(z) = f(z) − A/(z − 1) has a holomorphic extension to Re z ≥ 1. Then the series Σ a_n / n converges with sum f(1) − A.

---

## Hypotheses extracted (for the Wiener-Ikehara Theorem 4.2 reading)

| # | Hypothesis | Source |
|---|---|---|
| h_1 | **S(t) = 0 for t < 0** (function vanishes on the negative real axis). | line 442 |
| h_2 | **S(t) is nonnegative for t ≥ 0**. | line 442 |
| h_3 | **S(t) is nondecreasing for t ≥ 0**. | line 442 |
| h_4 | **Laplace-Stieltjes transform f(z) exists for Re z > 1**. | line 444 |
| h_5 | **g(z) = f(z) − A/(z − 1) extends to Re z = 1** with continuous (or uniform / L¹ on finite intervals) boundary function. | lines 451-456 |

Output: **e^{−t} S(t) → A** as t → ∞ — a 1-term asymptotic (no expansion to higher order).

---

## Hypotheses extracted (for the Newman-Korevaar Theorem 6.1 reading)

| # | Hypothesis | Source |
|---|---|---|
| h_1' | **{a_n} bounded:** |a_n| ≤ C. | Theorem 6.1 |
| h_2' | **f(z) = Σ a_n n^{−z} converges for Re z > 1**. | Theorem 6.1 |
| h_3' | **g(z) = f(z) − A/(z − 1) has holomorphic extension to Re z ≥ 1**. | Theorem 6.1 |

Output: **Σ a_n / n converges to f(1) − A**. Note: this is convergence of the series at z = 1, NOT an asymptotic expansion of partial sums.

---

## Notational mapping for our use case

To use Theorem 4.2: construct a nonnegative nondecreasing S(t) from inputs (1)-(4) whose Laplace-Stieltjes transform has the prescribed analytic-continuation property.

To use Theorem 6.1 / 8.1: construct a bounded (or one-sided-bounded) sequence (a_n) whose Dirichlet series f(z) = Σ a_n n^{−z} has the prescribed analytic-continuation property.

Possible target sequences:
  T1. a_n = ε_n (signed exact rationals).
  T2. a_n = |ε_n| · 2^n (positive normalized magnitudes).
  T3. partial sums S(n) = Σ_{k ≤ n} |ε_k| · 2^k.
  T4. |μ̂_n(ξ)|² partial sums.

The c=7/45 closure would need an explicit identification of the pole residue A = 7/45 (or a rational involving 7/45) in g(z) = f(z) − A/(z − 1).

---

## End of C HYPOTHESES extraction.
