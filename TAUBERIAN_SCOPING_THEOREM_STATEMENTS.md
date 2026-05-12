# TAUBERIAN_SCOPING_THEOREM_STATEMENTS — Phase 2: candidate Tauberian theorem hypotheses

**Date:** 2026-05-12. Wilson. Phase 2 of the Tauberian scoping probe.

## Purpose

For each primary candidate Tauberian theorem in Hank's corpus (`C:/Users/Nate/OneDrive/Documents/tauberian/`), articulate the theorem's hypothesis and conclusion precisely, with verbatim quotes where load-bearing. Note specifically the *form* of singularity each theorem requires (pole, branch, log, etc.) and the *integrability/regularity* hypotheses.

---

## (a) arxiv 2507.15394 (Chevalier 2025): Square-root singularity Tauberian theorem

**File:** `arxiv_2507.15394_Tauberian_Square_Root_Singularity.pdf`. Author: Guillaume Chevalier, July 22, 2025. Standalone arXiv preprint, intended as part of a triptych on random walks on free groups.

### Hypothesis (Theorem 1.12, first-order version)

> *"Let g be a continuous function on D̄ whose restriction to D is holomorphic. Let Σ a_n z^n be the power series expansion of g in the neighbourhood of z = 0. Suppose there exists a function h that is holomorphic in a neighbourhood of the closure of the set D(1,1)^{1/2} = {ω^{1/2} : |ω − 1| < 1} in C, denoted D̄(1,1)^{1/2}, such that ∀z ∈ D̄, g(z) = h(√(1−z))."*

The hypothesis form is **g(z) = h(√(1−z))**, with h holomorphic on a neighborhood of the half-disk image. This is a structural assumption that g has a square-root branch singularity at z=1, factoring through a degree-2 ramified cover.

### Conclusion (Theorem 1.12)

> *"a_n = −h'(0)/(2√π) · 1/n^{3/2} + o(1/n^{3/2})."*

### Full asymptotic expansion (Theorem 1.14)

> *"There exist a unique sequence of constants (c_l)_{l≥1} and a constant C such that for any positive integer K: a_n = (1/n^{3/2})·(C + c_1/n + c_2/n^2 + ... + c_K/n^K + O(1/n^{K+1})). And C = −h'(0)/(2√π)."*

### Generalization to α ∈ (0,1) (Remark 1.15)

> *"If we replace the exponent 1/2 by some real number α ∈ ]0,1[, in the previous tauberian theorems, that is to say, g is of the form g(z) = h((1−z)^α), then ... a_n = Σ_{j<K, j∉ℕ} c_j/n^j + O(1/n^K). And the first constant c_{α+1} equals −h'(0)α/Γ(1−α)."*

So the **general singularity-type → coefficient asymptotic rule** the paper uses is:

If g(z) = h((1 − z/z₀)^α) for α ∈ (0,1), h holomorphic with h'(0) ≠ 0, then
**a_n ~ (−h'(0)·α/Γ(1−α)) · z₀^{-n} · n^{-α-1}** as n → ∞.

For α = 1/2: a_n ~ z₀^{-n} · n^{-3/2} times a constant.

### Generalization to pole-type via meromorphic h (Theorem 1.16)

> *"Suppose there exists h_p a meromorphic function ... only possesses a pole at 0 with multiplicity M ≥ 1, and suppose that h_p satisfies ∀z ∈ D̄, g(z) = h_p(√(1−z)). Then b_n = D/n^{3/2−M} · (1 + ...). And if M = 1 then D = −Res(h_p, 0)/√π."*

So a meromorphic h with pole of order M at 0 (i.e., h_p(w) = D_M/w^M + ... near w=0) gives:

**b_n ~ z₀^{-n} · n^{M-3/2}** (note: positive exponent of n when M ≥ 2; n^{-1/2} when M=1).

### Key hypothesis fidelity points

1. **g must be continuous on the closed disk** D̄ and holomorphic on D. This requires the *only* singularity on the boundary |z|=z₀ to be the one at z=z₀ (the square-root branch). Other boundary singularities falsify the hypothesis.

2. **h must be holomorphic on a neighborhood of D̄(1,1)^{1/2}** — the "lemon-shaped" half-disk image. This is the cleanest statement; it implies the branch cut on (z₀, ∞) is the only singularity of g.

3. **h'(0) ≠ 0** for the leading constant to be nonzero. If h'(0) = 0 (so h has a higher-order zero), the leading term vanishes and one must look deeper.

4. **No non-negativity required on a_n.** This is unlike Newman-Zagier / Wiener-Ikehara. Chevalier's theorem applies to signed-coefficient power series.

5. **Ramified-cover extension (Corollary 4.1)** handles z₀^n · (root of unity)^n type oscillations via the Z/dZ-equivariance condition. This handles cases where multiple equivalent boundary singularities exist as A-orbit.

---

## (b) Flajolet-Sedgewick Analytic Combinatorics Ch. VI: Singularity analysis

**File:** `Flajolet_Sedgewick_Analytic_Combinatorics.pdf`. The 12MB textbook is the standard reference; I could not page-read it in this session (PDF too large for the harness's PDF read tool), but Chevalier 2025 §Introduction confirms its content and explicitly cites Theorem VII.6 in subsection VII.6.3 as the closest FS analog.

### Core FS singularity analysis result (the "standard table")

For a function f(z) = Σ a_n z^n with a unique dominant singularity at z = ρ on its circle of convergence, where f(z) has a *singular expansion*:

f(z) ~ Σ_k c_k · (1 − z/ρ)^{α_k}, as z → ρ in a Δ-domain (sectorial domain avoiding the cut along (ρ, ∞)),

then by FS Theorem VI.1 (Transfer theorem):

**[z^n] (1 − z/ρ)^{−α} ~ (1/Γ(α)) · ρ^{-n} · n^{α−1}**, for α ∉ {0, −1, −2, ...}.

The standard FS "singularity type → coefficient asymptotic" table:

| Singularity type at z=ρ | (1−z/ρ)^β | Coefficient asymptotic [z^n] |
|---|---|---|
| Simple pole | (1−z/ρ)^{−1} | ρ^{-n} · 1 (constant, leading) |
| Pole of order m | (1−z/ρ)^{−m} | ρ^{-n} · n^{m-1}/Γ(m) |
| Square-root branch | (1−z/ρ)^{1/2} | ρ^{-n} · (−1/(2√π)) · n^{-3/2} |
| (1−z/ρ)^{−1/2} branch | (1−z/ρ)^{−1/2} | ρ^{-n} · (1/√π) · n^{-1/2} |
| Generic α-branch | (1−z/ρ)^{α} | ρ^{-n} · (−1/Γ(−α)) · n^{−α-1}, α ∉ ℕ |
| Logarithmic | log(1/(1−z/ρ)) | ρ^{-n} · (1/n) |
| Power × Log | (1−z/ρ)^{α}·log(1/(1−z/ρ))^k | ρ^{-n} · n^{−α−1} · (log n)^k / Γ(−α) |

### Hypothesis fidelity for FS singularity analysis

1. **Δ-domain analyticity:** f must be analyzable in a sectorial domain Δ(ρ, η, φ) ⊂ {z : |z| < ρ + η, |arg(z−ρ)| > φ} for some η > 0, φ < π/2. This is **weaker** than Chevalier's hypothesis (which requires g continuous on closed disk).

2. **Singular expansion at z=ρ:** f(z) = σ(z) + R(z), where σ is the *singular part* (matching the table entries) and R is "smaller" in the Δ-domain (typically O((1−z/ρ)^A) for some A larger than the leading exponent).

3. **No non-negativity required.**

4. **Multiple dominant singularities** on |z|=ρ: handled by summing the contributions (each gives a ρ^{-n}·ζ^n term where ζ is a root of unity for symmetric multiple singularities).

FS Ch. VI is the **modern operational toolkit** for extracting a_n from local singularity behavior. Chevalier's paper (a) uses a different contour-integral proof but is equivalent in conclusion to the FS α=1/2 entry plus the pole-corrected version (M ≥ 1).

---

## (c) Newman 1980 + Zagier 1997: Pole-based Tauberian (PNT-style)

**Files:** `Newman_1980_Simple_Analytic_Proof_PNT.pdf`, `Zagier_Newmans_Short_Proof_PNT.pdf`.

### Newman's convergence theorem (verbatim from Newman 1980 p.693)

> *"THEOREM. Suppose |a_n| ≤ 1 and form the series Σ a_n n^{-z} which clearly converges to an analytic function F(z) for Re z > 1. If, in fact, F(z) is analytic throughout Re z ≥ 1, then Σ a_n n^{-z} converges throughout Re z ≥ 1."*

This is a **Dirichlet series** statement, not a power series statement. Hypothesis:
- a_n bounded (|a_n| ≤ 1).
- F(z) = Σ a_n n^{-z} analytic for Re z > 1 (convergent there).
- F(z) extends analytically (no singularity) to Re z ≥ 1.

Conclusion: Σ a_n n^{-z} converges throughout Re z ≥ 1.

### Application to PNT (Newman §"Second Proof")

Newman applies this to f(z) = Σ_p (log p)/(p^z − 1) = Σ a_n/n^z where a_n = Σ_{p<n}(log p/p) − log n − c. By showing F(z) = f(z) + ζ'(z) − cζ(z) extends analytically across Re z = 1 (the only obstacle being the double pole of f at z=1, which is exactly compensated by ζ'(z) − cζ(z)), the convergence theorem yields PNT.

### Hypothesis fidelity for E(z) probe

**Newman-Zagier is NOT applicable to E(z) = Σ ε_n z^n directly.** Reasons:

1. Newman-Zagier is a **Dirichlet series** theorem, indexed by 1/n^z. The setup is fundamentally different from a power series Σ a_n z^n.

2. The hypothesis requires the function to **extend analytically across the natural boundary line**, with poles isolated and removable by adding analytic auxiliary series. R77.6 finds E(z) has a **branch cut** at z=2 — not a removable singularity. The branch cut is *not* a pole and cannot be removed by adding an analytic correction term.

3. The Dirichlet series analog of E(z) would be A(s) = Σ ε_n / n^s, which is a different object than E(z) entirely. Computing A(s)'s analytic structure would be a separate substantial probe.

**Disposition:** Newman-Zagier matches when the singularity is a pole (and one can remove it by analytic continuation). It does NOT match R77.6's branch-cut detection. **Excluded from the candidate set** for E(z)'s singularity.

---

## (d) arxiv 2504.16233 (Pierce-Turnage-Butterbaugh-Zaman 2025): Guide to Tauberian theorems

**File:** `arxiv_2504.16233_Guide_Tauberian_Arithmetic_Apps.pdf`. The "best single modern entry point" per Hank.

### Theorem A (Hypothesis A, weakest version, simple pole or pole of order m)

Hypothesis: A(s) = Σ a_n / λ_n^s (general Dirichlet series, **non-negative coefficients**) converges for Re(s) > α, and in this region A(s) = g(s)/(s − α)^m + h(s) where g is a polynomial of degree at most m − 1 with g(α) ≠ 0 and h is holomorphic in Re(s) > α and continuous in Re(s) ≥ α.

Conclusion: Σ_{λ_n ≤ x} a_n ~ c·x^α·(log x)^{m−1} as x → ∞.

### Theorem B (Hypothesis B, with power-saving remainder)

Adds growth condition |A(s)| ≤ C(1 + |Im(s)|)^κ along Re(s) = α − δ for some δ > 0, κ ≥ 0.

Conclusion: Σ_{λ_n ≤ x} a_n = Res_{s=α}[A(s)x^s/s] + O(x^{α − δ/(κ+1)}·(log x)^{m−1}).

### Pole of real order γ ∉ ℤ (the guide §14.1.1, Delange's generalization)

Hypothesis (verbatim from §14.1.1): *"Suppose that f(s) = (s−a)^{−γ} g(s) + h(s) for Re(s) > α, for a fixed real number γ with γ ∉ {0, −1, −2, ...}, and functions g, h that are holomorphic for Re(s) ≥ α, with g(α) ≠ 0. Then a(t) ~ g(α)/Γ(γ) · e^{αt} · t^{γ-1}."*

For the Dirichlet series setting: Σ_{n≤x} a_n ~ g(α)/(α·Γ(γ)) · x^α · (log x)^{γ−1}.

This is Delange's theorem (Delange 1954, Thm. III).

### Singularities of logarithmic type (the guide §14.1.3)

Hypothesis (verbatim): *"f(s) = 1/(s − α)^γ · Σ_{j=0}^k g_j(s) · (log(1/(s−α)))^j + h(s)"* where g_j, h are holomorphic on Re(s) ≥ α, g_k(α) ≠ 0.

For γ ∉ {0, −1, −2, ...}: a(t) ~ g_k(α)/Γ(γ) · e^{αt} · t^{γ−1} · (log t)^k.

For γ = −m (integer ≥ 0): different leading-order form involving (log t)^{k-1}.

### Key hypothesis fidelity points for E(z) probe

The guide is **Dirichlet series**-centric (Σ a_n / λ_n^s), not power-series. **But** the underlying machinery (Wiener-Ikehara, Delange, Ikehara-Delange) is the canonical operator-free way to extract coefficient asymptotics from analytic information about the generating function.

For a power series Σ a_n z^n, the analog is FS Ch. VI singularity analysis or Chevalier's contour-integral approach. The guide §1.4 explicitly notes that for power-series problems, FS Ch. VI is the standard tool.

**The guide is NOT the matching framework** for E(z)'s power-series setting. It is useful for **understanding the general structure** (singularity types → asymptotic shape) but not as a directly-applicable theorem.

---

## Secondary candidates briefly inspected (PDFs available in corpus)

### (e) arxiv 1910.03036 (Ingham-Euler-Maclaurin)

"Ingham's Tauberian theorem applied to generating functions with exponential growth near radius of convergence — modular forms / partition-type counting." Hank's note. The setup is for series like p(n) (partition function) where log a_n ~ c√n type growth, with exponential singularity at radius. **Not applicable** to E(z) which has polynomial-rate (1/2)^n decay.

### (f) arxiv 2511.15928 (Selberg-Delange / L-functions)

Extends Tenenbaum's ζ(s)^α · G(s) setup to products of L-functions. **Dirichlet series**, requires Euler product structure or close analog. **Not applicable** to E(z) (no Euler product structure for Collatz ε_n).

### (g) arxiv 2508.20814 (Explicit Tauberian with averaged inputs)

Removes pointwise-bound input requirement (Landau's finite-differencing). Counts Abelian number fields. **Dirichlet series**, arithmetic-structure-specific. **Not applicable** to E(z).

### (h) Karamata regular variation (foundational, in Haggstrom + Borwein)

Karamata's theorem: if f(z) = Σ a_n z^n with a_n ≥ 0, and f(z) ~ c/(1−z)^β as z → 1^-, then Σ_{n≤N} a_n ~ c·N^β/Γ(β+1). Requires **non-negativity** of a_n. R76 §11 gives ε_n alternating signs (positive at n=2, negative thereafter), so Karamata's non-negativity hypothesis **fails**. Karamata-Hardy-Littlewood family **not directly applicable**.

---

## Summary Phase 2 table

| Theorem | Singularity type covered | Power series or Dirichlet? | Non-negativity required? | Output |
|---|---|---|---|---|
| Chevalier 2507.15394 Thm 1.14 | (1−z/z₀)^α for α ∈ (0,1) including α=1/2 | Power series | No | a_n ~ C·z₀^{-n}·n^{-α-1} |
| Chevalier Thm 1.16 | meromorphic h, pole of order M at 0 | Power series | No | b_n ~ D·z₀^{-n}·n^{M-3/2} |
| FS Ch. VI | Full table: poles, branches, logs, products | Power series | No | Singularity-type → asymptotic from table |
| Newman-Zagier 1980 | Pole only (or removable boundary singularity) | Dirichlet series | No (|a_n|≤1) | Convergence on closed half-plane |
| PTBZ Thm A | Pole of integer order m at α | General Dirichlet | YES | Σ_{n≤x} a_n ~ c·x^α(log x)^{m-1} |
| PTBZ Thm B | Pole of integer order m + growth condition | General Dirichlet | YES | Σ_{n≤x} a_n = Res + power-saving error |
| Delange 1954 (PTBZ §14.1.1) | Pole of real order γ | Dirichlet | YES | Σ_{n≤x} a_n ~ c·x^α(log x)^{γ-1} |
| Delange §14.1.3 | (s−α)^{-γ}·(log(1/(s−α)))^k | Dirichlet | YES | (log t)^k factor |
| Karamata | (1−z)^{-β}, a_n ≥ 0 | Power series | YES | Partial sum asymptotic |
| Ingham-Euler-Maclaurin | Exponential singularity at radius | Generating function | ~ Yes | log a_n ~ √n type |

**Candidate theorems applicable to R77.6's setup (E(z) = power series, branch-cut at z=2, signed coefficients):**

1. **Chevalier 2025 (arxiv 2507.15394)** — directly covers (1−z/z₀)^α singularity with α ∈ (0,1), including α=1/2.
2. **Flajolet-Sedgewick Ch. VI** — broader singularity-analysis framework, covering all branch types and log factors via the standard table.

Phase 3 will match these against R77.6's empirical detection and the cached ε_n values.
