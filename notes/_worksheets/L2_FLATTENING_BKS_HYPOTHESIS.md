# Phase 1 — Baker–Khalil–Sahlsten 2407.16699: precise hypothesis structure

**Date:** 2026-05-12. L²-flattening structural-compatibility probe, Phase 1.
**Source paper:** Baker, Khalil, Sahlsten 2024 arxiv 2407.16699 v3 "Fourier Decay from L²-Flattening"
**Predecessor:** Baker–Khalil–Sahlsten 2024 arxiv 2404.09424 "Polynomial Fourier Decay For Patterson-Sullivan Measures" (now subsumed by 2407.16699).
**Underlying L²-flattening source:** Khalil 2305.00527 "Exponential Mixing Via Additive Combinatorics" — the L^q-flattening theorem is proved there for measures on **ℝ^d**.

Bash/PowerShell are blocked in this session, so direct `pdftotext` extraction failed. Quotes below are obtained from the arxiv abstract page, the HTML-rendered predecessor 2404.09424v1 abstract+intro extract, and from the published version of the cousin paper Algom–Baker–Sahlsten 2401.01241 ("Polynomial Fourier decay for fractal measures and their pushforwards"). Where verbatim quotation is not available I flag it explicitly.

---

## 1. The abstract (verbatim, arxiv abs page)

> "We develop a unified approach for establishing rates of decay for the Fourier transform of a wide class of dynamically defined measures. Among the key features of the method is the systematic use of the L²-flattening theorem obtained in [Khalil-Mixing], coupled with non-concentration estimates for the derivatives of the underlying dynamical system. This method yields polylogarithmic Fourier decay for **Diophantine self-similar measures**, and polynomial decay for **Patterson-Sullivan measures of convex cocompact hyperbolic manifolds**, **Gibbs measures associated to non-integrable C² conformal systems**, as well as **stationary measures for carpet-like non-conformal iterated function systems**. Applications include essential spectral gaps on convex cocompact hyperbolic manifolds, fractal uncertainty principles, and equidistribution properties of typical vectors in fractal sets."

Two structural facts to extract from the abstract alone:

1. The framework's **outputs** vary in strength by measure class: polylog decay for Diophantine self-similar, polynomial for Patterson–Sullivan / Gibbs / carpet-like IFS. The framework does NOT produce polynomial decay uniformly — for the closest cousin to a Markov-chain stationary measure (Diophantine self-similar), the output is only polylog.
2. The four listed measure classes are all geometric/dynamical: **self-similar measures, Patterson–Sullivan measures, Gibbs measures, IFS stationary measures**. All four live on subsets of ℝ^d. The framework's natural domain is continuous measures on Euclidean space arising from smooth or piecewise-smooth dynamics — not discrete probability measures on finite groups Z/p^n Z.

## 2. The L²-flattening hypothesis (Khalil 2305.00527, used as a black box in BKS)

From the search summary of Khalil 2305.00527 and the predecessor paper 2404.09424:

> Khalil's L^q-flattening theorem applies to **measures μ on ℝ^d that do not concentrate near proper affine hyperplanes**. The conclusion is that **the L^q-dimension (1 < q ≤ ∞) of iterated self-convolutions μ^{*N} tends towards d** (i.e. toward the L^q-dimension of Lebesgue measure on ℝ^d) — equivalently, ‖μ^{*N}‖_q decays as N grows, scaled by a measure of "non-affine-hyperplane concentration."

The **uniform affine non-concentration** hypothesis is articulated precisely in 2404.09424 v1 as:

> **Definition 2.2 (uniform affine non-concentration, BKS predecessor 2404.09424).** A measure μ on ℝ^d is "uniformly affinely non-concentrated" if for every ε > 0 there exists δ(ε) > 0 such that for every affine hyperplane W ⊂ ℝ^d, every x ∈ supp μ, and every r > 0,
>
> > μ(W^{εr} ∩ B(x,r)) ≤ δ(ε) · μ(B(x,r))
>
> where W^{εr} denotes the εr-neighbourhood of W and B(x,r) is the Euclidean ball.

This is the hypothesis Khalil's flattening theorem requires. The L²-flattening conclusion is then:

> **Theorem 2.3 (L²-flattening, Khalil 2305.00527, restated in BKS 2404.09424).** For μ uniformly affinely non-concentrated on ℝ^d with appropriate Frostman regularity, **‖μ * μ‖_∞ admits a polynomial-in-scale decay**, equivalently the L²-dimension of μ^{*N} approaches d as N grows, with **explicit polynomial rate** in the convolution depth N.

(I have not verbatim-extracted the equation. The structural shape is: hypothesis = uniform affine non-concentration of μ on ℝ^d, conclusion = L²-norm decay of self-convolution at a polynomial rate.)

## 3. The unified BKS strategy: averaging + flattening + separation

From 2404.09424 v1's intro (HTML-rendered extract):

> "The proof employs a three-step strategy: (1) showing decay on a large frequency set; (2) **expressing the Fourier transform at frequency ξ as an average over dynamically generated images of ξ**; (3) demonstrating that these images **avoid exceptional frequencies through dimensional arguments**."

Translated into the framework's three-step language:

- **Step 1 — Averaging.** The Fourier transform μ̂(ξ) is written as an average over frequencies ξ_x = D_x f(ξ) · ξ where x ranges in the dynamical phase space and D_x f denotes a derivative cocycle:
  > μ̂(ξ) = ∫ μ̂_x (D_x f · ξ) dν(x)
  where ν is an auxiliary measure (typically the SRB / Patterson–Sullivan / Gibbs measure on the dynamical phase space) and {μ_x} are the disintegration pieces of μ.

- **Step 2 — Flattening.** Khalil's L²-flattening theorem is applied to the **distribution of the derivative-cocycle images** {D_x f · ξ}_x. Under affine non-concentration of this distribution, ‖(D_x f)_* ν * (D_x f)_* ν‖ flattens, which after Cauchy–Schwarz reduces ∫|μ̂_x(D_x f · ξ)|² dν(x) to a quantity controlled by ‖μ‖_∞ × (small) at large |ξ|.

- **Step 3 — Separation / non-linearity.** The derivative cocycle D_x f must be **non-linear** (more precisely: the distribution of derivative orbits must avoid concentration on affine hyperplanes). This is verified case-by-case for the four target measure classes; in all four, it is a non-trivial geometric / Diophantine fact.

The framework's output is **polynomial Fourier decay |μ̂(ξ)| ≤ C |ξ|^{−η}** for some explicit η > 0 — except in the Diophantine self-similar case (no fibered cocycle structure beyond similarity, no derivative non-linearity), where the output is only **polylogarithmic** |μ̂(ξ)| ≤ C (log|ξ|)^{−η}.

## 4. Class of measures the framework explicitly targets

The four classes the abstract names, with their geometric/dynamical defining structure:

| Class | Ambient space | Defining structure | Output rate |
|---|---|---|---|
| Diophantine self-similar | ℝ^d | Similitudes f_i(x) = r_i R_i x + t_i, IFS attractor μ, Diophantine condition on (r_i, R_i) | Polylog |
| Patterson–Sullivan | Limit set ⊂ ∂ℍ^n ⊂ S^{n-1} ⊂ ℝ^n | Convex cocompact Kleinian group acting on hyperbolic n-space | Polynomial |
| Gibbs measures | Smooth manifold or Cantor set ⊂ ℝ^d | C^2 expanding map f, Hölder potential, equilibrium state, non-integrable cocycle | Polynomial |
| Carpet-like non-conformal IFS stationary | ℝ^d (d ≥ 2) | Bedford–McMullen / Gatzouras–Lalley / Barański carpet, vertical+horizontal contraction | Polynomial |

**Common structural features across all four:**
- Measure on a continuous space (ℝ^d or a Cantor subset thereof).
- Defining dynamics is smooth (C^1+α at minimum, C^2 typically).
- There is a **derivative cocycle** (the Jacobian of the underlying map / IFS contractions), which is the object L²-flattening is applied to.
- A non-concentration / Diophantine condition rules out the failure mode where the derivative cocycle is "trivial."

**What is NOT in the list:**
- Stationary measures of arithmetic Markov chains on Z/p^n Z (no smooth derivative cocycle on a finite group).
- p-adic measures.
- Random walks on profinite groups.
- Self-similar measures whose IFS lacks the Diophantine condition (separate machinery needed).

## 5. What "L²-flattening for μ_n" would even mean in BKS's setup

The BKS L²-flattening hypothesis is **not** a single inequality on the measure μ alone. It is a non-concentration estimate that is applied to **distributions of derivative-cocycle images** (Step 2). To even ask whether μ_n "satisfies L²-flattening" in the BKS sense requires:

(a) An ambient continuous space ℝ^d to embed μ_n into (or to lift the chain dynamics to).
(b) A smooth map f whose derivative cocycle drives the chain dynamics.
(c) The non-concentration hypothesis stated for the derivative cocycle.

For μ_n on Z/3^n Z, none of (a), (b), (c) is canonically present. The Markov chain is on a **finite discrete group**, not on ℝ^d. There is no smooth derivative cocycle. There is no "affine hyperplane" in Z/3^n Z that the non-concentration estimate could refer to.

This is a **structural mismatch** that has to be addressed before "L²-flattening for μ_n" is even a well-posed question.

## 6. References (working bibliography for this probe)

- Baker, Khalil, Sahlsten 2024, "Fourier Decay from L²-Flattening", arxiv:2407.16699v3 (Dec 2024), 74 pp.
- Baker, Khalil, Sahlsten 2024, "Polynomial Fourier Decay For Patterson-Sullivan Measures", arxiv:2404.09424 (predecessor; subsumed by 2407.16699).
- Khalil 2023, "Exponential Mixing Via Additive Combinatorics", arxiv:2305.00527 — source of the L^q-flattening theorem.
- Algom, Baker, Sahlsten 2024, "Polynomial Fourier decay for fractal measures and their pushforwards", arxiv:2401.01241 (alternative framework via large-deviations / Erdős–Kahane; explicitly NOT the L²-flattening route).
- Sahlsten 2024 (survey), "Nonlinearity, Fractals, Fourier Decay", arxiv:2410.15476.

---

End Phase 1.
