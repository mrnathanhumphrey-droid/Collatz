# FG Candidate G — Saloff-Coste 2004 (Random walks on finite groups, mixing time)

**PDF:** Saloff_Coste_2004_Random_Walks_Finite_Groups.pdf.
**Extracted text:** `C:/tmp/fg/saloff_coste_2004.txt`.

---

## Saloff-Coste 2004 Theorem 2.1 (VERBATIM, p. 267)

> "Let K be an irreducible Markov kernel on a finite state space X. Then K admits a unique invariant distribution π and ∀ x, y ∈ X, lim_{t→∞} H_t(x,y) = π(y). Assume further that K is aperiodic. Then the chain is ergodic, that is, ∀ x, y ∈ X, lim_{n→∞} K^n(x,y) = π(y)."

### Hypotheses (typed):

- h_G.2.1.state: X finite state space. [TYPE (i)]
- h_G.2.1.K: K irreducible Markov kernel. [TYPE (ii)]
- h_G.2.1.aperiodic: K aperiodic (for ergodicity conclusion). [TYPE (ii)]

### Conclusion:

- Existence & uniqueness of π; convergence of K^n(x, ·) to π.

This is **qualitative ergodicity**, not Fourier decay rate.

---

## Saloff-Coste 2004 Reversibility framework (eq 2.5, p. 268)

> "When (K, π) satisfies (2.5) [∀ x, y ∈ X, π(x) K(x, y) = π(y) K(y, x)], one says that K is reversible with respect to π and that π is a reversible measure for K. Equation (2.5) is also called the detailed balance condition in the statistical mechanics literature."

The rest of the survey (spectral gap λ_1, log-Sobolev, Nash inequality, mixing time bounds) **assumes reversibility (detailed balance) throughout**.

### Hypotheses for spectral-gap mixing bounds (e.g. Theorem 3.1, p. 285):

- h_G.3.1.K_finite: K on finite group G. [TYPE (i)]
- h_G.3.1.reversible: detailed balance (2.5). [TYPE (ii)]
- h_G.3.1.spec_gap: spectral gap λ_1 = 1 − β_1 of K-operator on L^2(π). [TYPE (ii)]

### Conclusion:

- Mixing time T_TV(G, p) = O((1/λ_1) log |G|), etc.

---

## Phase 1 — hypothesis × input matrix

| Hypothesis | (1)-(4) Disposition |
|---|---|
| h_G.state: X finite | SATISFIED at each n: (Z/3^n)* is finite (size 2·3^{n-1}). |
| h_G.K_irreducible | SATISFIED: Tao's chain K_n on (Z/3^n)* is irreducible (verified empirically + ergodicity stated in Tao §6). |
| h_G.aperiodic | SATISFIED: the step distribution Geom(2) has positive probability at each v ∈ N+1, and 2 generates (Z/3^n)*, so the chain is aperiodic. |
| h_G.reversible (detailed balance) | **FAILED.** The Syracuse chain step is x → x · 2^{-v} (with v ~ Geom(2)) which is **directed**: the reverse step would be x → x · 2^{+v} (mod 3^n), but the chain only does 2^{-v} → 2^{-v'} type transitions. The stationary distribution π_n is not reversible w.r.t. the Tao kernel; it is a directed Markov chain. This is the same problem as Probe R2 ("transient vs stationary"). FAILED. |
| h_G.spec_gap λ_1 well-defined | For reversible chains, λ_1 is the L^2(π) operator-norm gap. For non-reversible Syracuse, the L^2 operator K has complex spectrum; the "spectral gap" is interpreted as gap to 1 in absolute value, which exists for ergodic chains. But the Saloff-Coste machinery (Cheeger / Nash / log-Sobolev / functional inequalities) **assumes reversibility** to convert spectral gap to mixing time. FAILED at the load-bearing assumption. |

**Phase 1 disposition: NO_FIT** on reversibility hypothesis.

---

## Phase 2 — conclusion shape

Saloff-Coste's results bound **mixing time of K^n(x, ·) to π** in TV / L^2 / L^∞ distance. This is *not* Fourier decay of π. The Fourier coefficients of π are the closure target; mixing-time bounds tell you how fast the chain approaches π, not how regular π is. **CONCLUSION_SHAPE_MISMATCH.**

(Side note: there's a *vague* connection — fast L^2 mixing implies the chain "smooths" Fourier coefficients quickly — but this is the wrong direction. L^2 mixing depends on |π̂_n(ξ)| being small for the eigendecomposition of K, not on |π̂(ξ)| being small as a *property* of π.)

---

## Phase 3 — profinite extension

Reversibility is the load-bearing assumption. The Syracuse chain is fundamentally non-reversible: it implements the dynamics T(x) → x' with the 2-adic Geom(2) step, which has no time-reversal symmetry.

There exist non-reversible-Markov-chain mixing-time results (Diaconis-Saloff-Coste 1996, Wilson 2004) but these typically give *worse* bounds, and still don't deliver Fourier decay.

**Phase 3 disposition: STRUCTURALLY_BLOCKED at reversibility.** The Saloff-Coste framework is fundamentally inapplicable to non-reversible chains.

---

## Disposition G: **NO_FIT** (categorical at reversibility, with conclusion-shape mismatch).

- Phase 1: reversibility fails categorically.
- Phase 2: conclusion shape is mixing-time, not Fourier decay.
- Phase 3: structural block (non-reversible chains outside the framework).
