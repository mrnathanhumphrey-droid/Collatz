# FG Candidate A — Furstenberg 1963 (Lyapunov / non-commuting random products)

**PDF:** `furstenberg_grav/Furstenberg_1963_Noncommuting_Random_Products.pdf`
**Extracted text:** `C:/tmp/fg/furstenberg_1963.txt`
**Statements located:** Theorem 1.1 (p.9), Theorem 1.2 (p.10), Theorem 8.5 (p.48), Theorem 8.6 (p.50).

The Furstenberg-Kesten 1960 statement is NOT in corpus (paywalled); Theorems 8.5 + 8.6 are the closest structural Lyapunov statements in the 1963 paper.

---

## Theorem 1.1 (VERBATIM, p. 384–385):

> "If G is any locally compact group, p a probability measure on G, then G possesses a stationary measure for p only if p has its support in a compact subgroup of G."

### Hypotheses (typed):

- h_A.1.1.group: G is a **locally compact group**. [TYPE (i)]
- h_A.1.1.walk: p is a probability measure on G. [TYPE (ii)]
- h_A.1.1.stationary: existence of a stationary measure for p on G. [TYPE (iii)]

### Conclusion:

- C_A.1.1: support of p ⊂ compact subgroup of G. [TYPE (iv) — structural negative]

---

## Theorem 1.2 (VERBATIM, p. 386):

> "Let p be a probability distribution on the set of m × m unimodular matrices and let G be the smallest closed subgroup of SL(m,R) containing the support of p. Then p admits a stationary measure on the G-space R^m − {0} only if G is either compact or reducible."

### Hypotheses (typed):

- h_A.1.2.group: G < **SL(m,R)** (unimodular matrices). [TYPE (i)]
- h_A.1.2.walk: p probability measure on SL(m,R). [TYPE (ii)]
- h_A.1.2.stationary: stationary measure for p on the G-space R^m − {0}. [TYPE (iii)]

### Conclusion:

- C_A.1.2: G is compact OR G is reducible. [TYPE (iv) — structural]

---

## Theorem 8.5 (VERBATIM, p. 424):

> "If G is irreducible and ∫ log ||g|| dp(g) < ∞, the expression ∫∫ p_1(g,ξ) dp(g) dπ(ξ) for a stationary measure π for p is independent of the stationary measure π. Denoting the common value by α_p(p_1), we have, with probability 1, n^{-1} log ||X_n…X_1 u|| → α_p(p_1) for all nonzero vectors u ∈ R^m."

### Hypotheses (typed):

- h_A.8.5.group: G < SL(m,R), G **irreducible**. [TYPE (i)]
- h_A.8.5.walk: p probability measure on G with **∫ log||g|| dp(g) < ∞**. [TYPE (ii)]
- h_A.8.5.stationary: stationary measure π for p on P^{m-1}. [TYPE (iii)]

### Conclusion:

- C_A.8.5: Lyapunov exponent α_p(p_1) := ∫∫ p_1(g,ξ) dp dπ is well-defined (independent of π) and equals n^{-1} log ||X_n…X_1 u|| almost surely. [TYPE (iv) — Lyapunov exponent]

---

## Theorem 8.6 (VERBATIM, p. 426):

> "Let G be a noncompact subgroup of SL(m,R) such that no subgroup of G of finite index is reducible. Then α_p(p_1) > 0 and with probability 1, ||X_n…X_1 u|| grows exponentially as n → ∞ for all u ∈ R^m − {0}."

### Hypotheses (typed):

- h_A.8.6.group: G < SL(m,R), **noncompact**, no reducible finite-index subgroup. [TYPE (i)]
- h_A.8.6.walk: p with finite log-moment (inherited from 8.5). [TYPE (ii)]

### Conclusion:

- C_A.8.6: positive Lyapunov α_p > 0; exponential growth of ||X_n…X_1 u||. [TYPE (iv) — strict positivity of Lyapunov]

---

## Phase 1 — hypothesis × input matrix for Theorem 8.5+8.6 (the quantitative-rate Lyapunov statement)

| Hypothesis | (1) Tao form | (2) C1 renewal | (3) R75/76/77 | (4) eps_k k=1..8 |
|---|---|---|---|---|
| h_A.group: G < SL(m,R) | **FAILED**: Syracuse chain on (Z/3^n)*; no matrix-product representation; (Z/3^n)* is abelian (cyclic of order 2·3^{n-1}), no analog of SL(m,R) structure | FAILED | FAILED | FAILED |
| h_A.group: G irreducible | N/A (no matrix structure) | N/A | N/A | N/A |
| h_A.walk: p prob on SL(m,R) | **FAILED**: Syracuse step is 2-adic Geom(2); domain is N+1 (geometric on positive integers), maps to (Z/3^n)* via 2^{-a}; NOT a probability on a matrix group | FAILED | FAILED | FAILED |
| h_A.walk: ∫ log||g|| dp finite | N/A (no matrix norm) | N/A | N/A | N/A |
| h_A.stat: stationary π exists | SATISFIED (Tao's π_n is the chain's stationary distribution on (Z/3^n)*) | SATISFIED | SATISFIED | SATISFIED |

**Phase 1 disposition for Theorem 8.5/8.6:** NO_FIT. h_A.group fails on every input — Syracuse is on a **profinite abelian** group, Furstenberg-Kesten-Lyapunov machinery is for *matrix products* in SL(m,R) (or analog non-abelian Lie groups). The Lyapunov exponent itself measures growth of ||X_n…X_1 u|| — Syracuse has no "u" vector to multiply, and the chain step does not act linearly on a vector space.

## Phase 1 — Theorem 1.1/1.2 (existence-of-stationary, structural negative)

Theorem 1.1 is qualitative: stationary measure exists only if support compact. Syracuse satisfies this trivially (support of step is in compact subgroup of the action, but action is on (Z/3^n)* via x → x · 2^{-a} which is a contraction structure, not a group element acting on a homogeneous space).

But: Theorem 1.1 CONCLUSION is "support of p is in a compact subgroup" — which is not a Fourier-decay statement at all. **Conclusion-shape mismatch — NO_FIT for closure target.**

## Phase 2 — conclusion shape

For 8.5/8.6: even if hypotheses were satisfied, the conclusion is a **Lyapunov exponent** (rate of exponential growth of vector norm under random matrix product). This does NOT bound |mu_n(xi)| in any direct sense. Lyapunov exponents give *asymptotic norms*, not Fourier decay. CONCLUSION_SHAPE_MISMATCH.

## Phase 3 — profinite extension

The proofs of 8.5/8.6 rely on:
- Submultiplicativity of matrix norms ||g_n g_{n-1}|| ≤ ||g_n|| · ||g_{n-1}||.
- Cocycle structure σ(gh, ξ) = σ(g, hξ) σ(h, ξ) on the projective space P^{m-1}.
- Compactness of P^{m-1} (a smooth manifold).

In the profinite setting (Z/3^n)*: no projective space, no submultiplicative norm of multiplication, no cocycle structure of matrix type. The whole Lyapunov framework is built on the **non-commutative matrix structure** of SL(m,R) — the Syracuse multiplicative group (Z/3^n)* is abelian, so any "Lyapunov exponent" is trivially the drift of log|x_n|, which doesn't exist on a *finite* group (no escape to infinity). **Extension is STRUCTURALLY_BLOCKED.**

---

## Disposition A: **NO_FIT** (and structurally blocked at extension).

Rationale:
- Group hypothesis fails (matrix-product structure absent).
- Conclusion shape mismatch (Lyapunov rate, not Fourier decay).
- Profinite extension structurally blocked (abelian → trivial Lyapunov).

Not MODE_H_CIRCULAR — the hypothesis is on the matrix structure, not on Fourier decay. Just categorically wrong tool.
