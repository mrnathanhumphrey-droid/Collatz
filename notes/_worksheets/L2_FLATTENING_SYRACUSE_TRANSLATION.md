# Phase 2 — L²-flattening hypothesis articulated for Syracuse μ_n

**Date:** 2026-05-12. L²-flattening structural-compatibility probe, Phase 2.

---

## 1. The Syracuse stationary measure μ_n: precise definition

The relevant chain in Tao 2019 (arxiv 1909.03562) and in this project (see `s_infinity_exact.py`, R76, R77) is the **Syracuse Markov chain on (Z/3^n Z)\*** = {r ∈ {0, …, 3^n − 1} : r ≢ 0 (mod 3)}.

**Transition kernel.** With v drawn from Geom(1/2) conditioned to v ∈ {1, …, M} (where M = 2·3^{n−1}, the order of 2 mod 3^n),
> K_n[r → s] = P((3r + 1) · 2^{−v} ≡ s (mod 3^n))

where 2^{−v} is inverted mod 3^n. Concretely, with normalisation Z_v = 1 − 2^{−M},
> K_n[r → s] = Σ_{v ≥ 1, (3r+1)·2^{−v} ≡ s} (2^{−v} / Z_v).

**Stationary distribution.** μ_n = π_n is the (unique) left eigenvector of K_n with eigenvalue 1, normalized to a probability distribution on (Z/3^n Z)\*:
> π_n K_n = π_n, Σ_r π_n(r) = 1, π_n(r) = 0 for r ≡ 0 (mod 3).

|supp π_n| = 2·3^{n−1}; ambient group order |Z/3^n Z| = 3^n.

**Algebraic content of μ̂_n.** The Fourier coefficient is
> μ̂_n(ξ) = Σ_r π_n(r) e^{2πi r ξ / 3^n}, ξ ∈ Z/3^n Z (or in {0, …, 3^n − 1}).

This is the function Tao Prop 1.17 bounds with `C_A`. The polynomial-in-A Fourier bound is exactly a statement about |μ̂_n(ξ)| for "most" or "all sufficiently high-frequency" ξ.

## 2. The structural mismatch with BKS

The BKS hypothesis requires the measure to live on ℝ^d with a smooth derivative cocycle (Phase 1, §3 and §5). The Syracuse μ_n lives on a **finite discrete group** Z/3^n Z. There is no canonical:

- Ambient ℝ^d the measure embeds into.
- Smooth map f whose derivative cocycle generates the dynamics.
- Affine hyperplane structure for the non-concentration hypothesis to refer to.

Three plausible translations of "L²-flattening for μ_n" can be articulated, each of which has its own structural issue.

## 3. Translation A — direct discrete-group L²-flattening

Restate L²-flattening as a property of μ_n on the finite abelian group G_n = Z/3^n Z, in analogy with Bourgain-Konyagin-style "additive flattening" on F_p:

> **(A-hypothesis).** ‖π_n * π_n‖_2 ≤ N^{−γ} · ‖π_n‖_2 for some γ > 0 independent of n,
> where N = |G_n| = 3^n and ‖·‖_2 is the L²(G_n, counting) norm.

This is the discrete analogue. The connection to Fourier decay is via Plancherel: ‖π_n * π_n‖_2² = Σ_ξ |π̂_n(ξ)|^4, so L²-flattening of π_n controls the **fourth moment** of |π̂_n|, which in turn controls a Fourier-tail bound by Cauchy–Schwarz.

**Problem with Translation A.** The Bourgain–Konyagin / additive-combinatorics L²-flattening on F_p relies on **sum-product** in F_p (the Erdős–Szemerédi inequality |A·A| + |A+A| ≥ |A|^{1+δ}). The Syracuse chain on Z/3^n Z is not a sum-product setting: the chain uses **multiplication by 3 and division by 2^v mod 3^n** — multiplication by 3 is multiplication by 0 (since 3 ≡ 0 mod 3^n is the divisor of the modulus), so the standard sum-product on (Z/3^n Z)\* doesn't structure the chain. Furthermore, BKS's L²-flattening is the Khalil 2305.00527 version, which is for ℝ^d, not for discrete groups.

**Status of Translation A:** the hypothesis is well-posed but the BKS framework's machinery doesn't apply to it. The relevant body of theory for Translation A is the Bourgain–Konyagin → Bourgain–Gamburd → Varjú lineage of L²-flattening on **discrete groups** (already in your literature corpus under "Bourgain-Konyagin/").

## 4. Translation B — embed μ_n into ℝ via the 3-adic Cantor set

Identify Z/3^n Z with the n-digit base-3 expansions {0, 1, 2}^n, and embed into [0, 1] via
> ι: r ↦ Σ_{j=0}^{n−1} d_j(r) · 3^{−j−1}, where r = Σ d_j(r) 3^j.

Then μ_n pushes forward to a probability measure ι_* μ_n on the **3-adic Cantor set** C_3 ⊂ [0, 1]. As n → ∞, the projective limit μ_∞ = lim μ_n is a measure on the full 3-adic integers Z_3, or equivalently (via ι) a probability measure on C_3 ⊂ [0, 1].

In this translation:

- The ambient space is ℝ (so d = 1).
- The Cantor set C_3 is self-similar (it IS an IFS attractor: the IFS f_0(x) = x/3, f_1(x) = (x+1)/3, f_2(x) = (x+2)/3).
- The stationary measure μ_∞ on C_3 is the pushforward of the Markov chain's invariant measure.

**Problem with Translation B (Diophantine condition).** For the BKS framework to deliver **polynomial** decay (rather than only polylog), the Diophantine condition on the IFS contractions is required. The 3-adic Cantor set IFS has contractions (1/3, 1/3, 1/3) — a single contraction ratio. This is the **trivial** Diophantine case: log(1/3) = log(1/3), no Diophantine non-resonance. The BKS framework here delivers only **polylogarithmic** decay (the abstract's worst case). The decay rate is **(log|ξ|)^{−η}**, which is far too slow to beat Tao's effective `C_A`.

**Furthermore**, the BKS framework requires the measure on the IFS attractor to be the **self-similar measure with respect to weights p_i** (i.e., μ = Σ p_i (f_i)_* μ as a fixed-point equation). The Syracuse μ_∞ is **NOT** self-similar in this sense: its 3-adic digit distribution is **non-i.i.d.** (the Markov chain has non-trivial transitions, the digits are correlated through the Syracuse dynamics; see e.g. Matthews 1989 in the bundle's Collatz-specific section).

**Status of Translation B:** the embedding is well-defined but (i) the Diophantine condition is trivially violated by the homogeneous IFS, and (ii) the Syracuse measure is not the Bernoulli self-similar measure on C_3.

## 5. Translation C — L²-flattening as a Markov-chain mixing statement

Reinterpret L²-flattening as the L²-mixing rate of the Syracuse chain:
> **(C-hypothesis).** ‖π_n − u_n‖_2 → 0 as n → ∞, where u_n is the uniform measure on the chain's state space, with a polynomial decay rate in some parameter (mixing time? chain depth?).

This is the Foster–Lyapunov / spectral-gap angle (the "drift-condition fallback" in INDEX). The Syracuse chain on (Z/3^n Z)\* has a finite state space of size 2·3^{n−1}, hence a finite L²-gap; the question is whether this gap is uniform in n.

**Problem with Translation C.** This is a statement about the **spectrum of the transition operator K_n**, not about the **Fourier transform of the stationary measure π_n**. They are related (Plancherel + spectral decomposition of K_n on L²(G_n)) but not identical. More importantly, this translation **moves the question outside the BKS framework**: BKS doesn't use Markov-chain mixing as its L²-flattening hypothesis; it uses Khalil's affine non-concentration of derivative-cocycle distributions. The two are conceptually different L²-flattening statements.

**Status of Translation C:** well-defined, but it's the **Foster–Lyapunov / drift-condition** route (the parallel "Probe 2" or "fallback" in INDEX), not the BKS L²-flattening route.

## 6. Connection to existing project work

The project already has substantial structural understanding of K_n:

- R74–R77 (this session's earlier results): T_lead = (1/45)·[[7,9],[28,36]] is a 2×2 reduction of the chain's dynamics on the class-resolved Plancherel-squared mass; its eigenvalue 43/45 governs the rate-(1/2) deviation envelope.
- R76 §11, R77 §1: T_diag = (1/5)·[[1,1],[4,4]] rigorously controls the diagonal of the recursion.
- R76: the **conservation law** Σ_j M_{n+1}(η_0 + j·3^n) = 0 is a structural identity on |μ̂_n|^2 sums.
- The ε_n cache (`experiments_output/result_77_7_eps_exact_through_k7_v2.json`) gives exact rationals through n = 7.

These are all statements about the **spectral structure of K_n** and the **Plancherel-sum of |μ̂_n|**, not about ‖μ_n‖_2 or ‖μ_n * μ_n‖_2 in the BKS-Khalil sense. The project has thus already mapped out Translation C structurally (Foster–Lyapunov / Markov-chain-mixing flavored work), but Translations A and B are unexplored and structurally blocked for the reasons above.

## 7. Summary — what would need to be true for BKS to apply

For Baker–Khalil–Sahlsten 2407.16699 to **deliver polynomial decay on μ_n via L²-flattening**, all of:

1. μ_n would need to embed naturally into ℝ^d (some d ≥ 1) as a **continuous** measure (or as a measure on a smooth-IFS attractor).
2. The embedding would need to come with a **smooth derivative cocycle** D_x f generated by the chain dynamics (not the discrete Markov kernel directly).
3. The derivative-cocycle distribution would need to satisfy **uniform affine non-concentration** (Khalil 2305.00527 hypothesis).
4. The IFS / cocycle would need to satisfy a **Diophantine non-resonance** condition strong enough to upgrade polylog → polynomial decay.

The Syracuse μ_n on Z/3^n Z satisfies **none** of these out of the box:

- (1) is satisfied only in the Translation-B sense (embedding via 3-adic digits), and only as a measure on the trivially-Diophantine 3-adic Cantor set.
- (2) fails: the chain dynamics generate **discrete digit transitions** on a finite group, not a smooth cocycle.
- (3) is not even well-posed in Z/3^n Z without further structure.
- (4) is the worst case: homogeneous contraction ratios = trivially Diophantine, BKS framework gives only polylog.

The structural mismatch is fundamental to the framework's design, not a bookkeeping or technical issue.

---

End Phase 2.
