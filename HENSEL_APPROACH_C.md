# HENSEL Phase 2 Approach C — Vinogradov mean-value direct

**Date:** 2026-05-11. Analyst: Wilson. Reports to: Nathan.

## Disposition: NOT_TRIGGERED (Approach A succeeded provisionally)

Approach C was queued as fallback if Approaches A and B failed. Approach A produced a structurally exact Hensel-lifted closed form. C is not needed for closure at family level.

However, this document records what Approach C would look like, for completeness and as comparison with the literature.

## VMV reformulation

The original Approach C asked: write `|Σ 1̂(p·a) · F̂_p(p·a)|²` as a 6th-moment-like object and apply Bourgain-Demeter-Guth 2016 (cubic VMV).

The squared bilinear:
> `|S_p|² = Σ_a Σ_b 1̂(p·a) · conj(1̂(p·b)) · F̂_p(p·a) · conj(F̂_p(p·b))`

Substituting `F̂_p(p·a) = p·e_q(c)·G_p(a)` and Approach A's closed form `G_p(a) = √q · η_p · e_q(P_a(s*(r)))`:

> `|S_p|² = p²·q·|η_p|² · Σ_a Σ_b 1̂(p·a) · conj(1̂(p·b)) · e_q(P_a(s*_a) − P_b(s*_b))`

With |η_p| = 1 and substituting `P_a(s*(r))` polynomial in s_a = s*(r)(C_a) of degree r:

> `|S_p|² = p²·q · Σ_{(s_a, s_b)} 1̂(p·a(s_a))·conj(1̂(p·a(s_b))) · e_q(Polynomial(s_a) − Polynomial(s_b))`

The phase difference `Polynomial(s_a) − Polynomial(s_b)` is a polynomial of degree r in (s_a, s_b) jointly. After expanding mod p^{r+1}:

- Leading order p²: `-(s_a² − s_b²)/2 = -(s_a + s_b)·(s_a − s_b)/2`
- Order p³: `(s_a³ − s_b³)/6 − (s_a·s_a_2 − s_b·s_b_2)·something`
- ...

This is the **6-th moment object** in VMV language: the squared sum has a "lattice point counting" interpretation where (s_a, s_b) range over Z/p^{r-1}.

## BDG 2016 applicability

Bourgain-Demeter-Guth's 2016 Annals paper proved the sharp Vinogradov mean-value exponent for cubic phases on integer intervals `[1, N]`:
> `||Σ_{n=1}^N e(α_1·n + α_2·n² + α_3·n³)||_{L^6} ≤ N^{3/2+ε}`

This bounds 6-moment objects with cubic polynomials of degree 3 in n on the integer interval.

**Our object is different:**
- Modulus p^{r+1} (finite group), not integer interval.
- Polynomial degree r (not 3) — degrees 2, 3 only at r=3, but higher at r ≥ 4.
- s_a ranges over a multidimensional digit space, not a 1-D interval.

**Literature analog:** the "Heath-Brown hybrid" framework treats cubic exponential sums over prime power moduli, which IS our setup at r=3 with the cubic phase. R79b walks back the cubic-character framing because at r ≥ 4 the phase isn't strictly cubic in `a`.

For r ≥ 4 with the Hensel-corrected closed form's polynomial of degree r in s*(r), Heath-Brown's machinery wouldn't directly apply — the phase is HIGHER-degree.

## BDG explicit-constants question

For Approach C to produce a useful bound, BDG's exponent must be combined with EXPLICIT CONSTANTS. BDG 2016 proves the right exponent (sharp Vinogradov conjecture) but the implicit constants in the bound are NOT determined in their paper — they depend on the polynomial structure.

To close the bilinear bound to strict 2√N at r ≥ 4 via Approach C, we'd need:
1. A specific BDG-style decoupling estimate applied to the prime-power modulus + multi-digit polynomial.
2. Explicit constants matching the family-level required uniformity.

**This is a research-level paper of its own.** Not achievable in a single session even if Approach A had failed.

## Verdict

> **APPROACH_C_NOT_TRIGGERED.** Approach A's direct method succeeded; Approach C's VMV route is unnecessarily complex and would only produce structural-match-only verdicts (per Milicevic-Banks analog) without explicit constants. If Approach A's numerical verification fails AND Approach B fails to reveal the issue, Approach C should be referred to as VMV-LITERATURE-LEVEL work, not session-level.

## Adversarial check (A4) for VMV literature

Even if Approach C were triggered, BDG 2016 gives:
- **Right exponent (Vinogradov conjecture):** sharp, established.
- **Explicit constants:** NOT established at the level needed for family-level uniform bilinear closure.
- **Setup:** integer intervals + cubic, not prime-power + multi-digit polynomial.

Per the brief's adversarial check (A4): VMV-LITERATURE-LEVEL work would be required, beyond competent careful manipulation. **Flag for INCONCLUSIVE if Approach C were the only option.**
