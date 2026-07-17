# Pascadi 2025 (arXiv:2511.08445) — Extraction disposition

**Date:** 2026-05-11. Reader: extraction agent. Paper: Pascadi, "Non-abelian amplification and bilinear forms with Kloosterman sums."

## Disposition

> **EXTRACTION_FAILS_STRUCTURAL**

## One-paragraph rationale

Pascadi's main bilinear theorem bounds `Σ_{m, n} α_m β_n · S(am, n; c)` — a sum of the **classical Kloosterman sum** across two free integer interval variables, with arbitrary complex weights, for composite moduli `c`. The saving comes from a **non-abelian amplifier on SL₂(Z/cZ)** built over the normal subgroup `Γ_c(d)`, combined with combinatorial counting of word equations in SL₂(Z/p^k). Our object — `Σ_{a ≡ 1 mod p in Z/p^r} 1̂(p·a) · F̂(p·a)` — has (i) **one** coset variable, not two; (ii) **no Kloosterman sum** anywhere; (iii) its summation group is the **abelian** quotient (Z/p^r)^×, in which Pascadi's amplifier becomes trivial by Pascadi's own §2.3 remark (the abelian analog weights `χ̄'(ℓ)χ(ℓ)`, which is zero when only `{I}` is available). The lit-scan "method-shape match" was real at the keyword level (bilinear, amplification, composite moduli, Fourier) but false at the mechanism level: Pascadi's saving is harvested from the non-commutative arithmetic of SL₂(Z/cZ) on a Kloosterman kernel — neither of which is present in our setting.

## Comparison against GY (the prior structural failure)

Pascadi **resolves all three** of the obstructions that killed GY:

| Obstruction | GY | Pascadi |
|---|---|---|
| (A1) Smooth-weight requirement | Required Ŵ-rapid-decay; killed by our Dirichlet kernel | **No smoothness** required (arbitrary `(α_m), (β_n)`) |
| (A2) Scope `r ≥ 3` | `d² | q | d³` forced `r ∈ {1, 2}`; killed our r = 8..20 | **Admits all prime powers** via factorization c = dd'e with e = 1 |
| (A3) AFE / L-function input | Required AFE-supplied `1/√(mn)` amplitudes | **No AFE, no L-function** input; saving purely combinatorial + character-theoretic |

But Pascadi introduces a **new, more fundamental obstruction** at the level of object shape:

| Obstruction | Status |
|---|---|
| (A4) Cardinality vs amplitude×phase | **New shape**: Pascadi's saving is **combinatorial counting in SL₂(Z/p^k Z)**, neither cardinality nor amplitude×phase. The amplifier `Σ_{ℓ ∈ Γ_c(d)} ρ'(ℓ)*⊗ρ(ℓ)` lives on a **NORMAL subgroup of a NON-ABELIAN group**. Our problem is abelian; Pascadi's §2.3 explicitly notes the abelian analog of the amplifier is trivial. |
| Object shape | **Fundamental mismatch**: Pascadi takes two interval variables `m, n` and the classical Kloosterman sum `S(am, n; c)` across both. We have one coset variable `a` and no Kloosterman sum. |
| (A5) Magnitude / order-of-magnitude | Pascadi's saving is δ-improvement over `c^{1+o(1)}` (sub-Polya-Vinogradov); our target is √N · √q saturation. Different sides of "trivial bound", not comparable rates. |

## Where the mechanism diverges (specific citation)

Pascadi §2.3 — exact quote from extracted HTML:

> "To the best of our knowledge, [this amplifier construction is] the first instance of such a construction [in the non-abelian setting]. For abelian (Dirichlet) characters, the analogous amplifier weights `χ̄'(ℓ)·χ(ℓ)`, which is trivial when only `{I}` is available."

This is the load-bearing comment. Our problem is summed over the **abelian** group `(Z/p^r)^×`, restricted to the principal-unit coset `{a ≡ 1 (mod p)}`. The image of this coset under the diagonal embedding `(Z/p^r)^× → SL₂(Z/qZ)` is **NOT a normal subgroup of SL₂(Z/qZ)** (the only normal subgroups in this range are the congruence kernels `Γ_q(d)`, which contain off-diagonal elements). So Pascadi's amplifier construction does not apply to our coset.

Furthermore Pascadi §2.2 (Kloosterman matrix Fourier transform): the **first step** of the proof takes the c×c Kloosterman matrix `K = (S(m,n;c))` and applies unitary Fourier transforms in `m, n`. This step has no analog for our object, which is **not a matrix** indexed by two variables — it is a single sum over a coset of a single F̂·1̂ product.

## Refs

- [PASCADI_TRANSLATION.md](PASCADI_TRANSLATION.md) — notation correspondence (Pascadi vs ours), side-by-side bilinear sum shapes, scope check showing `N = p^{r-1}` exceeds Pascadi's `≪ √c` allowed range for r ≥ 4.
- [PASCADI_MECHANISM.md](PASCADI_MECHANISM.md) — step-by-step trace of where the `c^{-1/12}` saving emerges in Pascadi (§2.2 → §2.3 → §3.3 → §6 → §2.5), the three input requirements (Kloosterman matrix, two interval variables, range `≪ √c`), and the three things the mechanism does not do.
- [PASCADI_TRANSLATION_ATTEMPT.md](PASCADI_TRANSLATION_ATTEMPT.md) — four direct-translation attempts (all fail at attempt 1.3a-1.3d before reaching a starting line), adversarial checks A1-A6 (A1/A2/A3 resolve, A4 introduces new structural obstruction, A5 wrong order of magnitude, A6 honest scope says "not translation, new mathematics").

## What this leaves open / does not rule out

This disposition rules out **direct or nearly-direct translation of Pascadi's published method** to our problem. It is silent on:

- Whether some entirely different abelian-amplifier construction tailored to (Z/p^r)^× and the principal-unit coset can give √N saturation. The DFI / Heath-Brown abelian amplification is a natural precursor, already not known to deliver in the depth-r regime.
- Whether the explicit Cochrane-log closed form `G(a) = √q · e_q(P_a(s*(C_a)))` (Theorem 78.6) admits direct evaluation of `Σ_a 1̂(p·a) · G(a)` without external bilinear-bound machinery. This is a "do the computation" path, not a "cite a theorem" path.
- Whether a future generalization of Pascadi's non-abelian amplifier to non-normal subgroups, or to non-Kloosterman kernels, could be developed. Currently not visible in the published method.

After Garcia-Young and Pascadi (the two freshest method-shape candidates from the lit scan), **the published-literature space of nearly-direct bilinear extractions for our √N target appears to be exhausted at the "method-shape match" level**. The remaining paths are either: (i) parallel construction / new theorem, or (ii) direct evaluation using the R78.4-78.6 explicit closed form. The "translate an existing theorem" path is closed by this disposition + the GY disposition together.
