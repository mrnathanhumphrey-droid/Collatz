# Probe R10 — the character ledger (the spectral form) — **FULL PASS (A gate + B/C welds + D/E)**

**Date:** 2026-07-21  Exact rationals via roots-of-unity traces (no cyclotomic field needed — everything closes in
ℚ). Probe `probes/probe_charledger_R10.py` (reuses R7 build_mu/mu1/cram). Gates Wilson's spectral session: **S is a
partial sum of a FROZEN character series.**

**Objects.** Multiplicative group G_r = (1+3ℤ)/(1+3^{r+1}ℤ) = ⟨4⟩, cyclic order 3^r, dlog₄ the iso to ℤ/3^r;
ν = pushforward of Syrac(ℤ/3^r)=μ_r under s↦1+3s; χ_k(X)=ζ^{k·dlog₄(X)}, ζ=e(1/3^r). Spectral form
S_{n+1}=2Σ_{χ∈Ĝ_n}|ν̂(χ)|²/(4χ(4)−1); frozen layer **Λ_r := Σ_{ord(χ)=3^r}|ν̂(χ)|²/(4χ(4)−1)**.

**Exact mechanism (character side only).** |ν̂(χ_k)|² = Σ_u g_r(u)ζ^{ku}, g_r = dlog-domain autocorrelation of ν_r
(rational, fresh from Syrac+dlog). The geometric weight closes as a rational trace:
**A_N(u) := Σ_{ωᴺ=1} ωᵘ/(4ω−1) = N·4^{N−i₀}/(4ᴺ−1)**, i₀=((u−1) mod N)+1. Primitive-k restriction (order exactly
3^r) by inclusion–exclusion: **Λ_r = Σ_u g_r(u)[A_{3^r}(u) − A_{3^{r−1}}(u)]**. Uses dlog table + Syrac(ℤ/3^r)
**only** — never the C-tables, γ-tables, or level r+1 data.

## R10-A — Λ-LEDGER GATE: **PASS**
2Λ_r = OffDiag_{r+1} = S_{r+1} − S_r **exact, r = 1…5**, Λ from the character side alone.

| r | Λ_r | 2Λ_r | OffDiag_{r+1} = S_{r+1}−S_r | verdict |
|---|---|---|---|---|
| 1 | −2/21 | −4/21 | −4/21 | ✅ |
| 2 | **−1490/203889** | −2980/203889 | −2980/203889 | ✅ |
| 3 | +2849957897648150/2159281421340253987 | +5699915795296300/… | = | ✅ |
| 4 | +347914009142242452048362…/… | +695828018284484904096725…/… | = | ✅ |
| 5 | +307926927625750707504348…/… | +615853855251501415008697…/… | = | ✅ |

Anchors hit: Λ₁ = −2/21, 2Λ₂ = −2980/203889 (the banked increment is now the frozen statement Λ₂ = −1490/203889).
**Each ledger increment IS one character layer.** Σ_{r≥1} Λ_r = ½·Σ_{k≥2} OffDiag_k = ½·(−1/5) = **−1/10** (with
Λ₀ = 1/(4−1) = 1/3), so S_∞ = 7/15 ⟺ Σ_{r≥1} Λ_r = −1/10. Walk-back #31 **not** incurred; the spectral form lives.

## R10-B — LAYER MASS (multiplicative cross-thread weld): **PASS**
Σ_{ord(χ)=3^r}|ν̂(χ)|² = S_r **exact, r = 1…5** (2/3, 10/21, 31370/67963, …). The **multiplicative order-3^r
character mass equals the shell mass** — computed as Σ_u g_r(u)·c_{3^r}(u) (Ramanujan), it equals the additive
primitive Plancherel mass of R6. The two interlocked self-similarities (additive-primitive and
multiplicative-order) meet on the same numbers, level by level. A genuine weld, unlike R9-D's two vacated numerals.

## R10-C — TRACE IDENTITY (measure-free): **PASS**
Σ_{χ∈Ĝ_r} 1/(4χ(4)−1) = A_{3^r}(0) = 3^r/(4^{3^r}−1) **exact, r = 1…4**: 1/21, 1/29127, 1/667199944795629,
1/72172920362019897195243695442779408378070742863. **The (4ᴺ−1) Mersenne denominators are literally the trace of
the weight function 1/(4ω−1) over the full character group** — the steer bottomed out.

## R10-D — WITHIN-LAYER PROFILE (measurement, NO fit, NO verdict)
Λ_r actual vs the uniform baseline Λ_r^unif = S_r·⟨w⟩_r (⟨w⟩_r = total primitive weight ÷ 2·3^{r−1}), each column
labeled with what it compares against; palindrome |ν̂(χ_k)|²=|ν̂(χ_{−k})|² **exact at every r**.

| r | Λ_r [actual] | Λ_r^unif [uniform-ν̂ baseline] | excess | note |
|---|---|---|---|---|
| 2 | −0.007308 | −0.003777 | **1.935×** | matches Wilson's pre-reg 1.93× exactly |
| 3 | +0.001320 | −0.0000012 | (baseline ≈0) | sign flip; Λ carried entirely by non-equidist. |
| 4 | +0.000650 | −4e−17 | (baseline dead) | " |
| 5 | +0.000327 | −2e−49 | (baseline dead) | " |

**r=2 is the last level where uniformity is a meaningful baseline (1.935× excess, matching the hand-derivation).**
From r=3 on the uniform baseline ⟨w⟩_r is doubly-exponentially small (the weight-trace collapses as 1/(4^{3^r}−1)),
so **the layer value is carried entirely by within-layer NON-equidistribution of |ν̂(χ)|² over the angles χ(4)** —
R8-A's excess-over-uniform, one representation deeper. (Uniform ν kills every layer r≥1 identically, inheriting
the R8-A uniform kill.) The |ν̂|² profiles: r=2 three distinct values {0.0459, 0.0493, 0.1428}; r=3..5 fan out to
9/27/79 distinct values in bands [0.0012,0.0517] / [0.0003,0.0200] / [0.00004,0.0078]. Λ sign flips −,−,+,+,+
(= the OffDiag −,−,+,+ pattern, since 2Λ_r = OffDiag_{r+1}).

## R10-E — DECAY PROFILE (measurement; raw material, no law fitted)
| r | #prim | max\|ν̂\| | ℓ¹ layer mass Σ\|ν̂\| |
|---|---|---|---|
| 1 | 2 | 0.577350 | 1.154701 |
| 2 | 6 | 0.377924 | 1.628719 |
| 3 | 18 | 0.227389 | 2.674073 |
| 4 | 54 | 0.141569 | 4.664744 |
| 5 | 162 | 0.088201 | 7.938319 |

max|ν̂| contracts ~0.62/layer; ℓ¹ mass grows ~1.7/layer. Raw material for the multiplicative-1.17 question; no rate
claimed.

## Status
**R10 FULL PASS.** The spectral form is verified: **2Λ_r = OffDiag_{r+1} exact r=1…5** from the character side
alone (R10-A gate), so S is a partial sum of the frozen layer series and S_∞=7/15 ⟺ Σ_{r≥1}Λ_r = −1/10. The
multiplicative order-3^r character mass welds to the shell mass S_r (R10-B); the Mersenne (4ᴺ−1) denominators are
the character-group trace of the weight (R10-C); and the within-layer measurement (R10-D) reproduces the
pre-registered 1.93× excess at r=2, then shows the uniform baseline dying (r≥3 carried purely by
non-equidistribution). **Still owed (pen):** the frozen layer sum Σ_{r≥1}Λ_r = −1/10 in closed form — i.e. the
within-layer distribution of |ν̂(χ)|² over χ(4)-angles that produces each Λ_r (equivalently γ_∞ / C̄_∞(j)). No
fitting; exact rationals; the individual |ν̂|² floats and the ℓ¹/max decay are labeled measurement.
