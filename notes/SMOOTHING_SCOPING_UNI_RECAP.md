# Phase 1: ARHW UNI condition — exact form needed for the smoothing scoping probe

Source: arxiv:2306.01275v2 extracted via pypdf to C:/tmp/arhw_full.txt. All quotations preserve OCR artifacts verbatim.

## (a) UNI condition statement — what the cocycle needs to satisfy

The load-bearing form used by Theorem 2.8 is **condition (10), p. 10**, on the induced IFS Φ_N supplied by Claim 2.1:

> "**There exists m', m > 0 such that:** For every x ∈ [0, 1], for both i = 1, 3,
>
>     m ≤ |d/dx (log f'_i − log f'_{i+1})(x)| ≤ m', and m − 2·C̃·sup_{f∈Φ} ||f'||_∞ > 0.    (10)"

This is the UNI condition in its operative form. Theorem 2.4 part (5) lifts it to "UNI in all parts" of the random model:

> "**5. (UNI in all parts)** There exist m', m > 0 and N_0 ≥ 0 such that for all N ≥ N_0, for every ω ∈ Ω there exist α_1^N, α_2^N ∈ X_N^{(ω)} such that
>
>     m ≤ |d/dx (log f'_{α_1^N} − log f'_{α_2^N})(x)| ≤ m', for all x ∈ [0, 1]."

The equivalence with non-conjugate-to-linear (Claim 2.2 proof, p. 7):

> "if for all σ-periodic ξ, ζ ∈ A^N and x ∈ K we have
>
>     lim_n d/dx log f'_{ξ|n}(x) = lim_n d/dx log f'_{ζ|n}(x),    (6)
>
> then Φ is C² conjugate to a linear IFS."

**Therefore UNI ⟺ ∃ two σ-periodic codings ξ, ζ with distinct limit-derivatives of log f'.** Verifying UNI for a candidate δ_a means exhibiting two such σ-periodic codings and computing the limit derivative explicitly.

## (b) Quantitative dependence — does α(UNI strength) → 0 as UNI strength → 0?

ARHW Theorem 2.8 (p. 17, verbatim):

> "Then there exist C, γ, ε, R > 0 and some 0 < α < 1 such that for all |b| > R, a ∈ R with |a| < ε, and n ∈ N
>
>     ||P^n_{a+ib}||_{C¹} ≤ C · |b|^{1+γ} · α^n.    (21)"

The 0 < α < 1 in Theorem 2.8 is the **L²-contraction rate of the transfer operator's pieces** (Proposition 2.9). Tracing the proof chain:

- Proposition 2.9 contraction α depends on Naud-style oscillatory-integral estimates in the L² norm.
- Those estimates depend on **the UNI lower bound m** (the lower bound in condition (10)).
- Specifically, in Naud's framework the rate of L² decay degrades polynomially in 1/m: as m → 0 the cancellation in ∫ e^{2π·i·b·c(I,x)} g ◦ f_I dp(I) becomes weaker, and α(m) → 1.

**ARHW does NOT extract an explicit functional form α(m).** The paper's exponent α in the final Fourier-decay bound |F_q(ν)| = O(|q|^{−α}) is implicit through (γ, α_contraction, ε, R) — no closed-form dependence on m, m', ρ, ρ_min is given.

**Critical implication for smoothing-route limit:** If a Syracuse-derived δ_a has UNI strength m(ε) → 0 as the smoothing parameter ε → 0, then α(ε) → 0, and the Fourier-decay exponent vanishes in the limit. This is the H_DELTA_EXISTS_BUT_UNI_DEGENERATE failure mode the pre-registration anticipates.

The paper's empirical claim (p. 2): "easily verifiable conditions when Theorem 1.1 may be applied" — conditions (9) and (10) — are existence conditions, not quantitative-strength conditions. No theorem in the paper says "α ≥ f(m) for explicit f."

## (c) C² regularity requirement: piecewise or global?

p. 1–2, verbatim:

> "let Φ = {f_1,...,f_n} be a finite set of strict contractions of a compact interval I ⊆ R (an IFS - Iterated Function System), such that every f_i is differentiable. We say that Φ is C^α smooth if every f_i is at least C^α smooth for some α ≥ 1."

p. 10, on building the induced IFS Φ_N:

> "Fix a C² IFS Φ as in Section 2.1, and let us retain the other assumptions and notations from that Section."

The framework requires **each f_i ∈ Φ to be C² globally on the compact interval I = [0,1]** (or equivalently, C² on a neighborhood containing [0,1]). The "piecewise" framing is not directly supported: each f_i is a single C² map on the whole interval, and the IFS structure builds the attractor by iteration of these (globally-smooth) maps.

This rules out δ_a constructions that are C² only off a Cantor set or only piecewise on [0,1] with C⁰ or C¹ joins.

## Summary

| Aspect | Requirement | Comment |
|---|---|---|
| UNI form | (a) Two σ-periodic codings ξ, ζ with distinct lim d/dx log f' OR (b) condition (10): m ≤ |d/dx(log f'_i − log f'_{i+1})| ≤ m' on [0,1] with m > 2·C̃·ρ | Equivalent up to passing to induced IFS Φ_N |
| α dependence on UNI strength | α(ν) > 0 but **not extracted in closed form**; degrades as m → 0 | Limit ε → 0 would send α → 0 if m(ε) → 0 |
| C² regularity | Each f_i ∈ C²([0,1]) globally | Piecewise constructions excluded |
| Uniform contraction | 0 < inf|f'| ≤ sup|f'| < 1 on [0,1] | Non-negotiable |
| Disjointness of cylinders | Condition (9): at least one disjoint triple of 1st-gen cylinders | Achievable for small δ_a via Claim 2.1 induced IFS |
