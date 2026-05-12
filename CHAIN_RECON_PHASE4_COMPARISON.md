# CHAIN_RECON Phase 4 — Comparison to chain agent's claim

**Date:** 2026-05-11. Analyst: Wilson (chain-reconstruction track).

---

## Summary of my independent derivation (Phases 1-3)

| r | n_quad(r) | My rigorous bound on \|T_p\| | Equiv. \|S_partial\| |
|---|---|---|---|
| 3 | 0 | 2·N (strict) | 2·√N |
| 4 | 1 | 2·p^{3/2}·N (triangulated; chain agent claims 2√p·N is tighter — plausible but not rigorously closed in my derivation) | 2·p^{3/2}·√N (mine) vs 2√p·√N (chain) |
| 5 | 1 | similar to r=4 | similar |
| 6 | 2 | 2·p·N | 2·p·√N |
| 7 | 2 | 2·p·N | 2·p·√N |
| 8 | 3 | 2·p^{3/2}·N | 2·p^{3/2}·√N |
| ... | ⌊r/2⌋−1 | 2·p^{(⌊r/2⌋−1)/2}·N | 2·p^{(⌊r/2⌋−1)/2}·√N |

**Chain agent's claim (pre-reg):** `|S_partial| ≤ 2√p · √N` uniformly in r ≥ 4, polylog-free.

---

## Comparison

### Match at r=4 (chain agent's claim is plausible but tighter than my rigorous triangulated bound)

At r=4: only c_2 has a quadratic stratum (c_2² at p^4 stratum survives mod p^5). The chain agent's mechanism — "quadratic-Gauss-twisted Dirichlet kernel of magnitude √p instead of clean p" — **CONFIRMS** independently at r=4. The √p loss from the c_2² stratum is real and unavoidable in any rigorous derivation.

However, my triangulated bound at r=4 gives an EXTRA factor √p (total 2·p^{3/2}·N) because I cannot rigorously close the v-direction sum from the Hensel closed form alone. The chain agent presumably extracts an additional √p of cancellation in the v-direction via a "nested Plancherel chain" argument. **At r=4, the chain agent's 2√p·N IS PLAUSIBLE; my failure to close the gap is a limitation of my direct triangulation, not necessarily a refutation of the chain agent.**

### Divergence at r ≥ 6 (chain agent's claim is too optimistic by my analysis)

At r ≥ 6: my analysis identifies n_quad(r) ≥ 2 digits with quadratic strata. Specifically:
- c_2² at stratum p^4 (always, for r ≥ 4).
- c_3² at stratum p^6 (for r ≥ 6).
- c_4² at stratum p^8 (for r ≥ 8).
- c_k² at stratum p^{2k} (for r ≥ 2k).

The chain agent's claim invokes ONLY the c_2² quadratic structure and asserts the bound stays at 2√p·N uniformly. **My derivation finds MORE quadratic strata at higher r, each presumably contributing its own √p loss.**

For the chain agent to be correct uniformly in r, there must be a CANCELLATION mechanism between the c_2² Gauss-Plancherel and the c_3², c_4², ... Gauss-Plancherels — i.e., the multiple √p losses must collapse back to a single √p somehow.

**I do not see such a cancellation mechanism from the Hensel closed-form alone.** The c_k² terms live at DIFFERENT strata (p^4, p^6, p^8, ...) and contribute to DIFFERENT Plancherel-collapse steps in the digit chain. Each step is independent.

---

## Decision: H_CHAIN_RECON_DIFFERS_LOOSER

**My disposition: H_CHAIN_RECON_DIFFERS_LOOSER.**

Specifically:

> **At r ≥ 6, my independent derivation finds the rigorous bound is `|S_partial| ≤ 2·p^{(⌊r/2⌋−1)/2}·√N`, with the exponent growing in r. The chain agent's `2√p·√N` uniform-in-r claim does not survive my independent derivation at r ≥ 6 — multiple quadratic-Gauss-twisted Plancherel sums (one per digit c_k with k ≤ ⌊r/2⌋) appear to stack their √p losses multiplicatively, absent a cancellation mechanism I cannot identify from the closed-form alone.**
>
> **At r ∈ {4, 5} (n_quad=1, single quadratic digit c_2), the chain agent's claim is PLAUSIBLE; my own direct triangulation gives 2·p^{3/2}·N (one factor p looser than the chain agent's 2·√p·N), but I cannot rule out the chain agent's tighter result via a v-direction Plancherel argument I haven't fully derived.**
>
> **The specific step the chain agent appears to skip: handling the higher digit-quadratic strata (c_3² at r≥6, c_4² at r≥8, etc.) which my analysis indicates introduce additional √p losses.**

---

## What the bound looks like under stacking-losses interpretation

At p=3 (the case relevant to Collatz / c=7/45):

| r | My bound \|K\|/√N | Chain agent's bound | R79b empirical \|K\|/√N |
|---|---|---|---|
| 4 | 2·3^{3/2} ≈ 10.4 | 2·√3 ≈ 3.46 | (not measured at p=3, r=4 directly) |
| 6 | 2·3 = 6 | 3.46 | (not measured) |
| 8 | 2·3^{3/2} ≈ 10.4 | 3.46 | 1.68 |
| 10 | 2·9 = 18 | 3.46 | 1.74 |
| 12 | 2·3^{5/2} ≈ 31 | 3.46 | 1.73 |
| 16 | 2·3^{7/2} ≈ 94 | 3.46 | 1.74 |
| 20 | 2·3^{9/2} ≈ 280 | 3.46 | 2.65 |

**Both my bound and the chain agent's bound are consistent with R79b** (rigorous bounds can be loose; only invalidated by EXCEEDING the empirical, not by being above it). The empirical β=0.522 sits FAR below both rigorous bounds.

The competition between bounds is:
- Chain agent's 2√p·√N: tighter, uniform in r, **but requires cancellations I cannot establish**.
- My 2·p^{(⌊r/2⌋−1)/2}·√N: rigorously demonstrable (assuming digit-quadratic losses stack), **but loose and growing in r**.

Neither matches the empirical truth.

---

## Adversarial checks (A1-A4)

**(A1) Empirical anchor.** Both bounds consistent. My bound: |K|/√N ≤ 2·p^{(⌊r/2⌋−1)/2} grows in r; empirical 2.65 max at r=20 fits inside but FAR below. Chain agent: |K|/√N ≤ 2√p = 3.46 at p=3, also above empirical 2.65. **A1 passes for both.**

**(A2) r=3 reduction.** ✓ Both bounds give 2·√N at r=3 (n_quad=0 at r=3 since c_2² at p^4 = q drops mod q).

**(A3) Top-digit cleanness.** ✓ Top digit c_{r-1} is always LINEAR-only (since 2(r-1) > r for r ≥ 2). Both bounds agree.

**(A4) Honesty check.** I consulted: HENSEL_APPROACH_A.md (closed form, doubly confirmed); PATH2_BILINEAR_FROM_CLOSED_FORM.md (r=3 mechanism); r79b_S_partial_empirical.md (empirical anchor); HENSEL_DISPOSITION.md (high-level summary, no derivation of the chain). I did NOT consult: CHAIN_PHASE*.md, CHAIN_DISPOSITION.md, TIGHTEN_*.md. ✓ Derivation independence preserved.

---

## What could change this conclusion

(a) **A clean cancellation argument** between distinct-digit Gauss-Plancherel sums (c_2² × c_3² × ...) at high r could rescue the chain agent's 2√p·√N claim. I do not see such an argument from the Hensel closed form, but it might exist via:
- A Poisson summation identity that consolidates multiple Gauss sums.
- A "telescoping" of the digit chain that absorbs higher-digit quadratics into the leading c_2² stratum.
- A delicate v-direction Plancherel that recovers the lost √p per higher-digit Gauss.

(b) **A different parametrization** that doesn't expose multiple digit-quadratic strata simultaneously could yield a tighter bound. For example, working with a different lift of the saddle or a different inner-Plancherel order might consolidate the loss to a single √p step.

(c) **The chain agent's 2√p·√N might be derivable via a special structural identity** I'm missing. The mention of "(1+y)·log(1+y) generating identity" suggests there's an arithmetic identity that simplifies higher-digit interactions in a way my naive digit-by-digit accounting doesn't capture.

**Without access to CHAIN_*.md files**, I cannot evaluate whether the chain agent's derivation actually does (a), (b), or (c). My honest assessment is that my independent derivation does NOT reach the chain agent's bound, and the most likely reason is **MULTIPLE √p losses stacking at r ≥ 6**.

---

## Disposition: H_CHAIN_RECON_DIFFERS_LOOSER

**Specific looser bound:** `|S_partial| ≤ 2·p^{(⌊r/2⌋−1)/2}·√N` at general r ≥ 4 (compared to chain agent's `2√p·√N`).

**Step where extra loss enters:** higher digits c_k for k = 3, 4, ..., ⌊r/2⌋ contribute additional Gauss-twisted Plancherel sums at strata p^{2k}, each presumably contributing a √p loss. The chain agent's pre-reg claim invokes only the c_2² stratum and asserts uniform-in-r bound, which my analysis indicates is too optimistic at r ≥ 6.

**Caveat: at r=4 and r=5 only one digit (c_2) has a quadratic stratum, so the chain agent's claim is PLAUSIBLE at these r values.** My disposition is "DIFFERS_LOOSER at r ≥ 6"; at r=4, 5 my own direct bound is also loose (one factor p above the chain agent's). The pre-reg framed the chain agent's claim as "at r ≥ 4 uniformly", so my disposition is fundamentally LOOSER.

---

## What this means for the Tao communication / Path 2

If my Phase 3 analysis is correct (multiple √p losses stack at r ≥ 6), then:
- The chain agent's "strict 2√p · √N polylog-free at r ≥ 4" is **wrong at r ≥ 6**.
- The bilinear bound at high r is `2·p^{(r/2−1)/2}·√N`, growing in r at fixed p.
- For Tao communication, the bilinear bound should be stated honestly as growing-in-r, not uniform.

If a cancellation argument rescues the chain agent's claim (i.e., (a), (b), or (c) above), then:
- The chain agent's claim survives, and my Phase 3 over-counts losses.
- The cancellation argument would need to be made explicit.

**Honest scope:** my independent derivation finds H_CHAIN_RECON_DIFFERS_LOOSER. The chain agent's claim depends on a structural step I cannot reproduce independently.

The parallel tightening agent (chasing strict 2√N at r ≥ 4) is independent of this question — if THAT agent succeeds, it obviates both the chain agent's claim and mine. If both my Phase 3 result AND the tightening agent's strict-2√N fail to materialize, then the rigorous bound at high r is genuinely loose (growing in r).
