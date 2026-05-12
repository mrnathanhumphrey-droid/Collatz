# CHAIN_RECON Disposition

**Date:** 2026-05-11. Analyst: Wilson (chain-reconstruction adversarial track).

---

## Disposition: **H_CHAIN_RECON_DIFFERS_LOOSER**

> Independent re-derivation of the nested inner-Plancherel chain at r ≥ 4 produces a LOOSER bound than the chain agent's claim of `|S_partial| ≤ 2√p · √N` uniform in r.
>
> Specifically, my Phase 3 analysis identifies n_quad(r) = ⌊r/2⌋ − 1 digit-quadratic strata at general r ≥ 4 (the digits c_k for k = 2, 3, ..., ⌊r/2⌋ each have a c_k² term surviving at stratum p^{2k} mod q=p^{r+1}). Each digit-quadratic produces a √p Gauss-twisted Plancherel sum (the same mechanism the chain agent invokes for c_2). Under the natural assumption that these losses STACK MULTIPLICATIVELY (no cancellation between distinct-digit Gauss sums), the bound is:
>
> > **|S_partial| ≤ 2·p^{(⌊r/2⌋ − 1)/2}·√N at r ≥ 4**
>
> growing polynomially in r at fixed p.
>
> The chain agent's `2√p·√N` (uniform in r) matches my analysis only at r ∈ {4, 5} (where n_quad=1, single quadratic digit c_2). At r ≥ 6, my analysis identifies additional quadratic digits (c_3 at r ≥ 6, c_4 at r ≥ 8, etc.), each presumably contributing its own √p loss. Absent a cancellation mechanism I cannot identify from the Hensel closed-form alone, the chain agent's uniform-in-r claim is too optimistic.

---

## Confirmation: I did NOT read the contaminated files

Confirmed: I did NOT read any of:
- CHAIN_PHASE1_R3_RECAP.md
- CHAIN_PHASE2_DIGIT_EXPANSION.md
- CHAIN_PHASE3_HIGHER_DIGIT_COLLAPSES.md
- CHAIN_PHASE4_NESTED_CHAIN.md
- CHAIN_PHASE5_NUMERICAL.md
- CHAIN_DISPOSITION.md
- Any TIGHTEN_*.md files

I DID consult (per task permissions): HENSEL_APPROACH_A.md (doubly-confirmed closed form), PATH2_BILINEAR_FROM_CLOSED_FORM.md (r=3 base case mechanism), r79b_S_partial_empirical.md (empirical anchor), HENSEL_DISPOSITION.md (top-level summary, no derivation), result_78_FINAL.md (foundational theorems).

---

## Precise looser bound and step that introduces extra loss

**Precise bound:** `|S_partial| ≤ 2 · p^{(⌊r/2⌋ − 1)/2} · √N` at general r ≥ 4.

| r | n_quad(r) | My bound shape | Chain agent's claim |
|---|---|---|---|
| 3 | 0 | 2·√N (strict) | 2·√N (agrees) |
| 4 | 1 | 2·√p·√N (plausibly; my direct triangle gives 2·p^{3/2}·√N) | 2·√p·√N |
| 5 | 1 | 2·√p·√N | 2·√p·√N |
| 6 | 2 | 2·p·√N | 2·√p·√N (LOOSER by factor √p) |
| 7 | 2 | 2·p·√N | 2·√p·√N (LOOSER by √p) |
| 8 | 3 | 2·p^{3/2}·√N | 2·√p·√N (LOOSER by p) |
| 20 | 9 | 2·p^{9/2}·√N | 2·√p·√N (LOOSER by p^4) |

**Step introducing extra loss:** the digit c_k for k = 3, 4, ..., ⌊r/2⌋ each have c_k² entering P_a(s*(r)) at stratum p^{2k}, with that stratum surviving mod p^{r+1} when 2k ≤ r. Each such digit produces a quadratic-Gauss-twisted Plancherel sum of magnitude √p (analog of the chain agent's c_2² mechanism at p^4 stratum), and these losses appear to stack independently across digits since they live at different strata of the phase decomposition.

---

## Why the chain agent's claim might still be right (caveats)

My disposition is LOOSER, but there are scenarios where the chain agent could be correct:

**(a) Cancellation between distinct-digit Gauss-Plancherels.** The c_2² and c_3² (and higher) Gauss-twisted sums might interfere with each other, with the combined magnitude saturating at single √p rather than stacking. I do NOT see such an interference from the Hensel closed-form structure, but it might exist via a "telescoping" identity or Poisson summation that consolidates multiple Gauss sums.

**(b) Tighter v-direction Plancherel.** At each digit step in the chain, my direct triangulation loses a factor of p in the "u-direction" sum (the v variable parametrizing the restricted u-residue-class). A more careful Plancherel against the e_{p^{r-k}} Dirichlet kernel could recover this factor, in which case the chain agent's claim becomes 2√p · √N strict.

**(c) "(1+y)·log(1+y) generating identity" providing structural simplification.** The chain agent's pre-reg mentions this identity as enabling the uniform-in-r save. If the identity produces an algebraic consolidation of higher-digit quadratics into a single effective c_2² stratum, the chain agent is correct. I do not see such consolidation from my direct expansion.

**My honest assessment:** absent the chain agent's specific argument for (a), (b), or (c), my independent derivation finds the bound growing in r at fixed p. The chain agent's uniform-in-r claim depends on a structural step I cannot reproduce.

---

## What would resolve the disagreement

To settle whether the chain agent's claim is correct OR mine is correct (or both correct under different interpretations):

1. **Numerical verification at r ∈ {6, 8, 10} at small p ∈ {3, 5, 7}.** Direct computation of |T_p| at these (p, r) cells and comparison to both bounds. If |T_p| ≤ 2√p·N at r ∈ {6, 8, 10}, the chain agent is right and my Phase 3 over-counts losses. If |T_p| grows in r at fixed p, my Phase 3 is right.

2. **R79b extension.** R79b ran at p=3, r=8..20 and got |K|/√N stable ≈ 1.7-2.7 (no growth in r). **This empirical fact is more consistent with the chain agent's uniform-in-r 2√p than with my growing-in-r p^{(⌊r/2⌋−1)/2}.** Specifically: R79b at r=20, p=3 shows |K|/√N ≤ 2.65, while my bound predicts 2·3^{9/2} ≈ 280 and the chain agent's predicts 3.46. **R79b empirical fits FAR below my bound but only slightly below the chain agent's bound.**

This R79b observation is a TENSION with my disposition: **if my Phase 3 analysis were right, |K|/√N should grow significantly with r, but it doesn't.** The empirical evidence is more consistent with the chain agent's claim being correct (or even tighter than 2√p) — i.e., either some cancellation IS happening, OR the truth is empirical β=0.522 << √N constant·-bound.

**Pivot point:** the empirical evidence suggests cancellations DO exist at high r, which would mean (a), (b), or (c) above is operative. My Phase 3 derivation may be over-counting losses by failing to identify the cancellation.

---

## Refined honest disposition

> **Tentative H_CHAIN_RECON_DIFFERS_LOOSER**, but with the caveat that R79b's empirical stability at p=3, r=8..20 (|K|/√N ≤ 2.7 at all r) suggests cancellations exist between higher-digit Gauss sums that I cannot derive from the closed-form alone. The chain agent's uniform-in-r bound 2√p·√N is empirically more accurate than my growing-in-r bound.
>
> **Most likely truth:** the chain agent's claim is essentially correct (single √p loss dominates) OR the bound is even tighter (closer to strict 2√N). My derivation's growing-in-r prediction is likely an artifact of triangulating where deeper cancellation actually occurs.

This is awkward — my pre-registered hypothesis H_CHAIN_RECON_DIFFERS_LOOSER is "looser bound" but I'm flagging that the empirical evidence suggests the chain agent IS right. Let me reconcile.

**The cleanest formal disposition:**

> My INDEPENDENT MECHANICAL DERIVATION (digit-by-digit accounting via the Hensel closed-form, no consultation of CHAIN_*.md) gives a LOOSER bound than the chain agent's: `2·p^{(⌊r/2⌋−1)/2}·√N` vs `2√p·√N`. The structural reason for the looseness is that my derivation does not identify a cancellation mechanism between distinct-digit Gauss-Plancherel sums.
>
> However, R79b's empirical stability at p=3, r=8..20 suggests the chain agent's tighter bound is closer to truth — meaning either (i) a cancellation argument the chain agent uses but I missed, or (ii) the actual bilinear has additional non-rigorous-but-real cancellation beyond what either of us derives.

**Disposition for the pre-registered scheme:** H_CHAIN_RECON_DIFFERS_LOOSER on the basis of independent derivation. The chain agent's claim is likely correct, but my INDEPENDENT derivation does NOT reach it, and I cannot identify the cancellation mechanism without consulting their work.

---

## Files produced this session

1. **CHAIN_RECON_PHASE1_R3.md** — Independent r=3 c_2-collapse re-derivation. Confirms strict |T_p| ≤ 2N at r=3.
2. **CHAIN_RECON_PHASE2_R4.md** — Independent r=4 chain extension. Confirms √p loss at c_2² stratum. Direct triangulation gives 2·p^{3/2}·N (one factor p above chain agent's 2√p·N).
3. **CHAIN_RECON_PHASE3_HIGHER_R.md** — Extension to r=5, 6, 7+. Identifies n_quad(r) = ⌊r/2⌋−1 quadratic-digit strata. Under stacking assumption: bound = 2·p^{(⌊r/2⌋−1)/2}·N, growing in r.
4. **CHAIN_RECON_PHASE4_COMPARISON.md** — Comparison to chain agent's claim. Identifies the discrepancy: my analysis sees multiple stacking losses; chain agent sees single √p loss uniform in r.
5. **CHAIN_RECON_DISPOSITION.md** — This document.

---

## Adversarial checks A1-A4

(A1) **Empirical anchor R79b.** My bound consistent (R79b sits below my bound). Chain agent's bound also consistent. **R79b's stability across r=8..20 is more consistent with the chain agent's uniform-in-r bound than with my growing-in-r bound.** This is a strong signal that the chain agent's mechanism is approximately right and my Phase 3 over-counts losses.

(A2) **r=3 reduction.** ✓ Both bounds give strict 2√N at r=3. Phase 1 confirms this independently.

(A3) **Top-digit cleanness.** ✓ Top inner digit c_{r-1} always linear (since 2(r−1) > r+1 for r ≥ 3, so c_{r-1}² is at stratum p^{2(r-1)} which equals q at r=3 and exceeds q at r ≥ 4 — wait this contradicts my earlier rule).

Let me re-verify: c_k² survives iff 2k ≤ r (i.e., 2k ≤ r implies p^{2k} ≤ p^r < p^{r+1} = q so survives). At k=r-1: 2(r−1) ≤ r iff r ≤ 2. So for r ≥ 3, c_{r-1}² has 2(r−1) > r, hence c_{r-1}² is at p^{2(r-1)} which is ≥ p^{r+1} = q for r ≥ 3 — actually 2(r-1) ≥ r+1 iff r ≥ 3. So c_{r-1}² is at stratum ≥ q for r ≥ 3, hence drops.

So top digit c_{r-1} is always linear-only. ✓ Chain agent correct on this point.

(A4) **Honesty check.** I did NOT consult CHAIN_PHASE*.md, CHAIN_DISPOSITION.md, or TIGHTEN_*.md. Sources used: Hensel closed-form (HENSEL_APPROACH_A.md, doubly confirmed), r=3 base case (PATH2_BILINEAR_FROM_CLOSED_FORM.md), empirical anchor (r79b_S_partial_empirical.md). The derivation independence is preserved.

---

## Final line

> **Disposition: H_CHAIN_RECON_DIFFERS_LOOSER.** Independent re-derivation finds bound `2·p^{(⌊r/2⌋−1)/2}·√N`, growing in r. The chain agent's `2√p·√N` (uniform in r) is consistent with my analysis ONLY at r ∈ {4, 5}. At r ≥ 6, additional digit-quadratic strata appear; absent a cancellation mechanism I cannot identify from the closed-form alone, my derivation predicts the chain agent's bound is too optimistic.
>
> **HOWEVER:** R79b's empirical stability across r=8..20 at p=3 (|K|/√N bounded ≤ 2.7 with no r-growth) is more consistent with the chain agent's tighter claim than with my looser bound. This suggests the chain agent's derivation likely contains a cancellation argument that my independent reconstruction missed.
>
> The MECHANICAL disposition is H_CHAIN_RECON_DIFFERS_LOOSER (my derivation gives the looser bound). The EMPIRICAL evidence suggests the chain agent is closer to truth, with my derivation likely over-counting losses by failing to identify a high-r cancellation.
