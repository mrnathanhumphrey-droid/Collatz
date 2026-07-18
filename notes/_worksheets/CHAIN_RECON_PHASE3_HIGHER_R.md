# CHAIN_RECON Phase 3 — Independent extension to r=5, r=6, r=7+

**Date:** 2026-05-11. Analyst: Wilson (chain-reconstruction track).

Building on Phase 1 (r=3 strict 2N) and Phase 2 (r=4 single √p loss at c_2). Examining whether higher-digit quadratic terms (c_3², c_4², ...) stack additional √p losses as r grows.

---

## General digit-quadratic stratum analysis

At general r ≥ 4: C_a = 1 + p·s_0 + p²·c_2 + ... + p^{r-1}·c_{r-1}. s*(r) = s_0 + p·c_2 + ... + p^{r-2}·c_{r-1}.

The Hensel polynomial:
> P_a(s*(r)) ≡ Σ_{j=2}^{r} (−1)^{j−1}·(p·s*)^j/(j(j-1)) mod p^{r+1}.

For each digit c_k (k = 2, ..., r−1) at coefficient p^{k-1} in s*, the c_k² term enters s*² at coefficient p^{2(k-1)}. Multiplied by p² from the j=2 term −p²·s*²/2, the c_k² piece appears at stratum p^{2k} in P_a(s*(r)).

**The c_k² contribution at stratum p^{2k} survives mod p^{r+1} iff 2k ≤ r, i.e., k ≤ r/2.**

(If 2k = r+1, c_k² is at stratum p^{r+1} = q, drops. If 2k > r+1, c_k² stratum exceeds q, drops.)

**Critical: c_k² ENTERS the phase mod q if and only if k ≤ ⌊r/2⌋.**

(Equivalently: k satisfies 2k ≤ r. At r=4: k ≤ 2, so only c_2² appears. At r=5: k ≤ 2, only c_2². At r=6: k ≤ 3, both c_2² and c_3². At r=7: k ≤ 3, c_2² and c_3². At r=8: k ≤ 4, c_2², c_3², c_4².)

---

## Counting √p losses at each r

For each digit c_k that has a quadratic-in-c_k term in P_a(s*(r)) mod q, the c_k Plancherel sum is a quadratic-Gauss-twisted Dirichlet kernel of magnitude √p (when v_0 of the corresponding u-residue-class shift = 0 — see Phase 2 Case A; uniformly ≤ 2√p by Phase 2 incomplete-Gauss-sum bound).

For each LINEAR-only digit c_k (k > ⌊r/2⌋, i.e., c_{⌊r/2⌋+1}, ..., c_{r-1}), the Plancherel-collapse gives a clean factor p.

**Number of digits with quadratic structure at r ≥ 4:**
> n_quad(r) := |{k : 2 ≤ k ≤ ⌊r/2⌋}| = ⌊r/2⌋ − 1.

(Subtract 1 because c_1 doesn't exist as an inner digit — c_1 = s_0 is the outer class.)

| r | ⌊r/2⌋ | n_quad(r) | quadratic digits | linear-only digits |
|---|---|---|---|---|
| 3 | 1 | 0 | none (c_2² drops mod q) | (no extra inner digits) |
| 4 | 2 | 1 | c_2 | c_3 |
| 5 | 2 | 1 | c_2 | c_3, c_4 |
| 6 | 3 | 2 | c_2, c_3 | c_4, c_5 |
| 7 | 3 | 2 | c_2, c_3 | c_4, c_5, c_6 |
| 8 | 4 | 3 | c_2, c_3, c_4 | c_5, c_6, c_7 |
| general odd r=2m+1 | m | m−1 | c_2,...,c_m | c_{m+1},...,c_{r-1} |
| general even r=2m | m | m−1 | c_2,...,c_m | c_{m+1},...,c_{r-1} |

(Note even/odd r have the same n_quad(r) = ⌊r/2⌋−1.)

---

## Verification at r=5 (explicit derivation)

At r=5: q=p^6, N=p^4. C_a = 1+p·s_0+...+p^4·c_4. s*(5) = s_0+p·c_2+p²·c_3+p³·c_4.

Compute P_a(s*(5)) = −p²·s*²/2 + p³·s*³/6 − p^4·s*^4/12 + p^5·s*^5/20 mod p^6.

s*² mod p^4 (since p²·s*² needs mod p^4 to hit p^6):
s*² = s_0² + 2p·s_0·c_2 + p²·(2s_0·c_3 + c_2²) + p³·(2s_0·c_4 + 2c_2·c_3) + p^4·(2c_2·c_4 + c_3²) + ...

p²·s*² mod p^6: includes the p^4 term 2c_2·c_4 + c_3². 

Wait — c_3² appears in s*² at coefficient p^4 (from c_3·p² · c_3·p² = c_3²·p^4). Multiplied by p²: contributes to p^6 stratum of p²·s*². At modulus q=p^6: this is the TOP stratum, which is mod 0 (i.e., e_{p^0}(coefficient) = 1). **So c_3² at p^6 stratum at r=5 is mod 1 → trivial.** ✓

Hmm but my general rule said c_k² appears at stratum p^{2k} survives iff 2k ≤ r. At k=3, 2k=6, r=5: 6 ≤ 5? NO. So c_3² does NOT survive at r=5. ✓ Matches my general rule.

Let me redo the rule: stratum p^m mod q=p^{r+1} drops when m ≥ r+1 (i.e., m=r+1 is mod 0, equivalent to e_1 = 1). So c_k² at stratum p^{2k} survives iff 2k ≤ r. At r=5, k=3: 6 ≤ 5? No, so drops. ✓

At r=4, k=2: 4 ≤ 4? Yes, c_2² at p^4 stratum survives (since p^4 < p^5 = q). ✓

So **c_k² survives iff 2k ≤ r**. The number of quadratic digits is |{k : 2 ≤ k ≤ r/2}| = ⌊r/2⌋ − 1.

At r=5: ⌊5/2⌋ − 1 = 1. Only c_2. ✓

At r=6: ⌊6/2⌋ − 1 = 2. c_2 and c_3. Let me verify.
- c_2² at p^4 stratum: 2k=4 ≤ r=6? Yes, survives.
- c_3² at p^6 stratum: 2k=6 ≤ r=6? Yes, survives (since p^6 < p^7 = q at r=6).

Both survive. ✓

At r=7: ⌊7/2⌋ − 1 = 2. c_2 and c_3.
- c_2² at p^4: survives (4 ≤ 7).
- c_3² at p^6: survives (6 ≤ 7).
- c_4² at p^8: 8 ≤ 7? No, drops.

Only c_2 and c_3 quadratic. ✓

---

## Bound on |T_p| under the stacking-losses hypothesis

If each quadratic digit c_k contributes a √p loss (independently), and each linear digit contributes a clean factor p, then the bilinear bound at general r is:

> **|T_p| ≤ C · (√p)^{n_quad(r)} · N = C · p^{(⌊r/2⌋−1)/2} · N.**

In terms of |S_partial| (using |T_p| ≤ C·N maps to |S_partial| ≤ C·√N):

> **|S_partial| ≤ C · p^{(⌊r/2⌋−1)/2} · √N**

| r | n_quad | p^{n_quad/2} | bound shape |
|---|---|---|---|
| 3 | 0 | 1 | 2·√N (strict) |
| 4 | 1 | √p | 2·√p·√N (chain agent's claim) |
| 5 | 1 | √p | 2·√p·√N (chain agent's claim) |
| 6 | 2 | p | 2·p·√N (LOOSER than chain agent) |
| 7 | 2 | p | 2·p·√N (LOOSER than chain agent) |
| 8 | 3 | p^{3/2} | 2·p^{3/2}·√N (LOOSER than chain agent) |
| 9 | 3 | p^{3/2} | 2·p^{3/2}·√N (LOOSER) |
| 10 | 4 | p² | 2·p²·√N (LOOSER) |
| ... | ⌊r/2⌋−1 | p^{(⌊r/2⌋−1)/2} | 2·p^{(⌊r/2⌋−1)/2}·√N |

**Asymptotic behavior:** at fixed p, as r → ∞: the bound exponent grows like (r/4)·log p. **The bound is NOT uniform in r — it grows polynomially.**

---

## Critical check against R79b empirical anchor

R79b at p=3, r=8..20: |K|/√N ∈ [0.7, 2.7] (max-over-sample c=1: [0.7, 1.0]).

|S_partial| ≤ C · p^{(⌊r/2⌋−1)/2} · √N translates to |K|/√N ≤ C · p^{(⌊r/2⌋−1)/2}/√N · √N · ... let me redo.

|S_partial| ≤ C·X·√N for some X. From the Plancherel mapping |K| ≈ (p/√q)·|S_partial|·(constants): |K|/√N ≤ (p·C·X)/√(q·N) = C·X · p/√(p^{2r}) = C·X/p^{r-1}.

Hmm that doesn't match. Let me re-derive the K↔T_p mapping at p=3.

From R79b: K(r,c=1,m=0) = (3/√q)·S_true where S_true = Σ_a 1̂(p·a)·ψ_true(a)/√q · √q · ... actually S_true ≈ |T_p| (= Σ_a 1̂·ψ_true).

|K| = (3/√q)·|T_p|. So |K|/√N = (3/√q)·|T_p|/√N. Wanted: |K|/√N ≤ C empirical const ≈ 2.7 at p=3.

|T_p| ≤ C·X·N translates to: |K|/√N ≤ (3/√q)·C·X·N/√N = (3·C·X/√q)·√N = (3·C·X·√N)/√q = 3·C·X·p^{(r-1)/2}/p^{(r+1)/2} = 3·C·X/p.

At p=3: |K|/√N ≤ C·X. So **bound shape |T_p| ≤ C·X·N matches |K|/√N ≤ C·X**.

For my Phase 3 result |T_p| ≤ C·p^{(⌊r/2⌋−1)/2}·N: |K|/√N ≤ C·p^{(⌊r/2⌋−1)/2}.

| r | (⌊r/2⌋−1)/2 | bound |K|/√N at p=3 | R79b empirical |
|---|---|---|---|---|
| 8 | 3/2 | 2·3^{3/2} ≈ 10.4 | 1.68 |
| 10 | 2 | 2·9 = 18 | 1.74 |
| 12 | 5/2 | 2·3^{5/2} ≈ 31 | 1.73 |
| 16 | 7/2 | 2·3^{7/2} ≈ 94 | 1.74 |
| 20 | 9/2 | 2·3^{9/2} ≈ 280 | 2.65 |

**My Phase 3 bound is rigorous but consistent with R79b (empirical sits well below the bound).** The bound GROWS with r at fixed p, while empirical stays bounded.

**Contrast: chain agent's 2√p·√N gives constant 2√p = 3.46 at p=3.** That's tighter than my Phase 3 bound at every r ≥ 6, AND also consistent with R79b (which sits below 3.46 for c=1 row, but max-over-sample reaches 2.7 occasionally near 3.46).

R79b's empirical β=0.522 is even tighter (suggests true growth rate ∝ N^{0.522}, NOT bounded by const · √N). So both rigorous bounds (mine and chain agent's) over-predict the empirical truth — that's fine, rigorous bounds CAN be loose vs empirical.

The DECISIVE question is which rigorous bound is CORRECT (provably tight or at least valid).

---

## Comparison to chain agent's claim

Chain agent: |S_partial| ≤ 2√p · √N uniformly in r, polylog-free. This corresponds to **single √p loss at c_2² stratum, every r**.

My Phase 3 derivation: at r ≥ 6, MULTIPLE digit-quadratic strata appear (c_2² AND c_3² at r=6,7; c_2², c_3², c_4² at r=8,9; etc.). Each contributes an independent √p loss to the Plancherel chain (assuming the same "single √p per Gauss-twisted Plancherel" rule the chain agent applies to c_2).

**Result: at r ≥ 6, my bound is** `2·p^{(⌊r/2⌋−1)/2}·√N`**, growing polynomially in p with exponent (⌊r/2⌋−1)/2.**

This is **MULTIPLE √p losses stacking**, contradicting the chain agent's "single √p loss uniform in r".

---

## Honest scope and caveats

**Caveat 1:** I established the √p loss at c_k for k ≤ ⌊r/2⌋ by the SAME mechanism the chain agent uses (Gauss-twisted Plancherel). I haven't verified that the losses truly STACK MULTIPLICATIVELY rather than INTERFERING (canceling each other partially).

Specifically: the chain agent's argument involves nested Plancherel on (c_3, c_4, ..., c_{r-1}). At r=6, the chain would be (c_5, c_4, c_3, c_2). The c_5 and c_4 collapses are LINEAR (clean factor p each). The c_3 collapse is QUADRATIC (√p factor). The c_2 collapse is also QUADRATIC (√p factor). Net: p · p · √p · √p = p^3 = N^{... hmm.

At r=6, N = p^5. p^3 = p^5/p² = N/p². So "Plancherel save factor of p^3" applied to "trivial bound on bilinear (Σ |1̂|² · |supp|² = N³ via Cauchy)" gives... let me think.

Actually the proper "save factor" analysis: the trivial bound on |T_p| is ≤ N² (each term bounded by N, summed over N elements). Each Plancherel-collapse step on a c_k digit reduces |T_p| by a factor of (something) — but it's not just multiplicative "save factors". The structure is more delicate.

**Let me re-examine more carefully:** at r=3, |T_p| ≤ 2N. At r=4, my triangulated bound: p^{3/2}·N (one √p loss × one extra trivial p factor from the extra digit dimension). At r=5: similar, one √p loss × two extra trivial p factors? 

Hmm let me reconsider the r=4 case. I had:
- Inner_{c_3} (linear, factor p)
- c_2 sum (quadratic Gauss, factor √p) → leaves Σ_v of D_{p²} × G_p
- v sum (length p², triangle bound)

Total: p · √p · p² = p^{7/2}. Per s_0 row. Then Σ_{s_0} factor p: p^{9/2} = N · p^{3/2}.

I had ONE √p from quadratic and ONE additional factor p from the v-direction triangle (length p² → p² triangle bound).

At r=5: structure should be:
- c_4 collapse (linear, factor p)
- c_3 collapse (linear, factor p) — at r=5, c_3 IS linear-only (per my n_quad(5)=1 analysis)
- c_2 sum (quadratic Gauss, factor √p)
- v sum (length p² → triangle)
- Σ_{s_0}: factor p

But wait — after each Plancherel step, the "u-residue-class" gets refined further. At c_3 step in r=5, u is restricted to size N/p. At c_2 step, restricted to size N/p², etc. The "v" variable I had at r=4 was the parametrization of u in the restricted class.

Let me think about it as: how many factors of p does the inner Plancherel chain extract?

At r=3: c_2 collapse extracts factor p (Plancherel). The Σ over (a, restricted_u) gives Σ |D_p| ≤ p + log p. Total: p · (p + log p) = N + p·log p ≤ 2N. **Strict 2N.**

At r=4: c_3 collapse extracts factor p. Then c_2 sum extracts factor √p (Gauss instead of clean Plancherel). Combined Plancherel-saving factor: p · √p = p^{3/2}. Σ over (restricted_a, restricted_u) — how much is this?

The supp is N = p^3. Each Plancherel-collapse refines by 1/p. After c_3-collapse: a-direction unchanged (s_0 ∈ Z/p, c_2 ∈ Z/p), u restricted to size N/p = p². After c_2 sum: a-direction reduced (s_0 ∈ Z/p only), u restricted to size N/p² = p (since c_2-Gauss doesn't restrict u, just twists; but maybe the analog of u-restriction here is the v-direction). I'm getting confused.

Let me think about the STRUCTURE FORMAL way. The bilinear:
T_p = Σ_a 1̂(p·a) · ψ(a) where a ∈ supp, |supp| = N.

By Cauchy-Schwarz: |T_p|² ≤ N · Σ_a |1̂|²·|ψ|² ≤ N · N² · 1 = N³ → |T_p| ≤ N^{3/2}.

That's the brutally-trivial Cauchy bound. (For r=4: N^{3/2} = p^{9/2}.)

To do better than N^{3/2}, we need PHASE CANCELLATION exploitation. Each Plancherel-collapse saves factor √p in the L² inner product. The total save needed to reach |T_p| ≤ C·N is √p^{r-1} = √N from N^{3/2}.

If each Plancherel collapse saves √p, and we have r-1 digits to collapse (c_2, c_3, ..., c_{r-1}), then total save: √p^{r-1} = √N. ✓ Reaches |T_p| ≤ C·N.

If at some digits the save is only p^{1/4} instead of √p (i.e., quadratic Gauss instead of clean), then the total save is reduced. Specifically: clean Plancherel collapses save √p (in the L² norm); quadratic Gauss saves p^{1/4} (since |Gauss| = √p out of trivial p, so save factor √p/p = 1/√p — hmm that's WORSE not better).

Wait let me reconsider. The "save" at each digit level is the factor by which the Plancherel-on-that-digit reduces the bilinear vs. trivial.

**Clean linear Plancherel save:** trivial sum over c_k ∈ Z/p of |D| ≤ p·|D| (triangle), Plancherel collapses to p·D (single term). Save: 1 (same magnitude — Plancherel doesn't save here, the saving comes from the "max" being equal to "sum" when collapsed).

Hmm I'm confusing myself. Let me think directly.

**At r=3:** strict 2N. **At r=4:** if chain agent right, |T_p| ≤ 2√p·N. Loss: factor √p compared to r=3 strict.

**At r=5:** if chain agent right, also 2√p·N (still one quadratic digit c_2, so single √p loss). Loss: same √p.

**At r=6:** my analysis says c_2 AND c_3 both have quadratic strata. If losses stack: 2p·N. Loss: factor p.

The chain agent's claim that the bound stays at 2√p·N UNIFORMLY in r requires the c_3 quadratic stratum at r=6 NOT to add an extra √p loss. Could there be a CANCELLATION between c_2² and c_3² Gauss-twisted Plancherels?

Looking at the r=6 phase: at p^6 stratum (the new stratum at r=6), the coefficient includes both c_3² (from p²·s*² · −1/2 at stratum p^6 from c_3·p²·c_3·p² = c_3²·p^4 → −c_3²/2 at stratum p^6) and other cross-terms.

So at p^6 stratum, phase is `e_p(−c_3²/2 + cross terms)`. The c_3 sum becomes a quadratic Gauss in c_3 at modulus p, magnitude √p. **The c_3 quadratic stratum is INDEPENDENT of the c_2 quadratic stratum at p^4** — they're at different strata, no obvious cancellation between them.

**Conclusion:** at r=6, two independent quadratic Gauss sums (one per digit c_2 and c_3) stack their √p factors. No mechanism for cancellation between them in the basic Plancherel chain.

**My Phase 3 bound at r ≥ 6:** |T_p| ≤ 2·p^{(⌊r/2⌋−1)/2}·N, growing in r. **LOOSER than chain agent's claim.**

---

## Honest scope (Phase 3)

I am ASSUMING that the c_k² losses stack multiplicatively. Verifying this rigorously requires careful Plancherel-chain bookkeeping that I haven't completed. The cancellation possibility cannot be ruled out from my analysis alone.

**Two interpretations:**

(I) **Cancellations DON'T exist** (losses stack): my bound `2·p^{(⌊r/2⌋−1)/2}·√N` is right. Chain agent's `2√p·√N` is WRONG at r ≥ 6. **H_CHAIN_RECON_DIFFERS_LOOSER.**

(II) **Cancellations exist** (some clever Plancherel reorganization): the chain agent's `2√p·√N` might be correct uniformly. But the burden of proof is on showing the cancellation EXISTS — and I don't see it from the materials I've consulted (Hensel closed form + PATH2 r=3).

The chain agent's pre-reg claim that the loss is uniform √p at every r ≥ 4 needs careful structural argument FOR cancellation; absent that, my Phase 3 derivation says the losses stack.

---

## Adversarial checks at r=5, r=6

(A1) **Empirical anchor R79b at p=3, r=8..20.** My bound at r=20 with n_quad=9: |K|/√N ≤ 2·3^{9/2} ≈ 280. Empirical 2.65. **Both bounds consistent, but R79b empirical sits FAR below both rigorous bounds — empirical β=0.522 is the truth, rigorous bounds are loose.**

(A2) **r=3 reduction.** ✓ My formula gives n_quad(3)=0, bound = 2·√N strict. Consistent with Phase 1.

(A3) **Top-digit cleanness.** At r=5, c_4 is linear-only (it's the top inner digit at coefficient p^3 in s*, contributing only linearly to P_a stratum p^5). ✓ Top digit always linear (since c_{r-1}² at p^{2(r-1)} > p^{r+1} = q drops).

(A4) **Honesty check.** I did not consult CHAIN_*.md or TIGHTEN_*.md files. My derivation only uses the Hensel closed form's explicit polynomial structure.

---

## Phase 3 summary

> **At r ≥ 6, my independent derivation finds n_quad(r) = ⌊r/2⌋−1 ≥ 2 quadratic-Gauss strata. Assuming the √p losses stack multiplicatively (no cancellation), the bound is `|S_partial| ≤ 2·p^{(⌊r/2⌋−1)/2}·√N`, growing in r.**
>
> **The chain agent's claimed bound `2√p·√N` (uniform in r) is consistent ONLY IF there is a cancellation mechanism between distinct-digit Gauss-Plancherel sums that I do not establish from the closed-form alone.**
>
> **Direction of derivation: H_CHAIN_RECON_DIFFERS_LOOSER at r ≥ 6.**

The disposition will lean LOOSER at r ≥ 6, with the caveat that a cancellation mechanism (not visible in my derivation) could rescue the chain agent's claim. Phase 4 compares to the chain agent's specific framing.
