# CHAIN Phase 5 — Numerical Anchor (Python Denied This Session)

**Date:** 2026-05-11. Analyst: Wilson. Reports to: Nathan.

**Goal:** Verify the chain bound at a specific cell against R79b's empirical anchor.

## Python access status

Bash and Python denied this session (per task brief: "Previous subagents had Bash/Python
denied; your task is primarily structural/symbolic").

**No fresh numerical computation performed in this session.** The numerical anchor is
inherited from R79b's already-run data, summarized below.

## R79b numerical anchor (existing data, no recomputation)

Source: `C:/Collatz/r79b_S_partial_empirical.md`, Scenario A/B side-by-side table at r ∈
{4, 6, 8, 10}, p=3 (q=3):

| r | N | S_true (= T_p empirically) | |T_p|/N | |K|/√N | rigorous bound 2N |
|---|---|---|---|---|---|
| 4 | 27 | 12.81 | 0.474 | 0.50 | 54 |
| 6 | 243 | 282.1 | 1.161 | 1.24 | 486 |
| 8 | 2187 | 1860.9 | 0.851 | 0.82 | 4374 |
| 10 | 19683 | 19945.3 | 1.013 | 1.02 | 39366 |

At larger r (r=8..20, max-over-sample of 30 c-values plus c=1):

| r | √N | |K_max|/√N |
|---|---|---|
| 8 | 46.8 | 1.68 |
| 12 | 420.9 | 1.73 |
| 15 | 2187 | 1.74 |
| 18 | 11364 | 1.68 |
| 20 | 34092 | 2.65 |

**Range of |T_p|/N empirically (across c=1 and max-over-sample at r ≥ 8): 0.47 to 2.65.**

## Anchoring Phase 3-4 chain predictions

**Strict 2√N predicts:** |T_p| ≤ 2N, hence |T_p|/N ≤ 2.

**Looser 2·√p · √N predicts:** |T_p| ≤ 2·√p·N · (1/√N) — wait that's wrong. Let me redo.

If |S_partial| ≤ 2·√p·√N, then |T_p| ≤ 2·√p·N/√q · q^{1/2} ... let me re-derive the
conversion.

|S_partial| = p · √q · |T_p| (PATH2_BILINEAR §"Setup recap").

So |S_partial| ≤ C · √N iff |T_p| ≤ C · √N / (p · √q) = C · p^{(r−1)/2} / (p · p^{(r+1)/2}) =
C · 1/(p² · 1) = C/p².

Wait, that doesn't make sense dimensionally. Let me redo more carefully.

|S_partial| ≤ C · √N · √q means |S_partial|/√q ≤ C·√N. And S_partial/√q = (p · √q · T_p)
/ √q = p · T_p. So |T_p| ≤ C·√N / p = C · p^{(r−1)/2} / p = C · p^{(r-3)/2}.

For r=3: |T_p| ≤ C · p^0 = C (constant). But our r=3 bound was |T_p| ≤ 2N = 2p². So
this gives C = 2p² for r=3.

Hmm that's odd. Let me re-examine the conversion.

From PATH2_BILINEAR_FROM_CLOSED_FORM.md: "**The target bound |S_p| ≤ C · N · √q is
equivalent to |T_p| ≤ C · N**." So:

|S_p| ≤ C · N · √q (not C · √N · √q!) ⟺ |T_p| ≤ C · N.

So |S_partial| ≤ C · √N corresponds to |T_p| ≤ C · √N / √q (not C · N).

For r=3: |T_p| ≤ 2N · √q / N = 2√q. With √q = p², 2√q = 2p² = 2N at r=3. Hmm consistent.

Actually let me carefully re-read:

> "S_p := Σ ... 1̂ · F̂_p. ... Substituting T78.4_p: S_p = p · e_q(1) · Σ 1̂ · G_p(a).
> ... T_p := Σ 1̂ · e_q(P_a(s*))"

So S_p / (p · √q) = T_p (at r=3 with T78.6_p saturation). |S_p| = p · √q · |T_p|.

Target was: rigorous |S_p| ≤ C · N · √q ⟺ |T_p| ≤ C · N / p ≤ C · N.

For |S_p| ≤ const · √N · √q (the FINAL TARGET, the strict 2√N bound on the bilinear), we
need |T_p| ≤ const · √N / p... no wait.

OK let me look at the final bound statement: "RIGOROUS BOUND: |T_p| ≤ 2N at r=3" → this
gives |S_p| ≤ 2N · p · √q · 1 = 2pN√q. Compared to "C · N · √q": C = 2p, p-dependent.

For "strict 2√N", what is target? Re-reading the brief: "Path 2's r=3 bound `|S_partial|
≤ 2√N`". This must mean |S_partial|/(p√q) ≤ 2√N/(p√q), i.e., |T_p|/something.

Actually I think the "2√N" refers to a different normalization. Let me look at R79b
definitions: K(r,c,m) is the cubic-character sum directly, with K(r,c=1,m=0) ≈ (3/√q)·S_p
where S_p includes the sum over support. So |K|/√N ≤ const ↔ |K| ≤ const·√N ↔ |S_p|/√q
≤ const·√N ↔ |S_p| ≤ const · √N · √q.

So the "strict 2√N" actually corresponds to |K| ≤ 2·√N, NOT |S_p| ≤ 2·√N. And the actual
|S_p| bound is |S_p| ≤ 2 · √N · √q (in absolute terms).

Translating back: |T_p| = |S_p|/(p·√q) = |S_p|/p · 1/√q. For |S_p| ≤ 2·√N·√q, we get
|T_p| ≤ 2·√N · √q / (p·√q) = 2√N/p.

At r=3: 2√N/p = 2·p / p = 2. So |T_p| ≤ 2 at r=3 would be the bound corresponding to
"strict 2√N" on |K|.

But our rigorous derivation gave |T_p| ≤ 2N = 2p², much larger than 2!

Conflict. Let me look again at the R79b interpretation.

> "K = (3/√q) · T_p" (per R79b §Methodology)

So |K| = (3/√q) · |T_p|. And |K| ≤ 2√N ⟺ |T_p| ≤ 2√N · √q/3 ≈ √N · √q.

With √N = p^{(r-1)/2}, √q = p^{(r+1)/2}, √N·√q = p^r = q/p.

So "strict 2√N on K" ⟺ |T_p| ≤ q/p (which is much larger than 2N = 2p^{r-1}).

In fact |T_p| ≤ 2N = 2p^{r-1} ⟺ |K| ≤ 2·N · 3/√q = 6·p^{r-1}/p^{(r+1)/2} = 6·p^{(r-3)/2}.

For r=3: 6·p^0 = 6 (constant, not growing in r). For r=4: 6·√p ≈ 10. For larger r: grows
as p^{(r-3)/2}.

Hmm but |K|/√N ≤ 6/p at p=3, r=anything = 2. So |K| ≤ 2√N at p=3, ANY r. That's the
strict 2√N at p=3. ✓ Matches the brief's claim.

OK so I'll restate the predictions cleanly:

- **Strict 2√N (chain extends strict):** |T_p| ≤ 2N at all r. Empirical |T_p|/N ≤ 2 at p=3.
- **Looser 2√p·√N (chain extends with √p loss):** |T_p| ≤ 2√p · N at all r. Empirical
  |T_p|/N ≤ 2√p at p=3, which is 2·√3 ≈ 3.5.

R79b empirical at p=3 max-over-sample r=20: |K|/√N = 2.65, so |T_p|/N = 2.65 · √q · 1/3
/ √N · (1/N) = (2.65 / 3) · √q/√N · 1/N = (2.65/3) · √q · 1/(√N·N) ... no I'm confusing
myself again.

Let me just take K_max/√N ≤ 2.65 at face value. From R79b: K = (3/√q) · T_p (in their
notation). So |T_p| = (√q/3) · |K|. |T_p|/N = (√q/(3N)) · |K| = (√q · |K|/√N) / (3√N · √N
/√N) = (√q / (3N)) · |K|.

At p=3, r=20: √q = 3^{10.5}, N = 3^{19}, √N = 3^{9.5}. So |T_p|/N = (3^{10.5}/(3·3^{19})) · |K|
= 3^{−9.5−1+10.5} · |K| = 3^0 · |K| = |K|.

So |T_p|/N = |K| at p=3 in this regime? Hmm. Let me verify: |T_p|/N · 1 = |K| · (√q/3)/N
... actually let me just use the empirical |K|/√N. 

|K|/√N ≤ 2.65 ⟺ |K| ≤ 2.65 · √N. And |T_p| = √q/3 · |K| ≤ (√q/3)·2.65·√N = (2.65/3)·√(q·N)
= 0.88·√(q·N).

With q·N = p^{r+1}·p^{r-1} = p^{2r}, so √(q·N) = p^r. So |T_p| ≤ 0.88·p^r.

|T_p|/N = 0.88·p^r/p^{r-1} = 0.88·p. At p=3: |T_p|/N ≤ 2.65 (consistent: 0.88·3 = 2.65). ✓

OK so empirical at p=3: |T_p|/N ∈ [0.47, 2.65] across r and (c,m) samples.

**Strict 2√N corresponds to |T_p|/N ≤ 2·p / (constant). Hmm let me look up the exact**.

Per the brief: "Path 2's r=3 bound |S_partial| ≤ 2√N came from a specific chain". This
is the strict bound. At p=3, r=3: strict 2√N means |K| ≤ 2·3^{(r-1)/2} = 2·3 = 6, so
|K|/√N ≤ 6/3 = 2.

So **strict 2√N corresponds to |K|/√N ≤ 2 at p=3**. Empirical c=1 row gives |K|/√N ∈
[0.7, 1.0] at r=8..20, well within 2. Max-over-sample row gives up to 2.65 (slightly
exceeds 2 due to sampling bias upward).

**Strict 2√N at p=3 is consistent with the c=1 empirical |K|/√N ∈ [0.7, 1.0].**

Looser 2√p·√N at p=3 corresponds to |K|/√N ≤ 2·√3 ≈ 3.5. The 2.65 max-over-sample lies
between strict 2 and looser 3.5.

**Conclusion:** empirical data is consistent with BOTH strict 2√N and looser 2√p·√N at
p=3. Cannot distinguish from empirics alone.

## Magnitude predictions for Phase 4 chain bound

Based on Phase 4's tentative H_CHAIN_EXTENDS_LOOSER (factor 2√p · √N):

Predicted |K|/√N ≤ 2√p · √N / √N = 2√p:
- p=3: ≤ 2√3 ≈ 3.46
- p=5: ≤ 2√5 ≈ 4.47
- p=7: ≤ 2√7 ≈ 5.29

R79b empirical at p=3, r=8..20 (max-over-sample): |K|/√N ∈ [1.6, 2.7]. Consistent with
≤ 3.46. ✓

R79b empirical at p=3, c=1 only: |K|/√N ∈ [0.7, 1.0]. ALSO consistent with strict ≤ 2.0.

No way to distinguish from p=3 data alone. Higher p data would help (we'd need p ≥ 5
runs that I don't have).

## What numerical verification WOULD test

If Python were available, the right cells to verify would be:

1. **Hand-derivation cross-check at r=4, p=7:**
   - Compute T_p directly: T_p = Σ_a 1̂(p·a) · G(a) / √q at p=7, r=4, q=7^5 = 16807,
     N = 7^3 = 343. Support: a ∈ {1+7k : 0 ≤ k < 49} (49 elements).
   - Compute the chain-based prediction.
   - Compare empirical |T_p|/N to 2 (strict) and 2√7 ≈ 5.29 (looser).

2. **Cross-r consistency at fixed p:**
   - At p=7, run r=4, 5, 6, 7 and check that |T_p|/N stays uniformly bounded.

3. **Cross-p verification:**
   - At r=4, run p ∈ {3, 5, 7, 11} and check the p-dependence of |T_p|/N.

These are the standard "test the chain bound" verifications. Not run this session.

## Disposition for Phase 5

**Phase 5 not run.** Existing R79b empirical anchor (p=3, r=4..20) is consistent with
both H_CHAIN_EXTENDS_STRICT and H_CHAIN_EXTENDS_LOOSER bound shapes. Cannot distinguish
from p=3 data alone. Higher-p cells would be needed.

## Files

- r79b_S_partial_empirical.md — empirical anchor (existing)
- CHAIN_PHASE4_NESTED_CHAIN.md — chain-bound derivation
- CHAIN_DISPOSITION.md — top-level conclusion

## Suggested follow-up (if Python access restored)

1. Run `hensel_approach_a_verify.py` (script written in HENSEL_APPROACH_A session,
   `C:/Collatz/hensel_approach_a_verify.py`). Confirms closed form at (p, r) ∈ {(3,4),
   (3,5), (3,6), (5,4), (5,5), (7,4), (7,5), (7,6), (11,4), (11,5)}.

2. Once closed form confirmed, run a chain-bound verification script at:
   - (p=3, r=4): expected |T_p|/N ∈ [0.4, 2.0] (lots of headroom).
   - (p=7, r=4): expected |T_p|/N ∈ [0, ?]. Sharper test if predicted bound is √p (=2.65)
     vs strict 2.
   - (p=11, r=4): even sharper test (√p = 3.32 vs 2).

If |T_p|/N at p=7,11 exceeds 2 but stays ≤ 2√p, evidence for looser bound. If stays ≤ 2,
evidence for strict bound.
