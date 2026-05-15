# W2.B — Monotone partition combinatorics: single-2-block patterns

**Date:** 2026-05-14
**Task:** Track A wrinkle 2, step 2. Count monotone partitions on [n] with
exactly one 2-block and (n−2) singletons, verbatim from Hasebe monograph
Defn 3.23 + Prop 3.25 + Thm 3.26.
**Mode E:** verbatim citations.

---

## 1. Verbatim Hasebe monograph Defn 3.23 (p. 39, from Deliverable A §3)

> "Let T be a totally ordered finite set. An ordered set partition `π = (ρ, ≤)`
> of T is called a **monotone set partition** if
> - ρ is a noncrossing set partition,
> - if B, B' ∈ ρ satisfies B ⪯ B' then B ≤ B'.
>
> The set of monotone set partitions of T is denoted by `M(T)`."

Here `B ⪯ B'` is the covering relation: B covers B' if min B ≤ i ≤ max B
for all i ∈ B'. So **inner blocks come higher in linear order than outer
blocks**.

## 2. Verbatim Hasebe monograph Prop 3.25 (p. 39, from Deliverable A §3)

> "The cardinality of M(T) is `(|T| + 1)! / 2`."

Check: `|M(1)| = 1`, `|M(2)| = 3`, `|M(3)| = 12`, `|M(4)| = 60`. ✓

## 3. Verbatim Hasebe monograph Thm 3.26 / HS 2011 Thm 6.1 (verbatim Deliverable A §3)

> "On any nc-probability space (A, φ) and for any x ∈ A, we have
> `φ(x^n) = Σ_{π ∈ M(n)} (1/|π|!) κ_π(x),  n ∈ N.`  (3.13)"

where `κ_π := κ_{|B_1|} κ_{|B_2|} ⋯ κ_{|B_k|}` and `|π|` = number of blocks.

---

## 4. Counting one-2-block, (n−2)-singleton monotone partitions on [n]

**Step 1: Underlying set partition `ρ`.**

For `ρ` to be **noncrossing** with one 2-block + (n−2) singletons, the 2-block
must be an **interval** `{i, i+1}` for some `i ∈ {1, ..., n−1}`. (If the 2-block
were `{i, j}` with `j > i+1`, any singleton `{k}` with `i < k < j` would cross
it — but singletons don't actually cross in the strict noncrossing-partition
sense; however, the noncrossing constraint for arbitrary 2-blocks `{i, j}`
allows non-interval blocks. Let me reconsider.)

**Correction:** A 2-block `{i, j}` with `j > i+1` is **noncrossing** with all
singletons (singletons trivially don't cross anything). So **any** 2-block
`{i, j}` with `1 ≤ i < j ≤ n` qualifies as the underlying noncrossing
partition. Position count: `C(n, 2) = n(n−1)/2`.

**Hmm — but in Deliverable C §4 the count was given as (n−1). Let me check
against Hasebe Defn 3.23 once more.**

Defn 3.23 first bullet: ρ is noncrossing. A 2-block `{i, j}` with the
remaining elements as singletons is **always** noncrossing (noncrossing is
about pairs of blocks intersecting in a "crossing" way; a single 2-block plus
singletons has at most one nontrivial block and singletons can't cross it).

So the **noncrossing constraint alone** gives `C(n, 2) = n(n−1)/2` choices.

**Step 2: Apply monotone ordering (second bullet of Defn 3.23).**

For each underlying ρ, we then order the blocks linearly subject to: if
`B ⪯ B'` (B covers B' in the nesting sense), then `B ≤ B'` in the linear
order.

For our setup with one 2-block `{i, j}` and (n−2) singletons:
- The 2-block `{i, j}` covers each singleton `{k}` with `i ≤ k ≤ j`, which is
  `j − i + 1 − 2 = j − i − 1` singletons (the ones strictly between i and j).
- The remaining `(n−2) − (j − i − 1) = n − j + i − 1` singletons are outside
  the 2-block; the 2-block does not cover them.
- Within "covered" singletons, they are pairwise incomparable in ⪯ (singletons
  cover no one).
- Within "outside" singletons, they are pairwise incomparable in ⪯.
- The 2-block IS COVERED BY... nothing (it's not nested inside any other
  block in this configuration). So the 2-block is at the top of the ⪯ partial
  order.

**Wait — the convention in Defn 3.23 says inner blocks come HIGHER in the
linear order than outer blocks.** So if 2-block `{i, j}` is "outer" and the
singletons it covers are "inner"... no, the singletons covered ARE inner.

Re-reading carefully: "if B ⪯ B' then B ≤ B'". B ⪯ B' means B covers B'
(B is bigger / outer). The condition `B ≤ B'` means outer comes before
inner in the linear order. So the 2-block (outer) comes BEFORE the
singletons it covers (inner) in the linear ordering.

Number of valid linear orderings of `{2-block, n−2 singletons}` subject to:
- 2-block ≤ each of its (j − i − 1) covered singletons.

For a 2-block at position {i, j}, with `k := j − i − 1` covered singletons and
`m := n − 2 − k` outside singletons:
- Total elements to order: n − 1 (1 block + n−2 singletons)
- Constraint: 2-block must come before its k covered singletons in the ordering
- The outside m singletons + 2-block + k covered singletons: total n−1
- Number of linear orderings = (n−1)! / (k+1) [the 2-block must be the
  FIRST among {2-block, k covered singletons}, which has (k+1) elements;
  probability 1/(k+1)]

**Total count of monotone partitions with this fixed underlying ρ:**

`(n − 1)! / (k + 1) = (n − 1)! / (j − i)`

**Step 3: Sum over choices of 2-block {i, j}.**

For each `j − i = d ∈ {1, 2, ..., n−1}`, the number of 2-blocks with that
gap is `(n − d)`. Total monotone partitions M(n) with one 2-block:

`Σ_{d=1}^{n−1} (n − d) · (n − 1)! / d  =  (n−1)! · Σ_{d=1}^{n−1} (n − d)/d`
                                       =  (n−1)! · [n · Σ_{d=1}^{n−1} 1/d − (n−1)]
                                       =  (n−1)! · [n · H_{n−1} − (n−1)]

where `H_{n−1}` is the (n−1)-th harmonic number.

**Sanity check at n = 2:**
- Only 2-block is {1,2}, k=0 covered singletons.
- Count = 1 · (2-1)!/1 = 1.
- Plus the all-singletons partition (no 2-block) = 1.
- M(2) total = 1 + 1 = 2.
- But Prop 3.25: |M(2)| = (2+1)!/2 = 3.

**Discrepancy.** Let me redo the M(2) count by hand. T = {1, 2}. Noncrossing
partitions of {1,2}: (a) `{1,2}` as one 2-block, or (b) `{{1},{2}}` as two
singletons. For (a), the linear order on the single block is trivial: 1 way.
For (b), the linear order on `{{1}, {2}}` has 2! = 2 orderings, both
trivially monotone (singletons cover nothing, so no constraint). Total
M(2) = 1 + 2 = 3. ✓

**So one 2-block contributes 1 to M(2), and the (n−2)=0 singletons-only case
contributes 0 (vacuous). My formula's `n=2` case `(n−1)!·[n·H_{n−1} − (n−1)]`
= `1·[2·1 − 1] = 1`. ✓ matches.**

**Sanity check at n = 3:**
- 2-blocks: {1,2}, {1,3}, {2,3}.
- For {1,2}: covered singletons = 0, outside = 1 ({3}). Count = (3-1)!/1 = 2.
- For {1,3}: covered singletons = 1 ({2}), outside = 0. Count = 2!/2 = 1.
- For {2,3}: covered singletons = 0, outside = 1 ({1}). Count = 2!/1 = 2.
- Total: 2 + 1 + 2 = 5.

My formula at n=3: `(n−1)!·[n·H_{n−1} − (n−1)] = 2·[3·(1+1/2) − 2] = 2·[4.5 − 2] = 5`. ✓

**Step 4: Coefficient in moment-cumulant formula (Thm 3.26).**

The contribution to E_B(X^n) from all monotone partitions with one 2-block
(of any underlying noncrossing shape) and (n−2) singletons:

`Σ_{π ∈ M(n), one 2-block} (1/|π|!) · κ_2^B · (κ_1^B)^{n−2}`
   = (1/(n−1)!) · [number of such monotone partitions] · κ_2^B · (κ_1^B)^{n−2}
   = (1/(n−1)!) · (n−1)! · [n · H_{n−1} − (n−1)] · κ_2^B · (κ_1^B)^{n−2}
   = [n · H_{n−1} − (n−1)] · κ_2^B · (κ_1^B)^{n−2}

(since `|π|! = (n−1)!` for any partition with `n−1` blocks).

**This contribution grows like `n · log(n) · κ_2^B · (κ_1^B)^{n−2}`** in n —
NOT like `(1/2)^n`.

---

## 5. Reconciliation with Deliverable C §4

MONOTONE_CUMULANTS_C_ASYMPTOTIC.md §4 stated:

> "From Hasebe monograph Defn 3.23: a monotone set partition is non-crossing
> with inner blocks higher in linear order. The 2-block can sit at any
> position, but the noncrossing constraint forces it to be an interval block
> {i, i+1}. The number of choices is (n−1) (positions for the interval)."

**This was incomplete.** The noncrossing constraint does NOT force the
2-block to be an interval (singletons can sit anywhere without violating
noncrossing). The correct count is `n · H_{n−1} − (n−1) ~ n · ln(n)`, not
(n−1).

**Implication for the closed-form derivation:** the leading contribution
to E_B(X^n) from one-2-block monotone partitions is `n · ln(n) · κ_2^B ·
(κ_1^B)^{n−2}`, not `(n − 1) · κ_2^B · (κ_1^B)^{n−2} / (n−1)!`.

This grows even FASTER than the (n−1) estimate in n — still polynomial in n,
still no `(1/2)^n` decay.

**The rate `(1/2)^n` cannot come from the combinatorial structure of monotone
partitions. It MUST come from the B-valued spectral structure of κ_2^B itself
— specifically from the rate-(1/2) eigenvalue of the bilinear pair operator
T_M acting on κ_2^B's domain.**

---

## 6. Disposition of W2.B

**What's rigorous (verbatim Hasebe):**
- Count of monotone partitions on [n] with exactly one 2-block + (n−2)
  singletons is `n · H_{n−1} − (n−1)`.
- Coefficient in Thm 3.26 moment-cumulant formula: `[n·H_{n−1} − (n−1)] /
  (n−1)! · κ_2^B · (κ_1^B)^{n−2}`.
- The combinatorial factor is polynomial-times-log in n, not exponential.

**What this implies:**
The factor `14 = 2·7` does NOT arise from a monotone-partition combinatorial
count. The Hasebe Defn 3.23 enumeration produces growth `n · H_{n−1}`, not a
constant 14. The `(1/2)^n` rate of the subdominant must come from the
**spectral structure of κ_2^B** (the rate-(1/2) eigenvalue of the bilinear
pair operator, R77 §3), NOT from the monotone-partition combinatorics.

**Cleanest statement.** The coefficient `1/30` of the `(1/2)^n` subdominant is
**not derivable from the Hasebe Defn 3.23 / Thm 3.26 monotone-partition
formula alone**. The framework supplies the mechanism (κ_2^B subdominant) and
the per-step structure, but the closed-form amplitude requires the bilinear
pair operator T_M's λ_2 eigenvalue and its eigenvector amplitude — exactly
the open step in R75 §8 / R76 §6 / R77 §6.

---

## Files

- MONOTONE_CUMULANTS_A_VERBATIM.md (verbatim Hasebe Defn 3.23, Thm 3.26)
- MONOTONE_CUMULANTS_C_ASYMPTOTIC.md §4 (prior reasoning, with the count
  error noted)
- W2_KAPPA2_CALC.md (κ_2^B amplitude calculation on (1, 4))
