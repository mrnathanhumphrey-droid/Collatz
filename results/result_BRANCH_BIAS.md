# RESULT — BRANCH-BIAS: the inverse-tree conditional 0.792 = (√21−3)/2 exactly, growth rate λ=(3+√21)/6 (2026-07-27)

**Probe:** `probes/probe_branch_bias.py`. Closes inverse-tree Phase-5 open question #1 (`inverse_tree_findings.md`):
*"why is P((n−1)/3 odd | n≡1 mod 3) = 0.792?"* — the single conditional that is the entire source of the tree's
branching bias (branching ⟺ n≡4 mod 6, and 0.264 = 0.333 × 0.792). Answer: **it is 1/λ where λ is the tree's growth rate**,
and both are exact algebraic numbers.

## Setup
Inverse Collatz tree (`inverse_tree/build_tree.py`): node `n` has children `left = 2n` (always) and
`right = (n−1)/3` (only when n≡1 mod 3 **and** (n−1)/3 odd ⟺ **n≡4 mod 6**). The "one number" is the node-measure
conditional `P((n−1)/3 odd | n≡1 mod 3) = P(n even | n≡1 mod 3) = P(n≡4 mod6)/P(n≡1 mod3)`.

## Structural key — the source of the asymmetry
**`left = 2n` is ALWAYS even, so odd nodes arise ONLY from the `(n−1)/3` branch of the ≡4-mod-6 nodes.** That is the
entire reason the conditional exceeds ½: even residues are fed by the always-on doubling; odd residues only by the
sparse branch.

## Derivation (Perron eigenvector of the mod-6 offspring matrix)
The tree's residue-count vector obeys `x_{d+1} = A x_d`, `A` = mean-offspring matrix on residues mod 6
(`left`: 0→0,1→2,2→4,3→0,4→2,5→4; `right`: parent 4 → {1,3,5} each ⅓, from the exact three lifts mod 18). The stationary
residue distribution is the Perron eigenvector `v`. By the always-even/only-branch structure, **all even residues share one
value E, all odd residues share one value O**, collapsing `A v = λ v` to two equations:
```
E + O = λ E        (even children: doubling of the previous residue)
(1/3) E = λ O      (odd children: the branch, mass 1/3 from a ≡4 node)
```
⟹ `O = E/(3λ)`, substitute: `1 + 1/(3λ) = λ` ⟹ **`3λ² − 3λ − 1 = 0`** ⟹ **`λ = (3+√21)/6 = 1.2637626`**.

Then `P(n≡1 mod3) = O + E = λE`, so the conditional is
> **`P((n−1)/3 odd | n≡1 mod3) = E/(O+E) = E/(λE) = 1/λ = (√21 − 3)/2 = 0.7912878` — EXACT.**

## Gates (all pass)
- **Exact ⟺ mod 6.** The full mod-`2·3^k` offspring operator (right child's extra 3-adic digit handled by the exact three
  lifts, no equidistribution assumption) gives the *identical* conditional 0.7912878 and `λ=1.263763` at **every k=1…7** —
  mod-6 is exact, no finer correction.
- **Empirical.** `tree_d50.parquet` (379,600 nodes): P(even|≡1mod3) → 0.7905–0.7922 at depths 40–50, cumulative depth≥30 =
  **0.79143**, converging to (√21−3)/2.
- **Closed form vs numerical Perron:** match 1.1e-16.

## The full closed-form picture (three constants, one root)
- **Growth rate / mean offspring / a★_6 λ_max:** `λ = (3+√21)/6 = 1.263763` (root of `3λ²−3λ−1=0`). The old a★_6 chain had
  this only numerically (1.2638, eigvec-match 3e-4); it is now exact.
- **Branching ratio** `P(n≡4 mod6) = 1/(3λ) = (√21−3)/6 = 0.263763`.
- **`P(n≡1 mod 3) = 1/3` exactly** (the mod-3 marginal is uniform even though the full mod-6 is not).
- **The conditional** `= 1/λ = (√21−3)/2 = 0.791288`. Elegant reading: *the fraction of ≡1-mod-3 nodes that branch is exactly
  the reciprocal of the tree's growth rate.*
- **Consistency:** `P(≡1mod3)·conditional = (1/3)·(√21−3)/2 = (√21−3)/6 = branching ratio` ✓.

## Bonus correction (inverse-tree growth exponent)
The tree grows like `λ^d`, exponent **`log λ = log((3+√21)/6) = 0.234093`** — **not** `log(4/3)=0.287682`. The empirical
layer-count slope 0.2343 (`phase1_checkpoint.md`, window d=10..50, R²=0.9998) had already reached the *true* asymptote; the
checkpoint's "finite-window slope expected lower than log(4/3)" note was a misread — 0.234093 **is** the asymptote.

## Net
The 0.792 conditional, the 0.264 branching ratio, and the 1.2638 growth rate are **one algebraic object**, all determined by
`3λ²−3λ−1=0`, `λ=(3+√21)/6`. Closes inverse-tree Phase-5 OQ#1 exactly and pins the a★_6 `λ_max` in closed form. This is the
2-adic prefix arc's structure, disjoint from the 3-adic 7/15 cascade (per `PROP_TEST_RESULTS.md`); it does not bear on S_∞.
**Not at stake:** anything in the current 7/15 / SOLSTICE / GARSIA program.
