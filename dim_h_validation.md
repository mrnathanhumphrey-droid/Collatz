# Bowen-Sullivan-Pollicott dimension validation — outcome (γ)

**Status.** Decisive. The claimed match in Result 57 between
`2·log(λ_max(M_closed))/log(2) = 0.6755` and Chang's "Hausdorff dimension ≈ 0.68"
is numerical happenstance against a rounded heuristic, not structural
identification. Chang's rigorous figure is **dim_H(C) = log(φ)/log(2) ≈ 0.69424**
(exact closed form, golden-ratio dimension), which differs from our 0.6755 by
**0.019**, not 0.005 as Result 57 stated.

The factor of 2 in our formula has no rigorous derivation from any natural
Bowen pressure equation on M_closed. The natural Furstenberg-Hutchinson
formula `log(λ_max)/log(2) = 0.3378` describes M_closed's branching entropy
in 2-adic norm — a structurally meaningful but DIFFERENT quantity from
Chang's survivor-set dim_H.

**Action for v3.6 / framework synthesis:** walk back the central "same
fractal" claim. The two operators describe different dynamics on different
sets; their leading eigenvalues encode different invariants.

## Step 1 — M_closed setting

Result 23's `M_closed` is the natural-density transition matrix of the
**backward inverse-Collatz tree** on residues mod 2^k:

  M[r, 2r mod 2^k]  += 1                  (doubling, ALWAYS)
  M[r, child_r]     += 1/3                (inverse-3, only EVEN r)

This is a **multitype Galton-Watson branching matrix**, not a probability
transition matrix. Row sums:
- ODD rows: 1 (only doubling)
- EVEN rows: 4/3 (doubling + 1/3 inverse-3)

Leading eigenvalue λ_max sits in (1, 4/3).

**Empirically (k ∈ {5, 6, 8, 10}):** λ_max = 1.263763 invariant to 6 decimals
across all k. Sub-dominant 0.263763, gap |λ_2/λ_1| = 0.2087. Real positive,
isolated, irreducibly Perron-Frobenius.

## Step 2 — Chang's actual dim_H computation

Chang 2603.11066v6, Section 11 (line 7724-7742). The 2-adic Hausdorff
dimension of the survivor set C ⊂ Z_2 is computed via the **Bowen pressure
equation** P(s) = log_2(λ(T(s))) = 0, where:

```
T(s) = [[ 2^(-2s),  2^(-s),   0      ],
        [ 2^(-2s),  0,        2^(-s) ],
        [ 2^(-2s),  0,        2^(-s) ]]
```

indexed by **safe forward Syracuse classes {1, 3, 7} mod 8**, with weights
2^(-v·s) for valuation cost v ∈ {1, 2}.

**Closed-form solution.** All rows of T(s) sum to 2^(-2s) + 2^(-s). Row
stochasticity gives the eigenvalue λ = row sum. Setting u = 2^(-s):

  λ(T(s)) = 1 ⟺ u² + u = 1 ⟺ u = (-1 + √5)/2 = 1/φ

where **φ = (1+√5)/2 is the golden ratio**. Therefore:

  **dim_H(C) = s = log(φ)/log(2) = 0.6942419136...** (closed form, exact)

Chang's paper reports 0.6942 to 4 decimals; my replication gives 0.694241913
to 9 decimals (`exp 80`). The rounded "0.68" cited in some sections is a
heuristic, not the rigorous value.

## Step 3 — Our claimed match: numerical comparison

| quantity | value | source |
|---|---:|---|
| Chang dim_H (rigorous) | 0.6942419136 | T(s) Bowen pressure |
| Chang heuristic (rounded) | 0.68 | informal references |
| Our 2·log(λ_max)/log(2) | 0.6754519 | Result 57 |
| log(φ)/log(2) | 0.6942419136 | exact closed form |

| comparison | gap |
|---|---:|
| Our 0.6755 vs rigorous 0.6942 | **0.0188** (NOT within 0.005) |
| Our 0.6755 vs heuristic 0.68 | 0.0045 (within 0.005) |

**Result 57 compared to the heuristic, not the rigorous value.** This was
the source of the "match within 0.005" claim. Against the rigorous value,
the gap is 4× larger.

## Step 4 — Bowen pressure applied to M_closed: tested, no match

Tested three natural parametrizations of M_closed(s) where transitions are
weighted by 2^(-cost·s), solving λ(M(s)) = 1 for s:

| (cost_double, cost_inv3) | rationale | s |
|---|---|---:|
| (1, 0) | natural 2-adic: doubling = +1 bit, inv-3 = unit | 0.4150 |
| (1, log_2(3)) | doubling = 2-adic, inv-3 = 3-adic | 0.3077 |
| (1, 1) | both contribute equally | 0.3377 |
| target | s = 0.6942 needs scale = exp(log_λ/0.6942) | 1.4010 |

**None of these natural parametrizations land at Chang's 0.6942.** Reverse-
engineering: for log(λ_max)/log(scale) = 0.6942, the scale needs to be
1.4010 — not a clean physical quantity (between √2 = 1.4142 and 4/3 = 1.3333).

## Step 5 — Candidate dimension formulas

| formula | value | gap vs rigorous | gap vs heuristic |
|---|---:|---:|---:|
| log(λ)/log(2)  [Furstenberg-Hutchinson] | 0.3377 | −0.357 | −0.342 |
| 2·log(λ)/log(2)  [our claim] | 0.6755 | **−0.019** | −0.005 |
| log(λ)/log(4/3)  [K_h scale] | 0.8137 | +0.120 | +0.134 |
| log(λ)/log(3) | 0.2131 | −0.481 | −0.467 |
| log(λ)/log(3/2) | 0.5773 | −0.117 | −0.103 |
| 1 − log(λ)/log(2)  [codim] | 0.6623 | −0.032 | −0.018 |

The closest natural formula is `1 − log(λ)/log(2) = 0.6623`, the
**density codimension** of inverse-tree-reachable integers in [1, N] —
since #(reachable at depth d ≤ N) ~ λ^d / 2^d ~ N^(log(λ)/log(2) − 1) gives
density ~ N^(−0.6623). But this is NOT what Chang computes either.

## Step 6 — Sensitivity: where would M_closed need to be?

For 2·log(λ_max)/log(2) = 0.6942:
  λ_max needs to be √φ ≈ 1.272020
  Our value: 1.263763
  Gap: −0.0083

For log(λ_max)/log(2) = 0.6942:
  λ_max needs to be φ ≈ 1.618034
  Tested perturbation: weight inverse-3 step at **1** (not 1/3): λ_max → exactly φ.

The φ-eigenvalue appears when inverse-3 acts on every EVEN residue with
weight 1 — but that's NOT M_closed's structure (inverse-3 has natural-density
weight 1/3 from the 1-in-3 lift selection). The natural M_closed sits at
λ_max = 1.2638, structurally independent of φ.

## Step 7 — Why the operators describe different things

**M_closed (backward, Galton-Watson):**
- Acts on residues mod 2^k in inverse-Collatz tree
- λ_max = mean offspring number per generation
- Limit set = closure of "all integers reaching 1" → all of Z_2 conjecturally (full dim 1)
- log(λ)/log(2) = 0.3378 = 2-adic Furstenberg dim of branching IFS
  (heuristic; not the actual H-dim of any sharp set)

**Chang's T(s) (forward, parametrized transfer):**
- Acts on safe Syracuse classes {1, 3, 7} mod 8 with valuation weights 2^(-v·s)
- Bowen pressure eq solves for s = dim_H(C) of survivor set
- C = divergent-orbit set in Z_2, Cantor-type, dim_H = log(φ)/log(2)
- Real Hausdorff dimension of a real fractal subset

These operators are NOT in the same family. They describe complementary
phenomena (backward growth vs forward survival) on the same integer support
but their leading eigenvalues encode unrelated invariants.

## Verdict

**Outcome (γ).** Numerical coincidence against a rounded heuristic.

The 0.6755 ≈ 0.68 match in Result 57 came from comparing our value to the
rounded "0.68" heuristic that appears in Chang's introduction. Against the
rigorous 0.6942 (= log(φ)/log(2), exact closed form), the gap is 0.019 and
the "factor of 2" in our formula has no derivation.

The framework synthesis chapter's claim that "M_closed and Chang's P_K live
on the same fractal invariant set Λ" needs revision before external
correspondence. The two frameworks compute different quantities; their
agreement at 4-decimal precision does not extend to a structural
identification.

## What IS structurally meaningful from M_closed

- **λ_max = 1.263763 invariant across k ∈ {5..11}**: Result 23's principal
  finding stands. The inverse tree's natural-density branching matrix has a
  k-stable Perron eigenvalue.
- **log(λ_max)/log(2) = 0.3378**: Furstenberg-Hutchinson dimension of the
  inverse-tree branching IFS in 2-adic norm. This is the entropy rate per
  backward step; each backward step adds 0.338 bits of which-integer info.
- **Density codim = 1 − 0.3378 = 0.6622**: thinning rate of inverse-tree-
  reachable integers as N grows. Real and structural.

**What is NOT meaningful from M_closed:** any direct identification with
Chang's survivor-set dim_H. The two are different invariants of different
operators describing different dynamics.

## Salvageable structural connection (if any)

Both Chang's dim and M_closed's λ_max involve specific algebraic numbers
(golden ratio for Chang; an algebraic number near 1.26376 for M_closed).
Whether these arise from a deeper common structure — e.g., both being
eigenvalue-equations on Markov chains over {1, 3, 7}-class Syracuse safe
states — is open. Currently no derivation links them.

The relationship dim_forward + dim_backward = 1 (Hausdorff-codimension
duality) gives a near-miss: 0.6942 + 0.3378 = 1.0320, off by 0.032.
Suggestive but not exact.

## For framework synthesis chapter

Replace "the dimension match strongly suggests both frameworks live on the
same fractal" with:

> Result 23's M_closed has Perron eigenvalue λ_max = 1.263763 invariant
> across moduli mod 2^k for k ∈ {5..11}. The implied 2-adic Furstenberg
> branching dimension log(λ_max)/log(2) = 0.3378 is structurally distinct
> from Chang's rigorous survivor-set dim_H(C) = log(φ)/log(2) = 0.6942.
> The two encode complementary forward/backward Collatz dynamics; no
> direct identification holds.

This is honest and preserves Result 23 / Chang's individual findings without
making claims the validation does not support.

## Files

- `experiments/80_dim_h_validation.py` — main validation, all candidate formulas
- `experiments/81_dim_h_sensitivity.py` — perturbations, golden-ratio surfacing
- `experiments_output/dim_h_calculations.csv` — formula table

Compute: <2 seconds total (linear algebra on small matrices).
