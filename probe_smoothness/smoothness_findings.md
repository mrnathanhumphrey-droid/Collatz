# Smoothness partition probe — Plancherel mass localization

**Date:** 2026-05-06.
**Brief:** test whether B-smooth vs B-rough residues (classified by max prime
in forward Collatz orbit of smallest positive lift, depth-200) carry
distinguishably different shares of Plancherel mass at the c=7/45 character
group. Tested k ∈ {5, 6, 7}, B ∈ {7, 50, 100, 1000, 10000}.

## Verdict — Outcome A (refined): mass DOES localize, but not on smoothness alone

> **Smooth residues carry systematically less Plancherel mass than their
> count share** across every (k, B) tested where the smooth set is non-trivial
> (mass_ratio_smooth ∈ [0.24, 0.92], not 1.0). Equivalently, rough residues
> carry up to 2× their fair share at moderate B (k=5 B=1000: rough is 31% of
> count but 65% of mass, ratio 2.06).
>
> **But the mechanism isn't "smoothness controls π"** — Pearson correlation
> r(π, max_prime) is +0.30 at k=5 and decays to −0.02 at k=7. The aggregate
> mass-localization signal masks substantial within-partition heterogeneity:
> at k=5 some 7-smooth residues (e.g., r=20, max_prime=5, orbit length 7)
> have high π = 0.019 while some 1000-rough residues have low π. The
> localization is a coarse partition effect, not a fine-grained correlation.
>
> **The cleanest sub-finding** is the "1619 club": at k=5, 8 of the 10
> highest-π residues share max_prime = 1619 with orbit lengths 90–100. At
> k=6, 9 of 10 highest-π residues share max_prime = 1619. This identifies
> a specific orbital family carrying disproportionate Plancherel mass —
> residues whose forward orbits pass through prime 1619 with characteristic
> long-orbit length.
>
> **Cross-level S_k_smooth vs S_k_rough** is methodologically ill-defined
> here: the smooth/rough partition uses each level's own smallest-positive
> lift, so partitions at k=5 and k=6 aren't marginally consistent. The
> additive split S_k = S_k_smooth + S_k_rough sums to S_k_total = 7/15 +
> ε_k always (identity), but the individual values are dominated by partition
> mismatch across levels (S_k_smooth often negative and large).

## Setup

For each r ∈ (Z/3^k)*:
1. Smallest positive lift n = r.
2. Run forward Collatz (n → n/2 if even, 3n+1 if odd) until n=1 or 200 steps,
   capped at max_value 10¹².
3. Track maximum prime factor across all integers visited.
4. Classify B-smooth if max_prime ≤ B and orbit reached 1; else B-rough.

For each (k, B):
```
X_k         = 3^k · Σ_r π_k(r)²
X_k_smooth  = 3^k · Σ_{r smooth} π_k(r)²
X_k_rough   = X_k − X_k_smooth
mass_ratio_smooth = (X_k_smooth / X_k) / (|smooth| / |coprime|)
                  = 1.0 if π² mass uniform across smoothness;
                  ≠ 1.0 indicates localization.
```

All k=5,6,7 residues reached 1 within 200 steps; max orbit length at k=7
was 179. No divergent or capped orbits.

## Headline table — partition geometry

| k | B | n_smooth | n_rough | count_share_sm | mass_share_sm | mass_ratio_sm | mass_ratio_rough |
|---|---|---|---|---|---|---|---|
| 5 | 7 | 14/162 | 148 | 0.086 | 0.062 | **0.72** | 1.03 |
| 5 | 50 | 63 | 99 | 0.389 | 0.233 | **0.60** | 1.26 |
| 5 | 100 | 75 | 87 | 0.463 | 0.240 | **0.52** | 1.42 |
| 5 | 1000 | 111 | 51 | 0.685 | 0.353 | **0.51** | **2.06** |
| 5 | 10000 | 162 | 0 | 1.000 | 1.000 | 1.00 | (n/a) |
| 6 | 7 | 18/486 | 468 | 0.037 | 0.012 | **0.34** | 1.03 |
| 6 | 50 | 95 | 391 | 0.196 | 0.067 | **0.34** | 1.16 |
| 6 | 100 | 122 | 364 | 0.251 | 0.082 | **0.33** | 1.23 |
| 6 | 1000 | 280 | 206 | 0.576 | 0.219 | **0.38** | **1.84** |
| 6 | 10000 | 483 | 3 | 0.994 | 1.000 | 1.01 | 0.02 |
| 7 | 7 | 21/1458 | 1437 | 0.014 | 0.003 | **0.24** | 1.01 |
| 7 | 50 | 130 | 1328 | 0.089 | 0.082 | 0.92 | 1.01 |
| 7 | 100 | 177 | 1281 | 0.121 | 0.091 | 0.75 | 1.03 |
| 7 | 1000 | 579 | 879 | 0.397 | 0.255 | **0.64** | 1.24 |
| 7 | 10000 | 1394 | 64 | 0.956 | 0.984 | 1.03 | 0.36 |

The dominant pattern: **mass_ratio_smooth < 1 robustly**. Rough residues
over-represent mass at every (k, B) tested (except B=10000 trivial case
where smooth ≈ full).

## Refined mechanism — the "1619 club"

Pearson correlation r(π, max_prime):

| k | r(π, max_prime) | r(π, log₁₀ max_prime) | spearman(π, max_prime) |
|---|---|---|---|
| 5 | +0.299 | +0.222 | +0.153 |
| 6 | +0.021 | +0.166 | +0.093 |
| 7 | −0.018 | +0.035 | +0.022 |

The correlation is weak and decays with k. **So the bulk-mass-localization
signal is NOT a smooth monotone "rougher → larger π" relation.** It's a
coarse partition effect: rough residues have a *higher mean π²*, but the
within-partition variance is huge.

What drives the partition-level signal: a small set of residues with
**characteristic long Collatz orbits passing through specific medium-large
primes**. Top-10 highest-π residues at each k:

**k = 5** (top 10 carry 27% of total mass with 6% of count):

| r | π | max_prime | orbit_length |
|---|---|---|---|
| 242 | 0.0452 | **1619** | 96 |
| 152 | 0.0339 | 29 | 23 |
| 107 | 0.0339 | **1619** | 100 |
| 182 | 0.0337 | **1619** | 93 |
| 161 | 0.0336 | **1619** | 98 |
| 121 | 0.0226 | **1619** | 95 |
| 188 | 0.0198 | **1619** | 106 |
| 20 | 0.0194 | **5** | **7** |
| 233 | 0.0194 | **1619** | 83 |
| 206 | 0.0192 | **1619** | 88 |

8 of 10 share max_prime = 1619 with orbit lengths 90–106. The
"1619 club" is a coherent orbital family.

**k = 6** (top 10): 9 of 10 share max_prime = 1619.
**k = 7** (top 10): 5 of 10 share max_prime = 1619; the rest scatter across
{1093, 1367, 41, 577}. Less concentrated at higher k, but 1619 is still the
modal max_prime.

**The exception — r=20 at k=5:** orbit 20→10→5→16→8→4→2→1, length 7,
max_prime = 5 (the smallest possible nontrivial). This residue is 7-smooth
(B=7) but has high π = 0.019. **So 7-smooth ≠ low-π universally.**

This is what produces the within-partition heterogeneity: most 7-smooth
residues have low π, but r=20 is an outlier with high π. The aggregate
mass_ratio_smooth = 0.72 captures the bulk effect; individual residue
behavior is more heterogeneous.

## Cross-level eps reading (caveat: not a Plancherel restriction)

For each B, computed S_k_partition = X_k_partition − X_{k−1}_partition
where each level's smooth/rough partition is built from its own
smallest-positive-lift orbits. **This is not a Plancherel restriction**
because the partitions at k and k-1 are inconsistent (the smooth residues
at level 6 are not the marginal projection of smooth residues at level 7).

| k | B | S_k_smooth | S_k_rough | S_k_total | ε_k |
|---|---|---|---|---|---|
| 6 | 7 | −0.170 | +0.636 | 0.466 | −5.0e-4 |
| 6 | 50 | −0.556 | +1.022 | 0.466 | −5.0e-4 |
| 6 | 100 | −0.518 | +0.984 | 0.466 | −5.0e-4 |
| 6 | 1000 | −0.372 | +0.838 | 0.466 | −5.0e-4 |
| 7 | 7 | −0.034 | +0.500 | 0.465 | −1.2e-3 |
| 7 | 50 | +0.100 | +0.366 | 0.465 | −1.2e-3 |
| 7 | 100 | +0.078 | +0.388 | 0.465 | −1.2e-3 |
| 7 | 1000 | +0.263 | +0.202 | 0.465 | −1.2e-3 |

S_k_total recovers 0.466 (= 7/15 + ε_k) at every B by construction (additive
identity). **The individual S_k_smooth and S_k_rough vary wildly** across B
and even flip signs (S_6_smooth < 0 at every B; S_7_smooth flips positive).
This isn't structural information — it's the partition mismatch dominating.

**Interpretation:** the brief's pre-registered Q1 ("does S_k_smooth converge
to a different value than 7/15") is methodologically ill-posed under the
chosen partition rule. To answer cleanly, the partition would need to be
**marginally consistent** across k: classify at the finest tested k=7 and
project the partition down to k=6, k=5 via mod-3^{k-1} projection. Then
S_k_smooth and S_k_rough trace coherent sequences whose limits could be
compared to 7/15. Deferred.

## Pre-registered outcome reconciliation

| Outcome | Status |
|---|---|
| **A: smooth/rough partition gives distinguishable convergence behavior** | **PRIMARY (refined)** — partition-level mass localization confirmed at every (k, B) tested with non-trivial partition; smooth residues consistently under-weighted in mass; mechanism is partition-coarse, not fine-grained |
| **B: S_k_smooth = S_k_rough = 7/15 trivially** | **REJECTED** at the partition-mass level (smooth carries 0.24–0.92× its share); but cross-level S_k_smooth convergence question is ill-posed under chosen methodology, so B not directly testable |
| **C: conditional measures don't converge** | **NOT TESTED** — would require marginally-consistent partition |
| **D: partition degenerate at all B** | **REJECTED** — partitions are non-trivial at B ∈ {50, 100, 1000} across all tested k |

## What this means for the framework

**Substantive partition-level finding.** The chain's stationary measure π_k —
a 3-adic object computed from Tao-Syracuse iteration — exhibits a non-uniform
distribution across smoothness classes defined by integer Collatz orbits
(a multi-prime, integer-arithmetic object). The chain doesn't "see" the
integer orbit's prime factorization at all; the correlation is emergent
from how (Z/3^k)* residues map to integer-orbit families under the smallest-
positive lift.

**The "1619 club" is the cleanest structural feature.** Residues whose
forward orbits visit prime 1619 with orbit length 90–100 carry
disproportionate Plancherel mass. This is a specific arithmetic family,
not just a partition statistic.

**Why might 1619 specifically?** 1619 is prime. In the Collatz orbit space,
trajectories that hit 1619 share a common "branch" of the inverse Collatz
tree. The residues lifting to such trajectories cluster in (Z/3^k)* in a
specific 3-adic pattern. R58's inverse-tree subtree-size measure
(Pearson 0.86 with D_emp) is the natural framework for understanding why
specific orbital branches carry disproportionate mass; the "1619 club"
finding here is consistent with that picture (R58 says trajectory mass
concentrates on specific inverse-tree subtrees; smoothness classification
just happens to pick out one such subtree at k=5,6).

**The framework's structural constant 7/15 is NOT carried by smooth or
rough residues separately** — it's the aggregate. The partition-level
localization concerns where π²-mass lives (so X_k = 3^k · ‖π‖² is
biased toward rough), but S_k = X_k - X_{k-1} = 7/15 + ε_k holds for
the whole chain, not for smooth/rough separately under this methodology.

## Suggested follow-ups

Ranked by tractability:

1. **Marginally-consistent partition.** Classify at k=7 and project down via
   mod-3^{k-1} reduction. Then S_k_smooth, S_k_rough trace coherent sequences;
   can ask whether they converge to 7/15·f_smooth and 7/15·f_rough or
   different limits. Cheap (~minutes); directly addresses Q1, Q2.
2. **Investigate the "1619 club" structurally.** Why is 1619 the modal
   max_prime among high-π residues? Trace the orbit branches; relate to
   inverse-tree (R58) subtree structure. Could surface a specific
   trajectory family that captures most of the Plancherel mass.
3. **Compare to alternative partitions.** Smoothness is one partition;
   others (orbit length, max value reached, v_2 pattern, R58 subtree size)
   could give cleaner localization signals. Run at k=5,6,7 with each
   partition; compare mass_ratio across partitions. Identifies which
   partition feature most cleanly localizes π² mass.
4. **Per-residue scatter diagnostics.** Plot π(r) vs (max_prime, orbit_length,
   max_value) across all (Z/3^k)* residues. Could reveal hidden 2D structure
   that the 1D max_prime correlation misses.

The cheapest decisive next probe is **(1)** — converts the ill-posed cross-level
question into a well-posed one and answers Q1/Q2 directly.

## Files

- [smoothness_probe.py](smoothness_probe.py) — main probe script
- [smoothness_supplement.py](smoothness_supplement.py) — correlation + extreme-residue analysis
- [smoothness_probe.log](smoothness_probe.log) / [smoothness_supplement.log](smoothness_supplement.log) — full stdout
- [result_smooth_rough_partition.csv](result_smooth_rough_partition.csv) — per-residue (k, r, max_prime, orbit_length, π)
- [result_S_k_conditional.csv](result_S_k_conditional.csv) — per-(k, B) partition geometry
- [result_eps_conditional.csv](result_eps_conditional.csv) — cross-level S_k_smooth/S_k_rough (caveat: partitions inconsistent across k)

## Honest framing

**The probe answers a coarse question (does smoothness localize mass)
positively** but the underlying mechanism is NOT "smoothness drives π" —
it's "specific orbital families carry disproportionate mass, and these
happen to fall outside the strict-smooth set at moderate B." The
"1619 club" finding is the substantive structural pull from this probe;
the smoothness framing is the entry point that surfaced it but isn't
itself the explanatory variable.

The cross-level S_k_smooth question (Q1) is ill-posed under the chosen
partition methodology and would need a marginally-consistent partition
to answer. The cleanest interpretation is partition-mass localization at
fixed k (mass_ratio_smooth < 1 robustly), not "smoothness defines a
sub-chain whose stationary converges to a different limit."

Per the brief's note ("Outcome B (no smoothness signal) would be
informative — would tell us the framework's structural information is
uniform across smoothness-defined subsets, consistent with the structure
being algebraic rather than analytic"): we got *not-quite-B*. Smoothness
alone isn't the right partition (correlation r decays with k), but
*some* orbit-derived partition is non-uniform. The structure has analytic
signatures, but the analytic feature isn't smoothness — it's the specific
orbital family ("1619 club") that smoothness happens to partially detect
at small k.
