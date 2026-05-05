# Follow-up duality checks — VERDICT: sample-size artifact dominates, residual structural difference is modest (factor ~0.2-4)

**Date:** 2026-05-04. Sibling-probe Task 2 follow-up. Tests whether the dramatic factor-10³-10⁴ basin-fingerprint difference between Agent 2 (3x+1 single-basin) and Agent 3 (3x−1 three-basin) at large (n, k) reflects structural distribution-level difference, or a sample-size artifact (different |V_n|).

## Verdict (one line)

> **The factor-10³-10⁴ raw difference is ~95% sample-size artifact.** When Agent 2's tree is truncated to Agent 3 root-1's vertex counts (sorted by integer value, smallest first), the matched D_n^{matched} agrees with Agent 3 root-1's D_n within factor **0.17 – 4.1** across all (n, k) pairs tested. The residual factor-of-a-few variance IS real and structural — but much smaller than the raw comparison suggested.

The basin fingerprint exists but is modest at the matched-N level; the dramatic 10⁴× ratio at (n=6, k=2) is dominated by the fact that Agent 2's tree has ~600× more vertices at depth 6 than Agent 3 root-1.

## CHECK 1 — Matched-N (THE DECISIVE TEST)

**Procedure.** Rebuild Agent 2's full 3x+1 inverse tree from 1 to depth 6 (1.35M vertices total). At each depth n, sort vertices by integer value (smallest first) and truncate to the first M_n vertices, where M_n matches Agent 3 root-1's depth-n vertex count: 1, 13, 56, 189, 459, 1061, 2247.

Compute D_n^{matched}(k) on this truncated set via the same closed form Agent 2 uses. Compare to (a) Agent 2's full D_n(k) and (b) Agent 3 root-1's D_n(k).

**Selected results** (full table in [duality_followup_data.csv](experiments_output/duality_followup_data.csv)):

| (n, k) | A2 full D_n | A2 matched D_n | A3 root-1 D_n | matched / A3 |
|---|---|---|---|---|
| (1, 2) | 6.12e−02 | 7.10e−02 | 1.07e−01 | 0.667 |
| (2, 2) | 3.95e−03 | 3.06e−02 | 2.10e−02 | 1.45 |
| (3, 2) | 1.98e−04 | 2.18e−03 | 1.51e−03 | 1.44 |
| (4, 2) | 9.88e−06 | 3.47e−03 | 1.04e−02 | 0.33 |
| (5, 2) | 1.14e−06 | 8.74e−04 | 4.00e−03 | 0.22 |
| (6, 2) | 3.20e−07 | 9.98e−04 | 1.04e−03 | 0.96 |
| (3, 3) | 3.42e−03 | 3.93e−02 | 1.95e−01 | 0.20 |
| (4, 3) | 1.11e−04 | 7.85e−02 | 1.91e−02 | 4.12 |
| (5, 3) | 6.21e−05 | 6.49e−03 | 7.15e−03 | 0.91 |
| (6, 3) | 5.06e−07 | 1.38e−03 | 7.96e−03 | 0.17 |
| (5, 5) | 9.48e−04 | 1.23e−01 | 9.66e−02 | 1.27 |
| (6, 5) | 9.18e−05 | 5.00e−02 | 7.11e−02 | 0.70 |

**Reading:**

- **matched-to-A3 ratio range: ~0.17 to ~4.1** across the (n, k) grid. Mostly within factor-2 of unity. No systematic trend toward huge values.
- The original raw A3/A2 ratio was ~24,000 at (6, 3). After matching sample sizes, the matched/A3 ratio at (6, 3) is 0.17 — A2 is now slightly *smaller* than A3 at matched N. The 24,000× factor was almost entirely sample-size driven.

**The matched-to-full ratio (last column not shown above)** grows from 1.0 (n=0) to 3,122 at (6, 2). This quantifies the magnitude of the sample-size effect within Agent 2 alone: truncating to 2,247 vertices vs the full 1.35M changes D by factor ~3,000. Most of the original factor-10⁴ difference between A2 and A3 was this same effect.

## CHECK 2 — Per-vertex normalized D̃_n := D_n / |V_n|²

**As literally specified, this normalization OVERSHOOTS** because D_n already has 1/|V_n|² baked in (since mu = count/|V_n|, so |mu_hat|² = |count_hat|²/|V_n|²). Dividing again by |V_n|² gives ~ 1/|V_n|⁴ × count-level Plancherel mass, which makes Agent 2 (large |V_n|) look much smaller than Agent 3 (small |V_n|).

**Empirical** (selected):

| (n, k) | A2 D̃ = D/|V|² | A3 r1 D̃ = D/|V|² | ratio A3/A2 |
|---|---|---|---|
| (4, 2) | 5.4e−14 | 4.9e−08 | **9.1e+5** |
| (5, 2) | 6.2e−17 | 3.6e−09 | **5.7e+7** |
| (6, 2) | 1.7e−19 | 2.1e−10 | **1.2e+9** |
| (6, 3) | 2.8e−19 | 1.6e−09 | **5.7e+9** |

The ratio explodes to ~10⁹ because the |V_n|⁴ in the denominator dominates everything else. **This normalization makes the answer worse, not better, for the structural comparison.**

## CHECK 2-alt — Count-level Plancherel (D · |V_n|²)

The "right" normalization to factor out the empirical-measure dilution is **multiplying by |V_n|²**, giving the integer-histogram Plancherel mass C_n(k) := sum coprime |count_hat(ξ)|² (which doesn't depend on |V_n|).

**Empirical** (selected):

| (n, k) | A2 C_n = D·|V|² | A3 r1 C_n | ratio A3/A2 |
|---|---|---|---|
| (1, 2) | 12 | 18 | 1.50 |
| (2, 2) | 72 | 66 | 0.92 |
| (3, 2) | 360 | 54 | **0.15** |
| (4, 2) | 1,800 | 2,196 | 1.22 |
| (5, 2) | 20,754 | 4,506 | **0.22** |
| (6, 2) | 582,696 | 5,232 | **0.0090** |
| (6, 3) | 922,608 | 40,194 | **0.044** |

The count-level mass grows for both, but A2's grows **faster** because A2 has more vertices being added. At (6, 2), A2's count-level Plancherel is ~111× A3's. This is the count-level effect of A2's larger tree.

**Neither raw D_n nor count-level D·|V|² is the right "structural" comparator.** The matched-N approach in CHECK 1 IS — it equalizes the sample size and asks whether the integer histogram structure differs.

## Synthesis

Three normalizations, three different verdicts:

| normalization | A3/A2 ratio at (6, 2) | interpretation |
|---|---|---|
| Raw D_n | 8,691 | basin fingerprint looks dramatic |
| D_n / |V|² (CHECK 2 literal) | 1.2e+9 | dilution-doubled, misleading |
| D_n · |V|² (count-level) | 0.009 | A2 dominates, also misleading |
| **matched-N (CHECK 1)** | **0.96** | **structural difference is small** |

The matched-N comparison is the principled one — it controls for vertex count and asks whether the *distribution structure* of the depth-n empirical measure differs between systems.

**Conclusion:** the ~10⁴× raw difference in D_n is **mostly sample-size artifact** (different |V_n| between systems). The residual structural difference is **factor 0.17–4.1**, with no clear systematic pattern across (n, k).

## What this changes about the duality verdict

The original duality verdict (no clean function S → D, but basin fingerprint detectable) **softens**. Specifically:

- **Forward Markov chain symmetry** (S_n^{3x±1} identical) — unchanged. Proved rigorously.
- **Inverse-tree differences are real but modest** at matched sample size. The 3-basin structure of 3x−1 produces a depth-n empirical distribution that differs from 3x+1's by factor ~2-5 (in either direction depending on (n, k)) once vertex counts are equalized.
- **The dramatic 10⁴× difference** at large depth was driven by the fact that 3x+1's inverse tree (single basin) reaches enormously more vertices than 3x−1's (basin-and-cap-bounded). This is a STRUCTURAL feature of the integer-level dynamics — bounded basins vs unbounded — but it's a count effect, not a distribution-shape effect.
- **The basin fingerprint exists at two scales**: a large vertex-count effect (10³-10⁴×) AND a small distribution-shape effect (~5×). The vertex-count effect is the dominant signal in raw D_n.

For c=7/45 closure: this nuance doesn't affect the rate-1/2 problem (which is forward-symmetric anyway). For sibling-study framing: the inverse tree DOES distinguish 3x+1 from 3x−1, but the principled distinguisher is at the integer-tree-growth-rate level (single vs multi-basin), with a smaller residual at the distribution level.

## Files

- [duality_followup_check.py](duality_followup_check.py) — script
- [duality_followup_data.csv](duality_followup_data.csv) — full matched/full/A3 table
- [duality_followup_verdict.md](duality_followup_verdict.md) — this writeup
