# PADE_EXTENSION_TABLE — Phase 1: extended Padé table

**Date:** 2026-05-12. Wilson. Phase 1 of the Padé-extension probe.

---

## Headline finding

**With ε_n cached through k=6 (5 input coefficients of the auxiliary series f̃(z) = (E(z) − ε_1 z)/z²), no NEW Padé approximants beyond R77.6's set are computable.** The constraint m + n ≤ #(coefficients) − 1 = 4 is already saturated by R77.6's enumeration:

| m+n | Approximants R77.6 enumerated |
|---|---|
| 2 | (1,1) |
| 3 | (2,1), (1,2) |
| 4 | (3,1), (2,2), (1,3), (4,0), (0,4) |

The (m+n=4) row is complete — every (m,n) with m+n=4 and m,n ≥ 0 is in R77.6's set. Padé approximants with m+n=5 would need ε_7 (k=7 Markov chain compute, ~hours).

So Phase 1's "extended table" task collapses to **re-examination of the existing eight approximants for full pole structure (secondary poles, complex poles, cluster patterns) beyond just closest-to-z=2**.

---

## Full pole structure per approximant

Floating-point pole locations (numpy.roots on Q(z) with float-converted exact-rational coefficients; precision ≤ 10⁻¹⁰ for polynomial degree ≤ 3, ≤ 10⁻⁸ for degree 4 — well below the |z−2|-distance scale we report).

### [1/1] (lowest-order diagonal)

| pole idx | z (complex) | \|z\| | \|z−2\| | role |
|---|---|---|---|---|
| 0 | +2.0764 +0.0000j | 2.076 | 0.0764 | PRIMARY |

Single pole only (Q is linear). Pure real, positive, just above z=2.

### [2/1] (sub-diagonal m+n=3)

| pole idx | z (complex) | \|z\| | \|z−2\| | role |
|---|---|---|---|---|
| 0 | +2.1292 +0.0000j | 2.129 | 0.1292 | PRIMARY |

Single pole (Q linear). Pure real, positive, above z=2. **Closer to z=2 than [3/1] but FURTHER than [1/1]** — sub-diagonal pushes pole away from z=2 relative to diagonal.

### [1/2] (super-diagonal m+n=3)

| pole idx | z (complex) | \|z\| | \|z−2\| | role |
|---|---|---|---|---|
| 0 | +2.1299 +0.0000j | 2.130 | 0.1299 | PRIMARY |
| 1 | +155.4074 +0.0000j | 155.41 | 153.4074 | spurious |

Two poles (Q quadratic). Primary near z=2.13, spurious at z≈155. The spurious pole is a classical Padé artifact: when the linear system for q_1, q_2 forces a near-cancellation, one root flies to infinity. Doesn't affect the primary reading.

### [3/1] (sub-diagonal m+n=4)

| pole idx | z (complex) | \|z\| | \|z−2\| | role |
|---|---|---|---|---|
| 0 | +2.3132 +0.0000j | 2.313 | 0.3132 | PRIMARY |

Single pole. Pure real, FURTHEST from z=2 in the primary-cluster set (0.313 vs 0.05-0.13 for others). High m+low n approximants push the closest pole AWAY from z=2 — explained by branch-cut structure (taking many Taylor terms but few poles forces the pole to "stand in for" a much extended cut).

### [2/2] (diagonal m+n=4)

| pole idx | z (complex) | \|z\| | \|z−2\| | role |
|---|---|---|---|---|
| 0 | +2.0513 +0.0000j | 2.051 | 0.0513 | PRIMARY |
| 1 | +0.6878 +0.0000j | 0.688 | 1.3122 | secondary |

Two poles. PRIMARY at z=2.0513 (CLOSEST to z=2 across all approximants). Secondary at z=+0.6878 — pure real, INSIDE the unit disk-ish. **Reading of the secondary:**

- NOT consistent with R76 §10's predicted next singularity at z=4 (would expect secondary near z=4).
- NOT a complex secondary (pure real).
- Position z=0.69 has no immediate structural meaning for the Markov chain (no obvious connection to 1/2, 1/3, or 7/15 reference points).
- Most likely a Padé artifact at low order: when the diagonal [2/2] tries to capture both the leading rate-1/2 (z=2) and the next-order structure with only two poles, the secondary "stands in" for any structure not accounted for by the primary.

Treat as **likely artifact pending higher-order diagonal extension** (would need ε_7, ε_8 for [3/3], [4/4]).

### [1/3] (super-diagonal m+n=4)

| pole idx | z (complex) | \|z\| | \|z−2\| | role |
|---|---|---|---|---|
| 0 | +2.3485 +0.0000j | 2.348 | 0.3485 | PRIMARY |
| 1 | +7.9778 +0.0000j | 7.978 | 5.9778 | spurious |
| 2 | −11.9889 +0.0000j | 11.989 | 13.9889 | spurious |

Three poles. Primary near z=2.35, two large-magnitude spurious poles. Like [1/2], the super-diagonal with high-degree denominator generates spurious poles at large |z|. All real, no complex.

### [4/0] (Taylor truncation)

No poles (degree-0 Q). Sanity-check entry — confirms the Taylor truncation framework. No new information.

### [0/4] (all-pole approximant)

| pole idx | z (complex) | \|z\| | \|z−2\| | role |
|---|---|---|---|---|
| 0 | −0.9514 +0.7648j | 1.221 | 3.0488 | cc-pair |
| 1 | −0.9514 −0.7648j | 1.221 | 3.0488 | cc-pair |
| 2 | +0.4526 +1.0072j | 1.104 | 1.8464 | cc-pair |
| 3 | +0.4526 −1.0072j | 1.104 | 1.8464 | cc-pair |

Two complex-conjugate pairs. Both on circles inside |z|=2:
- |z|≈1.22 pair at arg ≈ 141°
- |z|≈1.10 pair at arg ≈ 66°

**Reading of [0/4]:** [0/n] is the "all-pole" approximant — numerator forced to constant. For functions with branch structure at z=ρ, the [0/n] approximant places poles around circles of radius slightly less than ρ (Padé "Froissart"-style behavior). For ρ=2, finding cc-pairs at |z|=1.10 and |z|=1.22 is consistent with this artifact pattern. **These are NOT real secondary singularities of E(z).** The reading from [0/4] does not contradict the readings from the diagonal/near-diagonal approximants.

---

## Summary table (closest pole only, sorted by |z−2|)

| (m,n) | role | closest pole | \|z−2\| |
|---|---|---|---|
| (2,2) | diagonal | +2.0513 +0.0000j | **0.0513** |
| (1,1) | diagonal | +2.0764 +0.0000j | **0.0764** |
| (2,1) | sub-diag | +2.1292 +0.0000j | 0.1292 |
| (1,2) | super-diag | +2.1299 +0.0000j | 0.1299 |
| (3,1) | sub-diag | +2.3132 +0.0000j | 0.3132 |
| (1,3) | super-diag | +2.3485 +0.0000j | 0.3485 |
| (0,4) | all-pole | +0.4526 +1.0072j | 1.8464 (cc-pair off-axis) |

Diagonal sequence is tight at the closest-to-z=2 end (0.05, 0.08). Off-diagonals fan outward in a clear pattern: more numerator degree pushes the pole further from z=2 with these N=5 coefficients.

---

## What's observed about secondary structure

1. **No primary-pole complex-conjugate pairs** in any of the (1,1), (2,1), (1,2), (3,1), (2,2), (1,3) approximants. The [0/4] cc-pairs are isolated to that specific approximant and behave like classical branch-cut artifacts of the all-pole construction.

2. **One real secondary pole** appears in (2,2) at z=0.6878 — no consistent partner in any other approximant. Likely artifact.

3. **Two large-magnitude spurious poles** appear in (1,3) at z≈+8 and z≈−12, and one in (1,2) at z≈155 — all classical Padé super-diagonal artifacts.

4. **No clustering of multiple poles near z=2** suggesting branch-cut densification — would need higher-order approximants to see (diagonals beyond [2/2]).

---

## Verdict (Phase 1 only)

The extended pole table CONFIRMS R77.6's branch-cut reading and **rules out a complex-conjugate primary structure** (which would have been H_COMPLEX_SECONDARY). It does not advance the branch-order disambiguation — no new approximants are computable. The phase shifts to Phase 2 (pattern analysis across the existing approximants) and Phase 3 (ratio diagnostic) for further information.
