# RESULT — PROBE MAXMODE2: the decisive saturation check + U² identity (2026-07-26)

**Probe:** `probe_maxmode2.py`. Wilson's hard flag: "max mode slowing, 4% margin" = the 0.89 configuration (rate still
drifting, read as an asymptote). Deparity + two-step the max-mode rate to k=16 before any lemma; if climbing, the lemma
is false as stated, not merely unproved.

## (1) DECISIVE — the max-mode rate SATURATES (the good outcome)

Deparitied / two-step amplitude rate of `P_max(k)=max_{3∤a}|μ̂_k(a)|²`:
- **early (k=8–10) mean 0.6544 ; late (k=14–16) mean 0.6607 ; drift +0.0063 over 8 levels.**
- **linear slope +0.0013/level — flat.**

The rate genuinely **saturates at ~0.655 amplitude**, it does not climb. R66's "slowing" (k=3→6: 0.668→0.743) was the
*transient*; by k=8–16 it has plateaued. **This is the opposite of the 0.89 trap** (which was still drifting at the
edge). And the 0.76 spikes that gave the "4% margin" were the raw period-2 artifact — deparitied, the true saturated
rate is 0.655, **comfortably under the 0.79 ceiling**. So a uniform bound can exist; the lemma is not false-as-stated.

## (2) The max MIGRATES — it's a sup over a moving frequency, not a fixed eigenvector

`a_max/3^k`: k=3: 0.259, k=5: 0.342, k=7: 0.127, k=9: 0.246, k=11: 0.124, k=13: 0.094, k=15: 0.106. The maximizing
frequency is **not fixed** — it drifts, tending toward `a/3^k ≈ 0.09–0.13` (or its conjugate 0.87–0.91) at large k.
Consequence for the lemma: the target is `sup_a |μ̂_k(a)|²` over a **migrating** location, not a single tracked
eigenvalue — harder than an eigenvector bound, and it means the lemma is genuinely a sup-norm control on R66's chain.

## (3) The U² identity — verified exactly, and U² is INVARIANT

`mean_k|γ_k(m)−1|² = Σ_{a≠0}|μ̂_k(a)|⁴` (= Gowers U² norm of ρ_k, minus the trivial mode), machine-exact (rel 1e-16
at k=12,14,16). **U²(ρ_k) = 0.29754, constant across levels** (per-level ratio 1.0000) — an invariant, the ℓ⁴ analog
of R66's ℓ² primitive-sum invariant (7/15). Decomposition: the between-class variance from the pinned means is
**exactly 2/9** (`(1/3)(5/3−1)²+(2/3)(2/3−1)² = 2/9 = 0.2222`); the remaining 0.0753 is within-class. So the channels'
mean-square spread around the pinned mean is level-invariant, with the between-class part exact.

This confirms Wilson's shelf move: the aggregate `Σ|μ̂|² = X−1 → ∞` is the wrong norm (triangle inequality hopeless,
phase cancellation does all the work); the native object is **ℓ⁴ = U²**, additive combinatorics, and it's invariant.
Honest limit (Wilson's): U² controls the family in **mean square**, not k=3,4 individually — so the max-mode target
survives for those two specific channels.

## (4) R66 answer (b): the ~2^{−k} max decay is FITTED, not derived
R66 §9 (on disk): "Max over primitive a: empirically follows (1/2)^{k−c} … NOT yet clean closed form." So **only the
AVERAGE 3^{−(k−1)} is derived** (Parseval + the invariant S_∞); the max law is fitted. ⟹ the lemma cannot cite R66's
max law — it must **derive** a sup bound on `|μ̂_k(a)|²` from R66's explicit stationary π_r.

## Net for the pen — green light, with the object named honestly
- **The rate saturates (~0.655, flat), not climbs** — there is a lemma to write, with room (0.655 vs 0.79 ceiling; the
  4% margin was a deparity artifact).
- **The object is `sup_a |μ̂_k(a)|²` on R66's explicit chain** — a migrating sup, derived (R66's max law is only
  fitted), NOT R74's aggregate operator (average, too fast) and NOT a fixed eigenvector.
- **The U² identity is the paper's mean-square statement** (invariant 0.2975, between-class part exactly 2/9); it
  controls the family but not k=3,4 individually, so the per-channel max-mode bound is still the binding lemma.

**Not at stake:** MEAN1, HIERARCHY, CHANNEL_ID, R1–R30. Cheap (cached ρ + build_nu(11), 8.5s).
