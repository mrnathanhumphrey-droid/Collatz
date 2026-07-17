# Result 59: Z_2 measure framework — variant b same as a (Pearson 0.86), measure is multifractal not Sullivan-conformal, σ-band match has clean mechanistic explanation

**Date:** 2026-05-03. Tests the Z_2 measure framework beyond Result 58. Goal: variant (b) branching-density weighting; Hausdorff dim of measure; σ-band conditional analog; conformality check.

**Verdict: outcome (β/γ) refined.**
- **(β) Pearson +0.86 holds** under variant (b) branching-density too — no improvement over Result 58's variant (a).
- **(γ) Sullivan-conformality FAILS** — image/preimage mass ratios span 0.92 (r=5) to 2.34 (r=21) under one Syracuse step. No constant Hausdorff dim δ fits.
- **(α) σ-band mechanism confirmed:** q1 (low-depth, near-root) reproduces D_emp at t=90 with Pearson +0.86; q2-q4 anticorrelate. Predicted by the structural mechanism: D_emp at survivor-time t = residue distribution over orbits with σ slightly > t whose position at iteration t is at low inverse-tree depth.
- **Multifractal:** mass dimension dim_q2 declines from 0.83 (k=7, coarse) to 0.54 (k=15, fine). Crosses Chang's H-dim 0.68 around k=12. NOT a uniformly-fractal measure.

Code: `zadic_measure_framework.py`. Compute: ~3s at N=2^22.

---

## 1. Setup

Inverse Collatz tree from m=1 at value-truncation N=2^22 (~1.25M odd nodes), per Result 58. Tested:

| Variant | Weight w(m) |
|---|---|
| (b) accumulator | Σ over descendants of n_branches(d) |
| (b) direct | n_branches(m) only, no descendants |
| (b) combined | subtree-size(m) × (1 + n_branches(m)) |

n_branches(m) = number of inverse predecessors of m within tree.

## 2. Variant (b) branching-density: marginal improvement, same magnitude as (a)

| Variant | Pearson(•, D_emp) | MAE |
|---|---|---|
| (a) subtree-size [Result 58] | **+0.857** | 0.118 |
| (b) accumulator | **+0.864** | 0.117 |
| (b) direct | +0.171 | 0.319 |
| (b) combined | +0.866 | 0.120 |

Branching density adds essentially nothing. The subtree-size already captures the relevant structural signal. Direct (just n_branches per node, no descendant accumulation) gives weak +0.17, comparable to the uniform-on-tree baseline (Result 58's variant f).

**Conclusion:** the structural object is the subtree-size weighting; branching multiplicity per node provides no additional signal beyond what subtree-size already encodes.

## 3. Hausdorff dim: multifractal, mass-dim declines with scale

Box-counting dim and mass-dim (Renyi q=2) at modulus 2^k:

| k | 2^k | #active cyls | box-dim | mass-dim_q2 | avg local-dim |
|---|---|---|---|---|---|
| 5 | 32 | 16 | 0.800 | 0.793 | 0.797 |
| 6 | 64 | 32 | 0.833 | 0.816 | 0.825 |
| 7 | 128 | 64 | 0.857 | **0.832** ← max | 0.845 |
| 8 | 256 | 128 | 0.875 | 0.823 | 0.851 |
| 9 | 512 | 256 | 0.889 | 0.797 | 0.848 |
| 10 | 1024 | 512 | 0.900 | 0.756 | 0.838 |
| 11 | 2048 | 1024 | 0.909 | 0.720 | 0.826 |
| 12 | 4096 | 2048 | 0.917 | **0.671** ≈ Chang 0.68 | 0.809 |
| 13 | 8192 | 4096 | 0.923 | 0.625 | 0.790 |
| 14 | 16384 | 8192 | 0.929 | 0.582 | 0.769 |
| 15 | 32768 | 16375 | 0.933 | **0.544** | 0.747 |

**Mass dim declines monotonically from k=8 to k=15.** This is the signature of a **multifractal** measure — local Hausdorff dim varies across the support. The measure concentrates on a fractal subset whose dim depends on the scale at which we measure it.

The measure crosses Chang's reported H-dim 0.68 around **k=12**, which corresponds to scale 2^-12 ≈ 1/4096 in 2-adic. This is suggestive — Chang's "divergent set" H-dim 0.68 may be the asymptotic mass dim of OUR measure at fine scales. But the multifractal structure means there's no single δ; the conformal-measure framework (Sullivan, Pollicott-Urbański) which assumes constant δ does not apply directly.

**Box-dim asymptotes to 1.0 as k → ∞**, consistent with the support being all of Z_2 (every odd integer reachable from 1 in finite inverse-Syracuse steps, modulo conjecture).

## 4. σ-band conditional: q1 (near-root) reproduces D_emp; mechanism clarified

Partition tree by orbit-σ quartile (= inverse-tree depth quartile):

| Band | Depth range | n_nodes | Pearson vs D_emp at t=90 |
|---|---|---|---|
| q1 | 0–29 | 336,987 | **+0.863** |
| q2 | 29–43 | 319,236 | −0.164 |
| q3 | 43–60 | 326,627 | +0.439 |
| q4 | 60–152 | 327,285 | −0.109 |

Only q1 matches D_emp at t=90. q2 and q4 give near-zero or negative correlations. q3 weak positive.

**Cross-table: predicted band vs empirical D_t at multiple t values:**

| Pred band ↓ vs t → | t=10 | t=30 | t=50 | t=70 | t=90 | t=110 |
|---|---|---|---|---|---|---|
| q1 (0-25%, low σ) | −0.737 | +0.279 | +0.452 | +0.710 | **+0.863** | **+0.877** |
| q2 (25-50%) | +0.096 | +0.476 | +0.525 | +0.215 | −0.164 | −0.157 |
| q3 (50-75%) | −0.505 | **+0.655** | +0.525 | +0.435 | +0.439 | +0.487 |
| q4 (75-100%, high σ) | +0.148 | +0.355 | +0.460 | +0.332 | −0.109 | −0.067 |

**As t grows, the matching predictor band shifts from q3 (mid σ) at t=30 → q1 (low σ) at t=90/110.** This is exactly what the structural mechanism predicts:

### Mechanism

For an orbit starting at integer n with σ(n) = depth in inverse tree:
- Position at iteration t = node at depth (σ(n) − t) from n along path to 1 = node at depth (σ(n) − t) in inverse tree from m=1
- Surviving at t means σ(n) > t

For surviving orbits at t=90, most have σ slightly > 90, so position is at low inverse-tree depth (close to m=1, where the m_j cylinder lives). Hence q1 (low depth in inverse tree) reproduces D_emp at t=90.

For surviving orbits at t=30, surviving requires σ > 30; their position at t=30 has depth σ−30, which spans more variation. Best match is q3 (mid depth).

**This is the cleanest mechanistic interpretation of D_emp's structural origin: the trajectory measure at survivor-time t is the residue distribution over inverse-tree nodes at low depth.** Result 58's Pearson +0.86 falls out of this — variant (a) is dominated by low-depth nodes due to subtree-size weighting (most mass near root), so it converges to the q1 distribution.

## 5. Conformality check: REJECTED — mass ratios vary 0.92 to 2.34

For the measure to be **Sullivan-conformal** with Hausdorff dim δ, we need μ(T(A))/μ(A) ≈ |T'|^δ_2 where |T'|_2 = 2^(−v_2(3m+1)) for one Syracuse step.

Empirical mass ratios on cylinders mod 64:

| r | mass(r) | #integers | Syracuse(r) mod 64 | mass(image) | ratio |
|---|---|---|---|---|---|
| 5 | 3,415,302 | 63,153 | 1 | 3,158,618 | **0.925** |
| 7 | 1,922,334 | 24,623 | 11 | 2,081,780 | 1.083 |
| 13 | 1,917,981 | 63,160 | 5 | 3,415,302 | 1.781 |
| 21 | 1,351,308 | 64,879 | 1 | 3,158,618 | **2.337** |
| 27 | 2,086,058 | 15,410 | 41 | 1,840,522 | 0.882 |
| 31 | 1,549,431 | 10,278 | 47 | 1,594,407 | 1.029 |
| 37 | 1,569,409 | 58,128 | 7 | 1,922,334 | 1.225 |
| 53 | 1,777,240 | 63,238 | 5 | 3,415,302 | 1.922 |

Ratios range **0.92 to 2.34**. If Sullivan-conformal with constant δ, all ratios on cylinders with the same v_2(3m+1) value should match.

| r | v_2(3r+1) | mass ratio |
|---|---|---|
| 5 | 4 | 0.925 |
| 37 | 4 | 1.225 |
| 53 | 5 | 1.922 |
| 21 | 6 | 2.337 |
| 7 | 1 | 1.083 |
| 27 | 1 | 0.882 |
| 31 | 1 | 1.029 |

Within v_2=1 cylinders alone: ratios vary 0.88-1.08 (suggestive of constant). But across v_2: 0.93 (v=4) to 2.34 (v=6) — not consistent with 2^{−vδ} for any single δ.

**Sullivan-style conformal measure with constant Hausdorff dim δ: REJECTED.** The measure has finer multifractal structure than constant-δ conformality permits.

## 6. Verdict per brief outcomes

| Outcome | Status |
|---|---|
| (α) Z_2 measure variant reproduces D_avg within bootstrap | **PARTIAL** — Pearson +0.86 (Result 58 + variant b confirmed), residual ~30% at QSD-tilt extremes (r=5, r=23, r=13) |
| (β) Captures shape but quantitative gaps | **PRIMARY** — exactly this. Shape matched; quantitative gaps at QSD-tilt extremes |
| (γ) No variant reproduces D_avg | **REJECTED** — variant (a/b) gives strong structural match |
| Sullivan-conformal identification | **REJECTED** — measure is multifractal, not constant-δ conformal |

## 7. What this opens vs closes

**Closed:**
- Z_2 measure framework gives Pearson +0.86 with D_emp via subtree-size weighting (Result 58 + 59 variant b)
- σ-band conditional structure has clean mechanistic explanation (Section 4 here)
- Sullivan-conformal identification: REJECTED. Multifractal measure beyond Sullivan/Pollicott-Urbański constant-δ machinery
- Result 23 leading eigvec is the depth-asymptote density (different from value-truncation marginal D_emp)

**Open:**
- Esscher-tilt closure of the +0.86 → +0.95+ residual gap (Result 22 σ-quartile machinery)
- Multifractal spectrum f(α) characterization (full set of local dimensions across the measure's support)
- Whether multifractal nature relates to Chang's H-dim 0.68 as ONE point on the multifractal spectrum (mass-dim crosses 0.68 around k=12)

## 8. Connection to literature

The multifractal measure ≈ Bernoulli convolution analog on Z_2: a self-similar measure on the 2-adic integers with non-constant local dimension. Classical literature (Erdős, Solomyak, Peres-Schlag-Solomyak) on Bernoulli convolutions on R has exact analog questions. Whether the Collatz inverse-tree measure on Z_2 belongs to a parametric family (singular vs absolutely continuous, depending on parameters) is an open analytical question pointed at by this framework.

## 9. What this rules out

- Sullivan/Pollicott-Urbański conformal measure on Z_2 with constant δ: **REJECTED.**
- Any framework predicting constant Hausdorff dim across the support: **REJECTED.**
- The measure being captured by Result 23's leading eigvec: **REJECTED** (depth-asymptote object, different from value-truncation marginal).

## 10. What this concludes about the framework synthesis

The trajectory measure's structural identity:

> D_empirical(r) at survivor-time t ≈ residue distribution over nodes at depth (σ-t) in the inverse Collatz tree from m=1, marginalized to mod 32, plus an Esscher tilt at the QSD extremes (r=5, r=23, r=13).

This is the **integer-level inverse-tree marginal at value-truncation matching the empirical sampling regime**. The residue-only frameworks (Chang cylinder, Result 51; depth-truncated tree, Result 52) fail because they project away the value information. The Z_2 measure framework, when restricted appropriately to the value-uniform regime (this work), captures the structure.

The remaining gap (~30% at r=5/23/13) is the survivor-conditioning Esscher tilt, analytically tractable via Result 22's σ-quartile machinery.

## 11. Files

- `zadic_measure_framework.py` — variant (b), Hausdorff, σ-band, conformality
- `experiments_output/zadic_hausdorff.csv` — dim across k
- `experiments_output/zadic_measure_predictions.csv` — D_pred per variant
- `experiments_output/zadic_measure_log.txt` — full log
- `zadic_measure_framework.md` — this document (Result 59)

## 12. Concrete next moves

1. **Esscher tilt closure** (highest leverage): w(m) = subtree-size(m) × tilt(σ(m), q≈0.72). Test if Pearson lifts +0.86 → +0.95+. This is the analytical missing piece.
2. **Multifractal spectrum f(α)** computation. Sweep q ∈ [-2, 2] in Renyi dims; characterize the spectrum width and shape.
3. **σ-band fine partitioning**. Run with 8 or 16 σ bands rather than 4; verify the q1-at-large-t pattern across finer band granularity.
4. **N=2^26 confirmation**. Push to confirm Pearson plateau and the multifractal mass-dim curve are both stable.
