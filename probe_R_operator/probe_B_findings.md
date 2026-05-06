# Probe B — composition of R_k forcing operators

Composing R_(m) = R_7 ∘ ι_6 ∘ ... ∘ ι_4 ∘ R_4: V_4 → W_8, for m = 1..4.

ι_k: W_{k+1} ↪ V_{k+1} is the orthonormal embedding (P_{W_{k+1}}^T). Each composition step adds one renormalization level. The geometric mean σ_1(R_(m))^(1/m) is the candidate Lyapunov rate per level.

## Per-m results

| m | output level | shape | σ_1 | σ_1^(1/m) | σ_2 | rank |
|---|---|---|---|---|---|---|
| 1 | W_5 | 108×54 | 0.670003 | 0.670003 | 0.670003 | 36 |
| 2 | W_6 | 324×54 | 0.400491 | 0.632843 | 0.396706 | 36 |
| 3 | W_7 | 972×54 | 0.240309 | 0.621713 | 0.237720 | 36 |
| 4 | W_8 | 2916×54 | 0.143209 | 0.615166 | 0.141655 | 36 |

## Lyapunov rate of the composition

If the composition has a stable per-level multiplicative rate λ, then σ_1(R_(m))^(1/m) → λ as m grows.

| m | σ_1^(1/m) | closest reference | distance |
|---|---|---|---|
| 1 | 0.670003 | 0.6706 (single-level σ_1) (0.6706) | 0.0006 |
| 2 | 0.632843 | 0.6706 (single-level σ_1) (0.6706) | 0.0378 |
| 3 | 0.621713 | 0.6706 (single-level σ_1) (0.6706) | 0.0489 |
| 4 | 0.615166 | 0.6706 (single-level σ_1) (0.6706) | 0.0554 |

## Decay rate σ_1(R_(m+1)) / σ_1(R_(m))

Per-step σ_1 ratio (multiplicative gain when adding one level to the composition).

| m → m+1 | σ_1(R_(m)) | σ_1(R_(m+1)) | ratio |
|---|---|---|---|
| 1 → 2 | 0.670003 | 0.400491 | 0.597745 |
| 2 → 3 | 0.400491 | 0.240309 | 0.600035 |
| 3 → 4 | 0.240309 | 0.143209 | 0.595936 |

## Verdict

σ_1^(1/m) sequence: 0.6700, 0.6328, 0.6217, 0.6152 (decreasing toward asymptote)
σ_1(R_(m+1)) / σ_1(R_(m)) ratio: 0.5977, 0.6000, 0.5959 (per-step multiplicative gain stable at ~0.60)

**Brief's third walk-back gate fires.** Composition does NOT reveal ρ_slow ≈ 0.83.

What the data shows:

1. **Per-step multiplicative gain stabilizes at ≈ 0.60.** Each new R added to
   the composition multiplies σ_1 by ~0.60. So σ_1(R_(m)) → σ_1(R_(1)) · 0.60^(m-1)
   asymptotically, and σ_1^(1/m) → 0.60 as m → ∞.

2. **Per-level rate ≈ 0.60, not 0.83.** The asymptotic Lyapunov rate of the
   single-direction composition is around 0.60, an entirely different number from
   the order-3 recurrence's ρ_slow ≈ 0.827.

3. **Rank stays at 36 = 2·n_4/3** through composition. Adding more R_k's to
   the chain doesn't expand the effective subspace — the composition operates
   on the same 36-dim non-kernel of R_4, with subsequent R's projecting it
   smaller in W-magnitude but not changing its dimension.

4. **No new structure from composition.** σ_2/σ_1 ratio is essentially 1 at
   every m (e.g., m=4: σ_2 = 0.1417 vs σ_1 = 0.1432, ratio 0.989) — the
   2-fold pairing seen at single-level persists through composition.

**Implication.** The order-3 recurrence's ρ_slow ≈ 0.827 (R²=0.797 fit on
ε_2..ε_11) is NOT in σ_1 of R^k composition. Three possibilities:

- ρ_slow lives in a non-σ_1 SVD direction, i.e., ε_k corresponds to a
  specific input/output functional that's not the dominant singular pair
  but a different mode whose decay is slower (0.83 > 0.60).
- ρ_slow is a finite-k recurrence-fit artifact (recall R² = 0.797 on
  10 points, not high), and the true asymptotic rate is closer to 0.60
  matching this composition probe.
- ε_k decay involves an additional structure outside R-composition (e.g.,
  forcing from a specific initial condition aligned with V → W^c or
  a non-linear effect).

**Walk-back gate disposition (full chain of probes):**

- Walk-back gate "σ_1 ≈ 1": NOT triggered at any single k or composed m.
- Walk-back gate "σ_1 ≈ 0": NOT triggered (rank stays 36, σ_1 stays > 0.14).
- "Cluster but no match": TRIGGERED at single-level (σ_1 ≈ 0.67) AND at
  composed level (per-level rate ≈ 0.60). Two stable clusters, neither
  matching ρ_slow.

The probe establishes that R_k composition is well-defined and produces
asymptotic per-level rate ~0.60, but does not identify the operator
structure responsible for ρ_slow ≈ 0.83.

## Files

- `probe_B_composition.csv` — per-m sigma values
- `probe_B_findings.md` — this writeup