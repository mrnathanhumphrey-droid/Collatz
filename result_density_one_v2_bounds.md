# Density-1 v_2 bound tests

Quantitative density-1 confirmation tests connecting the eps_k Syracuse-Markov framework to the Lagarias / Tao 2-adic equidistribution prediction.

**Data**: `data/v_seq_N8388608.parquet` — Syracuse trajectories for starting integers in [3, 8388607]. Filtered to odd starts coprime to 3 ⇒ 2,796,202 qualifying trajectories.

**Trajectory length**: max 248, median 51, mean 53.32. (Per-trajectory k=1000/10000 tests from the brief unreachable; adapted to length-binned density-1 convergence.)

**log_2(3)** = 1.5849625007

**Mean of v across all trajectories**: 2.102161 (geometric prediction: 2.0)

## TEST A: density of v >= k vs geometric null 2^{-(k-1)}

Geometric null: under Lagarias' 2-adic equidistribution, v_2(3 n_i + 1) is asymptotically Geom(1/2) on {1, 2, ...}, so P(v >= k) = 2^{-(k-1)}. Empirical density of v >= k along each trajectory is computed and compared.

| k | null = 2^{-(k-1)} | emp mean | emp median | frac >= null | frac > null |
|---|------------------:|---------:|-----------:|-------------:|------------:|
| 1 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 0.000000 |
| 2 | 0.500000 | 0.526301 | 0.507937 | 0.568757 | 0.512438 |
| 3 | 0.250000 | 0.287930 | 0.265306 | 0.608913 | 0.578633 |
| 4 | 0.125000 | 0.155074 | 0.138462 | 0.642110 | 0.619087 |
| 5 | 0.062500 | 0.065202 | 0.054054 | 0.400843 | 0.385648 |
| 6 | 0.031250 | 0.031350 | 0.023256 | 0.371930 | 0.362452 |
| 7 | 0.015625 | 0.016533 | 0.009804 | 0.400354 | 0.395704 |
| 8 | 0.007812 | 0.009187 | 0.000000 | 0.278616 | 0.278550 |
| 9 | 0.003906 | 0.005060 | 0.000000 | 0.155537 | 0.155537 |
| 10 | 0.001953 | 0.003226 | 0.000000 | 0.092712 | 0.092712 |

### TEST A by trajectory length bin

Fraction of trajectories with empirical density of v >= k exceeding the geometric null 2^{-(k-1)}:

| length bin | n_traj | k=1 | k=2 | k=3 | k=4 | k=5 | k=6 | k=8 | k=10 |
|---|---|---|---|---|---|---|---|---|---|
| [1,5) | 312 | 1.0000 | 0.9968 | 1.0000 | 1.0000 | 0.9872 | 0.9647 | 0.9327 | 0.7596 |
| [5,10) | 6,943 | 1.0000 | 0.9912 | 0.9996 | 0.9991 | 0.9967 | 0.9806 | 0.8541 | 0.5218 |
| [10,20) | 127,902 | 1.0000 | 0.9852 | 0.9984 | 0.9978 | 0.9503 | 0.9343 | 0.6102 | 0.2800 |
| [20,30) | 349,940 | 1.0000 | 0.9659 | 0.9887 | 0.9820 | 0.8642 | 0.8644 | 0.4129 | 0.1635 |
| [30,50) | 846,831 | 1.0000 | 0.8115 | 0.9022 | 0.8613 | 0.5654 | 0.4308 | 0.2917 | 0.1021 |
| [50,100) | 1,354,878 | 1.0000 | 0.3176 | 0.3376 | 0.4321 | 0.1550 | 0.1805 | 0.2127 | 0.0536 |
| [100,200) | 109,334 | 1.0000 | 0.0149 | 0.0030 | 0.0193 | 0.0075 | 0.0139 | 0.1379 | 0.0292 |
| [200,1000) | 62 | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1774 | 0.3548 |

## TEST B: trajectory mean v vs log_2(3)

Density-1 claim: for density-1 of starting integers, (1/L) Σ v_2(3 n_i + 1) > log_2(3) for L large. Per-trajectory mean v is computed, and fraction exceeding log_2(3) reported by trajectory length bin.

| length bin | n_traj | mean(mean_v) | std(mean_v) | frac > log_2(3) | mean excess | median excess |
|---|---|---|---|---|---|---|
| [1,5) | 312 | 6.887821 | 2.144668 | 1.000000 | +5.302858 | +5.165037 |
| [5,10) | 6,943 | 4.172558 | 0.581923 | 1.000000 | +2.587596 | +2.526149 |
| [10,20) | 127,902 | 2.939443 | 0.269971 | 1.000000 | +1.354480 | +1.297390 |
| [20,30) | 349,940 | 2.455911 | 0.116200 | 1.000000 | +0.870948 | +0.859482 |
| [30,50) | 846,831 | 2.149077 | 0.089537 | 1.000000 | +0.564114 | +0.561379 |
| [50,100) | 1,354,878 | 1.916532 | 0.061301 | 1.000000 | +0.331570 | +0.328081 |
| [100,200) | 109,334 | 1.782506 | 0.021242 | 1.000000 | +0.197544 | +0.200752 |
| [200,1000) | 62 | 1.693428 | 0.004667 | 1.000000 | +0.108465 | +0.109212 |

## Verdict

### TEST B: tautological for terminating trajectories — not an empirical confirmation

The 100% pass rate is **mathematically forced** for any Syracuse trajectory
that reaches 1. Direct algebra:

For trajectory n_0 → n_1 → ... → n_L = 1 with each step
n_{i+1} = (3 n_i + 1) / 2^{v_i}:

  log_2(n_{i+1}) − log_2(n_i) = log_2(3 + 1/n_i) − v_i
                              = log_2(3) − v_i + log_2(1 + 1/(3 n_i))

Summing over the trajectory:

  −log_2(n_0) = L · log_2(3) − L · mean_v + Σᵢ log_2(1 + 1/(3 n_i))

  ⇒  **mean_v − log_2(3) = (1/L) · [log_2(n_0) + Σᵢ log_2(1 + 1/(3 n_i))]**

Both terms on the right are non-negative for n_0 ≥ 3 and n_i ≥ 1. So
mean_v > log_2(3) **identically** for every Syracuse trajectory that
terminates at 1.

The 100% rate of 2,796,202 trajectories passing TEST B is therefore a
restatement of "all 2,796,202 trajectories reached 1" — already established
on this range from prior Collatz verification. It is not new empirical
content about Lagarias's prediction.

### What IS empirically informative from this probe

1. **Unconditional ensemble mean v_2 = 2.102** vs geometric Geom(1/2)
   prediction 2.0. The 5% deviation is modest and consistent with the
   asymptotic geometric marginal (Tao 2019) at this finite ensemble size.

2. **Mean excess decays as ~ log_2(n_0)/L by length bin.** Empirical
   excess: +1.35 at L∈[10,20), +0.56 at L∈[30,50), +0.20 at L∈[100,200),
   +0.11 at L∈[200,1000). With typical n_0 ≲ 2²³, log_2(n_0)/L ≈ 22/L,
   matching the by-bin numbers within the small Σᵢ log_2(1 + 1/(3 n_i))
   correction. This is the algebra above working out — also not new
   empirical content beyond the algebraic identity.

### TEST A: per-trajectory geometric null fails on long trajectories — selection effect

The pointwise statement "density of {i : v_i ≥ k} ≥ 2^{-(k-1)}" **does
not hold** in the strict per-trajectory sense once trajectories are long.
Pass rates at k=2 by length bin:

| length bin | n_traj | frac(density(v≥2) ≥ 1/2) |
|---|---|---|
| [10,20) | 127,902 | 0.985 |
| [30,50) | 846,831 | 0.812 |
| [50,100) | 1,354,878 | 0.318 |
| [100,200) | 109,334 | 0.015 |
| [200,1000) | 62 | 0.000 |

Mean_v by bin (TEST B table) drops from 2.94 at L∈[10,20) to 1.69 at
L∈[200,1000) — approaching log_2(3) from above. Selection on "trajectory
takes L >> 1 steps to reach 1" forces the v-distribution toward small
values: from the algebraic identity, mean_v − log_2(3) = (log_2(n_0) +
small)/L, so mean v is necessarily near log_2(3) when L is large relative
to log_2(n_0). The conditional v-distribution given large L is not
Geom(1/2); its mean is structurally constrained, and cumulative density
P(v ≥ k) drops below 2^{-(k-1)} for k ≥ 2.

### Positioning vs Lagarias

The empirically informative content is the unconditional ensemble mean
v = 2.102 ≈ 2 (Geom(1/2)). This is **consistent with** Tao 2019's
measure-theoretic Geom(1/2) prediction at strong sample resolution. It
is not stronger evidence than the measure-theoretic statement — the
finite-sample empirical confirmation is what one expects given the
asymptotic theorem. The tautological TEST B is **not** a quantitative
density-1 result strengthening the connection to Lagarias.

## Notes

- Per-trajectory CSV (50,000-row uniform random sample) at result_density_one_v2_bounds.csv.
- Diagnostic file documents failure-mode analysis (residue patterns, length distribution of failures).
- Trajectory length is the natural per-trajectory k; the brief's k=10000 sample point is unreachable per-trajectory in this data (max length = 248).
