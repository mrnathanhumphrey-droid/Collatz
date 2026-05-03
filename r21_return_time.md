# Return-time distribution to r ≡ 21 mod 32 (Result 47)

**Status.** Mixed outcome (b)+(d). Renewal structure holds at the autocorrelation
level (ρ_lag1 ≈ −0.005), but two additional structures matter:

1. **Band-mixing**: marginal P(G) is mixture of band-conditional geometric-like
   distributions with rates λ(q) varying 3× across σ-bands.
2. **Gap → cylinder coupling**: ρ(G_n, V_n) = −0.139. Longer waits select
   different m mod 2^k residues at the next r=21 visit.

The cylinder formulation is structurally close-but-not-quite-sufficient. The
trajectory measure decomposes as:

  μ_traj = (band selection) × ∏_n (band-conditional gap) × (gap-conditional cylinder draw)

Three local components, all well-defined. The Lagarias-class open piece
sharpens to "characterize gap-conditional cylinder draw at r=21 visits."

## Step 2: Marginal return-time distribution

250K orbits at N=2³⁴, 1.12M r=21 visits, 875K inter-visit gaps:

| statistic | value |
|---|---:|
| ⟨G⟩ | 16.37 |
| SD(G) | 15.42 |
| median | 11 |
| skew | 1.48 |
| excess kurt | 2.57 |

Linear fit log P(G=g) = −2.655 − 0.0743·g, **R² = 0.988** (geometric tail).
Implied λ = 0.074 → predicted Geom mean = 13.5 (off from empirical 16.37 by 21%).

Var/Mean² empirical = 0.888 vs Geom prediction 0.939 → **under-dispersed** vs
pure geometric. KS to Geom(1/⟨G⟩) = **0.033** (fails 0.02 strict threshold).

The marginal is geometric-like in tail decay but quantitatively distinct —
mixture of geometrics with band-dependent rates.

## Step 3: Per-σ-band conditional distributions

| band | n | ⟨G\|band⟩ | SD | λ_band ≈ 1/⟨G⟩ |
|------|---:|---------:|------:|-----------:|
| 0–25 | 198,351 |  8.92 |  7.94 | 0.112 |
| 25–50 | 207,495 | 13.93 | 12.46 | 0.072 |
| 50–75 | 220,740 | 17.97 | 15.36 | 0.056 |
| 75–95 | 192,172 | 22.02 | 17.97 | 0.045 |
| 95–100 | 55,770 | 26.09 | 21.13 | 0.038 |

**3× variation in rate across bands.** Low-σ orbits visit r=21 frequently;
high-σ orbits spend more steps in the deterministic 31-residue spectrum
between visits. Mechanism: high-σ orbits accumulate more steps per
log(n)-decay rate, more total steps total, more steps between r=21 visits.

Within each band, var/mean² ranges 0.66–0.79 — under-dispersed relative
to geometric prediction (0.89–0.96). Gaps are TIGHTER than pure geometric
within bands, suggesting some structural mechanism beyond i.i.d. memoryless
returns.

## Step 4: Autocorrelation — renewal structure intact

| pair | n | correlation |
|------|---:|---:|
| (G_n, G_{n+1}) | 627,272 | **−0.0054** |
| (G_n, G_{n+2}) | 394,706 | +0.0154 |

**Gaps are essentially independent** (lag-1 |ρ| < 0.01). No memory across
visits. Renewal structure holds at the gap level.

## Step 5: Gap-cylinder coupling — moderate

| pair | n | correlation |
|------|---:|---:|
| (G_n, V_n) | 874,528 | **−0.1385** |
| (V_n, G_{n+1}) | 627,272 | +0.0221 |

**Asymmetric coupling.** The gap leading TO visit n correlates with v_2 AT
visit n (ρ = −0.14). Longer waits → smaller v_2. But v_2 at current visit
does NOT predict the next gap (ρ ≈ 0). One-directional information flow.

Mechanism: long gap means orbit traversed many deterministic residues, accumulating
log m drift via known deterministic v values. The specific mod-2^k bit pattern
when arriving at r=21 depends on the gap length (which determines the cumulative
m mod 2^k transition through the deterministic chain). Long gaps → more cumulative
mixing → m's higher bits more "uniform-like" at arrival → more likely to land
on residues where v_2 is small (the typical case).

This is consistent with Result 42's finding that orbit history barely matters
beyond m mod 2^k — but the gap LENGTH gates which m mod 2^k residue you arrive at.

## Step 6: First-passage T_1 by starting residue

| r_0 mod 64 | ⟨T_1⟩ | median | P(T_1=0) | group |
|----------:|------:|-------:|---------:|-------|
| 21 |  0.00 |  0 | 1.000 | r=21 (start) |
| 53 |  0.00 |  0 | 1.000 | r=21 (start, 53 mod 32 = 21) |
| 7  | 16.72 | 12 | 0.000 | Chang I_2 |
| 27 | 18.82 | 14 | 0.000 | Chang I_2 |
| 31 | 17.87 | 13 | 0.000 | Chang I_2 |
| 59 | 17.34 | 12 | 0.000 | Chang I_2 |
| 63 | 19.85 | 15 | 0.000 | Chang I_2 |
|  3 | 16.91 | 12 | 0.000 | non-I_2 v=1 |

Chang I_2 residues have T_1 ≈ 16.7–19.9, similar to non-I_2 residues
(T_1 ≈ 16.9 for r=3). No clear spectral contraction visible at first-passage
time. Chang's I_2 contraction (ρ(B̃_2^ext) ≤ 5/32) operates on different
observable than first-passage to r=21.

## Verdict — mixed (b) + (d)

**(b) confirmed**: Renewal within bands but band-mixing in marginal.
λ(q) varies 3× across σ-bands. Marginal P(G) is mixture, not pure geometric.

**(d) flagged**: G→V coupling ρ = −0.139 means cylinder draw at r=21 is
gap-conditional, not gap-independent.

**Renewal property partially holds**: gaps are i.i.d. across visits within
band (autocorrelation = −0.005). But the cylinder draw is coupled to
incoming gap, so the trajectory measure is NOT (band × renewal × independent cylinder).

## Sharpened decomposition

```
μ_traj = (band selection P(q|N))
       × ∏_n (band-conditional gap G_n ~ near-geometric(λ(q)))
       × (gap-conditional cylinder draw V_n | G_n at r=21)
```

Three local components. Each well-defined and testable.

The Lagarias-class open piece: characterize the **gap-conditional cylinder
draw V_n | G_n** at r=21 visits. This is more specific than Result 42's
"visit-frequency on m mod 2^k" formulation — the visit-frequency is itself
G-conditional.

## What this means for the open piece

The decomposition gives a more granular formulation:

> Open piece = joint distribution of (G_n, V_n) at r=21 visits, given σ-band.
> G_n ~ near-geometric(λ(band)) with band-specific rate.
> V_n given G_n = arithmetic v_2 of (3·m+1)/2 where m's higher bits depend on
> the residue path traversed during the gap.

The cylinder formulation (Result 42) captures the marginal V_n distribution.
The renewal-coupled formulation (this Result) captures the joint (G_n, V_n)
distribution. Both are valid; the joint is finer-grained.

## For Chang correspondence

> Empirical renewal structure at r=21 visits: gaps are independent across
> visits (ρ_lag1 ≈ −0.005), with band-conditional near-geometric distribution
> (rate λ(q) varies 3× across σ-bands). G→V coupling ρ = −0.139 (asymmetric:
> long gaps select smaller v_2 at next visit). Chang I_2 residues show no
> distinctive first-passage behavior (T_1 ≈ 17–20 like other non-21 residues).
> The cylinder formulation decomposes cleanly into (band × renewal × G-conditional
> cylinder draw) — three local components.

## Files

- `experiments/74_r21_return_time.py`
- `experiments_output/74_r21_return_time_log.txt`
- `experiments_output/74_r21_return_time_marginal.csv`

Compute: 1.2s walking + analysis (250K orbits, 1.12M visits, 875K gaps).
