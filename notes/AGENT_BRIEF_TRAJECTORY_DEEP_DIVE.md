# Collatz Project — Context Dossier (2026-05-02)

This is a state-of-the-project snapshot. Not a task assignment. Not a brief
to a specific agent. The user is the conductor; this file exists only as
context they may choose to share.

## Project in one paragraph

`C:\Collatz\` contains a research project on the residue-class structure of
the Collatz total stopping time σ(n) (the 3x+1 problem) and a generalization
to qx+1 for q ≥ 5. The 3x+1 thread has produced the **prefix decomposition
theorem** (defined below). The qx+1 thread has produced an empirical Cramér
large-deviation law for convergence rates with θ(q) the root of
`q^(−θ) = 2^(1−θ) − 1`, matching at q=5 to 0.01%. The project is in
correspondence with N. Bonacorsi (Columbia) on a potential joint paper.

## Definitions and conventions

### Prefix decomposition (operational)

For odd n, write `n = 2^k · m + r` where `r = n mod 2^k` is an odd integer in
`[1, 2^k − 1]` and `m = (n − r) / 2^k` is the integer "tail". Track the
Collatz map symbolically on the state `state = a·m + c`, starting from
`(a, c) = (2^k, r)`. Apply the rules:

- If `a` is even: parity of state is determined by `c` (independent of `m`).
  - If `c` is even: state is even. Halve. `(a, c) → (a/2, c/2)`.
  - If `c` is odd:  state is odd. Apply 3n+1. `(a, c) → (3a, 3c + 1)`.
- If `a` is odd: parity of state depends on `m`. **Terminate.**

The number of iterations is `prefix_steps(r)` (typically 7–16 for k=6).
Final state is `a_final(r) · m + c_final(r)`. By construction `a_final(r) = 3^j`
for some `j ∈ {1, ..., k}` — `j` is the count of "3a" rule applications in
the prefix.

The prefix is fully deterministic per residue class. After it terminates,
all further dynamics depend on `m`.

The heuristic prediction for the per-class intercept of σ vs ln(n) is
```
α_det(r) = prefix_steps(r) + K · ln(a_final(r) / 2^k)
```
where `K = 3 / (ln 4 − ln 3) = 10.4282...` is the odd-n drift constant.

### Noise-floor ratio (operational)

Let `α_actual(r)` be the per-class OLS intercept from regressing σ on ln(n)
on data restricted to n ≡ r mod 2^k. Let `SE_α(r)` be the standard error of
that intercept estimate. Fit `α_actual(r) ~ a + b · α_det(r)` across
all classes, get residuals `resid(r) = α_actual(r) − (a + b · α_det(r))`.

```
noise_floor_ratio = SD(resid) / mean(SE_α)
```

Ratio ~1: residuals are at sampling-noise floor (decomposition explains all
detectable structure). Ratio ≫ 1: real residual structure beyond prefix
prediction. Ratio < 1: prefix decomposition over-explains relative to
sampling noise (k=4 in the existing data).

### Trajectory measure on v

The natural-density measure on odd integers gives `v = ν₂(3m+1) ~ Geom(1/2)`,
i.e., `P(v = k) = 2^(−k)` for k = 1, 2, 3, ....

The trajectory measure pools v values along iterates of the Syracuse map
`T(m) = (3m+1) / 2^v`. The trajectory measure differs from natural density:
v=4 and v=10 are over-represented (~+23%, +24%); v=6, 7, 8, 9, 11, 12, 13, 14
are under-represented (graded sags).

For v=k *exactly*, the residue `r_v = m mod 2^(k+1)` is a single fixed value
solving `3m + 1 ≡ 2^k (mod 2^(k+1))`. Concretely:
- v=4: m ≡ 5 mod 32 (= 2^5)
- v=10: m ≡ 341 mod 2048 (= 2^11)
- v=16: m ≡ ? mod 131072 (= 2^17)  [computed by `m ≡ 3⁻¹ · (2^16 − 1) mod 2^17`]

### Cramér large-deviation law (qx+1 thread)

For q ≥ 5, the qx+1 map is transient (drifts away from 1). Per-residue-class
convergence rate decays as
```
conv_rate(j; q) ≈ A(q) · exp(−θ(q) · log(q) · j)
```
where j is the prefix odd-step count and θ(q) is the unique positive root of
```
q^(−θ) = 2^(1−θ) − 1
```

The derivation uses Cramér's theorem applied to a random walk in
log-coordinates with steps `X = log(q) − v · log(2)` and v ~ Geom(1/2). The
quantity that enters Cramér's theorem is the moment generating function
```
M(s) = E[2^(s·v)]
```
evaluated at `s = θ − 1`. Under Geom(1/2):
```
M(s) = sum_{k=1}^∞ 2^(s·k) · 2^(−k) = 1 / (2^(1−s) − 1)
```
The Cramér root condition reduces to `M(θ − 1) · q^(−θ) = 1`, equivalent to
the quoted equation `q^(−θ) = 2^(1−θ) − 1`.

**This is what the trajectory-measure deviation could perturb.** Variance is
*not* the right test — what matters for Cramér's theorem is whether the
empirical trajectory MGF `M_traj(s) = Σ_v 2^(s·v) · P_traj(v)` at `s = θ − 1`
matches the Geom(1/2) MGF at the same s. A first-moment-skewed deviation can
preserve variance but shift the MGF.

## What's been done this session

### Experiment 24 — k-sweep alpha decomposition (done)

Sweep of the prefix decomposition across `k ∈ {4, 5, 6, 7, 8, 9, 10, 11, 12}`
on N=2^27 data (67M odd values). The noise-floor ratio sits in [0.90, 0.99]
across all k except k=4 (ratio = 0.70 — prefix over-explains at coarse
resolution). a_final levels = k exactly at every tested k (i.e., 2^(k−1)
odd classes collapse to k distinct distributions). Result confirms the
universality claim at finer modular grids than the original writeup tested.

CSV: `experiments_output/24_k_sweep_alpha_decomposition_N134217728.csv`

### Experiment 25 — trajectory measure of v (done)

N_start = 10^8 odd starts, T = 200 Syracuse steps each, 7×10^9 v-samples
pooled. First-step v matches Geom(1/2) to 4 decimal places (sanity ✓ on
natural density). Trajectory v deviates structurally:

| v | traj_ratio |
|---|---|
| 1 | 1.0001 |
| 2 | 0.972 |
| 3 | 0.996 |
| 4 | **1.232** |
| 5 | 0.966 |
| 6 | 0.823 |
| 7 | 0.802 |
| 8 | 0.754 |
| 9 | 0.720 |
| 10 | **1.237** |
| 11–14 | 0.61, 0.67, 0.59, 0.46 |
| 15 | 0.93 |
| 16 | 0.879 |
| 22 | 0.389 |
| 28 | 0.422 (high SE) |

The hypothesis "every v ≡ 4 mod 6 spikes" was rejected — only the first two
members of that residue class are spikes; v=16 onward sits inside the general
decay band. Mechanism for the v=4 / v=10 specificity is open.

CSV: `experiments_output/25_trajectory_measure_Nstart100000000_T200.csv`

### Experiment 26 — μ_β at large N (in progress)

Streaming OLS of σ vs ln(n) on odd n ∈ [3, N] for many N values. Result:

| N | log₂(N) | β | gap from heuristic |
|---|---|---|---|
| 2^20 | 20 | 10.3723 | +0.0559 |
| 2^22 | 22 | 10.3816 | +0.0466 |
| 2^23 | 23 | 10.3845 | +0.0437 |
| 2^24 | 24 | 10.4044 | +0.0238 |
| 2^25 | 25 | 10.4191 | +0.0091 |
| 2^26 | 26 | 10.4192 | +0.0090 |
| 2^27 | 27 | 10.4293 | **−0.0011 (crossed)** |
| 2^28 | 28 | 10.4298 | −0.0016 |
| 2^29 | 29 | 10.4252 | +0.0030 |
| 2^30 | 30 | 10.4236 | +0.0045 |
| 2^31 | 31 | 10.4213 | +0.0069 |
| 2^32 | 32 | 10.4187 | +0.0095 |

β is **non-monotone**. Approaches heuristic from below for N ≤ 2^26, jumps
above between 2^26 and 2^27, peaks at 2^28, then drifts back below the
heuristic and re-opens the gap to ~+0.010 by N = 2^32. The "monotone
approach from below" claim in the current writeup is wrong; the writeup A1
needs revision to reflect oscillatory behavior.

Cross-check: re-running streaming OLS at N=2^25 reproduces the existing
writeup's β = 10.4191 to 5 decimals — the non-monotone behavior is not a
method artifact.

## Files of record

- `writeup.md` — main 3x+1 writeup. Currently being edited by the primary
  agent (post-review fixes); A1 wording needs further revision in light of
  the non-monotone N-trend.
- `findings.md` — append-only audit log (469+ lines, latest entries on B1
  Pathfinder gap and qx+1 closure).
- `literature_check.md` — prior-art audit identifying Terras 1976 Lemma 4 as
  the asymptotic predecessor of the prefix decomposition.
- `project_collatz.md` and `project_collatz_qx1.md` (in user's auto-memory
  at `C:/Users/Nate/.claude/projects/c--As-Above-So-Below-Master/memory/`).

## Available data

- `data/main_N{N}.parquet` for `N ∈ {2^20, 2^22, 2^23, 2^24, 2^25, 2^27}`.
- `data/q_main_q{q}_N{N}.parquet` for q ∈ {3, 5, 7, 9, 11, 13} at varying N
  up to 10^9 for q=7 and q=11.
- σ caches at N ∈ {2^26, 2^28, 2^29, 2^30, 2^31, 2^32} were built in-memory
  by experiment 26 but not persisted to parquet.

## Open questions (phrased as questions, not assignments)

These are unresolved threads. Order is not implied. None of these has been
assigned to anyone.

1. **MGF transfer for the qx+1 Cramér result.** Compute the empirical
   trajectory MGF `M_traj(s) = Σ_v 2^(s·v) · P_traj(v)` from the experiment
   25 CSV at `s = θ(q) − 1` for q=5, 7, 9, 11. Compare to Geom(1/2) MGF.
   Magnitude of the deviation determines whether the q=5 0.01% Cramér match
   is fully protected or needs a small correction. (Variance check is a
   weaker form of this and can mask first-moment skewed deviations.)

2. **Mechanism of the v=4, v=10 spikes.** The trajectory measure pushforward
   on m mod 2^(k+1) tests whether specific residue classes are
   over-represented:
   - v=4 → m ≡ 5 mod 32. If trajectory measure on m mod 32 over-weights
     residue 5 by ~1.23×, the v=4 spike is mechanistically explained as a
     consequence of Syracuse's invariant measure favoring this class.
   - v=10 → m ≡ 341 mod 2048. Same test at finer resolution.
   - v=16 → m ≡ r mod 131072, where r solves 3m+1 ≡ 2^16 mod 2^17. If the
     over-weight DECAYS at higher v (consistent with the empirical
     ratio 0.879 at v=16), that explains why the spike pattern truncates.

3. **q=5 trajectory v-distribution shape.** Whether the spike pattern is
   q=3-specific or has a q-shifted analogue. Pooling rule matters: long
   convergent orbits at q=5 are rare; pooling all v's weights long
   trajectories more, while per-trajectory equal weighting or fixed-T
   truncation gives different measures. The choice should be specified
   before running.

4. **Cause of the β oscillation between N=2^26 and N=2^32.** Three candidates:
   (i) record-σ outliers in [2^26, 2^27] pulling β up disproportionately,
   testable by truncating top-K σ values per range; (ii) a log(log(n))
   correction term in `E[σ | n]`, testable by fitting
   `σ ~ a + b · ln(n) + c · ln(ln(n))`; (iii) connection to the v=4, v=10
   trajectory deviation if its magnitude is N-dependent. The first is the
   cheapest test.

5. **Neighbor coincidence at q=3.** From a separate measurement at k=10,
   ~25% of consecutive m's in the same residue class produce identical σ.
   Reduces by prefix algebra to: σ(post-state) = σ(post-state + a_final) for
   ~25% of consecutive m's. Constraint on the stochastic remainder S(n).
   Not characterized.

## House conventions

- Numbered experiments (`NN_description.py`) following the existing pattern
- chunk-parallel numba for histogram/reduction (the naive
  `tid = i % n_threads` pattern has a race condition; see experiment 25 for
  the working chunked template)
- Sanity-check protocol gate before any "structural" / "anomaly" claim:
  uniform-sampling baseline, definition recheck, finite-N guard, parity /
  off-by-one check, numerical-precision check (in
  `feedback_sanity_check_protocol.md`)
- Audit-trail-complete documentation: every empirical finding gets recorded
  in `findings.md` with the boring explanations that were ruled out
- No editorializing in writeup or one-sheet text — numbers and what they
  mean
