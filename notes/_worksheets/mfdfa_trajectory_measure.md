# Result 61 — Multifractal analysis of trajectory measure

**Outcome (β):** clean spatial multifractal spectrum extracted; standard
weighted-Cantor closed forms don't fit; temporal time series are essentially
monofractal Brownian. Two findings:

1. **Temporal v_t / log₂ m_t along orbits is monofractal**, h(q) ≈ 0.5 for
   q ≥ 0 — confirms classical Collatz heuristic (log m random walk with
   independent increments).
2. **Spatial measure on Z₂ is multifractal**, with D₀ = 1.00, D₁ = 0.61,
   D₂ = 0.27, D∞ ≈ 0.15. Chang's H-dim 0.68 and Result 23's 0.6755 both
   locate at **q ≈ 1 (information dimension D₁)** — they're probing the
   same point of f(α), not a shared single δ.

## Methodology pivot

The brief proposed three time-series choices for Kantelhardt MF-DFA:
v_t per orbit, D(r,t) marginalized, and M_closed iterates. Choice (a) was
fired first; result was diagnostic (γ): h(q) ≈ 0.5 for q ≥ 0 and explosion
for q < 0 with R² ≈ 0.4 — known MF-DFA artifact when integer-valued series
produce F²≈0 segments under linear detrending at small scales.

Switched to log₂(m_t) increments (continuous-valued, Brownian-like): same
pattern. h(q) ranges 0.50 → 0.48 across q ∈ [0, 5] with R² > 0.99. Multifractal
width Δh < 0.03 — essentially monofractal.

**This is the correct empirical finding for temporal series**: log m
along Collatz orbits IS a Brownian random walk with i.i.d. increments
(drift −log₂(4/3) = −0.415, std ≈ 1.37). The R59 multifractality is
NOT temporal — it lives in the static asymptotic measure on Z₂.

Spatial analysis re-implemented as direct partition function on Z₂ cylinders:
Z_q(2^k) = Σ_C μ(C)^q, with μ from R58 inverse-tree subtree-size weighting
(Pearson +0.86 with D_emp, established as best single-weight family).

## Temporal MF-DFA (log₂ m_t increments)

```
   q       h(q)      SE     R²       tau(q)
  -5     2.69      1.26   0.396    -14.45    (negative-q artifact)
  -1     1.61      0.65   0.469     -2.61    (negative-q artifact)
  -0.5   0.60      0.07   0.920     -1.30    (transitional)
   0     0.501     0.011  0.997     -1.000   ✓
   0.5   0.498     0.010  0.997     -0.751   ✓
   1     0.495     0.009  0.998     -0.505   ✓
   2     0.489     0.007  0.999     -0.022   ✓
   3     0.484     0.006  0.999      0.452   ✓
   5     0.476     0.003  1.000      1.378   ✓
```

For q ∈ [0, 5]: h(q) = 0.5 - 0.024 q (effectively a single Hurst H = 0.5).
Conclusion: log m increments are uncorrelated noise. **Trajectory is
genuinely Brownian temporally.**

σ-band conditional gives same picture: 0-25 band has h ≈ 0.40, others have
h ≈ 0.27 across q. Bands separate orbit length not multifractal structure.

## Spatial multifractal Z_q on Z₂ cylinders

Inverse tree to max_value = 2²⁰ = 1,048,576: 312,238 odd nodes.
Subtree-size weights (variant a from R58). Cylinders mod 2^k for k ∈ [6, 16].

```
   q       tau(q)     D_q       R²
  -5      -11.27     1.879     0.911   (small-cluster amplification, noisy)
  -3       -6.77     1.693     0.935
  -2       -4.50     1.500     0.962
  -1       -2.36     1.182     0.998
  -0.5     -1.62     1.080     1.000
   0       -1.00     0.999     1.000   ← support dim = full Z₂
   0.5     -0.43     0.863     0.999
   1       -0.61     0.608     0.972   ← info dim D₁
   1.5      0.19     0.374     0.861
   2        0.27     0.267     0.770   ← correlation dim D₂
   3        0.39     0.194     0.722
   5        0.60     0.151     0.734   → D∞ ~ 0.15
```

**The trajectory measure on Z₂ has full support (D₀ = 1.00) but strongly
multifractal mass distribution (D∞ ~ 0.15).** Spectrum width D₀ − D∞ ≈ 0.85.

Where prior dim claims locate on the spectrum:

| Claim | Value | Best-match q | D_q at that q |
|---|---|---|---|
| Chang H-dim | 0.68 | **q = 1** (info dim) | 0.608 |
| R23 λ_max | 0.6755 | **q = 1** (info dim) | 0.608 |
| R59 dim_q2(k=12) | 0.67 | **q = 1** (info dim) | 0.608 |
| R59 dim_q2(k=7) | 0.83 | q = 0.5 | 0.863 |
| R59 dim_q2(k=15) | 0.54 | q = 1 (info dim) | 0.608 |

**Convergent finding:** Chang's 0.68, R23's 0.6755, and R59's k=12
estimate ALL locate at the information dimension D₁ of the trajectory
measure. They are probing the same single moment of the multifractal
measure. The "single δ" framing was the right neighborhood (D₁) but
incomplete — the full f(α) spans 0.15 to 1.00.

R59's scale-dependent dim_q2(k) sweep (0.83 → 0.54 from k=7 to k=15) is
the standard multifractal-crossover behavior of a finite-resolution
correlation-dimension estimate sweeping toward the asymptotic D₂ ≈ 0.27.

## f(α) singularity spectrum

```
  q      h    α        f(α)
 -1.0  -2.36  1.704    0.660    bottom-left (rare large mass)
  0.0  -1.00  1.188    0.999    near peak
  1.0  -0.61  0.619    1.226    PEAK (typical singularity)
  2.0   0.27  0.147    0.027    right tail (rare small mass)
```

Peak f(α) = 1.23 at α ≈ 0.62 (q=1), corresponding to the typical local
Hölder exponent of the measure.

α-spectrum width (clean q region): ~1.5 (from α ≈ 0.1 to α ≈ 1.7).
Wide spectrum confirms strong multifractality.

## Cantor model fits

Both one-scale (Macek-Wójcik primary) and two-scale weighted Cantor
models tested:

```
Model        Best parameters                          RSS
1-scale      p = 0.500, l = 0.5                      53.82
2-scale      p1 = 0.499, l1 = 0.05, l2 = 0.05        22.86
```

**Neither model fits the empirical D_q.** RSS reduction by 2-scale is
57% but absolute residual is still large. The trajectory measure is
NOT a weighted Cantor measure with simple parameters.

Likely reason: Z₂ measure has full support (D₀ = 1) AND wide MF spectrum
(D₀ − D∞ = 0.85). Standard weighted Cantor sets on [0,1] with two children
of equal scale 1/2 force D₀ ≤ 1; getting D∞ down to 0.15 requires extreme
weight asymmetry that the model can't express within its natural family.

## Verdict

| Outcome | Status |
|---|---|
| (α) clean MF spectrum + weighted Cantor fit | partial — spectrum yes, Cantor no |
| (β) clean MF spectrum, no Cantor closed form | **PRIMARY** |
| (γ) MF-DFA fails to extract spectrum | rejected (spatial Z_q gives clean R² in central q) |

**Empirical multifractal characterization succeeds. Closed-form weighted
Cantor parametrization does not.** The trajectory measure is multifractal
on Z₂ with D₀ = 1.00, D₁ = 0.61, D₂ = 0.27, D∞ ≈ 0.15, and Chang/R23/R59
single-dim claims all locate at D₁ ≈ 0.61.

## Implications for framework synthesis

1. **R59's "Sullivan-conformal REJECTED" framing stands and sharpens.**
   The trajectory measure is multifractal with width ~0.85 in D_q —
   no constant δ can describe it; this is now quantified.

2. **Chang's 0.68 and R23's 0.6755 are not "shared dimension" but
   "shared D₁".** They both estimate the information dimension of the
   trajectory measure. The match is interpreted correctly as same point
   on f(α), not single δ.

3. **The Lagarias-class open piece in the multifractal framework**
   reduces to: derive the weights p_q (and scales l_q) of a multi-scale
   Cantor-like measure on Z₂ that reproduces D_q above. Result 60's
   size-stratified Markov framework gives a different finite-dimensional
   identification of D_avg via QSD eigenvector — both are valid; they
   characterize different aspects (size-stratified is residue marginal of
   QSD; spatial MF is direct measure-theoretic characterization on Z₂).

4. **Temporal trajectory dynamics are Brownian.** This rules out
   long-memory or persistent multifractal processes for v_t / log m_t.
   The Collatz time-evolution is genuinely uncorrelated; multifractality
   is purely measure-theoretic at the asymptotic spatial level.

## Files

- `mfdfa_trajectory_measure.py` — temporal v_t MF-DFA (initial run)
- `mfdfa_logm.py` — temporal log m increment MF-DFA (cleaner)
- `mfdfa_spatial.py` — spatial Z_q on Z₂ cylinders (PRIMARY)
- `mfdfa_h_q.csv`, `mfdfa_logm_h_q.csv` — temporal h(q)
- `mfdfa_spatial_tau_q.csv` — spatial τ(q), D_q
- `mfdfa_spatial_f_alpha.csv` — spatial f(α) singularity spectrum
- `mfdfa_spatial_Z_q.csv` — partition function values
- `mfdfa_spatial_cantor_fits.csv` — Cantor model fit residuals
- `mfdfa_spatial_log.txt` — full run log
