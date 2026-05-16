# Route B' — Direct fit to ε_k data underconstrains slow-mode parameters

**Date:** 2026-05-16 (after Route B R1+R2 negative on magnitude residual). **Verdict: the empirical (ρ=0.984, period=9.2) characterization is an EXTRAPOLATION from PADE Hadamard at n=13 — direct fitting to k=2..13 ε_k data does NOT robustly recover these parameters. To validate the asymptotic, ε_14, ε_15, ... computation is required (~3-10hr each per the project's compute estimate).**

## Setup

Per Route B disposition: the magnitude 0.984 lives in the ε_k bilinear pair-form generating function (NOT in the state-Markov-chain spectrum). Direct route to validate: fit a CC-pair ansatz to ε_k data and see if ρ = 0.984, θ = 2π/9.2 = 0.683 rad fall out.

Two probes:
1. `phase_routeB_prime_eps_fit.py`: fit `ε_k = A·ρ^k cos(θ k + φ)` (and 2-mode variant) to k=7..13, 5..13, 4..13, 3..13.
2. `phase_routeB_prime_eps_fit_v2.py`: fit `ε_k = A·(1/r_trans)^k + B·(1/r_slow)^k cos(θ k + φ)` (transient + slow-mode) to k=2..13, with differential evolution + 20 random seeds.

ε_k data: 13 points from `PADE_NUMERICAL_DATA.md` (exact rationals at k=1..6, numerical at k=7..13).

## Result

### Single-CC fit on k=7..13 (post-transient, 7 points)

| Model | ρ | θ | Period | Loss |
|---|---|---|---|---|
| CC | 1.075 | 0.299 | 21.0 | 4.7e-9 |
| 2-mode | 0.988 | 0.169 | 37.2 | 2.1e-9 |

Neither fit matches (0.984, 9.2). The 2-mode fit hits ρ=0.988 (close to 0.984) but with period 37 (far from 9.2). With only 7 data points and ~4 parameters, the fit is under-constrained.

### Transient + slow-mode fit on k=2..13 (12 points)

| k range | r_trans | r_slow | θ | Period | Loss |
|---|---|---|---|---|---|
| 2..13 | 1.56 | **2.50 (upper bound)** | 1.88 | 3.3 | 1.9e-5 |
| 3..13 | 4.00 (upper bound) | **1.00 (lower bound)** | 0.05 (lower bound) | 125.7 | 2.2e-6 |
| 4..13 | 1.50 (lower bound) | **1.00 (lower bound)** | 0.05 (lower bound) | 125.7 | 1.7e-6 |

The optimizer consistently hits boundary values for `r_slow` and `θ`. **The data does not robustly constrain the slow-mode parameters to (0.984, 9.2)** — the fit prefers either ρ=1.0 (no decay, with very slow oscillation) or ρ ≈ 1.2 (slow growth) or radius=2.5 (different structure entirely).

## Why this happens

The empirical PADE Hadamard at n=13 is 1.57 (= |ε_13|^{1/13}), and the asymptote z=1.016 was *predicted* in `STATE.md` based on different analysis. **At k≤13 the data is still in the transient regime** — the slow-mode at z=1.016 has not yet dominated. Per `PADE_NUMERICAL_DISPOSITION.md`:

> "The slow-mode singularity at z ≈ 1.016 (predicted from STATE.md's ρ ≈ 0.984 in k-space) is **NOT YET supported by the data.** The Hadamard radius at n=13 is 1.57, not 1.016. The slow-mode is the TRUE asymptotic but n=13 is still in the transient regime."

Direct fitting to k=2..13 cannot distinguish (0.984, 9.2) from a variety of other (ρ, θ) combinations that fit equally well.

## What this means for Route B'

The slow-mode CC pair (0.984, period 9.2) is a **PREDICTION** from STATE.md's analysis (slow-mode singularity at z=1.016 plus sign-pattern evidence of period-9 oscillation). It is **NOT directly measurable** from the k=2..13 ε_k data.

To validate Route B' (= identify the magnitude 0.984 as a specific operator eigenvalue):
- **Option (a):** compute ε_14, ε_15, ε_16, ... (each ~3-10 hrs compute) to extend Padé/Hadamard convergence beyond the transient. With ε_20 or so, the slow-mode should dominate enough to be fit robustly.
- **Option (b):** derive the asymptotic ρ analytically from the ε_k recursion structure (R76/R77 framework). This is conceptual work without compute, but requires deep familiarity with the recursion.
- **Option (c):** accept the structural reading "(0.984, 9.2) is the slow-mode asymptotic predicted by upstream analysis; direct ε_k confirmation requires k > 20" and document this as the boundary.

## Decision

Route B' is **data-limited** at k=13. The numerical PADE / Hadamard / direct-fitting work cannot validate (0.984, 9.2) without additional ε_k data. Pursuing Option (a) is a compute project (~30-100 hrs total for 5-10 more coefficients).

This bounds Route B's reach within current data. The structural identification (period-9 mechanism via cyclic-Z_18 symmetry of σ_{-1} on (Z/27)*) stands — but the asymptotic magnitude closure (validation of 0.984) is gated by compute.

## Files

- `phase_routeB_prime_eps_fit.py` + `experiments_output/phase_routeB_prime_eps_fit.json`
- `phase_routeB_prime_eps_fit_v2.py` + `experiments_output/phase_routeB_prime_eps_fit_v2.json`
- `ROUTEB_PRIME_EPS_FIT_RESULT.md` (this writeup)

## Net Route B disposition (after R1 + R2 + B' direct fit)

**Period-9 mechanism: STRUCTURALLY IDENTIFIED.** Three independent views (Phase 4 channel, (c, m) Markov chain at M=18, full state chain at n=5, 6) all give period ≈ 9.5 — within 8% of empirical 9.2. The cyclic-Z_{2·3^{n-1}} symmetry of σ_{-1} on (Z/3^n)* is the structural source.

**Magnitude 0.984: DATA-LIMITED.** Predicted asymptote from STATE.md upstream analysis. Direct fitting to k=2..13 ε_k underconstrains. Requires ε_14+ computation (~3-10hr each) OR analytic derivation from ε_k recursion (R76/R77 deep dive) OR acceptance of the predicted asymptote.

Route B is **structurally complete on the period-9 side, data-bounded on the magnitude side**. Net partial closure stands.
