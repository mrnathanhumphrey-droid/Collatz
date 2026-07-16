# Result 41 (qx+1 paper) — our r_q=1 boundary IS Siegel's σ_H=1 (α_H(0)=1) degenerate case: the spectral gap closes exactly at the dynamically MARGINAL map. (Non-singularity was the wrong object.)

**Date:** 2026-07-16. **Verdicts: ✗ H_NONSING REFUTED (my initial guess — non-singularity holds ∀q, pre-committed to lose) / ★ H_DEGEN CONFIRMED (α_H(0)=1 ⟺ σ_H=1 ⟺ q=3 ⟺ r_q=1, exact) / H_QUANT: r_q < 2^{1−σ_H} (flagged, no fit).**

**Headline: reading Siegel's dissertation (arXiv 2412.02902) named our phase boundary. Our r_q=1 / d=2 / exceptional-point boundary (q=3) = Siegel's `α_H(0)=1 ⟺ σ_H=1` DEGENERATE case, which he explicitly flags as a "degenerate measure" (eq 4.144–4.145, p277). σ_H = log_2((q+1)/2) is the growth exponent; σ_H=1 (q=3) is MARGINAL, σ_H>1 (q≥5) is expanding. So the spectral gap closes exactly at the dynamically marginal map — which is why 3x+1 is the critical case and 5x+1, 7x+1 (σ_H>1, known divergent orbits) are not. My initial "non-singularity" guess was WRONG (holds ∀q); the marginal-exponent degeneracy is the right object.**

Probe: `probe_41_siegel_alpha_boundary.py`. Log: `result_41_siegel_alpha_boundary_log.txt`. Runtime: instant (exact rationals).

## Context — the Siegel read

Siegel's 2022 USC dissertation *"(p,q)-adic Analysis and the Collatz Conjecture"* (arXiv 2412.02902, posted Dec 2024) is our anchor reference — it supplies the Parseval identity (his eq 2.180 = our `‖π_k‖²` identity = our `Parseval.lean`), the μ̂ recursion (Prop 2.18), the Fourier symbol, and the circulant eigenvalue theorem (Ch3). His Chapter 4 studies `χ_H` (the Collatz numen = our self-similar measure) via its Fourier transform. Lead 1 from that read: locate our phase boundary in his framework.

## Setup — our qx+1 map as a 2-Hydra map (p=2)

Branch 0 (even) `H_0(x)=x/2`: `H'_0(0)=1/2`, `H_0(0)=0`. Branch 1 (odd) `H_1(x)=(qx+1)/2`: `H'_1(0)=q/2`, `H_1(0)=1/2`. Siegel's Ch4 quantities:
```
    α_H(0)   = (1/p)Σ_j H'_j(0) = (1+q)/4         β_H(1/2) = -1/4
    α_H(1/2) = (1-q)/4                             γ_H(1/2) = 1/(q-1)
    σ_H      = log_p(Σ_j H'_j(0)) = log_2((q+1)/2)   (σ_H=1+log_2 α_H(0), α_H(0)=2^{σ_H-1})
```

## Results (exact rationals)

| q | α_H(0)=(1+q)/4 | α_H(1/2) | σ_H | γ_H(1/2) | r_q | 2^{1−σ_H}=1/α_H(0) |
|---|---|---|---|---|---|---|
| **3** | **1** | −1/2 | **1.0000** | 1/2 | **1.00** | 1.000 |
| 5 | 3/2 | −1 | 1.585 | 1/4 | 0.62 | 0.667 |
| 7 | 2 | −3/2 | 2.000 | 1/6 | 0.39 | 0.500 |
| 11 | 3 | −5/2 | 2.585 | 1/10 | — | 0.333 |
| 13 | 7/2 | −3 | 2.807 | 1/12 | — | 0.286 |

## H_NONSING — my initial guess, REFUTED (pre-committed to lose)

I first proposed the boundary was Siegel **non-singularity** (Def 4.2: `α_H(j/p)=0`). Computing: `α_H(0)=(1+q)/4` and `α_H(1/2)=(1−q)/4` are **never 0** for odd q. So non-singularity holds for **all** q and does not single out q=3. **My "non-singularity = boundary" framing was wrong** — reported as a loss, per the pre-registration.

## ★ H_DEGEN — the corrected identification, CONFIRMED

`{q : α_H(0)=1} = {3} = {q : r_q=1}`. Exactly:
```
    α_H(0) = (1+q)/4 = 1  ⟺  q=3  ⟺  σ_H = log_2((q+1)/2) = 1  ⟺  r_3 = 1 (gap closed).
```
And Siegel **explicitly** flags `α_H(0)=1` at p=2 as a **"degenerate measure"** (eq 4.144–4.145): "χ̂'_H(t) and χ̂_H(t) differ by a factor of γ_H(1/2)Â_H(t), which is a degenerate measure because α_H(0)=1." **That degenerate case is exactly our critical point / exceptional point (R39).** For q≥5, `α_H(0)>1` (`σ_H>1`) and the gap opens (`r_q<1`).

## The dynamical reading — the gap closes at the MARGINAL map

`σ_H = log_2((q+1)/2)` is the map's **growth exponent** (`σ_H = log_p Σ_j H'_j(0)`, the log-average branch multiplier). The three regimes:
- **q=3 (3x+1): σ_H = 1 — MARGINAL** (average multiplier `(q+1)/2 = 2 = p`, growth balances halving). Spectral gap **closed** (`r_3=1`). This is the conjectural/critical Collatz case.
- **q≥5 (5x+1, 7x+1, …): σ_H > 1 — EXPANDING.** These maps have known divergent orbits. Spectral gap **open** (`r_q<1`).

So the spectral-gap closure is not arithmetic happenstance — it coincides with the **dynamical marginality** of 3x+1. The gap opens precisely for the expanding (σ_H>1) siblings. This gives the boundary (R6's `(q−1)/d=1`, R39's exceptional point) a **dynamical meaning** and a **citable home** in Siegel.

## H_QUANT — bonus observation (flagged, NOT claimed)

`r_q` vs `2^{1−σ_H} = 1/α_H(0)`: q=5 → r=0.62 vs 0.667; q=7 → r=0.39 vs 0.500. So `r_q < 2^{1−σ_H}` at both — the true gap is **strictly stronger** than the marginal-exponent rate. This could be a real upper bound `r_q < 1/α_H(0)` (gap at least the marginal rate), but only 2 points and R28 killed elementary closed forms for r_q — so it is an observation to chase, not a result.

## What this does / does not do for L3

- **Does:** names the boundary. The R6/R39 phase boundary is Siegel's `σ_H=1` (`α_H(0)=1`) marginal-exponent degeneracy — a dynamical criterion, citable, and it explains *why* q=3 is the critical q (marginal growth), consistent with R40's finding that q=3 is uniquely critical.
- **Does NOT:** close L3. The **decay RATE / gap SIZE** (`r_q`, the q-adic correlation decay) is still the open piece — Siegel supplies existence/quasi-integrability of χ̂_H and the symbol, but the q-adic decay rate is exactly what he leaves open (the ‖π_k‖² spectral bound). So the division of labor holds: Siegel = foundation + boundary name; L3 = the rate, ours.

## Not at stake
R1–R40. This identifies the boundary in Siegel's framework; it does not change any r_q value or the d=2 characterization.

_Reporting discipline: my initial "non-singularity" identification was pre-registered AS my guess and AS predicted-to-lose (α_H never vanishes for odd q) — reported as a loss, not quietly swapped. The corrected H_DEGEN is an exact set-equality ({q:α_H(0)=1}={3}), not a fit. The r_q < 2^{1−σ_H} pattern is disclosed as a 2-point observation, not claimed (R28 discipline). All Siegel quantities computed as exact rationals from his Def 4.1/4.4._
