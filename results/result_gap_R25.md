# Probe R25 — the spectral gap (subcritical) — **the gatekeeper measured: healthy away from criticality, narrowing toward it (unresolved)**

**Date:** 2026-07-22  Reuses the renewal builder, now deep (r→14 via truncated v-loop). Probe
`probes/probe_gap_R25.py`. Wilson's redirect: drop X_∞; use **C(λ):=lim_r S_r/ρ^r** (per-r, no tail), where
ρ=3(1−λ)/(1+λ)=q·Σ_v p_v² is the exact **leading eigenvalue** of the transfer operator. Since ρ(½)=1, **C(½)=S_∞**,
so the theorem is the **boundary value C(½)=7/15** of a function analytic on (½,1) — *provided the spectral gap
survives*, i.e. **|λ₂|/ρ stays bounded < 1 as ε→0.** That is the gatekeeper, and R25-C is the first look at it.

## R25-A — DEFINITION AUDIT: the R24 inconsistency resolved
`1−ρ = 8ε/(3+2ε)` holds exactly (verified all ε). The R24 discrepancy (amplitude/(ε·X_∞) off by 71%→6%) was
**definitional**: the amplitude column used Y_∞ = 1+Σ_r S_r while the ε·X_∞ column used X_∞ = Σ_r S_r — the `X₀=1`
term. With the consistent amplitude `(1−ρ)·X_∞`, the ratio is `(1−ρ)/ε = 8/(3+2ε)` exactly. (And note
`(1−ρ)X_∞ = Cρ + subdominant corrections`, not clean C — which is *why* we drop X_∞ and use the plateau.)

## R25-B — THE PLATEAU C(λ)=lim S_r/ρ^r: the deep values point at 7/15, but the plateau isn't converged
`p_r = S_r/ρ^r` at the deepest reachable depth:

| ε | ρ | p₁₄ | drift near r=14 |
|---|---|---|---|
| 0.10 | 0.750 | 0.2751 | decreasing to C |
| 0.05 | 0.871 | 0.3741 | decreasing to C |
| 0.02 | 0.947 | 0.4328 | ~flat / turning |
| 0.01 | 0.974 | 0.4518 | **increasing, not converged** |
| 0.005 | 0.987 | 0.4611 | **increasing, not converged** |
| 0.002 | 0.995 | **0.4666** | **increasing, not converged** |

**`p₁₄(ε)` climbs monotonically to 7/15 = 0.46667 as ε→0** (0.275 → 0.4666), consistent with C(½)=7/15. But the
per-ε plateau is *not* converged at r=14 — `p_r` still drifts, and the direction flips (decreasing for large ε where
the gap is healthy, increasing for small ε where it is not). So this supports 7/15 but does not lock C(λ) at any
single ε.

## R25-C — THE GATEKEEPER |λ₂|/ρ: healthy away from criticality, **narrowing toward it**
From `d_r = p_{r+1}−p_r ~ C₂(λ₂/ρ)^r(λ₂/ρ−1)`:

| ε | ρ | `d_r` behavior (r→14) | `|λ₂|/ρ` |
|---|---|---|---|
| 0.10 | 0.750 | shrinking cleanly, real | **0.69** |
| 0.05 | 0.871 | shrinking cleanly, real | **0.57** |
| 0.02 | 0.947 | oscillatory transient, marginal | ~0.8 (noisy) |
| 0.01 | 0.974 | **`d_r` ≈ constant — not yet shrinking** | unmeasurable (>1 = artifact) |
| 0.005 | 0.987 | not shrinking | unmeasurable |
| 0.002 | 0.995 | not shrinking | unmeasurable |

- **Large ε: the gap is healthy** — `|λ₂|/ρ ≈ 0.57–0.69`, real subdominant, plateau forms.
- **Small ε (≤0.02): the plateau does not form by r=14** — `d_r` stops shrinking (nearly constant at ε=0.01), so the
  `|λ₂|/ρ > 1` readings are **artifacts of non-asymptotic depth, not real values.** Their *meaning* is that
  `|λ₂|/ρ` is close to 1 — the gap is narrowing as ε→0 — with an oscillatory (complex-pair) transient emerging
  around ε≈0.02.

**This is the live risk Wilson named, now observed.** At criticality the subdominant is the period-9 complex pair on
the unit circle (the ≈0.98-envelope oscillation); subcritically it sits inside, and the data **leans toward
`|λ₂|/ρ → 1` — the gap shutting exactly at the point of interest.** Not proven: r=14 is too shallow at small ε to
resolve the complex pair's rate. But the qualitative trend (fast-shrinking `d_r` at large ε, non-shrinking at small
ε) is the concerning one, not the reassuring one.

## R25-D — Richardson: unreliable here
With the per-ε plateau unconverged (esp. small ε), Richardson/Aitken on C(λ) is not trustworthy (the C estimates
scatter 0.27–0.52). The clean statement is the direct one: **p₁₄(ε) → 7/15 as ε→0.**

## Status
**R25: the gatekeeper is measured — and it is the crux, not a formality.** **A** resolves the R24 inconsistency
(Y_∞ vs X_∞). **B** the deep values `p₁₄(ε) → 7/15` support the boundary value C(½)=7/15, but the per-ε plateau is
not converged. **C (the theorem's actual gatekeeper, looked at for the first time):** the spectral gap is **healthy
away from criticality** (`|λ₂|/ρ ≈ 0.6` at ε=0.05–0.1, real) but **narrows toward criticality** — by ε≤0.01 the
plateau hasn't formed at r=14, an oscillatory complex pair emerges (ε≈0.02), and the trend leans toward
**`|λ₂|/ρ → 1`, the gap shutting at λ=½** (the period-9 mechanism). **D** Richardson unreliable.

**Consequence for the crux (owed to the pen):** the whole subcritical route's viability reduces to one measurable
number — **does the subdominant `|λ₂|/ρ` stay bounded below 1 as ε→0, or tend to 1?** R25 shows it is healthy at
ε≥0.05 and *appears to narrow* as ε→0, but **r=14 is too shallow to resolve the small-ε limit** (build wall:
support 2·3^{r−1}, r≈16 max even with the truncated-v build). Two ways to resolve it, both cleaner than S_r
asymptotics: (i) push the depth (r→16) and re-measure the complex-pair rate; (ii) **compute the transfer operator's
second eigenvalue directly** (a finite-matrix eigenproblem at fixed level, no r→∞ extrapolation) — the definitive
gatekeeper measurement. If `|λ₂|/ρ` is bounded < 1, the route closes and C(½)=7/15 follows by eigenvalue
perturbation; if it → 1, the subcritical approach fails for a *nameable* reason (the gap shuts, and the same complex
pair *is* the critical period-9 oscillation). No fitting; exact ρ, labeled numeric plateau/gap; the gatekeeper's
small-ε limit reported as unresolved at reachable depth.
