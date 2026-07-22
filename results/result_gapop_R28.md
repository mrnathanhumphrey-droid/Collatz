# Probe R28 — the gap-correlation operator — **A GATE PASS: R27-E operator certified (odd-kill exact); B/C bound is mechanism-level; D failed**

**Date:** 2026-07-22  Certified renewal `X'=1+3·2⁻ᵛX` on `1+3ℤ/3^{r+1}`. Probe `probes/probe_gapop_R28.py`. Gates
Wilson's R27-E derived operator: `ν̂_r(ξ)=e(ξ/3^{r+1})·E_v[ν̂_{r−1}(ξ2⁻ᵛ)]`, gap correlation
`R_r(d)=Σ_η ν̂_r(η)·conj(ν̂_r(η·2^d))`.

## R28-A — GATE: **PASS** — the operator closes exactly, odd-d killed exactly
The closing identity holds in the normalization `R_{r−1}(0)=3·X_{r−1}` (so `X_r = Σ_d P(d)R_{r−1}(d)` — equivalently
Wilson's `S_r = 3Σ_d P(d)R_{r−1}(d)` with his R = mine/3), `P(d)=(1−λ)λ^{|d|}/(1+λ)`, at λ=½:

| r | Σ_d P(d)R_{r−1}(d) | X_r (bank) | ratio | d=0 diag P(0)R(0) | X_{r−1} | increment | S_r | max\|odd-d R\| |
|---|---|---|---|---|---|---|---|---|
| 2 | 2.142857 | 2.142857 | 1.000000 | 1.666667 | 1.666667 | — | 0.476190 | 3e−16 |
| 4 | 3.068646 | 3.068646 | 1.000000 | 2.604432 | 2.604432 | 0.464214 | 0.464214 | 7e−16 |
| 7 | 4.465821 | 4.465821 | 1.000000 | 4.000330 | 4.000330 | 0.465491 | 0.465491 | 9e−16 |

**Three derived facts confirmed exactly:**
- **The operator closes:** `X_r = Σ_d P(d)R_{r−1}(d)` to ratio 1.000000, r=2…7; the increments `X_r−X_{r−1}` **are**
  `S_r` exactly.
- **Diagonal Flatness:** the d=0 term `P(0)R_{r−1}(0) = X_{r−1}` exactly (`3P(0)=1` at λ=½) — the diagonal channel
  is the previous level, recovered as a single channel.
- **Odd-gap kill, DERIVED:** all odd-d correlations vanish to **machine zero** (`max|odd-d R| ~ 1e−15`). The
  conjugate-kill Wilson derived (`ord₃(2)=2` ⟹ `3∤2^d−1` for odd d ⟹ maximal-order twist ⟹ total cancellation)
  now has its source, giving R7's odd-gap vanishing a mechanism.

**R27-E is certified. #43 not incurred.** The second-moment recursion has a clean derived operator form, and the even
channels `R(d=2m)` carry stratum index `v₃(m)` (R8's W_j classes), also as Wilson derived.

## R28-B — κ measured: the bound *direction* holds, but R(2)/R(0) is not the eigenvalue ratio
The m=±1 twisted channel over the untwisted, `|R(2)/R(0)|` per r:

| λ | r=3 | r=5 | r=7 | 2λ² |
|---|---|---|---|---|
| 0.50 | 0.3238 | 0.2304 | 0.1786 | 0.500 |
| 0.55 | 0.3435 | 0.2639 | 0.2273 | 0.605 |
| 0.60 | 0.3656 | 0.3065 | 0.2861 | 0.720 |

**The twisted channel is strictly smaller than the untwisted at every r** (`|R(2)/R(0)| < 2λ²` throughout) —
confirming Wilson's `κ < 1` *direction* (nontrivial character sum < untwisted ⟹ `|λ₂|/ρ < 2λ²`). But `|R(2)/R(0)|`
**declines with r** (0.37→0.18 at λ=0.6) and does **not** equal the measured eigenvalue ratio (`|λ₂|/ρ ≈ 0.6–0.68`
from R26/R27) — so the raw channel ratio is *not* κ as a clean eigenvalue observable; the eigenvalue link runs
through the operator's full spectral action, not a single channel.

## R28-C — channel mixing: LARGE ⟹ the bound is mechanism-level, not rigorous
| r | \|R(4)/R(2)\| (m=±2) | \|R(6)/R(2)\| (m=±3) |
|---|---|---|
| 3 | 0.706 | 2.059 |
| 5 | 0.698 | 1.874 |
| 7 | 0.675 | 1.744 |

**The higher channels are comparable to or LARGER than m=±1** (`|R(6)/R(2)| ≈ 1.7–2.1` — the m=±3 channel, a
*coarser stratum* `v₃(3)=1`, exceeds m=±1). So isolating the m=±1 channel is **not free** — the convolution couples
gaps strongly, exactly Wilson's stated caveat. The bound `|λ₂|/ρ < 2λ²` is therefore **mechanism-level, not
rigorous**: closing it needs control of channel mixing, and the mixing is measured to be O(1), not small.

## R28-D — real+pair fit: **FAILED** (deflation ill-conditioned)
Deflating by the fixed real mode 0.5 (`μ_r = Λ_{r+1} − 0.5·Λ_r`) gives `μ = +0.040, +0.005, −1e−5, +2e−6, −5e−4,
+4e−4` — the middle values `μ_3, μ_4 ≈ 1e−6` are near-zero (r=3,4 is almost exactly 0.5-geometric, so deflation
annihilates it), leaving an **ill-conditioned** sequence dominated by the transient (`μ_1, μ_2`) and a late
perturbation (`μ_5, μ_6`). The 2-term solves on μ give garbage (roots −283, +6.9). **Fixing the real mode did NOT
lift the ≥3-mode underdetermination** — contrary to the expectation that D would pay. The period remains unresolved.

## Status
**R28: the operator is certified; the gap bound is mechanism-level; the period is still open.** **A GATE PASS** —
Wilson's R27-E gap-correlation operator closes exactly (`X_r = Σ_d P(d)R_{r−1}(d)`, r=2…7), with **Diagonal Flatness
exact** and the **odd-gap conjugate-kill derived and confirmed to machine zero** (R7's kill now has a source);
even-d channels carry the `v₃(m)` strata. **#43 not incurred.** **B** — the twisted channel is `< 2λ²` at every r
(bound *direction* confirmed) but declines and isn't the raw eigenvalue ratio. **C** — channel mixing is **O(1)**
(m=±3 channel *larger* than m=±1), so the bound `|λ₂|/ρ < 2λ²` is **mechanism-level, not rigorous** — the
channel-mixing gap is real and measured. **D** — deflation-by-0.5 failed (ill-conditioned); period unresolved.

**Consequence for the crux (owed to the pen):** the gap-survival conclusion (`|λ₂| ≈ 0.5 < 1`, from R26/R27/R18-A)
now has a **certified operator foundation** — the second moment closes through Wilson's gap operator, with the
odd-kill and strata derived rather than observed. What R28 does *not* deliver: a rigorous bound (channel mixing is
O(1), so the clean `|λ₂|/ρ < 2λ²` argument needs mixing control) or the period (D failed). So the honest state is:
**the operator is real and certified; gap-survival is strongly supported (three numerical lines + the operator's
structure) but not yet proven; the remaining analytic work is control of the gap-channel mixing** — the one place
the bound leaks. No fitting; exact operator-closure and machine-zero odd-kill as gates, labeled channel ratios; the
non-rigor of the bound and the failure of D reported plainly.
