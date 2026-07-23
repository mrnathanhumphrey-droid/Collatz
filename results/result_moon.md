# Probe MOON — phase-diversity joint fit — **diversity is PRESENT but structured (bore lives in the sub-dominant oscillatory channels); the joint fit RAILS to the boundary — phase diversity does not rescue the period. Integral channel closed.**

**Date:** 2026-07-23. Probe `probes/probe_moon.py`. The m-channels `A_r(m) = C_{r+1}(m)/3 = γ_r(τ_m) − γ_{r−1}(τ_m)`
(so `Λ_r = Σ_m 4^{−m} A_r(m)`) are different observables of the same bore operator, sharing (ρ,θ) with per-channel
(A_m, φ_m). If the channel phases differ, phase diversity substitutes for the r-coverage we can't get. Exact A_r(m)
r≤7 (character ledger, cross-checked), float r=8..12 (build_nu). **Model-adequacy diagnostic, NOT a period
measurement (R27-A: L(z) not rational ⟹ Λ_r is not a finite sum of exponentials).**

## M-A — the diversity check: **present, and structured**
Sign sequences of `A_r(m)`, r=2..12:

| m | sign sequence | turnovers |
|---|---|---|
| 1 | `+ + + + + + + + + + +` | **none (monotone +)** |
| 2 | `− + − + − − − − − + −` | r=3,4,5,6,11,12 (oscillatory) |
| 3 | `− − − − − − − − − − −` | **none (monotone −)** |
| 4 | `+ + − − + − − − − + +` | r=4,6,7,11 (oscillatory) |

The channels **split**: `m=1` (monotone +), `m=3` (monotone −) carry **no oscillation**; `m=2, m=4` carry the
oscillation, turning near r=6 (consistent with R13-C's "common turnover ~r=6" — but only for the oscillatory
channels). **Structural finding:** the **dominant** channel `m=1` (weight `4⁻¹` in Λ_r) is monotone; the **bore
oscillation lives in the sub-dominant channels `m=2,4`** (weights `4⁻²,4⁻⁴`). So `Λ_r = (monotone dominant) + (small
oscillatory)` — the bore is a *weak* oscillation riding a *strong* monotone decay. That is precisely why it is so
hard to measure, and why every fit in this campaign has been pulled toward long/flat periods.

Diversity is present (turnovers spread r=3..12), so the kill condition is not met — M-B runs.

## M-B / M-C — the joint fit **RAILS to the boundary** (free and ρ-fixed, all subsets identically)
Joint fit `A_r(m) = A_m ρ^r cos(rθ + φ_m)`, shared (ρ,θ), r=6..12, admitted channels {1,2,3,4,9,27} (M-D: all
exceed transient contamination at r≥10):

- **Free (ρ,θ): rails to ρ = 0.940 (floor), P = 30.0 (ceiling)** — no interior optimum. φ-spread is large (227°) but
  that reflects the channels' disparate behavior, not a resolved shared mode.
- **ρ fixed at 0.984 (the known rate): P rails to 30.0** for all six channels **and** for both disjoint subsets:

| subset | best P (ρ=0.984) |
|---|---|
| all six | 30.0 (ceiling) |
| {1,3} (monotone) | 30.0 (ceiling) |
| {2,4} (oscillatory) | 30.0 (ceiling) |

Both subsets rail to the **same boundary** — this is **not** genuine agreement (adequacy) and **not** disagreement
(continuum finding); it is **identical degenerate failure**: the two-mode model wants an ever-longer period because
the r=6..12 window can't close a cycle, and the weak oscillatory channels can't pull P down against the monotone
ones. Extending the grid would only rail further. **Phase diversity is present but insufficient — it does not
substitute for the missing time coverage.**

## Scope / verdict
**MOON: diversity present, but the joint fit is under-determined — the integral channel is closed.** **M-A** —
diversity is real and structured: the bore oscillation is carried by the sub-dominant channels `m=2,4`, while the
dominant `m=1` and `m=3` are monotone. **M-B/M-C** — the joint fit rails to the boundary (long P, low ρ) whether free
or ρ-fixed, all subsets identically; no interior optimum, so the two-mode model does not resolve a shared (ρ,θ) on
r=6..12. This is **doubly consistent**: with the resolution wall (Probe D, SAT — the window is <1 period), and with
**R27-A's proven continuum** (L(z) not rational ⟹ no single discrete mode for the fit to converge to). The rail-to-
long-period agrees in direction with SAT-D (~17.6): the period is long and unresolvably so from r≤12.

**Consequence for the crux (owed to the pen).** Phase diversity was the last untried trick in the integral channel,
and it is spent: the channels ARE phase-diverse, but the oscillatory signal (in the weight-suppressed sub-dominant
channels) is too weak, over a <1-period window, to jointly pin (ρ,θ) — every fit rails long. The one genuinely new
result is structural, not spectral: **the bore is a weak oscillation on the sub-dominant m-channels, not on the
dominant one** — which both explains the measurement difficulty and localizes the mode (echoing F1-C: the mode lives
in the intermediate/sub-dominant structure, not the dominant part). Every integral-channel probe — Probe D (short
baseline), SAT (infinite baseline), MOON (phase diversity) — now returns the same verdict from a different angle:
**a long, non-integer, unpinnable period, requiring exact `S_k` to k≈17–20 to measure.** No fitting on r<6; the
boundary railing reported as such (not smoothed into a period); the scope (model-adequacy, not measurement) honored.
Not at stake: R1–R30, R80–R82, Probe D, SAT.
