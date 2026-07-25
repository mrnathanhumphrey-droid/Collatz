# Probe MODES — the n=1/n=2/n=3 decomposition of the coupling — **Re δ̂(1) (the dominant fluctuation mode) is single-signed POSITIVE for all r=2..12 and does NOT flip; the g_r crossings (r=2, r=6) are n=2 / n=3 EXCURSIONS overpowering a small-but-positive n=1 — a single-signed dominant mode with a bounded, spiking (lengthening-period) perturbation. This corrects `result_fiber.md`'s "no crossing / 0.477" and reduces the sign question to a ratio bound `|Re δ̂(2)| < 4|Re δ̂(1)|` (+ the n=3 bound).**

**Date:** 2026-07-25. Probe `probes/probe_modes.py`. `g_r = ⟨δ_r,Re w⟩ = Σ_{n≥1} 4^{−n} Re δ̂_r(n)`,
`δ̂_r(n)=Σ_k δ_r(k)e^{−2πink/3^r}` = lag-n autocorrelation of ν in the orbit coordinate. Re w's weights fall ×¼ per
mode; `3|n` = fiber-mean, `3∤n` = fiber-fluctuation. So n=1 carries ~78% of the fluctuation weight, n=2 ~19%.

## Check 1 — closure (PASS)
`¼Re δ̂(1) + 1/16 Re δ̂(2)` reproduces the fluctuation coupling `Σ_{3∤n}4^{−n}Re δ̂(n)` to **~0.3% at large r**
(r=12: +7.293e−4 vs +7.316e−4). And `g_r = fluctuation(3∤n) + fiber-mean(3|n)` exactly.

## Check 2 — THE crossings: n=1 never flips
| r | Re δ̂(1) | Re δ̂(2) | Re δ̂(3) | \|d2\|/\|d1\| | g_r |
|---|---|---|---|---|---|
| **2** | **+5.71e−2** | **−3.71e−1** | +0 | **6.50** | **−** |
| 3 | +1.98e−2 | +1.15e−2 | −2.04e−1 | 0.58 | + |
| 4 | +8.64e−3 | −3.42e−3 | −2.13e−2 | 0.40 | + |
| 5 | +8.04e−3 | +1.88e−3 | −9.48e−2 | 0.23 | + |
| **6** | **+7.74e−3** | −2.55e−2 | **−7.50e−2** | 3.29 | **−** |
| 7 | +5.79e−3 | −1.39e−2 | −7.02e−3 | 2.40 | + |
| 10 | +3.60e−3 | −5.00e−4 | −1.46e−3 | 0.14 | + |
| 12 | +2.96e−3 | −1.85e−4 | −8.43e−4 | 0.06 | + |

**`Re δ̂(1) > 0` for every r=2..12 — the dominant mode is single-signed positive.** The two g_r sign changes are
*not* n=1 flips:
- **r=2:** an **n=2 excursion** — `Re δ̂(2)=−0.371` (|d2|/|d1|=6.5) overpowers a small positive n=1.
- **r=6:** the **fluctuation stays positive** (+4.4e−4); the **n=3 fiber-mean** (`Re δ̂(3)=−0.075`) tips g₆ negative.

Both crossings are the dominant positive mode being *overpowered*, not *reversed*. This is the stronger structural
position: a single-signed dominant mode with a bounded perturbation.

## Check 3 — rates and the perturbation envelope
- **`Re δ̂(1)` decays at ~0.9** (late window r≥8: 0.924, 0.892, 0.923) — the ladder rate.
- **`|Re δ̂(2)|/|Re δ̂(1)|` spikes exactly at the crossings** (6.50 at r=2, 3.29 at r=6) and is small between
  (0.06–0.72 at large r). The spike envelope is **decreasing** (6.5 → 3.3), with the crossings at gaps 2, ~4 — the
  **lengthening-period** signature. The n=3 fiber-mean plays the same overpowering role at r=6.

## What this means for 7/15 vs 0.477 — and what it does NOT settle
The sign of g_r is `sign(dominant positive n=1 + n=2/n=3 perturbation)`. The dominant mode is **positive throughout**
and decays at 0.9; the perturbation **spikes periodically** (lengthening) and its envelope is **decreasing**. So:
- If the spike envelope continues below the crossing threshold, g_r eventually stays positive ⟹ **S_∞≈0.477**.
- If the (lengthening) spikes keep crossing the threshold, g_r keeps crossing ⟹ 7/15 stays live via a very-late
  turnover. The next spike is due past r≈16, unobserved.
- **Two spikes = two points (#30):** the decreasing envelope (6.5→3.3) is real but cannot be extrapolated to
  "eventually below threshold" from two data. So this **leans 0.477** (dominant mode positive, envelope decreasing)
  but does **not** establish it — the lengthening-crossing mechanism is alive, now identified as the
  n=1-vs-(n=2,n=3) competition.

## The reduced question (pen)
`Re δ̂_r(1)` is a **coefficient sign** — the lag-1 autocorrelation of ν in the orbit coordinate — and Bochner /
positivity **cannot** deliver it (a measure near x=½ has `Re δ̂(1)<0`); it needs an analytic lower bound from the
construction of ν_r. But given `Re δ̂(1)>0`, the crossings are governed by the **ratio** `|Re δ̂(2)| < 4|Re δ̂(1)|`
(and the n=3 analogue) — a bound positivity arguments *can* address. The two sub-questions:
1. **Is `Re δ̂(1)>0` provable** (single-signed dominant mode)? Analytic estimate on ν_r, not a positivity/pen step.
2. **Does the perturbation envelope `|Re δ̂(2)|/|Re δ̂(1)|` stay below 4 for large r** (bounded, decreasing)? A ratio
   bound — the first thing in this arc positivity can attack.

## Status
**MODES:** `g_r = Σ_{n≥1}4^{−n}Re δ̂_r(n)`; closure ✓ (¼d1+1/16 d2 = fluctuation to 0.3% late). **`Re δ̂(1)>0` for
all r=2..12 — single-signed dominant mode, decays ~0.9.** The g_r crossings (r=2, r=6) are **n=2 / n=3 excursions**
overpowering a small positive n=1 (|d2|/|d1| spikes 6.5, 3.3 at the crossings; lengthening gaps 2, ~4) — **NOT n=1
flips.** Corrects `result_fiber.md` ("no crossing / 0.477 single channel" retracted: g_r DOES cross, fiber-mean
dominates at r=6). Reduces 7/15-vs-0.477 to: (1) is `Re δ̂(1)>0` provable (analytic, needs ν_r construction — Bochner
insufficient), and (2) does `|Re δ̂(2)|<4|Re δ̂(1)|` hold for large r (ratio bound, positivity-addressable). Leans
0.477 (dominant mode positive + decreasing spike envelope) but does NOT establish it (2 spikes = #30; lengthening
mechanism alive). Not at stake: R1–R30, R80–R82, g_r ladder to r=12, Λ^unif closed form. The fiber ladder and S2
splitting are ONE method extended (S2 continues g_r past r=16), not two.
