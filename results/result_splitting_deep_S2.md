# Gate SPLITTING (S2) — the DEEP curve to r=40 — **a STATISTICAL ESTIMATE, not a proof (the exact ladder still stops at r=16): the log₃ ×3-crossing at [27,31] is REFUTED, and for the first time the deepest data does not point at 7/15 — ε̂ sits at the 0.477 no-turnover level, with only a sub-3σ late-rollover hint keeping the constant alive**

**Date:** 2026-07-25. Probe `probes/gate_splitting_S2_deep.py`, log `logs/S2_deep_r40_run.log`.
**This document reports a splitting ESTIMATE with an error budget — it is NOT welded to the exact ε ladder.**
Bigint accumulator (past the int64 wall r≈18) + ESS/level monitor. Stage A certified the bigint path
**bit-identical** to the S1-validated int64 estimator at r≤16 (`max|Δγ|=0.00e+00`). Stage B: N=4×10⁵, 20 replicas,
channels m∈{1,2}, r→40; **a curve, not an endpoint** — the shape is the observable (Wilson's revised spec).

## The ε̂ curve (telescoped: ε̂_r = ε_16 + 2Σ_m 4⁻ᵐ(γ_r(m)−γ_16(m)), anchor ε_16=4.685×10⁻³ exact/ν-validated)

| r | 20 | 25 | 28 | 31 | 35 | 40 |
|---|---|---|---|---|---|---|
| ε̂ (×10⁻³) | 6.37 | 7.86 | 8.41 | 8.99 | 9.18 | 8.37 |
| SE (×10⁻³) | 0.42 | 0.61 | 0.77 | 0.84 | 1.03 | 1.10 |
| z vs anchor | 4.00 | 5.24 | 4.85 | 5.14 | 4.37 | 3.36 |

Peak of ε̂ at **r≈36** (9.20×10⁻³). Drop peak→40 = +8.34×10⁻⁴ ± 3.66×10⁻⁴ (**+2.3σ**). Rise anchor→40 =
+3.68×10⁻³ ± 1.10×10⁻³ (**+3.4σ**). **Decision rule (pre-committed): INCONCLUSIVE** — neither the drop-from-peak
(2.3σ) nor the still-rising-at-40 branch cleanly ≥3σ. But the *shape* is decisive where it counts (below).

## Result 1 — the log₃ ×3-crossing at [27,31] is REFUTED (solid, tail-independent)
ε̂ climbs **monotonically and significantly (z≈5)** straight through the predicted turnover window — at r=28 and
r=31 it is still rising, near its maximum. **There is no peak in [27,31].** The log-periodic ×3 prediction
(crossing #3 at 3³≈27–31), which was Wilson's/the ×3-lead's and had already been downgraded to "imposed, not
measured" after the waveform retraction, is now **refuted by direct deep test.** This kill is clean, high-
significance, and does not depend on the noisy tail. Wilson owns the kill.

## Result 2 — 7/15 is under real pressure (stated without a hedge)
ε̂ sits at **~9.2×10⁻³ = the no-turnover extrapolation** of `result_logperiodic.md` (ε_∞≈+9.9×10⁻³ ⟺
**S_∞≈0.477, not 0.467**). This is an **independent method** (splitting, r=17–40) **agreeing with and extending the
exact ν-Λ trend, continuous at the r=16 seam**, +3.4σ above anchor. **For the first time in the campaign, the
deepest available data does not point at the constant** — it favors 0.477 over 7/15.

## Result 3 — the only lifeline is a sub-3σ late rollover (7/15's sole surviving path)
ε̂ peaks at r≈36 and dips to r=40, a **+2.3σ** drop-from-peak — **short of the 3σ bar, and it leans partly on
r=40**, the one level where ESS dropped (uniq-frac 0.167 vs the ~0.223 plateau). If real, this rollover **relocates
any turnover to r≥36, far outside [27,31] — so it too kills log₃** — and is the *only* path back to 7/15.

## Guardrail 1 — ESS holds (bars trustworthy through r≈38)
uniq-frac (state diversity after resample; single-level ceiling q(1−e⁻³)≈0.317) is a **flat plateau ~0.223 from
r=17 through r=38** (survival q_s≈0.3335 throughout), dropping only at the last levels (0.217 at r=39, 0.167 at
r=40). So the deep bars are reliable to ~r=38; only r=40 itself is the noisiest datum — which further weakens the
already-sub-threshold rollover hint that partly leans on it.

## The lengthening-transient seam — predicts this shape, and is therefore FLAGGED AS SUSPECT
Wilson's proposed mechanism `Λ_r = S_r·⟨δ_r, Re w⟩` (δ_r the localized deviation field): δ localizes near x≈0 where
the transport barely contracts (`|D(ξ)|→1`), so δ gets **pinned near the trivial character and held** — making the
crossings a **finite transient of a slow approach, not a fixed-period oscillation.** Each successive half-cycle
lengthens **without bound** (not ×3, not any fixed ratio) because the closer δ gets to x=0 the slower it can move —
asymptotically frozen. This predicts **exactly the measured shape**: no turnover in [27,31] ✓, a very late shallow
rollover (the r≈36 hint) ✓, and lengthening that makes *every* finite-period fit fail ✓ (the whole MOON/SAT/×3
history). **⚠️ This is precisely why it is distrusted right now** (Wilson's own flag): a model that survives a
strong negative by *predicting* the negative is either deep or unfalsifiable, and the numerics cannot yet tell
which — it has the **same two-points-of-support, retro-fits-everything shape the ×3 story had three turns before it
was refuted.** Recorded as the thing to be suspicious of, not the thing to believe.

## Why NOT to push to r≈50 (the depth strategy is the same wall's third hat)
If the transport genuinely pins δ near x≈0, **the rollover recedes as you chase it** — each level of depth buys less
because `|D|→1` slows the approach; one could burn r=50 and still get a ~2σ rollover hint at r≈47, forever. The
lengthening that killed the period fits (Probe D/SAT/MOON) **kills the deep-run strategy too.** Chasing the turnover
by depth is structurally the same mistake as chasing the period by depth. **Deeper compute is not the right spend.**

## What actually decides it — the positivity claim (scalar, dodges the growing-state wall) — USER to derive
The whole question is now **"does ⟨δ_r, Re w⟩ become eventually single-signed?"** — simultaneously (a) whether 7/15
holds, (b) why the rollover is late, and (c) whether depth will *ever* resolve it. It is **scalar**, so it dodges the
`3^{r+1}` growing-state wall entirely. The splitting run did not answer it but did something better: it showed **the
claim is on a knife's edge** — ε̂ at 0.477 with a sub-threshold rollover means `⟨δ_r, Re w⟩` is either *just barely*
turning single-signed very late, or not at all. That is exactly the regime where a proof is worth attempting because
the numerics cannot decide it and will not. **Three decisive outcomes:** (i) it closes ⟹ 7/15 proved, rollover
explained as the last sign change, depth-chasing correctly abandoned; (ii) it provably *fails* single-signedness ⟹
7/15 is wrong, ε_∞≈0.477 is real, and the constant was a low-r illusion (a genuine discovery); (iii) genuinely
undecided ⟹ the exact crux inequality is exhibited, a sharper open problem than "run it deeper." **Wilson takes
this (pen); Claude gates the derivation when posted.**

## Status
**S2 DEEP (statistical estimate, NOT proof):** bigint path certified (bit-identical to validated int64); ESS holds
to r≈38. **(1) log₃ ×3-crossing at [27,31] REFUTED** (ε̂ rises z≈5 through the window, no peak there) — clean,
tail-independent, kills Wilson's own ×3 prediction. **(2) 7/15 UNDER PRESSURE:** ε̂ at the **0.477 no-turnover
level** (+3.4σ above anchor), an independent method extending the exact ν-Λ trend — **first time the deepest data
does not point at the constant.** **(3)** only a **sub-3σ (2.3σ) late-rollover at r≈36**, partly leaning on the
noisiest level, keeps 7/15 alive — and if real relocates the turnover past 36, itself refuting log₃. Decision rule:
**INCONCLUSIVE**; shape decisive on the [27,31] refutation. **Next is NOT deeper compute** (rollover recedes as
chased if δ is pinned — the lengthening wall's third hat); it is the **scalar positivity claim `⟨δ_r, Re w⟩`
eventually single-signed**, which decides 7/15, the late rollover, and the depth question at once — **USER to
derive, Claude to gate.** The lengthening-transient seam predicts the observed shape and is therefore **flagged
suspect** (two-points, retro-fits, same shape ×3 had). Not at stake: R1–R30, R80–R82 (exact identities, M-reality).
At stake and now under empirical pressure: the asymptotic value **7/15** (deep data favors 0.477; conditional on an
unconfirmed r≥36 turnover). Epistemic guard held: this is a splitting estimate with an error budget, **never welded
to the exact ε ladder.**
