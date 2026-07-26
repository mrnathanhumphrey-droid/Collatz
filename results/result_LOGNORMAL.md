# RESULT — LOGNORMAL check: REJECTED; log-variance SATURATES (not linear); Σ|π̂|⁴ ≠ U² (2026-07-26)

**Probe:** `probe_lognormal.py`. Wilson's overdetermined test: model `L = ln|π̂(a)| ~ N(μ, σ²)` with `σ² ~ ck`; fit
(μ,σ²) from the direct mean+variance of L over the units, then predict `E|π̂|²=e^{2μ+2σ²}`, `E|π̂|⁴=e^{4μ+8σ²}` and
check against measured. Pass ⟹ 2-parameter closed-form magnitude distribution; miss ⟹ log-normality wrong.

## Log-normality REJECTED (3 independent ways)
| k | μ | σ²=Var(L) | skew(L) | exkurt(L) | E2meas/E2pred | E4meas/E4pred | Σ\|π̂\|⁴ |
|---|---|-----|------|------|------|------|------|
| 3 | −1.93 | 0.107 | −0.18 | −0.71 | 0.99 | 0.89 | 0.0166 |
| 4 | −2.56 | 0.233 | −0.63 | −0.40 | 0.91 | 0.53 | 0.0065 |
| 6 | −3.78 | 0.446 | −1.36 | +3.48 | 0.76 | 0.23 | 0.0011 |
| 9 | −5.52 | 0.468 | −0.68 | +1.28 | 0.86 | 0.59 | 0.00009 |
| 13 | −7.83 | 0.563 | −0.57 | +1.44 | 0.90 | 1.75 | 0.0000 |

1. **The moments do not close.** `E4meas/E4pred` swings 0.89 → 0.23 → 1.75 across k, never stable at 1; even
   `E2meas/E2pred ≈ 0.87` (not 1). The overdetermined check MISSES — log-normality is wrong. (Fails already at the
   *untruncated* k=3,4 where E4rat = 0.89, 0.53 — not a stationary_trunc artifact.)
2. **L is not Gaussian.** Persistent **negative skew ≈ −0.65** and **positive excess kurtosis ≈ +1.4**, stable from k≈8
   onward (k=6,7 show a transient spike to skew −1.4, exkurt +3.5). The log-magnitude spectrum converges to a **stable
   non-Gaussian shape** — a heavy LEFT tail (anomalously suppressed frequencies), consistent with VALPROFILE's
   low-tail-dominated geometric mean.
3. **σ² SATURATES — it is NOT linear in k.** Var(L) = 0.11, 0.23, 0.30, 0.45, 0.47, 0.44, 0.47, 0.49, 0.52, 0.54, 0.56
   — rises k=3..6 then **plateaus ~0.5** (linear fit `σ²=0.039k+0.10` only R²=0.82). ⚠️ **Correction to VALPROFILE:**
   the "log-variance grows ~linearly in k" was an artifact of comparing μ to the *model* typical `√k·3^{−k/2}`, not the
   actual quadratic mean. `μ ≈ −0.58k` and `½ln E2 ≈ −0.55k`, so `μ − ln(√k·3^{−k/2}) ≈ −0.03k` grew linearly — but
   that gap is μ-minus-a-model, NOT the true variance. The **directly measured Var(ln|π̂|) plateaus at ~0.5.** So the
   random-walk-over-levels premise is not supported: the per-level variance increments shrink to zero, they don't
   accumulate linearly.

## Bonus catch — Σ|π̂|⁴ ≠ U² (object mismatch avoided)
`Σ_{3∤a}|π̂(a)|⁴` = 0.0166 → 0.0000 (decaying ~geometrically), **NOT** the level-invariant `U²=0.29754`. So the
identification `E|π̂|⁴ = U²/φ` is WRONG: **U²=0.29754 is the ρ (dlog/channel) spectrum's ℓ⁴, not π̂'s** — the
additive(π̂)-vs-multiplicative(ρ̂) seam again. Measuring it directly (rather than assuming) caught the mismatch.

## Verdict
Wilson's pre-authorized fallback fires: **log-normality is wrong.** And the log-variance fact it rested on is weaker
than stated — `Var(ln|π̂|)` **saturates ~0.5**, it does not grow linearly. What survives is a genuine (if modest)
distributional fact: **the π̂ log-magnitude spectrum converges to a stable non-Gaussian shape** (skew ≈ −0.65,
exkurt ≈ +1.4, Var ≈ 0.5). This constrains the spectrum but with a shape, not two parameters — and it does **not**
bridge to the channels (the ℓ⁴ that defines γ_r(k) is the ρ-spectrum U², a different object from π̂'s decaying ℓ⁴). The
norm→channel bridge Wilson hoped for does not form via log-normality. Not at stake: VALPROFILE (denominator theorem
stands — that's the day's keeper), NORMCHECK, RECENTER, CHANNEL_ID, R1–R30. Cheap (8.5s).
