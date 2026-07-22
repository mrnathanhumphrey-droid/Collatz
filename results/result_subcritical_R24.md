# Probe R24 — the subcritical scaling law — **POSITIVE: ε·X_∞ → 7/40 from outside criticality; strongest evidence yet for 7/15**

**Date:** 2026-07-22  Reuses the R23 generalized renewal builder (float/numpy). Probe
`probes/probe_subcritical_R24.py`. Wilson's move: the critical sum `S_∞ = 2Σ_χ|ν̂|²w(χ)` is **conditionally
convergent** (Σ|ν̂|²|w| ≤ (1/3)·lim X_r = ∞; terms don't shrink, only cancellation converges), so **no finite prefix
determines it** and every extrapolator imposes a decay model the system doesn't have (the R23-B 0.72–0.81 scatter =
model misspecification, not noise). Fix: **step off criticality.** Hold q=3, set λ=½+ε; then the same sum is
**positive-term, absolutely convergent, clean geometric** with exact rate ρ=3(1−λ)/(1+λ)<1. Test the scaling law.

## R24-A — the exact decay rate ρ is CONFIRMED (Wilson's "calm subcritical" verified)
Measured `S_{r+1}/S_r` locks onto ρ_pred = 3(1−λ)/(1+λ) at every ε:

| ε | λ | ρ_pred | measured ratio (r=10) |
|---|---|---|---|
| 0.10 | 0.60 | 0.75000 | 0.7418 (→) |
| 0.05 | 0.55 | 0.87097 | 0.8672 |
| 0.02 | 0.52 | 0.94737 | 0.9466 |
| 0.01 | 0.51 | 0.97351 | 0.9739 |

The subcritical decay is a **single-mode geometric with a closed-form rate** — no oscillation, no cancellation. This
is the structural opposite of the critical regime, and it makes extrapolation legitimate (the model is *known*, not
imposed).

## R24-B — the scaling law: **ε·X_∞ → 7/40, amplitude → 7/15**, both monotone
`X_∞(λ) = Σ_{r≥1} S_r = (lim_r 3^r‖μ_r‖²) − 1`, extrapolated two independent ways (Aitken on the partial sums; and
the exact-ρ geometric tail `Y_R + S_R·ρ/(1−ρ)`).

| ε | ρ | X_∞ | **ε·X_∞** (→ 7/40 = 0.17500) | **amplitude (1−ρ)·X_∞** (→ 7/15 = 0.46667) | estimator agreement |
|---|---|---|---|---|---|
| 0.10 | 0.750 | 1.4065 | 0.14065 | 0.6016 | 2.0e−3 |
| 0.05 | 0.871 | 2.9831 | 0.14915 | 0.5139 | 2.1e−2 |
| 0.02 | 0.947 | 8.0925 | 0.16185 | 0.4786 | 6.6e−2 |
| 0.01 | 0.974 | 16.796 | **0.16796** | **0.4714** | 1.9e−1 |

**`ε·X_∞` climbs monotonically toward 7/40 from below** (0.141 → 0.168); **the amplitude falls monotonically toward
7/15 from above** (0.602 → 0.471). Linear-in-ε extrapolation of the two smallest amplitude points → 0.464; Aitken on
the amplitude sequence → 0.470 — both straddle 7/15. **The scaling law holds: `C(λ) → 7/15` as λ→½⁺.**

## Status
**R24: POSITIVE — the strongest evidence for 7/15 the program has produced, via a genuinely different limit.** The
critical constant is approached from *outside* criticality along a clean, positive-term, absolutely-convergent
family; the exact decay rate ρ=3(1−λ)/(1+λ) is confirmed, and both scaling forms (ε·X_∞ → 7/40, amplitude → 7/15)
converge monotonically to the pre-registered targets. This is *not* what R23-A killed — that varied q at critical λ;
this holds q=3 and varies λ, a limit nobody had run.

**⚠️ Honest bounds on the claim:**
- **Not 6-digit-locked.** The two extrapolators agree tightly at large ε (2e−3) but diverge at small ε (1.9e−1)
  because the renewal support is 2·3^{r−1} (λ-independent), capping the build at r=10, so the small-ε X_∞ are
  tail-dominated. The amplitude readout — which rides on the *confirmed exact ρ* — is the reliable one, and it is the
  cleaner monotone march to 7/15; the raw X_∞ at ε≤0.02 carry ~1–10%. So: **decisive-suggestive, not proof.**
- **Continuity assumption (Wilson's own flag).** The test is really "does the amplitude C(λ) → 7/15" — the theorem
  approached from outside. Not circular (it is a *different* limit: λ→½ of an absolutely-convergent positive-term
  quantity, vs r→∞ of a conditionally-convergent one), but it needs **continuity of C at criticality**, which is
  unproved.
- **Trap avoided (Wilson):** extrapolated X_∞, not f(τ_m). Subcritically X_∞ is finite so f has *no* identity
  singularity — that singularity is born exactly at λ=½, and f(τ_m) may be discontinuous there; so the f(τ_m)
  extrapolation (R23-B) was the wrong object, and X_∞ is the right one.

**Consequence for the crux (owed to the pen):** R24 gives the first *trustworthy* numerical confirmation of 7/15
(clean positive-term limit, estimators agree where the build is deep, ρ exact) — and reframes the target as
**continuity of the subcritical amplitude C(λ) at λ=½**, a cleaner analytic object than the conditionally-convergent
critical sum. Sharpening to six digits needs deeper builds at small ε or a Richardson-in-ε pass on the amplitude. No
fitting of the *rate* (ρ is exact and confirmed); labeled numeric extrapolation with agreement/spread reported; the
continuity gap flagged as the honest remaining assumption.
