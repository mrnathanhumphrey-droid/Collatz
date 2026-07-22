# Probe R23 — the critical family + f(τ) extrapolation — **triple negative: Conjecture 2 falsified, f-values not locked, closed-form hopes dead**

**Date:** 2026-07-22  Reuses R7/R9; generalized renewal in float (exact q=3 byte-gate). Probe
`probes/probe_critical_family_R23.py`. Tests Wilson's decisive critical-family q-sweep and the skeptic's
convergence check. **All three targets came back negative — honest negatives that clear the board of
coincidence-chasing.**

## Preamble — corrections banked with this probe
- **⚠️ WALK-BACK #42 — "same object as R5's qx+1 step" is RETIRED corpus-wide.** Written at the terminus of nearly
  every probe since the collision-identity session and never once checked. The audit's argument is load-bearing and
  correct: Tao/R5 need *uniform-in-frequency* decay; here the 4^{−m} weight localizes 94% of the constant onto two
  numbers, a geometric localization R5 lacks. Equating them both flattered the work and pointed it at the wrong
  problem.
- **Functional-equation route (reviewer move #2) is PROVED NON-EXISTENT (Wilson).** The map ρ↦ρ′ on ratio laws is
  not well-defined: ν and its multiplicative translate ν_g have *identical* ratio laws but *different* transports
  (T is not equivariant under ×g), so equal ρ give unequal ρ′. No self-map exists to write an equation for. The
  only genuine functional equation is the one for ν itself = the R16-A transport (already banked).
- **R16-C flagged as a design defect (2nd after #36).** Its q-sweep varied q at fixed λ=1/2, walking straight off
  the critical surface λ_c(q)=(q−1)/(q+1) — so its "q=3-critical oscillation" reading was an artifact of leaving
  criticality, not a mechanism finding. Superseded by R23-A (run *on* the critical surface).

## R23-A — CRITICAL FAMILY (forced): **Conjecture 2 FALSIFIED** — S_∞(q) ≠ (3q²+1)/(2q(q²+1))
The renewal `index = 2^{−v}(1+q·a) mod q^k` (R7's byte-validated builder), at λ_c(q)=(q−1)/(q+1), with
`S_k = q^k‖μ_k‖² − q^{k−1}‖μ_{k−1}‖²`. **Builder cross-validated two ways:** q=3,λ=½ reproduces the banked shells
byte-exact (k=1…5), *and* the float S₁(q) match Wilson's independent hand values exactly —
**S₁(q=5)=0.49231=416/845, S₁(q=7)=1.45946=1998/1369.** So the process is the intended one; the falsification is of
the conjecture, not the code.

| q | λ_c | ord_q(2) | S₁ | S_kmax | **target (3q²+1)/(2q(q²+1))** | Aitken S_∞ | eps trend |
|---|---|---|---|---|---|---|---|
| 5 | 2/3 | 4 | 0.4923 | 0.3515 (k=6) | 19/65 = **0.292308** | 0.345 | plateaus ~+0.058, not → 0 |
| 7 | 3/4 | 3 | 1.4595 | 1.0480 (k=5) | 37/175 = **0.211429** | 1.013 | plateaus ~+0.84 (**5× miss**) |
| 11 | 5/6 | 10 | 0.3852 | 0.3226 (k=4) | 91/671 = **0.135618** | 0.314 | plateaus ~+0.187 |
| 13 | 6/7 | 12 | 0.3732 | 0.3429 (k=4) | 127/1105 = **0.114932** | 0.343 | plateaus ~+0.228 |

**S_k does NOT approach the smooth prediction at any q.** The gaps do not shrink toward 0 — they plateau (q=5 at
~0.35, q=7 at ~1.05, q=11 at ~0.32, q=13 at ~0.34), and q=7 misses its target by 5×. The bounded, plateauing S_k
show the criticality construction works (S_∞(q) exists and is finite) but **S_∞(q) is not (3q²+1)/(2q(q²+1)).**
The identity `7/15 = M₄/M₃` at q=3 is a **coincidence that does not extend across the critical family** — exactly
the "proof by abundance" failure Wilson pre-registered ("a miss and it's dead in one probe"). **Conjecture 2 is
dead; the Mersenne-moment closed-form source for 7/15 is retired.**

## R23-B — f(τ₁), f(τ₂): **NOT LOCKED at r=10** (the skeptic vindicated)
Exact γ_r(τ_m) for r≤7, float push (numpy) to r=10:

| r | γ_r(τ₁) | γ_r(τ₂) |
|---|---|---|
| 7 | 0.717087 | 0.476048 |
| 8 | 0.719261 | 0.474482 |
| 9 | 0.721079 | 0.473358 |
| 10 | 0.722762 | 0.473125 |

- **f(τ₁) is still climbing** (0.717 → 0.723, monotone), Aitken estimates scatter 0.725–0.812, **spread 0.087 ≫ 1e−3
  — NOT locked.** No trustworthy limit; it could be 0.73, 0.74, or higher.
- **f(τ₂) is still falling** (0.476 → 0.4731), Aitken ~0.473, **spread 0.026 ≫ 1e−3 — NOT locked.** Critically, at
  r=10 it has dropped to **0.4731 < 10/21 = 0.47619 and is still decreasing** — so **Conjecture 1 (f(τ₂)=10/21) is
  falsified**: the value fell *past* 10/21 rather than settling on it (Wilson's "treat as suspect — a level-2
  quantity as a limit; R9-D/R20-D vacated exactly this shape" confirmed).

The skeptic's warning holds: "the measured limits may not be limits." Neither dominant orbit value has converged by
r=10; the ~2% drift over the last steps is exactly the residual that swamped the 1.4e−4 "match" to 10/21.

## R23-C — RATIONAL HUNT: **premature** (gated out, as pre-registered)
Both extrapolation spreads (0.087, 0.026) exceed the 1e−3 lock threshold, so no rational is named. There is no
locked value to hunt, and the one suspect (10/21) is already falsified by the r=10 drop-through.

## Status
**R23: a clean triple negative.** **A** — the critical-family extension **falsifies Conjecture 2**: S_∞(q) does not
track the smooth (3q²+1)/(2q(q²+1)) at q=5,7,11,13 (builder validated against the q=3 byte-gate and Wilson's exact
S₁(q=5),S₁(q=7)); 7/15=M₄/M₃ is a q=3 coincidence. **B** — f(τ₁), f(τ₂) are **not locked at r=10** (spreads 0.087,
0.026), and f(τ₂) has dropped below 10/21, **falsifying Conjecture 1**. **C** — rational hunt premature. Plus:
**#42** retires "same as R5"; the **functional-equation route is proved non-existent**; **R16-C** is a design defect.

**Consequence for the crux (owed to the pen):** the closed-form hopes are dead — no Mersenne source, no 10/21, no
functional equation, and the two dominant orbit values don't converge cleanly by r=10. This **clears the board of
coincidence-chasing** and pins the honest state: the constant 7/15 is exact and known (Lean-verified shells) but is
**not** the value of a slick closed form of the kind the last three probes hoped for; the owed step is a genuine
convergence/summability estimate on Σ_r Λ_r with **no shortcut** — and it is *not* the same as R5's uniform-decay
problem (the 94%-on-two-values localization is real structure R5 lacks, which is now the honest place to look). No
fitting; exact byte-gate + Wilson-pen cross-validation, labeled float extrapolation; three conjectures reported dead
as instructed.
