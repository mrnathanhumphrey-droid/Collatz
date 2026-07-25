# Probe POINCARE — Wilson's pen reduction verified, and the dilation mechanism tested — **the Poincaré form `Re δ̂_r(1) = 1 − ‖Δq‖²/(2‖q‖²)` (q = fiber-fluctuation profile Pρ) is VERIFIED to ~1e-13 at r=12..16; the Rayleigh quotient R=‖Δq‖²/‖q‖² creeps UP toward the white threshold 2 (1.9941→1.9961), margin 2−R=2·d1 shrinking ~0.90/level = MARGINAL BY CONSTRUCTION (no slack, the threshold sits exactly where a flat spectrum lands). Wilson's parity falsifier is DEGENERATE: Num_even=Num_odd=Num/2 to machine zero, forced by the m↦N−m conjugate symmetry (flips parity, preserves Ĉ·cos). The non-degenerate dilation-band cut (m'=2⁻¹m) is NOT the tight-cancellation organizer (residue/band 0.33→0.54, weakening, vs the cos-sign band's 165×) — but the big-modulation (low-m') band is net positive and dominates, so the dilation weakly tilts the SIGN the right way. The tight cancellation is intrinsically the cos-frequency near-symmetry.**

**Date:** 2026-07-25. Probe `probes/probe_poincare.py`, log `logs/poincare_run.log`. No build_nu (uses `scratchpad/rho_12..16.npy`). Verifies Wilson's pen reduction and runs his cheap falsifier + its non-degenerate correction.

## (1) The Poincaré form — VERIFIED (the solid deliverable)
Since P (kill 3|m modes) is a Fourier multiplier by an indicator, `1_{3∤m}|ρ̂|² = |1_{3∤m}ρ̂|²`, so the projected autocorrelation IS an autocorrelation: `PC = q⋆q̃`, `q = Pρ` (fiber-fluctuation profile). Target `(q⋆q̃)(1) = Σ_k q(k)q(k+1) > 0`. With `Σq(k)q(k+1) = ‖q‖² − ½‖Δq‖²`:
$$\operatorname{Re}\hat\delta_r(1) = 1 - \frac{\|\Delta q\|^2}{2\|q\|^2}, \qquad \operatorname{Re}\hat\delta_r(1)>0 \iff \|\Delta q\|^2 < 2\|q\|^2.$$

| r | ‖q‖² | ‖Δq‖² | Rayleigh R | 1−R/2 | banked d1 | rel | 2−R = 2·d1 |
|---|---|---|---|---|---|---|---|
| 12 | 8.824e−7 | 1.760e−6 | 1.994073 | 2.963565e−3 | 2.963565e−3 | 1.5e−13 | 5.927e−3 |
| 14 | 9.832e−8 | 1.962e−7 | 1.995118 | 2.440865e−3 | 2.440865e−3 | 2.4e−13 | 4.882e−3 |
| 16 | 1.095e−8 | 2.186e−8 | 1.996122 | 1.939225e−3 | 1.939225e−3 | 5.6e−13 | 3.878e−3 |

Both `1−R/2` and `lag1/‖q‖²` reproduce the banked d1 to ~1e-13. **A discrete Poincaré inequality on the fiber-fluctuation subspace.** The Rayleigh quotient ranges over [0,4]; **the value 2 is exactly the white / median-at-N/4 configuration**, and R sits just under it, **creeping UP** (1.9941→1.9961) as the margin `2−R = 2·d1` decays ~0.90/level. **MARGINAL BY CONSTRUCTION** — not a spectral gap with room, a threshold sitting precisely where a flat spectrum lands; the 0.19%/d1 excess is how far under 2 we are. **Any proof strategy that hopes for slack is ruled out by this** — a hard constraint to bank on future attempts.

## (2) Parity falsifier — DEGENERATE (forced by conjugate symmetry)
`Num_even = Num_odd = Num/2` to machine zero (r=16: even−odd = 0.0e+00; r=12: 1.4e−17). **Proof:** `m↦N−m` is a bijection of {3∤m} that flips parity (N=3^r odd), preserves `Ĉ(m)=Ĉ(N−m)` (even spectrum) and `cos(2πm/N)` (even) — so it maps the even-m sum onto the odd-m sum term-by-term. **Any parity cut is pinned to Num/2.** Parity cannot organize the cancellation; Wilson's parity falsifier is degenerate (balanced by symmetry, not by dynamics).

## (3) Dilation-band cut (non-degenerate) — NOT the tight-cancellation mechanism, but weakly sign-correct
The transport modulation `|ĝ(m)|² = 1/(5−4cos(2π m'/N))`, `m'=2⁻¹m mod N`, low-passes in the **dilated** coordinate m' (dlog(−2)=2⁻¹=(3^r+1)/2, the measured translation). Split Num by `|m'|<N/4` (modulation large) vs `>N/4` (small) — a conjugation-**invariant** cut (m'↦−m'), so *not* symmetry-forced:

| r | Num | big-mod (\|m'\|<N/4) | small-mod (\|m'\|>N/4) | \|resid\|/\|big\| | mass(big) |
|---|---|---|---|---|---|
| 12 | +1.390e−3 | +4.184e−3 | −2.794e−3 | 0.332 | 0.5047 |
| 14 | +1.148e−3 | +2.584e−3 | −1.436e−3 | 0.444 | 0.5028 |
| 16 | +0.914e−3 | +1.705e−3 | −0.791e−3 | 0.536 | 0.5016 |

Residue is 33–54% of the band and **weakening** — nowhere near the cos-sign band's 165× (`|resid|/|pos|=0.006`, SPECTILT). By Wilson's criterion (halves comparable to sum → route dead), **the dilation does NOT organize the tight cancellation.** But the big-modulation (low-m') band is **net positive and dominates** (so the residue = d1 > 0): the modulation tilts the sign the right way — a weak, reproducible signature, not the mechanism. The tight 165× cancellation is intrinsically the **cos-frequency near-symmetry** (⟨cos⟩≈0, the nearly-N/4-symmetric spectrum of SPECTILT); the dilation only imprints weakly on which way the residue falls.

## Status
**POINCARE (verify + falsify):** ⭐**Poincaré form VERIFIED** — `Re δ̂_r(1)=1−‖Δq‖²/(2‖q‖²)`, q=Pρ, matches banked d1 to ~1e-13 (r=12..16); target is the lag-1 autocorrelation of the fiber-fluctuation profile, `Σ_k q(k)q(k+1)>0`. **Rayleigh R creeps UP toward 2 (1.9941→1.9961); margin 2−R=2·d1 shrinks ~0.90 — MARGINAL BY CONSTRUCTION (threshold = white config, no slack). Banked as a hard constraint on future proofs.** ⚠️**Parity falsifier DEGENERATE:** Num_even=Num_odd=Num/2 to machine zero, forced by m↦N−m (flips parity, preserves Ĉ·cos) — parity can't organize the cancellation. ⚠️**Dilation-band cut (non-degenerate) NOT the tight-cancellation organizer:** residue/band 0.33→0.54 (weakening) vs cos-sign 165×; but big-modulation band net-positive dominates ⟹ dilation weakly tilts the SIGN correctly. Tight cancellation = cos-frequency near-symmetry, dilation = weak sign imprint. Both mechanism-cuts (parity, dilation) ruled out as the *organizer*; the sign obstruction remains the cos near-symmetry (SPECTILT), now equivalently the Poincaré threshold R<2. Not at stake: R1–R30, R80–R82, d1 ladder to r=16, MODES/RATIO-2/AC-LAGS/SPECTILT. rho_12..16 in scratchpad. commit pending.
