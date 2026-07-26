# Probe CONE — is the positive sign self-generated? — **NO. The flubber/self-generation frame is NOT supported. Lead result (family 5): Φ = C_q(1)/C_q(0) is EVEN in q, so "flip the fluctuation sign" and "reflect" both leave Φ = +d₁ EXACTLY — family 5 starts INSIDE 𝒦, not outside (the sharpest test can't be built from a sign-flip; the falsifier must disrupt the correlation STRUCTURE, not the sign). For the genuine adversarial starts (Φ₀ from −1.0 to 0): every one converges to +d₁, but NOT by entering 𝒦 and remaining — 𝒦 is NOT forward-invariant (interior starts exit; the transport drives the neutral δ_{X=1} to Φ = −0.21, exact-confirmed), and convergence is oscillatory/overshooting, crossing the boundary repeatedly. So the positive sign is carried by CONVERGENCE to the fixed point ν, not manufactured by the transport — exactly the CONE-D vacuity. The sign is a property of ν_r; there is no forward-invariant cone to prove.**

**Date:** 2026-07-25. Probe `probes/probe_cone.py`, log `logs/cone_run.log`. Certified R16-A renewal `T` (build_nu's step, reused): `a → a' = (2⁻ᵛ(1+3a)) mod 3^r`, `v~Geom(½)`; fixed point ν_r (R8 self-similarity). `P` (kill 3|m) = subtract the (N/3)-periodic part (exact, no FFT). `Φ(μ)=C_q(1)/C_q(0)`, `q=Pρ_μ`. Exact rational at r=4,5 (gate, constructions, first steps); float64 for r=4..8 trajectories (20 exact iterations blow up denominators; sign-crossings −1→+0.002 robust to 1e-12; **exact vs float validated to 9 digits**).

## CONE-A — GATE (pass)
| r | Φ(ν_r) vs d₁_r | ‖Tν − ν‖₁ | Φ(uniform) |
|---|---|---|---|
| 4 (exact) | rel 5.2e−14 | **6.0e−37** | q≡0 ✓ |
| 5 (exact) | rel 5.0e−15 | 2.9e−17 | q≡0 ✓ |
| 4..8 (float) | rel ≤1e−13 | ≤4.5e−16 | — |
**ν is the fixed point of T to machine precision; Φ(ν)=d₁; P correct (uniform → q≡0).** Exact-vs-float trajectory (δ_{X=1}, r=5): identical to 9 digits — float licensed.

## LEAD (CONE-B family 5): the sign of q is INVISIBLE to Φ
`Φ = C_q(1)/C_q(0)` is a ratio of quadratics in q ⟹ **even in q**. So:
- **F5a reflected** `ρ_ν(−s)`: q → q(−s), `Φ₀ = +d₁ EXACTLY` (r=4,6,8: +8.643e−3, +7.742e−3, +4.666e−3 = d₁).
- **F5b flip** `ρ_ν − 2q_ν`: q → −q, `Φ₀ = +d₁ EXACTLY`.

**Family 5 — "same envelope, wrong sign" — starts INSIDE 𝒦, not outside.** The "wrong sign" intuition does not map to a Φ<0 start: flipping or reflecting the fluctuation is invisible to the lag-1 autocorrelation. A falsifier can only be built by disrupting the correlation STRUCTURE (F1–F4), never the sign. (This is the informative single result requested to lead — and it says the sign question is intrinsically about *structure of q*, not *sign of q*.)

## CONE-C/D/E — the genuine adversarial starts
| family | Φ₀ (r=8) | settle | stays +? | L1@first-cross | Φ@first-cross |
|---|---|---|---|---|---|
| F1 extremal (high-freq mode) | **−1.0000** | +d₁ | **NO** | 0.66 | +0.024 (overshoot) |
| F2 alt-comb | −0.9998 | +d₁ | **NO** | 0.80 | +0.0035 |
| F3 ×4-anticorrelated | −0.7710 | +d₁ | **NO** (r8) | 0.56 | +0.043 |
| F4a point-mass generic | 0 | +d₁ | **NO** | 1.99 | +3e−19 |
| F4b point-mass δ_{X=1} | 0 | +d₁ | **NO** (r8) | 1.94 | +1e−7 |

Three decisive facts, all pointing the same way:
1. **The transport does NOT manufacture positivity.** The neutral start δ_{X=1} (Φ=0) is driven to **Φ = −0.21** by T (exact-confirmed: steps 0→4 = 0, +3e−6, −0.013, −0.016, −0.213). If the sign were self-generated, T would push toward +; it pushes strongly −.
2. **𝒦 is NOT forward-invariant.** F5a starts at Φ=+d₁ (inside) yet `stays+? = NO` — it exits 𝒦 and returns. The deep-negative starts, once they first touch positive, dip back below 0 (`stays+? = NO` for F1, F2, F3@r8, F4). Trajectories cross the boundary repeatedly.
3. **Convergence is oscillatory/overshooting, not self-restoring.** Φ@first-cross overshoots above d₁ (F1: +0.19, +0.05, +0.024 ≫ d₁≈0.008) then undershoots below 0 and spirals into +d₁ — damped oscillation about ν (the arc's subdominant oscillatory mode), which crosses zero, not a floor-restoration.

## Verdict — the flubber frame is NOT supported (CONE-D vacuity)
Every start converges to `Φ = +d₁` **because ν is the attracting fixed point** (mixing), and every ergodic start inherits `Φ(ν) > 0` asymptotically. That is trivial ergodicity, **not** self-generation. Wilson's CONE-D test — "does Φ turn positive while FAR from ν (manufacture) or only near ν (convergence)?" — resolves as **convergence**: the transport drives neutral/point starts *negative*, keeps 𝒦 non-invariant, and recovers Φ>0 only by mixing to ν, with oscillatory boundary-crossings along the way. The strict death-condition (a permanent Φ≤0 survivor) did not fire — but neither did self-generation. **There is no forward-invariant cone to prove invariant; the positive sign is a property of the stationary measure ν_r, not a dynamically-manufactured/self-restored attractor.** The three deaths (R29, PRODFORM, decay-shelf) plus this one now say: the sign lives in the fixed point, and the target is unchanged — prove `d₁(ν_r) = C_{q_ν}(1)/C_{q_ν}(0) > 0` for the specific ν_r (the AC-LAGS autocorrelation inequality, all terms ≥0 via ρ≥0, margin 0.19%). The self-consistency reframe does not open a cone route.

## Status
**CONE (falsifier, self-generation NOT supported):** GATE pass (ν = fixed point to 1e-17, Φ(ν)=d₁, exact=float to 9 digits). **LEAD: Φ even in q ⟹ family 5 (reflect/flip) starts at +d₁ INSIDE 𝒦** — the "wrong sign" test is void; falsifiers must disrupt q's structure not its sign. Genuine adversarial starts (F1 Φ₀=−1.0, F2 −1.0, F3 −0.77, F4 0) **all converge to +d₁ but NONE enter-and-remain**: (1) T drives neutral δ_{X=1} to Φ=−0.21 (does NOT manufacture +); (2) 𝒦 NOT forward-invariant (interior starts exit); (3) convergence oscillatory/overshooting, crosses zero repeatedly. ⟹ positivity is **carried by convergence to ν, not self-generated** (CONE-D vacuity). **No forward-invariant cone; the sign is a property of ν_r.** Flubber/self-generation frame CLOSED as a proof route; target unchanged = `d₁(ν_r)>0` (AC-LAGS autocorr inequality). Not at stake: R1–R30, R80–R82, all Thread-3 probes. commit pending.
