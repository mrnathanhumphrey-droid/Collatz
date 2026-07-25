# Probe RATIO-2 — is the oscillating component of Re δ̂_r(1) SUBDOMINANT? — **YES at the r≤14 test: |D_14|=0.0103 lands in the pre-registered subdominant band [−0.0104,−0.0026] and is 3.5× below the not-subdominant floor 0.0363, so the "not-subdominant / 7/15-live" branch is REJECTED; the oscillation damps at μ=(D_14/D_10)^¼=0.73 < the dominant rate ~0.90, and the per-period trend collapsed 0.044→0.001. Re δ̂(1)>0 through r=16. This SUPERSEDES RATIO's "7/15 live again" (RATIO's rule conflated a subdominant pair with a dominant one — Wilson withdrew it). Two honest misses: the pre-registered sign pattern +,+,−,+ FAILED (observed −,−,−,−, ρ_r monotone-decreasing r=11→15 = a broad damped trough, not a period-4 dip), and the damping is DECELERATING (μ rose 0.615→0.73), so asymptotic subdominance is supported-but-not-locked. Provability untouched. Leans 0.477 more firmly than RATIO.**

**Date:** 2026-07-25. Probe `probes/probe_ratio2.py`, log `logs/ratio2_run.log`. Exact-side only (build_nu→dlog→|FFT|² profile path, validated at G0). Tests whether the oscillation RATIO found is **dying relative to the dominant mode** — the only version of the question that touches the sign. **A subdominant complex pair produces exactly what RATIO measured (a strictly positive sequence with oscillatory ratio ripples); only a *dominant* complex eigenvalue can flip Re δ̂(1).**

## Primary observable (fit-free)
`D_r := ln ρ_r − ln ρ_{r−1} = ln d1_{r+1} − 2 ln d1_r + ln d1_{r−1}` (2nd difference of ln d1_r), `d1_r := Re δ̂_r(1)`, `ρ_r := d1_{r+1}/d1_r`. `D_r` annihilates the dominant geometric rate exactly (no fit, no window) — it **is** the oscillating component.

## Exact closed form (derived; gives a rational reference and cross-checks the FFT path for free)
`Re δ̂_r(1) = [ 2C(1) − (C(N/3−1)+C(N/3+1)) ] / [ 2(C(0) − C(N/3)) ]`, `N=3^r`, `C(d)=Σ_s ρ(s)ρ((s+d) mod N)`.
(Because `Σ_{k∈prim} cos(2πka/N) = N·[N|a] − (N/3)·[(N/3)|a]`, only lags {0,1,N/3−1,N/3,N/3+1} survive.) Reproduces every banked d1_r exactly.

## R2-A — PRECISION FLOOR (gate) — **PASS by 13 orders of magnitude**
- **A0** profile-path (P) vs closed-form (F), float, r=2..12: `max rel|P−F| = 5.8×10⁻¹⁴` (both paths agree; formula validated).
- **A1** exact-rational reference r≤7: `|float(P) − exact|/|exact|` = **6×10⁻¹⁵ … 2×10⁻¹⁴**, flat in r (no growth).
- **A2** working end: tol(1e-18 vs 1e-26) rel = 2×10⁻¹⁵ (truncation negligible); P-vs-F rel = 6×10⁻¹⁵.
- **Floor `δD ~ 4σ_rel = 9.4×10⁻¹⁴` vs half-signal 0.0026.** GATE PASS — this is a real high-precision result, not a null-from-noise. (The #32/#40 failure mode is excluded.)

## R2-B — the ladder, extended to r=16 (build_nu reached r=16; auditable)
| r | 11 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|
| d1_r (×10⁻³) | 3.2111 | 2.9636 | 2.6962 | 2.4409 | 2.1871 | 1.9392 |
**Re δ̂(1) > 0 at every r=2..16** (P-vs-F agreement ≤6×10⁻¹⁴ at each new level). `d2,d3` also banked (see log; needed for shape).

## R2-C — the dip test (the decider)
| r | ρ_r | D_r |
|---|---|---|
| 11 | 0.9229 | +0.0343 |
| 12 | 0.9098 | −0.0143 |
| 13 | 0.9053 | −0.0050 |
| 14 | 0.8960 | **−0.0103** |
| 15 | 0.8867 | −0.0105 |

**Pre-registered predictions vs observed:**
| quantity | subdominant/dying | not-subdominant | **observed** | verdict |
|---|---|---|---|---|
| D_14 | −0.0052, band [−0.0104,−0.0026] | \|D_14\|≥0.0363 | **−0.0103** | **IN subdominant band** (upper-mag edge); 3.5× below not-sub floor ⟹ **not-subdominant REJECTED** |
| signs D_12..D_15 | +,+,−,+ | — | **−,−,−,−** | **MISS** (see below) |
| dip location | r=14 | — | broad trough r=12–15 | period lengthened |

`μ = (D_14/D_10)^¼ = 0.73`, `ρ_c = μ·ρ₁ ≈ 0.67` — the oscillation decays at ~0.73/level, **faster than the dominant mode ~0.90** ⟹ **subdominant.** But `μ` ROSE from 0.615 (D_6→D_10) to 0.73 (D_10→D_14): **the damping is decelerating.**

## R2-D — shape: is ρ converging? — **YES (trend collapsed)**
Per-period mean of D: **r=7..10 = +0.0442** (banked +0.0444 ✓), **r=11..14 = +0.0012.** The upward trend that inflated the ratios is gone ⟹ the earlier 3:1 sign asymmetry was **trend+oscillation**, and the Aitken ~0.90 limit is meaningful (ρ settling from above, not still climbing). Tail factor `1/(1−ρ)≈10` is **not** understated (Wilson's R2-D worry not realized — ρ is decreasing, not rising).

## R2-E — record check (no compute): the banked ρ≈0.984 is NOT a competing rate
- `result_epsilon_11.md` (line 18): model `ε_k ≈ scale·ρ^k·cos(kθ+φ)`, ρ=0.984 — **fitted to ε_k as a DECAYING quantity ⟹ presupposes ε_∞=0 = 7/15.** Never an independent rate.
- `result_F1_period9.md` F1-A: where 0.984 appears as a Λ-ratio (0.987, 0.989) it is measured in the **flat transient window r=8–11** (Λ near its maximum). The actual deep Λ decay (`result_crossing.md`, r=11–16) is **~0.90 — "Λ decays FASTER than 0.984^r"** — which is exactly RATIO-2's d1 rate.
- **⟹ No 0.984-vs-0.90 rate tension. RATIO-2's ~0.90 IS the real rate**, and the crossing Λ-rate, unified. The 0.984 either presupposed 7/15 or was a transient-window artifact.

## What this settles, and what it does NOT
- **Settles (the question asked):** the oscillating component is **subdominant** at r≤14 — three independent reads: (i) |D_14|=0.0103 in the subdominant band and 3.5× below the not-subdominant floor; (ii) μ=0.73 < dominant 0.90; (iii) trend collapse 0.044→0.001. A subdominant pair **cannot flip** the dominant positive mode ⟹ **RATIO's "7/15 live again" is superseded** (RATIO's withdrawn rule conflated *a* complex pair with the *dominant* one). Re δ̂(1)>0 holds through r=16, and 7/15 now needs a mechanism nobody has produced (a NEW dominant-scale negative excursion past r=16 — zero evidence, and now a dying-amplitude argument against).
- **Does NOT settle (honest):**
  1. **Sign-pattern miss.** Predicted +,+,−,+; observed −,−,−,− (ρ_r monotone-decreasing r=11→15). The oscillation is **not a period-4 isolated dip but a broad, shallow, damped trough** — the period lengthened. This does not rescue 7/15 (amplitude dying) but the specific "period-4 complex pair" model is wrong; the shape is a **freezing/lengthening transient** (consistent with the lengthening-transient seam: damping + lengthening = finite, asymptotically frozen).
  2. **Decelerating damping.** μ rose 0.615→0.73. If μ keeps creeping toward ρ₁≈0.90, subdominance becomes marginal asymptotically. Two μ points can't tell (#30). Subdominance is **locked at r≤14, supported-but-not-locked asymptotically.**
  3. **Provability untouched.** Even a clean dying-oscillation leaves the real burden: an **analytic lower bound on Re δ̂_r(1) from the construction of ν_r.** Bochner/positivity still insufficient. RATIO-2 is evidence, not proof.

## Status
**RATIO-2 (exact-side, precision-gated):** floor 9.4×10⁻¹⁴ ≪ signal 0.0052 (exact-rational reference agrees to ~1e-14; build_nu to r=16). **The oscillating component of Re δ̂(1) is SUBDOMINANT at r≤14** — |D_14|=0.0103 in the pre-registered subdominant band, 3.5× below the not-subdominant floor 0.0363 (**REJECTED**); damping μ=0.73 < dominant ~0.90; per-period trend collapsed 0.044→0.001 (**ρ converged, R2-D**). Re δ̂(1)>0 through r=16. **SUPERSEDES RATIO's "7/15 live"** — a subdominant pair cannot flip the dominant positive mode. **Two misses:** pre-registered sign pattern +,+,−,+ FAILED (observed −,−,−,− = broad damped trough, period lengthened, freezing-transient shape); damping decelerating (μ 0.615→0.73) ⟹ asymptotic subdominance supported-but-not-locked. **R2-E:** banked 0.984 presupposed 7/15 (fit to ε_k) or was a transient-window Λ-ratio; RATIO-2's ~0.90 is the real rate (= crossing deep Λ-rate), **no rate tension.** Leans 0.477 more firmly than RATIO/MODES. **Provability untouched** (analytic lower bound on Re δ̂_r(1) from ν_r construction remains the burden). Not at stake: R1–R30, R80–R82, g_r/d1 ladder now to r=16, Λ^unif closed form, MODES' Re δ̂(1)>0.
