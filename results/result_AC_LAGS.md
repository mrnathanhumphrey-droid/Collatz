# Probe AC-LAGS — the five-lag structure of the dominant mode Re δ̂_r(1) — **the G−K reframing's PREMISE is wrong (G,K are not in (0,1) — they are large NEGATIVE and diverge linearly at −0.652/level, from −1.03 at r=2 to −10.25 at r=16, because C(1)<C(N/3) and the denominator C(0)−C(N/3)→0). The split is exact but ill-conditioned and does NOT isolate the signal. The real structure: the autocorrelation has THREE near-identical peaks at k=0, N/3, 2N/3 (the 3-fold fiber structure), whose local bump shapes are identical except the center (C(0)=1 vs C(N/3)=0.919 at r=16). Since C(0)−C(N/3)>0 strictly, the inequality the pen must prove is `C(1) > ½[C(N/3−1)+C(N/3+1)]` — the lag-1 autocorrelation exceeds the mean of the two lags flanking the secondary peak — which holds at r=16 by a relative margin of only 0.19%, tight and tightening.**

**Date:** 2026-07-25. Probe `probes/probe_ac_lags.py`, log `logs/ac_lags_run.log`. Pure measurement, r=2..16, on the validated build_nu/dlog/profile path. No new operator, no new depth. `C(d)=Σ_s ρ(s)ρ((s+d) mod N)`, N=3^r. Purpose (Wilson): give the pen the **shape of the inequality**, not decide anything. rho_12..16 + shape grids dumped to scratchpad for reuse.

## AC-A — GATE (pass)
5-lag closed form vs profile-path d1, r=2..16 (exact rationals r≤7): **worst rel = 7.8×10⁻¹³ [PASS]**. **C(0) > C(N/3) strictly at every r** (denominator positive; no index-3-coset degeneracy).

## AC-B — the normalized table, and the premise failure
`Re δ̂_r(1) = G_r − K_r`, `G_r=[C(1)−C(N/3)]/[C(0)−C(N/3)]`, `K_r=Δ²C(N/3)/(2[C(0)−C(N/3)])`. Selected rows (C·/C0):

| r | C(1)/C0 | C(N/3)/C0 | C(N/3−1)/C0 | C(N/3+1)/C0 | G_r | K_r | G−K = d1 |
|---|---|---|---|---|---|---|---|
| 2 | 0.32381 | 0.66667 | 0.22857 | 0.38095 | −1.0286 | −1.0857 | +5.714e−2 |
| 8 | 0.14584 | 0.85829 | 0.13963 | 0.15074 | −5.0275 | −5.0321 | +4.666e−3 |
| 12 | 0.10667 | 0.89660 | 0.10510 | 0.10762 | −7.6397 | −7.6426 | +2.964e−3 |
| 16 | 0.08406 | 0.91859 | 0.08359 | 0.08421 | −10.2512 | −10.2532 | +1.939e−3 |

**Wilson's premise ("G_r should sit in (0,1)") is FALSE.** It assumed C(1) > C(N/3); measured, **C(1) < C(N/3)** (near lag far below the fiber lag). So the numerator C(1)−C(N/3) < 0 and G < 0. And C(N/3)/C0 **→ 1** monotonically (0.667→0.919), so the denominator C(0)−C(N/3) → 0 and both G, K → −∞. d1 = G−K is a tiny positive difference of two diverging negatives.

## AC-C — cancellation depth (deepening, but not yet exhausting r=16)
"Digits of G and K that coincide" grows monotonically **1.28 (r=2) → 3.72 (r=16)**, ~+0.08/level, still climbing. Float floor on G−K reaches **5.3×10⁻¹³ relative** at r=16 — the exact side still has ~13 digits on G−K, so resolution is not the problem; the ill-conditioning is *why the split is the wrong object*, not a numerical wall.

## AC-D — which term moves: NEITHER converges (pre-registered outcome 3)
`dG` and `dK` are **constant at ≈ −0.652/level** from r=4 on (−0.645, −0.652, −0.656, …, −0.652, −0.653, −0.654). **G and K both drift linearly to −∞**, locked together (G−K = d1, positive, decaying ~0.90). This is outcome 3: **the G−K split isn't the right decomposition** — both terms diverge and the signal is their vanishing difference. The bound cannot be attempted on G or K individually.

## AC-E — fiber-lag antisymmetry: alternating, decaying, and a diagnostic (not a driver)
`A_r = [C(N/3+1) − C(N/3−1)]/C(0)` **alternates sign** and decays at **~0.68/level**: +0.152, −0.092, +0.054, −0.037, +0.024, …, +6.2×10⁻⁴ (r=16); relative to C(1) it shrinks **0.47 → 0.0073**. C is *not* locally symmetric about N/3 but is becoming so. **A_r does NOT enter the d1 numerator** (which uses the *sum* C(N/3−1)+C(N/3+1)), so it is a diagnostic, not a driver — but its alternating ~0.68 decay is plausibly the same object as RATIO-2's subdominant oscillation (μ≈0.73).

## AC-F — the full shape: THREE replica peaks (the structural key)
Fine windows (r=16, C(k)/C0, offsets −5..+5):
- **@ k=0:** 0.088, 0.099, 0.142, 0.055, 0.084, **1.000**, 0.084, 0.055, 0.142, 0.099, 0.088
- **@ k=N/3:** 0.088, 0.099, 0.142, 0.055, 0.0836, **0.919**, 0.0842, 0.055, 0.142, 0.099, 0.088
- **@ k=2N/3:** (mirror of N/3)

**The local bump shape is IDENTICAL at all three peaks (flanks match to ~4 digits); only the CENTER differs** — C(0)=1.000 vs C(N/3)=C(2N/3)=0.919. So `C(k) ≈ [3-fold quasi-periodic envelope with peaks at 0, N/3, 2N/3] + [a self-correlation spike at k=0 only]`. This is the level-lift fiber structure (shift by N/3=3^{r−1} maps fibers to fibers) imprinted on the autocorrelation, and C(N/3)/C(0)→1 is the profile becoming fiber-uniform. Lag 1 and the N/3-neighborhood are **not** generic points on a smooth decay — they are the +1 offsets of the primary and secondary peaks of a 3-peak structure. (Full 200-pt grids saved: `scratchpad/ac_shape_r{12,16}.npy`.)

## The inequality, restated for the pen
Because C(0)−C(N/3) > 0 strictly, **`sign(Re δ̂_r(1)) = sign( C(1) − ½[C(N/3−1)+C(N/3+1)] )`** = sign of *(offset-1 flank of the primary peak) − (mean offset-1 flank of the secondary peak).* Given the peaks are replicas, this is a **second-order comparison of the same local bump at the primary vs secondary peak** — the two differ only because the primary carries the k=0 self-spike. At r=16 the inequality holds by numerator `2C(1)−C(N/3−1)−C(N/3+1) = +3.16×10⁻⁴·C0`, i.e. **C(1) exceeds the flank-mean by a relative margin of only 0.19%** (0.084055 vs 0.083897), and the absolute margin decays ~0.90/level. **This replaces the divergent G−K split with a bounded, structural, tight inequality.** Both numerator (the margin) and denominator (C0−C(N/3)) → 0 at ~0.90; the sign lives entirely in whether the margin stays positive.

## Status
**AC-LAGS (measurement, decides nothing by design):** GATE pass (7.8e-13; C(0)>C(N/3) strict). **Wilson's G−K reframing premise FALSIFIED** — G,K are large negative, diverge linearly at −0.652/level (−1.03→−10.25, r=2→16), because C(1)<C(N/3) and C(N/3)/C0→1; the split is exact but ill-conditioned and isolates nothing (AC-D outcome 3, neither term converges). **Real structure (AC-F): the autocorrelation has 3 replica peaks at 0,N/3,2N/3** (fiber structure), identical local bumps, centers 1.000/0.919/0.919; lag 1 and N/3±1 are the +1 offsets of the primary/secondary peaks. **The inequality for the pen: `C(1) > ½[C(N/3−1)+C(N/3+1)]`**, holding at r=16 by a 0.19% relative margin (tight, tightening ~0.90/level). AC-E: fiber-lag antisymmetry A_r alternates, decays ~0.68 (diagnostic, not in the numerator; likely RATIO-2's oscillation). AC-C: cancellation deepens 1.3→3.7 digits, resolution ample to r=16. Not at stake: R1–R30, R80–R82, d1 ladder to r=16, MODES' Re δ̂(1)>0, RATIO-2's subdominance. This hands the pen a bounded structural inequality in place of the divergent G−K near-cancellation.
