# Probe R1 — the renewal/loss spine (thread 3) — raw dumps on the frozen instrument

**Date:** 2026-07-18  CPU. Probe `probes/probe_thread3_R1.py`, log `logs/probe_thread3_R1_log.txt`,
TSVs `outputs/thread3_renewal_ledger.tsv` + `outputs/thread3_mass_sequences.tsv`. ONE frozen instrument
(`build_M_gen`, the pair operator), iterated. No fitting anywhere; raw dumps are the deliverable.

**Renewal variable (pinned from R7/R16):** Y_k = 3^k·‖π_k‖² (the "3" is universal in q), a_k = 1^T Mᵏ v₀,
v₀=δ(1,1,0). c_k = Y_k − Y_{k−1} (per-scale loss). X_k (accumulation) = Y_k = Σ_{j≤k} c_j.

## R1-D — bookkeeping welds (do this first, per the worksheet): **CONFIRMED**
- **(i) Cross-era weld EXACT.** Exact-rational c_k on the frozen build: **c₁ = 2/3, c₂ = 10/21** — identical to
  the Era-6 R70 sequence (`c₁==2/3: True`, `c₂==10/21: True`). c₃ = 31370/67963 = 0.46157 (the exact
  continuation). The frozen instrument reproduces Era-6's constants byte-for-byte. Float L=3 iteration agrees.
- **(ii) 7/45 normalization (frozen def, transcribed verbatim from `c_seven_forty_fifth_derivation.py` / Paper 5):**
  > ‖d_{k+1}‖² = Σ_{r'} π_{k+1}(r')² − (1/3)·Σ_r π_k(r)²   [R74 identity];   ‖d_{k+1}‖² ~ c·(1/3)^k, **c = 7/45**.

  **Bridge (zero weight; the pen derives):** 3c = 3·(7/45) = **7/15** = the flat per-scale level of c_k. Verified
  exact as Fractions. The 7/45 is the Fourier/Plancherel-side price tag of the 7/15 real-side flat level.

## R1-A / R1-B — accumulation vs saturation: **the boundary contrast holds** (with two honest deviations)
| q (L) | X_k behavior | X_30 | c_k behavior | verdict |
|---|---|---|---|---|
| **3** (L=3) | **linear accumulation** | **14.578** | oscillates in 0.44–0.46 band around 7/15 | **driftless accumulation** |
| 5 (L=2) | **saturates** | 1.234 (flat) | → 0 | geometric saturation |
| 7 (L=2) | **saturates** | 1.370 (flat) | → 0 | geometric saturation |

**The boundary statement is clean:** at q=3 the corpus recursion loses contraction — X_k grows without bound
(14.6 by k=30, and climbing linearly); at q=5,7 it saturates (X_k flat, per-scale c_k → 0). Flatness at q=3 IS
the boundary — exactly R16's "diverges linearly at q=3, converges at q≥5," now on the one frozen instrument.

**Deviation 1 (q=3 asymptotic level = L-truncation, not a fit).** The frozen L=3 c_k tail sits at ~0.437–0.442,
**−2.5% below 7/15** (X_k−(7/15)k drifts 1.20 → 0.578 over k=1..30, i.e. frozen-L=3 slope ≈ 0.447 < 7/15). This
is the depth-L truncation: for k>L the frozen operator underestimates deep collisions, so it accumulates
slightly slower than the true measure. The **exact** constants (k≤3: 2/3, 10/21) are exact; the asymptotic level
approaches 7/15 as L→∞ at the R73 rate — the convergence the pen's bridge quantifies (L=4 full-M available to
tighten it; not run here, per "re-weld constants first"). Reported as a deviation, not smoothed.

**Deviation 2 (q=5,7 per-scale ratio 3/q not cleanly resolved at L=2).** The pre-registered c_{k+1}/c_k → 3/q is
**contaminated at L=2**: only k≤2 is exact, and c_k drops into the O(1e-3) band-oscillation noise floor by k≈8
(c_k even changes sign — the L=2 band ringing around the saturated X_∞). The **saturation itself is unambiguous**
(X_k flat to 3 digits); the precise 3/q rate would need larger L. The clean signal is saturation, not the ratio.

## R1-C — mass sequences m_n = ‖M_tower^n · uniform‖₁, n=1..200 — **RAW DUMP (shape sealed for the pen)**
Delivered raw to `outputs/thread3_mass_sequences.tsv` (4 sequences × 200 steps). **No exponent extraction, no
regression, no log-log** — per the guard. What the raw sequences show at the surface:
- **q3_L3_tower / q3_L4_tower / q3_L3_full:** dominated by the partner geometric ρ_L ≈ 1/3 (successive ratio
  ≈ 3.0). Wilson's Prediction R (ringing period 3^{L−1} = 9 at L=3 / 27 at L=4, damping lifetime ≈ 8 / 74 steps)
  lives in the **subleading residual** m_n − c·ρ_L^n — present as a small oscillation in the successive ratios
  (band |λ|≈0.30 < partner 0.3335, so it decays ~0.90^n ≈ one 9-step cycle at L=3), **left for the pen to extract**.
- **q7_L2_tower (control):** decays at ≈ 1/6.5 per step (partner ≈ 0.154, far below q=3's 1/3) with **no slow
  ringing** — fast-damped, gapped band, exactly the control's committed signature.
- **q3_L3_full** (kinematic + tower) provided for the two-sector interference version the bridge will need.

**Deviation 3 (control at L=2, not L=3).** q=7 L=3 = 7.4M states — build infeasible on CPU. The control ran at
L=2; the qualitative signature (gapped, fast-damped, no slow ringing) is L-robust. Noted, not hidden.

## Status
Frozen-instrument spine laid: **exact constants welded (2/3, 10/21, 3·7/45=7/15)**; boundary contrast confirmed
(q=3 accumulates linearly, q≥5 saturate); mass sequences dumped raw for the pen's sealed Prediction-R shape
judgment. Three deviations reported as deviations (q=3 slope = L-truncation −2.5%; q=5,7 ratio contaminated at
L=2; control at L=2). The pen's session-two bridge (7/15 flat level ⟷ Tao's per-scale O(1) loss, with 7/45 as
the Plancherel price tag) rides on this spine and the frozen 7/45 definition transcribed above. No fitting; raw
dumps in the TSVs.
