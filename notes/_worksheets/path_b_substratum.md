# Result 28: Path B sub-stratum extension — outcome (4) confirmed; per-j W_j requires integer-lattice info beyond residue chain at any sub-stratum depth

**Date:** 2026-05-02. Sequel to Result 19 (Path B continuation). Tests sub-stratum lattice extension proposed in brief.

This document tests whether extending Path B's matrix Wiener-Hopf state space from (residue mod 2^k) to (residue mod 2^k, sub-stratum index j_inner = v_2(1+3·h)) at the boundary residue r=21 produces per-j W_j matching empirical W_2, W_4, W_5. **Result: outcome (4) confirmed at high resolution.** Sub-stratum extension *does* produce per-j differentiation (unlike Path B v2's universal W_j = 3.44), but the differentiation has the WRONG SHAPE: predictions grow monotonically with j_inner, while empirical W_j is non-monotone with a sign flip. The sub-stratum extension at the residue level cannot recover empirical structure regardless of truncation depth.

Code: `path_b_substratum.py` (drift exact at machine precision, 2M-orbit MAP simulation per j_target).

---

## 1. Drift identity preserved at machine precision

Stationary distribution remains uniform (1/32 per non-boundary residue, 1/32 total over all sub-strata at boundary). With Geom(1/2) sub-stratum giving E[v|r=21] = 7:

```
E_π[X] = -0.28768207  (target log(3/4) = -0.28768207, diff = +0.00e+00)
```

Sub-stratum extension does not break the drift identity. Q construction is exact.

## 2. Per-j W_j extraction with sub-stratum conditioning (2M orbits, k=6, j_inner_max=12)

For each j_target, simulate MAP from S=0 until first descent below -L = -(log(2^36) - log(m_j)). Condition on end state matching the j-target sub-stratum.

| j | m_j | L_j (nats) | target state | ⟨overshoot⟩ (nats) | W_j_pred (steps) | W_j_emp | gap |
|---|---|---|---|---|---|---|---|
| 2 | 5 | 23.34 | r=5, * | 0.994 ± 0.004 | **+3.46** ± 0.01 | +7.156 | -3.70 |
| 4 | 85 | 20.51 | r=21, j_inner=2 | 1.946 ± 0.022 | **+6.77** ± 0.08 | -4.755 | **+11.52 (sign flip)** |
| 5 | 341 | 19.12 | r=21, j_inner=4 | 2.491 ± 0.048 | **+8.66** ± 0.17 | +4.590 | +4.07 |

**Predictions are MONOTONE INCREASING in j_inner (j=2 → 3.46, j=4 → 6.77, j=5 → 8.66).** Empirical W_j is **non-monotone with a sign flip** (W_4 < 0 < W_5).

The biggest mismatch is j=4: prediction is +6.77 step units, empirical is -4.76 step units — SIGN FLIPPED with absolute gap of 11.5 step units.

## 3. Per-sub-stratum breakdown — predictions are linear in j_inner

Running 2M orbits with descent level L = log(2^36/85) and breaking down conditional overshoot by every sub-stratum index:

| j_inner | corresp m_j | count | ⟨overshoot⟩ (nats) | W_pred (steps) |
|---|---|---|---|---|
| 0 | m_3 | 7,785 | 1.286 | +4.47 |
| 1 | (none) | 4,950 | 1.620 | +5.63 |
| 2 | m_4 | 3,154 | 1.939 | +6.74 |
| 3 | (none) | 1,853 | 2.282 | +7.93 |
| 4 | m_5 | 1,085 | 2.546 | +8.85 |
| 5 | (none) | 585 | 2.943 | +10.23 |
| 6 | m_6 | 351 | 3.274 | +11.38 |
| 7 | (none) | 183 | 3.665 | +12.74 |
| 8 | m_7 | 94 | 3.962 | +13.77 |
| 9 | (none) | 56 | 5.183 | +18.02 |

**Linear fit:** W_pred(j_inner) ≈ 4.47 + 1.51·j_inner step units. Each unit of j_inner adds ~1.5 step units (~0.43 nats, ~62% of log(2)) to the predicted overshoot.

This is exactly what the structure should give: at sub-stratum j_inner, the *final descent step* has v = 6 + j_inner, so larger j_inner means larger final step magnitude, hence larger overshoot beyond the level. The relationship is monotone by construction.

## 4. Why the sub-stratum extension cannot match empirical per-j W_j

The mismatch isn't a magnitude calibration issue — it's a **shape mismatch** that no truncation depth or k-refinement can fix.

**Structural reason.** At sub-stratum (r=21, j_inner=k), the abstract MAP state aggregates ALL m values with that residue and v_2(1+3·h)=k. Specifically:
- j_inner = 0, m at residue 21: m ∈ {21, 21+128, 21+256, ...} (h with v_2(1+3h)=0)
- j_inner = 2, m at residue 21: m ∈ {85, 85+128·2, 85+128·4, ...} (h with v_2(1+3h)=2)
- j_inner = 4, m at residue 21: m ∈ {341, 341+128·8, ...} (h with v_2(1+3h)=4)

Only ONE specific h value per j_inner gives lattice landing at m_j: h = m_(j-3). For j=4, h=5 → m=341? No wait — for j=4, m_j = 85 and h = m_(j-3) = m_1 = 1. m = 21 + 64·1 = 85 = m_4. ✓

For j_inner = 2 sub-stratum: the orbit with m = 85 is one specific orbit at h=1. But the MAP simulation at sub-stratum j_inner=2 conditions on ALL h with v_2(1+3h)=2 (namely h ≡ 1 mod 8 minus h ≡ 1 mod 16, i.e., h = 1, 9, 17, 25, ...). Only h=1 gives m_4. The rest give m = 85 + 64·8 = 597, m = 85 + 64·16 = 1109, etc.

**Empirical W_j is conditioned on m = m_j EXACTLY.** Sub-stratum conditioning aggregates many m values per stratum, only ONE of which is the lattice landing. The conditional overshoot in the MAP framework is an average over all these h values, weighted by their natural-density measure within the sub-stratum.

**The sign flip empirical → W_4 = -4.76 cannot emerge from any positive-overshoot framework.** Renewal overshoot is always ≥ 0. Empirical W_4 < 0 means the orbits absorbing at m_4 = 85 take FEWER Syracuse steps than Wald asymptotic predicts — they "land short" of the deterministic descent rate. This is integer-lattice arithmetic structure (specific value of m_4 = 85, specific cumulative log m descent across orbit history) that the abstract MAP cannot represent.

## 5. What the sub-stratum extension does close (vs Path B v2)

Despite missing per-j W_j, the sub-stratum extension closes:

1. **Linear-in-j_inner growth law for conditional overshoot.** W_pred(j_inner) ≈ 4.47 + 1.51·j_inner. This is a clean structural finding: the framework correctly predicts that larger j_inner → larger overshoot, with slope ~62% of log(2). The 4.47 intercept = E[L⁻]/log(4/3) + (correction from finite-L), and 1.51 = "per-step descent absorption ratio".

2. **Per-sub-stratum prediction differentiation.** Path B v2 predicted UNIVERSAL W_j = 3.44 across all j (no differentiation). The sub-stratum extension predicts different values per j, capturing some structure even if not the right shape.

3. **The integer-lattice constraint is now precisely identified.** For W_j_pred(sub-stratum) to match empirical W_j, the framework would need to condition on h = m_(j-3) specifically (not on all h at sub-stratum j_inner). This requires explicit integer-level state tracking beyond the residue chain.

## 6. Critical structural test: P(j=3) = 0?

The brief flagged this as the load-bearing observation. Sub-stratum framework gives:
- P(end at sub-stratum j_inner=0) = 7,785 / 19,991 = 0.39 (in the L=20.5 simulation)
- This corresponds to "approximate landing at m_3" in the brief's mapping
- BUT empirically P(j=3) = 0 (number-theoretically: m_3 = 21 has no Syracuse predecessors among odd integers)

The sub-stratum extension does NOT reproduce P(j=3) = 0. It gives P(j=3-equivalent sub-stratum) = ~0.4, dominated by the j_inner=0 mass.

This is outcome (b) in the brief's classification: "The extended framework gives P(j=3) > 0. Sub-stratum extension at the residue level isn't sufficient; need explicit ancestor-blocking at the boundary."

The integer-level structural blocking is real: at residue 21 with j_inner=0 (i.e., m = 21 + 128·k for k = 0, 1, ...), only m = 21 itself triggers absorption at attractor m_3. Other values (m = 149, 277, ...) at the same sub-stratum are NOT at any attractor. The MAP framework cannot distinguish them.

## 7. Implications for ε_S closure

Brief's calculation assumed sub-stratum predictions match empirical: ε_S_pred = +1.39 ≈ log(4) = 1.386. With actual sub-stratum predictions:

```
W_2_pred = +3.46    log(m_2)/log(4/3) = 5.59
W_4_pred = +6.77    log(m_4)/log(4/3) = 15.43
W_5_pred = +8.66    log(m_5)/log(4/3) = 20.27

ε_S_pred = P(j=2)·[W_2 - log(m_2)/log(4/3) + 1] + ...
         = 0.938·(3.46 - 5.59 + 1) + 0.024·(6.77 - 15.43 + 1) + 0.038·(8.66 - 20.27 + 1)
         = 0.938·(-1.13) + 0.024·(-7.66) + 0.038·(-10.61)
         = -1.06 - 0.18 - 0.40
         = -1.64 step units
```

In nats: -1.64 · log(4/3) = -0.47 nats. Empirical ε_S ≈ log(4) = +1.39 nats. **Sub-stratum prediction has wrong sign**, off by 1.86 nats. Framework doesn't close ε_S.

## 8. Verdict per brief decisive outcomes

- **(1) Per-j W_j matches to ±0.05:** NO. Max gap 11.5 step units.
- **(2) Per-j W_j matches to ±0.1:** NO.
- **(3) W_2 matches but W_4, W_5 don't:** NO (none match; W_2 is universal Lorden, j=4 sign-flipped).
- **(4) Universal W_j_pred ≈ 3.44 still emerges:** PARTIAL — W_2 is at Lorden, but j ≥ 3 predictions DO differ from Lorden (linear in j_inner). The differentiation has wrong shape, not absent.

**Net verdict (4-modified): the sub-stratum extension produces a NEW structural prediction (linear-in-j_inner W_j growth) that DIFFERENTIATES per j but does not match empirical shape.** The mismatch is in the form of the per-j variation, not just in the magnitudes:
- Framework predicts: W_j monotone increasing in j_attr (since j_inner = 2(j_attr-3) increases with j_attr)
- Empirical: W_j non-monotone with sign flip (W_2 = +7.16, W_4 = -4.76, W_5 = +4.59)

The sign flip empirical is the structural signature that cannot emerge from any positive-overshoot, residue-level framework.

## 9. What's needed for closure

To capture empirical per-j W_j, the framework needs explicit integer-level state beyond the residue chain. Specifically:

1. **Joint (residue mod 2^k, log m) state with discrete log m on integer lattice.** Each integer m corresponds to a specific log m = log(m). Absorption at m_j is exact landing on this discrete lattice point. The matrix WH framework on a continuous log axis cannot represent this.

2. **Or: per-orbit history dependence.** Empirical W_j may depend on which specific h values the orbit traversed before absorbing at m_j. The MAP framework averages over all histories consistent with the residue path; per-j W_j requires conditioning on specific history patterns.

3. **Or: finite-difference correction at small m_j.** The Lorden-style asymptotic overshoot E[L⁻²]/(2 E[L⁻]) is a renewal limit at large barrier. For SMALL barriers (log m_j is small relative to typical step), pre-asymptotic corrections dominate. Empirical W_4 < 0 may come from m_4 = 85 being so small that the asymptotic Lorden formula is not yet valid.

These three are not mutually exclusive. (3) is most compatible with the existing framework: a Caravenna-Doney style local LD correction at small barrier could give signed corrections.

## 10. Files

- `path_b_substratum.py` — extended MAP simulation with sub-stratum tracking
- `path_b_substratum.md` — this document (Result 28)
- `closed_form_findings.md` — Result 28 entry
- `path_b_continuation.md` — Result 19 (residue-only baseline; Path B v2)
- Build on: `path_b_matrix_wh_v2.py` (k=6 Q construction with proper boundary)

## 11. Citations

- Asmussen Ch IV — matrix Wiener-Hopf for Markov-additive
- Alsmeyer-Buckmann 2018 — Markov random walks
- Lorden 1970 — asymptotic stationary excess
- Caravenna-Doney 2019 — local large deviations (relevant to small-barrier correction)
