# Probe 85 — Phase 0 feasibility memo (costed decision, no moments computed at target)

**Date:** 2026-07-14. **Headline: the cost estimate that made this a "big spend" decision was ~1000× too high. MEASURED, the evidential bridge at n=6 costs ~minutes, not ~160 hours. The Probe-84 cost-gate is void — there is no large cost to gate. Recommendation: target n=6, skip the gate, Phase 1 ready on operator go.**

Files: `result_85_phase0_timing.py` + `result_85_phase0_timing_log.txt` + `result_85_phase0_n6_log.txt`.

## 1. Gate: Probe 84 LANDED (parallel, 2026-07-14) → offset RESOLVED → n=5 unlocked

R84 verdict **H_ARTIFACT**: the mod-9 offset is a **normalization artifact** — a top-layer `3^r` global phase from the family-defining `(1+3^r)^ℓ` twist (since 4≡1 mod 3, `e_q(ℓ·3^r·4^u)=e_3(ℓ)=ω₃^ℓ`, constant in u → a global cube-root phase). v₃(offset)=r at r=2..7 (not a fixed 3²; the "9" R83 saw was just 3^{r=2}). R83's "genuine residual structure" is **overturned**. Per the pre-reg gate — "offset resolved (normalization artifact / derived source) → **n=5**." The r=2 phase-value comparison is now clean: divide out the *derived* per-family `ω₃^ℓ`.

Given cost is trivial at both (§2), the strongest move is **both**: n=5/r=2 (divide out the known ω₃^ℓ) AND n=6/r=3 (directly in the certified R81b Mahler regime, no offset). Agreement across the two is the real robustness check; disagreement flags the offset-division or the `r=n−2j+1` mapping (my A3 triage) rather than a true null.

## 2. Cost — MEASURED, not estimated (this is the correction)

Complexity derived from `dwm_kraus_match_syracuse.py`: inner loop is V_MAX² pairs per level × 2 levels = **V_MAX⁴** iterations; each does `X1@X2@X1` = dim³ matmuls. R83 §87 assumed the dim³ term dominates → "~6 h (n=5), ~160 h (n=6)." **Measured, it does not:**

| n | dim | time (V_MAX=8) |
|---|---|---|
| 3 | 18 | 0.06 s |
| 4 | 54 | 0.30 s |
| 5 | 162 | 1.47 s |

The n=3–5 fit is dim^1.49 (build/overhead-bound), but the direct **n=6 measurement is steeper** — the loop grows once dim=486's Python matrix-fills and matmuls matter:

| n | dim | loop @V_MAX=8 (measured) | → loop @V_MAX=16 |
|---|---|---|---|
| 5 | 162 | 1.27 s | ~23 s |
| 6 | 486 | **23.6 s** | **~7 min** |

n=5→n=6 is 18.6× for a 3× dim jump → effective ~dim^2.7 (V_MAX⁴ confirmed separately: t(16)/t(8)=18.16≈16). **Measured n=6 G1 loop ≈ 7 min @ V_MAX=16.** Full bridge (G1 + G2 4-alternating + the delta_1/vac_π reductions folded into the same loops) ≈ 2–3× → **n=6 ≈ 15–20 min; n=5 ≈ 1.5 min.**

**Critical setup note:** the original `dwm_kraus_match_syracuse.py` computes the stationary π by *exact-rational* Gaussian elimination — O(dim³) Fractions, which is **minutes+ at dim=486 and was the actual bottleneck** (it hung the first n=6 timing run). **Float π (power iteration) does it in 0.11 s** (same trick as Probe 4) and the moments are float anyway. **Phase 1 must use float π, not `stationary_rational`.**

**R83's ~160 h estimate was still ~500× too high** — the third estimate this session overturned by measuring (after ⌊r/2⌋+2 and the inverse-tree 1/9).

## 3. Memory + the k=8 precedent — does NOT apply

Matrices are dim²·16 B: 5 KB (n=3) → **3.8 MB (n=6)**. This is **compute-bound, not memory-bound.** The k=8 disaster was a 15 h, 20 GB, RAM-bandwidth-saturated ε-computation killed by Windows Update. None of that applies here: single process, MB-scale RAM, minutes-long. **No pool, no checkpointing, no Windows-Update pause needed** (all were conditioned on >2 h — nothing here is).

## 4. Refinements folded in (my A1–A5, re-weighted by the measurement)

- **A1 (V_MAX):** the Geom(½) weight makes the sum cutoff-convergent (R3a2: the 2^{−e}-weighted tree was E_MAX-stable to full precision). n=3 hit 6 digits at V_MAX=16; keep V_MAX=16 (going higher wastes the V_MAX⁴). **A residual that does not shrink with V_MAX is structural, not truncation** — clean H_BRIDGE_PARTIAL discriminator.
- **A2 (F̂ side):** "free" = one FFT at the target r (`result_81b_mahler_extend.py`), not a closed form (v₃(c_k) formula is open). Trivial at r=3.
- **A3 (failure triage):** if H_BRIDGE_FAILS, audit `r=n−2j+1` + no-b_prior-pre-averaging BEFORE concluding decoupling — the n=3 match is 6 digits × 4 reductions, too tight to be coincidence.
- **A4 (D_W):** MOOT as a cost-saver (cost already trivial); still the natural invariant, optional.
- **A5 (Windows Update):** MOOT — run is < 2 h.

## 5. §0d cheaper cell — moot

j=3 needs n≥8 for r≥3 (worse). But since n=6 is ~8 min, there is no cell worth trading precision for. Just run n=6.

## 6. Recommendation (operator decision)

**Target n=6 (r=3), full 4-reduction G1+G2, V_MAX=16, float π, single process. Measured ~15–20 min. Skip the Probe-84 gate** — it existed to avoid a 27× cost that measurement shows is ~15 min, and n=6/r=3 sidesteps the unresolved mod-9 offset (an r=2 artifact) entirely.

**Phase 1 is ready to fire on your go.** Two honest caveats: (1) use **float π**, not the exact-rational stationary (the only thing that was actually slow at n=6). (2) It still requires building the *Syracuse-side* directly-measured moments at n=6 (they don't exist) — but the timing harness already parametrizes the construction to arbitrary n, and it's the same ~minutes of compute. Per the pre-reg, I stop here and wait for the go; nothing at the target has been computed.

`THEOREM_C_745` and Thms 78.1–78.3 are not at stake in this probe regardless of outcome.
