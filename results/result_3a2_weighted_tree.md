# R3a2 — REFORMULATION: 2^{−e} edge-weighted inverse tree. Verdict: H_STABLE_NONCLEAN

**Date:** 2026-07-14. **The correct (cutoff-free) object is well-defined but has NO clean decay rate. 1/9 is dead in BOTH formulations. Sibling quantitative-rate thread CLOSED.**

Probe `result_3a2_weighted_tree.py`. R3a killed the uniform-count inverse measure (E_MAX-cutoff-dependent by orders of magnitude). This gives the idea its best shot: the natural Syracuse harmonic/pushforward measure — each inverse edge `y→g_-(y;e)` carries weight `2^{−e}` (the forward map's halving probability), cutoff-free by Σ_e 2^{−e} < ∞.

## 1. Validation + the methodological win

- **Validation PASS:** with unit edge weight the tower reproduces all **35/35** committed uniform D_n(k), n≤6.
- **Cutoff-free CONFIRMED.** Weighted D_n(k=4) across E_MAX ∈ {20,30,40}:

| n | E_MAX=20 | E_MAX=30 | E_MAX=40 |
|---|---|---|---|
| 6 | 2.708963 | 2.708919 | 2.708919 |
| 10 | 1.493132 | 1.493080 | 1.493080 |
| 14 | 0.431350 | 0.431330 | 0.431330 |

E30 = E40 to all shown digits (tail 2^{−40} ≈ 1e-12). **The weighting cures R3a's cutoff pathology — the object is genuinely well-defined.** This is the right fix, and it works.

## 2. But no clean rate

Weighted D_n(k=4), n=1..14: `30.5, 1.87, 10.4, 12.4, 10.0, 2.71, 3.97, 1.43, 2.00, 1.49, 0.822, 0.383, 0.460, 0.431` — **non-monotone, oscillating, not geometric.** Ratios D_{n+1}/D_n swing 0.27–1.47 (some > 1: D_n *increases* at steps), nowhere near 1/9 or any constant. Same wild picture at k=2,3,5. **H_STABLE_NONCLEAN.**

## 3. Verdict — the thread ends here, correctly

- **1/9 is dead in both formulations:** uniform-count (R3a: cutoff artifact) AND the correct 2^{−e}-weighted harmonic measure (R3a2: well-defined but no clean rate). This is the *strong* negative — the idea got its best shot (the right measure) and the inverse-tree Plancherel mass simply has no clean geometric decay for the 3x+1 single basin.
- **NOT a measure-choice artifact.** R3a2 rules out "we just used the wrong measure." The quantitative sibling *decay-rate* thread is **CLOSED**.
- **Untouched:** the qualitative sibling *paradox* (forward K₋=σK₊σ symmetric, inverse trees differ) needs no rate and stands. 3b (cycle-count detector) as a *qualitative* fingerprint could still be posed, but any *quantitative* version inherits the no-clean-rate finding — do not pursue a numeric rate.

## 4. Process

R3a2 is what a reformulation should do: take the falsification's diagnosis (wrong, cutoff-dependent object), build the principled fix (harmonic 2^{−e} measure), verify the fix addresses the diagnosis (E_MAX-stable — it does), and honestly report what the corrected object shows (no clean rate). The fix succeeded; the hoped-for result did not survive it. Files: `result_3a2_weighted_tree.py` + `result_3a2_log.txt`.
