# Probe B1 — THE BLINDNESS THEOREM (3x−1 spectral identity) — **CLAIM CONFIRMED; MECHANISM WALKED BACK**

**Date:** 2026-07-18  CPU, exact. Probe `probes/probe_blindness_B1.py`, log `logs/probe_blindness_B1_log.txt`.
M₋ built **directly** from the 3x−1 map (additive constant −1 ⟹ Syracuse variable π_k → −π_k ⟹ difference
coupling T=2^{−S}−2^{−S'} → −T), with **no use of Σ or of M₊** — the conjugacy is discovered, not baked.

## Verdict in one line
**The second-moment spectrum is sign-blind — spec(M₋) = spec(M₊) exactly — but the intertwiner is the
pair-swap P, not the committed negation Σ.** Theorem 7 (Blindness) stands; its lemma changes.

## B1-2 — SPECTRAL IDENTITY: **CONFIRMED to machine precision**
| quantity | M₊ | M₋ | \|diff\| |
|---|---|---|---|
| **partner (tower Perron), L=2** | 0.34682666 | 0.34682666 | 9.99e-16 |
| **partner, L=3** | 0.33323630 | 0.33323630 | 7.22e-16 |
| **doublet-top, L=3** | 0.23763996+0.18303042j | 0.23763996+0.18303042j | 1.03e-15 |
| **kinematic c_k family** (36 modes L=2 / 324 L=3) | — | — | **0.00e+00 exact** |
| full spectrum L=2 (324, nearest-neighbor) | — | — | 7.21e-9 |
| full spectrum L=3 (8424) | — | — | **exact by permutation similarity** |

The partner values reproduce Wilson's pre-registered targets **0.34682666 / 0.33323630 exactly**; the
doublet is the banked member #1. **No spectral difference exists at machine precision — the KILL condition
did not fire, the frozen construction is clean.** The partner, doublet, braid, crossing, and the entire
kinematic family are identical between the map with three nontrivial cycles and the map conjectured to have one.

## B1-1 — THE MECHANISM: negation FAILS, pair-swap is EXACT (walk-back #25)
Entry-by-entry, both L, exact rationals:

| intertwiner | test | max\|diff\| | verdict |
|---|---|---|---|
| **pair-swap** P:(a,b,γ)↦(b,a,γ) | M₋ = P M₊ P | **0.000e+00** | **PASS — this is the intertwiner** |
| negation Σ:(a,b,γ)↦(−a,−b,−γ) | M₋ = Σ M₊ Σ (Wilson) | 0.258 (L2) / 0.250 (L3) | **FAILS — not the intertwiner** |
| negation Σ | Σ M₊ Σ = M₊ (symmetry?) | 0.258 / 0.250 | **also fails — Σ is not even a symmetry** |

**Why Σ dies, and why P lives (the correction, for the pen's two-line lemma):**
- The pair operator sees only the **difference** T = ap − bp. Exchanging the two copies (P) sends T → −T with
  the carry untouched — a pure **integer relabeling** — so P M₊ P is exactly the −T operator = M₋. Clean.
- Negation Σ flips T → −T **and** γ → −γ and forces the carry γ′ = ⌊(γ+T)/q⌋ through modular negation. The
  carry-floor arithmetic does **not** commute with modular negation (⌊·⌋ vs −mod), so Σ fails by exactly the
  0.25–0.26 residual — the **same breakage that killed the J-involution (walk-back #14)**: "any invariance
  must act on the carry as an INTEGER map, not modular." Σ violates that; P (integer relabel) satisfies it.
- **The physics is unchanged:** T → −T *is* the 3x+1 ↔ 3x−1 additive-constant flip. The pen just realizes it as
  copy-exchange, not sign-negation. The blindness is if anything **cleaner** via P — an exact permutation
  similarity (0.000e+00), no arithmetic subtlety.

## B1-3 — CONTRAST COLUMN (documented; the operator provably cannot see this)
| | 3x+1 | 3x−1 |
|---|---|---|
| **spectrum** (B1-1/B1-2) | \| — identical to 1e-15 — \| |
| **positive cycles** | **1**: {1,4,2} | **3**: {1,2}, {5,14,7,20,10}, {17,50,25,74,37,110,55,164,82,41,122,61,182,91,272,136,68,34} |
| **inverse trees** | single basin | three basins (one per cycle) |

Cycle census computed directly (full 3n±1 map, starts ≤ 3000) — matches Wilson's committed list exactly.
The inverse-tree density asymmetry (banked 10³–10⁴×; three distinct 3x−1 basins vs one for 3x+1) is the
**one B1-3 sub-item cited, not re-run** — banked code path `probes/agent3_inverse_tree_3xm1_Dn.py`, records
`result_3a_inverse_tree_ninth.md` / `result_inverse_tree_residue.md`. Re-running the Fourier-tree infra for a
documented contrast column would be disproportionate; the crisp contrast (identical spectra vs 3-vs-1 cycles)
is self-contained above.

## For THEOREM 7 (the pen)
The operator through which every distributional approach factors assigns **identical spectra** to a map with
three nontrivial cycles and a map conjectured to have none — the distributional-to-pointwise barrier, located
and **proven** (exact permutation similarity M₋ = P M₊ P), not lamented. **Lemma correction:** the intertwiner
is the pair-exchange P, not the negation Σ (which fails on carry-floor/modular non-commutation, cf. walk-back
#14). Pre-registered magnitudes (partner 0.34682666 / 0.33323630) reproduced exactly; cycle list exact.
