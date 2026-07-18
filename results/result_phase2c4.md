# Result — PROBE 2c4: the destination-τ matrix is τ-uniform by row (L=3 exactly [1/3,1/3,1/3]) — no source-τ correlation (A). The g_diag = (−13,+2,+11)/49 DIRECTION is confirmed (the v₀-flattening g is parallel, cos 0.9993) but the single v₀ trit cannot fully flatten v₀ (floor 0.083) and only shrinks the bracket 1.14× — v≥1 is now binding and gets perturbed, so the PAIR of trits (v₀+v₁) is needed (B). The v₁ mod-9 gate tables are delivered (C).

**Date:** 2026-07-17. Rung-2 dressing (A destination-τ matrix, B the search) + the v₁ front tables (C). Direct/exact + a 2-parameter search, no eigen-solves. Probe `probes/probe_phase2c4.py`, log `logs/probe_phase2c4_log.txt`, dump `outputs/v1_gate_mod9_q3.tsv`. τ = r-rung at β*=3/5 (τ=0↔12/49, 1↔17/49, 2↔20/49).

**Headline: (A) the destination-τ matrix is τ-uniform by row at both L — the destination τ-distribution is independent of source τ (no correlation); at L=3 it is the FULLY uniform [1/3, 1/3, 1/3], at L=2 it is [1/6, 1/2, 1/3] (boundary). (B) The reference direction g_diag = (−13,+2,+11)/49 is CONFIRMED — the v₀-spread-minimizing g is parallel to it (cos 0.9993, both L) — but the single v₀ trit CANNOT fully flatten v₀ (floor spread 0.083, not 0) because of the destination back-reaction and the undressed DEEP branch's dilution; and the full Collatz–Wielandt bracket only shrinks 1.14× (9/49 → 0.1605) because v≥1 is now the binding constraint AND gets perturbed by the v₀ dressing. The residual carrier is v≥1 — so the v₀ trit alone is nearly spent, and the PAIR of trits (v₀ + v₁, dressed jointly) is the real lever. (C) The v₁ mod-9 gate/carry tables (even-e′ channels, LTE t-classes) are delivered for the blind v₁ derivation.**

## A — the destination-τ matrix (rung-3 raw material + diagonal-approx adjudicator)
For each source τ-level on v₀, the mass-distribution of destination τ among its v′=0 targets (3×3, exact):
| L | src τ=0 (12/49) | src τ=1 (17/49) | src τ=2 (20/49) | verdict |
|---|---|---|---|---|
| 3 | [1/3, 1/3, 1/3] | [1/3, 1/3, 1/3] | [1/3, 1/3, 1/3] | **τ-uniform by row, fully uniform** |
| 2 | [1/6, 1/2, 1/3] | [1/6, 1/2, 1/3] | [1/6, 1/2, 1/3] | **τ-uniform by row** (dist [1/6,1/2,1/3]) |

- **All three source rows are IDENTICAL at both L — no source-τ → destination-τ correlation.** At L=3 the row is the fully uniform [1/3,1/3,1/3]; at L=2 it is [1/6,1/2,1/3] (a small-system boundary form). This is the clean rung-3 transfer material and the adjudicator for the diagonal approximation.

## B — the search: g_diag direction confirmed, magnitude/flatten limited
Corrector `h₂ = h_β·(1 + g(τ))` on v₀, g mean-zero. Reference g_diag = (−13,+2,+11)/49 (= the source-factor that would send r/(1+g_τ) → 1/3 if destinations didn't react).
| readout | L=2 | L=3 |
|---|---|---|
| g_diag: v₀ residual spread | 0.1195 (NOT flat) | 0.1208 (NOT flat) |
| **v₀-spread-min g, direction cos vs g_diag** | 0.9990 | **0.9993** |
| v₀-spread-min: g values | (−0.181,+0.036,+0.145) | (−0.184,+0.035,+0.149) |
| v₀ floor spread (best achievable by the trit) | 0.0804 | 0.0834 |
| full-bracket min g* | (−0.100,−0.013,+0.113) | (−0.109, 0, +0.109) |
| full bracket width at g* | 0.1605 | 0.1605 (rung-1 was 9/49 = 0.1837; **shrink 1.14×**) |
| residual carrier | v≥1 (O,v1) | v≥1 (O,v1) |

- **g_diag DIRECTION confirmed:** the g that best flattens v₀ is parallel to (−13,+2,+11) to cos 0.9993 — Wilson's diagonal optimum is the correct leading direction (consistent with A being τ-uniform, no correlation).
- **But the single v₀ trit cannot fully flatten v₀:** the best achievable v₀ spread is 0.083 (not 0), and g_diag itself leaves 0.12. The reason is the **destination back-reaction** — dressing v₀ changes the h at the UP/DOWN destinations (which are themselves v₀, weighted by h_β), so the correction is not the clean `r → r/(1+g_τ)`; plus the **undressed DEEP branch** (v′≥1 targets) dilutes. So the diagonal approx captures the direction but overshoots the magnitude (optimal is ~0.7× g_diag).
- **The bracket shrinks only 1.14× (9/49 → 0.1605)** because **v≥1 is now the binding constraint** and, worse, the v₀ dressing PERTURBS v≥1 (its spread grows from ~1/7 to ~0.16). So the full-bracket optimum uses a small g (spares v≥1) rather than the larger v₀-flattening g. **The v₀ trit alone is nearly spent.**
- **⟹ the PAIR of trits is required:** v₀ and v₁ must be dressed JOINTLY, because the v₀-only dressing trades against v≥1. This is exactly Wilson's "rung two = a pair of trits" — and C provides the v₁ half.

## C — the v₁ analog tables (the pen's next blind target)
mod-9 gate/carry for v≥1 sources (γ≡0 mod 3), even-e′ channels, dump `outputs/v1_gate_mod9_q3.tsv`:
- **e′≡0 mod 6** (c=0, LTE t=2 — the digit-shift channel), **e′≡2** (c=6, t=1), **e′≡4** (c=3, t=1). All six units mod 9 {1,2,4,5,7,8} pass (the discriminating structure is in the (u,γ) mod9 → γ′ mod3 carry map, tabulated in the dump).
- R₀(s) mass weights (L=3): {262145, 131074, 65540, 32776, 16400, 8224}/786429.
- Same protocol as v₀: tables first (here), Wilson derives the v₁ trit blind, judge after. The v₁ 1/7 residual needs its own trit; the two are then dressed together.

## Adjudication
| part | verdict |
|---|---|
| A destination-τ | ✅ τ-uniform by row (L=3 [1/3,1/3,1/3] fully uniform, L=2 [1/6,1/2,1/3]); no source-τ correlation. |
| B g_diag direction | ✅ confirmed (v₀-flatten g ∥ g_diag, cos 0.9993). |
| B single-trit flatten | ⚠️ LIMITED — v₀ floor 0.083 (back-reaction + DEEP dilution); bracket shrinks only 1.14×; v≥1 binding + perturbed ⟹ pair-of-trits needed. |
| C v₁ tables | delivered (even-e′ mod-9 gate/carry + R₀(s)). |

**⟹ The destination-τ matrix is clean and τ-uniform (A), and the g_diag direction is confirmed — Wilson's diagonal optimum is the right leading correction. But the honest finding is that the v₀ trit alone is nearly exhausted: it cannot fully flatten v₀ (destination back-reaction) and, more decisively, v≥1 is now the binding residual and gets perturbed by the v₀ dressing, so the single-trit bracket shrinks only 1.14×. The rung is a PAIR of trits (v₀+v₁) that must be dressed jointly, and C hands over the v₁ material for the blind derivation. Next: Wilson derives the v₁ trit; the joint 2-trit search is the real rung-2 shrink.**

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, H, J, L-A, L-B, Real-T1, R1-STOP, Probe P/C/E/W/F2/G0/G1/G2, 2c(0+1), 2c0, 2c1, 2c2, 2c3(A+B), 2c3-gate. No `r_q` value changes; no rate-law fit. τ, β*=3/5, and the ladder {12,17,20}/49 carry over from 2c2/2c3 unchanged.

_Reporting discipline: A is reported exactly (τ-uniform, with the L=2 vs L=3 distribution difference stated). B is an HONEST MIXED result — the g_diag DIRECTION is confirmed (cos 0.9993) but the "g_diag flattens v₀" reading is REFUTED (spread 0.12, not 0), and the reason is named (destination back-reaction + DEEP dilution), not smoothed; the modest 1.14× bracket shrink and the v≥1-binding/perturbation are reported as the reason the v₀ trit alone is nearly spent. The corrector form used (h_β·(1+g(τ))) is stated explicitly so Wilson can reconcile if his intended form differs. The pair-of-trits conclusion is drawn from the data (v≥1 binding), not assumed. C tables are exact._
