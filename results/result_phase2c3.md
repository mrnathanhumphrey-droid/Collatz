# Result — PROBE 2c3 (A+B): the rung-2 JUDGE PASSES. At β*=3/5 the v₃=0 ratio is a well-defined function of (a mod9, γ mod9, e mod6) taking EXACTLY three values {12/49, 17/49, 20/49} in a balanced 72/72/72 split (both L) — the pre-registered SHAPE. τ exists at mod-9 resolution and is genuinely 3-variable (it does NOT reduce to a+γ mod9 or any coarser key). Mod-9 gate/carry tables + R₀(s) mass weights delivered (A). The v≥1 residual is a separate 5-value structure on (a mod9, γ mod27, e mod6), spread 1/7 — its own dressing. C held for Wilson's τ.

**Date:** 2026-07-16. Rung-2 corrector inputs (A) + level-set judge (B); C (2-parameter line search) held until Wilson sends τ. Direct/exact, no eigen-solves. Probe `probes/probe_phase2c3.py`, log `logs/probe_phase2c3_log.txt`, dumps `outputs/rung2_{gate_mod9,judge_v0_q3_L{2,3}}.tsv`. Claude gates; τ is Wilson's to derive blind.

**Headline: the pre-registered SHAPE holds exactly. The Collatz–Wielandt ratio r at β*=3/5, restricted to v₃=0 tower states, is WELL-DEFINED on the mod-9 key (a mod9, γ mod9, e mod6) and takes EXACTLY three values {12/49, 17/49, 20/49}, split 72/72/72 across the 216 keys — identical at L=2 and L=3. It does not factor through any coarser key I tested ((a+γ) mod9, (a·γ) mod9, (a mod9, γ mod9) without e) — τ is a genuine three-variable label. The level-set table is dumped for Wilson's τ to be judged entry-by-entry. A delivers the mod-9 gate/carry tables (all six units mod 9 can pass for every e′ mod 6; the structure lives in the (u,γ) mod9 → γ′ mod3 carry map) and the exact R₀(s) mass weights. The v≥1 residual is a distinct, richer object — 5 values on (a mod9, γ mod27, e mod6), spread 1/7 — confirming Wilson's read that rung two is a pair of trits, one per v-class.**

## B — THE JUDGE: v₀ ratio level sets at β*=3/5
| L | key | well-defined | # values | values | split |
|---|---|---|---|---|---|
| 2 | (a mod9, γ mod9, e mod6) | ✅ True | **3** | {12/49, 17/49, 20/49} | 72 / 72 / 72 |
| 3 | (a mod9, γ mod9, e mod6) | ✅ True | **3** | {12/49, 17/49, 20/49} | 72 / 72 / 72 |

- **Pre-registered SHAPE MET: exactly three populations on v₀, r constant on the level sets of (a mod9, γ mod9, e mod6).** τ exists at mod-9 resolution; the three level sets are the balanced 72/72/72 partition of the 216-key domain (6 units_a × 6 units_γ × 6 phases_e). Values {12/49, 17/49, 20/49} are the 2c2 residual ladder (mean exactly 1/3).
- **τ is genuinely 3-variable — it does NOT reduce.** r does not factor through (a+γ) mod9, (a+γ) mod9 with e, (a·γ) mod9 with e, or (a mod9, γ mod9) without e — all fail well-definedness. All three coordinates (a, γ to mod 9; e to mod 6) are load-bearing and combine non-trivially. (Structural glimpse from the dump: for fixed (a,γ) the value cycles through the three as e mod6 varies, but with a γ-dependent grouping — e.g. (a,γ)=(1,1): e∈{0,1}→12/49, {2,3}→17/49, {4,5}→20/49; (1,2) shifts the grouping. This is the τ signature Wilson must reproduce.)
- **Level-set table dumped** `outputs/rung2_judge_v0_q3_L{2,3}.tsv` (216 keys → exact fraction) — the adjudicator for Wilson's derived τ.

## A — mod-9 gate/carry tables + mass weights (rung-2 inputs)
- **Gate/carry mod 9 (v₀ sources), dump `outputs/rung2_gate_mod9_q3.tsv`:** for each e′ mod 6 (which sets c = 1−2^{e′} mod 9 and the LTE class), the (u mod9, γ mod9) → (gate pass, γ′ mod 3) map. **All six units mod 9 {1,2,4,5,7,8} can pass for every e′ mod 6** (the gate is mod 3; for any u some γ passes) — the discriminating structure is entirely in the carry γ′ mod 3 as a function of (u, γ) mod 9. c by class: e′≡0→c=0 (t=2, the digit-shift), e′≡1,3,5→c=8,2,5 (t=0, odd), e′≡2,4→c=6,3 (t=1).
- **R₀(s) mass weights (exact):** L=2: {65/189, 34/189, 20/189, 16/189, 20/189, 34/189} (s=0..5, symmetric s↔D−s); L=3: {262145, 131074, 65540, 32776, 16400, 8224}/786429. These weight the shift-s contributions to the ratio — the move-algebra half of τ.

## The v≥1 residual (rung-2, one level up)
- The v₃≥1 ratio at β*=3/5 is **well-defined on (a mod9, γ mod27, e mod6)** with **5 values {2/7, 11/35, 5/14, 2/5, 3/7}** (spread 3/7−2/7 = **1/7**, matching Wilson's flagged 1/7), both L. Richer than v₀'s 3-value ladder and needing γ to mod 27 (one 3-adic level deeper — the v≥1 carry sits deeper). **Confirms rung two is a pair of trits, one per v-class**: v₀ gets a 3-value dressing on mod-9, v≥1 a finer one on mod-27.

## Adjudication
| part | verdict |
|---|---|
| B judge (v₀) | ✅ SHAPE MET — 3 values {12/49,17/49,20/49}, well-defined on (a,γ mod9, e mod6), 72/72/72, both L; τ genuine-3-variable, dumped. |
| A gate tables | delivered (mod-9 gate/carry + R₀(s) exact) — the rung-2 τ inputs. |
| v≥1 residual | 5 values on (a mod9, γ mod27, e mod6), spread 1/7 — separate dressing (pair-of-trits confirmed). |
| C | HELD for Wilson's τ (2-parameter line search β=3/5 frozen, β₂ on τ). |

**⟹ The judge passes on the nose: the bad-cell residual at β*=3/5 is an exact three-valued function of the mod-9 data, balanced and L-invariant — the SHAPE Wilson pre-registered. The level sets are dumped as the blind-derivation adjudicator, and the mod-9 gate/carry + R₀(s) tables are handed over as τ's raw material. τ is genuinely three-variable (no coarse reduction), so the derivation is a real combinatorial object, not a relabeling. The v≥1 sector is a parallel, deeper (mod-27) 5-value residual — the second trit of the rung. C fires once τ lands.**

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, H, J, L-A, L-B, Real-T1, R1-STOP, Probe P/C/E/W/F2/G0/G1/G2, 2c(0+1), 2c0, 2c1, 2c2. No `r_q` value changes; no rate-law fit (r is contained in the bracket; nothing fitted; τ is Wilson's derivation, the level sets only judge it). β*=3/5 and the ladder {12/49,17/49,20/49} carry over from 2c2 unchanged.

_Reporting discipline: the judge is reported as SHAPE-MET (3 populations, well-defined on the pre-registered key), with the non-reduction to coarser keys stated so τ is not understated as a simple sum. The structural glimpse (e-cycling with γ-dependent grouping) is offered as an OBSERVATION for Wilson's derivation to match, not as a derived τ (his to build blind). The v≥1 residual is reported at its actual resolution (mod-27, 5 values), not forced into the v₀ mold. All values exact rationals. C is explicitly held._
