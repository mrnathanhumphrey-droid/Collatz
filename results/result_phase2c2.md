# Result — PROBE 2c2: FIRST CORRECTOR WORKS. G1/G2 gate the survival dichotomy exactly ({1/9,4/9} with 4/9 ⟺ a+γ≡0 mod 3; (even,v₃=0)=2/9). The parity-sector dressing h_β with β* = 3/5 (exact) shrinks the bracket 0.4857 → 9/49 ≈ 0.1837 (2.64×), KILLS the (odd,v₃=0) survival dichotomy, and leaves a finer, parity-INDEPENDENT residual {12/49, 17/49, 20/49} on the two v₃=0 cells — the next corrector's target.

**Date:** 2026-07-16. First-order corrector (2c-2): gate the algebraic survival claims, then line-search the parity dressing. Direct/entry algebra + scalar line search, no eigen-solves. Direct at q=3. Probe `probes/probe_phase2c2.py`, log `logs/probe_phase2c2_log.txt`. Claude gates; ALGEBRAIC/STATISTICAL labeled. Structural pre-registration only (ledger 0-for-20; no magnitude pre-commit on β* or width).

**Headline: the two algebraic claims pass exactly — every (odd,v₃=0) tower state has survival ∈ {1/9, 4/9} with 4/9 ⟺ a+γ≡0 mod 3 (both L, machine precision), and every (even,v₃=0) state has survival exactly 2/9. The corrector h_β(x) = h⁰(cell)·(1 + β·σ·1_{(odd,v₀)}), σ=+1 ⟺ a+γ≡0 mod 3, has optimal β* = 3/5 EXACTLY: the Collatz–Wielandt bracket shrinks from width 0.485714 (17/35) to 9/49 ≈ 0.183673 — a 2.644× reduction — still containing 1/3 and the partner. Both structural pre-registrations are met: the width strictly shrinks, and the (odd,v₃=0) survival dichotomy is KILLED (σ=+1 and σ=−1 states now share the identical residual range). The residual is carried by BOTH v₃=0 cells with an identical, parity-INDEPENDENT 3-value structure {12/49, 17/49, 20/49} (mean exactly 1/3) — no longer the survival dichotomy, and the next corrector's target. L-invariant (identical L=2,3).**

## Gates — the survival dichotomy, exact
| gate | claim | L=2 | L=3 | verdict |
|---|---|---|---|---|
| **G1** | (odd,v₃=0) survival ∈ {1/9, 4/9}; 4/9 ⟺ a+γ≡0 mod 3 | ✅ (values {1/9,4/9}, dichotomy exact) | ✅ | **PASS** |
| **G2** | (even,v₃=0) survival = 2/9 exactly | ✅ (max\|·−2/9\| = 5.6e-17) | ✅ (8.3e-17) | **PASS** |

- **G1 confirmed ALGEBRAIC:** each (odd,v₃=0) state survives with mass exactly 1/9 or 4/9, and the 4/9 states are EXACTLY those with a+γ≡0 mod 3 (the σ=+1 set); exactly half the cell is at 4/9 (fraction 0.5000). This is the residue-forced dichotomy that made (odd,v₃=0) the sole statistical survival cell (2c-0a) — now resolved as a clean a+γ mod 3 split.
- **G2 confirmed ALGEBRAIC:** (even,v₃=0) survival is 2/9 with no per-state variation (the "cancellation" Wilson flagged) — the sign structure of the even sector cancels to a constant.

## Corrector — β* = 3/5, width 0.4857 → 9/49
`h_β(x) = h⁰(cell)·(1 + β·σ(x)·1_{(odd,v₀)}(x))`, σ = +1 iff a+γ≡0 mod 3, operator Mᵀ (flow convention).
| L | baseline β=0 | β* | bracket at β* | width at β* | shrink | 1/3 in | partner in |
|---|---|---|---|---|---|---|---|
| 2 | [0.118367, 0.604082], w=0.485714 | **+3/5** | [12/49, 3/7] = [0.244898, 0.428571] | **9/49 = 0.183673** | **2.644×** | ✅ | ✅ |
| 3 | same | **+3/5** | same | **9/49** | 2.644× | ✅ | ✅ |

- **β* = 3/5 exactly** (line-search minimum, recognized as 3/5 to machine resolution). `width(β)`: −0.5→1.195, −0.25→0.690, 0→0.486, +0.25→0.344, +0.5→0.223, **+0.6→0.184** — a clean single minimum at 3/5.
- **The bracket shrinks 2.644× and stays exact-rational:** [12/49, 21/49], width **9/49**. (Baseline was on denominator 245 = 5·49; β=3/5 clears the factor of 5.) The partner and 1/3 remain inside.
- **The (odd,v₃=0) survival dichotomy is KILLED:** at β*, the σ=−1 and σ=+1 sub-populations both span [12/49, 20/49] — they OVERLAP completely. The variance β was designed to absorb (the survival split) is gone, exactly as pre-registered.

## Residual — the next corrector's target
Per-cell ratio structure at β* (L=3):
| cell | range | spread | distinct ratios |
|---|---|---|---|
| (even, v₃=0) | [12/49, 20/49] | 8/49 | **{12/49, 17/49, 20/49}** |
| (odd, v₃=0) | [12/49, 20/49] | 8/49 | **{12/49, 17/49, 20/49}** |
| (even, v₃≥1) | [2/7, 3/7] | 5/35 | {2/7, 11/35, 2/5, 3/7} |
| (odd, v₃≥1) | [2/7, 3/7] | 5/35 | {2/7, 5/14, 3/7} |

- **The residual max-spread is carried by BOTH v₃=0 cells, with an IDENTICAL, parity-INDEPENDENT 3-value structure {12/49, 17/49, 20/49}** (mean = 49/3 / 49 = 1/3, still centered). The dressing collapsed the even/odd distinction in the v₀ sector: what was a parity-split survival dichotomy is now a single parity-blind 3-value ladder. **This is the next corrector's target** — a finer (deeper-digit) structure on the v₀ cells, no longer the survival dichotomy.
- The v₃≥1 cells carry a slightly smaller residual (spread 5/35 ≈ 0.143) with their own 3–4 value structure.
- **L-invariant:** β*, the width, and the residual {12/49,17/49,20/49} are identical at L=2 and L=3 — the first corrector operates entirely on the L-invariant cell/survival structure. The partner's actual L-flow toward 1/3 lives in still-finer structure that this order does not yet resolve (the bracket 9/49 is L-independent, so it contains the partner but does not track it — the higher-order / depth corrections must supply that, per the 2c0 selection rule × cascade tax).

## Adjudication
| item | verdict |
|---|---|
| G1 dichotomy | ✅ PASS — survival ∈ {1/9,4/9}, 4/9 ⟺ a+γ≡0 mod 3, both L (ALGEBRAIC). |
| G2 cancellation | ✅ PASS — (even,v₃=0) = 2/9 exact (ALGEBRAIC). |
| corrector shrink | ✅ β*=3/5, width 0.4857 → 9/49, 2.644× — strictly shrinks (pre-reg met). |
| dichotomy killed | ✅ σ=±1 overlap at β* — the survival variance is gone (pre-reg met). |
| residual carrier | v₀ cells, parity-blind {12/49,17/49,20/49} — NO longer the dichotomy; next target. |

**⟹ The first corrector works exactly as Wilson designed: a single rational parameter β*=3/5 dressing the bad cell shrinks the Collatz–Wielandt bracket 2.6× and eliminates the (odd,v₃=0) survival dichotomy, leaving a cleaner, parity-independent 3-value residual {12/49,17/49,20/49} on the v₃=0 cells. Both structural pre-registrations hold. The residual is the corrector chain's next rung — a finer-digit structure, and (being L-invariant) still a shallow one; the partner's L-dependence awaits the depth-graded orders. β* and every bracket value are exact rationals (denominator 49 = 7²), no fitting.**

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, H, J, L-A, L-B, Real-T1, R1-STOP, Probe P/C/E/W/F2/G0/G1/G2, 2c(0+1), 2c0, 2c1. No `r_q` value changes; no rate-law fit (β* minimizes bracket WIDTH, not fitted to the partner; the partner is contained, not targeted). The scale-ledger negative (2c-1C) and the selection rule (2c0) are the guides for the deeper orders.

_Reporting discipline: β*=3/5 is reported as the exact line-search minimum (recognized rational, not asserted a priori — no magnitude pre-commit per the 0-for-20 ledger). The width shrink is reported with exact rationals (9/49). The dichotomy-killed claim is verified by the σ=±1 overlap, not assumed. The residual carrier is reported honestly: the max-spread CELL is still (odd,v₀) [tied with (even,v₀)], but the STRUCTURE is no longer the σ-dichotomy — it is the parity-blind {12/49,17/49,20/49} ladder; both facts stated. L-invariance is flagged (the corrector doesn't yet reach the partner's L-flow)._
