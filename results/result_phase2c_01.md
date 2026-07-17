# Result — PHASE 2c stages 0+1: the defect object is FROZEN (reconstruction exact, amplitude = |c_k|/c₀) and the ZEROTH BRACKET is [0.1184, 0.6041] (width 0.486, L-invariant, ∋ 1/3 & partner). The named 2c-3 obstruction is CONFIRMED numerically: defect amplitudes do NOT decay (0.974 at D=18→1) and the defect Frobenius mass GROWS relative to the mean-field (2.3→9.8) — amplitude decay cannot carry the contraction; depth is the only route.

**Date:** 2026-07-16. Phase 2c stages 0 (freeze defect object) + 1 (zeroth Collatz–Wielandt bracket). Claude derives / code gates. Direct/exact at q=3 (INSTRUMENT LAW). Probe `probes/probe_phase2c_01.py`, log `logs/probe_phase2c_01_log.txt`. No proof authored, no rate fit. Per-state claims labeled ALGEBRAIC (residue-forced, exact) vs STATISTICAL (class-averaged only).

**Headline: (2c-0) the defect object freezes cleanly — M̄ = gauge-average closure (the (e_ρ,γ)-compressed chain lifted, = k=0 unit-character component), D = M_tower − M̄ character-expands over the unit group with Parseval reconstruction EXACT (residual 1e-13), and the per-k amplitude spectrum ‖D_k‖/‖D_0‖ EQUALS |c_k|/c₀ (the twisted-autocorrelation ratio) — exact at D=18. (2c-1) the zeroth bracket [min,max] (Mᵀh⁰)/h⁰ = [0.11836735, 0.60408163], width 0.486, L-INVARIANT, contains both 1/3 and the partner — the baseline every correction order must beat. The named 2c-3 obstruction is CONFIRMED in two independent measures: |c_1|/c₀ = 0.832 (D=6) → 0.974 (D=18) → 1 (non-decaying), AND the defect/mean-field Frobenius-mass ratio GROWS 2.34 (L=2) → 9.80 (L=3). Amplitude decay is dead as a contraction mechanism; depth-weighting is the only remaining route.**

## 2c-0a — scope correction on G0-2 (per-STATE survival). Only ONE cell is statistical.
G0-2's "spread 0" was the CLASS-averaged (compressed-chain) survival. Per-STATE (individual (a,b,γ)):
| cell | per-state survival | verdict |
|---|---|---|
| (even, v₃=0) | 2/9 exactly (spread 0) | **ALGEBRAIC** |
| (even, v₃≥1) | 5/9 exactly (spread 0) | **ALGEBRAIC** |
| (odd, v₃≥1) | 4/9 exactly (spread 0) | **ALGEBRAIC** |
| **(odd, v₃=0)** | **[1/9, 4/9], mean 5/18 (spread 1/3)** | **STATISTICAL** |

- **Three of four cells have algebraically-constant per-state survival; only (odd, v₃=0) varies per state** — range [1/9, 4/9], class-mean 5/18, half-range **1/6 = 0.1667 ≈ Wilson's flagged ≈0.17** (full spread 1/3). So G0-2's cell survivals {2/9, 5/18, 5/9, 4/9} are exact per-state for 3 cells and a class-average for the 4th. The STANDING RULE (algebraic vs statistical) is now anchored: (odd, v₃=0) survival is the sole statistical row-sum, and must never be used per-state.

## 2c-0b — the defect object, frozen. Reconstruction exact; amplitude = |c_k|/c₀.
`M̄[dst,src] = (1/D)Σ_{s∈⟨2⟩} M[dst, s·src]` (gauge s: (a,b,γ)→(sa,sb,γ); orbit = fixed (e_ρ,γ), varies e_a). This is the (e_ρ,γ)-compressed chain lifted = the k=0 unit-group character component. `D = M_tower − M̄ = Σ_{k≠0} D_k`, `D_k` the twist-k character component, `‖D_k‖²_F = (1/D)Σ_orbit‖S_orbit v_k‖²`.
| L (D) | Parseval `Σ_k‖D_k‖² = ‖M_tower‖²_F` | residual | ‖M̄‖²=‖D_0‖² | defect Σ_{k≠0}‖D_k‖² | defect/mean-field |
|---|---|---|---|---|---|
| 2 (6) | 10.64499874 = 10.64499874 | 1.8e-15 | 3.19 | 7.46 | **2.34** |
| 3 (18) | 306.00466923 = 306.00466923 | 1.1e-13 | 28.34 | 277.67 | **9.80** |

- **Reconstruction EXACT (residual ≤1e-13) — M̄ + Σ_{k≠0}D_k = M_tower to machine precision, both L.** The defect object is frozen and complete.
- **Amplitude formula CONFIRMED: the per-k amplitude spectrum √(‖D_k‖²/‖D_0‖²) equals |c_k|/c₀ = |Σ_δ w_δ² ω^{kδ}|/Σ_δ w_δ² — EXACTLY at D=18** (k=1..6: 0.974, 0.910, 0.832, 0.759, 0.700, 0.655, matching |c_k|/c₀ to the digit). At D=6 the amplitude ratios (0.789, 0.621, 0.569, …) track |c_k|/c₀ (0.832, 0.655, 0.600) more loosely — small-tower boundary effect. **The defect amplitude per k IS the twisted-autocorrelation ratio** (the R_k family, Real-T1's object, reappearing as the tower defect amplitude).
- **The defect DOMINATES and GROWS: defect/mean-field Frobenius² = 2.34 (L=2) → 9.80 (L=3).** The mean-field M̄ is a minority of M_tower, and its share shrinks with L. This is the same fact as |c_k|/c₀→1, seen in the operator norm.

## 2c-1 — the zeroth Collatz–Wielandt bracket. [0.1184, 0.6041], the baseline.
h⁰ = coarse Perron right-eigenvector lifted by cell: (even,v₃=0)→2/3, (even,v₃≥1)→5/3, (odd,v₃=0)→5/6, (odd,v₃≥1)→4/3 (verified `T h⁰ = (1/3)h⁰` on G1's 4×4). Operator = Mᵀ (matching G1's `T[src,dst]=flow(src→dst)` convention; `M h⁰` would be inflow-weighted with degenerate zero rows).
| L | bracket [m₀, M₀] | width | 1/3 inside | partner inside |
|---|---|---|---|---|
| 2 | [0.11836735, 0.60408163] | 0.4857 | ✅ | ✅ |
| 3 | [0.11836735, 0.60408163] | 0.4857 | ✅ | ✅ |

- **Pre-registration met: 1/3 is inside, and so is the partner (ρ_L), at both L.** Collatz–Wielandt guarantees ρ_L ∈ [m₀, M₀]; the width 0.486 is the zeroth defect scale — **the baseline every correction order must beat.**
- **The bracket is L-INVARIANT** (identical [0.11836735, 0.60408163] at L=2 and L=3, shrink factor 1.00×) — the extreme ratios come from an L-invariant local configuration (min at (even,v₃≥1), max at (even,v₃=0)). **The zeroth bracket does NOT shrink with L**; the shrink must be produced by the corrector (2c-2), and — given 2c-0b — the corrector's contraction cannot come from amplitude decay.

## Adjudication
| stage | verdict |
|---|---|
| 2c-0a scope | 3 cells ALGEBRAIC (exact per-state survival), (odd,v₃=0) STATISTICAL (±0.167 ≈ Wilson's 0.17). |
| 2c-0b defect | frozen; reconstruction EXACT (1e-13); amplitude = \|c_k\|/c₀ (exact @D18); defect dominates & grows (2.3→9.8). |
| 2c-0b obstruction | CONFIRMED — amplitudes non-decaying (0.974→1) in TWO measures; amplitude route dead, depth is the only mechanism. |
| 2c-1 bracket | [0.1184, 0.6041], width 0.486, L-invariant, ∋ 1/3 & partner. Baseline for 2c-2. |

**⟹ Stages 0 and 1 land as designed and produce banked value independent of the crux. The defect object is a clean, complete, gauge-character decomposition whose amplitudes are exactly the R_k twisted-autocorrelation ratios; the zeroth bracket is the L-invariant baseline (width 0.486). Most importantly, the named 2c-3 obstruction is now numerically nailed from two sides — the character amplitudes |c_k|/c₀ → 1 AND the defect Frobenius mass grows relative to the mean-field — so 2c-2's corrector and 2c-3's contraction estimate CANNOT rely on amplitude decay. The depth-weighted norm (a defect at digit depth j diluted by the cascade's 3⁻ʲ, per the confirmed D2 digit law) is the only route left, exactly as Wilson anticipated.**

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, H, J, L-A, L-B, Real-T1, R1-STOP, Probe P/C/E/W/F2/G0/G1/G2. No `r_q` value changes; no rate-law fit (the 2.9e-3, 1.0e-4 sequence untouched; the bracket contains the partner but is not fitted to it). G0-2's cell survivals stand as class-averages; 2c-0a refines them to per-state (3 algebraic + 1 statistical). D2/U0 (G2) feed 2c-2/2c-3 unchanged.

_Reporting discipline: the bracket convention error (M vs Mᵀ) was caught and fixed before banking — the degenerate [0, 2.28] was a wrong-direction artifact, the correct bracket is [0.118, 0.604]. The scope correction is reported as a genuine refinement of G0-2 (only 1 cell statistical), with Wilson's ≈0.17 identified as the half-range. The amplitude=|c_k|/c₀ match is reported exact at D=18 and loose at D=6 (not overclaimed as exact at both). The obstruction is reported as CONFIRMED (numerically, two measures), framed as it kills the amplitude route — a load-bearing negative for 2c-3, not a failure. Algebraic/statistical labels applied per the standing rule._
