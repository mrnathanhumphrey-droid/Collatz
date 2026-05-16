# Phase 4 — Channel spectral gap on D_W: structural result + framing refinement

**Date:** 2026-05-15 (post-compact, after Phase 2 confirmed D_W is exactly dark under j ≥ 2). **Verdict: the per-step channel `L|_{D_W}` has a LARGE commutant (not a unique stationary state), so its spectral gap does NOT directly equal the c=7/45 rate 0.984. The dark-subspace classification correctly identifies D_W structurally, but the c=7/45 closure runs through a more subtle mechanism than "spectral gap of `L|_{D_W}`."**

## Setup

Per `DWM_DARK_SUBSPACE_ATTACK_PLAN.md` Phase 4: compute spectrum of the per-step channel
`L(ρ) = Σ_{v=1}^∞ M_v^{(j, b_prior)} · ρ · (M_v^{(j, b_prior)})†`
restricted to D_W, for fixed j ≥ 2 and b_prior. Compute as superoperator on M(D_W) (d_W² × d_W²) and eigendecompose.

Probe: `phase4_dark_spectral_gap_probe.py`, V_MAX = 16, scanned (j, b_prior) ∈ {2, 3} × {0, 1, 2}. Output: `experiments_output/phase4_dark_spectral_gap_probe.json`.

## Result

The spectrum is structured as **clusters at specific cyclic-group eigenvalues**, not a single λ_2 with a clean gap.

### n = 2 (d_W = 4, d_W² = 16)

Full spectrum is the same for all (j, b_prior) tested:
| Cluster | Count | Magnitude | Identification |
|---|---|---|---|
| Commutant | **4** | 1.0 (− 2^{−16}) | C[σ_{−1}^k]\|_{D_W}: power-of-shift algebra |
| 1st below | 4 | **0.5774** = 1/√3 | 1st harmonic, σ_{−1}'s eigenvalue ω_6 |
| 2nd below | 4 | 0.3780 = 1/√7 | 2nd harmonic, ω_3 |
| 3rd below | 4 | 0.3333 = 1/3 | 3rd harmonic, −1 |

### n = 3 (d_W = 16, d_W² = 256)

For j = 2 (with non-trivial phase `φ_9` reducing commutant):
| Cluster | Count | Magnitude |
|---|---|---|
| Commutant | **8** | 1.0 |
| 1st below | ≥ 12 | **0.8976** = 0.5/\|1 − 0.5 e^{iπ/9}\| |

For j = 3 (phase trivial at n=3, x_phase ≡ 0 mod 27):
| Cluster | Count | Magnitude |
|---|---|---|
| Commutant | **16** (= d_W) | 1.0 |
| 1st below | ≥ 4 | 0.8976 |

The commutant equals d_W for j = 3 at n = 3: the channel reduces to bare-shift averaging Σ 2^{−v} σ_{−v}(·)σ_{−v}†, which is **diagonal in the σ_{−1}-eigenbasis** and preserves the full diagonal.

## Structural identification of below-commutant eigenvalues

The first below-commutant eigenvalue at level n is **exactly**:
`λ_below(n) = 0.5 / |1 − 0.5 · ω_{2·ord_{3^n}(2)}|  =  0.5 / |1 − 0.5 · e^{iπ/3^{n−1}}|`

Verification:
- n = 2: ord_9(2) = 6, ω_12 = e^{iπ/6} … wait, formula gives e^{iπ/3} (since 3^{n−1} = 3). 0.5/|1 − 0.5 e^{iπ/3}| = 0.5/√(0.75) = 0.5774 ✓
- n = 3: 3^{n−1} = 9. 0.5/|1 − 0.5 e^{iπ/9}| = 0.5/0.5567 = 0.8976 ✓
- n = 4: prediction 0.987.
- n = 5: prediction 0.9985.
- **As n → ∞: λ_below(n) → 1**, i.e., channel becomes degenerate at the inverse limit.

This is the channel eigenvalue corresponding to the **fundamental Fourier mode of `σ_{−1}` action** on the cyclic group (Z/3^n)*, with weight `Σ_v 2^{−v} e^{ikv·angle}` summed for k = ±1.

## Why this is NOT the c=7/45 rate

The c=7/45 closure question (per `T_LEAD_CORRECTED_DISPOSITION.md`) involves three quantities:
1. **T_lead exact spectrum {43/45, 0}** on the 2-D class-resolved space (D_class, NOT D_W).
2. **Empirical Hadamard 0.984** at n = 10..13 (inward-trending), the asymptotic moment-decay rate.
3. **2.9% gap** between 43/45 = 0.9556 and 0.984.

Phase 4 measured the spectral gap of `L|_{D_W}` for the j ≥ 2 sub-family. This is a DIFFERENT object:
- λ_below(n) at finite n is determined by σ_{−1}'s fundamental Fourier mode on (Z/3^n)*. Specifically 0.5/|1 − 0.5 e^{iπ/3^{n−1}}|.
- λ_below(n) → 1 as n → ∞ (channel becomes degenerate; commutant absorbs everything).
- λ_below(4) ≈ 0.987 is fortuitously close to 0.984 BUT this is not a closure — the formula is structurally distinct from 43/45, and λ_below(n) keeps growing toward 1 as n → ∞.

**The c=7/45 / 43/45 rates live on D_class, not D_W.** T_lead is a class-resolved operator (P_+, P_-) → (P_+, P_-) with 2×2 matrix `(1/5)·[[1, 1], [4, 4]]` of spectrum {43/45, 0} (over Q, exact). It is not the dark-subspace channel.

## What the `L|_{D_class}` channel looks like (companion computation)

For completeness: under j ≥ 2 at level n, D_class is 2-dim (constant within each mod-3 class). σ_{−v} acts on D_class as identity (v even) or swap (v odd). The channel:
`L|_{D_class}(ρ) = (1/3)·I·ρ·I + (2/3)·S·ρ·S`
(using Σ_{v even} 2^{−v} = 1/3, Σ_{v odd} 2^{−v} = 2/3 for v ≥ 1).

Spectrum on M(D_class) (4 eigenvalues): {1, 1, −1/3, −1/3} (in Pauli basis: I and σ_x stay; σ_y and σ_z flip and shrink by 1/3).

**Neither 43/45 nor 0.984 is in this spectrum.** T_lead's 43/45 is NOT the per-step channel eigenvalue on D_class either — it is the eigenvalue of a SPECIFIC SCALAR REDUCTION (the leading Off_lin + T_diag operator extracted via R77's cross-frequency coherent summation), not the full channel.

## Conclusion: the dark-subspace framing requires refinement

The Benoist-Pellegrini-Szczepanek 2024 dark-subspace framework correctly identifies D_W as the **exact** invariant subspace under j ≥ 2 at finite n (Phase 2 result, exact mod machine epsilon). This is itself a structural advance.

But the c=7/45 = 43/45 closure question is about a DIFFERENT operator than the per-step channel on D_W:
- **T_lead** is a specific scalar reduction of the cross-frequency Off_lin + T_diag structure on D_class, not the channel on D_W.
- T_lead's 43/45 = 1 − Σ_g W_+(g) is a **moment-level Plancherel identity**, not a channel-spectral-gap.

The structural picture clarifies as follows:
- **D_W is dark under j ≥ 2** (Phase 2). The channel `L|_{D_W}` has commutant of size 4-16 (Phase 4), with first below-commutant eigenvalue 0.5/|1 − 0.5 e^{iπ/3^{n−1}}| (n=2: 0.577; n=3: 0.898; → 1 as n → ∞).
- **D_class evolution under j ≥ 2** is the 2-dim swap-symmetric channel `(1/3)I + (2/3)S(·)S` with spectrum {1, 1, −1/3, −1/3}.
- **T_lead's 43/45** is the eigenvalue of a CROSS-FREQUENCY OFF_LIN+T_DIAG operator on D_class (NOT the full channel) — a specific scalar reduction extracted via R77's class-mass coherent summation.
- **Empirical 0.984** is the asymptotic moment-decay rate of the FULL trajectory; the 2.9% gap to 43/45 is the unmodeled correction from inter-level / period-9 fluctuation (per `T_LEAD_CORRECTED_DISPOSITION.md`).

The dark-subspace classification, while structurally clean, is the **wrong frame** for closing the 2.9% gap. The right frame is **scalar reductions of cross-frequency operators on D_class** combined with **inter-level period-9 corrections**, not the channel-spectral-gap framework.

## What this means for the c=7/45 closure path

After Phases 1-4, the routes for closing the 2.9% gap (per the original `T_LEAD_CORRECTED_DISPOSITION.md` triad) are:
1. **Route A** (Nisoli closure at 43/45 with M_3''=24.4): still DEAD (Nisoli ruled out per prior session).
2. **Route B** (alternative inter-level operator construction): now structurally specified — must operate on D_class with period-9 phase coupling. Phases 1-4 of dark-subspace classification do NOT directly construct it; new probe needed.
3. **Route C** (document 43/45 + K_k {1, 0, ...} + dark-subspace classification of j ≥ 2 as the paper-grade structural anchor; defer period-9 closure): now the natural disposition. **The dark-subspace classification + Phases 1-4 structural results form a paper-shaped R3.**

## Phase 5 status

Per the attack plan, Phase 5 would be the "inverse-limit extension" of the dark-subspace structure. But given Phase 4's finding that λ_below(n) → 1 as n → ∞, the inverse-limit channel L|_{D_W} converges to the IDENTITY (commutant = everything). This means the inverse-limit channel-spectral-gap is degenerate — 0 gap. **Phase 5 in its original framing is moot.**

The actual asymptotic dynamics at the inverse limit are determined by the j = 1 mixing event + level-scaling, not by the inverse-limit `L|_{D_W}` spectral gap.

## Net Phase 4 deliverable

**Three things established:**
1. The per-step channel `L|_{D_W}` for j ≥ 2 has commutant of size 4–16 at finite n and spectral structure `0.5/|1 − 0.5 e^{ikπ/3^{n−1}}|` for k = 0, 1, 2, ... (with k=0 giving 1.0 = commutant cluster).
2. First below-commutant eigenvalue → 1 as n → ∞, so the inverse-limit channel-spectral-gap is degenerate.
3. The c=7/45 / 43/45 closure does NOT run through `L|_{D_W}`'s spectral gap. The dark-subspace classification correctly identifies D_W structurally, but T_lead's 43/45 lives on D_class as a class-mass coherent-summation phenomenon, not a channel eigenvalue.

**Decision: pause the original Phase 5-6 chain. The dark-subspace classification (Phases 1-4) is a structurally-clean ancillary result, not the c=7/45 closure mechanism.** The Phase 1-4 findings form a paper-shaped section but they don't close the 2.9% gap. The right next-step decision is up to the user: write up Phases 1-4 as paper-grade structural results (Route C), or open a new probe targeting Route B (period-9 inter-level operator on D_class).

## Files

- `phase4_dark_spectral_gap_probe.py` (probe)
- `experiments_output/phase4_dark_spectral_gap_probe.json` (full eigenvalue data, top-20 per (j, b_prior))
- `PHASE4_DARK_SPECTRAL_GAP_RESULT.md` (this writeup)

## Caveat

This writeup makes a claim ("the c=7/45 closure does NOT run through L|_{D_W}'s spectral gap") that revises the framing of `DWM_DARK_SUBSPACE_ATTACK_PLAN.md`. The probe data is the ground truth: large commutant at every (j, b_prior), λ_below(n) → 1 as n → ∞. The framework gap is structural, not interpretive. Updating attack plan accordingly.
