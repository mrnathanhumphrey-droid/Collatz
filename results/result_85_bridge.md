# Result 85 rungs 1–2 — DWM operator-DFT bridge + cyclotomic-"7" test (disposition)

**Date:** 2026-07-14. **Verdicts: rung 1 = operator-level chirp identity CONFIRMED (single-r, unextended); rung 2 = cyclotomic-"7" REFUTED.**

Probes: `probe_85_operator_dft.py`, `probe_85_orbit_structure.py`, `probe_85_rung2_eisenstein.py` (axis A), `probe_85_rung2_axisB.py` (axis B), `probe_85_rung2_dinf.py` (D_inf). Logs: `result_85_operator_dft_log.txt`, `result_85_rung2*_log.txt`. Builds on R82 (H_BRIDGE_PARTIAL), R81/R81b (certified F̂ phase), R85 Phase 0 (feasibility).

## Setup — the placeholder was invalid; measure fixed it

R85 Phase 1's first evidential comparison used a placeholder DWM chirp multiplier (=1) and compared pointwise on ⟨4⟩. `probe_85_orbit_structure.py` measured the *real* DWM j=2 chirp argument `Δ = 2^{−b}(2^{−v}−2^{−v'})` and found it spans the **whole ring** Z/3^{r+1} (34.6% ⟨4⟩, 34.6% 2⟨4⟩, 30.8% non-units, identical at n=5,6) — so a pointwise chirp-vs-F̂ swap is **category-invalid** (2/3 of the DWM arguments aren't even in F̂'s ⟨4⟩ support). The valid test is a Fourier one.

Key algebraic fact: mod q=3^{r+1} the unit group is Z/2 × Z/3^r, and ⟨4⟩=⟨2²⟩ is the *entire* Sylow-3 part; R81's family c_{ℓ,ε}=2^ε(1+3^r)^ℓ collapses to **2 chirp-classes** (ε→coset), the ℓ-index being only the R84 global ω₃^ℓ phase. So the pointwise chirp-identity is *forced* by R82 step-1 (same character species) — the earlier "same for ε=0 / differs for ε=1" pattern was detecting the coset bit, not the bridge.

## Rung 1 — operator-DFT bridge (n=6, r=3)

The real DWM j=2 single-Kraus chirp, on its ⟨4⟩ sub-orbit (even v=2k), is `e_q(c'·4^{−k})`, c'=−xi·2^{−b}. DFT'd and run through R81's own cyclotomic certification (`result_81_fhat_phase_profile.py`) as an independent oracle:
- magnitude **3√q=27 flat** on a size-9 support (Th 78.3 holds for the DWM chirp),
- phase certifies, **degree 3** (R81's growing-degree law, same value),
- 1st-difference J₄ multiset matches R81's family reference on both cosets.

**⇒ the DWM operator's chirp is a bona-fide R78/R81 certified chirp** — R82's formula-level common-form lifted to an operator-level, independently-certified statement. **Caveat:** single-r (n=6 only); checks degree + 1st-diff multiset, not the full J₄ index. Every unit multiplier gives this structure, so it pins *species*, not discriminating family-membership. Same small-window profile the arc killed 3× (H_QUAD, ⌊r/2⌋+2, mod-9 offset) — **owed an r=5 (n=8) extension before trusted.**

## Rung 2 — is the "7" the cyclotomic norm N(2−ω)? REFUTED

The discriminating content is the measure-assembly. Tested whether the "7" in 7/45 is N(2−ω)=Φ_p(2)=2^p−1 by generalizing the multiplier 3→p (base-2 halving fixed; ω→primitive p-th root):
- **Axis A** (vary halving base m): S_k(m) diverges off m=2 — that axis breaks S_k convergence (test-design flaw, not a verdict).
- **Axis B / D_inf** (vary multiplier p): `S_k(p)~(p/3)^k`, so the convergent observable is `D_inf(p)=lim(3/p)^k S_k` (anchored D_inf(3)=7/15). Back-out D(p)=N_p/D_inf:

| p | D(p) | 2^{p+1}−1 | ratio |
|---|---|---|---|
| 3 | 15.03 | 15 | 1.00 |
| 5 | 63.4 | 63 | 1.01 |
| 7 | 162.7 | 255 | **0.64** |
| 11 | 2809 | 4095 | **0.69** |
| 13 | 10641 | 16383 | **0.65** |

The D=15,63 (=2⁴−1,2⁶−1) fit at p=3,5 is a **2-point coincidence**; p=7,11,13 break it ~0.65×. **7=Φ_3(2) is numerical coincidence, not mechanism** — the "7" stays R77's T_diag (1,4)-projection mass (7=1−8/15). Same death as inverse-tree 1/9 and ⌊r/2⌋+2: small-window fit that dies on extension.

**Distinct from R1B:** R1B confirmed the *measure*-robustness of 7=N(2−ω) (H_ROBUST, invariant marginal exact to 1e-7). Rung 2 refuted the *multiplier*-generalization 7=Φ_p(2). Different claims — the mechanism holds against the true measure; it does not generalize across multipliers.

## Byproduct — the rate 1/3 is multiplier-independent (handed to parallel)

`S_k(p)~(p/3)^k` ⇒ `‖d_k(p)‖²~3^{−k}` at rate **1/3 for every multiplier p** (p=3,5,7,11,13). This confirms STATE's two-month-old open "why q/3? — publishable theorem candidate" from an independent probe. The *rate* generalizes (base-2 halving vs 3-adic modulus); the *numerator* does not (per rung 2). **Derivation (separation-of-variables) handed to the parallel agent.**

## Not at stake

THEOREM_C_745 (c=7/45), Th 78.1–78.3, R81/R81b certification — all untouched. Rung 1 is a positive on the R82 bridge; rung 2 is a clean negative on a spine-decoration; neither bears on the rigorous constant.

_Reporting discipline: rung-1 caveat (single-r, species-not-family) stated up front, not buried; rung-2 refutation stated as a refutation; the p=3,5 near-fit reported AND flagged as a 2-point coincidence, not softened._
