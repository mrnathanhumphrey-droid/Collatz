# Result — MICRO-PROBE G0: PARTNER-CHAR CONFIRMED — the dynamical partner IS the Perron root of the carry-tower principal submatrix, EXACTLY (8e-16 / 2e-14 at L=2/3). G0-2: C2's survival range reconciles cleanly (parity-averaging), but the pre-registered {4/27,5/27,4/9,5/9} is REFUTED — the true exact survivals are {2/9, 5/18, 4/9, 5/9} with the v₃-dependence REVERSED.

**Date:** 2026-07-16. Cheap decisive gate. Direct/dense at q=3 (INSTRUMENT LAW). Probe `probes/probe_phase2b_G0.py`, log `logs/probe_phase2b_G0_log.txt`. No proof authored, no rate fit.

**Headline: the partner is fully in hand as a spectral-radius. Deleting the γ=0 sector from the full operator and taking the top eigenvalue of the remaining principal submatrix gives the true partner to machine precision — 0.34682666 (L=2, diff 8.3e-16) and 0.33323630 (L=3, diff 2.3e-14) — and it is genuinely the Perron (next tower eigenvalues are complex, 0.022 / 0.238). The partner = ρ(M_tower), no subspace extraction needed. G0-2 separately: C2's 0.25–0.50 survival range is EXACTLY the (e mod 3) parity-average of the fine survivals (0.25=(2/9+5/18)/2, 0.50=(5/9+4/9)/2) — a clean grouping effect, no normalization mismatch. But the pre-registered fine values {4/27,5/27,4/9,5/9} are REFUTED: the true values (exact, spread 0, no O(2⁻ᴰ) correction) are {2/9, 5/18, 4/9, 5/9}, and the v₃-direction is reversed — v₃(γ)≥1 states survive MORE, not less.**

## G0-1 — THE theorem gate: partner = Perron(M_tower). CONFIRMED.
`M_tower` = the principal submatrix of the full pair operator on the γ≠0 (carry-tower) states; direct dense eig; ρ = top eigenvalue by modulus.

| L | dim(full) | dim(tower) | c₀ | true partner (full, subspace) | ρ(M_tower) | \|ρ − partner\| | is the Perron |
|---|---|---|---|---|---|---|---|
| 2 | 324 | 288 | 0.34391534 | 0.34682666 | **0.34682666** | **8.33e-16** | yes (next: 0.022±0.24j) |
| 3 | 8748 | 8424 | 0.33333588 | 0.33323630 | **0.33323630** | **2.31e-14** | yes (next: 0.238±0.18j) |

- **PRE-REGISTERED EXACT — CONFIRMED at both L, to machine precision.** ρ(M_tower) equals the true partner to 8e-16 (L=2) / 2e-14 (L=3), and it is the dominant eigenvalue of the tower block (the next tower eigenvalues are complex and well-separated: |·|≈0.022 at L=2, ≈0.238 at L=3).
- **Why it holds (F2-4 + P, now a clean spectral statement):** the partner's right eigenvector is zero on γ=0 (Probe P) and c₀'s kinematic mode never feeds the tower (F2-4: B[kin,tow]=0). So `M·r_p = λ_p·r_p` with `r_p|_{γ=0}=0` restricts exactly to `M_tower·(r_p|_tower) = λ_p·(r_p|_tower)` — the partner is an eigenvalue of the tower principal submatrix, and G0-1 confirms it is the TOP one.
- **⟹ PARTNER-CHAR CONFIRMED: the dynamical partner IS ρ(M_tower).** The last object of the entrance exam is now a spectral radius of a concrete, γ=0-deleted principal submatrix — "the one thing near 1/3 that is not an autocorrelation" is the Perron root of the carry tower. No EP-fragile extraction is needed to define it (though the tower block is still defective-adjacent, so extraction near it stays direct/LU per the instrument law).

## G0-2 — tower survival row-sums by (e_ρ parity, v₃(γ)). C2 RECONCILED; pre-reg VALUES refuted.
Row-sums (survival = total per-state out-weight passing the gate) of the uniform compressed chain, tower classes (γ≠0), grouped by (e_ρ mod 2, v₃(γ)=0 vs ≥1). **Spread within each group ≈ 1e-16 ⟹ exact rationals, no finite-L correction.**

| e_ρ parity | v₃(γ) | survival (exact) | pre-registered |
|---|---|---|---|
| even | 0 | **2/9** = 0.22222 | (was 4/27) |
| odd | 0 | **5/18** = 0.27778 | (was 5/27) |
| odd | ≥1 | **4/9** = 0.44444 | (was 4/9, wrong cell) |
| even | ≥1 | **5/9** = 0.55556 | (was 5/9, wrong cell) |

- **C2 reconciliation — CLEAN (bookkeeping, as hoped).** C2 grouped by (θ=e mod 3, γ); each θ-class lumps 3 even + 3 odd e_ρ values (e_ρ ∈ {r, r+3, r+6, …}), averaging the even/odd survival alternation. The result is EXACT: C2's 0.25 = (2/9 + 5/18)/2 (the v₃=0 parity-average), C2's 0.50 = (5/9 + 4/9)/2 (the v₃≥1 parity-average). Verified to machine precision. **C2's survival is the SAME row-sum object — no normalization difference; the 0.25–0.50 range is purely the (e mod 3) parity-washing.** The fine (parity, v₃) split resolves the alternation C2 averaged away.
- **Pre-registered SET {4/27, 5/27, 4/9, 5/9} — REFUTED (real error in the alternation derivation, NOT bookkeeping).** Two things are wrong: (i) the low-cell values are **2/9 and 5/18, not 4/27 and 5/27** (the predicted 1/3-suppression of the v₃=0 cells is absent); (ii) the **v₃-dependence is REVERSED** — v₃(γ)≥1 states survive MORE (4/9, 5/9), v₃(γ)=0 survive LESS (2/9, 5/18). Mechanism: at 3|γ the gate `(γ+T)≡0 mod 3` forces `T≡0 mod 3` = same-parity moves = the heavy-weight channel, so 3|γ states retain more mass, not less. The alternation derivation's cell assignment and v₃-direction need correction (flagged for Nathan). The **correct exact target set is {2/9, 5/18, 4/9, 5/9}.**
- **Bonus structural fact:** survival depends only on (parity, whether 3|γ) — v₃=1 and v₃=2 classes give identical survival (the ≥1 bucket has spread 0). Depth beyond the first carry level does not change the survival, only the parity and the 3|γ indicator.

## Adjudication
| item | verdict |
|---|---|
| **G0-1 partner = Perron(M_tower)** | **CONFIRMED** — exact to 8e-16/2e-14, and it is the top of the tower block. Partner fully characterized as a spectral radius; object in hand. |
| G0-2 C2 reconciliation | CLEAN — C2's 0.25/0.50 = exact (e mod 3) parity-averages of the fine survivals; no normalization difference (bookkeeping/grouping). |
| G0-2 pre-reg values {4/27,5/27,4/9,5/9} | REFUTED — true exact set is {2/9, 5/18, 4/9, 5/9}; v₃-dependence reversed (3\|γ survives more). Real alternation-derivation error, flagged. |

**⟹ The partner now has a clean operational definition — ρ(M_tower), the Perron root of the γ=0-deleted principal submatrix — confirmed exact at L=2,3. The compressed-chain survival structure is pinned to exact rationals {2/9, 5/18, 4/9, 5/9} (parity × 3|γ), which both explains C2's coarse 0.25–0.50 (parity-average) and corrects the pre-registered alternation values. Wilson's alternation derivation should be re-cast on {2/9, 5/18, 4/9, 5/9} with v₃(γ)≥1 as the high-survival cell.**

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, H, J, L-A, L-B, Real-T1, R1-STOP, Probe P/C/E/W/F2. No `r_q` value changes; no rate-law fit (the true 2.9e-3, 1.0e-4 sequence untouched). E-FORM (F2-1) unaffected — the survivals here are row-sums of the same gate-confirmed entries.

_Reporting discipline: G0-1's pass is reported with the exact residual AND the is-Perron check (top-of-block, not merely present). G0-2 is reported as a SPLIT outcome — C2 reconciliation clean (exact parity-average), but the pre-registered values REFUTED and named as a real error (not smoothed into "consistent with"); the corrected exact set and the reversed v₃-direction are stated plainly with the gate mechanism. The exact-rational status (spread 0) is verified, so the values are not called "approximate up to O(2⁻ᴰ)" — they are exact at this grouping. c₀-masquerade / instrument law respected throughout._
