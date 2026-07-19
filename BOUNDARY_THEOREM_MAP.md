# BOUNDARY THEOREM — theorem → gate bridge

Maps the master write-up (`D:/Resolve Research/Collatz Documents/Papers/boundary_theorem_master.md`,
Thm 1–6) to this repo's verification record. Wilson's pen writes the theorems; these are the probes/gates that
judge them. Status key: **PROVEN** (transcribe from source) · **OWED** (pen derivation remaining) ·
**NUMERICAL** (pre-registered, commit-stamped). Live chronological log: `STATE.md`.

| Thm | Statement (short) | Status | Gate / probe | Result doc |
|---|---|---|---|---|
| **1** | Kinematic spectrum = {c_k}, c_k = Σw²ω^{kδ}; closed-form ℓ_k; rest nilpotent | PROVEN | Real-T1, gate 18/18 (commit `1fe5a64`) | result_phase2b_H / LALB |
| **2** | spec M = {c_k} ∪ {0} ∪ spec M_tower; partner = ρ(M_tower) | PROVEN | PARTNER-CHAR, G0 @ 1e-14 both L | result_phase2c0_* |
| **3** | Crossing: c₀(λ)=(1−λ)/(1+λ) moves, ρ(S)=1/3 weight-free ⟹ collide at λ=½ | PROVEN | G1 arc + **G-D0.3** gate (probe_trackD0_gate, λ=0.4/0.5/0.6) | result_trackD0_gate |
| | └ numerical companion: partner pinned 1/3 within 6e-4 while c₀ tracks | NUMERICAL | D-1 crossing rider (probe_trackD1_fork) | result_trackD_fork |
| **4** | Ladder **labeling** total {±1,±2,±4}+DC (626/626 L3); triangular selection 3^{L−1−j}\|κ; **NOT a reduction** (purity 0.07–0.32) | PROVEN (totality+rule); scope from measured purity | **J1** (probe_judge_J1) + 2c0-G3 | **result_judge.md §J1** |
| **5** | σ(θ)=(1/3)((1+e^{iθ})/2)²; bulk=½-flux, edge=Lebesgue restriction; dominant mode = d within couplings ≤3.6e-5 | PROVEN (sketch rigor); **Lemma D OWED** (closed-form d + march) | D2-a ½-flux (exact), D2-b/c edge, **D2-e** ladder matrices, **J2** couplings | result_trackD2a / result_trackD2e / result_judge.md §J2 |
| **6** | Coalescence at (3,½): condensation 0.681→0.900→0.987, phase→2π/3^{L−1}, gap (1/3)sin²(θ/2) @2% by L4; braiding; overlap→1/defect 17→189 | NUMERICAL spine complete; **rides on Lemma D** | G4 braid, D-1 EP recon, W4 witnesses, **J2** arm + **rider** (12-digit L4 doublet) | result_phase2d_G4 / result_trackD1 / result_W4 / result_judge.md §J2 |
| **7** | **Blindness:** spec(M₋)=spec(M₊) exactly — sign-blind; the distributional→pointwise barrier, located & proven | PROVEN (exact perm similarity); **lemma corrected** | **B1** — M₋=P M₊ P (pair-swap) EXACT 0.0; partner 0.34682666/0.33323630 reproduced | **result_blindness_B1.md** |

**Thm 7 mechanism note (walk-back #25):** the intertwiner is the **pair-swap** P:(a,b,γ)↦(b,a,γ), **not** the
committed negation Σ:(a,b,γ)↦(−a,−b,−γ). Σ fails (M₋≠ΣM₊Σ, 0.25–0.26) on carry-floor ⌊·⌋ vs modular-negation
non-commutation — the same breakage as the J-involution (walk-back #14); the durable constraint "invariance must
act on the carry as an INTEGER map, not modular" selects P (an integer relabel that sends T→−T = the +↔− flip).
The claim is unchanged and cleaner (exact permutation similarity).

## OWED (Wilson's pen — the two remaining derivations)
- **Lemma D** — closed-form d(gf, L) (exact lattice sum) with σ as its limit and the correction law reproducing
  the phase march 0.705 → 0.940 → 0.9928. Feeds Thm 5 limit + Thm 6 condensation.
- **Lemma B** — diagonal-dominance bound for dominant modes (empirically 2.5e-7 → 3.6e-5) from the selection
  structure. Feeds Thm 5's "dominant mode = d" claim.

## Judge (Probe J) — the verdict, unsealed 2026-07-18 (commit `d7692ab`)
- **J1 completeness (L=2,3, dense):** totality ✅ (all 626 modes → {±1,±2,±4}+DC, no 4th family); jury modes
  captured 1e-3→1e-2; **strong pre-reg FAILS** for subdominant (purity 0.07–0.32, resid→mode-scale) ⟹ ladder =
  complete labeling + dominant-mode theory, **not** a spectral reduction. (Thm 4's scope clause is this result.)
- **J2 arm (L=4, SpMV, no dense eig):** PASS ×3 — (i) |d(1,4)|/σ(θ₁)=1.00059 ∈(1,1.0201), dressing
  +2.01%→+0.059%; (ii) phase 0.9928 ∈(0.9434,1); (iii) couplings 3.6e-5 « L3. Cross-check: reproduces banked
  block-6 L4 doublet to 1.6e-4.
- **Rider (12-digit L=4 doublet split):** RUNNING (within-block subspace iteration; near-EP floor-limited to
  ~10 digits). Feeds §6 block-splitting ledger (target split ratio ×0.076). *(append on landing)*

## Verification ethos (§8)
Derive-blind / machine-gate / adjudicate; magnitudes pre-registered before the runs that judge them; the
walk-back ledger (24 logged retractions, public); instrument law near the EP (dense/direct or within-block
power iteration; NO ARPACK/shift-invert; exact-rational; no rate-fitting). Honest negatives ARE the win —
the 27^{−L} law (44× miss), the ±4-seat withdrawal, and J1's strong form are all in the record.
