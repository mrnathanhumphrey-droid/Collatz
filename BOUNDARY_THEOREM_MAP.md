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
- **Rider (12-digit L=4 doublet split):** LANDED, **12 digits** (converged it 259, Δ=7.5e-14). Members
  0.320422712770+0.075242317692j & 0.320222549235+0.075251807019j; **split 2.00388e-4, ratio-to-L3 0.0758 ≈
  ×0.076** (block-splitting ledger confirmed). Partner 0.333499901322 = g4 ρ₄ to 12 digits. §6 entry complete.

## THREAD 3 — the constant (7/15, 7/45) and the loss ledger (Theorem F + the off-diagonal derivation)
The paper's constant, welded to the frozen instrument and its derivation reduced to a finitely-computable ledger.
Shell sequence S_k (= corpus per-scale L² cost, backward convention S_k = 3^k(a_k − a_{k−1}/3) = A_k − A_{k−1}).

| item | statement | status | gate / probe | result doc |
|---|---|---|---|---|
| **Theorem F (flat level)** | S_k → S∞ = 7/15 = 3·(7/45); three faces (spectral ρ=1/3, renewal driftless, Fourier non-decaying shells) | **welded** (S₁=2/3, S₂=10/21 exact; asymptote L-truncated at finite L) | **R1** re-weld + boundary contrast; **R2** A∞=(3/2)S∞=7/10 | result_thread3_R1 / result_thread3_R2 |
| A∞ = 7/10 (amplitude) | c₀-mode overlap of the independent pair | value holds; **ℓ₀-route DEAD** (2531/4095≠7/10) → it's the Jordan coupling (3/2)·β, not isolated ℓ₀ | **R2-B** (walk-back on crown route) | result_thread3_R2 §R2-B |
| **Theorem S (secular)** | S∞ = 3·g·φ_tow·ψ_kin; diverging g/Δ × vanishing Δ cancel | **mechanism gated** (Δ-cancellation L=3 ratio 0.974, L=4 0.9916); **L→∞ closed form OWED** | **R3** (L=4 product discarded: near-EP underconvergence) | result_thread3_R3 |
| edge-density law | 7/15 = L→∞ 1/θ edge residue (Dirichlet integral) | **NOT numerically confirmable** (finite-L band discrete, EP wall) — must close symbolically (Ĝ) | **R4** (scaling laws refuted at finite L) | result_thread3_R4 |
| dark-state selection | some band modes unread by the agreement functional | **readout zero** ⟨1\|r⟩=0 (DC-free), NOT a symmetry (no involution commutes; P dies #29) | **D1** (walk-back #29) | result_darkstate_D1 |
| deviation law | S_k = 7/15 + d_k | candidate 1/(5·21^{k−1}) **DEAD** at k=3 (wrong sign); true d_k signs +,+,−,−,−,− (overshoot) | **R5** (exact S₁..S₆ from Basic.lean) | result_deviation_R5 |
| **off-diagonal ledger** | S_k = S_{k−1} + OffDiag_k (diagonal replicates via 3×⅓=1); Σ_{k≥2} OffDiag = −1/5 | **R6-A GATE PASS**: OffDiag₂ = −4/21 derived from v≠v′ Ramanujan sums; diagonal replicates exactly k=2,3 | **R6** | result_offdiag_R6 |
| **channel engine (⟨4⟩-orbit law)** | OffDiag_k = (2/3)Σ_m 4^{−m}C_k(m), C_k(m) = twisted ⟨4ᵐ⟩-orbit character sum of μ_{k−1}, period 3^{k−1} | **R7 FULL GATE PASS k=2,3,4,5** (engine = frozen); sign-flip DERIVED (C_k(1) crosses zero); odd gaps ≡0; Mersenne denoms 4^P−1 | **R7** | result_engine_R7 |
| **uniform kill + strata** | correlation lives in μ non-uniformity (uniform μ ⟹ C_k≡0); OffDiag_k = (2/3)Σ_j W_j C̄_k(j), W_j k-independent | **R8 ALL PASS**: uniform kill C_k≡0 exact; band count m-independent (3^{k−1}, 2·3^{k−1}); ledger↔deviation weld exact; W_j closed form; overshoot = C̄_k(0) crosses zero k=3→4 | **R8** | result_strata_R8 |
| **collision identity (γ on ℤ₃)** | S_K = 2Σ_m 4^{−m}γ_{K−1}(τ_m), γ_n(τ)=collision density, γ_n(0)=X_n; whole campaign = one stationary γ_∞ | **R9 FULL PASS**: identity exact K=2..6 (from μ tables); γ_n(0)=X_n weld (τ=0 = qx+1 corpus); DC self-similarity C_k(DC)=3S_{k−1}; off-DC bounded; P_n=S_n, ⟨γ_n⟩=3S_{n+1}/2 | **R9** | result_gamma_R9 |

**Thread-3 net:** the constant's derivation is reduced to **Σ_{k≥2} OffDiag_k = −1/5**, a finitely-computable
exact-rational ledger. The **channel engine (R7)** now derives every term through k=5 from first principles:
OffDiag_k = (2/3)Σ_m 4^{−m}C_k(m), with C_k the twisted ⟨4ᵐ⟩-orbit character sum of μ_{k−1} (period 3^{k−1}),
engine == frozen for k=2,3,4,5. The two-sign tail mechanism is derived, not observed: the gap-2 sign flip
(−1/6 → +2/147) is **C_k(1) crossing zero** (C₂(1)=−1, C₃(1)=+4/49); odd-gap channels vanish identically
(conjugate-kill, verified); a positive DC/self-orbit class of collapsing Mersenne weight 1/(4^{3^{k−1}}−1) races
the negative bulk to produce the −,−,+,+ overshoot. The **uniform kill + strata (R8)** sharpen this: the whole
correlation lives in μ's non-uniformity (uniform μ ⟹ C_k ≡ 0 exact), the affine-band pair-count is m-independent
(3^{k−1}, 2·3^{k−1}), and OffDiag_k = (2/3)Σ_j W_j C̄_k(j) with **W_j k-independent** (closed form x/(1−x)−x³/(1−x³),
x=4^{−3^j}) — so the limit law is entirely C̄_∞(j), and the −,−,+,+ overshoot is **C̄_k(0) (the dominant W₀=20/63
stratum) crossing zero between k=3 and k=4**. The **collision identity (R9)** completes the reduction: the whole
campaign is **one stationary function γ_∞ on ℤ₃** — S_K = 2Σ_m 4^{−m}γ_{K−1}(τ_m) exact K=2…6, the τ=0 line welds
to the qx+1 corpus (γ_n(0)=X_n), DC self-similarity is literal (C_k(DC)=3S_{k−1}), and off-DC columns are bounded
(divergence confined to τ=0). **Still owed (pen):** γ_∞(τ) closed form and Σ_m 4^{−m}γ_∞(τ_m) = 7/30 (⟺ mean
twisted-collision density 7/10 ⟺ S_∞ = 7/15) — equivalently the stationary C̄_∞(j) and Theorem S's L→∞ limit.
Walk-backs this arc: #25 (P not Σ), #26/#27 (shell convention — killed, backward pinned), #29 (dark = readout not
symmetry); crown ℓ₀-route retracted; deviation candidate 1/(5·21^{k−1}) killed.

## Verification ethos (§8)
Derive-blind / machine-gate / adjudicate; magnitudes pre-registered before the runs that judge them; the
walk-back ledger (~29 logged retractions, public); instrument law near the EP (dense/direct or within-block
power iteration; NO ARPACK/shift-invert; exact-rational; no rate-fitting). Honest negatives ARE the win —
the 27^{−L} law (44× miss), the ±4-seat withdrawal, J1's strong form, the crown ℓ₀-route (2531/4095), the
deviation candidate 1/(5·21^{k−1}), and R4's finite-L edge-density scalings are all in the record.
