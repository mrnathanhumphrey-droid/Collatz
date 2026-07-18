# RESEARCH_ARC.md — the whole C:/Collatz arc, scoped

**Audit date:** 2026-07-17. Purpose: a single top-level map of the *entire* research program — every result (proven / open / killed / superseded), the chronology, the mathematical fields, the canonical spine, and the repo mechanics. Indexed three ways: **by finding · by month · by field of math**. Built from a full read-only audit of all 1,749 tracked files.

> **Read this first — a naming collision.** Three different "Phase N" numberings coexist and are NOT the same thing:
> 1. **LIVE campaign** = Phase 0 → 1 → 2a → 2b → 2c (July 2026, the L3 spectral-gap proof). This is the frontier.
> 2. **Proof scaffolds** = the `PHASE3_*.md` root docs (pen-and-paper briefs for the L3 crux — part of the live campaign).
> 3. **SUPERSEDED (May)** = `PHASE1_DARK_SUBSPACE_RESULT.md`, `PHASE2_APPROX_DARK_RESULT.md`, `PHASE4_DARK_SPECTRAL_GAP_RESULT.md` — the May "dark-subspace / DWM" ancillary, a q=3-only prior effort. Same word "Phase," unrelated campaign.

---

## 0. What this repo is (one paragraph)

A research program on the **statistical and spectral properties of the Syracuse (3x+1) map and its qx+1 family**. It has produced several finished results and is currently one step from a standalone paper. The through-line: the Collatz multiplier `3` sits at a *critical point* of its own family, and the program's job is to prove the exact senses in which that is true. The current deliverable is the **qx+1 universal-rate paper** (`‖π_k‖² ~ C_q·3^{−k}`), whose final step is a spectral-gap theorem (**L3**) reduced to a **limit theorem** `ρ(M_tower,L) → 1/3` at q=3 — the live frontier (Phase 2c).

---

## 1. BY FINDING — the results ledger

### ✅ PROVEN / finished
| Finding | Statement | Home |
|---|---|---|
| **Prefix decomposition** (founding, Paper 2) | The `2^{k−1}` odd residue classes mod `2^k` collapse to `k` conditional σ-distributions via a deterministic prefix ending at `a_final ∈ {3^j}`; count `C(k−1,j−1)`. | `main.tex` (untracked), `notes/writeup.md` |
| **⟨α_det⟩ = log6/log(4/3)** exactly, ∀k | Closed-form class-mean invariant; `α_det(r)=prefix_steps(r)+K_h·log(a_final/2^k)`, `K_h=3/log(4/3)`. | `notes/closed_form_findings.md`, `_paper2_corrections.docx` |
| **Tao-2022 bridge** | `s_mean(r;f) ≈ α_det(r) + K_h·log(N/f(N)) + ε`; slope on α_det ∈ [0.994,1.001] at textbook K_h, no fit. Verified over 40 cells (k∈{8,10,12,14}×5 observables×2 scales), N up to 2³². | `notes/tao_bridge_findings.md` |
| **c = 7/45** (q=3, RIGOROUS UNCONDITIONAL) | Leading coeff of the Syracuse Plancherel mass; `S_∞ = 7/15`. Via R74+R75(Plancherel)+R76(conservation)+R77+R64.B. **Lean-verified** (`__lean_check/`). | `THEOREM_C_745.md`, `notes/c_seven_forty_fifth.md`, `results/result_76_conservation_law.md` |
| **3x+1 ↔ 3x−1 conjugacy** | The two Markov chains are conjugate by negation. | `notes/sibling_3x_minus_1_symmetry_verdict.md` |
| **Syracuse = DWM quantum trajectory** | Davies-Wiseman-Milburn identification, verified to 6 sig figs; dark-subspace classification (D_W dark for j≥2). | `notes/FRAMEWORK_IDENTIFICATION.md`, `R3_DARK_SUBSPACE_STRUCTURAL.md` |
| **Universal rate — at mechanism** | `‖π_k‖² ~ C_q·3^{−k}`, rate q-independent; `3 = 1/Σ_v 4^{−v}` (halving second-moment, q-blind); `D₂ = log3/log q`, so **`D₂=1 ⟺ q=3`**. Survived every adversarial-q falsifier. | `QX1_UNIVERSAL_RATE_WRITEUP_2026_07_14.md`, `results/result_5_universal_rate.md` |
| **Constant** `c̃_q = (q−3)/q` | Same factorization; bracketed leading term is an exact rational identity. | `results/result_4_ctilde_ord2.md`, `notes/c_tilde_structure_verdict.md` |
| **Object identity (R7)** `X_k−X_{k−1}=M_k(1)` | One-line Pythagoras; `M_k(1)=q^k‖d_k‖²` is the clean primitive (carries the rate at every q incl. q=3), `X_k` degenerates to linear at q=3. | `results/result_7_object_identification.md` |
| **THEOREM D1 (toy gap)** | For `M(q,−1,λ)`: `r(λ) = (1−λ²)/(1+λ²)`; maximality via nilpotence of the e=−1 carry graph (LEMMA D1-MAX). Derived, 5/5 vs sweep. | `BRIEF_D1_TOY_GAP.md`, `results/result_phase2b_{s2,F,Dmax}.md` |
| **THEOREM Real-T1** | Real q=3 eigenvalues = twisted autocorrelations `c_k=Σ_δ w_δ²ω^δ`; closed-form left eigenvectors `ℓ_k=ω^{−e_a}R_k(e_ρ)/R_k(0)[γ=0]`. Gate 18/18 @L=3. Closes the **kinematic half** of the q=3 boundary. | `results/result_phase2b_T1.md` |
| **Partner-char (G0)** | The dynamical partner IS `ρ(M_tower)` exactly (Perron of the γ≠0 principal submatrix). | `results/result_phase2b_G0.md` |
| **Corrector chain rung-1 / rung-2** | β*=3/5 (bracket 0.486→9/49, 2.64×); rung-2 trit τ: `r = 4/9 − κ(pop)·W⁻(τ)`, three populations merge (κ·W⁻ pop-independent). Wilson's blind derivation gate-vindicated (5.6e-16). | `results/result_phase2c2.md`, `result_phase2c3_gate.md` |
| **R78 bilinear bound, r≤3** | `|S_partial| ≤ 2√N`, strict, p-uniform. | `notes/PATH2_DISPOSITION.md` |

### 🔶 OPEN (live or standing)
- **L3 spectral-gap crux (the contraction).** Prove `ρ(M_tower,L) → 1/3`, i.e. the corrector chain contracts. Amplitude-decay routes RULED OUT (`|c_k|/c₀→1`; defect mass grows 2.3→9.8; raw mass doesn't discriminate q3/q7). Only surviving mechanism: the depth selection-rule `3^{L−1−j}|k` × cascade-tax `3^{−j}`. **The paper's genuine novel core.** *(Phase 2c, live frontier.)*
- **Universal-rate collision bound.** Off-diagonal collision mass `= O(3^{−k})` on `(Z/q^k)*` (counting problem; forced `C_q≥1` is a theorem, the `O` is conjectured). This is what L3 formalizes.
- **Subdominant rate ε_k** (the c=7/45 *sub*-leading resonance: 0.984 / period-9 / branch-cut). Terminal obstruction after ~75 probes: needs a **rigorous profinite analytic transfer-operator theory** (WATSON+FAURE+BGT triple-PARTIAL all point here). Standing, not on the current critical path.
- **R78 bilinear r≥4** (2×polylog, PARTIAL/walkback); **ε constant −2.35 & ½-slope** closed form; **qx+1 first-passage analog** of α_det.

### ⛔ FALSIFIED / dead ends (honest record — do not retry)
cascade/bridge Φ: F̂_p→μ̂_n (**NO_BRIDGE_FOUND**) · T_3 3×3 companion matrix (spectrum an artifact) · Candidate A (FALSIFIES_F2, no cross-level) · W4 Faure √3 identification (FALSIFICATION_FINAL) · dark-subspace = subdominant-rate mechanism (MOOT, λ→1) · SL₂ embedding (DOESN'T EXIST) · Atkinson/Rota-Baxter operator (only trivial projectors) · Ayyer-Singla / Diaconis-Graham framework (the "+1 breaks group-walk") · Hecke/grössencharacter L-value match (May-31 tangent, nowhere) · most operator-algebra imports (AFL, Belavkin, QSC, amalgamated-freeness, Bruhat-Tits) · θ=e mod3 compression · swap-involution J & rotation S (floor-carry ≠ modular).

### ⚠️ RETRACTIONS (corrected, load-bearing claim survived)
`δ_q ≈ 0.82/ord_q(2)` REFUTED by 2.55×10¹³ → `δ_q = 2^{1−ord}·(q−3)/q` · nilpotency index "=2L" → `≤2L` (true {2,4,5}) · toy "r=0.831 flat" was an ESPRIT mass-artifact → true 0.60 · pre-reg fine survivals {4/27,5/27,4/9,5/9} → {2/9,5/18,4/9,5/9} · "one line to a theorem: generalize R76 conservation" DEAD (conservation ports free, so can't be the missing step; Thm 76.3 leading-mode collapse needs `(q−1)/2=1 ⟺ q=3`).

---

## 2. BY MONTH — the chronology (164 commits, 2026-05-01 → 2026-07-17)

| Month | Commits | What happened |
|---|---|---|
| **2026-05** | 71 | **Founding + broad exploration burst.** Prefix decomposition + Tao bridge (the founding result); c=7/45 thread (R74–R79, Plancherel/Hensel/bilinear); qx+1 generalization + Cramér law; ε_6…ε_16; DWM dark-subspace identification (verified 05-15); Tauberian/Padé arc opened (05-12); late-May Hecke/gchar PARI tangents (05-30/31). |
| **2026-06** | 3 | **Packaging interlude — no new math.** JNT submission bundle; remove copyright PDFs from tracking; add the Lean 4 `__lean_check/` verification project. |
| **2026-07** | 90 | **The L3 spectral-gap campaign + the big reorg.** R6–R44 (object-fixing, pillar-3 rewrite, collision-count identity, Konyagin/Chang/Siegel edges); `r_q` pinned as subdominant eigenvalue; q=3 recast as an **order-2 exceptional point**; the 6-phase falsifier-first campaign (Phase 0→2c) delivering THEOREM D1 + Real-T1, then the Phase-2c limit-theorem / corrector chain (through 07-17, Probe 2c4). The 07-16 "Reorganize into type folders" commit re-committed everything (this is why git mtimes are all July — use commit *subjects*, not mtimes, to date work). |

**The 5 pivots:** (1) ~05-12 empirical → operator-spectral + Tauberian; (2) 05-15 DWM quantum-trajectory reframing + "structural boundary" mapped; (3) June research → publication/Lean; (4) early-July → the standalone qx+1 paper as primary deliverable (old 0.82/ord route killed); (5) mid-late July → the L3 6-phase campaign, q=3 as an exceptional point (the current governing picture).

---

## 3. BY FIELD OF MATH — the taxonomy

The program is genuinely cross-disciplinary. Fields, with representative artifacts and active months:

1. **Number theory — Collatz stopping times & prefix decomposition** (founding). `notes/writeup.md`, `experiments/01_*`. *May.*
2. **Analytic number theory — exponential/character sums & Plancherel.** `THEOREM_C_745.md`, `results/result_76_*`, R78 bilinear, Konyagin/Burgess. *May (R74–79), revived July (collision counts).*
3. **Arithmetic dynamics / transfer operators** — the spine of the L3 campaign. `probe_25_transfer_operator_Aprime.py`, `M_tower`, `probe_R_operator/`. *May → July.*
4. **Dynamical-systems functional analysis — anisotropic Banach / Ruelle resonances** (the mislabeled "`fluid_dynamics/`" corpus: Liverani/Baladi/Tsujii/Faure — **not** actual turbulence). *May, scoping only.*
5. **Spectral theory / exceptional points** — q=3 as an order-2 EP; Jordan blocks; ESPRIT/Krylov. `probe_39_exceptional_point.py`, Phase 2b/2c. *Heavily July.*
6. **p-adic & profinite methods** — profinite transfer operator, Hensel, LTE, 3-adic. `probe_profinite/` (live ρ_slow≈0.83), `zadic_measure_framework.py`. *Late May + July (Chang 2-adic).*
7. **Harmonic analysis / Fourier decay of self-similar measures (Bernoulli-convolution flavor)** — r_q as an L² Fourier-decay rate; Varjú lit. `probe_self_similarity/`, `probe_bernoulli_char_decomp_*`. *May-31 + July.*
8. **Representation theory / Hecke L-functions** — grössencharacters on Q(i), PARI/GP. The four root `.gp` scripts. *End of May (dead tangent).*
9. **Analytic combinatorics / Tauberian theory** — E(z)=Σε_n z^n, Padé/singularity, Flajolet-Sedgewick. ~30 `notes/TAUBERIAN_*`. *05-12 → July (fully disposed).*
10. **Probability / Markov-chain convergence** — stationary S_k, ε_k = S_k−7/15. `probe_epsilon_12..16/`. *May, reused July.*
11. **Open quantum systems / quantum trajectories (DWM)** — Kraus family, dark subspaces. `notes/DWM_*`, `PHASE1/2/4_*DARK*`. *Mid-May.*
12. **Hierarchical Bayesian statistics** — NB2 GLM for σ vs log n, Stan. `model.stan`, `experiments/04_*`, `fits/`. *May, then dormant.*
13. **Formal verification (Lean 4 + mathlib)** — machine-checked c=7/45. `__lean_check/collatz_verify/`. *June.*
14. **qx+1 generalization (cross-cutting)** — the q·x+1 family threading through all the above; the current paper. `QX1_UNIVERSAL_RATE_WRITEUP_*`. *May → July (primary deliverable).*

---

## 4. THE CANONICAL SPINE — what to read, and the replication path

**Three distinct results live here (not one paper):**
- **Paper 1 (finished, publishable):** the prefix decomposition + Tao bridge — `notes/writeup.md`, `main.tex` (currently *untracked*).
- **The c=7/45 result (finished, Lean-verified, q=3):** `THEOREM_C_745.md` + `__lean_check/`.
- **The qx+1 universal-rate paper (current, one step from done):** `QX1_UNIVERSAL_RATE_WRITEUP_2026_07_14.md` + the L3 campaign.

**L3-campaign replication path** (run probes from `probes/`, dense/direct only — iterative eigensolvers FAIL on the defective q=3 Jordan operator; the whole 2c chain is exact-rational):
`probe_25_transfer_operator_Aprime.py` → `probe_phase2a_q2b_q6.py` (the `build_M_gen` builder) + `probe_phase2a_recon.py` → D1 via `probe_phase2b_F.py`/`Fcor.py` → Real-T1 via `probe_phase2b_H.py`→`LALB.py`→`T1.py` → partner via `probe_phase2b_P.py`→`G0.py` → compression `E.py`→`W.py` → defect freeze `F2.py`→`probe_phase2c0.py` → corrector chain `probe_phase2c1→2→3→3_gate→4.py`. Live state: `STATE.md` (line 3 = the 07-17 summary).

**Canonical synthesis docs (authoritative, cite these not the worksheets):** `STATE.md` (live log), `README.md`, and within `notes/`: `logical_chain_findings_to_c745.md`, `findings.md`, `framework_cohesion.md`, the per-probe `*_DISPOSITION.md` capstones.

---

## 5. REPO MECHANICS

- **1,749 tracked files, 180 MB.** ~125k *untracked* files on disk = gitignored literature PDFs + extracted page images + paper-bundle zips + regeneratable numeric outputs (`experiments_output/`, `outputs/`, `fits/`) + visualization renders + local-only author drafts (`_Author_Papers/`, the finished `A_Symbolic_Prefix_Decomposition…` paper, `JNT Submission/`) + a local-only "Hateley framework" branch.
- **Weight is concentrated:** `npz` = **131 MB** of the 180 MB, in **three files** — `probe_self_similarity/pi_15_truncated.npz` (86.8 MB), `pi_14` (28.9 MB), `pi_13` (9.6 MB). All regeneratable (their own gitignore excludes pi_16/17/18 but pi_13/14/15 slipped through tracked). Then tsv 17 MB, csv 12 MB, md 6.8 MB (628 files), py 5.2 MB (565 files), wav 4.1 MB (3 audio renders), 2 DOCX.
- **Tracked by top dir:** notes 452, probes 403, outputs 187, logs 180, results 120, experiments 112, root 25, + ~30 scattered May side-thread dirs (~340 files).

### Audit flags (replication + clarity blockers)
1. **Heavy tracked binaries** — 3 npz files ≈ 125 MB are regeneratable; the single biggest bloat.
2. **Thin, mislocated dependency manifest** — only `logs/requirements.txt` (4 lines: numpy/numba/polars/matplotlib); no scipy/mpmath/sympy pin despite heavy use; no `pyproject.toml`/`environment.yml`.
3. **Hardcoded absolute paths** — **301 of 565** `.py` files contain literal `C:\`/`D:\` paths; won't run off this machine.
4. **The "Phase N" naming collision** (see top warning) — three unrelated campaigns share the word.
5. **~300 of 452 notes/ files are subordinate worksheets** collapsible to ~40 `*_DISPOSITION.md` capstones (nothing should be *deleted* — it's the honest falsification record — but it can be collapsed for presentation).
6. **The finished Paper-1 (`main.tex`) is untracked** — the replication surface doesn't contain the project's most publishable result.
7. **Scattered May side-dirs** — most are dead tangents (Hecke/gchar, Ayyer-Singla, Atkinson, sonification `audio_data/`, one-off `probe_*` folders); a handful are live (`experiments/`, `data/`, `probe_profinite/`, `__lean_check/`, `probe_self_similarity/`, `inverse_tree/`).

---

_This file is a navigational map, not a result. All statements are sourced from a 2026-07-17 read-only audit of the tracked tree; findings' authoritative statements live in the files cited. For the live frontier see `STATE.md`; for the founding result see `notes/writeup.md`; for the current paper see `QX1_UNIVERSAL_RATE_WRITEUP_2026_07_14.md`._
