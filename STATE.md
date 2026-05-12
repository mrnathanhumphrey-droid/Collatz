# Live State — Collatz framework synthesis

**Last updated:** 2026-05-12 (**seven-probe spectral trajectory mapped + Tauberian framework arc opened**; c=7/45 Nisoli closure recharacterized as **structurally inapplicable** — not just unresolved; **cross-frequency closure on V_M = span{M_n^{ab}(g, c)} found positive** but **V_M does not close under iteration** (phase + parity obstructions); **Tauberian scoping landed H_AMBIGUOUS** with Flajolet-Sedgewick Ch. VI / Chevalier 2507.15394 Thm 1.16 as candidate framework, single-theorem selection pending ε_7 exact-rational compute. **Post-consolidation post-mortem:** Delta diagnostic + Padé extension probes both confirm n=2..6 window is DOUBLY pre-asymptotic — combining with prior-session numerical ε_7..ε_13 measurements showing 4.7× envelope jump at n=7 + slow oscillating mode at ρ ≈ 0.984, **the true leading singularity may not be at z=2 but at z ≈ 1.016 (or complex-conjugate pair)**; R77.6's branch-cut at z=2 may be SECONDARY structure visible only in fast-transient window. **R77.7 v2 new solver fired** (modular CRT + rational reconstruction, targeting <2hr at k=7 vs original 8.5hr-killed Fraction Gauss elimination); **Padé numerical extension fired** using numerical ε_7..ε_13 to test the singularity shift hypothesis this session.

**Both landed:** Padé numerical extension → H_TWO_SINGULARITIES_VISIBLE (z=2 REFUTED as leading singularity; Hadamard radius at n=10..13: 2.06 → 1.81 → 1.66 → 1.57, monotone inward, complex-conjugate-pair structure consistent with period-9 sign pattern). R77.7 v2 → **ε_7 EXACT-RATIONAL COMPUTED in 39 min, 151 primes, 0 reconstruction failures, witness-verified. ε_7 = -1.175236830374320×10⁻³ (Fraction ~4485 bits); |ε_7|·2^7 = 0.1504 confirms 4.7× envelope jump at exact precision; R76 §11's "(1/30)·(1/2)^n + O((1/4)^n)" REFUTED algebraically.** ~13× speedup vs original R77.7 killed.

**Framework-reopening map:** T_lead Nisoli bypass (R77.3) is LIVE again — cross-freq derivation today provides the machinery, just need to retarget Off_lin's contribution to corrected rate (ρ ≈ 0.984) instead of 1/2. Joint 2-3-adic Bohr conditionally live at k-iteration timescale (retirement was at v-bracket scale, different timescale). R58 Esscher and R77.4 K_k stay dead. **Tao 2020 measure-theoretic framework stays the strategic anchor; the broader (q-sweep, σ-records, trajectory measure, bilinear bound, F̂_p) are unchanged.**). Earlier 2026-05-06 / 2026-05-05 entries retained below for continuity.

**2026-05-12 entry: seven-probe spectral boundary + Tauberian shift.** Continuation session after 2026-05-11's Path 2 + Hensel bilinear bound delivery. Trajectory: (1) Reading A scoping → Candidate A (locally constant functions on Ẑ_3^×) tractable; (2) **Candidate A construction landed H_CANDIDATE_A_FALSIFIES_F2** — c_{n, k} ∈ Q computed at n=1..6, 15 of 21 are exactly 0/1; only diagonal k=n−1 nonzero; c_{n, n−1} = S_n → 7/15 (does NOT decay); structural diagnosis: K_n(d) supported on d ≡ 0 mod 3^{n−1} forces φ_n ∈ W_{n−1} **by construction** — rate-1/2 does NOT live in W_k filtration; (3) R76 §11 verification landed INCONCLUSIVE — (1,4)/eigenvalue-1/2 is §11's own "structural conjecture / open" with T_diag = (1/5)·[[1,1],[4,4]] having spectrum {0, 1} (eigenvector (1, 4) at λ=**1**, not 1/2); (4) T_N construction landed H_OFF_LIN_UNDERSPECIFIED — R77 sketch §5 articulates Off_lin as procedure-plus-claim, not 2x2 matrix; cross-frequency closure asserted not derived; (5) **Cross-frequency derivation landed H_CROSS_CLOSES_ON_ENLARGED_SPAN** — first positive structural advance after five negatives: closure exists on V_M = span{M_n^{ab}(g, c)} parameterized by g = v' − v; g=0 slice IS span{P_n^{ab}(c)} (T_diag rigorous); g≥2 genuinely new; **mixed-parity (v_3(d)=0) vanishing upgrades R76 §11's empirical P^{+−}=0 to rigorous algebraic identity** via lift-fiber orthogonality; (1,4)-direction preservation falls out for all g via W_−(g)/W_+(g)=4 uniform in g; **R77 sketch §5's "quadratic forms in {P_n^{ab}(c)}" assertion is FALSE as stated** (correct statement on V_M); empirical verification at n=2,3 confirms augmented rank 6/7 vs P-only rank 1; (6) **T_V spectrum probe landed H_M_RECURSION_UNDERSPECIFIED** — V_M does NOT close under iteration n→n+1; two new obstructions: F1 phase offset θ_{v,g}=2^v·ẽ_g/3 generically not expressible as ẽ_{G''}−ẽ_G (worked example g=2,v=2,v'=3: phase 5/8 has no integer G); F2 parity: incoming g∈{2,4} produces only ODD outgoing G, but V_M = span{g∈{0,2,4,...}} doesn't contain odd-G moments; g=0 closure holds (T_diag); empirical t_v_compute.py confirms both; (7) **Tauberian scoping probe landed H_AMBIGUOUS / INCONCLUSIVE** — framework right (Flajolet-Sedgewick Ch. VI singularity analysis), single-theorem selection pending ε_7 data; Chevalier 2507.15394 Thm 1.14 (pure √-singularity) FALSIFIED at leading order (predicted n^{-3/2} factor gives growing product 0.108 → 0.468 across n=2..6, not constant); Newman-Zagier excluded (Dirichlet-series, not power-series, requires pole R77.6 rules out); Chevalier Thm 1.16 (meromorphic h with pole of order M at 0 → n^{M − 3/2}) is cleanest single-theorem candidate but M not determined at N=5; **subleading δ_n := |ε_n|·2^n − 1/30 non-monotone, changes sign between n=5 and n=6** — inconsistent with single-term ansatz, possibly multiple competing subleading terms / second-sheet secondary singularities; reconciliation: leading |ε_n|·2^n ≈ 1/30 simple-pole-like + subleading branch-cut consistent with R77.6 Padé pattern (Padé picks up branch part) and R77.4 envelope (Jordan ruled out, log/power tied). **Main-thread digs:** (a) (1,4)-direction structurally forced by R64.B's squared class mass (1/3)² : (2/3)² = 1:4; T_diag's λ=1 eigenvector IS (1, 4) — threads "what (1, 4) is" and "what T_diag conserves" collapse into single finding; sharp structural target Off_lin · (1, 4) = scalar −1/2 for the conjectured eigenvalue 1/2; (b) R77.6 thorough re-read confirms branch-cut at z=2 multi-mode consistent with prior-session ε_7..ε_13 oscillation (fast transient k=2..6 rate ~1/2 + slow oscillating mode k=7+ with ρ≈0.984); (c) Bohr probe history confirmed retired (STATE item 11 — descent-funnel artifact at v≤100, not joint structure at large v); R77.7 status file's "Bohr supersedes" framing itself outdated. **c=7/45 closure landscape recharacterized:** three-obstruction landscape from 2026-05-11 ("Tao C_A INFEASIBLE / Bilinear |K| DELIVERED / Spectral M_3 NOT RESOLVED") becomes **two-plus-structural-boundary**: spectral M_3 is now **STRUCTURALLY INAPPLICABLE** — Nisoli requires discrete eigenvalue of resolvent at rate 1/2, seven probes + R77.6 + R77.4 + ε_k oscillation all point to branch-cut endpoint / multi-mode, no discrete eigenvalue exists. **Tauberian framework arc** is the live direction — Flajolet-Sedgewick Ch. VI is the right abstraction, ε_7 exact-rational compute (R77.7 re-fire; previously killed at 8.5hr) gates single-theorem selection. **Literature bundle assembled** at `C:/Users/Nate/Documents/burgess/literature/` — 73 PDFs across 7 math-field lots (tauberian/burgess_bound/bilinear_character_sums/multiplicative_subgroups/heilbronn_sums/taos_machinery/lagarias_collatz); master math-field INDEX.md at `burgess/literature/INDEX.md`. **Updated artifacts:** `SESSION_DISPOSITIONS_2026_05_12.md`, `README.md`, memory file `project_collatz_r78_bilinear_cracked.md`, `burgess.zip`. **All deliverables this session:** READING_A_SCOPING_*, CANDIDATE_A_*, R76_S11_*, T_N_*, CROSS_FREQ_*, T_V_*, TAUBERIAN_SCOPING_* + verification scripts (`candidate_a_compute.py`, `cross_freq_compute.py`, `t_v_compute.py`, `tauberian_verify.py`).

(legacy entry from 2026-05-06 follows below for continuity)

**Earlier last-updated:** 2026-05-06 (**ε_13 measured = +2.948e-3 — order-3 recurrence on ε_k WALKED BACK as window-unstable**; Δ_k entropy-deficit Outcome B/C confirmed at k=12→13 = 0.879 with widening gap to ρ_slow → distinct modes; ρ_slow ≈ 0.83 reliably identified ONLY as L¹/TV inverse-limit rate from probe_profinite, NOT from any order-N recurrence fit on ε_k). Earlier 2026-05-05 (Bohr signal RETIRED + R58 closure attempt #2 failed + two-mode ε_k decomposition FALSIFIED at k=6 + **rate-1/2 envelope claim FALSIFIED at k=7, non-monotone S_k** + **q-spectrum probe: \|λ_2\|^(q) is q-universally near-zero, not 1/2 — R77.4 erratum extended to q ∈ {5,7,11,13}** + **ε_10 measurement (k=8,9,10 extension) confirms oscillation in ε_k beyond ambiguity, slow-mode envelope rate ρ ≈ 0.984, period ≈ 9.2 in k-space, rate-1/2 \|ε_n\|·2^n bound empirically refuted in extrapolation** + **low-v_2 residue admissibility probe returned 18/18 cells null — bit-budget contradiction mechanism ruled out**).

**2026-05-06 entry: ε_13 + recurrence walk-back.** ε_13 = +2.948e-3 measured via FFT on cached `pi_13_truncated` (truncation error sub-machine-precision). Continues the post-zero-crossing rising sequence |ε_10..13| = 0.72, 1.50, 2.27, 2.95 × 10⁻³ with decelerating ratios 2.08 → 1.51 → 1.30 (approaching a peak; another sign flip plausibly imminent). **Order-3 recurrence on ε_k has been walked back as not a structural feature**: refitting on ε_2..ε_12 gives dominant root +1.030 (not 0.827); on ε_2..ε_13 gives +1.115. Order-3 prediction for ε_13 is 9.2% off measured. The "ρ_slow = 0.826934 from order-3 recurrence" entry in `result_renormalization_recurrence_fits.csv` was fitted on ε_2..ε_10 in an earlier window before the post-zero-crossing growth (ε_10..13 all positive and growing) destabilized the fit. **The reliable identification of ρ_slow ≈ 0.83 is the L¹/TV inverse-limit rate from probe_profinite (R² = 0.97), not the recurrence root.** Documents updated: `framework_cohesion.md`, `probe_profinite/profinite_findings.md`, `probe_self_similarity/self_similarity_findings.md`. New: `probe_epsilon_13/epsilon_13_findings.md`. **Δ_k entropy-deficit Outcome confirmed at k=12→13 = 0.879** (continues monotone rise; OLS gap to ρ_slow widens from 2.7% to 5.8% as fit window extends — distinct-modes signature). Candidate analytic value for ρ_Δ: 7/8 = 0.875 (suggestive, no derivation).

(legacy entry from 2026-05-05 follows below for continuity)
 Today's update: the joint 2-3-adic Bohr empirical positive cited yesterday has been deflated — bracket-stratification probe ([result_bohr_probe_strat.md](result_bohr_probe_strat.md)) shows per-bracket chi²/df at v > 10⁶ is 0.94–0.95 (z ≈ −1, statistically CRT-independent within ±2σ); the original aggregate z=16.5 was the v ≤ 100 descent funnel, not joint structure at scales relevant to D_emp. Also: second R58 → D_emp closure attempt via Esscher tilt at log(R58/R60) per R69's mechanism FAILED (Pearson 0.857 → 0.867, only +0.010 improvement; r=5/13/23 improved but r=1/r=21 broke). **Two R58 Esscher attempts now both rejected** — closure requires non-uniform tilt or a sign-aligned observable; uniform Esscher cannot fix opposite-sign residuals. **Convergence-shape ε_k = A·(1/2)^k + B·(1/3)^k two-mode hit (logged earlier today) FALSIFIED at k=6** ([result_epsilon_6.md](result_epsilon_6.md)): float64 power iteration on K_6 (486 states) gives ε_6 = −4.98×10⁻⁴ vs predicted −5.86×10⁻³; off by ~10×. **Rate-1/2 envelope FALSIFIED at k=7** ([result_epsilon_7.md](result_epsilon_7.md)): power iteration on K_7 (1458 states) gives ε_7 = −1.18×10⁻³, |ε_7/ε_6| = **2.36** — ratio reversed up, sequence non-monotone. S_k has local maximum at k=6 (S_6 = 0.46617 closest to 7/15) then backs off at k=7 (S_7 = 0.46549). |ε_n|·2^n envelope: 0.038, 0.041, 0.039, 0.037, 0.032, **0.150** at n=2..7 — supposed-stable-near-0.04 was visible portion of longer-period oscillation, NOT a true rate-1/2 envelope. Cross-validated: power-iteration and scipy.eigs agree to 1e-15 at k=6 and k=7; at k=5 float matches exact rational to 1e-15. **Asymptotic-rate question is wide open.** Plausible candidates (NO REFIT per brief): damped oscillation with complex-conjugate rates ρ·e^±iθ, longer-period structure, non-elementary shape. **ε_10 measurement (2026-05-05) RESOLVES the candidate ambiguity:** ε_10 = +7.21×10⁻⁴ confirms oscillation in ε_k trajectory beyond ambiguity — ε_10 sign-flipped from k=9 (first positive since k=2), magnitude rebounded to ~|ε_8| (|ε_10/ε_9| = 95.8). Two-mode model fit: fast transient mode dominating k=2..6 (rate ~1/2), slow oscillating mode emerging k=7+ with envelope rate ρ ≈ 0.984 per k-step, period ≈ 9.2 in k-space. Three methods (power iter, scipy.eigs, FFT) agree to 4×10⁻¹⁵ at k=10. The slow-mode rate kills the |ε_n|·2^n bounded-envelope reading: at ρ ≈ 0.984, |ε_n|·2^n grows as (1.968)^n. **Rate-1/2 conjecture as |ε_n|·2^n bound is empirically refuted in extrapolation.** Within-level K_k spectral gap (~0.998) is consistent — within-level mixing is fast; the slow oscillating mode must live in **inter-level renormalization**. Confirms R77.4 erratum and q-spectrum probe finding: convergence rate is a tower-renormalization phenomenon, not a within-level operator phenomenon. Predictions for ε_11–15 documented in [result_epsilon_10.md](result_epsilon_10.md) for falsification testing. **Low-v_2 admissibility probe (2026-05-05) returned 18/18 cells null** ([result_low_v2_residue_admissibility.md](result_low_v2_residue_admissibility.md)). Bit-budget contradiction mechanism ruled out: low-v_2 trajectories show LARGER admissibility horizons (mean h ≈ m/E[v_2_low] with E[v_2_low] ≈ 1.73), not smaller. Bit-budget consumption is deterministic at rate Σ v_i per step; slower consumption produces longer admissibility, not contradiction. Ensemble structural numbers consistent with bit-budget theory: mean v_2 = 2.10, mean horizon ≈ m/2, all 3.3M trajectories converged within 256 steps. GPT's "contradiction-mod-2^m" framing not supported at simplest interpretation. Reframings (joint admissibility, residue counting, stuck residues) flagged as potential future directions but not pursued — each requires substantial analytical work and drifts toward closure-attempt territory. **Yesterday's landmarks (still standing):** R77.3 falsified the 3-mode geometric ansatz for ε_n over Q; R77.4 envelope fits gave (M) — Jordan ruled out, log/power tied at N=5; R77.4 erratum showed K_k itself has no eigenvalue near 1/2 (not the rate operator); R77.6 generating-function probe found branch-cut signature at z=2 (type indeterminate); R77.7 (k=7 ε extension) killed at ~8.5 hr (cache retained for future re-fire; the original "superseded by Bohr" framing is now also retired). Path-C subroutes all closed (R79b: ‖S_partial‖ saturated β=0.522; saddle-class: structural but not closure path; C2/BGK: random-like, multiplicative energy ≈ N³; C3/band-l¹: ‖ĥ‖_{ℓ¹(D)} saturates at trivial bound; band-spectral: lf_mass→0.25, smooth-completion weakened). **Q-sweep test 2:** literal hypothesis c_q = S_∞^{(q)}/q falsified for q ≥ 5 (S_k diverges); renormalized c̃_q := lim S_k^{(q)}/(q/3)^k exists universally; **c̃_q = (q−3)/q confirmed at q=11, 13, 17 within 1%** (q=17 ruled out non-prim-root explanation for q=7's anomaly). **Sibling study (3x±1):** forward S_n^{3x−1} = S_n^{3x+1} via K_- = σK_+σ chain symmetry (proved at k=1..4); inverse-tree D_n differs structurally (basin fingerprint after matched-N analysis: factor ~0.2-4 residual, ~95% of raw 10⁴× difference is sample-size). **Archive:** `closed_form_findings.md` (79+ results). **Latest results:** `c_seven_forty_fifth.md` (R75), `result_76_conservation_law.md` (R76), `result_77_T_lead_spectrum.md` (R77), `result_77_3_nisoli_bypass.md` (R77.3), `result_77_4_operator_shape.md` + `result_77_4_K_spectrum_erratum.md` (R77.4), `result_77_6_generating_function.md` (R77.6), `result_77_7_status.md` (R77.7 NOT COMPLETED), `result_78.md` (R78), `result_79.md` (R79), `r79b_S_partial_empirical.md` (R79b), `saddle_class_subsum_analysis.md`, `band_l1_analysis.md` (C3), `bk_moments_analysis.md` (C2), `band_spectral_decomposition.md`, `result_q_sweep_test_2_c_q.md`, `c_tilde_structure_verdict.md` + `c_tilde_q17_probe.py`, `sibling_3x_minus_1_symmetry_verdict.md`, `duality_S_vs_D_verdict.md` + `duality_followup_verdict.md`, `result_bohr_probe.md` + verification chain (RETIRED — bracket strat deflated), `esscher_tilt_r58_closure_v2_verdict.md` (second failed R58 closure), `result_epsilon_8.md` + `result_epsilon_9.md` + `result_epsilon_10.md` + `result_epsilon_11.md` (ε_k oscillation confirmed beyond ambiguity; ε_11 = +1.50e-3 grew further), `result_renormalization_spectrum.md` (recurrence-fit slow-mode rate ρ_slow ≈ 0.83 real, R_k operator on W_k spectrum is 0 — W_k is forcing-only, not propagating), `result_low_v2_residue_admissibility.md` (bit-budget contradiction mechanism null at 18/18 cells).

---

## Active framework synthesis

**The trajectory measure on Z_2** (= the survivor-conditioned residue distribution D_avg(r) for r mod 32) has the following multi-layer characterization:

1. **Integer-level identification** (R58): D_avg = mod-32 marginal of the inverse Collatz tree from m=1 weighted by subtree-size, value-truncated at N=2^22. Pearson 0.86, MAE 0.118 in mean-1 units. Stable across N = 2^16 to 2^22.

2. **3-adic Fourier closure** (R63 → R66 → R70 → R73/R74):
   - |μ̂(1/3)|² has closed form = ½[(a−b)² + (b−c)² + (a−c)²] / (a+b+c)² with (a,b,c) = mass-fractions at residues (0,1,2) mod 3
   - First-principles: (a, b, c) = (1, D+2, 2D−3)/(3D) where D = ⟨inverse-tree path length⟩; → (0, 1/3, 2/3) as D → ∞ (R64.B)
   - At higher 3-adic levels: |μ̂(a/3^k)|² has closed form via Markov chain on (Z/3^k Z)* with v ~ Geom(½) heuristic (R66)
   - Asymptotic decay: average |μ̂|² over primitive a → S_∞ / (2·3^(k-1)) (R66)
   - **S_∞ = 7/15** strongly evidenced (10⁻⁴ extrapolation, R70). Equivalent algebraic identities: S_{k+1} = 3·X_k·⟨ψ−1/3⟩_w (R73) = 3^(k+1)·||d_{k+1}||² (R74). Both give same 7/15 limit via different decompositions.

3. **Structural class** (R59, R62):
   - Multifractal Z_2 measure with wide D_q spectrum (D_0 = 1.00, D_∞ ≈ 0.15)
   - Fourier dimension σ = 0 (atomic-class spectrum, NOT Sullivan-conformal)
   - Resonances at 3-adic rationals reflect (3m+1) arithmetic structure
   - Right literature home: multiplicative number theory measures on Z_p (NOT Bernoulli convolutions or Erdős lacunary)

4. **Mechanism** (R45-R47, R64.B):
   - r ≡ 0 mod 3 are LEAVES in inverse tree (no Syracuse predecessors)
   - r ≡ 2 mod 3 has smaller smallest-pred (≈2m/3) than r ≡ 1 (≈4m/3) → 2× more mass on r ≡ 2
   - Forward Syracuse v-parity → next-residue rule: v even → 1 mod 3, v odd → 2 mod 3 (verified 100% on 1M pairs)

---

## Closed forms LOCKED

| Quantity | Value | Source |
|---|---|---|
| ⟨α_det⟩ | log(6)/log(4/3) | R1 |
| Asymptotic (a, b, c) mass fractions mod 3 | (0, 1/3, 2/3) | R64.B |
| \|μ̂(1/3)\|² (D → ∞) | 1/3 | R64.B |
| \|μ̂(1/3)\|² (finite D) | (D²−4D+7)/(3D²) | R64.B |
| \|μ̂(1/2)\|² | 1 | R63 (trivial, all m odd) |
| ψ at k=1 (sub-cell purity, lifting k=1→2) | 3/7 | R70/R73 |
| S_1 = primitive Fourier sum at level 1 | 2/3 | R70 |
| S_2 | 10/21 | R70 |
| S_∞ = invariant primitive Fourier sum | 7/15 | R70 (extrapolation, rigorous proof open) |
| Asymptotic invariant X_k·⟨ψ−1/3⟩_w | 7/45 | R73 |
| ⟨\|μ̂(a/3^k)\|²⟩_a (avg over primitive) | 7/30 · 3^(−(k−1)) | R66+R70 |
| ‖d_{k+1}‖² leading rate (sets 7/15 constant) | 1/3 per level | R74 |
| \|S_{k+1} − 7/15\| subleading rate (convergence) | 1/2 per level | R73 |
| Decomposition: ‖d‖² = (7/45)·(1/3)^k + ε_k/3^(k+1) | algebraic identity | reconciliation script |
| P(V=k \| r=21 cylinder, uniform m) | shifted Geom(½): 2^(−(k−4)) | R52.A |
| **Universal qx+1 ratio** S_{k+1}^{(q)} / S_k^{(q)} → q/3 | empirical (4-sig-fig at q=11,13 by k=2) | q-sweep test 2 |
| **c̃_q := lim S_k^{(q)} / (q/3)^k exists** ∀ tested q (no closed form yet) | (q,c̃_q) ≈ (3, 7/15), (5, 0.488), (7, 0.78), (11, 0.7288), (13, 0.7698) | q-sweep test 2 |

---

## Major supersessions

| Old | New | Reason |
|---|---|---|
| R52.B (inverse-tree miss, Family C) | **R58** (value-truncation gives Pearson 0.86) | Wrong truncation regime; depth-50 tree concentrates pathologically on m_j chain. Value-truncation matches D_avg's integer-uniform sampling. |
| R57 H-dim coincidence 2·log(λ_max)/log(2) ≈ 0.68 | R61.A walked back | Walk-back: the comparison was against a heuristic Chang value; multifractal analysis (R61.B) shows wide spectrum, single-δ doesn't apply. |
| R71 conjecture λ_2(K_k) = 1/2 | R71.B + R73 | K_k has rank 2 with λ_2 = 0; convergence rate from level-lifting structure, NOT chain spectrum. |
| R73 max\|d\| → 0 at rate 1/2 | R73 (revised) + R74 | max\|d\| only decays at rate 0.97/level. Right invariants: X_k·⟨ψ−1/3⟩_w (rate 1/2, R73) or ‖d‖² (rate 1/3, R74). Both give same S_∞ = 7/15. |
| Apparent R73 vs R74 inconsistency (1/2 vs 1/3) | reconciliation script (2026-05-03) | NOT contradictory: ‖d‖² = (7/45)·(1/3)^k + ε_k/3^(k+1) where ε_k decays at 1/2. R74 leading rate, R73 subleading rate, both empirically present. |
| R60 "size-stratified Markov derives D_avg from first principles" | R77 (γ) | R60's empirical K = K_dynamics × W_visit. Derived K (no visit-weights) recovers local dynamics (per-state Pearson 0.96) but Perron eigvec is uniform → D_pred ≈ 1.0, FAILS to recover D_avg. R60's identification works because empirical K bakes in visit weights. Reframe as "Perron of empirical kernel" not "first-principles derivation." |
| "Three independent confirmations of dim/D_avg" | independence_audit (Validation Task 2) | 0 of 7 cross-result claims are independent confirmations; 6 are consistent characterizations of same Collatz dynamics; 1 is propose-not-demonstrate. R58↔R60 cross-Pearson +0.92 (same identification two ways). |
| R23 / Chang dim sharing | dim_h_validation (γ) | Five distinct values in 0.07-wide window; no algebraic identity. Furstenberg branching 0.338 vs Chang exact log(φ)/log(2) = 0.694 vs spatial info dim D_1 = 0.608. |
| Lagarias-Sinai v ~ Geom(½) exact | R68 (γ) | Marginal moments E[v], P(v=1) hold within 1%; full distribution has 5-25% structural deviations at specific j. Affects K_h precision (~0.5%) and downstream closed forms. |
| Chang ↔ trajectory operator factorization | R69 (γ) | Distinct dynamical observables, no algebraic relation. |
| R65 conjecture 4^(−(k−1)) decay | R66 | Asymptotic factor is 3, not 4. Conjecture was finite-k artifact at k=1→2. |
| R77.2 conjecture: T_3 has spec {1/2, 1/4, 1/8} over Q (Nisoli bypass) | R77.3 (β) | 3-mode geometric ansatz fitted from {ε_1,ε_2,ε_3} over Q gives A = −157462/3058335 ≠ −1/30; predictions miss ε_4..ε_6 by 28-41% relative. Bypass FAILS; Tao Prop 1.17 effective C_A still required. |
| R77 implicit framing: "T's spectrum" shapes ε_n envelope | R77.4 erratum + q-spectrum probe (item 14) | K_k itself has \|λ_2\| ≈ 10⁻⁵..10⁻³ at k=3..6, NO eigenvalue in [0.3, 0.7]. K mixes in O(1) steps. The (1/2)^n envelope is INTER-LEVEL renormalization, not within-level spectral gap. **Now q-universal:** same spectral triviality holds at q ∈ {5, 7, 11, 13} per [result_qspectrum.md](result_qspectrum.md) — \|λ_2\|^(q) ∈ [10⁻⁶, 10⁻⁵] across all primes tested, never near 1/2. |
| Hypothesis: 7/45 is the q=3 instance of a closed-form c_q family | q-sweep test 2 (NO-PATTERN literal / UNIVERSAL-SHAPE renormalized) | S_∞^{(q)} doesn't exist for q ≥ 5 (sequence grows like (q/3)^k). 7/45 is q=3-specific. But renormalized c̃_q does exist universally — partial finding. |
| R63 {m_j} atomic decomposition for resonance | R63 (revised, full-population partition) | {m_j} chain accounts for only 0.15% of |μ̂(1/3)|²; resonance comes from full-population mod-3 mass asymmetry. |
| Sullivan-conformal measure framing (R59) | R62 | Multifractal with σ = 0 — not in Sullivan/Pollicott-Urbański constant-δ machinery. |
| Bernoulli convolution / Erdős-class lacunary framing (R62) | R63 | Resonance is population-level, not chain-level. Right home: multiplicative number theory on Z_p. |

---

## Open pieces (active path)

1. **Esscher-tilt closure for R58 → D_emp gap — TWO ATTEMPTS BOTH FAILED.**
   - Attempt #1 ([esscher_tilt_r58_closure.md](esscher_tilt_r58_closure.md), σ_orbit observable): best λ ≈ −0.01, +0.014 Pearson improvement. Rejected.
   - Attempt #2 ([esscher_tilt_r58_closure_v2_verdict.md](esscher_tilt_r58_closure_v2_verdict.md), log(R58/R60) observable per R69's mechanism, 2026-05-05): best λ = 0.95, +0.010 Pearson improvement (0.857 → 0.867). r=5/13/23 improved but r=1/r=21 broke. Rejected.
   - **Structural takeaway:** a uniform Esscher tilt observable cannot carry per-residue sign information when residuals at different residues demand opposite-direction corrections (r=5 needs +, r=13 needs −). Both attempts share this failure mode. R69's weight-ratio mechanism is partial and operationally insufficient.
   - Closure path requires either (a) non-uniform tilt with residue-conditional λ_r, (b) a sign-aligned observable across residues (neither σ_orbit nor log(R58/R60) qualifies), or (c) a closure mechanism outside the Esscher-tilt family. R58 stays at Pearson 0.857 as the best the inverse-tree subtree-size measure can achieve.

2. **Rigorous proof S_∞ = 7/15.** PARTIAL: c = 7/45 now algebraically anchored (R75/R76); rate-½ proof remains.
   - **Plancherel formula (R75, RIGOROUS):** S_k = Σ_{ξ ∈ Z/3^k, 3∤ξ} |μ̂_k(ξ)|². So c = (1/3)·lim Σ |μ̂|² over high-freq.
   - **Conservation law (R76, RIGOROUS):** Σ_{j=0,1,2} M_{n+1}(η_0 + j·3^n) = 0 where M_n(η) := Σ_ξ μ̂_n(ξ) μ̂_n*(ξη).
   - **Leading-mode identity (R76, RIGOROUS):** S_{n+1} = −2·M_{n+1}(1+3^n). Reduces rate question to scalar sequence R_n := M_n(1+3^{n−1}) → −7/30.
   - ~~**Empirical rate ½ verified through k=5** (|ε_n|·2^n stable at C ≈ 0.04 for n=2..5).~~ **WALKED BACK 2026-05-05 by ε_7 measurement** (see R77.x bullet below): |ε_n|·2^n envelope is non-stable — at n=7 it jumps to 0.150, ~4× the supposed plateau. The k=2..5 stability was a finite-window artifact of a longer-period non-monotone oscillation. Asymptotic rate is unknown; ratio-based extrapolation requires k≥8 before any refit.
   - ~~**Provisional certified bound (R75):** assuming rate ½, |c − S_k/3| ≤ 0.013·(1/2)^k. At k=5, bound 4.2×10⁻⁴, actual 3.8×10⁻⁴.~~ **CERTIFICATE NOT VALID 2026-05-05** under the new k=7 data: rate-½ assumption is the load-bearing hypothesis and it's now empirically contradicted. The bound at k=5 happened to hold by coincidence of finite-k phase; at k=7 the actual deviation 1.18×10⁻³ exceeds 0.013·(1/2)^7 = 1.0×10⁻⁴ by an order of magnitude. Bound and proof framework need either (a) a different rate hypothesis, (b) a non-monotone-tolerant bound, or (c) abandonment until structural form is identified.
   - **Comprehensive obstruction map** for the analytical closure (Kalafatelis eq 190) — six probes, all closed/weakened:
     - Cochrane Theorem 2 (R78): D = 0 obstruction sharp; trivial-bound only. ❌ closed.
     - Pólya-Vinogradov (R78): worse than trivial for r ≥ 3. ❌ closed.
     - van der Corput B=1, B=2 (R79): rigorous rate 0.73 / 0.81 — sub-trivial but well above empirical √N ≈ 0.5. ⊳ insufficient.
     - Empirical |S_partial(r)| sweep (R79b): direct measurement at r=8..20 gives β = 0.522 ± 0.008 (R² = 1.0000); ĥ Plancherel-saturated, no concentration to exploit. ⊳ confirms vdC stalls below empirical.
     - Saddle-class subsum (saddle_class_subsum_analysis.md): β_j ≈ 0.92, 1.06, 0.98 across j=0,1,2; partition is structural but each subsum scales linearly in n_j → not a closure path.
     - C2 / BGK on ⟨4⟩ (bk_moments_analysis.md): M_4 slope = 3.0059 (R²=0.9999), random-like multiplicative arc. BGK applies to primal Kalafatelis sum directly. ⊳ partial closure of rate-1/2; eq 190 still requires more.
     - C3 / direct band-l¹ (band_l1_analysis.md): N_r^{-1/2} ‖ĥ‖_{ℓ¹(D_{r,t}(η))} grows exactly N_r^{1.0} across 36 (ℓ,t,η) cells. Saturates trivial bound; no inter-m cancellation. ❌ closed.
     - Band-spectral decomposition (band_spectral_decomposition.md): lf_mass → 0.25 (uniform-in-k baseline) at r=10..14; smooth-weight cancellation gives only constant factor. Smooth completion (R78 path 2) empirically WEAKENED. ⊳ closed for smooth-weight purposes.
     - **Remaining un-attacked:** C1 (5x+1 sibling-attack reframing — multi-day rebuild). Other Path-C routes exhausted.
   - **R77.x operator-shape probes** (per-strand status):
     - R77.3 (β): finite-mode geometric ansatz for ε_n over Q FALSIFIED. Cleanly rules out simple-rational spectrum.
     - R77.4 (M): envelope curve fits at N=5 (n=2..6) inconclusive between Jordan / log / power-law; Jordan ruled out by direction (b<0); H2 ≈ H3 tied at ΔAIC = 0.23.
     - R77.4 erratum: K_k spectrum has nothing near 1/2 — the rate operator is INTER-LEVEL renormalization, not within-level mixing.
     - R77.6 (G-branch-cut, type indeterminate): generating function Σ ε_n z^n's Padé approximants place poles in [2.05, 2.35] real-axis; diagonal [n/n] converges 2.076 → 2.051 from above; consistent with branch cut at z=2, NOT simple pole. Power vs log unresolved at N=5.
     - **R77.7 (NOT COMPLETED):** k=7 Markov chain extension killed at ~8.5 hr by user direction. Original framing claimed "substantially superseded by joint 2-3-adic Bohr empirical positive" — that successor framing is **also retired** (Bohr signal deflated 2026-05-05). R77.7 now stands as: killed before completion, no current empirical successor in place; cache retained for future re-fire if a deeper Padé probe of E(z)'s singularity at z=2 becomes load-bearing. See [result_77_7_status.md](result_77_7_status.md).
     - R77.5 (inter-level renormalization residual operator): NOT YET RUN. Identified as the natural follow-up after R77.4 erratum. Would give the operator's spectrum directly.
     - **Convergence shape probe 2026-05-05 ([result_convergence_shape.md](result_convergence_shape.md)) — SUPERSEDED by ε_6 falsification below:** ε_k = S_k − 7/15 was fit with **(1/2)^k + (1/3)^k decomposition**, amplitudes (A, B) ≈ (−0.49, +1.33), AICc beat alternatives by 25 points on k=1..5. Predicted ε_6 ≈ −0.00583 from this fit. Sign pattern (+,+,−,−,−) confirmed (one flip, monotone from below). Algebraic obs: ε_1 = 1/5, ε_2 = 1/(3·5·7) clean; ε_3 onward non-elementary at current resolution. **The two-mode form turned out to be a finite-k OLS artifact** — see ε_6 entry below.
     - **ε_6 prediction test 2026-05-05 ([result_epsilon_6.md](result_epsilon_6.md)):** computed S_6 via float64 power iteration on K_6 (486 states, 7 iterations, residual 2.6e-16). Result: **ε_6 = −4.98×10⁻⁴**, vs two-mode prediction −5.86×10⁻³ — **off by ~10×; two-mode FALSIFIED at k=6** per the brief's decision rule (ε_6 ∈ [−0.005, +0.001] → "TWO-MODE INCOMPLETE / amplitudes wrong"). |ε_{k+1}/ε_k| ratios at k=2..6 were 0.535, 0.482, 0.470, **0.432** — appeared monotone decreasing toward < 1/2; that read was ALSO wrong, see ε_7 entry below.
     - **ε_7 measurement 2026-05-05 ([result_epsilon_7.md](result_epsilon_7.md), [result_epsilon_7_verify.py](result_epsilon_7_verify.py)):** float64 power iteration on K_7 (1458 states, 8 iterations, residual 2.8e-16). Result: **ε_7 = −1.18×10⁻³**, **|ε_7/ε_6| = 2.36 — ratio REVERSED UP, non-monotone trajectory.** S_k has a *local maximum* at k=6 (S_6 = 0.46617, closest to 7/15 yet seen) and backs off at k=7 (S_7 = 0.46549 ≈ S_5). Cross-validated by scipy.eigs: power-iteration and scipy.eigs agree on π_6 and π_7 to L1 distance 1e-15. **Major collateral falsification: the "|ε_n|·2^n stable near 0.04 for n=2..5" rate-1/2 envelope is a finite-k illusion.** Updated envelope at n=2..7: 0.038, 0.041, 0.039, 0.037, 0.032, **0.150** — k=7 jumps 4× the supposed envelope. The "stable near 0.04" was the visible portion of a longer-period oscillation, not a tight envelope. Asymptotic-rate question wide open. Per brief: do not refit until at least k=8. Plausible candidate forms: damped oscillation with complex-conjugate rates ρ·e^±iθ, longer-period structure, non-elementary shape.

11. **Joint 2-3-adic Bohr empirical signal — RETIRED 2026-05-05.**
    - Original claim ([result_bohr_probe.md](result_bohr_probe.md), 2026-05-04): Syracuse iterates show non-CRT-independence on (Z/2^a)* × (Z/3^b)*, z=16.5 at k=20 (a=5, b=4).
    - **Deflated 2026-05-04/05** by bracket-stratification probe ([result_bohr_probe_strat.md](result_bohr_probe_strat.md)): per-bracket chi²/df at v ∈ (10⁶, 10⁹] is 0.95 (z = −0.99) and at v > 10⁹ is 0.94 (z = −1.18) — statistically CRT-independent within ±2σ. The original aggregate signal was driven by the v ≤ 100 bracket-A descent funnel (low-v trajectories transiting toward 1), not by joint structure at the scales relevant to D_emp.
    - **Status:** the joint 2-3-adic Bohr structure does NOT exist at large v. Whatever signal exists lives at v ≤ 10⁶ (descent regime). Item retired as a load-bearing structural object.
    - For audit: do not cite the Bohr finding as a closure path or as evidence of R58/R60 gap structure — it has been retired.

12. **Sibling study (3x+1 vs 3x-1):**
    - **Forward symmetry K_- = σK_+σ proved** ([sibling_3x_minus_1_symmetry_verdict.md](sibling_3x_minus_1_symmetry_verdict.md)): the q=3 Syracuse Markov chains for 3x+1 and 3x-1 are conjugate by negation σ(r) = -r. Implies S_n^{3x-1} = S_n^{3x+1} as exact rationals at every n (verified k=1..4); all R76/R77 derived quantities transfer; c=7/45 is automatic for 3x-1 by the same evidence chain.
    - **Inverse-tree D_n asymmetry** ([duality_S_vs_D_verdict.md](duality_S_vs_D_verdict.md), [duality_followup_verdict.md](duality_followup_verdict.md)): despite forward equality, integer-level inverse-tree Plancherel masses differ by 10³-10⁴× at large depth. After matched-N control (Agent 2 truncated to Agent 3 root-1's |V_n|), residual structural difference is factor ~0.2-4 — ~95% of raw difference is sample-size driven. Forward chain symmetry doesn't propagate to inverse-tree integer-level dynamics.
    - **No clean forward-backward duality** D = f(S) in any candidate form.

13. **q-sweep follow-up** (q_sweep test 2 + c̃_q structure test, partial finding):
    - Universal asymptotic ratio S_{k+1}^{(q)} / S_k^{(q)} → q/3 across all tested q ∈ {3,5,7,11,13,17}.
    - **c̃_q ≈ (q−3)/q** confirmed at q=11, 13, 17 within 1% (q=17 has ord(2 mod 17) = 8 ≠ q-1; non-prim-root status is NOT the differentiator).
    - q=3 separate regime (c̃_3 = 7/15 from forward limit; (q-3)/q form vanishes there).
    - q=5 mild deviation δ ≈ 0.09 (Aitken extrapolation suggests genuine, not finite-k); q=7 large deviation δ ≈ 0.21 (q-specific or finite-k, not the non-prim-root pattern q=17 disconfirmed).
    - Publishable theorem candidate independent of Collatz closure status.

9. **Markov-side first-principles derivation of K** (R77 + R78): local dynamics K_dynamics derivable (per-state Pearson 0.96), but Perron eigvec under uniform-within-state is uniform (D_pred ≈ 1.0). R78 tested path (b) — derive W_visit from R66's 3-adic Bohr π_4. **Result: marginal works (Pearson 0.987 with π_4), per-cell conditional varies wildly, K_full with π_4 weights still gives uniform Perron**. The earlier framing here claimed a "JOINT 2-3-adic Bohr structure that breaks CRT independence within (r mod 32, b) cells" — that claim was the source of the now-retired Bohr empirical positive (item #11). With the bracket-stratification deflation (joint structure CRT-independent at v > 10⁶), the joint stationary on (Z/2^j × Z/3^k)* hypothesis is no longer load-bearing as a closure path. Open question reframed as: what IS the structural object that breaks the per-cell uniform Perron, given that joint 2-3-adic isn't it at large v?

3. **Operator factorization Chang ↔ trajectory** (R69 REJECTED for explicit factorization). Different question: structural relation between σ_Chang ≈ 1 and σ_traj = 0 Fourier classes. Open.

4. **Per-a magnitude pattern** (R72 partial). Asymptotic distribution ≈ Exp(1) but no closed form for individual primitive a values.

5. **σ-band conditional D_avg** characterization beyond R59's mechanism observation. Verify the structural mechanism (D_emp at survivor-time t = inverse-tree depth-(σ−t) marginal) holds at all empirical t values (currently checked at t=10, 30, 50, 70, 90, 110).

6. **Lagarias-Sinai precision** (R68 outcome γ): v ~ Geom(½) heuristic deviates 0.5%-25% at specific j. Whether refining the v-distribution improves R66 closed forms is open. **Refined positive finding:** v_t given m_t mod 2^k is exactly arithmetic-deterministic; the Geom-like marginal arises from trajectory measure being non-uniform mod 2^k for k ≥ 3. Pinning down that mod-2^k profile would give exact closed forms.

7. **R66 4^(-k) decay law** (R74 implies it's wrong). Should re-test |μ̂(a/3^k)|² ~ const · 3^(-k) (rate 1/3, not 1/4) against R65 empirical 0.306, 0.114, 0.023 at k=1,2,3. Likely const = 7/30.

8. **Apply audit reframings** to `lagarias_framework_synthesis.docx` (per `independence_audit.md`). External-facing copy; needs explicit user go before changes.

10. **q-sweep follow-up** (q_sweep test 2 partial finding):
    - Universal asymptotic ratio S_{k+1}^{(q)} / S_k^{(q)} → q/3 across all tested q ∈ {3,5,7,11,13}. Empirically clean (q=11,13 hit 4 sig figs by k=2). **Why exactly q/3?** No analytic explanation yet — could be a clean derivation from Tao's Plancherel framework.
    - c̃_q := lim S_k^{(q)} / (q/3)^k exists for every q tested, but values don't reveal a closed form at 5 q-points. q=7 is the only non-(2-primitive-root) case in the set, possible structural distinguisher.
    - **Publishable theorem candidate** ("the universal q/3 ratio") independent of any Collatz closure attack. Would benefit from q ∈ {17, 19, 23, ...} extension (cheap at k=2; minutes per q).
    - **For c=7/45 closure: nothing changes.** 7/45 confirmed q=3-specific.

14. **Rate-question reframing 2026-05-05 — within-level vs inter-level operator** ([result_qspectrum.md](result_qspectrum.md), combined with [result_epsilon_6.md](result_epsilon_6.md), [result_epsilon_7.md](result_epsilon_7.md), and [result_77_4_K_spectrum_erratum.md](result_77_4_K_spectrum_erratum.md)):

    **q-spectrum probe across q ∈ {3, 5, 7, 11, 13}** returns universal triviality of K_k^(q) spectrum: |λ_1| = 1 isolated, all sub-leading eigenvalues clustered near zero (10⁻⁶ to 10⁻⁵), no q in the tested family has an eigenvalue near 1/2. q-universality holds at the algebraic level — spectrum structure is {1} ∪ near-zero cluster regardless of prime. q=7 anomaly is **not spectral** (q=7's |λ_2| = 4.5×10⁻⁵ sits in same band as q=5's 5.7×10⁻⁵ and q=3's 2.9×10⁻⁵). (q-3)/q closed form at q ∈ {11, 13, 17} **does not originate in K_k spectrum**; it must live in inter-level renormalization.

    **Combined with k=5,6,7 ε_k results and R77.4 erratum:** the rate-1/2 conjecture as a within-level eigenvalue claim is **universally false**. K_k for any prime q is too rapidly mixing to have a 1/2-rate eigenmode; ε_k's apparent rate-1/2 envelope at k=2..5 (now also walked back by the k=7 non-monotone bounce) was never reflected in K_k's spectrum.

    **The convergence rate question is now reframed:** what is the spectrum of the **inter-level renormalization operator** R_k that maps π_k → π_{k+1} (or its associated dynamic on the residual subspace), and is its dominant non-trivial eigenvalue 1/2? This operator was identified as the right object in R77.4 erratum but its spectrum has not been computed at any q. R77.5 ("inter-level renormalization residual operator") was identified as the natural next probe but never run.

    **What this opens (deferred — flagging only):**
    - Build R_k explicitly at q=3 for k=4..7. Compute its top eigenvalues. Test whether |λ_2(R_k)| ≈ 1/2 (rescuing the rate-1/2 claim at the right operator) or some other value (different rate, different structure).
    - Once R_k's spectrum is in hand at q=3, run the same probe at q ∈ {5, 7, 11, 13} for q-universality at the inter-level layer.
    - The non-monotone ε_k bounce at k=6→7 (|ε_7/ε_6| = 2.36) suggests R_k may have **complex eigenvalues** giving oscillation in the projection; this is the natural test.

    **Publication relevance:** the cleanly q-universal "K_k spectrum is trivial" finding is itself a substantive negative result — settles whether per-level Markov mixing is the rate-controlling object across primes. Pairs with the q-universal q/3 growth ratio (item 10) as two separate q-universal theorems independent of Collatz closure status.

    **Update 2026-05-05 (post-compact, new probes):**
    - **Direct K_k top-10 eigenvalue spectrum at q=3, k=5,6,7** ([result_eigenvalue_spectrum.md](result_eigenvalue_spectrum.md)): |λ_2| = 3.42e-4, 4.41e-4, 1.82e-3 — three orders of magnitude below 1/2 at every tested k. Top 10 magnitudes form a tight cluster within ~16% of |λ_2| (no isolated sub-leading mode). Extends R77.4 erratum to k=7; consistent with the q-universal triviality finding above. Three data points don't constrain k → ∞ behavior.
    - **R̃_k = L · K_{k+1}^m · P inter-level construction at q=3, k=4,5,6,7** ([result_R_operator_spectrum.md](result_R_operator_spectrum.md)): for both lifts (uniform A and conditional-from-stationary B) and both K_{k+1} powers (m=1, 2), R̃_k's leading non-trivial eigenvalue lives at 10⁻⁵ to 10⁻³ — algebraically because the construction reduces to K_k on the level-k space. This R̃ does NOT carry the inter-level rate. The probe's brief targeted ρ ≈ 0.984 (now Tier 5 / falsified at k=11), so the negative result is consistent with that walk-back rather than independent evidence. **Side finding (Tier 4):** under uniform lift, ‖L_A · π_k − π_{k+1}‖_∞ approximately halves per k-step at k=4..7 (ratios 0.4998, 0.4992, 0.5000). Striking 3-point pattern in the same k-regime where rate-1/2 was originally fit and later walked back at k=7..11; treat as candidate finite-k transient until reproduced at k=8..11.
    - **Lagarias / Tao 2-adic density-1 v_2 probe** ([result_density_one_v2_bounds.md](result_density_one_v2_bounds.md)) on `data/v_seq_N8388608.parquet` (2,796,202 odd, coprime-to-3 starts in [3, 8388607]): TEST B (mean v > log_2(3)) returned 100% pass rate, but **the result is mathematically tautological** — direct algebra gives mean_v − log_2(3) = (log_2(n_0) + Σᵢ log_2(1 + 1/(3 n_i)))/L > 0 for every Syracuse trajectory terminating at 1, so 100% on the ensemble simply restates "all 2.8M trajectories reached 1" (already known on this range). **The actual empirical content (Tier 3)** is unconditional ensemble mean v_2 = 2.102 (Geom(1/2) prediction 2.0; 5% deviation, consistent with Tao 2019 measure-theoretic asymptotic at finite resolution). TEST A: per-trajectory geometric null `density(v_i ≥ k) ≥ 2^{-(k-1)}` fails for long trajectories (selection effect — conditioning on large L squeezes mean v toward log_2(3)). Do NOT cite TEST B's 100% pass rate as a quantitative density-1 strengthening of the Lagarias connection — it is an algebraic identity, not new empirical evidence.
    - **Mode amplitudes v2 (2026-05-05)** ([mode_amplitudes_v2_findings.md](probe_mode_amplitudes_v2/mode_amplitudes_v2_findings.md)): δ_k = L_{k-1} π_{k-1} − π_k decomposition onto R_k singular vectors confirms R is involved in slow-rate dynamics but **action is band-collective, not single-direction**. σ_1 captures ≤ 3% of ‖δ_k‖² across k=5,6,7 (2.73%, 3.55%, 0% at k=7 where rank 8 leads). Top-20 capture decreases with k (17.86% → 5.22% → 1.65%) — top-20/dim_R shrinks fast as state space grows. ρ_slow ≈ 0.83 **not present in any single K_k or R_k mode at any k tested**. Top-20 R_k singular values cluster tightly in [0.658, 0.671] across all three k — near-degenerate band, no isolated dominant direction. Decomp A (K_k right eigvec inner-product per brief) captured ≈ 0% — confirmed structural (right eigvec of non-symmetric K is biorthogonal-dual to LEFT eigvec, not to itself; right Perron is constant vector and δ_k has sum 0). **Three remaining hypotheses for ρ_slow origin:** (a) **composition across levels** — slow rate emerges from product of R-actions, not single R eigenvalue; (b) **functional-projection averaging** — ε_k is a specific scalar projection of δ_k that selects the "slow envelope" of the band-collective action; (c) **finite-k recurrence-fit artifact** — true rate near 0.60 (composition asymptotic) and the 0.83 fit is a small-k transient. Distinguishing requires k=8..11 ε measurements and computing R_k composition explicitly.

---

## Inactive / parked

- R7 ε_S = log(4) suggestive but not decisive at 50M precision
- R12 piecewise body-Gaussian + tail-GPD — partial closure
- R15 Wiener-Hopf attempt — partial
- R16 Esscher-duality FALSIFIED
- R28 Path B sub-stratum — per-j W_j cannot emerge from residue chain alone
- R29 Edgeworth third-moment shape (R²=0.87, coefficient empirical)
- R31.B Edgeworth standardization wrong direction
- R32.A per-attractor inv-tree spectral bypass FAILS
- R33.B ΔK U-shape was baseline artifact
- R51 QSD framework REJECTED entirely
- R57 H-dim coincidence walked back
- R64.A R60 v2 finer binning — overfitting

---

## Pointers

- **Full archive:** `closed_form_findings.md` (74 results, ~6900 lines) with index at top
- **Detailed writeups by result:** `*.md` files in `C:\Collatz\` matching topic
- **Code + data:** `C:\Collatz\experiments_output\`, `C:\Collatz\inverse_tree\`, `C:\Collatz\data\`
- **Visualization assets:** `C:\Collatz\blender_residue_graph\` (Round 1 + Round 2 PNGs)
- **Reference papers:** `C:\Collatz\lagarias\` (14 PDFs incl. 2603.11066v6, AST_1990, log(φ)/log(2) Chang exact)
- **External-facing synthesis:** `lagarias_framework_synthesis.docx` (do NOT modify without explicit user go)

### Critical scripts (centerpiece, load-bearing)

- `size_stratified_markov.py` — R60 D_avg identification (α-result)
- `mj_resonance_full_partition.py` — R63 |μ̂(1/3)|² closed form
- `higher_q_partition_test.py` — R65 3-adic specificity
- `s_infinity_exact.py` — R70 S_∞ = 7/15 evidence
- `alpha_beta_gamma_decay.py` / `alpha_beta_gamma_weighted.py` — R73 weighted-product rate 1/2
- `lifting_operator_spectral.py` — R74 ‖d‖² rate 1/3 + L_k SVD
- `r73_r74_reconciliation.py` — single-source-of-truth on the two rates
- `kernel_first_principles_v2.py` — R77 derived K vs R60 empirical (γ at Perron)
- `w_visit_derivation.py` — R78 W_visit from π_4 (γ; marginal works, conditional needs joint 2-3-adic)
- `lagarias_sinai_validation.py` — R68 v ~ Geom(½) deviations
- `independence_audit_compute.py` — Validation Task 2 cross-checks

### Audit / validation documents

- `independence_audit.md` — 7 cross-result claims, 6 are "consistent" not "independent"
- `dim_h_validation.md` — R23 / Chang dim walk-back
- `lagarias_sinai_validation.md` — Geom(½) heuristic precision
- `r73_r74_reconciliation.md` — leading-vs-subleading rates

---

## Convention note for parallel agents

When adding a new result:
1. **Check this STATE.md** before claiming "rejected" / "superseded" — many results are revisions of older ones with different framings, not contradictions. Especially watch for "apparent inconsistency" between two rates / values — usually it's leading-vs-subleading or different observables of the same quantity.
2. **Result numbers have collided** (R23, R31, R32, R33, R34, R52, R53, R61, R62, R64, R71 all duplicated). When referencing, use suffix .A / .B AND a topic keyword. Eventually renumber chronologically for v3.7.
3. **Update this STATE.md** when adding a result that supersedes an active claim, locks a new closed form, or closes/opens an item in the active path.
4. **Append don't rewrite** the index at top of `closed_form_findings.md` when adding new results — keep existing entries even if the result is later revised; mark status with ↻ or ✗.
5. **Independence is rare.** Almost every cross-result agreement traces back to the single Collatz-map root. Frame "consistent characterization" not "independent confirmation" unless cross-validating against external input (Tao K_h, Chang's exact formulas, classical Cramér-Lundberg, etc.). See `independence_audit.md`.
6. **Honesty over polish.** This archive includes lots of walk-backs (R57 H-dim, R71 1/2 conjecture, Lagarias-Sinai marginal exactness). Documenting walk-backs is part of the rigor signal — don't quietly omit them from new writeups.
