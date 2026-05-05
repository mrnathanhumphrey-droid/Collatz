# Live State — Collatz framework synthesis

**Last updated:** 2026-05-04 (post-Bohr empirical positive / sibling-study completion / R77.7 superseded snapshot). Today's landmarks: R77.3 falsified the 3-mode geometric ansatz for ε_n over Q; R77.4 envelope fits gave (M) — Jordan ruled out, log/power tied at N=5; R77.4 erratum showed K_k itself has no eigenvalue near 1/2 (not the rate operator); R77.6 generating-function probe found branch-cut signature at z=2 (type indeterminate); **R77.7 (k=7 ε extension) NOT COMPLETED — superseded by Bohr empirical positive; killed at ~8.5 hr; cache retained for future re-fire.** Path-C subroutes all closed (R79b: ‖S_partial‖ saturated β=0.522; saddle-class: structural but not closure path; C2/BGK: random-like, multiplicative energy ≈ N³; C3/band-l¹: ‖ĥ‖_{ℓ¹(D)} saturates at trivial bound; band-spectral: lf_mass→0.25, smooth-completion weakened). **Q-sweep test 2:** literal hypothesis c_q = S_∞^{(q)}/q falsified for q ≥ 5 (S_k diverges); renormalized c̃_q := lim S_k^{(q)}/(q/3)^k exists universally; **c̃_q = (q−3)/q confirmed at q=11, 13, 17 within 1%** (q=17 ruled out non-prim-root explanation for q=7's anomaly). **Sibling study (3x±1):** forward S_n^{3x−1} = S_n^{3x+1} via K_- = σK_+σ chain symmetry (proved at k=1..4); inverse-tree D_n differs structurally (basin fingerprint after matched-N analysis: factor ~0.2-4 residual, ~95% of raw 10⁴× difference is sample-size). **Joint 2-3-adic Bohr empirical positive** (`result_bohr_probe.md`): Syracuse iterates show structured non-CRT-independence on (Z/2^a)* × (Z/3^b)*, growing to z=16.5 at k=20 — direct empirical view of the analytical obstruction R77.7 was reaching for indirectly. **Archive:** `closed_form_findings.md` (79+ results). **Latest results:** `c_seven_forty_fifth.md` (R75), `result_76_conservation_law.md` (R76), `result_77_T_lead_spectrum.md` (R77), `result_77_3_nisoli_bypass.md` (R77.3), `result_77_4_operator_shape.md` + `result_77_4_K_spectrum_erratum.md` (R77.4), `result_77_6_generating_function.md` (R77.6), `result_77_7_status.md` (R77.7 NOT COMPLETED), `result_78.md` (R78), `result_79.md` (R79), `r79b_S_partial_empirical.md` (R79b), `saddle_class_subsum_analysis.md`, `band_l1_analysis.md` (C3), `bk_moments_analysis.md` (C2), `band_spectral_decomposition.md`, `result_q_sweep_test_2_c_q.md`, `c_tilde_structure_verdict.md` + `c_tilde_q17_probe.py`, `sibling_3x_minus_1_symmetry_verdict.md`, `duality_S_vs_D_verdict.md` + `duality_followup_verdict.md`, `result_bohr_probe.md` + verification chain.

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
| R77 implicit framing: "T's spectrum" shapes ε_n envelope | R77.4 erratum | K_k itself has \|λ_2\| ≈ 10⁻⁵..10⁻³ at k=3..6, NO eigenvalue in [0.3, 0.7]. K mixes in O(1) steps. The (1/2)^n envelope is INTER-LEVEL renormalization, not within-level spectral gap. |
| Hypothesis: 7/45 is the q=3 instance of a closed-form c_q family | q-sweep test 2 (NO-PATTERN literal / UNIVERSAL-SHAPE renormalized) | S_∞^{(q)} doesn't exist for q ≥ 5 (sequence grows like (q/3)^k). 7/45 is q=3-specific. But renormalized c̃_q does exist universally — partial finding. |
| R63 {m_j} atomic decomposition for resonance | R63 (revised, full-population partition) | {m_j} chain accounts for only 0.15% of |μ̂(1/3)|²; resonance comes from full-population mod-3 mass asymmetry. |
| Sullivan-conformal measure framing (R59) | R62 | Multifractal with σ = 0 — not in Sullivan/Pollicott-Urbański constant-δ machinery. |
| Bernoulli convolution / Erdős-class lacunary framing (R62) | R63 | Resonance is population-level, not chain-level. Right home: multiplicative number theory on Z_p. |

---

## Open pieces (active path)

1. **Esscher-tilt closure** for the +0.86 → +0.95+ gap in R58. Reduces to Result 22's σ-quartile machinery applied with weight = subtree-size × Esscher-tilt(σ, q≈0.72). Most likely high-leverage move.

2. **Rigorous proof S_∞ = 7/15.** PARTIAL: c = 7/45 now algebraically anchored (R75/R76); rate-½ proof remains.
   - **Plancherel formula (R75, RIGOROUS):** S_k = Σ_{ξ ∈ Z/3^k, 3∤ξ} |μ̂_k(ξ)|². So c = (1/3)·lim Σ |μ̂|² over high-freq.
   - **Conservation law (R76, RIGOROUS):** Σ_{j=0,1,2} M_{n+1}(η_0 + j·3^n) = 0 where M_n(η) := Σ_ξ μ̂_n(ξ) μ̂_n*(ξη).
   - **Leading-mode identity (R76, RIGOROUS):** S_{n+1} = −2·M_{n+1}(1+3^n). Reduces rate question to scalar sequence R_n := M_n(1+3^{n−1}) → −7/30.
   - **Empirical rate ½ verified through k=5** (|ε_n|·2^n stable at C ≈ 0.04 for n=2..5).
   - **Provisional certified bound (R75):** assuming rate ½, |c − S_k/3| ≤ 0.013·(1/2)^k. At k=5, bound 4.2×10⁻⁴, actual 3.8×10⁻⁴.
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
     - **R77.7 (NOT COMPLETED):** k=7 Markov chain extension killed at ~8.5 hr by user direction. Original purpose (distinguish log vs power-law asymptotic for ε_n via deeper Padé) **substantially superseded by joint 2-3-adic Bohr empirical positive** (`result_bohr_probe.md`), which directly observes the analytical obstruction R77.7 was indirectly probing. Cache retained for future re-fire if Bohr finding requires ε_n-side corroboration. See [result_77_7_status.md](result_77_7_status.md).
     - R77.5 (inter-level renormalization residual operator): NOT YET RUN. Identified as the natural follow-up after R77.4 erratum. Would give the operator's spectrum directly.

11. **Joint 2-3-adic Bohr empirical positive** ([result_bohr_probe.md](result_bohr_probe.md), 2026-05-04):
    - Syracuse iterates of N=10⁷ uniform odd starts show structured non-CRT-independence between Z/2^a residues and Z/3^b residues, with chi²-departure z-score growing monotonically with iteration depth k (max z=16.5 at k=20, a=5, b=4).
    - **Direct empirical view of the analytical obstruction** the trajectory measure carries — same multiplicative joint structure that R77.x and R66/R74 closed-forms were probing indirectly via Plancherel sums and ε_n shapes.
    - Verification chain: `result_bohr_probe_verify.md`, `result_bohr_probe_check.md`, `result_bohr_probe_cliff.md`, `result_bohr_probe_asymptote.py`.
    - **Implication:** the joint 2-3-adic Bohr coupling is now the load-bearing structural object. Future closed-form work on c=7/45 should target this directly, not via ε_n's analytic structure (the indirect path).

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

9. **Markov-side first-principles derivation of K** (R77 + R78): local dynamics K_dynamics derivable (per-state Pearson 0.96), but Perron eigvec under uniform-within-state is uniform (D_pred ≈ 1.0). R78 tested path (b) — derive W_visit from R66's 3-adic Bohr π_4. **Result: marginal works (Pearson 0.987 with π_4), per-cell conditional varies wildly, K_full with π_4 weights still gives uniform Perron**. The trajectory measure has a JOINT 2-3-adic Bohr structure that breaks CRT independence within (r mod 32, b) cells. Closing this requires the joint stationary on (Z/2^j × Z/3^k)*, not just R66's mod 3^k alone. **Possibly the same joint structure also gives R74's c = 7/45 — one analytical breakthrough closes both.**

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
