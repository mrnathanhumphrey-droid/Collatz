# Live State — Collatz framework synthesis

**Last updated:** 2026-05-04 (post-R79 van der Corput attack on Kalafatelis eq 190; outcome (γ) overall — vdC stalls at rigorous rate ~0.73 (B=1) / ~0.81 (B=2), well above the empirical √N rate ~0.5; even ideal pointwise √N is INSUFFICIENT for eq 190 closure — band-l¹ cancellation in dangerous band is what's actually needed. R79 sub-trivial rate ~0.73 recorded as side product; combined R78 + R79 sharpens the obstruction map for c = 7/45's rigorous closure). **Archive:** `closed_form_findings.md` (79+ results). **Latest results:** `c_seven_forty_fifth.md` (R75), `result_76_conservation_law.md` (R76), `result_77_T_lead_spectrum.md` (R77), `result_78.md` (R78), `result_79.md` (R79).

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
   - **R78 + R79 obstruction map** for the analytical closure (Kalafatelis eq 190):
     - Cochrane Theorem 2 (R78): D = 0 obstruction sharp; trivial-bound only. ❌ closed.
     - Pólya-Vinogradov (R78): worse than trivial for r ≥ 3. ❌ closed.
     - van der Corput B=1 (R79): rigorous rate ~0.73 (sub-trivial side product). ⊳ insufficient for eq 190.
     - van der Corput B=2 (R79): rigorous rate ~0.81 (worse than B=1; iteration loss). ⊳ insufficient.
     - Even *ideal pointwise √N* would NOT close eq 190 (R79 Step 4): need band-l¹ cancellation between m-values in dangerous band, which differencing can't access.
     - Open paths: Bourgain-Konyagin sum-product on ⟨4⟩, direct band-l¹ analysis of ĥ on D_{r,t}(η), or smooth completion (R78 path 2). All substantial research projects.
   - **What's still open (R77 target):** spectral analysis of T_lead operator (acting on R_n sequence) to certify rate ½ rigorously via Nisoli framework. Mod-3 class decomposition μ̂_n = μ̂_n^+ + μ̂_n^- gives a 2D operator structure; build matrix + spectrum.

9. **Markov-side first-principles derivation of K** (R77 + R78): local dynamics K_dynamics derivable (per-state Pearson 0.96), but Perron eigvec under uniform-within-state is uniform (D_pred ≈ 1.0). R78 tested path (b) — derive W_visit from R66's 3-adic Bohr π_4. **Result: marginal works (Pearson 0.987 with π_4), per-cell conditional varies wildly, K_full with π_4 weights still gives uniform Perron**. The trajectory measure has a JOINT 2-3-adic Bohr structure that breaks CRT independence within (r mod 32, b) cells. Closing this requires the joint stationary on (Z/2^j × Z/3^k)*, not just R66's mod 3^k alone. **Possibly the same joint structure also gives R74's c = 7/45 — one analytical breakthrough closes both.**

3. **Operator factorization Chang ↔ trajectory** (R69 REJECTED for explicit factorization). Different question: structural relation between σ_Chang ≈ 1 and σ_traj = 0 Fourier classes. Open.

4. **Per-a magnitude pattern** (R72 partial). Asymptotic distribution ≈ Exp(1) but no closed form for individual primitive a values.

5. **σ-band conditional D_avg** characterization beyond R59's mechanism observation. Verify the structural mechanism (D_emp at survivor-time t = inverse-tree depth-(σ−t) marginal) holds at all empirical t values (currently checked at t=10, 30, 50, 70, 90, 110).

6. **Lagarias-Sinai precision** (R68 outcome γ): v ~ Geom(½) heuristic deviates 0.5%-25% at specific j. Whether refining the v-distribution improves R66 closed forms is open. **Refined positive finding:** v_t given m_t mod 2^k is exactly arithmetic-deterministic; the Geom-like marginal arises from trajectory measure being non-uniform mod 2^k for k ≥ 3. Pinning down that mod-2^k profile would give exact closed forms.

7. **R66 4^(-k) decay law** (R74 implies it's wrong). Should re-test |μ̂(a/3^k)|² ~ const · 3^(-k) (rate 1/3, not 1/4) against R65 empirical 0.306, 0.114, 0.023 at k=1,2,3. Likely const = 7/30.

8. **Apply audit reframings** to `lagarias_framework_synthesis.docx` (per `independence_audit.md`). External-facing copy; needs explicit user go before changes.

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
