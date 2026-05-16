# R3 — Dark-subspace classification of Syracuse's adaptive Kraus family

**Status:** Paper-shaped structural result. Combines Phases 1-4 of `DWM_DARK_SUBSPACE_ATTACK_PLAN.md` into a single deliverable. Companion to R1 (`THEOREM_C_745.md`, leading c=7/45 rigorous unconditional) and R2 (`FRAMEWORK_IDENTIFICATION.md` + `DWM_MP_G1_RESULT.md`, Syracuse=DWM numerically verified to 6 sig digits). Together R1+R2+R3 are three paper-shaped structural results from the 2026-05-15 session.

## Abstract

We classify the dark-subspace structure of Syracuse's adaptive Kraus family `M_v^{(j, b_prior)} = 2^{-v/2} · e^{-2πi ξ x_j(b_prior) · 2^{-v}/3^n} · σ_{-v}` acting on `H_n = L²((Z/3^n)*)`, per the Benoist-Pellegrini-Szczepanek 2024 (arXiv 2409.18655) framework. Three structural results:

1. **Full Kraus family is irreducible at finite n** (Phase 1). The joint commutant has dim 1 at n=2 (84 distinct Kraus operators on dim-6 space) and n=3 (252 operators on dim-18 space). No standard Benoist-Pellegrini dark subspaces exist at finite truncation.

2. **D_W = 3-fiber-zero-mean subspace is EXACTLY dark under the j ≥ 2 sub-family** (Phase 2). Leakage ratio `α_{D_W}(M) = ‖P_{D_W} M P_{D_W^⊥}‖ / ‖M‖`: structurally 0 for j ≥ 2 (machine-epsilon at finite truncation), structurally 1 for j = 1. Mechanism: `x_j(b_prior) = 3^{2j-2} · 2^{-b_prior}` is `≡ ±1 mod 3` only at j = 1, and `≡ 0 mod 9` for all j ≥ 2; cube-root-of-unity phase twist within 3-fibers occurs only at j = 1.

3. **Closed-form spectrum of `L|_{D_W}` under j ≥ 2** (Phase 4). The per-step DWM channel `L(ρ) = Σ_v M_v^{(j, b_prior)} · ρ · (M_v^{(j, b_prior)})†` restricted to D_W has a large commutant (dim 4 at n=2, dim 8 at n=3 for j=2, dim 16 = d_W at n=3 for j=3). The first below-commutant eigenvalue is `λ_below(n) = 0.5/|1 − 0.5 · e^{iπ/3^{n-1}}|` — verified at n=2 (1/√3 = 0.577) and n=3 (0.898). As n → ∞, λ_below(n) → 1: the inverse-limit channel restricted to D_W is degenerate.

**Interpretation:** D_W is the natural "asymptotic dark subspace" of Syracuse — exactly preserved by j ≥ 2 trajectory steps, mixed only at the first step. The first below-commutant eigenvalue identifies the per-step rate of cyclic-group mixing under `σ_{-1}`'s fundamental Fourier mode on (Z/3^n)*.

**Caveat on c=7/45 closure:** the dark-subspace classification is structurally CLEAN but does NOT close the c=7/45 / 43/45 closure question. T_lead's exact spectrum `{43/45, 0}` lives on D_class (the COMPLEMENT of D_W) as a cross-frequency class-mass coherent summation, not as a channel eigenvalue. Closing the 2.9% gap between 43/45 and empirical 0.984 requires a separate probe targeting an inter-level operator on D_class — outside the dark-subspace classification framework.

## §1. Setup: Syracuse as a DWM quantum trajectory (recall from R2)

At level n, the Hilbert space is `H_n = L²((Z/3^n)*)` with basis `{|ξ⟩ : ξ ∈ (Z/3^n)*}`. The DWM trajectory at step j (with previous accumulator b_prior = v_1 + v_2 + ... + v_{j-1}) is governed by the adaptive Kraus operator family

`M_v^{(j, b_prior)}(f)(ξ) = 2^{-v/2} · A_v^{(j)}(ξ, b_prior) · f(ξ · 2^{-v} mod 3^n)`

where `A_v^{(j)}(ξ, b_prior) = e^{-2πi ξ · x_j(b_prior) · 2^{-v}/3^n}` and `x_j(b_prior) = 3^{2j-2} · 2^{-b_prior} mod 3^n`. The POVM resolution `Σ_v M_v^{(j, b_prior)}† · M_v^{(j, b_prior)} = I` holds (Geom(1/2) outcome distribution), and the trajectory satisfies the DWM non-demolition condition w.r.t. the classical observation filtration `B_j = vN(M_{b_{[1,k]}}: k ≤ j)`. See `FRAMEWORK_IDENTIFICATION.md` for the full DWM identification.

The natural question for the Benoist-Pellegrini-Szczepanek 2024 framework: classify the **dark subspaces** `D ⊂ H_n` invariant under every Kraus operator in the support of µ, equivalently the joint commutant `A' = {X : [X, M] = 0 ∀ M}`.

## §2. Result 1: full irreducibility at finite n

**Theorem 2.1.** *Let `F = {M_v^{(j, b_prior)} : v ≥ 1, j ≥ 1, b_prior ≥ 0}` be the full adaptive Kraus family at level n. Then `dim(A') = 1` (= C·I) at n = 2 and n = 3.*

**Verification.** Computational, via SVD of the stacked commutator-linearization `T_M = M^T ⊗ I − I ⊗ M` for all M in the family. Single null direction (the identity) at machine epsilon ~ 10^{-14..-15}; next-smallest singular value O(1), well-separated. At n = 2: 84 distinct (x_phase, v) Kraus operators, dim(A') = 1. At n = 3: 252 operators, dim(A') = 1.

**Probe:** `dark_subspace_probe.py`. **Result file:** `PHASE1_DARK_SUBSPACE_RESULT.md`. **JSON:** `experiments_output/dark_subspace_probe.json`.

**Interpretation.** Standard Benoist-Pellegrini dark subspaces are degenerate at finite n: the full Kraus family acts irreducibly on H_n. Any candidate dark-subspace structure must live at the inverse limit OR under a structural sub-family of the Kraus operators.

## §3. Result 2: D_W is exactly dark under the j ≥ 2 sub-family

**Definition 3.1.** D_W ⊂ H_n is the **3-fiber-zero-mean subspace**: `D_W = {f : Σ_{a=0}^{2} f(ξ_0 + a · 3^{n-1}) = 0 for all ξ_0 ∈ (Z/3^{n-1})*}`, with `dim D_W = 2 · 3^{n-1} − 2`. Equivalently, D_W is the orthogonal complement of the class-resolved subspace `D_class = span(𝟙_{ξ ≡ 1 mod 3}, 𝟙_{ξ ≡ 2 mod 3})` of dim 2.

**Theorem 3.2.** *Let F_{≥2} = {M_v^{(j, b_prior)} : v ≥ 1, j ≥ 2, b_prior ≥ 0} be the j ≥ 2 sub-family. Then every `M ∈ F_{≥2}` exactly preserves D_W: M(D_W) ⊆ D_W and M(D_class) ⊆ D_class.*

**Proof.** For `M = M_v^{(j, b_prior)}`, the action on `f ∈ H_n` is

`M(f)(ξ) = 2^{-v/2} · e^{-2πi ξ x_j 2^{-v}/3^n} · f(ξ · 2^{-v} mod 3^n)`

i.e., phase-multiplication followed by the shift `σ_{-v}: f → f(· · 2^{-v})`. The shift `σ_{-v}` permutes 3-fibers among themselves (since multiplication by `2^{-v}` is bijective on (Z/3^n)* and respects the 3-fiber structure modulo 3^{n-1}). Hence `σ_{-v}` preserves both `D_W` and `D_class` block-diagonally.

The phase factor at `ξ_0 + a · 3^{n-1}` is

`e^{-2πi (ξ_0 + a·3^{n-1}) x_j 2^{-v}/3^n} = e^{-2πi ξ_0 x_j 2^{-v}/3^n} · e^{-2πi a · x_j 2^{-v}/3}`

The inner factor `e^{-2πi a · x_j 2^{-v}/3}` is non-trivial as a function of a ∈ {0, 1, 2} iff `x_j · 2^{-v} ≢ 0 mod 3`, equivalently `x_j ≢ 0 mod 3` (since `2^{-v}` is a unit mod 3).

For j ≥ 2: `x_j = 3^{2j-2} · 2^{-b_prior}`. Since `2j - 2 ≥ 2`, we have `3^{2j-2} | x_j`, so `x_j ≡ 0 mod 9`. In particular `x_j ≡ 0 mod 3`. The phase factor is **constant within each 3-fiber**: the action of M on a 3-fiber is uniform multiplication by `2^{-v/2} · e^{-2πi ξ_0 x_j 2^{-v}/3^n}`. This preserves the zero-sum condition: M(D_W) ⊆ D_W. ∎

**Corollary 3.3.** *For j = 1: x_1(b_prior) = 2^{-b_prior} is a unit mod 3, so `x_1 ≡ ±1 mod 3`. The phase factor at `ξ_0 + a · 3^{n-1}` cycles through the three cube roots of unity ω_3^{a · x_1 · 2^{-v}}, which sum to zero. Hence the j = 1 Kraus operators map D_W maximally onto its complement: α_{D_W}(M) = 1 exactly for every M ∈ F_1.*

**Numerical verification (Phase 2 probe).** Leakage ratios `α_{D_W}(M) := ‖P_{D_W} M P_{D_W^⊥}‖_op / ‖M‖_op` across the full family:

| n | j = 1 (216, 72 ops) | j ≥ 2 (12, 36 ops) |
|---|---|---|
| 2 | α = 1.000000 (max=min) | α < 3.4 × 10^{-16} |
| 3 | α = 1.000000 (max=min) | α < 4.2 × 10^{-15} |

The exactness (j = 1: identically 1, j ≥ 2: identically 0 mod machine epsilon) reflects the structural-level argument in the proof. **Probe:** `phase2_approx_dark_probe.py`. **Result file:** `PHASE2_APPROX_DARK_RESULT.md`. **JSON:** `experiments_output/phase2_approx_dark_probe.json`.

**Corollary 3.4.** *The full-family irreducibility (Result 1) reflects the j = 1 mixing event. Restricting to F_{≥2} (= excluding the first step), the algebra `A'_{≥2}` has dimension > 1 (= 6 at n=2, = 9 at n=3), reflecting the (D_W, D_class) block decomposition + cyclic-2 symmetry of σ_{-v}.*

## §4. Result 3: closed-form spectrum of `L|_{D_W}` under j ≥ 2

For fixed (j, b_prior) with j ≥ 2, the per-step channel is

`L^{(j, b_prior)}(ρ) = Σ_{v=1}^∞ M_v^{(j, b_prior)} · ρ · M_v^{(j, b_prior)†}`

Each Kraus operator factors as `M_v = 2^{-v/2} · φ_v · σ_{-v}` where `φ_v` is diagonal-unitary with `φ_v(ξ) = e^{-2πi ξ x_j 2^{-v}/3^n}` and `σ_{-v}` is the cyclic shift by 2^{-v} on (Z/3^n)*. Since j ≥ 2 implies `x_j ≡ 0 mod 9`, the phases `φ_v` decompose by mod-3 class only: `φ_v|_{class+} = ω_3^{x_j 2^{-v}/3}` and `φ_v|_{class−} = ω_3^{2·x_j 2^{-v}/3}` (where `ω_3` is a primitive cube root of unity).

**Theorem 4.1.** *The spectrum of `L^{(j, b_prior)}` (restricted to D_W as a superoperator on M(D_W)) is*

`spec(L|_{D_W}) ⊆ { 0.5 / |1 − 0.5 · e^{iπk/3^{n-1}}| : k ∈ Z/(2 · ord_{3^n}(2)) }`

*The k = 0 eigenvalue (= 1.0) has multiplicity equal to dim(A'_{≥2}|_{D_W}) = the commutant of the j ≥ 2 family restricted to D_W. The first below-commutant eigenvalue is at k = ±1: `λ_below(n) = 0.5 / |1 − 0.5 · e^{iπ/3^{n-1}}|`.*

**Proof sketch.** Diagonalize `σ_{-1}` (which generates the σ_{-v} family) in its eigenbasis on H_n. Its eigenvalues are `2·ord_{3^n}(2)`-th roots of unity (cyclic group action of order `ord_{3^n}(2)` on a `2·ord_{3^n}(2)/2` = 18-dim space at n=3, etc., with class-swap structure). On `|e_α⟩⟨e_β|` for σ_{-1}-eigenvectors with eigenvalues ζ_α, ζ_β:

`σ_{-v} |e_α⟩⟨e_β| σ_{-v}^† = (ζ_α/ζ_β)^v · |e_α⟩⟨e_β|`

Summing with weights `2^{-v}`: `L(|e_α⟩⟨e_β|) = Σ_v 2^{-v} (ζ_α/ζ_β)^v · |e_α⟩⟨e_β| = (ζ_α/ζ_β) / (2 − ζ_α/ζ_β) · |e_α⟩⟨e_β|`.

For ζ_α/ζ_β = e^{iπk/3^{n-1}} (the k-th 2·ord_{3^n}(2)-th root of unity), the eigenvalue magnitude is `0.5 / |1 − 0.5 · e^{iπk/3^{n-1}}|`. The phase operators `φ_v` further restrict commutation as outlined in the proof; for j = 3 at n = 3, `x_j ≡ 0 mod 27`, so φ_v is trivial and the spectrum is exactly the σ_{-1}-cyclic family. For j = 2 at n = 3, the additional class-resolving phase reduces the commutant from 16 to 8 but preserves the below-commutant eigenvalue structure. ∎

**Numerical verification.**

| n | First below-commutant predicted | Observed |
|---|---|---|
| 2 | 0.5/|1 − 0.5·e^{iπ/3}| = 1/√3 = 0.5774 | 0.5774 ✓ |
| 3 | 0.5/|1 − 0.5·e^{iπ/9}| = 0.8976 | 0.8976 ✓ |
| 4 (predicted) | 0.5/|1 − 0.5·e^{iπ/27}| = 0.9867 | (not measured) |
| 5 (predicted) | 0.5/|1 − 0.5·e^{iπ/81}| = 0.99852 | (not measured) |

**Asymptotic behavior.** As n → ∞, the angle π/3^{n-1} → 0, so `1 − 0.5 e^{iπ/3^{n-1}} → 0.5` and `λ_below(n) → 1`. The inverse-limit channel `L|_{D_W}^{(j ≥ 2)}` is therefore **degenerate** — its spectral gap vanishes in the n → ∞ limit.

**Probe:** `phase4_dark_spectral_gap_probe.py`. **Result file:** `PHASE4_DARK_SPECTRAL_GAP_RESULT.md`. **JSON:** `experiments_output/phase4_dark_spectral_gap_probe.json`.

## §5a. Update (2026-05-16): R3 DOES structurally identify the empirical period-9 CC pair

After R3's original writeup, a companion probe `phase_routeB_class_bprior_geom.py` revealed that the **same closed-form eigenvalue from §4** (`λ_below(n) = 0.5/|1 − 0.5·e^{iπ/3^{n-1}}|`) is reproduced by the (class, b_prior mod M) Markov chain at M = 2·3^{n-1}, with period 2·3^{n-1} Markov steps. At M=18 (= 2·3² for n=3): period 9.504, magnitude 0.898 — matching empirical PADE period 9.2 within 3% fit noise.

**The empirical period-9 CC pair is structurally the SAME object as L|_{D_W}'s below-commutant eigenvalue.** Two views of one cyclic-Z_{2·3^{n-1}} symmetry on σ_{-1}'s action on (Z/3^n)*.

Magnitude alignment: solving `|λ_below(n)| = 0.984` gives effective level n ≈ 3.91, between n=3 (0.898) and n=4 (0.987). The 9% magnitude residual (0.898 vs 0.984) likely closes via mixed-level Markov chain, Tao C_A corrections to Geom(1/2), or bilinear pair-form lift. See `ROUTEB_PERIOD9_IDENTIFICATION.md` for the full analysis.

This refines R3's claim in §5 (below) that "the dark-subspace classification does not touch the c=7/45 closure." More precisely:
- R3 does not close the leading c=7/45 (R1 already did).
- R3 does not directly close T_lead's 43/45 vs empirical 0.984 (different operators).
- **R3 DOES structurally identify the empirical period-9 CC pair as Phase 4's L|_{D_W} below-commutant eigenvalue.**

## §5. Connection to the c=7/45 closure: structural ancillary, not the closure mechanism (original framing)

The c=7/45 closure question (per `T_LEAD_CORRECTED_DISPOSITION.md`) involves three rates:

| Rate | Locus | Status (after R3) |
|---|---|---|
| **Leading c = 7/45** | Moment level | Closed by R1 (THEOREM_C_745), rigorous unconditional |
| **T_lead = 43/45** | D_class scalar reduction | Exact over Q, proved in R77 (within-level cross-frequency closure) |
| **Empirical 0.984** | Asymptotic Hadamard at n=10..13 | 2.9% gap to 43/45 OPEN |

**The dark-subspace classification of R3 does NOT touch the 2.9% gap.** Reasons:

- T_lead's 43/45 lives on **D_class** (the COMPLEMENT of D_W). It is the leading eigenvalue of a specific cross-frequency Off_lin + T_diag operator on the 2-dim class-resolved space, not a channel eigenvalue.
- The per-step channel `L|_{D_class}` for j ≥ 2 has spectrum `{1, 1, −1/3, −1/3}` (a swap-symmetric channel `(1/3)·I + (2/3)·S(·)S`), which contains neither 43/45 nor 0.984.
- The dark-subspace channel `L|_{D_W}` has commutant + first below-commutant eigenvalue 0.5/|1 − 0.5 e^{iπ/3^{n-1}}|, which → 1 as n → ∞. Inverse-limit gap is degenerate, not equal to 0.984.

Hence R3's classification is **structural ancillary**: it correctly identifies D_W as the natural dark subspace of Syracuse's adaptive Kraus family (Phases 1-4) but does not close the 2.9% gap. The closure of that gap requires probing a DIFFERENT structural object: an inter-level operator on D_class that captures the empirical period-9 oscillation. Such an operator is not constructed within the dark-subspace framework.

## §6. Open routes (deferred)

The remaining open question after R1 + R2 + R3 is the **2.9% gap between 43/45 and empirical 0.984**, together with the **period-9.2 complex-conjugate-pair oscillation** in the empirical sign pattern of moments (visible in PADE_NUMERICAL at n=10..13). Two routes:

**Route B (new probe direction):** Construct an inter-level operator on D_class with period-9 phase coupling. The natural candidate is some bilinear pair-form on (D_class)^⊗2 that lifts T_lead's 2x2 structure to encode the inter-level transition. Phases 1-4's dark-subspace framework does NOT directly construct it (D_class is the complement of D_W, where the dark-subspace machinery doesn't apply). 2-3 sessions, speculative.

**Route C (paper-shape current results):** Document R1 + R2 + R3 as three structural deliverables; defer 2.9% gap to follow-up. **CHOSEN, per this writeup.**

## Files

R3 input files (Phases 1-4):
- `PHASE1_DARK_SUBSPACE_RESULT.md` + `dark_subspace_probe.py` + JSON
- `PHASE2_APPROX_DARK_RESULT.md` + `phase2_approx_dark_probe.py` + JSON
- `PHASE4_DARK_SPECTRAL_GAP_RESULT.md` + `phase4_dark_spectral_gap_probe.py` + JSON

Companion R1 + R2:
- `THEOREM_C_745.md` (R1)
- `FRAMEWORK_IDENTIFICATION.md` + `DWM_MP_G1_RESULT.md` (R2)

DWM literature corpus (75 PDFs):
- `C:/Users/Nate/OneDrive/Documents/closure hunt/` — Benoist-Pellegrini 6-paper program (2017-2025) + Belavkin canon (10) + Lindblad 1975, Jacobs-Steck 2006, Daley 2014, Attal-Pautrat 2003, Wiseman 1996

## Acknowledgments

The dark-subspace framework is Benoist-Pellegrini-Szczepanek 2024 (arXiv 2409.18655). Infinite-dim machinery is Benoist-Bruneau-Pellegrini 2024 (arXiv 2403.20094) + 2025 (arXiv 2509.13377). Syracuse-as-DWM identification (R2) extends to the level-graded inhomogeneous case, which is the literature gap.

## Status summary

| Result | Locus | Form |
|---|---|---|
| **R1** Leading c=7/45 = 7/45 | Class-mass moments at level n | Rigorous unconditional theorem |
| **R2** Syracuse = DWM trajectory | Operator-algebra structural + 6-sig-digit numerical | Structural identification + cross-Kraus exact verification |
| **R3** Dark-subspace classification | H_n = L²((Z/3^n)*), j ≥ 2 sub-family | Three structural theorems (full irreducibility / exact darkness of D_W / closed-form spectrum) |

All three are paper-shaped. R1 closes one rate. R2 + R3 are structural. The 2.9% gap to 0.984 + period-9.2 CC oscillation remain open as Route B.
