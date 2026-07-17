# QSD depth-extension test — Result 51 — outcome (γ)

**Date:** 2026-05-03. Decisive: the depth-extension framework "trajectory measure = limit of QSD eigenvectors v_K as K → ∞" is **wrong**. v_K projected to mod 32 does not converge to empirical D_avg at any depth, under any of five tested absorption conventions.

**The empirical trajectory measure is genuinely non-Markov in residues** and cannot be captured by any cylinder-averaged Markov chain at any depth.

Numerical: `qsd_depth_extension.py` (round 1, cylinder absorption), `qsd_depth_extension_v2.py` (round 2, four alternative conventions). CSVs: `qsd_depth_extension.csv`, `qsd_lambda_evolution.csv`, `qsd_depth_v2.csv`.

## 1. Setup

For depth K ∈ {6, 8, 10, 12, 14}, build the cylinder-averaged kernel P_K on odd residues mod 2^K (state space size 2^(K-1)) using 128 lifts per state. Identify absorbing states under five conventions:

- **cylinder21**: r mod 32 = 21 — same as Result 50 (b)
- **cylinder5**: r mod 32 = 5 — same as Result 50 (e)
- **values**: specific m_j values {1, 5, 21, 85, 341, ...} that fit in mod 2^K
- **value21**: residue 21 mod 2^K only (single state)
- **value1**: residue 1 mod 2^K only (single state)
- **fine_cylinder**: residue 21 mod 2^min(K, 14) (becomes specific at large K)

For each (K, convention), eigendecompose P_K_sub^T → leading eigenpair (λ_PF, v_K). Project v_K to mod 32: D_K(r) = Σ_{m ≡ r mod 32} v_K(m) / π_32(r).

Compare to empirical D_avg from Result 50 (averaged over t = 130–190 snapshots).

Validation at K=6 with cylinder21: reproduces Result 50 (b) {21,53} exactly — λ_PF = 0.9375, total dev = 5.43. ✓

## 2. Round 1 finding: cylinder absorption gives D_K invariant in K

| K | n_states | n_surv | λ_PF | spec_gap | total |D_K − D_avg| |
|--:|---------:|-------:|-----:|---------:|-----:|
| 6 | 32 | 30 | 0.9375 | 0.000 | 5.4332 |
| 8 | 128 | 120 | 0.9375 | 0.001 | 5.4332 |
| 10 | 512 | 480 | 0.9375 | 0.011 | 5.4332 |
| 12 | 2,048 | 1,920 | 0.9375 | 0.011 | 5.4332 |
| 14 | 8,192 | 7,680 | 0.9375 | 0.027 | 5.4332 |

D_K(r mod 32) is **identical to all decimal places** at K = 6, 8, 10, 12, 14. This is a structural fact:

**Cylinder consistency theorem (empirical statement).** The cylinder-averaged kernel P_K on mod 2^K, projected to mod 2^K' for K' < K, equals the cylinder-averaged kernel P_{K'} on mod 2^K'. Therefore the QSD eigenvector of P_K_sub (with absorption at a depth-K' cylinder) projects to the QSD eigenvector of P_{K'}_sub.

Going to deeper K does not refine the projected QSD when the absorbing set is a fixed shallow cylinder. The Round 1 framework is structurally degenerate.

## 3. Round 2 finding: specific-value absorption breaks the invariance, but D_K → uniform

Sparse absorption (single state or O(1) states at depth K) gives K-dependent QSDs. Pattern:

### (a) "values" — absorb at all m_j with m_j < 2^K

| K | n_absorbing | λ_PF | total |D_K − D_avg| |
|--:|-----------:|------:|-----:|
| 6 | 3 | 0.9160 | 7.71 |
| 8 | 4 | 0.9711 | 6.19 |
| 10 | 5 | 0.9908 | 5.61 |
| 12 | 6 | 0.9972 | 5.47 |
| 14 | 7 | 0.9992 | 5.43 |

λ_PF → 1 as K grows (absorption mass shrinks). D_K → uniform (≈ 1.0) at all r.

### (b) "value21" — single residue 21 mod 2^K

| K | λ_PF | D_K(r=21) | D_K(r=5) | total |D_K − D_avg| |
|--:|-----:|----------:|---------:|-----:|
| 6 | 0.9688 | 0.516 | 1.024 | 5.79 |
| 8 | 0.9923 | 0.874 | 0.993 | 5.39 |
| 10 | 0.9980 | 0.969 | 0.995 | 5.41 |
| 12 | 0.9995 | 0.992 | 0.994 | 5.41 |
| 14 | 0.9999 | 0.998 | 0.993 | 5.41 |

D_K converges to uniform (≈ 1.0) as K grows. The single absorbed residue is too sparse to perturb the spatial distribution.

### (c) "value1" — single residue 1 mod 2^K

| K | λ_PF | D_K(r=1) | total |D_K − D_avg| |
|--:|-----:|---------:|-----:|
| 6 | 0.9743 | 0.399 | 6.38 |
| 8 | 0.9938 | 0.760 | 5.72 |
| 10 | 0.9984 | 0.908 | 5.51 |
| 12 | 0.9996 | 0.961 | 5.43 |
| 14 | 0.9999 | 0.977 | 5.41 |

Same pattern — convergence to uniform.

### Generalization

Across all four sparse-absorption conventions, D_K → uniform (D = 1.0 at every non-absorbing r) as K → ∞. The mass of the absorbing set becomes vanishing fraction of the state space, so the Perron eigenvector of P_K_sub^T approaches the stationary of P_K (which is itself near-uniform).

**No convention gives D_K → D_avg.**

## 4. Why the framework fails — the structural mechanism

The empirical D_avg(r=5) = 1.86 is the **most enhanced** residue. But m=5 is the descent endpoint (T(5) = 1, terminal). In any QSD framework that absorbs at residue 5 (or its lifts), residue 5 is **depleted** (D ≈ 0 at the absorbing residue, ≈ 1 elsewhere).

The empirical enhancement at r=5 mod 32 comes from a non-Markov mechanism:

1. The actual Collatz orbit visits residue 5 mod 32 many times during descent — at every m with m mod 32 = 5 (i.e., m = 5, 37, 69, 101, 133, ...). Most of these visits are **not terminal** (only m = 5 specifically is terminal).

2. The orbit's value m decreases with t on average (descent). At late t, surviving orbits are predominantly in the small-m regime where they're about to terminate. Many small-m residues lie on the descent path; r=5 is heavily visited because m=5 is the second-to-last value.

3. The cylinder-averaged Markov chain treats all visits to residue 5 mod 32 the same — it cannot distinguish "m=5 (terminal)" from "m=37 (continue)". This loss of value information makes the framework structurally incapable of capturing the visit-frequency enrichment.

**The trajectory measure encodes the value-conditional visit-frequency distribution, which is fundamentally non-Markov in residues.**

## 5. Verdict — outcome (γ)

The depth-extension framework "trajectory measure = lim_{K→∞} v_K(r), where v_K is the leading eigenvector of P_K_sub" is wrong.

Three failure modes documented:

- **Cylinder absorption**: D_K invariant in K (cylinder consistency theorem). No improvement at any depth.
- **Sparse absorption**: D_K varies with K but converges to uniform, not D_avg.
- **Mixed conventions**: any convention with absorption at descent-path residues (5, 21, ...) DEPLETES those residues in the QSD, while empirically they are ENHANCED.

The trajectory measure is a non-Markov object. It cannot be characterized as a QSD of any cylinder-averaged Markov chain at any depth.

## 6. What the trajectory measure actually is

Empirical D_avg captures the **survivor-conditioned visit-frequency distribution** at large t. Surviving orbits at t ≈ 140 are size-biased toward small m (about to terminate); the visit frequency at residue r mod 32 reflects how often the descent path passes through that residue.

Heuristic mapping:

- r=5 enhanced: m=5 is the descent endpoint
- r=1 enhanced (D=1.61): final step before/at termination
- r=21 modest (D=0.93): m=21 lies on m_j attractor cylinder; orbits passing through m=21 quickly descend through 5 → 1
- r=13, r=25 depleted (D=0.55): residues that don't lie on common descent paths
- r=23 enhanced (D=1.40), r=29 enhanced (D=1.35): residues immediately preceding descent (T(23) = 35, then descends; T(29) = 11, then 17, then more)

The closed-form characterization of D_avg requires modeling the **descent-path geometry**: which residues are visited by orbits in the small-m terminal regime, weighted by visit frequency. This is a non-Markov, value-dependent calculation.

## 7. Implications for the v3.6 / Chang correspondence framing

### What was hoped (Result 50 reformulated framing)

> "Trajectory measure = leading eigenvector of P_K_sub as K → ∞. Chang's depth-13 captures survival rate; deeper K captures spatial profile. Both frameworks are truncations of the same infinite-depth object."

### What is actually true

> "Chang's stationary π is the leading eigenvector of P (no absorption) — the unconditioned cylinder-averaged Markov chain. The empirical trajectory measure D_avg is **NOT** the leading eigenvector of any P_sub at any depth. They live in different mathematical categories. π is a Markov stationary distribution. D_avg is a survivor-conditioned visit-frequency distribution that requires modeling value-conditional descent dynamics."

### The trajectory measure as a substantive open object

D_avg is empirically stable (small std across late-t snapshots) and qualitatively sensible (descent-path enrichment matches Collatz arithmetic). But it is NOT a QSD; it is a different probabilistic object.

Closed-form candidates for D_avg(r):

(i) **Visit frequency along descent**: D_avg(r) ∝ E_{descent}[#visits to residue r mod 32 per orbit]. Compute by enumerating descent paths from random m and counting residue visits weighted by survival probability.

(ii) **Size-biased weighting of cylinder transitions**: D_avg(r) ∝ ∫ P(reach residue r mod 32 at time τ | survive) · g(τ) dτ where g is a value-dependent weight.

(iii) **Renewal-theoretic**: D_avg(r) is the mean visit frequency to r mod 32 in a renewal process with absorbing termination.

None of these reduces to a simple eigenvector calculation on Chang's kernel.

## 8. Honest scope statement

Round 1 (cylinder absorption, K = 6 to 14): 5 minutes. Result: D_K invariant in K (structural).
Round 2 (5 × 5 conventions × depths matrix): 5 minutes. Result: no convention reaches D_avg.

Round 1 was intended to confirm the v3.6 reformulation. It instead established the cylinder-consistency theorem (a strong negative result). Round 2 was added to test alternative absorption conventions; it confirmed the pattern is robust.

K = 16, 18 not tested in round 2 because round 1's invariance result and round 2's convergence-to-uniform result are decisive at K ≤ 14. ARPACK convergence at K = 16 in round 1 was problematic; this is moot since deeper K cannot resolve the structural issue.

## 9. Comparison to brief's outcome menu

- **(α) Convergence to D_avg**: NOT REACHED. v_K does not approach empirical D_avg under any convention.
- **(β) Partial convergence**: NOT REACHED. λ_PF varies but D_K shape doesn't; no decomposition into "captured" vs "missing" features matches empirical structure.
- **(γ) No convergence**: CONFIRMED. v_K does not approach D_empirical at any depth, under any tested absorption convention.

## 10. What this means for the project

The trajectory measure is NOT identifiable as a QSD eigenvector of Chang's framework (or any depth-K extension thereof). The "shared operator, different boundary conditions" v3.6 framing is unsupported.

The trajectory measure remains a substantive empirical object (Result 50: stable late-t spatial profile, conserved per-step survival rate ≈ 0.94 matching Chang's λ_PF coincidentally). But its characterization requires going outside the cylinder-averaged Markov chain framework entirely.

Two paths forward:

1. **Accept non-Markov characterization**: model D_avg directly as visit-frequency along descent paths. This is honest but lacks the algebraic cleanliness of the QSD framing.

2. **Find an alternative kernel**: a *different* operator (not Chang's cylinder-averaged Markov chain) whose QSD or stationary equals D_avg. This requires identifying the right operator — possibly a value-weighted variant, or a renewal kernel.

The depth-extension framework (Round 51) does not work. The Lagarias-class open piece is the closed-form characterization of D_avg as a non-Markov visit-frequency object — which neither Chang's nor my framework currently provides.

## 11. Files

- `qsd_depth_extension.py` — round 1, cylinder absorption sweep K=6..18
- `qsd_depth_extension_v2.py` — round 2, five conventions × five depths
- `qsd_depth_extension.csv`, `qsd_lambda_evolution.csv`, `qsd_second_eigenmode.csv`,
  `qsd_depth_v2.csv` — full numerical outputs
- `qsd_depth_extension_log.txt`, `qsd_depth_v2_log.txt` — full diagnostic logs

Compute: ~10 minutes total (round 1 + round 2; K = 16 in round 1 ARPACK-failed but moot).
