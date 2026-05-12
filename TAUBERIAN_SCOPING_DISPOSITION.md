# TAUBERIAN_SCOPING_DISPOSITION — top-level disposition

**Date:** 2026-05-12. Wilson (analyst) reporting to Nathan. Top-level disposition of the first Tauberian-framework probe (post seven-probe spectral trajectory).

---

## DISPOSITION: **H_AMBIGUOUS / INCONCLUSIVE**

> **N=5 data is insufficient to select among candidate Tauberian theorems for E(z)'s singularity at z=2. The Chevalier (2507.15394) square-root Tauberian theorem (Thm 1.14, pure α=1/2 branch) is EMPIRICALLY FALSIFIED at the leading-order level — its prediction ε_n ~ C·(1/2)^n·n^{-3/2} grows by 4.3× across n=2..6 instead of converging to a constant. The structural match between R77.6's branch-cut detection and the candidate Tauberian theorems passes only for Flajolet-Sedgewick Ch. VI's MIXED-singularity framework (leading simple-pole + subleading branch), and within that framework Chevalier Thm 1.16 (meromorphic h, pole of order M ≥ 1) is the closest single-theorem candidate. But the precise value of M and the precise meromorphic-h structure cannot be determined from N=5.**
>
> **The subleading correction δ_n := |ε_n|·2^n − 1/30 is non-monotone and changes sign between n=5 and n=6. This strongly suggests multiple competing subleading terms — possibly multiple singularities (z=2 branch endpoint + secondary singularities on the second sheet) — which 5 data points cannot fit cleanly. No single-term ansatz (n^{-3/2}, n^{-1}, n^{-2}, exponential) reproduces this behavior.**

---

## Why H_AMBIGUOUS and not the alternatives

- **H_SQUARE_ROOT_MATCHES_PLUS_EMPIRICAL** (Chevalier 2507.15394 Thm 1.14 prediction matches empirical to within tolerance): **REJECTED**. The prediction n^{-3/2} is off by 4.3× across n=2..6.

- **H_SQUARE_ROOT_MATCHES_BUT_EMPIRICAL_AMBIGUOUS** (singularity is √-branch but N=5 insufficient): **REJECTED**. Even at the *subleading* level, the data does not fit α=1/2. δ_n behavior is inconsistent with any single-α power.

- **H_GENERAL_BRANCH_MATCHES** (FS Ch. VI is right framework, branch type undetermined): **CONSISTENT but underdetermined**. FS Ch. VI is the right framework, but the framework alone does not pick out a single theorem; it's a toolkit.

- **H_NEWMAN_ZAGIER_MATCHES**: **REJECTED**. Newman-Zagier requires the singularity to be a pole (removable by analytic continuation). R77.6 explicitly rules out a simple pole at z=2 in the branch-cut sense. Newman-Zagier is Dirichlet-series anyway; not the right setting for power-series E(z).

- **H_NEEDS_DIFFERENT_TAUBERIAN**: **PARTIALLY**. The mixed-singularity reading (leading simple-pole + subleading branch) is closest to Chevalier Thm 1.16 (meromorphic h with pole of order M at 0), but exact M is indeterminate.

- **H_AMBIGUOUS / INCONCLUSIVE**: **CHOSEN**. Honest disposition: the Tauberian framework is the right tool, FS Ch. VI is the right framework, but no single theorem is selected by N=5 data.

---

## Phase-by-phase summary

**Phase 1** (TAUBERIAN_SCOPING_R77_6_REREAD.md): R77.6 detected branch-cut at z=2 (rules out simple pole), but cannot distinguish power-law vs log at N=5. Cached ε_n through k=6 shows |ε_n|·2^n ≈ 1/30 (matching R76 §10's leading coefficient) with small non-monotone deviations.

**Phase 2** (TAUBERIAN_SCOPING_THEOREM_STATEMENTS.md): Read the four primary Tauberian theorem statements. Chevalier 2507.15394 covers α-branch singularities for power series, no non-negativity needed. FS Ch. VI is the broader singularity-analysis framework. Newman-Zagier is Dirichlet-series and pole-specific (doesn't apply). PTBZ guide is Dirichlet-series-centric.

**Phase 3** (TAUBERIAN_SCOPING_MATCH.md): The matching reveals that Chevalier Thm 1.14's pure α=1/2 prediction is incompatible with the empirical leading behavior |ε_n|·2^n ≈ const. The data is qualitatively consistent with FS Ch. VI's mixed singularity (leading pole + subleading branch), with Chevalier Thm 1.16 (meromorphic h) as the closest single theorem.

**Phase 4** (TAUBERIAN_SCOPING_VERIFICATION.md): Empirical verification at n=2..6:
- Chevalier Thm 1.14 (α=1/2): FALSIFIED (n^{3/2} product grows 4.3×).
- Chevalier Thm 1.16 (M=1, predicting n^{-1/2}): WEAK SUPPORT (1.5× range, peak at n=5).
- Simple-pole leading + branch subleading: QUALITATIVELY CONSISTENT but precise subleading undetermined.
- δ_n = |ε_n|·2^n − 1/30 is non-monotone, sign-flipping. No single-term subleading model fits.

**Adversarial checks** (per A1-A5 in brief):
- A1: R77.6's N=5 limitation is HONORED — not over-claiming disambiguation power.
- A2: Theorem hypothesis fidelity: Newman-Zagier excluded because it's Dirichlet-series, not power-series. Chevalier theorems require power-series with specific singularity types; matching attempted carefully.
- A3: Six data points is acknowledged as small-sample. Strength-of-consistency reported honestly.
- A4: No Tauberian-on-Collatz prior literature — disposition is "this Tauberian theorem appears most plausible" not "this is proven."
- A5: For each candidate, what additional data resolves is documented.

---

## Trajectory placement

The seven-probe spectral trajectory (T_3 → R_k → Candidate A → R76 §11 2D → T_N → Cross-frequency → T_V) converged on H_M_RECURSION_UNDERSPECIFIED — every operator construction over Q at finite truncation hit a structural wall. The consistent reading was: rate-1/2 lives at a branch-cut endpoint of a continuous-spectrum operator on infinite-dimensional space, not a discrete eigenvalue.

**The Tauberian framework arc opens** by accepting the branch-cut reading directly and asking: which Tauberian theorem extracts the rate-1/2 from E(z)'s singularity structure? This first probe finds:

1. **Framework match: Flajolet-Sedgewick Ch. VI singularity analysis.** This is the canonical operator-free toolkit for power-series coefficient asymptotics from local singularity behavior. Chevalier 2507.15394 is a specific case (square-root branch with α=1/2) and Thm 1.16 (meromorphic h with pole of order M at 0) is its closest single-theorem candidate for the mixed structure E(z) appears to have.

2. **Single-theorem match: undetermined at N=5.** The empirical data falsifies the cleanest Tauberian theorem (Chevalier Thm 1.14, pure α=1/2). The closest match (Thm 1.16) requires knowledge of M (the pole order of the meromorphic h at 0), which is not determinable from N=5.

3. **No Tauberian-on-Collatz prior literature** (confirmed by Hank's INDEX). This work is novel framework application.

---

## Recommended next probe

The natural next step to push the Tauberian framework forward, in priority order:

### Route A: Compute ε_7 (and possibly ε_8) — HIGHEST PRIORITY

R77.6 estimates ε_7 (k=7 Markov chain) takes ~hours of compute (Gauss elimination over Q on a 1458-state chain). This single new data point would:
- Add a diagonal Padé [3/3] point: predicted z-coordinate ≈ 2.030-2.040 if power-law, ≈ similar but slower if log.
- Constrain the subleading δ_7 to confirm or refute the sign-flip pattern at n=5→6.
- Allow a 6-point log-log fit on δ_n with one extra degree of freedom.

With ε_7 + ε_8, the Tauberian framework's branch-type discrimination becomes feasible. **Recommended.**

### Route B: Fit a multi-term mixed-singularity ansatz to existing data

Without new ε_n data, attempt to fit a more sophisticated form:

ε_n = A·(1/2)^n + B·(1/2)^n/n^α + C·(1/2)^n·cos(γn + φ)/n^β

to the available 5 subleading data points. Six free parameters with 5 data points — underdetermined, but tests whether the non-monotonicity is consistent with secondary-singularity oscillation vs. transient.

This is a **partial-information probe** and provides indicators rather than discrimination. Lower priority than Route A.

### Route C: Differential approximants (D-Padé)

A D-Padé approximant fits the ODE that f̃(z) satisfies near its singularity, encoding the branch exponent directly rather than reading it off pole drift. Can sometimes give cleaner exponent extraction at lower N than ordinary Padé. Requires additional numerical infrastructure (Maple or specialized Python library). **Medium priority.**

### Route D: Direct contour-integral attack on Chevalier Thm 1.14's hypothesis verification

If we cannot verify Chevalier Thm 1.14's hypothesis empirically (g(z) = h(√(1−z/2)) form), perhaps we can verify it **structurally** by computing E(z) − leading-pole-part directly and asking whether the remainder is of the form h(√(1−z/2)). This requires the closed-form of E(z), which is not available for Collatz. **Low priority** unless a closed-form path appears.

### Route E: Document this scoping probe as Tauberian-framework groundwork; publish

The structural reading (Tauberian framework is the right operator-free tool; FS Ch. VI is the right level of abstraction; Chevalier Thm 1.16 is the cleanest candidate; N=5 is the limiting factor) is itself publishable as part of the seven-probe + Tauberian-arc trajectory. Combined with the bilinear-bound side (burgess.zip, 25-commit work), the publishable claim expands.

This is the **most honest framing without committing to multi-session work**. Could pair with Route A as "we found the framework, here's what we'd need next."

---

## Deliverables produced (in C:/Collatz/, all prefixed TAUBERIAN_SCOPING_*)

- **TAUBERIAN_SCOPING_R77_6_REREAD.md** (Phase 1, this probe's articulation of R77.6's empirical finding)
- **TAUBERIAN_SCOPING_THEOREM_STATEMENTS.md** (Phase 2, candidate theorems read)
- **TAUBERIAN_SCOPING_MATCH.md** (Phase 3, matching table)
- **TAUBERIAN_SCOPING_VERIFICATION.md** (Phase 4, empirical falsification of α=1/2)
- **TAUBERIAN_SCOPING_DISPOSITION.md** (this file)
- **tauberian_verify.py** (numerical verification script for main-thread execution)

## Files referenced

- Internal: `result_77_6_generating_function.md`, `result_76_conservation_law.md` §10 §11, `STATE.md` item 11, `experiments_output/result_77_7_eps_exact_through_k7.json`, the seven-probe DISPOSITION.md files (CROSS_FREQ, T_V, T_N, R76_S11, CANDIDATE_A, R_K, M3).
- External (Tauberian corpus): `arxiv_2507.15394_Tauberian_Square_Root_Singularity.pdf` (Chevalier 2025), `arxiv_2504.16233_Guide_Tauberian_Arithmetic_Apps.pdf` (Pierce-Turnage-Butterbaugh-Zaman 2025), `Newman_1980_Simple_Analytic_Proof_PNT.pdf`, `Flajolet_Sedgewick_Analytic_Combinatorics.pdf` (Ch. VI singularity analysis, not page-read but content referenced via Chevalier's introduction and the PTBZ guide). Hank's INDEX: `C:/Users/Nate/OneDrive/Documents/tauberian/INDEX.md`.

## Synopsis (one paragraph)

The first Tauberian-framework probe scopes which Tauberian theorem matches E(z)'s singularity at z=2 as detected by R77.6's Padé analysis. Phase 1 articulates that R77.6 detected a branch-cut at z=2 (rules out simple pole) but cannot disambiguate power-law from log at N=5. Phase 2 reads four candidate theorem statements: Chevalier 2507.15394 (square-root and α-branch power-series Tauberian), Flajolet-Sedgewick Ch. VI (broader singularity-analysis framework), Newman-Zagier (Dirichlet-series, pole-specific — excluded), and PTBZ guide (Dirichlet-series, framework). Phase 3 matches: Chevalier Thm 1.14's pure α=1/2 prediction ε_n ~ C·(1/2)^n·n^{-3/2} is empirically falsified (grows 4.3× across n=2..6). The FS Ch. VI framework with mixed singularity (leading simple-pole-like + subleading branch at z=2) qualitatively fits, with Chevalier Thm 1.16 (meromorphic h with pole of order M at 0) as the closest single-theorem candidate. Phase 4 confirms: the subleading correction δ_n = |ε_n|·2^n − 1/30 is non-monotone and sign-flips between n=5 and n=6, indicating multiple competing subleading terms — likely a leading branch + secondary singularity oscillation. No single-term ansatz fits. **Disposition: H_AMBIGUOUS / INCONCLUSIVE.** N=5 is the limiting factor. Recommended next probe: compute ε_7 (and ideally ε_8) via the k=7 Markov chain (estimated hours of compute), which adds a critical diagonal Padé point and lets the empirical fit distinguish branch type. The Tauberian framework is the right operator-free tool for E(z) — confirmed structurally. The specific Tauberian theorem within it that closes the asymptotic extraction is not yet pinpointed.
