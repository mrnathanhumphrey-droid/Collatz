# R3a — the 1/9 inverse-tree decay: FALSIFIED (small-window + cutoff artifact)

**Date:** 2026-07-14. **Verdict: H_NONCONVERGENT — the "D_n(k) decays at ~1/9" claim does NOT survive extension. It was a joint n≤6 + E_MAX=30 artifact.** The 1/9 is removed as a "3² spine" entry.

Probe `result_3a_inverse_tree_ninth.py` (exact precision-tower, extends `result_inverse_tree_residue.py`), `result_3a_emax_check.py` (E_MAX sensitivity).

## 1. Method + validation

`D_n(k)` needs only vertex counts mod 3^k. Child residue mod 3^{m−1} is a function of parent residue mod 3^m and e, so the exact counts propagate via a precision tower (track Z/3^{k+n−d} at depth d, reduce one power per forward step). Small (≤~50–60k states), exact rationals, same inverse map / same E_MAX=30 / same D-formula as the committed script.

**Validated:** the tower reproduces all **35/35** committed exact values (`result_inverse_tree_residue.csv`, n=0..6, k=1..5). The tower is correct; the findings below are real, not a bug.

## 2. Extension to n=13 — no 1/9

With the original E_MAX=30, exact ratios `ρ_n(k)=D_{n+1}/D_n` past the n≤6 window (target 1/9=0.1111):

| k | ρ at n=6→7 … 12→13 |
|---|---|
| 2 | 0.0057, 0.585, 0.113, 0.155, 0.868, 0.013, 0.269 |
| 3 | 0.579, 0.155, 0.131, 0.299, 0.012, 0.617, 0.023 |
| 4 | 0.108, 0.105, 0.289, 0.031, 0.158, 0.093, 0.245 |
| 5 | 0.053, 0.234, 0.035, 0.218, 0.058, 0.222, 0.089 |

**No convergence to 1/9** — the ratios oscillate over two orders of magnitude. **H_EXACT_NINTH refuted; H_OTHER_RATIO (clean constant) refuted; H_K_DEPENDENT refuted.** Only **k=1 is stable: ρ=1.000 exactly for n≥2** — the 2/9 fixed point (mod-3 equidistribution), a genuine but *separate* phenomenon, not a 3² decay.

## 3. The confound is real: D_n is E_MAX-cutoff-dependent

`D_n(k=4)` and its ratios at E_MAX ∈ {30,40,50,60}:

| n | D @E30 | D @E40 | D @E50 | D @E60 |
|---|---|---|---|---|
| 6 | 1.76e-5 | 1.78e-5 | 9.83e-7 | 1.20e-7 |
| 8 | 2.01e-7 | 2.03e-8 | 5.36e-9 | 2.79e-10 |
| 12 | 2.63e-11 | 4.85e-13 | 3.77e-14 | 7.14e-15 |

**The D_n values change by orders of magnitude with the cutoff**, and the ratios never approach 1/9 at any E_MAX. The uniform-inverse-tree-count measure with an e≤E_MAX cap is **not a well-defined asymptotic object** — its large-n behavior is dominated by the arbitrary cutoff, not by intrinsic dynamics.

## 4. Verdict and consequences

- **The 1/9 inverse-tree decay is FALSIFIED as an asymptotic rate.** It held only in the n≤6, E_MAX=30 window — a small-window coincidence exactly analogous to R81's ⌊r/2⌋+2. **Remove it from the "9 = 3²" spine** (`STATE.md` / spine table): it is not a genuine instance of the squared-class-mass 3².
- **The right object was never built.** Uniform count + hard e-cutoff is arbitrary. The natural Syracuse inverse measure is the **edge-weighted tree (weight 2^{−e} per predecessor)** — the actual pushforward / harmonic measure, cutoff-free by absolute convergence (Σ_e 2^{−e} < ∞). Any decay-rate question must be asked of THAT object. Deferred to a reformulated probe (own pre-reg).
- **3b (inverse-tree = cycle-count detector) inherits this fragility.** The 3x+1-vs-3x−1 basin fingerprint (`duality_S_vs_D_verdict.md`, already ~95% sample-size per matched-N control) rests on the same cutoff-sensitive object. Do NOT promote it to a quantitative claim without first re-deriving under the 2^{−e}-weighted measure and checking E_MAX/N-robustness.

## 5. Process

The n≤6 "~1/9" verdict in `result_inverse_tree_residue.md` was honest for its window but did not extend or vary the cutoff. This probe did both — validation-gated (35/35) then extended — and the claim did not survive either test. Files: `result_3a_inverse_tree_ninth.py` + `result_3a_dn_tables.csv` + `result_3a_log.txt` + `result_3a_emax_check.py` + `result_3a_emax_log.txt`.
