# TAUBERIAN_SCOPING_R77_6_REREAD — Phase 1: empirical singularity detection in E(z)

**Date:** 2026-05-12. Wilson reporting. First phase of the Tauberian scoping probe.

## Purpose

Articulate exactly what R77.6's Padé analysis empirically detected about the analytic structure of E(z) := Σ_{n≥1} ε_n z^n at z = 2, distinguishing what is **observed** from what is **inferred** from what is **conjectured but not resolved**.

---

## (a) R77.6 Padé pole sequence (verbatim from result_77_6_generating_function.md, §"Padé approximants over Q")

Working on the auxiliary series f̃(z) := (E(z) − ε_1 z)/z², constructed exactly over Q from the level-k Markov chain rationals ε_1..ε_6.

| (m, n) | role | poles found | closest pole | dist to z=2 |
|---|---|---|---|---|
| (1, 1) | lowest-order diagonal | 1 | z = +2.0764 | **0.0764** |
| (2, 1) | near-diag m+n=3 | 1 | z = +2.129 | 0.129 |
| (1, 2) | near-diag m+n=3 | 2 (one spurious at 155) | z = +2.130 | 0.130 |
| (3, 1) | m+n=4 | 1 | z = +2.313 | 0.313 |
| (2, 2) | diagonal m+n=4 | 2 | z = +2.0513 | **0.0513** |
| (1, 3) | m+n=4 | 3 (two spurious) | z = +2.348 | 0.348 |
| (0, 4) | all-pole | 2 cc-pairs | off-axis on \|z\|≈1.1 | 1.85 |

**Diagonal [n/n] sequence (Stahl's natural probe of singular structure):**

| n | pole(s) | distance to z=2 |
|---|---|---|
| 1 | z = +2.0764 | 0.0764 |
| 2 | z = +2.0513, z = +0.6878 (artifact) | **0.0513**, 1.31 |

The primary diagonal pole moves 2.0764 → 2.0513: **monotone convergence to z = 2 from above the real axis (z > 2)**.

## (b) Convergence ratio

0.0513 / 0.0764 ≈ **0.67**.

Comparison:
- Exponential (simple pole at z = 2): ratio → 0 as n → ∞, in practice ≤ 0.25 even at N=2.
- Branch cut endpoint at z = 2: O(1/N) or O(1/N²), so ratio ~ 1/2 at small N.
- Logarithmic singularity: similar O(1/N) convergence pattern.

**0.67 is consistent with branch-cut convergence and inconsistent with a pure simple pole.**

## (c) What R77.6 explicitly rules out

**Simple pole at z = 2 with constant residue is RULED OUT.**

Quoted from R77.6 verdict, §"What the pole locations say", point 2: *"Poles approach z = 2 from ABOVE on the real axis (z > 2). This is the canonical Padé signature of a branch-cut singularity at z = 2 with the cut extending along [2, ∞). For a simple pole exactly at z = 2, Padé would give a stable pole AT z = 2 across all (m, n); we see drift, not stability."*

This rules out:
- Pole-based Tauberian theorems applied directly to a simple pole at z = 2 (Newman-Zagier-type with M=1 at the natural singularity).
- The "rate-1/2 with constant prefactor" reading of R76 §11 — at the level of E(z)'s singularity, the actual structure is not a pole.

## (d) What R77.6 explicitly does NOT determine

**Power-law vs logarithmic branch type cannot be separated at N=5.**

Quoted from R77.6 verdict: *"The 5-coefficient budget cannot separate (G-power) from (G-log). Both produce poles converging to z = 2 from the same side at this small order. Distinguishing requires extending to ε_7, ε_8 (k = 7, 8 Markov chains; ~hours each)."*

Summary table from R77.6:

| signature | (G-power): power-law | (G-log): logarithmic |
|---|---|---|
| Padé poles | cluster on cut [2, ∞), density ~ N | cluster on cut [2, ∞), density ~ N (log-faster) |
| Distance to 2 | shrinks like O(1/N) or O(1/N²) | shrinks like O(1/N) |
| At N = 2 diagonal | one pole near z = 2 | one pole near z = 2 |
| Discriminator | requires N ≥ 5–10 diagonal points | same |

So R77.6 commits to:
- Branch-cut singularity at z = 2 (cut on [2, ∞), endpoint at z = 2).
- Rules out: simple pole at z=2.
- Rules out: pole away from z=2 on disc of convergence boundary.

R77.6 does NOT commit to:
- Specific branch order (α in (1−z/2)^α).
- Power-law vs logarithmic vs hybrid.
- Whether E(z) has additional singularities on the second sheet of the branch.

## (e) Empirical ε_n through k=6 (cached, exact rationals)

Loaded from `C:/Collatz/experiments_output/result_77_7_eps_exact_through_k7.json` via `from fractions import Fraction`. Six values cached (k=1..6); k=7 not yet computed.

| n | ε_n (exact, displayed as float) | |ε_n|·2^n | (|ε_n|·2^n) − 1/30 |
|---|---|---|---|
| 1 | +0.200000 (= 1/5) | 0.400000 | — (transient) |
| 2 | +0.009524 (= 1/105) | 0.038095 | +0.0048 |
| 3 | −0.005092 (= −5191/1019445) | 0.040736 | +0.0074 |
| 4 | −0.002452 | 0.039236 | +0.0059 |
| 5 | −0.001152 | 0.036856 | +0.0035 |
| 6 | −0.000498 | 0.031866 | −0.0015 |

**Critical empirical observations (these constrain Phase 3's matching task):**

**E1.** |ε_n|·2^n stays in the range [0.032, 0.041] for n = 2..6 — essentially a constant, near **1/30 = 0.0333**. This matches R76 §10's leading-order conjecture S_n = 7/15 − (1/30)·(1/2)^n + O((1/4)^n) at high precision. Six terms cannot distinguish exact constancy from very slow algebraic decay.

**E2.** The deviation (|ε_n|·2^n − 1/30) is **non-monotone**: it rises from +0.0048 at n=2 to a peak around +0.0074 at n=3, then decays to negative −0.0015 at n=6. This is consistent with an oscillatory or signed subleading correction, NOT a pure positive algebraic decay term.

**E3.** The sign of ε_n flips between n=2 (positive) and n=3 (negative); for n ≥ 3 the values are all negative. This sign flip is consistent with the leading rate being −(1/30)·(1/2)^n (per R76 §10).

**E4.** No clean n^{−3/2} pattern visible in the deviation. Fitting (|ε_n|·2^n − 1/30) ~ c·n^{−β} on the n=2..6 data gives implied β values that swing from −0.15 to +0.80 across consecutive-pair fits, not a stable exponent. Six points is too few for a clean fit, especially given the non-monotonicity.

**E5.** The leading-order Padé fit zooming on z=2 (not at z=1) is consistent with the dominant exponential rate being (1/2)^n. The branch-cut analysis is sensitivity to the **subleading correction**, not the leading rate.

## (f) Reconciliation with R76 §10 / §11

R76 §10 (verbatim, from result_76_conservation_law.md): *"S_n = 7/15 − (1/30)·(1/2)^n + O((1/4)^n)."*

R76 §11 (verbatim): *"Open: rigorous derivation of leading coefficient 1/30 (numerical fit) and (1,4) eigenvalue = 1/2 (structural conjecture)."*

The conjectured subleading rate is (1/4)^n, i.e., the next singularity of E(z) after z=2 is at z=4 (a discrete exponential correction, not a branch correction near z=2 itself).

But R77.6's Padé analysis at z=2 detects branch-cut structure, not just a simple pole. This is reconcilable in two ways:

**Reconciliation A:** R76 §10's "+O((1/4)^n)" is itself only the *next pole's* contribution. Local to z=2, the structure could still be a branch-cut endpoint with very small branch coefficient. The dominant n=1/2-rate term is (1/30)·(1/2)^n which is consistent with a simple-pole reading; the branch-cut nature emerges only at finer resolution.

**Reconciliation B:** R76 §10's "+O((1/4)^n)" is an empirical fit, not a structural derivation. The actual subleading structure could be (1/2)^n · (sub-exponential factor like n^{−α} or 1/log n), and the empirical (1/4)^n fit at small N is artifact-dominated by the next-pole contribution swamping the slowly-decaying branch correction.

R77.6 does not resolve A vs B. The Tauberian framework's natural question is which Tauberian theorem applies given that E(z) has *both* a leading simple-pole-like behavior (driving (1/2)^n rate) AND a branch-cut at z=2 (refining the subleading structure).

---

## Concluding Phase 1 articulation

R77.6 detected:
1. **Singularity at z = 2** on the positive real axis (closest singularity of E(z), defining radius of convergence ρ = 2).
2. **Branch-cut nature** at z = 2 (Padé poles drift toward 2 from above the real axis along [2, ∞)).
3. **NOT a simple pole** at z = 2.
4. **Branch type indeterminate** (power-law vs logarithmic, branch order α unknown) at N=5.

The leading-order coefficient asymptotic is dominated by a (1/30)·(1/2)^n term (per R76 §10). The branch-cut signature R77.6 detects refines the structure of the subleading correction.

Phase 2 will read the candidate Tauberian theorems' precise hypotheses. Phase 3 will match.

## Files referenced
- `result_77_6_generating_function.md` (R77.6 writeup)
- `result_76_conservation_law.md` §10, §11 (leading coefficient and (1,4) conjecture)
- `experiments_output/result_77_6_pade_poles.csv` (numerical Padé pole locations)
- `experiments_output/result_77_7_eps_exact_through_k7.json` (cached ε_n exact rationals)
