# Result 20 — PHASE 1c: ★ R13's collision law GENERALIZES (level 2, verified on 1.86M colliding pairs). ⚠️ v1's headline hypothesis was a TAUTOLOGY. My `d/(q−1)` level-rate prediction LOST; the truth is `1/q`.

**Date:** 2026-07-15. **Verdicts: ★ H_L2 ✓ CONFIRMED at scale / ⚠️ H_CASC (v1) was VACUOUS — dropped, not salvaged / ⚠️ H_LEVELRATE — my prediction lost.**

**Headline: the level-2 condition generalizes — `2^{−S'_2} ≡ 2^{−S_2} + j_1·s·2^{−S_1} (mod q)` with suffix sums `S_1,S_2` in place of R13's `v_2,A`. Verified on 1,861,921 colliding pairs to k=4 across four primes, zero failures. Levels ≥3 are NOT derived.**

Probe: `probe_20_phase1c_cascade.py` (v2). Log: `result_20_cascade_log.txt`. Runtime: **3.7 s**.

## ⚠️ FIRST: v1 of this probe was defective, and it is recorded not hidden

v1's headline hypothesis **H_CASC** ("collision ⟺ all k cascade levels pass") was a **TAUTOLOGY**. `cascade_levels` computes `R = Σ_m q^{m−1}T_m` and divides by `q` k times; "passed all k" ⟺ `q^k | R` ⟺ `R ≡ 0 mod q^k` ⟺ **collision, by definition.** Its `0` false-positive / `0` false-negative columns were guaranteed by arithmetic and said nothing about Collatz.

**This is the seventh pre-registration defect this arc and the worst kind** — the first six were mis-specified *thresholds*; this was a hypothesis with *no content*. v2 **drops it entirely** rather than salvaging it.

## ★ H_L2 — the real claim, confirmed at scale

**Derivation (before running).** `value = Σ_{m=1}^{k} q^{m−1}2^{−S_m} mod q^k`, `S_m = v_{k−m+1}+…+v_k`; collision ⟺ `Σ_m q^{m−1}T_m ≡ 0 mod q^k`, `T_m := 2^{−S_m} − 2^{−S'_m}`.
Level 1: `T_1 ≡ 0 mod q` ⟺ `S'_1 = S_1 + j_1 d`. Writing `2^d = 1+q·s` and expanding `2^{−j_1 d} = 1 − j_1 q s + O(q²)` gives `T_1 = q·U_1` with `U_1 ≡ j_1·s·2^{−S_1} (mod q)` — **the non-trivial step.** Dividing by q and reducing mod q:

> ## `2^{−S'_2} ≡ 2^{−S_2} + j_1·s·2^{−S_1} (mod q)`

**At k=2, `S_1 = v_2` and `S_2 = v_1+v_2 = A` — exactly R13.** So R13 is the k=2 shadow of a general law.

| q | k | cells | collisions | L2 bad (coll) | sample | L2 bad (non) |
|---|---|---|---|---|---|---|
| 3 | 4 | 11,664 | **1,253,880** | **0** | 200,000 | **0** |
| 5 | 3 | 8,000 | **316,000** | **0** | 200,000 | **0** |
| 7 | 3 | 9,261 | **287,091** | **0** | 200,000 | **0** |
| 11 | 2 | 1,100 | 4,950 | **0** | 200,000 | **0** |

Every colliding pair satisfies the explicit form; on sampled non-collisions the explicit form never disagrees with the cascade's level-2 test. **Exact set equality — no threshold to mis-specify** (deliberate, given six mis-specified thresholds this arc).

## ⚠️ H_LEVELRATE — my prediction LOST. The truth is `1/q`.

| q | k | d | `1/d` | `d/(q−1)` (my prediction) | measured P(pass m \| passed m−1) |
|---|---|---|---|---|---|
| 3 | 4 | 2 | 0.5000 | 1.0000 | L1: **0.5001** · L2: 0.3315 · L3: 0.3354 · L4: 0.3292 |
| 5 | 3 | 4 | 0.2500 | 1.0000 | L1: **0.2485** · L2: 0.1996 · L3: 0.1945 |
| 7 | 3 | 3 | 0.3333 | 0.5000 | L1: **0.3312** · L2: 0.1427 · L3: 0.1391 |
| 11 | 2 | 10 | 0.1000 | 1.0000 | L1: **0.0993** · L2: 0.0848 |

- **Level 1 lands on `1/d` almost exactly** ✓
- **Levels ≥2: my `d/(q−1)` is WRONG** (predicts 1.0 at q=3,5; measured 0.33, 0.20). The rates sit at **≈ `1/q`** and **stabilize** (`L2 ≈ L3 ≈ L4`).

| q | `1/q` | measured L≥2 |
|---|---|---|
| 3 | 0.3333 | 0.3315, 0.3354, 0.3292 |
| 5 | 0.2000 | 0.1996, 0.1945 |
| 7 | 0.1429 | 0.1427, 0.1391 |
| 11 | 0.0909 | 0.0848 |

**Why I was wrong:** each level asks whether a *specific* residue equation holds mod q — probability `1/q` — not whether something lands in `H` (probability `d/(q−1)`). **I had the wrong event.** *(Pattern: structural priors land, quantitative priors lose. Now ~20-for-31.)*

## Byproduct — flagged as a HYPOTHESIS, explicitly not a result

`P(collision) ≈ (1/d)·(1/q)^{k−1}` would give:

| q,k | predicted | measured | off |
|---|---|---|---|
| 3,4 | 1,259,600 | 1,253,880 | 0.5% |
| 5,3 | 320,000 | 316,000 | 1.3% |
| 7,3 | 291,700 | 287,091 | 1.6% |

**NOT CLAIMED.** It is post-hoc pattern-matching on five conditional rates — exactly what `feedback_r2_cannot_discriminate_monotone_fits` exists to prevent — and it is the **unweighted** count, whereas `c_k` is a **weighted** mass, so it is not even the right object. Needs its own pre-committed test at unseen (q,k) before it means anything.

## Phase 1c — honest scope

**DELIVERED:** the level-2 condition generalizes (R13's law with `S_1,S_2`), verified at scale.
**NOT DELIVERED:** the explicit form at levels m ≥ 3. Level 3 reads `W_2 + T_3 ≡ 0 mod q` with `W_2 := (U_1+T_2)/q` — a **definition, not a formula**. An explicit `W_2` in terms of `(j_1, j_2, s, S_*)` needs a second-order expansion of `2^{−jd} mod q³`, which I have not done.

## Plan status

| phase | status |
|---|---|
| 0, 1, 1b | **DONE**, exact |
| **1c** | **PARTIAL** — level 2 delivered + verified at scale; levels ≥3 NOT derived |
| **2 — literature (Konyagin–Shparlinski)** | **NEXT** — may hand us the general form rather than deriving it |
| 3 — the bound `c_k ≤ C_q·r_q^k` | the math |

## Not at stake
R10–R19. A refutation would kill only the Phase-1c level-2 derivation.

_Reporting discipline: v1's vacuous hypothesis is disclosed in this file AND in the probe's own docstring, and was dropped rather than salvaged — a tautology that "passes" is worse than a rule that fires wrong. H_L2 is an exact set equality across 1.86M colliding pairs. H_LEVELRATE was pre-committed as a MEASUREMENT with no verdict, which is the only reason the `1/q` finding is reported as a loss for my prediction rather than dressed as a discovery. The `(1/d)(1/q)^{k−1}` byproduct is explicitly withheld as post-hoc._
