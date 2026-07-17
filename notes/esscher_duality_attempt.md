# Result 16: Esscher-tilted ascending-ladder duality attempt for closed-form E[L⁻]

**Date:** 2026-05-02. Sequel to `path_c_derivation.md` (Result 15).

This document attempts the closed-form derivation of E[L⁻] (strict descending ladder mean of the iid Syracuse log-walk under P) via the Esscher-tilted measure P* and the conjectured ascending-ladder duality. Numerical verification: `esscher_duality.py`.

---

## 1. Setup

**Walks:**
- P: v ~ Geom(1/2). Step X = log(3) − v·log(2). Mean E[X] = −log(4/3) ≈ −0.288 (negative drift).
- P*: v ~ Geom(3/4) (Esscher tilt at w* = 1). Step X same form. Mean μ* = log(3·2^(−4/3)) ≈ +0.174 (positive drift).

**Empirical targets (Path C):**
- E_P[L⁻] = 1.00456 ± 0.00061 nats (10⁷-orbit simulation)
- q = P(σ⁺ < ∞) = 0.71355 ± 0.00028 (10⁷-orbit simulation, separate run)

**Identity established (Result 15, derivation in this document Section 2):**

> **E[L⁻] = μ / (1 − q)**, where μ = log(4/3), q = P(σ⁺ < ∞)

Verified to 10⁻⁴ precision: μ/(1−q) = 0.287682/0.286348 = 1.00466 vs simulated 1.00456.

## 2. Wiener-Hopf identity for E[L⁻] (proof)

WH factorization for the negative-drift walk:
> 1 − φ(θ) = (1 − κ⁻(θ)) · (1 − κ⁺_def(θ))

where κ⁻ proper (κ⁻(0) = 1), κ⁺_def defective (κ⁺_def(0) = q < 1).

Differentiating at θ = 0: 1 − φ(0) = 0, both factors share the trivial zero. Use L'Hôpital:

> (d/dθ)(1 − κ⁻(θ))|_{θ=0} = (d/dθ)[(1 − φ(θ)) / (1 − κ⁺_def(θ))]|_{θ=0}
>                          = [−φ'(0)·(1 − q) − 0·(−κ⁺_def'(0))] / (1 − q)²
>                          = −φ'(0) / (1 − q)
>                          = i·μ / (1 − q)

LHS = i·E[L⁻] (since κ⁻'(0) = i·E[S(σ⁻)] = −i·E[L⁻]).

Equating: **E[L⁻] = μ/(1 − q)**. ∎

## 3. The user's conjectured duality: E[L⁻] = E*[L⁺*]

The brief proposed: by Feller Vol II Ch XII duality, the descending ladder mean under P equals the ascending ladder mean under P*. If E*[L⁺*] has closed form via tilted-walk machinery, E[L⁻] follows.

**Direct simulation test (5×10⁶ orbits each, numba-accelerated):**

| Quantity | Simulated value | 95% CI |
|---|---|---|
| E_P[L⁻] (descending under P) | **1.00466** nats | ±0.00087 |
| E_P*[L⁺*] (ascending under P*) | **0.34426** nats | ±0.00010 |
| Difference | **+0.66040 nats** | — |
| z-score for difference = 0 | **1480σ** | — |

**Verdict on the duality identity: FALSIFIED at >1000σ.**

The descending ladder mean under P (1.0047 nats) and the ascending ladder mean under P* (0.3443 nats) are not equal. They differ by a factor of nearly 3.

**Why the duality fails:** the Feller Vol II duality I am familiar with operates via *time reversal*, not Esscher tilt. Time-reversed walk has steps −X = v·log(2) − log(3), which is genuinely different from Esscher-tilted walk (v ~ Geom(3/4) with same step formula). The descending ladder under P has the same distribution as the ascending ladder under the *time-reversed* walk (by Feller's duality lemma, Vol II Ch XII Lemma 1) — but the time-reversed walk is **not** the Esscher-tilted walk.

For our specific walk:
- Time-reversed: −X has values {log(4/3), log(8/3), log(16/3), ...} ∪ {−log(3/2)} with probs {1/4, 1/8, 1/16, ...} ∪ {1/2}
- Esscher-tilted: X has same support as original {log(3/2), log(3/4), log(3/8), ...}, but with v ~ Geom(3/4)

Different walks. The ascending ladder under each is a different random variable.

## 4. What the test DID confirm

Three structural identities are verified at high precision:

**(a)** E[L⁻] = μ/(1 − q) (negative-drift WH derivation): 1.00466 = μ/(1−q) ✓ [Section 2]

**(b)** E*[L⁺*] = μ*/(1 − q*) (positive-drift WH derivation, analogous): 
> 0.34426 (sim) vs μ*/(1 − q*) = 0.17442/0.50674 = 0.34419 ✓

**(c)** Esscher relation: q* = P*(σ⁻* < ∞) = E_P[e^(−L⁻)] at w* = 1.
> Direct sim: q* = 0.49326 ± 0.000218
> Esscher: E_P[e^(−L⁻)] = 0.49334 ± 0.0006
> diff = −7.3 × 10⁻⁵ (−0.33σ) ✓

These three identities form a **closed system** but with two free parameters (q and q*) related by the Esscher relation:

> 1 − q* = 1 − E_P[e^(−L⁻)]
> E[L⁻] = μ/(1 − q)
> E*[L⁺*] = μ*/(1 − q*)

To close E[L⁻] in closed form, we need q in closed form, which the Esscher tilt doesn't provide directly.

## 5. Sparre-Andersen formula for q — and why it's not closed-form

For any iid walk with E[X] < 0:

> 1 − q = exp(−Σ_{n=1}^∞ (1/n) · P(S_n > 0))

For our walk: P(S_n > 0) = P(T_n ≤ ⌊n · log₂ 3⌋), where T_n ~ NegBin(n, 1/2).

The truncation at ⌊n · log₂ 3⌋ involves the Gelfond-Schneider transcendental log₂ 3. By Weyl equidistribution, {n · log₂ 3} mod 1 is uniform on [0,1), so the truncation pattern is "random" in the Diophantine sense.

The series Σ (1/n) P(S_n > 0) can be computed numerically (truncated at large N), but no closed-form summation is available. This is the **same obstruction** identified in Path C.

## 6. Closed-form candidates for q — all rejected at 10⁷-orbit precision

Tested in `path_c_q_attempt.py` against empirical q = 0.71355 (SE = 0.00014):

| Candidate | Value | z-score | Verdict |
|---|---|---|---|
| 5/7 = 0.71429 | 0.71429 | 5.1σ | rejected |
| 1 − log(4/3) = 0.71232 | 0.71232 | 8.6σ | rejected |
| log(2)/(log 2 + log(4/3)) | 0.70670 | 48σ | rejected |
| 3/4 (= P*(v=1)) | 0.75000 | 255σ | rejected |
| 1 − 1/e | 0.63212 | 570σ | rejected |
| log(2) | 0.69315 | 143σ | rejected |

**No simple closed-form candidate matches q within sampling error.** 5/7 is closest at 5.1σ off — *suggestive* of underlying rational structure, but rejected at this resolution.

## 7. Verdict on Result 16

| Aim | Status |
|---|---|
| E[L⁻] = E*[L⁺*] (user's duality claim) | **FALSIFIED at 1480σ** |
| E[L⁻] = μ/(1−q) identity (Result 15 hypothesis) | **CONFIRMED** at 10⁻⁴ |
| q* = E_P[e^(−L⁻)] (Esscher relation) | **CONFIRMED** at 0.33σ |
| E[L⁻] in closed form via Esscher route | **NOT achieved** |
| Closed form for q via Esscher tilt | **NOT achieved** |

**Conclusion:** the Esscher-tilt route does NOT close E[L⁻]. The proposed duality identity (descending under P = ascending under P*) is empirically false. The structural identities involving Esscher are valid but they form a closed system with q and q* both unknowns.

The fundamental obstruction is the same as Path C: q (and equivalently 1 − q*) require summing P(S_n > 0) at irrational truncation indices ⌊n · log₂ 3⌋, which has no closed-form summation due to Gelfond-Schneider transcendence of log₂ 3.

**The Esscher-duality route is exhausted.** Per the user's brief decision criterion: "If no clean closed form emerges from the calculation: Esscher-duality route exhausted; E[L⁻] is genuinely not closed-form even with the Esscher machinery; Path B is the only remaining route." That is the verdict.

## 8. Path B (Result 17) is the next move

Why: Path B operates on the *Markov-modulated absorbing chain* on residues mod 2^k (algebraic, eigenvalues of rational-entry matrix Q), bypassing the iid characteristic-function transcendence that blocks Path C and now Path A's Esscher extension.

Path B is the operationally tractable route to W_j (and via aggregation to ε_S and ε(σ)). It requires:
1. Build Q on residues mod 2^k (k=6 → 32 odd states).
2. Compute matrix Wiener-Hopf factorization (Alsmeyer-Buckmann 2018).
3. Extract W_j from the factorization.

The matrix factorization handles the algebraic combination of transcendental scalar pieces cleanly; Q's eigenvalues are algebraic (roots of rational-entry characteristic polynomial), avoiding the log₂ 3 obstruction at the matrix level.

---

## Files

- `esscher_duality.py` — duality test (5M-orbit simulation, numba-accelerated)
- `path_c_q_attempt.py` — earlier focused q test (10M orbits, q = 0.71355)
- `esscher_duality_attempt.md` — this document
- `closed_form_findings.md` — Result 16 entry below

## Citations

- Spitzer 1956 — Spitzer-Andersen formula for q (1 − q = exp of partial-sum series)
- Kyprianou 2014 Ch 6 — Wiener-Hopf factorization, Esscher transform
- Feller Vol II Ch XII — duality lemma (time reversal, NOT Esscher tilt)
- Asmussen Ch VIII — Esscher transform for random walks, Cramer-Lundberg framework

## Honest scope statement

The user's brief proposed Feller XII as the source of the descending-vs-ascending duality identity. After re-checking: Feller's duality is via time reversal, which preserves step distribution for iid walks but is structurally different from Esscher tilt. The user's brief identified the wrong duality form. The correct Feller XII identity is the *path-level* duality for max/min of the same walk, not a P → P* relationship. The Esscher tilt is its own machinery (Asmussen Ch VIII) and does not provide the claimed equality of ladder means across measures.

This is not a "this didn't work" result — it's a "the proposed identity doesn't hold" result. The Esscher-tilt machinery is internally consistent (the three identities I confirmed are real), but it doesn't pin q and hence doesn't close E[L⁻].
