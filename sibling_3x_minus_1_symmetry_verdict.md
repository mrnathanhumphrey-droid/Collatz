# 3x+1 ↔ 3x−1 conjugation symmetry on Z/3^k — VERDICT: PROVED

**Date:** 2026-05-04. Sibling-probe Task 1 (Wilson-prompted).

## Verdict

> **PROVED.** The Syracuse Markov chains for 3x+1 and 3x−1 are conjugate by the negation permutation σ(r) = −r mod 3^k. Algebraically: T_-(r, v) = −T_+(−r, v) for every (r, v) ∈ (Z/3^k)* × {1, …, M}. Computationally: K_-(r, s) = K_+(−r, −s) entry-by-entry, π_-(r) = π_+(−r), and S_n^{3x−1} = S_n^{3x+1} as exact rationals at all tested n.
>
> **Therefore S_∞^{3x−1} = 7/15** (sharing all evidence with the 3x+1 case), **c = 7/45 is automatic for the 3x−1 forward Markov chain by symmetry**, and the +1 vs −1 forward-direction sibling probe is **closed**. No new compute was required beyond the existing q=3 framework — the result is zero-cost.

## Algebraic proof (Markov-chain heuristic)

In Tao's q=3 Syracuse framework, the chain on coprime residues r ∈ (Z/3^k)* uses transitions

  r ↦ ((3r ± 1) · 2^{−v}) mod 3^k,   v ~ Geom(1/2) truncated to [1, M], M = ord_{3^k}(2).

The key fact: **v is treated as an independent geometric random variable, not derived from v_2 of any specific lift.** This means the chain transition is purely modular arithmetic in (Z/3^k, ×, +).

Define σ: (Z/3^k)* → (Z/3^k)* by σ(r) = −r (well-defined since gcd(−r, 3) = gcd(r, 3) = 1). Compute T_+(σ(r), v):

  3·(−r) + 1 = 1 − 3r = −(3r − 1)   in Z/3^k

so

  T_+(−r, v) = (1 − 3r) · 2^{−v} = −(3r − 1) · 2^{−v} = −T_-(r, v)   in Z/3^k.

Equivalently: **T_-(r, v) = −T_+(−r, v)** for every (r, v).

For the chain matrices indexed by (initial, target) = (r, s):

  K_-(r, s) = Σ_v Pr(v) · 1[T_-(r, v) = s]
            = Σ_v Pr(v) · 1[−T_+(−r, v) = s]
            = Σ_v Pr(v) · 1[T_+(−r, v) = −s]
            = K_+(−r, −s).

So K_- = σ K_+ σ as matrices indexed by coprime residues. The two chains are unitarily conjugate by the permutation σ.

**Consequences:**

1. Spectra coincide: spec(K_-) = spec(K_+).
2. Stationary distributions are related by π_-(r) = π_+(−r).
3. ‖π_-‖² = Σ_r π_-(r)² = Σ_r π_+(−r)² = Σ_r π_+(r)² = ‖π_+‖² (just relabeling).
4. X_k^{(−)} = 3^k · ‖π_-‖² = 3^k · ‖π_+‖² = X_k^{(+)}.
5. S_n^{(−)} = X_n^{(−)} − X_{n-1}^{(−)} = X_n^{(+)} − X_{n-1}^{(+)} = S_n^{(+)} as exact rationals.
6. All Plancherel-side derived quantities (Tao's S_∞ = 7/15, R74's c = 7/45, R76's M_n family, R77.x's ε_n envelope) are identical between the two systems.

## Computational cross-check

Built K_+ and K_- as exact rational matrices over Q at k = 1, 2, 3. Solved stationaries via Gauss elimination on Fractions. Verified all five identities entry-by-entry:

| k | N | M | states | K_-(r,s) ?= K_+(−r,−s) | π_-(r) ?= π_+(−r) | X_k^{(−)} | S_k^{(−)} | matches known S_k^{(+)} |
|---|---|---|---|---|---|---|---|---|
| 1 | 3  | 2  | 2  | PASS | PASS | 5/3 | 2/3 | ✓ |
| 2 | 9  | 6  | 6  | PASS | PASS | 15/7 | 10/21 | ✓ |
| 3 | 27 | 18 | 18 | PASS | PASS | 177005/67963 | 31370/67963 | ✓ |
| 4 | 81 | 54 | 54 | PASS | PASS | 6626070796594781675/2159281421340253987 | 143195649659456490/308468774477179141 | ✓ |

Total runtime: <1 second. Algebraic identity holds with zero numerical noise (exact-rational equality).

## Why this works (the v_2 issue, addressed)

The brief noted: |·|_2 is sign-invariant, so v_2(3r+1) = v_2(−(3r+1)) = v_2(−3r−1). But on Z/3^k Z the negation operation gives 3·(−r) + 1 = 1 − 3r, an integer DIFFERENT from −(3r+1) when r ≠ 0. So v_2 of the actual lift is not the same.

**This concern doesn't apply** to Tao's Markov chain heuristic, because the chain doesn't take v from the lift's actual 2-adic valuation. It samples v geometrically with probability 2^{−v}/Z_v, independent of r. The map (r, v) ↦ ((3r ± 1)·2^{−v}) mod 3^k is purely modular arithmetic, where 1 − 3r = −(3r − 1) holds tautologically.

If one were studying the **integer-level** Syracuse map (where v IS v_2 of the lift), the symmetry would be more subtle — different lifts of the same residue can give different v's. But the Markov-chain framework in use here doesn't feel that subtlety.

## Implication for the sibling study

**Forward direction (3x+1 vs 3x−1) is closed by symmetry.** All R75/R76/R77.x results transfer with the relabeling σ.

**The genuinely new sibling work is on the inverse-tree / (x+1)/3 side**, where:

- 3x−1 has at least three integer cycles ({1,2}, {5,7,10,14}, {17, …, 34}) — the integer-level trajectory measure depends on which basin of attraction one starts in.
- The inverse tree of (x+1)/3 from 1 explores only one of those basins.
- D_avg^{3x−1} (mod 32) starting from each cycle is a separate object; the multifractal class might be the same or different.
- The Plancherel-side symmetry above does NOT extend to inverse-tree integer-level measures (which involve the actual lifts and their v_2's, not the modular Markov chain).

So the next probe — **inverse-tree (x+1)/3 D_avg per basin** — is independent compute, not subsumed.

## Files

- [sibling_3x_minus_1_symmetry_check.py](sibling_3x_minus_1_symmetry_check.py) — script
- [sibling_3x_minus_1_symmetry_verdict.md](sibling_3x_minus_1_symmetry_verdict.md) — this writeup

## STATE.md impact

Add to closed-form lock-ins: **K_- = σ K_+ σ on Z/3^k for σ = negation, hence S_n^{3x−1} = S_n^{3x+1} exactly, hence S_∞^{3x−1} = 7/15 and c^{3x−1} = 7/45 by the same evidence chain.**

Add to obstruction map: **the rate-1/2 rigorous-closure problem for c = 7/45 is identical between the two systems**; nothing in the +1/−1 split distinguishes them at the Plancherel-Markov level.

Reframe sibling-study direction: forward direction is symmetric. Real new science is in the inverse-tree side and the multi-basin question.
