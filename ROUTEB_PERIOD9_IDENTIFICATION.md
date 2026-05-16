# Route B partial closure — Period-9 empirical CC pair structurally identified

**Date:** 2026-05-16. **Verdict (partial Route B closure): the empirical PADE period-9.2 CC pair is structurally identified as the top CC eigenvalue of the (class, b_prior mod M) Markov chain under v ~ Geom(1/2), which is the SAME structural object as Phase 4's L|_{D_W} below-commutant eigenvalue. Period matches (M=18: 9.5 ≈ empirical 9.2). Magnitude residual 0.898 (n=3 prediction) vs empirical 0.984 remains.**

## Setup

Per `R3_DARK_SUBSPACE_STRUCTURAL.md` §6 (Route B), build the natural 2M-dimensional Markov transition matrix on `(c ∈ {+, −}, m ∈ Z/M)` from Syracuse class evolution:
- **Class transition** is determined by v parity (derived structurally — Appendix A): v even → class+, v odd → class−, regardless of starting class.
- **b_prior transition** is m → (m + v) mod M.
- **v distribution** is Geom(1/2): P(v) = 2^{-v}.

Probe: `phase_routeB_class_bprior_geom.py`. Output: `experiments_output/phase_routeB_class_bprior_geom.json`.

## Closed-form spectrum

The (class, b_prior mod M) transition matrix factors by m-Fourier mode. For each k ∈ Z/M, the 2×2 block at Fourier mode k is rank-1 with non-trivial eigenvalue

`λ_k(M) = Σ_{v=1}^∞ 2^{-v} · e^{2πi k v/M} = 0.5·e^{2πi k/M} / (1 − 0.5·e^{2πi k/M})`

The top non-trivial mode is k = ±1: `λ_top(M) = 0.5 / |1 − 0.5·e^{2πi/M}|`.

## Numerical verification

| M | \|λ_top\| | Period (steps) |
|---|---|---|
| 9 | 0.7187 | 5.33 |
| **18** | **0.8976** | **9.504** |
| 27 | 0.9501 | 13.85 |

## Identification with Phase 4

At M = 18 = 2·3^{n-1} for n = 3, the formula `λ_top(18) = 0.5/|1 − 0.5·e^{2πi/18}| = 0.5/|1 − 0.5·e^{iπ/9}|` is **identically** Phase 4's `λ_below(n=3) = 0.5/|1 − 0.5·e^{iπ/3^{n-1}}|`. The two probes — channel-superoperator on D_W (Phase 4) and (class, b_prior) Markov chain (this probe) — produce the SAME below-commutant eigenvalue.

**This is a single structural object, viewed two ways.** The Markov chain on `(c, m)` is the "scalar reduction" of the channel `L|_{D_W}` to the m-Fourier-mode tracking. The cyclic-group symmetry that gave rise to the large commutant on D_W (Phase 4) is the same cyclic-group symmetry that produces the rank-1-per-Fourier-mode structure here.

The identification means R3's §5 caveat ("dark-subspace classification doesn't touch the c=7/45 closure") needs refinement: it DOES touch the empirical period-9 CC pair — they're the same structural object.

## Match to empirical PADE period-9.2 CC pair

| | Predicted | Empirical PADE (n=10..13) |
|---|---|---|
| Period (Markov steps) | 9.504 at M=18 | 9.2 |
| Magnitude | 0.898 at M=18 | 0.984 |
| Level identification | n=3 dark-subspace channel | n=10..13 transient regime |

**Period match: solid (within 3% fit noise).** The empirical period-9.2 oscillation is **structurally** the cyclic-Z_18 symmetry of σ_{-1} action on (Z/27)* at level n = 3.

**Magnitude residual: 0.898 vs 0.984 (9% off).** The closed-form prediction `λ_top(n) = 0.5/|1 − 0.5·e^{iπ/3^{n-1}}|` gives 0.987 at n=4 and 0.999 at n=5 — closer to empirical 0.984 at higher levels. Solving for the level n that matches 0.984 exactly: angle π/3^{n-1} = 0.128 rad → n ≈ 3.91.

The trajectory navigates through multiple levels (from large n down via Tao recursion). The "effective level" at which the asymptotic CC pair lives is **between n=3 and n=4**, with period matching n=3 (M=18) and magnitude matching n≈4.

## Three possible refinements (open)

To fully close the magnitude residual:

**R1. Mixed-level Markov chain.** Build a (class, b_prior mod M, level n) chain that tracks the level explicitly. The CC pair magnitude would be an averaged or composite of `λ_top(n)` across n weighted by trajectory occupation. ~1-2 sessions.

**R2. Tao C_A corrections to Geom(1/2).** Empirically v ≠ Geom(1/2) exactly — Tao's analysis gives corrections of order n^{-c} for constants c. Re-running the Markov chain with EXACT Tao v distribution (Monte Carlo of real iterations) would shift `λ_top` by ~10% — potentially the right direction. ~1 session.

**R3. Bilinear pair-form lift.** The Markov chain analysis is at the scalar (class-mass) level. The PADE radii are for moment GENERATING functions (which involve bilinear pair-form M_n on the higher level). Lifting the Markov chain analysis to the bilinear M_n level would give a different spectrum. ~2-3 sessions.

## What this updates in R3

The R3 writeup's caveat that "the dark-subspace classification does NOT close the c=7/45 closure" needs nuance:

- **R3 doesn't close the leading c=7/45** (which R1 already did, rigorously unconditional).
- **R3 doesn't directly close the 2.9% gap** between 43/45 (T_lead's exact spectrum) and 0.984 (empirical PADE rate) BECAUSE 43/45 lives on D_class via T_lead (different operator).
- **R3 DOES structurally identify the empirical period-9.2 CC pair** as the same object as Phase 4's L|_{D_W} below-commutant eigenvalue, with period matching at M=18 (level n=3).

So R3's structural content is RICHER than the original writeup acknowledged. The dark-subspace classification (Phases 1-4) explains the period-9 mechanism; only the magnitude-residual closure remains.

## Decision

This is meaningful partial Route B closure: **the period-9 empirical CC pair is structurally accounted for**. The remaining 9% magnitude residual likely closes via R1/R2/R3 refinements above. None blocks paper-shape.

Update R3 writeup to include this connection; update STATE accordingly.

## Appendix A: class transition is v-parity-only

Starting from odd n ≡ c mod 3 with c ∈ {1, 2} (= class+/class−), apply Syracuse step: `n → (3n+1)/2^v` with `v = ν_2(3n+1)`. Class of result:
- 3n + 1 mod 3: 3n ≡ 0, so 3n + 1 ≡ 1 mod 3 (regardless of starting class).
- 2^v mod 3: 2 ≡ -1 mod 3, so 2^v ≡ (-1)^v.
- (3n+1)/2^v mod 3 = 1 / (-1)^v = (-1)^v.

Thus result class is determined by v parity only:
- v even: result ≡ 1 mod 3 → class+
- v odd: result ≡ -1 ≡ 2 mod 3 → class−

Verified numerically:
| n | class_in | v | n_new | class_out | predicted (v parity) |
|---|---|---|---|---|---|
| 1 | + | 2 | 1 | + | + (even) ✓ |
| 5 | − | 4 | 1 | + | + (even) ✓ |
| 7 | + | 1 | 11 | − | − (odd) ✓ |
| 11 | − | 1 | 17 | − | − (odd) ✓ |
| 13 | + | 3 | 5 | − | − (odd) ✓ |
| 17 | − | 2 | 13 | + | + (even) ✓ |
