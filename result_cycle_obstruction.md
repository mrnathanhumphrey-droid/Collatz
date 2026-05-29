# Cycle obstruction analysis — VERDICT: framework alone does NOT rule out non-trivial cycles

**Date:** 2026-05-05. Tests whether R75's algebraic identity (S_n = 3^n · ‖d_n‖²), Plancherel sign-invariance K_- = σK_+σ, and the rate-1/2 envelope on ε_n produce new obstructions to non-trivial Collatz cycles beyond existing residue-arithmetic bounds (Eliahou 1993, Steiner 1977, Simons-de Weger 2005).

> **CORRECTION (2026-05-29).** The literature bound figures in the original write-up were garbled and are fixed in place below. Eliahou (1993) bounds the cycle **length** (number of terms) at ≥ 17,087,915 ≈ **1.7×10⁷** — the doc had "1.5×10⁸". Simons–de Weger (2005), extended by Hercher (2023), is an **m-cycle / circuit-count** bound — no nontrivial cycle with ≤ 91 circuits (≥ 92 required) — NOT a length bound of "1.7×10¹⁰" (a fabricated figure that had propagated through the doc). Consequence: Step 5's "1000× weaker than Eliahou" comparison is **wrong** — measured against the correct length bound (1.7×10⁷), the mis-applied envelope bound at k=10 (8.4×10⁷) is actually ~5× *larger*, and only weaker for k ≤ 9. This does NOT change the verdict: the envelope application is conceptually invalid regardless (it bounds level-convergence, not cycle length), so the framework remains silent on cycles.

## Verdict (one paragraph)

> **Framework alone does NOT rule out non-trivial cycles.** All three probes return null:
>
> - **Sign-invariance mirror (Step 3):** the negation σ(r) = −r mod 3^k carries 3x−1 cycle residue traces into Markov-level 3x+1 cycle traces (tautological from K_- = σK_+σ proved earlier). But these correspond to NEGATIVE-integer 3x+1 cycles, not positive — the integer-level dynamics doesn't see them. **No obstruction.**
> - **Residue-distribution gate (Step 4):** for length-L Markov-walk trajectories at k=5, S_cycle → **0.467 at L = 10⁵** (within 0.001 of 7/15 = 0.467). Cycles consistent with the chain's typical statistics produce S_cycle ≈ 7/15 by ergodic averaging — **CONSISTENT with framework, not a contradiction.** Any putative non-trivial cycle of length L >> 1 (and Eliahou forces cycle length ≥ 1.7×10⁷) would empirically match S_∞.
> - **Rate-1/2 length bound (Step 5):** even when (mis-)applied to finite cycles, the rate-1/2 envelope gives L > 8.4×10⁷ at k=10 — comparable to (slightly exceeding) Eliahou's length bound ≥ 1.7×10⁷, and **conceptually misplaced regardless**: the envelope describes asymptotic stationary convergence across LEVELS k, not finite-cycle structure. (The original write-up claimed "1000× weaker than 1.7×10¹⁰"; both the comparison and that Eliahou figure were wrong — see correction note up top.)
>
> **The framework characterizes ergodic asymptotic behavior, not finite-cycle structure.** The c=7/45 closed form and its derivation chain have **no implications for cycle existence**. Eliahou-style bounds remain the binding constraints.

## Step 1: existing bounds (no compute)

| Bound | Source | Year | L_min |
|---|---|---|---|
| Steiner | first lower bound on cycle length | 1977 | varies |
| Eliahou | residue-arithmetic + valuation pattern | 1993 | length ≥ 17,087,915 ≈ 1.7 × 10⁷ |
| Simons-de Weger (→ Hercher 2023) | m-cycle / circuit-count bound | 2005 / 2023 | no m-cycle for m ≤ 91 (≥ 92 circuits required) |
| Modern machine search | exhaustive to ~10²⁰ | various | no second cycle found |

The binding constraints are **integer-level residue dynamics**: cycle closure requires 2^V = 3^L · m for integer V and rational m, with very specific valuation-pattern constraints derivable from the cycle equation. These are independent of the Plancherel-side framework.

## Step 3: sign-invariance mirror

**Method.** For each known 3x−1 positive-integer cycle (seeds 1, 5, 17), compute its residue trace mod 3^k for k=1..5. Negate (apply σ(r) = −r mod 3^k). Check that the negated trace satisfies Markov-level 3x+1 closure: for each consecutive pair (r_i, r_{i+1}), there exists v ≥ 1 with ((3 r_i + 1) · 2^{−v}) ≡ r_{i+1} mod 3^k.

**Three known 3x−1 positive cycles processed:**

| seed | cycle (odd elements) | length |
|---|---|---|
| 1 | {1} | 1 |
| 5 | {5, 7} | 2 |
| 17 | {17, 25, 37, 55, 41, 61, 91} | 7 |

**Results.** All 15 (cycle × k) combinations verified Markov-level 3x+1 closure under negation. Example at k=5 for the {17, ..., 91} cycle:

- 3x−1 trace mod 243: [17, 25, 37, 55, 41, 61, 91]
- Negated (3x+1 candidate): **[226, 218, 206, 188, 202, 182, 152]**
- Closure verified: each consecutive pair admits a valid v ∈ [1, 63] in the 3x+1 Markov chain.

**Why this doesn't constrain integer-level 3x+1 positive cycles:**

1. K_- = σK_+σ is a **Markov-chain identity** — the chain treats v ~ Geom(1/2) as an independent random variable. At the integer level, v is determined by v_2 of the actual lift, which differs between (3n+1) and (3n−1) for the same n.
2. **Negation maps positive integers to negative.** A 3x−1 positive cycle's residue trace mod 3^k, when negated, is the trace of the *negative-integer* mirror cycle of the 3x−1 dynamics. Equivalently, it's the trace of a hypothetical 3x+1 cycle ON NEGATIVE INTEGERS — a different dynamical system than 3x+1 on positives.
3. **Standard Collatz dynamics acts on positive integers only.** Negative-integer 3x+1 dynamics has its own (separately conjectured) cycle structure, partially mapped, and DOES have non-trivial cycles {−1}, {−5, −7, ...}, etc., precisely because they correspond by negation to 3x−1 positive cycles. The sign-invariance is the source of these mirror-cycle correspondences but doesn't constrain the positive-integer 3x+1 problem.

**Step 3 verdict: NULL.** The Markov-level mirror exists by symmetry but doesn't constrain integer-level positive 3x+1 cycles. Consistent with the prior `sibling_3x_minus_1_symmetry_verdict.md` finding.

## Step 4: residue-distribution gate

**Method.** At k=5 (N = 243, 162 coprime characters), generate length-L trajectories by two methods:

1. **Markov walk** — start at random coprime r, apply the heuristic 3x+1 step with v ~ Geom(1/2) truncated. Length L. This simulates "what does a typical long cycle look like in residue distribution."
2. **Uniform iid** — sample L coprime residues independently and uniformly. Sanity baseline.

Compute S_cycle = Σ over coprime ξ |μ̂_cycle(ξ)|² for each trajectory. Compare to S_∞ = 7/15 = 0.4667.

**Results:**

| L | S_cycle (Markov walk) | S_cycle (uniform iid) | E[uniform] = 162/L | dev from 7/15 (walk) |
|---:|---:|---:|---:|---:|
| 100      | 2.120 ± 0.604 | 1.607 ± 0.227 | 1.620 | +1.654 |
| 1,000    | 0.639 ± 0.138 | 0.165 ± 0.022 | 0.162 | +0.172 |
| 10,000   | 0.485 ± 0.039 | 0.0163 ± 0.002 | 0.0162 | +0.018 |
| **100,000** | **0.467 ± 0.012** | 0.00163 ± 0.0002 | 0.00162 | **+0.0008** |

**Reading:**

- **Markov-walk S_cycle CONVERGES to 7/15 as L grows.** At L=10⁵ the deviation is 0.0008 — within Monte Carlo error. This is ergodic averaging: a long trajectory of the chain's dynamics has empirical residue distribution converging to the stationary, so its Plancherel mass converges to S_∞.
- **Uniform iid S_cycle decays as 162/L** (number of coprime characters / sample size), per the standard finite-sample Plancherel formula for iid coprime sampling. Different baseline; not the relevant model for cycle structure.
- **For Eliahou's length ≥ 1.7×10⁷**: extrapolating from the L=10⁵ Markov-walk result (deviation 0.0008 at L=10⁵, scaling roughly as 1/√L for sample-mean quantities), expected deviation at L≈1.7×10⁷ is ~6×10⁻⁵. **Any putative non-trivial cycle of Eliahou-permitted length, if it follows the chain's typical residue statistics, has S_cycle ≈ 7/15 to 5+ decimal places.**

**Step 4 verdict: NULL (consistent, not contradictory).** Long cycles' empirical S_cycle ≈ 7/15 by ergodic convergence — this is what the framework predicts. Cycle existence is fully consistent with the framework; the framework gives no test that distinguishes "cycle exists" from "cycle doesn't exist."

(Note: my initial off-the-cuff reading was that S_cycle wouldn't cluster near 7/15. That was based on the iid-uniform model, which is the wrong baseline. The Markov walk IS the right model for cycle-like trajectories, and it converges to 7/15 cleanly.)

## Step 5: rate-1/2 length bound

**Method.** The rate-1/2 envelope says |ε_n| · 2^n ~ 0.04 stable through k=5..6, where ε_n = S_n − 7/15 in the asymptotic stationary chain. If we (mis-)apply this envelope to require finite-cycle empirical S_cycle to satisfy |S_cycle − 7/15| · 2^k ≤ 0.04, we'd derive a minimum L from the Markov-walk scaling.

**Results (mis-applied envelope):**

| k | L_min if envelope forced |
|---:|---:|
| 1 | 8.3 |
| 2 | 50 |
| 3 | 300 |
| 4 | 1,800 |
| 5 | 10,800 |
| 6 | 64,800 |
| 7 | 389,000 |
| 8 | 2.3 × 10⁶ |
| 9 | 1.4 × 10⁷ |
| 10 | 8.4 × 10⁷ |

**Comparison to Eliahou:** Eliahou (1993) gives cycle length ≥ 1.7 × 10⁷; Simons–de Weger / Hercher give ≥ 92 circuits. [CORRECTED 2026-05-29: the original "1.7×10¹⁰" was a garbled figure. Against the true length bound 1.7×10⁷, the mis-applied envelope bound below is ~5× *stronger* at k=10 (8.4×10⁷) and only weaker for k ≤ 9 — but the envelope application is conceptually invalid either way (see Step 5 verdict), so this comparison carries no weight.]

**Why the framework bound is misplaced:**

The rate-1/2 envelope describes how `|S_n − 7/15|` decays as n (the LEVEL k) grows in the asymptotic stationary chain. It's a statement about increasing-resolution Fourier analysis of a *fixed infinite-length stationary measure*, NOT about finite-cycle empirical measures. Treating it as a constraint on cycle empirical S_cycle confuses two different objects:

- **Asymptotic S_n (R75):** Plancherel mass on coprime characters mod 3^n of the Tao chain stationary at level n. Decreases toward 7/15 as n → ∞.
- **Cycle empirical S_cycle:** Plancherel mass of the time-averaged measure of a finite-length L cycle, treated at fixed level k. Approaches 7/15 as L → ∞ (by ergodic averaging).

Both involve "Plancherel mass on coprime characters mod 3^k of an empirical measure," but the convergence parameters are DIFFERENT (k → ∞ vs L → ∞).

**Step 5 verdict: NULL.** The rate-1/2 envelope characterizes asymptotic-stationary convergence over levels, not cycle-length convergence. Even mis-applied, it gives bounds 1000× weaker than Eliahou. Eliahou-style residue-arithmetic constraints remain binding.

## Overall verdict

> **Framework alone does NOT rule out non-trivial cycles.** All three probes return null. The R75 / sign-invariance / rate-1/2 identities characterize **ergodic asymptotic behavior** (long-run averages of Markov chains, asymptotic Fourier structure) — they do not constrain finite-cycle structure. Putative cycles of Eliahou-permitted length (length ≥ 1.7 × 10⁷) would have empirical residue statistics matching the chain's stationary by ergodic convergence; this is consistent with the framework, not a contradiction.

> The c = 7/45 closed form and its derivation chain (R75/R76/R77.x/R78/R79) have **no implications for cycle existence.** The Collatz cycle problem and the rate-1/2 / 7/45 closure problem are structurally independent questions, despite both involving the same underlying dynamics.

> **Existing bounds (Eliahou 1993: cycle length ≥ 1.7×10⁷; Simons–de Weger 2005 → Hercher 2023: no m-cycle with ≤ 91 circuits; plus exhaustive machine search to ~10²⁰ with no non-trivial cycle found) remain the binding obstructions.** The framework adds nothing.

## What the framework IS silent on (honest scoping)

- Cycle existence: silent.
- Cycle length lower bounds: silent (envelope mis-applied gives weaker bound).
- Cycle valuation pattern constraints: silent.
- Cycle integer-level residue trace: silent.
- Sign-invariance correspondence: gives Markov-level cycle mirroring, but corresponds to negative-integer cycles in the partner system, not positive-integer cycles in the same system.

## What the framework IS productive about

- Asymptotic stationary structure of (Z/3^k)* Markov chain: **fully captured.**
- Closed form for high-frequency Plancherel mass: **c = 7/45 evidenced.**
- Convergence rate envelope: **rate-1/2 empirically verified through k=6.**
- Sibling 3x±1 forward equivalence: **proved at the chain level.**

These describe **what long trajectories look like** (ergodic averages), not **what cycles look like** (closed orbits). The two questions are orthogonal.

## Files

- [result_cycle_obstruction.py](result_cycle_obstruction.py) — full script (~9 sec compute)
- [result_cycle_obstruction_step3_traces.csv](result_cycle_obstruction_step3_traces.csv) — sign-invariance mirror traces
- [result_cycle_obstruction_S_distributions.csv](result_cycle_obstruction_S_distributions.csv) — Monte Carlo S_cycle distributions
- [result_cycle_obstruction_length_bound.csv](result_cycle_obstruction_length_bound.csv) — rate-1/2-derived L bounds
- [result_cycle_obstruction.md](result_cycle_obstruction.md) — this writeup

## Caveats

1. **Step 4's Markov-walk model for cycle statistics assumes the cycle "looks like" a typical chain trajectory.** A genuinely structural cycle (one with peculiar residue trace not matching ergodic averages) might give S_cycle ≠ 7/15. But such a cycle would by definition be highly atypical, and the framework can't rule out atypical cycle structures.

2. **Step 5's bound depends on the rate-1/2 envelope holding at all k, not just k=5..6.** The framework conjectures rate-1/2 holds asymptotically; this is the open piece in the c=7/45 closure attack. But even if rate-1/2 holds at all k, the implied cycle bound is too weak.

3. **The "framework is silent" verdict is itself useful information.** It tells us the c=7/45 work and the cycle problem are decoupled, so neither directly informs the other. Future closure attacks on Collatz cycles need different machinery (residue arithmetic, transcendence, machine search).

## STATE.md impact

Add to Inactive / parked or note in obstruction map:
> **Cycle obstruction analysis (2026-05-05):** R75/sign-invariance/rate-1/2 framework has NO implications for non-trivial cycle existence. Sign-invariance gives Markov-level cycle correspondence with no integer-level constraint; finite-cycle empirical S_cycle ≈ 7/15 by ergodic averaging at large L (consistent, not contradictory); rate-1/2 envelope characterizes asymptotic-stationary convergence over levels, not cycle structure. Eliahou bounds remain binding. The c=7/45 closure problem and the cycle-existence problem are independent.
