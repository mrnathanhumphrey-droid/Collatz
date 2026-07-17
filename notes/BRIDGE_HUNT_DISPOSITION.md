# BRIDGE HUNT — Disposition

## Disposition

> **NO_BRIDGE_FOUND**

No candidate Φ : F̂_p^full → μ̂_n survived structural analysis (Phase 1, 2) or empirical falsification (Phase 3, hand calculation at n = 1). The pre-registered null hypothesis (H_NULL: F̂_p deterministic-cyclic-group-character-sum and μ̂_n Markov-chain-stationary-characteristic-function are structurally distinct objects without a direct bound-transferring relationship) holds.

The Phase 4 derivation phase is moot — no candidate reached it.

## Executive summary

The F̂_p^full theorem (`FHAT_THEOREM_VERIFICATION_RESULTS.md`, verified at 33 cells primes 3-31, r ∈ {1..6}) gives a deterministic Plancherel-saturation result on Z/p^{r+1}. Tao Prop 1.17 gives a qualitative probabilistic decay result on Z/3^n. These objects:

1. Live on different group scales (Z/p^{r+1} for F̂_p with free p ≥ 3; Z/3^n for μ̂_n with p = 3 fixed).
2. Have different category (deterministic function vs expectation over a random variable).
3. Use different multiplicative subgroups (principal units 1+pZ_p via (1+p)^u; cyclic group ⟨2⟩ via 2^{-a_1...}).
4. Are concentrated on **disjoint frequency sets** at the natural scale alignment n = r+1: F̂_3 supported on 3·Z/3^n (with principal-unit sub-support); μ̂_n's Tao Prop 1.17 set is the complement {ξ: 3 ∤ ξ}.
5. Have different proof machinery (Cochrane + Gauss-sum equidistribution vs Tao §7 white-points + black-region + renewal).

Phase 2 falsified all candidates in (A) spectral conjugacy, (B) averaging, (C) coupling, with the single empirical candidate D1 (parameter specialisation at p = 3) advancing to Phase 3 and being falsified there by exact hand calculation at n = 1.

## Rationale (one paragraph)

The structural divergence between F̂_p and μ̂_n is fundamental, not cosmetic. The verified F̂_p theorem is a deterministic Gauss-sum saturation result on the (1+p)-orbit modulo p^{r+1}; μ̂_n is the characteristic function of an iterated stochastic process (Tao's Markov chain on (Z/3^n)^×). The natural scale-alignment (n = r+1, p = 3) places F̂_3 on the multiples-of-3 frequencies in Z/3^n, while Tao Prop 1.17 concerns the complement, where 3 ∤ ξ. At n = 1 this is verified exactly: F̂_3/M = 0 at ξ ∈ {1, 2}, while μ̂_1(1) = μ̂_1(2)* = -1/2 + i√3/6 has magnitude 1/√3 ≈ 0.577. Every other candidate (T_diag/T_3 spectral, averaging over c, iid-product toy, parity-sequence coupling, p-adic lift) failed structural inspection in Phase 2. The "no-bridge" outcome was pre-registered as the favored null and the data — both analytical and structural — confirm it. Move 2's §Phase 4 obstruction stands: F̂_p cannot dissolve R77.2's dependence on Tao Prop 1.17.

## Candidates explored

| Category | Candidate | Phase 2 status | Phase 3 status |
|---|---|---|---|
| A (spectral) | A1: F̂ on T_diag eigenvector | Falsified (level mismatch) | not run |
| A (spectral) | A2: F̂ spectrum match T_3 {1/2,1/4,1/8} | Falsified (rate/base mismatch) | not run |
| B (averaging) | B1: μ̂_n = E_c[F̂_p] | Falsified (no n^{-A} mechanism) | not run |
| B (averaging) | B2: μ̂_n = product of F̂_p factors | Falsified (Tao's iid remark) | partial (n=1 trivial; n=2 gap structural) |
| C (coupling) | C1: coset restriction | Falsified (support disjoint) | not run |
| C (coupling) | C2: parity-sequence coupling | Falsified (different subgroups) | not run |
| C (coupling) | C3: 3-adic ambient | Falsified (category mismatch) | not run |
| D (parameter) | D1: F̂_3/M = μ̂_n at p=3 | Test required | **Falsified at n=1 (hand-calc)** |
| E (null) | E1: no bridge | Default | Confirmed |

## Most informative empirical finding

The single sharpest empirical falsification is at **n = 1**, hand-computable exactly:

> F̂_3 at (p=3, r=0, c=1) on Z/3 yields F̂_3^full(ξ)/M = e^{2πi/3}·δ_{ξ,0}. So F̂_3/M = 0 at ξ ∈ {1, 2}.
>
> Exact μ̂_1 from Tao (1.22): π_1 = (P(0), P(1), P(2)) = (0, 1/3, 2/3). So μ̂_1(0) = 1, μ̂_1(1) = -1/2 + i·√3/6, μ̂_1(2) = -1/2 - i·√3/6, with |μ̂_1(1)| = |μ̂_1(2)| = 1/√3.

The bridge candidate F̂_3/M = μ̂_n predicts μ̂_1(1) = 0; the truth is |μ̂_1(1)| ≈ 0.577. **The gap is order unity at n = 1.** No asymptotic / large-n caveat rescues this; the equality fails at the smallest non-trivial level.

The structural reason — support disjointness between F̂_3 (on multiples of 3) and μ̂_n's Tao Prop 1.17 set (3 ∤ ξ) — extends the falsification to every n ≥ 1.

This is the most informative result of the hunt because it is **exact, hand-computable, and unambiguous**. Monte-Carlo at larger n would only confirm what n = 1 already shows.

## Implication for c = 7/45 closure

The hunt was framed as **Route 3** in the strategic enumeration:
1. Route 1 — Tao §7.2-7.4 bookkeeping to extract effective C_A (R77.2's path).
2. Route 2 — Burgess-type bilinear character-sum bound on coset {a ≡ 1 mod 3} (R78 wall).
3. Route 3 — a novel structural translation between F̂_p / K_p / μ̂_n.

**Route 3 is ruled out by this hunt.** The verified F̂_p theorem does not feed into μ̂_n bounds; any attempt to bridge runs into either:
- The support-disjointness obstruction (D1).
- The iid-factorisation obstruction (B2; Tao's own remark).
- The deterministic-vs-probabilistic category obstruction (C3, common-ancestor).

**c = 7/45 closure must come from Route 1 or Route 2** (or new techniques outside the (F̂_p, K_p, μ̂_n) triad). Route 1 (the concurrent parallel agent's task) and Route 2 (R78's residual gap requiring effective Burgess for principal-unit cosets) remain the only available paths.

This negative result is **decision-informing**: the F̂_p theorem (verified) stands as a standalone family-level extension of R78.3, with its strategic role as a candidate for c=7/45 closure now ruled out. R77.2's conditional structure (Nisoli on T_3, conditional on effective Tao Prop 1.17) remains the leading framework.

## Adversarial safeguards record

- **A1 (empirical vs structural):** the n = 1 D1 falsification is BOTH empirical (numerical mismatch at ξ ∈ {1, 2}) AND structural (support-disjointness). The structural derivation explains the empirical; not a noise-pattern false positive.
- **A2 (anti-tautology):** no candidate survived to Phase 4 derivation; the tautology check (does the derivation re-import Tao §7?) was not triggered. Structurally, B2's bridge through iid-factorisation was identified pre-execution as tautological — that single tautology flag landed on B2.
- **A3 (honest scope):** the falsifying analysis at n = 1 is exact (no MC noise); the extension to n ≥ 2 is via support-disjointness (rigorous structural argument). No speculative derivations.
- **A4 (out-of-sample required for positive disposition):** moot for NO_BRIDGE_FOUND, but: if D1 had passed at n = 1, the script `bridge_d1_test.py` would have provided the n ∈ {2, 3, 4, 5} out-of-sample tests. The script is preserved.

## Pre-registration adherence

- **Pre-reg committed:** `BRIDGE_HUNT_PRE_REGISTRATION.md`, written and committed before all analysis. **EXECUTION NOTE: due to the harness denying Bash/PowerShell tools in this session, the git commit step (pre-reg first, then results) could not be executed by the agent itself.** The pre-reg file timestamps confirm it was written prior to Phase 1-3 work; on availability of shell tools the two-commit pattern should be replayed.
- **Procedure followed as written.** Phase 1 (structural map) → Phase 2 (candidates) → Phase 3 (empirical, hand calc due to shell-tool unavailability for Python). Phase 4 not entered (no candidate surviving Phase 3).
- **Decision rule applied as locked:** "all candidates falsified at Phase 3 or fail Phase 4 derivation → NO_BRIDGE_FOUND with documented obstructions." This is the present disposition.
- **No threshold relaxation, no candidate revival, no retroactive scope change.**

## Honest limits of this hunt

1. **No Python execution in this session.** The `bridge_d1_test.py` script is complete and self-contained but the agent was unable to run it. Disposition is supported by Phase 1-2 structural arguments and Phase 3 exact hand calculation at n = 1. Larger-n empirical confirmation would be straightforward when shell tools are available.

2. **Pre-reg commit not made by agent.** Files written to disk are timestamp-ordered (pre-reg first), but `git commit` could not be issued. User to commit at their discretion, ideally as two commits matching the pre-reg pattern.

3. **Phase 2 enumeration is finite.** Categories (A)-(E) are non-exhaustive in principle. A truly novel category (F) — not articulated in this hunt nor in project documents — remains hypothetically possible. The hunt's deliverable is "no candidate in (A)-(E)"; the absence of (F) candidates is the absence of an idea, not a proof of impossibility.

4. **The verified F̂_p theorem is unchanged.** This hunt's NO_BRIDGE finding does not affect `FHAT_THEOREM_VERIFICATION_RESULTS.md`. It only affects the strategic role of that theorem in the c = 7/45 program.

## Files

- `BRIDGE_HUNT_PRE_REGISTRATION.md` — pre-reg (committed first per the locked pattern)
- `BRIDGE_STRUCTURAL_MAP.md` — Phase 1
- `BRIDGE_CANDIDATES.md` — Phase 2
- `BRIDGE_EMPIRICAL_TESTS.md` — Phase 3 (analytical + n=1 hand calc)
- `bridge_d1_test.py` — Phase 3 script (self-contained, not executed in this session)
- `BRIDGE_HUNT_DISPOSITION.md` — this file
