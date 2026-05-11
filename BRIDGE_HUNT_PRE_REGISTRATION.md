# BRIDGE HUNT — Pre-registration

**Timestamp:** 2026-05-11
**Working directory:** C:/Collatz/
**Parent commit:** d41899c (post F̂_p verification)
**Concurrency note:** parallel agent writing TAO_*/BOOKKEEPING_* files; this hunt writes ONLY BRIDGE_* files.

## Question

Does there exist an explicit math relationship Φ between
- **F̂_p^full(ξ)** — full-period Fourier transform of f_p(u) = e_M(c·(1+p)^u) on Z/M, M = p^{r+1} (deterministic cyclic-group character sum, family-level theorem verified at 33 cells)
- **μ̂_n(ξ)** — characteristic function of Syrac(Z/3^n), the Markov-chain stationary on (Z/3^n)\* (Tao Prop 1.17's object)

such that a bound on |F̂_p(ξ)| transfers to a bound on |μ̂_n(ξ')| via Φ, with the transferred bound tight enough to feed Nisoli's η < 1 requirement in R77.2?

The QX1 Move 2 attempt (`QX1_FAMILY_THEOREM_ATTEMPT.md` §Phase 4) flagged this transfer as missing. The verified F̂_p theorem at family level (`FHAT_THEOREM_VERIFICATION_RESULTS.md`) does NOT by itself dissolve R77.2's reliance on Tao Prop 1.17. The present hunt is for that missing bridge.

## Hypotheses

- **H_BRIDGE (alternative):** Φ exists; |F̂_p(ξ)| ≤ B(p,r,ξ) implies |μ̂_n(ξ')| ≤ B'(n,ξ') via an explicit math derivation, with B' tight enough that Nisoli's η = C_A · M_3 condition fires for c = 7/45.

- **H_NULL (default, favored):** F̂_p (deterministic, cyclic-group, scale Z/p^{r+1}) and μ̂_n (probabilistic-expectation over Markov-chain stationary, scale Z/3^n) are structurally different objects; the only way to bound μ̂_n with F̂_p data is to recover Tao §7.2-7.4 (which loops back to the R77.2 status quo). The Move 2 §Phase 4 obstruction was correct.

**Pre-registered favoring NULL.** Override requires:
1. **Empirical agreement** at parameters tested in Phase 3, AND
2. **Structural derivation** of Φ from established math, AND
3. **Out-of-sample verification** at parameters NOT used to construct Φ.

A single one or two of these without the third is insufficient.

## Procedure (locked)

### Phase 1: Structural articulation → `BRIDGE_STRUCTURAL_MAP.md`

Enumerate structural commonalities and divergences between F̂_p and μ̂_n with technical precision. Sub-tasks:

1.1. List shared arithmetic substrate (groups, characters, lifting from Z/p^r to Z/p^{r+1}, conservation laws, principal-unit structure, etc.).
1.2. List divergences (deterministic vs expectation; cyclic-group character on Z/p^{r+1} vs Z/3^n; period p^r vs Markov-chain stationary measure; etc.).
1.3. Tabulate against Tao 2022 §1.4-1.5, R76, R77.2, R78 to anchor structural claims.

### Phase 2: Candidate bridges → `BRIDGE_CANDIDATES.md`

Enumerate candidate maps Φ in the categories:
- (A) **Spectral conjugacy** via R76 §10 / R77.2 T_diag / R77 T_3.
- (B) **Averaging/iteration** of F̂_p over c, r, ξ to mimic μ̂_n's expectation.
- (C) **Coupling** via parity-sequence / p-adic lifting.
- (D) **Parameter specialization** F̂_p at p=3 vs μ̂_n at corresponding n.
- (E) **No-bridge** baseline (objects parallel but not directly linked).

Each candidate requires: precise math statement, empirical test design at small parameters, falsifiable criterion.

### Phase 3: Empirical testing → `BRIDGE_EMPIRICAL_TESTS.md` + scripts `bridge_*.py`

For each non-trivial candidate from Phase 2, run a numerical check at small p (focus p=3) and small n.

Compute μ̂_n empirically using Tao 2022 (1.26):
> Syrac(Z/3^n) ≡ Σ_{i=1..n} 3^{i-1} · 2^{-a_{[1,i]}}  (mod 3^n)

with iid Geom(2) a_i (independent over i = 1..n). Estimate
> μ̂_n(ξ) ≈ (1/N_samp) Σ_{j=1}^{N_samp} exp(-2πi · ξ · S_j / 3^n)

over N_samp ≥ 200000 samples per (n, ξ) cell.

Per candidate, three possible outcomes:
- AGREES (within Monte-Carlo precision + theory tolerance): advance to Phase 4.
- DISAGREES (cleanly outside MC error): candidate falsified.
- AMBIGUOUS: document; reassess if any candidate else clears.

### Phase 4: Structural derivation → `BRIDGE_DERIVATION.md` (only if candidate survives Phase 3)

Required for declaring BRIDGE_FOUND_RIGOROUS:
- Explicit derivation of Φ from established math (cite specific theorems and project results).
- Explicit error bounds B'(n,ξ') as function of B(p,r,ξ).
- Verification at parameters NOT used in Phase 3 (out-of-sample, A4).
- The bridge MUST NOT loop through Tao §7.2-7.4 implicitly (A2 tautology check).

### Phase 5: Disposition → `BRIDGE_HUNT_DISPOSITION.md`

Document the disposition, with summary at top and per-candidate results.

## Decision rules (locked, no hedging at disposition time)

- **BRIDGE_FOUND_RIGOROUS:** survived Phase 3 + Phase 4 + A4 out-of-sample. Mechanism explicit, not tautological. Route 3 viable for c=7/45 closure.
- **BRIDGE_FOUND_EMPIRICAL:** survived Phase 3 empirically, Phase 4 incomplete or unable to derive mechanism. Documented as empirical, NOT claimed as bridge.
- **BRIDGE_PARTIAL:** holds in restricted regime only (e.g. r=1 boundary, or specific frequency subset). Documented with scope.
- **NO_BRIDGE_FOUND:** all candidates falsified at Phase 3 or fail Phase 4 derivation. Document obstructions precisely.
- **INCONCLUSIVE:** small-parameter data insufficient; what computation would resolve.

**Early termination:** if Phase 1 reveals objects too structurally different for any promising candidate, disposition lands at NO_BRIDGE_FOUND without exhaustive Phase 2-4 work. Pre-committed to honest early termination.

## Adversarial safeguards

- **A1 — correlation ≠ structure:** numerical match without derivation = noise that aligns. Marked EMPIRICAL_PATTERN_NO_DERIVATION.
- **A2 — anti-tautology:** Φ must improve on Tao Prop 1.17, not restate it. If Phase 4 derivation routes through Tao §7.2-7.4 implicitly, it is bookkeeping in disguise; downgrade to no-bridge.
- **A3 — honest scope:** if the derivation requires expertise the project lacks (e.g. effective Burgess for principal-unit cosets at family level), flag explicitly. No speculative derivations.
- **A4 — out-of-sample required:** BRIDGE_FOUND_* requires at minimum one cell not used to fit/calibrate Φ.

## Files this hunt will write (BRIDGE_ prefix only)

- `BRIDGE_HUNT_PRE_REGISTRATION.md` (this file, committed FIRST in its own commit)
- `BRIDGE_STRUCTURAL_MAP.md` (Phase 1)
- `BRIDGE_CANDIDATES.md` (Phase 2)
- `BRIDGE_EMPIRICAL_TESTS.md` (Phase 3) + scripts `bridge_*.py`
- `BRIDGE_DERIVATION.md` (Phase 4, conditional)
- `BRIDGE_HUNT_DISPOSITION.md` (Phase 5, disposition at top)

## Commit pattern

- Pre-reg commit (this file ONLY): explicit `git add BRIDGE_HUNT_PRE_REGISTRATION.md` then commit.
- Later results commit: explicit `git add BRIDGE_*.md bridge_*.py` then commit.
- **Never** `git add .` or `git add -A` — parallel agent writes TAO_*/BOOKKEEPING_* concurrently.

## Constraints

- Pre-register before any compute fires (this file goes in BEFORE bridge_*.py runs).
- Exploratory research. Bridge declared only when empirical + structural + out-of-sample all align.
- Honest negative results constrain the search and inform Route 1 (Tao bookkeeping) vs Route 2 (Burgess) choice for c=7/45 closure.

---

End pre-registration.
