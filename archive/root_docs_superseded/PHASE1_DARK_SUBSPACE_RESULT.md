# Phase 1 — Dark-subspace commutant probe result

**Date:** 2026-05-15 (post-compact, immediately after DWM_DARK_SUBSPACE_ATTACK_PLAN.md). **Verdict: fully irreducible at finite n — predicted outcome confirmed.**

## Setup

Per `DWM_DARK_SUBSPACE_ATTACK_PLAN.md` §"Phase 1": probe `dark_subspace_probe.py` builds all distinct adaptive Kraus operators
`M_v^{(j, b_prior)} = 2^{-v/2} · e^{-2πi · ξ · x_j(b_prior) · 2^{-v}/3^n} · σ_{-v}` (with `x_j(b_prior) = 3^{2j-2}·2^{-b_prior} mod 3^n`) acting on H_n = L²((Z/3^n)*), and computes the dimension of the joint commutant `A' = {X : [X, M] = 0 for all M}` via SVD of the stacked commutator-linearization `T_M = M^T ⊗ I − I ⊗ M`.

Configuration: V_MAX=12, J_MAX=3, b_prior cycled through 0..36 (sufficient to enumerate all distinct x_phase = x_j(b_prior)·2^{-v} mod 3^n).

## Result

| n | dim H_n | # distinct Kraus operators | dim(A') | Verdict |
|---|---|---|---|---|
| 2 | 6 | 84 | **1** | fully irreducible: A' = C·I |
| 3 | 18 | 252 | **1** | fully irreducible: A' = C·I |

SVD spectrum at both n: one isolated singular value ≈ 10^{-14..-15} (machine epsilon, the identity-line null vector), the next-smallest singular value is O(1) and well-separated. Single null direction = trivially the multiples-of-identity line.

Output JSON: `C:/Collatz/experiments_output/dark_subspace_probe.json`.

## Interpretation

**Syracuse's level-graded DWM trajectory is fully irreducible at finite n.** No standard Benoist-Pellegrini dark subspaces exist at finite truncation — the Kraus family acts irreducibly on the entire coprime-residue space.

This is the **predicted outcome** per the attack plan and is fully consistent with:
- **K_k spectrum** {1, 0, ..., 0} (K_STRUCTURE_RESULT.md): K_k mixes in EXACTLY k Markov steps, no invariant sub-rotations.
- **U_n → W_n exact** (INTERLEVEL_U_PROBE_RESULT.md): inter-level operator has no discrete-eigenvalue sub-structure.
- All five evening probes confirming continuous-spectrum reading at finite truncation.

The c=7/45 subdominant rate cannot live in a standard dark-subspace at finite n. It must live either (i) in the **inverse limit** (Phase 5 territory, Benoist-Bruneau-Pellegrini 2024 + 2509.13377), or (ii) in an **approximately-dark structural subspace** (Phase 2: R76 / T_diag (1,4) / W_{k-1}) that becomes exactly dark only in the inverse limit.

## Decision

Per attack plan §"Phase 1 → Decision tree":
- dim(A') = 1 at finite n → **route to Phase 2 (approximate-darkness of structural subspaces)**.
- Skip Phase 3 (explicit decomposition unneeded — there are no non-trivial central projections at finite n).

## Phase 2 specification (next, 1-2 sessions)

For each candidate structural-invariant subspace D ⊂ H_n, compute the "approximate-darkness" ratio
`α_D(M) := ‖P_D · M · P_{D^⊥}‖ / ‖M‖`
for each Kraus operator M in the family. Small α = approximately dark (Kraus operators almost preserve D); large α = strongly mixed.

Candidates (each carries known c=7/45-related spectral content):
- **D_R76**: R76 conservation subspace, `{M_n : Σ_j M_n(η_0 + j·3^k) = 0 ∀ η_0, k}` — verified exact for n ≥ 2 (R76 §11). The "slow mode" at the moment level.
- **D_T_diag**: T_lead's (1, 4) eigenspace lifted from class-resolved (P_+, P_-) to full V_n^M. Carries the 43/45 eigenvalue.
- **D_{W_{k-1}}**: 3-fiber-zero-mean subspace at level k. Contains K_k's entire kernel (this session's K_STRUCTURE_RESULT).

If any candidate has uniformly small α across the Kraus family, it becomes the "candidate dark subspace" for the inverse limit (where exact darkness recovers per Benoist-Bruneau-Pellegrini 2024 infinite-dim machinery).

## Effort

Phase 1: ~30 min total (probe was pre-written; ran in <1s, wrote result in this file). On schedule per the 10-15 session estimate.

## Caveat

The fully-irreducible-at-finite-n result is **NOT a closure of c=7/45**. It is a confirmation that the closure must come via inverse-limit machinery, not finite-n discrete-eigenvalue extraction. Phase 2-5 do the substantive work; Phase 1 just rules out the "easy path."

The interpretation is structurally clean: Syracuse mixes fully at each finite level, but the c=7/45 spectral feature is an **asymptotic** property accessible only at the inverse limit.
