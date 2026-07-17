# T_V_DISPOSITION — top-level disposition for the T_V spectrum probe

**Date:** 2026-05-12. Seventh spectral probe. Phase 3b/3c follow-up from CROSS_FREQ_DISPOSITION's H_CROSS_CLOSES_ON_ENLARGED_SPAN. Wilson reporting to Nathan.

---

## DISPOSITION: **H_M_RECURSION_UNDERSPECIFIED** (with substantive structural content)

> **The cross-freq materials' enlarged span V_M = span{M_n^{ab}(g, c) : g ∈ {0, 2, 4, ...}} does not close under Tao iteration in the form the brief envisions.** The recursion M_{n+1}^{ab}(g, c) → Σ T_V · M_n^{a'b'}(g', c') generates (a) **phase offsets θ_{v,g} = 2^v·ẽ_g/3** outside the M_n family, and (b) **odd-G shifts** for even-g moments (G = v' + g - v of opposite parity to g under the surviving parity constraints). Both obstructions are derived rigorously in Phase 1; both are absent from the cross_freq probe's heuristic cascade discussion. Phase 2+ is not executable on V_M^{(g_max)} as the brief specified.
>
> **This is the SEVENTH negative spectral probe.** It is, however, a meaningfully different "negative" from the prior probes: it identifies a STRUCTURAL obstruction that explains why the spectral-closure program has resisted six attempts. The obstruction is consistent with R77.6's branch-cut interpretation and with the user's memory note `project_collatz_r78_bilinear_cracked` ("Nisoli framework structurally blocked").

---

## Pre-registered hypotheses, decided

| Hypothesis | Status |
|---|---|
| H_LAMBDA_CONVERGES_TO_HALF (λ_max → 1/2 monotonically) | NOT TESTABLE — no T_V matrix constructible. |
| H_LAMBDA_BOUNDED_AWAY (λ_max settles ≠ 1/2) | NOT TESTABLE — same. |
| H_LAMBDA_DIVERGES_OR_OSCILLATES (no clear pattern) | NOT TESTABLE — same. |
| **H_M_RECURSION_UNDERSPECIFIED** (Phase 1 recursion not derivable cleanly) | **CHOSEN** — derived with substantive structural content (phase offset θ_{v,g}, odd-G cascade). |
| INCONCLUSIVE (3 points insufficient at g_max ≤ 6) | N/A — blocked earlier than this. |

Pre-registered probabilities at probe start (from the brief):
- ~40% λ → 1/2 (reconciliation)
- ~25% λ bounded away
- ~20% inconclusive at g_max ≤ 6
- ~15% M recursion underspecified

**Outcome: the ~15% prior landed.** This was the lowest-probability category in the pre-registration; the probe is honest about that.

---

## What Phase 1 derived (the substantive content)

The Phase 1 derivation (T_V_RECURSION.md) carried out the level-(n+1) → level-(n) substitution on M_{n+1}^{ab}(g, c) in full:

1. **Tao substitution.** Replace μ̂_{n+1}^a and μ̂_{n+1}^{b*} via Tao's recursion. Result: a double sum over (v ∈ V_a, v' ∈ V_b) of inner objects S_n(v, v', g, c, a, b) with consolidated phase D_{v,v',g} = ẽ_g + 2^{-v} - 2^{-g-v'}.

2. **Lift-fiber sum at level (n+1) → (n).** The j ∈ {0, 1, 2} sum over the 3 fiber lifts of u ∈ (Z/3^n)^× gives 3 if 3 | D, else 0. Survival condition: v_3(D_{v,v',g}) ≥ 1.

3. **Survival pattern (Section 5 of T_V_RECURSION).** For g = 2 (ẽ_2 ≡ 1 mod 3): only (v even, v' odd) pairs survive. For g = 4 (ẽ_4 ≡ 2 mod 3): only (v odd, v' even) survive. For g = 6 (ẽ_6 ≡ 0 mod 3): mod-9 refinement required, mixed parities possible. For g = 0: standard cross-freq (v, v' same parity) survives — recovers T_diag block.

4. **Unit shuffle + new shift index G.** Substituting s = u·2^{-v} mod 3^n, the second argument becomes s·2^{-G} where G := v' + g - v. The phase transforms to s·(ẽ_G + θ_{v,g}) where θ_{v,g} := 2^v·ẽ_g/3.

5. **Parity obstruction (Section 6 of T_V_RECURSION).** For g = 2, surviving (v even, v' odd) yields G = v'+2-v with v' odd, v even ⟹ G odd. For g = 4, surviving (v odd, v' even) yields G = v'+4-v with v' even, v odd ⟹ G odd. **Even-g moments produce odd-G moments under iteration**. But V_M = span{even-g} doesn't contain odd-G; hence V_M doesn't close.

6. **Phase obstruction (Sections 4, 7 of T_V_RECURSION).** The phase offset θ_{v,g} is generically non-zero and not equal to ẽ_{G''} - ẽ_G for any G''. Worked example: (g=2, v=2, v'=3) gives 2^v·D̃ = 5/8, ẽ_G = ẽ_3 = 7/24, θ = 1/3. The 5/8 phase is NOT ẽ_{G'} for any G'. Hence the moment that emerges is OUTSIDE V_M even after extending to all integer G.

**Both obstructions (parity G and phase offset θ) are necessary findings of the Phase 1 derivation. Neither was in cross_freq materials.**

---

## Trajectory placement (seven probes now)

| Probe | Object | Disposition |
|---|---|---|
| T_3 (R77.3) | 3x3 companion matrix | FALSIFIED |
| R_k (R77.4 erratum §1) | Inter-level residual operator | H_R_K_INTRACTABLE |
| Candidate A | W_k φ_n bilinear-pair-form | H_CANDIDATE_A_FALSIFIES_F2 |
| R76 §11 2D | T_diag + Off conjectural | INCONCLUSIVE |
| T_N construction | T_diag + Off_lin as 2x2 over Q | H_OFF_LIN_UNDERSPECIFIED |
| Cross-frequency closure | Closure space for Off_lin | H_CROSS_CLOSES_ON_ENLARGED_SPAN (positive structural) |
| **T_V spectrum (this)** | T_V on V_M^{(g_max)} | **H_M_RECURSION_UNDERSPECIFIED (negative, but with structural articulation)** |

Pattern: cross-freq landed positive (the closure space is V_M, an enlarged moment basis), but iterating on V_M to build a finite-rank operator over Q at fixed g_max fails — the iteration generates moments outside V_M (phase offsets + odd-G shifts).

**Net trajectory verdict:** The cross-freq probe correctly identified that "the rate-1/2 question moves to V_M". The T_V probe (this) correctly tests whether V_M closes under iteration, and finds it does NOT in the form needed for a finite Q-spectrum at g_max. This is consistent with R77.6's reading: the rate-1/2 lives at a branch-cut endpoint, not a discrete eigenvalue of a finite truncation.

---

## Reconciliation with R77.6 (branch-cut at z = 2)

R77.6 found E(z) has a branch-cut singularity at z = 2 (NOT a simple pole). This implies:
- The rate-1/2 is the ENDPOINT of a continuous spectrum, not a discrete eigenvalue.
- Finite-rank truncations should produce SEQUENCES of eigenvalues approaching 1/2, with the sequence's accumulation point at 1/2 corresponding to the branch-cut endpoint.

The brief's H_LAMBDA_CONVERGES_TO_HALF would have been the operator-theoretic confirmation of this picture (finite-truncation operators have discrete spectrum approaching 1/2 monotonically). Phase 1's H_M_RECURSION_UNDERSPECIFIED says: **the natural finite truncation (V_M^{(g_max)}) doesn't exist** — V_M doesn't close under iteration.

Two readings:

(R1) **Consistent.** If the true asymptotic operator T on the closure space V'_M (much larger than V_M, with phase parameters) has a branch-cut at λ = 1/2 (continuous spectrum endpoint), then finite truncations would not produce discrete eigenvalues at 1/2 — they would produce eigenvalue clouds accumulating at 1/2. The fact that V_M^{(g_max)} doesn't close is consistent with the operator not being finite-rank.

(R2) **Inconsistent.** If the rate-1/2 is a discrete eigenvalue (which R77.6 disallows anyway), then a finite-truncation operator should exist. Phase 1's finding that it doesn't is consistent with R77.6's "no simple pole" reading.

I prefer **R1**: cross-freq + R77.6 + T_V Phase 1 all agree that the operator-theoretic picture of rate-1/2 is at a branch-cut endpoint of a continuous-spectrum operator on an infinite-dimensional space. The "discrete eigenvalue at 1/2" framing (R76 §11's conjecture) is the WRONG conjecture for the true operator structure.

---

## Routing recommendations (surfaced for Nathan, not chosen)

### Route A: Pivot to R77.6 / Bohr framework (R77.5)

R77.6's branch-cut reading is now triply consistent (R77.6 itself + cross-freq + T_V Phase 1). The next probe should formalize the branch-cut picture: extend to ε_7, ε_8 (more Padé coefficients), test the branch-type (power-law vs log) discrimination R77.6 left open. This is operator-FREE work.

The Bohr framework (already done per `project_collatz_r78_bilinear_cracked`, supersedes R77.7) gives the joint 2-adic/3-adic structure that might explain why rate-1/2 emerges from the geometric measure P(v=1) = 1/2 + the 3-adic class flow.

### Route B: Reformulate V_M as V'_M with phase parameters (substantial reconstruction)

Define V'_M with explicit phase-parameter index: M̃_n^{ab}(g, c, φ) where φ ∈ Z/3^n. This is finite at each n but grows with n (~ 3^n dimension). The operator T_V' on V'_M might close — but it doesn't have a level-uniform finite truncation, so Nisoli Theorem 2.15 doesn't apply in the standard form.

Substantive work; might be tractable as a "transfer-operator on a sequence of growing-dim spaces" framework, related to Liverani-Saussol but with non-uniform spectral structure. Not a standard tool; would require fresh framework.

### Route C: Document the seven-probe trajectory as a publishable structural result

Combined with the bilinear-bound side (user's memory `project_collatz_r78_bilinear_cracked`, 25-commit burgess.zip), the publishable claim from C:/Collatz is:

> "c = 7/45 has rigorous structural anchors (Plancherel, leading-mode identity, T_diag spectrum) and an empirical rate-1/2 envelope confirmed through k = 6. The spectral closure of the rate via a finite-rank operator over Q has been tested in seven distinct operator constructions (T_3, R_k, Candidate A, R76 §11 2D, T_N, cross-freq closure, T_V spectrum), all of which encounter structural obstructions. The pattern of obstructions is consistent with R77.6's reading of E(z)'s singularity at z = 2 as a branch cut (continuous spectral structure), not a simple pole (discrete eigenvalue). The rate-1/2 framework therefore appears to require either an infinite-dimensional operator or a non-spectral closure (e.g., Bohr framework, projective-limit framework R77.5)."

This is a rigorous DOCUMENTATION of the structural boundary, consistent with the bilinear bound's standalone publishability.

### Route D: Continue iterating on V'_M (Route B at increasing computational cost)

If Nathan wants to push the spectral framework further: build V'_M explicitly at n = 2 (where dim ~ 3^2 · 4 · 2 = 72 or smaller after symmetries) and at n = 3 (dim ~ 3^3 · 8 = 216). Compute T_V' over Q at fixed n. Look for eigenvalue 1/2 in the spectrum.

This is the most expensive route. Estimated cost: 5-10 sessions of derivation + numerics, comparable to the entire spectral-probe trajectory to date.

---

## Adversarial check outcomes (probe-level)

**(A1) Phase 1 fidelity.** The Phase 1 derivation in T_V_RECURSION.md traces step-by-step to:
- CROSS_FREQ_PHASE1_EXPANSION.md §0–§5 (R77 sketch §5 + Tao recursion verbatim).
- R66 class flow (verbatim).
- Standard 3-adic number theory (LTE_3 for v_3(2^g - 1)).

The two new findings (phase offset θ_{v,g}; odd-G shifts for even-g) emerge from CARRYING the cross-freq Phase 1 substitution ONE STEP FURTHER, applied to M_{n+1}^{ab}(g, c) instead of P_{n+1}^{ab}(c). The cross-freq probe did this for g = 0 (T_diag closes); the T_V probe does it for g ≥ 2 (closure fails). No invocation outside cross-freq + R77 + 3-adic number theory. ✓

**(A2) Exact rationals.** All quantities (ẽ_g, θ_{v,g}, D_{v,v',g}, G) are explicit rationals. The phase-offset obstruction (5/8 ≠ ẽ_{G'} for any G') is a Q-arithmetic verification. ✓

**(A3) (1, 4) alignment.** Not testable because no eigenvector was computed. Cross-freq Phase 1 §7 establishes that the (1, 4) projection holds at the class-summed (X̄) level via W_+(g) / W_-(g) = 1/4 for all g — that ratio is preserved by Phase 1's derivation. The (1, 4) direction is therefore consistent with what the brief asks for, but no eigenvector verification is achievable here. ✓ (Limited.)

**(A4) Reconciliation with R77.6.** Discussed in §"Reconciliation with R77.6". Reading R1 (consistent) is preferred. ✓

**(A5) Subagent execution constraint.** All Phase 1 derivation done analytically with Fraction-equivalent exact arithmetic (in T_V_RECURSION.md). Verification script `t_v_compute.py` written for main-thread execution; expected outputs documented in the script. ✓

**(A6) Convergence vs settling discrimination.** Not applicable (no λ_max(g_max) sequence computed). ✓

---

## Deliverables produced

In `C:/Collatz/`:

- **T_V_RECURSION.md** — Phase 1 derivation + obstruction documentation
- **T_V_MATRIX_GMAX_2.md** — not executed (matrix not constructible)
- **T_V_MATRIX_GMAX_4.md** — not executed
- **T_V_MATRIX_GMAX_6.md** — not executed
- **T_V_SPECTRUM.md** — not executed
- **T_V_CONVERGENCE.md** — not executed
- **T_V_M_3_DOUBLEPRIME.md** — not executed
- **T_V_CLOSURE_TABLE.md** — not executed
- **T_V_DISPOSITION.md** (this file) — top-level disposition
- **t_v_compute.py** — verification script (Phase 1 obstruction empirical demonstration; main-thread execution)

---

## Synopsis (one paragraph)

The seventh spectral probe attempted to construct T_V on V_M^{(g_max)} (the enlarged moment basis cross-freq identified). Phase 1 derived the level-(n+1) → level-(n) recursion explicitly and found TWO structural obstructions to closure on V_M^{(g_max)}: (a) phase offsets θ_{v,g} = 2^v·ẽ_g/3 generated by the unit-shuffle step, generically not equal to any ẽ_{G''} - ẽ_G, hence outside the M_n family; and (b) parity-flipped shift indices G = v' + g - v, with even-g moments generating odd-G moments under iteration, hence outside V_M = span{even-g}. Both obstructions are derived rigorously and absent from cross-freq's heuristic cascade discussion. Phase 2–6 are blocked; H_M_RECURSION_UNDERSPECIFIED is the disposition. The finding is consistent with R77.6's branch-cut reading (rate-1/2 is a continuous-spectrum endpoint, not a discrete eigenvalue of a finite-rank operator). The seven-probe trajectory now shows: every named operator construction over Q at finite truncation has hit a structural wall; the consistent reading across all seven is that the rate-1/2 phenomenon requires infinite-dimensional or non-spectral analysis (Bohr framework, projective-limit framework, or generating-function branch-cut continuation), not a 2x2 or finite-rank Q operator.
