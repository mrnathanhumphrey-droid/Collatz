# CROSS_FREQ_DISPOSITION — top-level disposition

**Date:** 2026-05-12. Sixth spectral probe (T_3 → R_k → Candidate A → R76 §11 2D → T_N construction → **cross-frequency closure derivation**). Wilson (analyst) reporting to Nathan. Fork 1 from T_N_DISPOSITION's H_OFF_LIN_UNDERSPECIFIED.

---

## DISPOSITION: **H_CROSS_CLOSES_ON_ENLARGED_SPAN**

> **The cross-frequency bilinears Q_n^{ab}(c; v, v') for v ≠ v' do reduce to a closed family of moments — but the family is NOT span{P_n^{ab}(c)}. The closure space is V_M := span{M_n^{ab}(g, c) : g ∈ {0, 2, 4, 6, ...}, (a,b) ∈ {+,-}², c ∈ {1, 2}}, an enlarged moment basis parameterized by g (the difference of geometric-step indices v' − v). The g = 0 slice IS span{P_n^{ab}(c)}; the g ≥ 2 slices are genuinely new dimensions, each carrying a sublattice-constrained bilinear of π_n.**
>
> **R77 sketch §5's assertion "Substituting Tao's recursion expresses [P_{n+1}^{ab}(c)] as quadratic forms in {P_n^{ab}(c)}" is FALSE AS STATED for v ≠ v'. The correct statement is: as quadratic forms in {M_n^{ab}(g, c)} for g in an enlarged index set. The (P_+, P_−) 2x2 picture of R76 §11 + T_lead_2x2.py is a (1, 4)-PROJECTION of the full operator on V_M, not a faithful 2x2 over Q.**
>
> **The rate-1/2 eigenvalue (empirical) corresponds to a leading mode of the full operator T_V on V_M, projected onto (P_+, P_−). Spectrum computation over Q requires constructing T_V at a truncated g_max — this is Phase 3b/3c follow-up work of similar scale to the original T_diag derivation, NOT completed in this session.**

---

## Pre-registered hypotheses, decided

| Hypothesis | Status |
|---|---|
| H_CROSS_CLOSES_ON_SAME_FREQ (closure on {P_n^{ab}(c)}, 2x2 over Q is well-defined) | **REJECTED** — sublattice constraint r' ≡ 2^g·(r + ẽ_g) mod 3^{n-1} is genuinely different from r = r' (the P-diagonal) for g ≥ 2. |
| **H_CROSS_CLOSES_ON_ENLARGED_SPAN** (closure on {P_n^{ab}(c)} ∪ {M_n^{ab}(g ≥ 2, c)}, finite or countably-infinite) | **CHOSEN** — Phase 1 derivation establishes this is the correct closure space. |
| H_CROSS_PARTIAL_CLOSURE (some (v, v') close, others don't) | NO — all (v, v') with v ≠ v' reduce identically to a g-parameterized family. The same-parity pairs (v_3(d) ≥ 1) survive; mixed-parity pairs (v_3(d) = 0) vanish via lift-fiber orthogonality. This is uniform across (v, v'). |
| H_CROSS_DOESNT_CLOSE (no finite-dim closure on any span) | NO — V_M is a well-defined closure; the question is whether the cascade of Tao iteration on V_M is finite-dim (Phase 3b/3c). |
| INCONCLUSIVE | NO — Phases 1, 2 derive the structural answer; Phase 3a documents the cascade structure. Phase 3b/3c is identified but recognized as a separate multi-session probe. |

---

## Trajectory placement

| Probe | Object | Disposition | Why |
|---|---|---|---|
| T_3 (R77.3) | 3x3 companion matrix | FALSIFIED | Spectrum {1/2, 1/4, 1/8} doesn't describe ε_n |
| R_k (R77.4 erratum §1) | Inter-level residual operator | H_R_K_INTRACTABLE | c_k = 0 structural; doesn't transport |
| Candidate A | W_k φ_n bilinear-pair-form | H_CANDIDATE_A_FALSIFIES_F2 | φ_n in W_{n−1} only |
| R76 §11 2D | T_diag + Off conjectural | INCONCLUSIVE | Off unconstructed |
| T_N construction | T_diag + Off_lin as 2x2 over Q | H_OFF_LIN_UNDERSPECIFIED | §5 procedure doesn't specify Off_lin |
| **Cross-frequency closure (this)** | Closure space for Off_lin | **H_CROSS_CLOSES_ON_ENLARGED_SPAN** | **Closure exists, but on V_M, not on {P_n^{ab}(c)}.** Spectrum computation is Phase 3b/3c. |

The six-probe pattern now reads: **every named candidate operator on the same-frequency moment basis {P_n^{ab}(c)} fails its closure inequality or doesn't exist. The actual closure of Tao's bilinear recursion lives on an enlarged moment basis V_M parameterized by g ≥ 0 (with g = 0 recovering the same-frequency case).** This is a positive structural finding, even though no Q-spectrum has yet been extracted.

**Comparison to prior probes:** T_N construction's H_OFF_LIN_UNDERSPECIFIED was a "the work to derive Off_lin is not in the project" finding. This probe DID the work — and found the closure exists on a larger space than R77 sketch §5 claimed. This is meaningfully different from the prior probes:

- It's not "the operator doesn't exist" (cf. K_n, R_k single-level operator failures).
- It's not "the operator exists but has wrong spectrum" (cf. T_diag's eigenvalue 1 on (1, 4)).
- It is: "the operator exists on a NATURAL ENLARGED SPACE; the rate-1/2 question is now about that space's spectrum."

Structurally: this probe converts the rate-1/2 question from "find an operator" to "find the spectrum of an operator we now have, on V_M." That is a real advance even without the spectrum being computed.

---

## What's new vs. T_N_DISPOSITION.md (the prior H_OFF_LIN_UNDERSPECIFIED probe)

T_N_DISPOSITION.md's three routing options:
1. **Route 1: Derive the cross-frequency closure (substantive analytical work).** ← THIS PROBE EXECUTED ROUTE 1, partially.
2. Route 2: Pivot to R77.6 branch-cut / density-of-states framing.
3. Route 3: Recognize the pattern, document structural boundary as no-go.

This probe found Route 1 IS tractable — the cross-frequency closure exists, on V_M. But:
- The closure space is LARGER than R77 sketch §5 claimed (not {P_n^{ab}(c)}; rather {M_n^{ab}(g, c)} for g ∈ {0, 2, 4, ...}).
- Computing the spectrum on V_M (truncated at g_max) is the natural Phase 3b/3c follow-up; not done here.
- The (P_+, P_−) 2x2 reduction is correct only as a (1, 4)-projection, not as a faithful Q-operator.

The structural meaning is: **the rate-1/2 phenomenon, if spectral, lives on V_M as a leading eigenvalue of T_V; its (1, 4)-projection onto (P_+, P_−) is what T_lead_2x2.py empirically sees.**

This is the most positive routing outcome consistent with the six-probe trajectory: the spectral framework IS viable, on a slightly larger dimension than originally hypothesized.

---

## Key algebraic findings (rigorous, derived in CROSS_FREQ_PHASE1_EXPANSION.md)

1. **Lift-fiber orthogonality (3-adic):** For v ≠ v', let d_{v,v'} := 2^{-v} - 2^{-v'}. The sum over the 3 lifts of u ∈ (Z/3^n)^× into (Z/3^{n+1})^× of e^{-2πi j·d/3} is:
   - 3 if 3 | d (v_3(d) ≥ 1)
   - 0 if 3 ∤ d (v_3(d) = 0)
   
   This kills mixed-parity (v even, v' odd) contributions to P^{++}, P^{−−} automatically. P^{+−} also receives zero from cross-parity AND from impossible same-parity. Hence **P^{+−} = 0 is preserved by Tao's recursion structurally** — recovering R76 §11's empirical observation as algebraic identity.

2. **3-adic valuation of d_{v,v'} for v, v' same parity:** For g := v' - v ≥ 2 even, v_3(d_{v,v'}) = 1 + v_3(g/2) ≥ 1. Hence same-parity pairs always survive the lift-fiber orthogonality.

3. **g-reduction:** The cross-frequency object Q_n^{ab}(c; v, v') depends only on g = v' - v (not on v, v' separately), after the lift-fiber + unit-shuffle reduction. This collapses the 2D (v, v') sum to a 1D g sum.

4. **Sublattice constraint:** M_n^{ab}(g, c) reduces to a partial sum over (r, r') pairs with r' ≡ 2^g·(r + ẽ_g) - 2^g·3^{n-1}·m mod 3^n, for some m ∈ {0, 1, 2}, weighted by ω^{-c·m}. The constraint sublattice depends nontrivially on g.

5. **(1, 4)-direction preservation:** The total off-diagonal contribution to (Off_{n+1}^{++}, Off_{n+1}^{−−}) is along (1, 4) for ALL g (because W_−(g) / W_+(g) = 4 for all g). So Off projected to (P_+, P_−) is rank-1 in the (1, 4) direction.

6. **Weights of the g-sum:** W_+(g) = 2^{-g+1}/15 (for P^{++}). Each step g → g + 2 attenuates by 1/4. So the cross-frequency expansion is convergent in g; the leading g = 2 term dominates with weight 1/30, g = 4 contributes 1/120, etc.

---

## Open work (Phase 3b/3c) for full Q-spectrum

For the rate-1/2 question to be settled rigorously over Q, the following must be done:

(B1) **Derive the recursion of M_n^{ab}(g, c) → M_{n+1}^{ab}(g', c') under Tao's iteration.** Each μ̂_{n+1} expansion at the M-moment definition produces another cross-frequency object. By the same Phase 1 §3–§5 analysis, this reduces to a g'-parameterized family — but with g' potentially expanding the index set (cascade question).

(B2) **Identify whether the cascade closes at finite g_max.** If yes: V_M^{(g_max)} is a finite-dim space, T_V is a finite-rank operator over Q. If no: V_M is infinite-dimensional, spectral analysis requires Banach-space methods.

(B3) **At finite truncation g_max ∈ {2, 4}, build T_V matrix over Q and compute spectrum.** Check for eigenvalue 1/2; identify its eigenvector and project to (P_+, P_−); verify (1, 4)-image.

(B4) **Nisoli closure on V_M.** If T_V has spectrum {1, 1/2, ...} stably across truncations and the truncation errors decay, apply Nisoli Theorem 2.15 to certify rate 1/2 rigorously.

Effort estimate: each of B1-B4 is comparable to one of the prior R77.x probes. Total: multi-session, perhaps 3-5 sessions of focused analytical work.

---

## Routing recommendation (surfaced for Nathan, not chosen)

### Route A: Execute Phase 3b/3c

Derive M(g) recursion (B1), check cascade closure (B2), compute spectrum (B3), apply Nisoli (B4). If everything works, the rate-1/2 question is closed rigorously over Q for the first time. If the cascade doesn't close at finite g_max, this route requires Banach-space framework and pivots toward Route C.

Recommended if a positive spectral closure of c = 7/45's rate is the priority and Nathan is willing to invest 3-5 sessions of focused work.

### Route B: Pivot to R77.6 branch-cut framing

Even if the cascade doesn't close finite-dim, the GENERATING FUNCTION G(z) := Σ S_n z^n's analytic structure encodes the rate-1/2 directly as a singularity at z = 2. The branch-cut framework (Result 77.6 candidates) sidesteps operator existence and reads rate-1/2 off the analytic continuation of G(z).

This is a parallel-track approach that doesn't require V_M's spectral structure. It's also the user's `project_collatz_r78_bilinear_cracked` memory-noted "Move 3" alternative (R77.5 §7+§9 projective-limit framework).

Recommended if the V_M spectral analysis is anticipated to be hard (closed-cascade unclear) and a parallel approach is worth pursuing.

### Route C: Document the positive structural finding now

The current probe's findings — **cross-frequency closure exists on V_M; (P_+, P_−) is a (1, 4)-projection of T_V; structural collapse P^{+−} = 0 is algebraic, not empirical** — are themselves publishable as a structural result.

Combined with the bilinear bound side (user's memory `project_collatz_r78_bilinear_cracked`, 25-commit burgess.zip), the publishable output expands from "bilinear bound + no-go for spectral side" to "bilinear bound + spectral side has positive structural anchor in V_M (cascade closure open, rate-1/2 spectral question reformulated)."

This is the most honest immediate framing without committing to further multi-session work.

---

## Adversarial check outcomes (probe-level)

**(A1) R77 sketch §5 fidelity.** All derivations trace to §5's verbatim text plus Tao's recursion (c_seven_forty_fifth.md §3) plus R66 class flow plus standard 3-adic number theory (LTE_3). No invoked identity is outside this scope. ✓

**(A2) Phase orthogonality vs cancellation.** Phase 1 §3 distinguishes v_3(d) = 0 (kills mixed-parity contributions) from v_3(d) ≥ 1 (passes to level-n character). Both outcomes documented. The "phase is non-trivial → reduces to {P}" alternative was rejected: the level-n character sum does NOT reduce to {P_n^{ab}(c)}; it reduces to sublattice-constrained moments. ✓

**(A3) (P_+, P_−) basis convention.** R76 §11 defines P_+ := P^{++}(c=1) = P^{++}(c=2) etc. for n ≥ 2. Off_lin's effective 2x2 on this basis is the rank-1 (1, 4)-image of the full T_V on V_M. Q-arithmetic on the 2x2 is not faithful; the (M(g ≥ 2)) side information is required. ✓ — this matches T_N_OFF_LIN_SPEC.md Obstruction 2.

**(A4) Exact rationals throughout.** All derived quantities (ẽ_g, weights W_±(g), sublattice equations) are in Q (or 3-adic integers reducible to Z/3^n exactly). No floating-point dependence in the structural claim. The cross_freq_compute.py script demonstrates the structural claim numerically; the algebraic claim does not depend on the numerical check. ✓

**(A5) Spectrum-vs-no-spectrum dichotomy.** H_CROSS_CLOSES_ON_ENLARGED_SPAN means the closure exists on V_M. Whether T_V on V_M has eigenvalue 1/2 over Q (computable at finite truncation) is a Phase 3b/3c question. Possible outcomes:
- Eigenvalue 1/2 survives at all truncations → rate-1/2 spectral closure achieved.
- Eigenvalue near-1/2 but doesn't stabilize at exact 1/2 → cascade may need infinite truncation; pivot to Route B.
- Different leading eigenvalue → reframes c = 7/45's rate.

Each of these is possible; Phase 3b/3c would decide. ✓

**(A6) Conflict with R77.4.** R77.4 ruled out K_n (residue-probability Markov operator) as carrying eigenvalue 1/2. T_V on V_M is a fundamentally different object (bilinear moment operator on class-resolved + g-twisted moments). The two operators act on different spaces; no conflict. ✓ Note: T_V might have eigenvalue 1/2 even when K_n doesn't, since the operators capture different aspects of π_n's dynamics.

---

## Deliverables produced (in C:/Collatz/)

- **CROSS_FREQ_PHASE1_EXPANSION.md** — Phase 1 leading-pair derivation (rigorous, traces to R77 §5 + Tao recursion + LTE_3)
- **CROSS_FREQ_PHASE1_SPAN.md** — Phase 2 span identification (sublattice analysis, V_M definition)
- **CROSS_FREQ_HIGHER_PAIRS.md** — Phase 3a cascade structure sketch (without Phase 3b/3c)
- **CROSS_FREQ_DISPOSITION.md** (this file) — top-level disposition + structural meaning + routing
- **cross_freq_compute.py** — reproducible verification script (n = 2, 3; X_n^{ab}(c; g) for g ∈ {2, 4, 6}; computational rank check)

Phases not produced:
- OFF_LIN_MATRIX.md — only applicable under H_CROSS_CLOSES_ON_SAME_FREQ; not chosen
- T_N_SPECTRUM_FROM_DERIVATION.md — same, not applicable
- (Phase 3b/3c deliverables) — multi-session follow-up, not in scope

---

## Files referenced

- result_77_sketch.md §5 — primary source for the asserted closure
- result_77_T_lead_spectrum.md §1, §2, §3, §6 — T_diag rigorous; off-diagonal heuristic + ledger
- result_77_T_diagonal.py — rigorous T_diag derivation (v = v' diagonal collapse)
- T_lead_2x2.py — numerical 2x2 fit
- result_76_conservation_law.md §11 — structural collapse and (1, 4)-eigendirection
- c_seven_forty_fifth.md §3 — Tao's recursion verbatim
- r66_per_a_decay_corrected.md — R66 class flow rule and per-class decay
- R76_S11_VERIFICATION.md / R76_S11_DISPOSITION.md — prior probe context
- T_N_OFF_LIN_SPEC.md / T_N_DISPOSITION.md — prior probe identifying the open closure
- result_77_5_inter_level_residual.md — multi-resolution decomposition context

---

## Synopsis (one paragraph)

The cross-frequency closure of Tao's bilinear recursion exists. It lives on the enlarged moment basis V_M = span{M_n^{ab}(g, c) : g ∈ {0, 2, 4, ...}}, NOT on the same-frequency basis {P_n^{ab}(c)} as R77 sketch §5 asserted. The derivation reduces every (v ≠ v') pair to a single g-parameter (= v' − v) family, with sublattice-constrained bilinears of π_n at each g. The (P_+, P_−) 2x2 picture is the (1, 4)-projection of the full operator on V_M, structurally consistent but not Q-faithful as 2x2. The rate-1/2 question reformulates as: does T_V on V_M (truncated at some g_max) have eigenvalue 1/2 over Q? This is Phase 3b/3c work — a multi-session follow-up, comparable in scale to the original T_diag derivation, identified clearly but not executed in this session.
