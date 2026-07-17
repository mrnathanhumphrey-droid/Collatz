# R76_S11_DISPOSITION — top-level summary of the R76 §11 (P_+, P_−) 2D probe

**Date:** 2026-05-12. Fourth spectral probe in the M_3 sequence (T_3 → R_k → Candidate A → R76 §11 2D).
Wilson (analyst) reporting to Nathan.

---

## DISPOSITION: **INCONCLUSIVE**

> **R76 §11 articulates the 2D (P_+, P_−) class-resolved structure rigorously, articulates the (1, 4)-eigendirection rigorously, but articulates the eigenvalue 1/2 claim on that direction only as a structural conjecture — explicitly labeled "Open" in §11's closing paragraph. The rigorously-derived version of the 2D operator (T_diag = (1/5)·[[1, 1], [4, 4]]) has spectrum {0, 1} over Q, with eigenvalue 1 on (1, 4) — NOT eigenvalue 1/2.**
>
> **The operator that would actually carry eigenvalue 1/2 on (1, 4) is T_diag + Off_lin, where Off_lin is the linearization of the off-diagonal bilinear correction described in R77 sketch §5. That linearization is open work — R77 sketch §10 schedules it as "1–2 hours of focused implementation" but it has not been executed.**
>
> **Phase 2 construction of (P_+, P_−) cannot verify R76 §11's eigenvalue-1/2 claim because the claim is conjectural in §11 itself, and the rigorous component of the operator demonstrably gives a different eigenvalue. Phase 2 would have to RECONSTRUCT the (claimed) full operator from scratch — exactly the reconstruction-instead-of-verification pattern that the parent task warns against (the R_K probe's intractability source).**
>
> **The probe terminates at Phase 1 honestly. Phases 2–5 deliverables are not produced. The pre-registration's INCONCLUSIVE disposition fits this outcome: "R76 §11 articulation is too implicit to construct (P_+, P_−) explicitly, even though the section exists and gestures at the structure." Modifying "too implicit" to "articulated as conjecture, rigorous version delivers a different eigenvalue" — the substance is the same.**

---

## Why INCONCLUSIVE rather than the alternatives

Five other dispositions were pre-registered:

- **H_2D_CARRIES_RATE_AND_CLOSES** / **H_2D_CARRIES_RATE_NEEDS_LARGER_K** / **H_2D_CARRIES_RATE_BORDERLINE**: All three require an explicit 2D operator with eigenvalue 1/2 in hand. R76 §11 does not provide one. The rigorously-derived part of §11's operator (T_diag) has eigenvalue 1 on (1, 4), not 1/2. The off-diagonal extension that would deliver eigenvalue 1/2 is open. None of these three dispositions can be honestly chosen.

- **H_2D_DOESNT_CARRY_RATE_AT_ALL**: This would require establishing that the 2D structure cannot carry rate-1/2 even in principle. The probe cannot establish this — the off-diagonal extension might well deliver eigenvalue 1/2 if rigorously constructed; the structural ε_n = 10δ_+ identity and the (1, 4) deviation direction are both rigorous and compatible with a 1/2-rate-carrying operator existing. The probe cannot rule it out.

- **H_R76_S11_DOESNT_HAVE_CLAIMED_CONTENT**: This would mean R76 §11 doesn't articulate the (1, 4) + 1/2 structure at all. But it does — verbatim quoting in R76_S11_VERIFICATION.md §(a) shows R76 §11 explicitly names "(1, 4)" and "1/2" together. The Candidate A agent did not invent the claim. The agent's overinterpretation was treating R76 §11's "structural conjecture" as established content. That's a degraded version of "content present", not "content absent".

**INCONCLUSIVE** is the disposition that says "Phase 1 verification yields: content present in conjectural form; operator with claimed spectrum not constructed in project; constructing it is open work that R76 §11 itself declares open." The probe cannot proceed past Phase 1 without doing the open work — which is a separate probe, not this one.

---

## Structural meaning of this disposition

The R76 §11 2D probe was the **first positive candidate** the project's own documentation named. The previous three probes (T_3, R_k, Candidate A) each carved the "where rate-1/2 lives" boundary more sharply by ruling out a candidate operator/decomposition. This fourth probe was supposed to be different in kind — verifying a candidate the project itself had identified, not introducing a new one.

The probe's actual finding is structurally similar to the previous three: **the project's identified candidate is articulated as conjecture, and the rigorous operator-theoretic content does not deliver the claimed spectrum.** The "carving the boundary" pattern continues, with a slightly different shape:

- **T_3 (R77.3):** order-3 companion with conjectured spectrum {1/2, 1/4, 1/8}. **Falsified** — the spectrum doesn't fit.
- **R_k (R77.5):** inter-level residual operator with conjectured eigenvalue 1/2. **Intractable** — the framework is vector-valued, not operator-valued; M_3' uncharacterizable.
- **Candidate A / W_k via φ_n (this week):** L²(Ẑ_3^×) multiresolution decomposition with hypothesized cross-k structure carrying rate-1/2. **Falsified** — φ_n lives entirely in W_{n−1}, single-level dominance F2.
- **R76 §11 2D (this probe):** 2D (P_+, P_−) operator with conjectured eigenvalue 1/2 on (1, 4). **Inconclusive** — operator articulated rigorously in part (T_diag with spectrum {0, 1}), conjecturally in part (full T with spectrum claim 1/2 on (1, 4)). Rigorous part doesn't deliver 1/2; conjectural part is open work R76 §11 itself names open.

The pattern across all four probes: **every named candidate operator either fails its closure inequality, fails to exist as a well-defined operator, or has its rate-1/2-carrying spectrum left as open conjecture by the project's own documentation.** Four probes, four different shapes of "the operator does not exist in the project today".

This is meaningful information about the project's structural shape, not just a sequence of dead ends. It suggests one of:
1. **The rate-1/2 phenomenon is a branch-cut/density feature** (R77.6 framing), not an eigenvalue, and any spectral probe will keep returning "no spectral carrier". The W_k/Candidate A finding hints at this — F2 shows φ_n concentrates at the moving finest scale, consistent with a density-flow rather than a fixed-spectral-direction picture.
2. **The rate-1/2 carrier exists but requires constructing an operator the project has not yet built.** R77 sketch §5/§7 outlines the construction (off-diagonal bilinear linearization on (P_+, P_−)). It's "1–2 hours of focused implementation" in §10's estimate. Doing it would either verify R76 §11's conjecture (in which case Phase 2–5 of this probe become well-defined) or falsify it (in which case the carving continues).
3. **A different operator-theoretic framework is needed** (Kozyrev wavelets via Candidate B, or transfer operator on coherent extension via Candidate C, or something else entirely).

The probe doesn't decide among these three. It does establish that R76 §11 alone, as currently written, is not a verifiable spectral anchor — it's a well-articulated conjecture waiting for the open derivation §11 itself names.

---

## Phase 1 summary

Full content in R76_S11_VERIFICATION.md. Headlines:

1. **R76 §11 names the 2D (P_+, P_−) structure rigorously** — the structural collapse P^{+−}(c) = 0 + class-c-symmetry is proved (Theorem 76.x machinery from §11 + §10).
2. **R76 §11 names the (1, 4)-eigendirection rigorously as the asymptotic deviation direction** — empirical through k=6, structurally grounded in R64.B's class mass ratio (1/3, 2/3) squared = (1, 4).
3. **R76 §11 names the eigenvalue 1/2 on (1, 4) as a "structural conjecture"** — verbatim: "Open: rigorous derivation of leading coefficient 1/30 (numerical fit) and (1,4) eigenvalue = 1/2 (structural conjecture). Both reduce to algebraic identities from R66's chain dynamics."
4. **The rigorously-derived component T_diag = (1/5)·[[1, 1], [4, 4]] has char poly λ² − λ = 0 over Q, with spectrum {0, 1} and eigenvector (1, 4) at λ = 1.** The eigenvalue 1/2 claim is NOT delivered by T_diag.
5. **The full operator T = T_diag + Off_lin that R77 sketch §3 conjectures has eigenvalue 1/2 on (1, 4) is unconstructed.** R77 sketch §6 ledger places it on the open side. R77 sketch §10 schedules it as 1-2 hours of focused implementation but the implementation has not happened.

Phases 2–5 are not executed. The gating condition ("Phase 1 confirms R76 §11 has the content") is not cleanly met: content present as conjecture, not as established result.

---

## Routing recommendation

Three plausible next moves, surfaced for Nathan (not chosen by this probe):

### Route 1: Execute R77 sketch §7's implementation outline

The "1–2 hours of focused implementation" R77 sketch §10 names. Concretely:
1. Build T_N = T_diag + Off_lin as a 2x2 matrix over Q, where Off_lin is the linearization of the off-diagonal bilinear correction. This requires writing out the bilinear sum Σ_{v ≠ v'} 2^{−v−v'} A_v(ξ) A_{v'}*(ξ) ... over the (P_+, P_−) basis explicitly.
2. Compute spectrum of T_N over Q. Verify eigenvalue 1/2 on (1, 4) or refute it.
3. If verified: re-run Phases 2–5 of this probe with T_N as the operator. Compute M_3'' = ‖(I − T_N)^{−1}‖ and evaluate the Nisoli closure inequality.

If Route 1 succeeds, the rate-1/2 phenomenon has its first rigorous operator-theoretic anchor. The original R76 §11 + R77 sketch trajectory closes.

If Route 1 fails (constructed T_N's spectrum doesn't include 1/2 on (1, 4)), the carving continues — R76 §11's conjecture is refuted, and the project needs a different operator.

This is a probe-shaped task, well-scoped, and natural to attempt as the immediate next move. Recommended.

### Route 2: Route to Candidate B (Kozyrev wavelets) per CANDIDATE_A_DISPOSITION.md's primary recommendation

Sidesteps R76 §11 entirely. The Kozyrev framework provides a different orthonormal basis on L²(Ẑ_3^×) where translation and dilation act jointly, potentially localizing rate-1/2 to a specific wavelet index.

The CANDIDATE_A_DISPOSITION.md honest caveat applies: if Kozyrev wavelets also concentrate φ_n at the finest scale (as the W_k filtration did), the F2-shaped finding generalizes and the rate-1/2 is not a wave-localization feature.

This is the pre-registered next probe from CANDIDATE_A_DISPOSITION.md. Lower scoping cost than Route 1's "build T_N explicitly" (Kozyrev basis is well-documented in the literature), but higher than what this probe was supposed to be (verification, not construction).

### Route 3: Recognize the structural pattern and pivot to non-spectral framing

Four probes, no positive spectral carrier. R77.6's branch-cut framing (rate-1/2 at z = 2 as a spectral-density feature, not an eigenvalue) is consistent with all four negative outcomes. Pivoting to a generating-function / density-of-states / analytic-singularity framework — rather than searching for a fixed spectral subspace carrying rate-1/2 — may be the right structural move.

This is a research-direction pivot, not a probe. Costs more upfront but might resolve the "where does rate-1/2 live" question by reframing rather than continuing the search.

---

## Adversarial check outcomes (probe-level)

**(A1) R76 §11 fidelity** — verbatim quoting in R76_S11_VERIFICATION.md §(a); the "Open" / "structural conjecture" labels are R76 §11's own. ✓

**(A2) (1, 4)-eigendirection basis** — basis is (P_+ − 7/150, P_− − 14/75), articulated explicitly in §11. ✓

**(A3) Eigenvalue 1/2 over Q vs approximately** — T_diag's eigenvalue 1 on (1, 4) is exact over Q (char poly λ² − λ = 0). The 1/2 eigenvalue is numerical only (`T_lead_2x2.py` empirical fits and off-diagonal ratio ≈ 0.503 at k=5→6). No rigorous Fraction(1, 2) result on any project operator over Q. ✓

**(A4) Operator ↔ moment functional relationship** — the link ε_n = 10 δ_+ via Plancherel is rigorous. (1, 4)-mode = ε_n-carrying mode is rigorous. ✓

**(A5) Conflict with R77.3 / R77.4 / R77.5** — no inter-result conflict; the 2D operator is distinct from K_n / T_3 / R_k. The intra-R76 §11 gap (rigorous T_diag has eigenvalue 1, conjectural T has eigenvalue 1/2) is documented. ✓

---

## Deliverables produced

- **R76_S11_VERIFICATION.md** — Phase 1 content verification, verbatim R76 §11 quoting + diagnosis + decision rationale
- **R76_S11_DISPOSITION.md** (this file) — top-level disposition + structural meaning + routing

Phases 2 (P_PLUS_MINUS_CONSTRUCTION.md), 3 (M_3_DOUBLEPRIME.md), 4 (P_CLOSURE_TABLE.md) deliverables are **not produced**. The pre-registered gating condition ("Phase 1 confirms R76 §11 has the content") was not cleanly met; the parent task's pre-registration explicitly says "Don't proceed past it if the content is ambiguous." The content is present as conjecture, which falls between "confirmed" and "absent". Honest non-production is the appropriate response.

---

## Files referenced

- `result_76_conservation_law.md` §11 — primary source for R76 §11 content
- `result_77_sketch.md` — companion sketch with off-diagonal scheme and §10 schedule
- `result_77_T_lead_spectrum.md` — rigorous T_diag derivation + ledger of rigorous vs empirical
- `result_77_T_diagonal.py` — T_diag char poly computation
- `T_lead_2x2.py` — numerical 2x2 fitting (off-diagonal ratio empirical → 0.503)
- `result_77_4_K_spectrum_erratum.md` — prior probe's K_k spectrum falsification (no conflict)
- `result_77_2_nisoli_certification.md` — Nisoli closure framework (the inequality we did not get to evaluate)
- `M3_CLOSURE_TABLE.md` — parameterized closure framework reused if Phase 4 were to run
- `CANDIDATE_A_DISPOSITION.md` — prior probe + the "active anchor" remark that prompted this probe
- `CANDIDATE_A_PATTERN_MATCH.md` — prior probe Phase 4
