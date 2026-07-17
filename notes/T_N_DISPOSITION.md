# T_N_DISPOSITION — top-level disposition of the T_N = T_diag + Off_lin construction probe

**Date:** 2026-05-12. Fifth spectral probe in the M_3 sequence (T_3 → R_k → Candidate A → R76 §11 2D → **T_N construction**).
Wilson (analyst) reporting to Nathan. Fork 1 from R76_S11_DISPOSITION.md.

---

## DISPOSITION: **H_OFF_LIN_UNDERSPECIFIED**

> **R77 sketch §5 articulates Off_lin as a procedure (substitute Tao's recursion into the bilinear moment definition) plus a claim ("quadratic forms in {P_n^{ab}(c)}") — not as a 2x2 matrix over Q.**
>
> **Expanding the procedure produces, at the v ≠ v' level, cross-frequency bilinears Q_n^{ab}(c; v, v') = Σ_{ξ ≡ c, 3∤ξ} A_v(ξ) A_{v'}^*(ξ) μ̂_n(ξ·2^{−v}) μ̂_n^*(ξ·2^{−v'}). These are bilinears in μ̂_n at TWO DIFFERENT FREQUENCIES (ξ·2^{−v} ≠ ξ·2^{−v'}). They do NOT collapse onto span{P_n^{ab}(c)} (the same-frequency basis) under any closure R77 sketch §5–§9 identifies.**
>
> **The §5 claim that the result is "quadratic forms in {P_n^{ab}(c)}" is the open closure, not a derivation. The project's own ledger (T_lead_spectrum.md §6) lists "off-diagonal exact bilinear-sum analysis" as Open. T_lead_spectrum.md §3 attempts the closure heuristically and stops at two contradictory candidates ("λ_2 = 4·(1/4) = 1? — no, more careful analysis needed" and "the cleanest derivation: λ_2 = P(v=1) = 1/2"). T_lead_2x2.py provides a numerical fit from data only, not a Q derivation.**
>
> **Without the cross-frequency closure, Off_lin is not constructable as a 2x2 matrix over Q. The probe terminates at Phase 1 honestly. Phases 2–5 deliverables are not produced. This is the pre-registered H_OFF_LIN_UNDERSPECIFIED failure mode: "R77 sketch §5 describes Off_lin in terms that don't pin down a specific 2x2 matrix over Q. The probe blocks at Phase 1 because reconstructing Off_lin would mean inventing rather than verifying."**

---

## Why H_OFF_LIN_UNDERSPECIFIED and not the alternatives

The pre-registration laid out six dispositions:

- **H_T_N_CONFIRMS_RATE_HALF / H_T_N_DIFFERENT_EIGENVALUE / H_T_N_NO_1_4_EIGENVECTOR / H_T_N_SPECTRAL_BUT_M3_DOESNT_CLOSE**: All four require Off_lin to be a concrete 2x2 matrix over Q so its spectrum can be computed (over Q). Phase 1 finds that R77 sketch §5 does not supply this matrix. None of these four can be honestly chosen.

- **INCONCLUSIVE**: This would apply if some other phase-internal obstruction blocked progress. The actual obstruction is at Phase 1's gating condition — the spec doesn't pin down the matrix — which is precisely **H_OFF_LIN_UNDERSPECIFIED**, not generic INCONCLUSIVE.

- **H_OFF_LIN_UNDERSPECIFIED** is what fits. Phase 1 is gating; pre-registration said "Phase 1 is gating just like in the prior R76 §11 probe. Don't reconstruct what isn't specified; only execute what is." The §5 specification is incomplete in a load-bearing way; the probe terminates honestly at Phase 1.

---

## Structural meaning

The T_N probe was the project's own scheduled work item — R77 sketch §10's "1–2 hours of focused implementation." The R76 §11 probe terminated INCONCLUSIVE specifically because the operator that would carry eigenvalue 1/2 (T_diag + Off_lin) is unconstructed; this Fork 1 probe was to do the construction.

What this probe found: the "1–2 hours of focused implementation" R77 §10 schedules is the **IMPLEMENTATION** assuming the closure exists (i.e., assume Off_lin is a known 2x2 matrix → compute spectrum → apply Nisoli → write up). The estimate refers to the bookkeeping cost, not the DERIVATION cost. The derivation — closing cross-frequency bilinears onto same-frequency moments — is its own substantive analytical step, of the same scale as the original T_diag derivation in result_77_T_diagonal.py.

The five-probe trajectory now reads:

| Probe | Object | Disposition | Why |
|---|---|---|---|
| **T_3 (R77.3)** | 3x3 companion matrix | **FALSIFIED** | Spectrum {1/2, 1/4, 1/8} doesn't describe ε_n; the recursion fit is artifact |
| **R_k (R77.4 erratum §1)** | Inter-level residual operator | **H_R_K_INTRACTABLE** | Reading A out of scope; Reading B not a Nisoli object; c_k = 0 → doesn't transport |
| **Candidate A (W_k φ_n)** | Bilinear-pair-form moment in W_k filtration | **H_CANDIDATE_A_FALSIFIES_F2** | Exact-Q: φ_n lives entirely in W_{n−1}, no cross-level structure |
| **R76 §11 (2D P_+, P_−)** | T_diag + (conjectured) full T | **INCONCLUSIVE** | T_diag rigorous with spec {0, 1}; full T (with eigval 1/2) unconstructed |
| **T_N construction (this)** | T_N = T_diag + Off_lin | **H_OFF_LIN_UNDERSPECIFIED** | R77 sketch §5's procedure doesn't pin down Off_lin as a 2x2 over Q; cross-frequency closure is the open step |

The pattern across all five probes: **every named candidate operator either fails its closure inequality, fails to exist as a well-defined operator, or has its critical specification left open by the project's own documentation.** Five probes, five different shapes of "the operator does not exist over Q in the project today".

This is meaningful structural information about the project's framework, not five separate dead ends. It points strongly toward one of:

1. **The rate-1/2 phenomenon is a branch-cut / density-of-states feature** (R77.6 framing), not an eigenvalue of any spectral operator. The W_k/Candidate A finding (F2: single-finest-scale dominance) and the cross-frequency obstruction here (no closure for the off-diagonal sum) are both compatible with a density-flow picture rather than a fixed-eigenvalue picture.

2. **A different operator-theoretic framework is required** (Kozyrev wavelets via Candidate B, transfer operator on coherent extension via Candidate C, or operator-on-Riesz-extension à la R77.6 §7+§9 projective-limit framework that the user's `project_collatz_r78_bilinear_cracked` memory mentions as a candidate Move 3).

3. **The closure-as-derivation work — i.e., closing the cross-frequency bilinear sum onto same-frequency moments — is the load-bearing open analytical step**, and any future T_N construction must do this derivation, not assert it. This is a future, more ambitious probe than R77 §10's "1–2 hours of focused implementation" suggested.

The probe doesn't decide among these three. It does establish that the natural next attempt in the project's own scheduled work (R77 §10) collapses to an open derivation the project's own ledger marks as open.

---

## Phase 1 detail

Full content in `T_N_OFF_LIN_SPEC.md`. Headlines:

1. **R77 sketch §5 is procedural** — quotes verbatim. The §5 sentence "Substituting Tao's recursion expresses them as quadratic forms in {P_n^{ab}(c)}" is the asserted-but-not-derived closure.

2. **Diagonal piece is rigorously closable.** The v = v' contributions reduce μ̂_n(ξ·2^{−v}) μ̂_n^*(ξ·2^{−v}) = |μ̂_n(ξ·2^{−v})|² (same frequency), which after ξ-sum and unit-shuffle change of variable lies in span{P_n^{ab}(c)}. This is exactly what `result_77_T_diagonal.py` derives → T_diag = (1/5)·[[1, 1], [4, 4]].

3. **Off-diagonal piece is not rigorously closable from §5.** The v ≠ v' contributions involve μ̂_n(u) μ̂_n^*(u') with u = ξ·2^{−v} ≠ u' = ξ·2^{−v'}, plus a non-trivial phase factor e^{−2πi ξ (2^{−v} − 2^{−v'})/3^{n+1}}. These do NOT collapse onto span{P_n^{ab}(c)} without an additional Plancherel-type closure that §5 does not provide.

4. **The project's own documentation flags this open.** `result_77_T_lead_spectrum.md` §3 offers TWO contradictory heuristic derivations (4·(1/4)=1 vs P(v=1)=1/2) and explicitly admits "more careful analysis needed". §6 ledger places "off-diagonal exact bilinear-sum analysis" on the OPEN side.

5. **No data-driven bridge is admissible.** Pre-registration requires decisions over Q (exact arithmetic). Numerical fits from (P_+, P_−) data at k = 2..5 produce approximations, not Fraction values; selecting among them is reconstruction.

Phases 2–5 are not executed; the gating condition was not met.

---

## Routing recommendation

Three plausible next moves, surfaced for Nathan (not chosen by this probe):

### Route 1: Derive the cross-frequency closure (substantive analytical work)

Make the §5 assertion "quadratic forms in {P_n^{ab}(c)}" rigorous by:

1. Expanding Q_n^{ab}(c; v, v') for v ≠ v' explicitly. The non-trivial phase factor 2^{−v} − 2^{−v'} has computable 3-adic valuation (e.g., 2^{−1} − 2^{−3} = 3/8, valuation 1; 2^{−1} − 2^{−5} = (16−1)/32 = 15/32, valuation 1; etc.).

2. Identifying which 3-adic-character orthogonalities collapse the cross-frequency bilinear onto same-frequency moments. The leading case (v = 1, v' = 3) gives valuation 1, suggesting the phase character sums to a level-n character with non-trivial structure — but exactly what that produces in terms of {P_n^{ab}(c)} is the open derivation.

3. If the closure exists: write Off_lin as a 2x2 over Q; spectrum computation is mechanical. If it doesn't: §5's assertion is wrong, T does NOT act on the 2D (P_+, P_−) subspace as a 2x2, and the operator-theoretic framework needs reframing onto the larger 6-dim or 8-real-dim space.

Cost estimate: this is real analytical work, of the same scale as the original T_diag derivation. Not "1–2 hours of focused implementation"; more like a multi-session probe with phases like CROSS_FREQUENCY_PHASE1_EXPANSION, CROSS_FREQUENCY_PHASE2_CHARACTER_SUMS, etc.

Recommended if a positive operator-theoretic anchor for rate-1/2 is wanted.

### Route 2: Pivot to R77.6 branch-cut / density-of-states framing

If the cross-frequency closure does NOT exist (i.e., the off-diagonal sum genuinely doesn't reduce to same-frequency moments), the rate-1/2 phenomenon may not live in any spectrum. Pivoting to a generating-function picture (R77.6) would:

- Treat S_n's generating function G(z) = Σ S_n z^n as the analytic object.
- Locate rate-1/2 as a singularity at z = 2 (radius-of-convergence boundary) — a branch-cut or pole feature.
- Bypass the operator-existence question entirely.

Cost: research-direction pivot, multi-session. The user's `project_collatz_r78_bilinear_cracked` memory mentions "Bilinear bound standalone publishable" — this Route 2 reframing of the spectral side might also stand alone publishable as a "density-of-states picture of c = 7/45's rate".

Recommended if continuing to search for a spectral carrier is no longer the priority.

### Route 3: Recognize the pattern and document the structural boundary

Five probes, no positive spectral carrier. Each probe has independently demonstrated that the natural candidate operator either doesn't exist over Q in the project, or does exist but doesn't have the spectrum the rate-1/2 conjecture requires.

Publishing this as a "structural boundary of the c = 7/45 spectral program" — alongside the bilinear bound and the F̂_p theorem that ARE rigorous — may be the most honest framing. The user's memory `project_collatz_r78_bilinear_cracked` already notes this: "c=7/45 Nisoli NOT closable" via R77.3 / R77.4 erratum routes. T_N construction joins this list.

Cost: a synthesis writeup, single session. The negative finding is publishable as a no-go theorem of sorts: "Within the project's bilinear pair-form moment framework, no Q-constructable 2x2 operator carries eigenvalue 1/2 on the (1, 4) deviation direction."

Recommended if the project's main publishable output is to be the bilinear bound side (already 25-commit, 183KB, 34-file burgess.zip per user memory).

---

## Adversarial check outcomes (probe-level)

**(A1) R77 sketch §5 fidelity.** Verbatim quoting in T_N_OFF_LIN_SPEC.md §(a); the "quadratic forms in {P_n^{ab}(c)}" assertion is §5's own wording, identified as the open closure. ✓

**(A2) Spectrum over Q.** Not reached — Phase 2 not executed. ✓ (vacuous gating).

**(A3) (1, 4) eigenvector check.** Not reached — Phase 2 not executed. ✓ (vacuous gating).

**(A4) Conflict with prior negatives.** No new conflict introduced. The H_OFF_LIN_UNDERSPECIFIED disposition is structurally consistent with R_K's H_R_K_INTRACTABLE (operator unconstructable) and Candidate A's F2 (single-finest-scale dominance) and R76 §11's INCONCLUSIVE (operator's spectrum claim is conjecture). No probe-internal conflict. ✓

**(A5) Resolvent norm honesty.** Not reached — Phase 4 not executed. ✓ (vacuous gating).

---

## Deliverables produced

- **T_N_OFF_LIN_SPEC.md** — Phase 1 spec extraction, verbatim §5 quoting + diagnosis + decision rationale
- **T_N_DISPOSITION.md** (this file) — top-level disposition + structural meaning + routing

Phases 2 (T_N_CONSTRUCTION.md), 3 (T_N_SPECTRUM.md), 4 (M_3_DOUBLEPRIME.md), 5 (T_N_CLOSURE_TABLE.md) deliverables are **not produced**. The pre-registered gating condition ("Phase 1 confirms R77 sketch §5 has the content as a 2x2 matrix over Q") was not met; the parent task's pre-registration explicitly says: "If Off_lin isn't specified enough to build, return H_OFF_LIN_UNDERSPECIFIED. Don't reconstruct."

A reproducible script `t_n_compute.py` is also not produced, because no T_N matrix exists to compute over.

---

## Files referenced

- `result_77_sketch.md` §5 — primary source for Off_lin's procedural articulation
- `result_77_T_lead_spectrum.md` §1, §2, §3, §6 — T_diag rigorous (§1), off-diagonal empirical (§2), heuristic candidate derivations of λ_2 (§3), open ledger (§6)
- `result_77_T_diagonal.py` — rigorous T_diag derivation (v = v' diagonal collapse)
- `T_lead_2x2.py` — numerical 2x2 fit from (P_+, P_−) data (not Q-derivation)
- `result_77_2_T_N_construction.py` — 1-dim T_N kappa_N flavor (data-driven)
- `result_76_conservation_law.md` §11 — original R76 §11 source
- `R76_S11_VERIFICATION.md` / `R76_S11_DISPOSITION.md` — prior probe, this one Fork 1 from
- `R_K_DISPOSITION.md` — same reconstruction-vs-verification failure mode, different operator
- `CANDIDATE_A_DISPOSITION.md` — F2 single-finest-scale dominance finding
- `M3_CLOSURE_TABLE.md` — closure inequality framework (not reached this probe)
- `result_77_2_nisoli_certification.md` — Nisoli closure structure (not reached this probe)
