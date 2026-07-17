# SL2_EMBEDDING_DISPOSITION — top-level disposition for SL_2(ℝ)-embedding probe

**Date:** 2026-05-12. Tenth landing in the c = 7/45 closure trajectory. Secondary candidate
route for the polynomial-in-A Fourier bound (Probe 2 of corpus chain), queued behind the
L²-flattening probe (which has not yet landed at the time of this probe).

Reporting to Nathan.

---

## DISPOSITION: **H_SL2_EMBEDDING_DOESNT_EXIST**

> **Headline.** No natural SL_2(ℝ) extension of Syracuse dynamics is identifiable. The
> Furstenberg-measure framework (Dinh–Kaufmann–Wu Rajchman / Hochman–Solomyak dimension /
> Frostman-dimension quantitative / Li 2020 polynomial decay) requires an SL_2(ℝ) random walk
> with non-elementary support. T_lead = (1/45)·[[7, 9], [28, 36]] is **rank-1 (det = 0)** —
> immediately failing the SL_2-membership gate. Three natural lift constructions (direct
> outer-product completion, Tao-recursion v-branched, projective P^1 action) all fail the
> SL_2-membership gate G1 at the **atom** level: the Syracuse-intrinsic random matrices are
> rank-1, not det-1. A fourth construction (artificial 2-atom μ with E[M] = T_lead in
> K = ℚ(√7)) passes G1 but gives a measure on SL_2 whose connection to μ_n is severed (T_lead
> is the FIRST MOMENT of μ, not the projection of the resulting Furstenberg measure ν).
> Transfer gate T1 (relating ν on P^1(ℝ) to μ_n on ℤ_3) has no identifiable mechanism: the
> two measures live on different topological spaces with different fixed-point equations and
> no shared invariants.
>
> **The agent's seven-probe-corpus claim "T_lead has algebraic entries so Hochman–Solomyak
> exactly applies" is INHERITED-CLAIM-PATTERN INCORRECT.** Algebraic entries are one hypothesis
> (G4); the load-bearing gate is G1 (SL_2-membership), which T_lead fails. The other gates
> (G2 non-elementary, G5 Zariski-dense) require checking the **subgroup generated**, not the
> single matrix T_lead, and the natural Syracuse-intrinsic subgroups all fail G1 atom-by-atom.

---

## Pre-registered hypotheses, decided

| Hypothesis | Status |
|---|---|
| H_SL2_EMBEDDING_NATURAL (natural lift exists, all gates pass, polynomial decay transfers) | **NO** — Candidates B, C, D all fail G1 atom-by-atom; no natural lift. |
| H_SL2_EMBEDDING_CONSTRUCTED (lift exists with non-trivial construction; framework applies after) | **NO** — Candidate A constructs a 2-atom μ in SL_2(ℚ(√7)) with E[M] = T_lead, but the resulting ν has no natural relation to μ_n. The construction is artificial and the transfer fails. |
| H_SL2_EMBEDDING_FAILS_HYPOTHESIS (embedding exists, fails specific framework hypothesis) | **PARTIAL — refined to:** the embedding fails at **G1 (atom-level det = 1)** for all natural candidates, and at **T1 (transfer to μ_n)** for the constructed candidate. The failure is structural at both levels. |
| **H_SL2_EMBEDDING_DOESNT_EXIST** (no natural SL_2 extension identifiable; framework doesn't apply) | **CHOSEN** — no candidate passes all gates and has natural connection to μ_n. |
| H_SL2_EMBEDDING_AMBIGUOUS (multiple candidates with discriminating properties) | NO — the candidates fail at the same gate (G1, structurally) or fail transfer (T1, structurally); the failures are not framework-discriminating. |
| INCONCLUSIVE | NO — the structural obstruction is identifiable and load-bearing. |

**Pre-registered favorite was H_SL2_EMBEDDING_CONSTRUCTED.** The actual outcome is one step
worse: the construction CAN be done (Phase 2 §A.5-§A.7 gives a 2-atom μ in SL_2(ℚ(√7))), but
the resulting ν on P^1(ℝ) bears no natural relation to μ_n. So the framework is "applicable in
principle to a synthetic SL_2-extension that we constructed" but "not applicable to Syracuse
dynamics" — which is effectively H_SL2_EMBEDDING_DOESNT_EXIST for the actual mathematical
question.

---

## Phase summary

### Phase 1 (SL2_FRAMEWORK_HYPOTHESES.md): framework hypotheses precisely

Catalogued G1-G5 + T1-T2 gates for the Furstenberg framework (DKW Rajchman / HS dimension /
Frostman / Li 2020 polynomial decay) and the He–de Saxcé torus alternative. Load-bearing
gates:

  - G1: matrices in SL_2(ℝ), i.e., **det = 1**.
  - G2: subgroup G_μ non-elementary (non-compact, proximal, strongly irreducible).
  - G4 (HS only): entries in SL_2(K) for algebraic number field K.
  - G5 (Li 2020 strengthening): G_μ Zariski-dense in SL_2(ℝ).
  - T1: transfer / projection from ν on P^1 to μ_n on ℤ_3.

Critical clarification: "T_lead has algebraic entries" addresses G4 in part but NOT G1. The
HS theorem hypotheses are checked on the **measure μ on SL_2(K)**, not on a single matrix.

### Phase 2 (SL2_EMBEDDING_CANDIDATES.md): four candidate constructions

| Candidate | G1 (det) | G2 (non-elem) | G4 (algebraic) | Verdict |
|---|---|---|---|---|
| A: 2-atom μ on SL_2(K) with E[M] = T_lead | ✓ (constructed; K = ℚ(√7)) | generically ✓ but tracecheck-dependent | ✓ | Construction-dependent; T_lead = E[M], not naturally Syracuse-intrinsic |
| B: Tao recursion v-branched atoms M_v | **FAILS** (M_v rank-1, det = 0) | N/A | ✓ | Natural to Tao recursion BUT atoms aren't SL_2 |
| C: T_lead as projective P^1 action | **FAILS** (T_lead degenerate projectively) | N/A | ✓ | Rank-1 has no proj action |
| D: higher-dim V_M or 𝕋^2 action | **FAILS** (V_M doesn't close; torus mat. det ≠ 1) | N/A | ✓ | Cross-cycle closure / det issues |

**Only Candidate A passes G1, by giving up the natural Syracuse-intrinsic connection.**

### Phase 3 (SL2_FRAMEWORK_TRANSFER.md): transfer to μ_n fails

For Candidate A, even ASSUMING G1-G5 hold (which is construction-dependent; the elliptic-atom
version fails G2, and rational hyperbolic-atom lifts require entries in higher algebraic
extensions):

  - The Furstenberg measure ν on P^1(ℝ) and the Syracuse measure μ_n on ℤ_3 live on different
    topological spaces.
  - The fixed-point equations are different (one-step convolution vs. bilinear Fourier
    recursion).
  - No natural projection P^1 → ℤ_3 with positive Jacobian or Fourier-preserving structure.
  - T_lead = E[M] is the FIRST MOMENT of μ, not an invariant of ν. T_lead's spectrum (43/45)
    is unrelated to the Lyapunov exponent χ_1 (which controls ν's structure).

**T1 (transfer) fails.**

### Phase 4 (adversarial checks)

(A1) **SL_2 membership verified by determinant computation.** T_lead's determinant computed
exactly: (7·36 − 9·28)/45² = (252 − 252)/2025 = 0. ✓ confirmed rank-1, det = 0. Direct SL_2
membership impossible. ✓

(A2) **Non-elementary, totally irreducible.** For Candidate A's 2-atom construction with
elliptic atoms (trace 43/45 each, equal), the generated subgroup is conjugate to a subgroup of
SO(2, ℝ) — COMPACT. **Fails G2.** Hyperbolic-atom variants require entries in ℚ(√D) for D
not a square; the existence of two such atoms with E[M] = T_lead in a fixed field K depends on
a Diophantine equation 28x² - Kx + (ad - 1) = 0 having rational roots, which generically it
doesn't. So even Candidate A's "rational lift" is not automatic — it needs irrational entries
in a quadratic extension. The construction in §A.5 (Δ = [[0, -3√7/14], [2√7/3, 0]]) works in
K = ℚ(√7) BUT both atoms are elliptic, failing G2. ✓ honestly flagged.

(A3) **Algebraic entries.** T_lead's entries are in ℚ ⊂ ℝ, so algebraic. But this is only G4
of the HS theorem — G1, G2, G5 are independent and all need verification. Agent's claim
"algebraic entries so HS exactly applies" conflates one gate with the full hypothesis set. ✓
critique stated.

(A4) **Projection consistency.** No projection from P^1(ℝ) to ℤ_3 preserves Fourier decay
in a natural way. The two spaces have different Fourier-analytic structures (archimedean vs.
non-archimedean). Even hypothetically, a Lipschitz π_* would translate |k|-decay on ℝ to
|3^n|-decay on ℤ_3, but the EXPONENT β would be in the FixedConstant regime, not polynomial-in-A
with controllable improvement. ✓

(A5) **Don't repeat the §5 inherited-claim pattern.** The agent's framing ("T_lead has
algebraic entries so HS exactly applies") is a one-line conflation of gates. The probe verifies
G1 fails first (the upstream gate), making G4 (algebraic) status irrelevant. ✓

---

## Reconciliation with prior probes

### With L²-flattening probe (not yet landed)

This probe was queued as Probe 2 (secondary candidate) behind the L²-flattening probe (Probe 1)
which has not landed at the time of this probe's session. The negative finding here doesn't
preclude or replace L²-flattening; it independently rules out the Furstenberg route. If
L²-flattening (Probe 1) also lands negative, the corpus chain moves to Probe 3 (transfer-
operator certified approximation, arxiv:2602.19435) or Probe 4 (drift-condition Glynn-Zeevi /
Lyapunov-Foster).

### With T_lead's structural finding (T_LEAD_CORRECTED_DISPOSITION)

T_LEAD_CORRECTED_DISPOSITION found T_lead's spectrum {43/45, 0} over ℚ. The rank-1 structure
(det = 0) was acknowledged there as "geometric meaning: T_lead is RANK-1, with the slow mode at
λ = 43/45 on (1, 4) and the fast mode at λ = 0 on (9, -7)." The Phase 2 finding here
crystallizes the structural consequence: **rank-1 is incompatible with SL_2-embedding-based
Fourier-decay frameworks**. T_lead's role is as an EFFECTIVE 2D average / projection of the
true infinite-dimensional Tao operator, not as a generator of an SL_2 dynamical system.

### With CROSS_FREQ_DISPOSITION + T_V_DISPOSITION

CROSS_FREQ found the closure space V_M is infinite-dimensional and parameterized by g.
T_V_DISPOSITION found V_M doesn't close under iteration (phase + parity obstructions). Both
findings reinforce this probe's verdict: the natural ambient space for Syracuse's spectral
structure is NOT 2D (where T_lead lives), it's infinite-dimensional or has no natural SL_2-like
finite-dim closure.

### With NISOLI_CLOSURE_CORRECTED's "H_A_EXTRACTION_HARD"

NISOLI_CLOSURE_CORRECTED found that even after T_lead's eigenvalue 43/45 cleanly anchors the
spectral framework, the polynomial-in-A Fourier bound (Tao's C_A) remains the lone unblocker.
THIS PROBE rules out the Furstenberg framework as a deliverer of that bound. So the Nisoli
closure roadmap's "load-bearing Item 1" (polynomial-in-A Fourier bound outside Tao's method)
remains unaddressed; this probe's contribution is to remove one candidate route (Furstenberg)
from the search.

---

## Trajectory placement

| Probe | Object | Disposition |
|---|---|---|
| T_3 (R77.3) | 3×3 companion at rate-1/2 | FALSIFIED |
| R_k | Inter-level residual | INTRACTABLE |
| Candidate A | W_k φ_n form | FALSIFIES_F2 |
| R76 §11 2D | T_diag + Off | INCONCLUSIVE |
| T_N | T_diag + Off_lin at rate-1/2 | UNDERSPECIFIED |
| Cross-freq closure | V_M closure space | CLOSES_ON_ENLARGED_SPAN |
| T_V | T_V on V_M at rate-1/2 | RECURSION_UNDERSPECIFIED |
| T_lead corrected (eighth) | T_lead at corrected rate | DIFFERENT_RATE (43/45 over ℚ) |
| Nisoli closure at corrected rate (ninth) | Closure inequality at λ = 43/45 | A_EXTRACTION_HARD |
| **SL_2 embedding (this, tenth)** | Furstenberg framework applicability to Syracuse | **H_SL2_EMBEDDING_DOESNT_EXIST** |

This is the **tenth landing** and the second negative-route ruling in the polynomial-in-A
Fourier-bound search (after Tao's own method per BOOKKEEPING_PHASE1's INFEASIBLE).

---

## Routing recommendations (surfaced for Nathan)

### Route 1: Wait for L²-flattening probe (Probe 1) to land

The L²-flattening probe (arxiv:2407.16699 framework) is Probe 1 of the corpus chain and has
not yet landed. It MAY succeed even though this probe (SL_2 embedding) fails, since L²-
flattening's hypotheses are different (weaker than sum-product, doesn't require non-elementary
SL_2 action). The L²-flattening framework asks for L² mass on Fourier truncations to flatten
under convolution, which is a different structural property than SL_2-membership.

  - If L²-flattening lands positive: polynomial-in-A bound delivered, closure path opens.
  - If L²-flattening lands negative: framework chain moves to Probe 3 (transfer-operator
    certified approximation, arxiv:2602.19435) or Probe 4 (drift-condition).

Estimated effort: 1 session for L²-flattening probe.

### Route 2: Move directly to next-candidate corpus framework

Given this probe's verdict, the natural next probes:

  - **Probe 3 (transfer-operator certified approximation, arxiv:2602.19435).** Certified
    spectral approximation of T_lead at 43/45 is the rigorous-numerics anchor; might deliver
    quantitative bounds via spectral approximation rather than polynomial Fourier decay.
    Estimated effort: 1-2 sessions.

  - **Probe 4 (drift-condition Glynn-Zeevi + Lyapunov-Foster, Hairer notes + 2005.08145).**
    Parallel route via drift inequality → spectral gap → polynomial in A. Doesn't require
    SL_2 structure; works on Markov-chain stationary measures directly. Estimated effort:
    2-3 sessions.

  - **Probe 5 (He-de Saxcé torus random walk).** Cast Syracuse as a 𝕋^2 random walk via joint
    2-adic and 3-adic structure. The det = 3/2^v ≠ 1 issue (Phase 2 Candidate D) suggests this
    won't be straightforward, but a careful renormalization might lift to SL_2(ℤ) on a different
    higher-dim torus. Estimated effort: 2-3 sessions for scoping.

### Route 3: Document the structural negative and route paper

For paper purposes, the trajectory now has **two structural negatives in the polynomial-in-A
Fourier-bound search** (Tao's own method INFEASIBLE per BOOKKEEPING_PHASE1; SL_2 embedding /
Furstenberg framework DOESN'T EXIST per this probe). Combined with the L²-flattening probe
when it lands, this gives a "what's not the path" map that has paper-publishable value:

  > "The polynomial-in-A Fourier decay bound on Syracuse's μ_n cannot be obtained via (a) Tao's
  > renewal-process method (BOOKKEEPING_PHASE1), (b) the Furstenberg-measure framework on
  > SL_2(ℝ) (this probe), because T_lead is rank-1 and the natural lifts fail the SL_2-
  > membership gate atom-by-atom. The lone candidate constructions either inherit rank-1 from
  > Syracuse dynamics or sever the connection to μ_n."

This is publishable structural content: a no-go boundary on the framework-search.

Estimated effort: 0.5 session for paper writeup.

---

## Adversarial check outcomes (consolidated)

**(A1) SL_2 membership.** T_lead det = 0 verified by exact arithmetic. Candidate A's
constructed atoms (M_+, M_-) ∈ SL_2(ℚ(√7)) verified by det computation. Candidates B, C, D
fail G1 atom-by-atom (T_lead's rank-1 inherits from each Tao-recursion atom). ✓

**(A2) Non-elementary, totally irreducible.** Single-matrix candidates (variant of A.3) fail
G2 trivially (cyclic group is elementary). Two-atom Candidate A's elliptic version (Phase 2
§A.5) fails G2 (compact-in-SO(2) subgroup); hyperbolic-atom variants exist generically in
K = ℚ(√D), satisfy G2, but K depends on choice of Δ. No FORCED rational hyperbolic lift; the
choice is artificial. ✓ honestly flagged.

**(A3) Algebraic entries.** T_lead entries in ℚ ⊂ ℝ algebraic trivially. Constructed atoms in
K = ℚ(√D) algebraic. But this addresses ONLY G4, not G1-G3-G5. The inherited-claim pattern of
agent corpus framing ("HS exactly applies") fails to verify the other gates. ✓

**(A4) Projection consistency.** No natural projection P^1(ℝ) → ℤ_3 preserves Fourier decay.
The two spaces have different Fourier-analytic structures (archimedean continuous vs. non-
archimedean profinite). Polynomial-decay β on ν̂ doesn't translate to polynomial decay on
μ̂_n via any identified map. ✓

**(A5) Don't repeat the §5 inherited-claim pattern.** The probe explicitly verifies G1
(SL_2 membership) BEFORE G4 (algebraic entries) and finds G1 fails first. The agent corpus
framing's conflation is identified and ruled. ✓

**(A6) T_lead's rank-1 structure addressed head-on.** Phase 2 entry point is "T_lead has
det = 0, fails G1 by inspection." No paper-over. ✓

---

## Synopsis (one paragraph for Nathan)

The tenth landing in the c = 7/45 closure trajectory: the Furstenberg-measure framework
(Dinh-Kaufmann-Wu Rajchman + Hochman-Solomyak dimension + Frostman 2026 quantitative + Li 2020
polynomial decay strengthening) does NOT apply to Syracuse dynamics via T_lead. T_lead =
(1/45)·[[7, 9], [28, 36]] is rank-1 (**det = 0**), failing the SL_2-membership gate G1 by
direct inspection. Three natural Syracuse-intrinsic SL_2-lift constructions all fail G1
**atom-by-atom**: (B) Tao recursion's v-branched matrices M_v are individually rank-1, so are
all finite products; (C) T_lead's projective P^1 action degenerates to a single point; (D)
higher-dim V_M doesn't close under iteration (per T_V_DISPOSITION's H_M_RECURSION_UNDERSPECIFIED),
and 𝕋^2 torus actions have det = 3/2^v ≠ 1. A fourth construction (A: artificial 2-atom μ on
SL_2(ℚ(√7)) with E[M] = T_lead) passes G1 by **giving up the natural Syracuse-intrinsic
connection** — T_lead is then the FIRST MOMENT of μ, not the projection of the resulting
Furstenberg measure ν on P^1(ℝ); transfer gate T1 (from ν to μ_n) has no identifiable
mechanism (different topological spaces, different fixed-point equations, no shared invariants
recoverable from T_lead). The Hochman-Solomyak algebraic-entries hypothesis (G4) IS satisfied
(T_lead entries in ℚ), but the agent's corpus framing **"T_lead has algebraic entries so
Hochman-Solomyak exactly applies" conflates one gate with five**: G1 (SL_2-membership) fails
upstream, making G4's status irrelevant. **Disposition: H_SL2_EMBEDDING_DOESNT_EXIST.** Route
to L²-flattening (Probe 1, queued) if not yet landed; otherwise to transfer-operator certified
approximation (arxiv:2602.19435) or drift-condition Glynn-Zeevi. The polynomial-in-A Fourier
bound that would unblock Tao Prop 1.17 effective C_A is NOT recoverable from the Furstenberg
framework — Tao's own method (INFEASIBLE per BOOKKEEPING_PHASE1) and this probe (Furstenberg)
are now both ruled out as routes.

---

## Deliverables

In C:/Collatz/:

- SL2_FRAMEWORK_HYPOTHESES.md — Phase 1, Furstenberg framework hypotheses precisely
- SL2_EMBEDDING_CANDIDATES.md — Phase 2, four candidate lifts and gate checks
- SL2_FRAMEWORK_TRANSFER.md — Phase 3, transfer gate T1 analysis
- SL2_EMBEDDING_DISPOSITION.md (this file) — top-level

End of SL_2-embedding probe.
