# R76_S11_VERIFICATION — Phase 1 content verification

**Date:** 2026-05-12. Phase 1 of the R76 §11 (P_+, P_−) 2D class-resolved recursion probe.
Wilson (analyst) reporting to Nathan.

---

## Headline

> **R76 §11 NAMES the (1, 4)-eigendirection with eigenvalue 1/2 — but as an OPEN STRUCTURAL CONJECTURE, not as a derived result.** The 2D class-resolved (P_+, P_−) operator IS articulated explicitly. The (1, 4)-eigendirection IS the asymptotic deviation direction (verified in §11 and in the companion R77 sketch). The eigenvalue 1/2 claim on that direction is **explicitly flagged "Open" / "structural conjecture"** in R76 §11's own closing paragraph. The rigorously-derived version of the 2D operator (T_diag from result_77_T_lead_spectrum.md §1) has spectrum {0, 1} on (1, 4), NOT spectrum {1/2, ...} on (1, 4).
>
> The Candidate A agent's lateral suggestion is therefore **half-correct**: the document names the 2D structure and the (1, 4)-eigendirection rigorously, but it names the 1/2-eigenvalue claim only conjecturally, and the rigorous version of the operator demonstrably does not deliver eigenvalue 1/2 on (1, 4) — it delivers eigenvalue 1.

This puts the probe at a structural boundary: not a clean "content present, proceed to Phase 2" (because the eigenvalue 1/2 claim is conjectural in R76 §11 itself), and not a clean "content not present" (because the document does articulate the structure and conjecture). The cleanest reading is **INCONCLUSIVE with respect to the pre-registered phrasing** ("R76 §11 names the (1, 4)-eigendirection with eigenvalue 1/2 as the rate carrier") — yes, it names it; no, it doesn't establish it.

The probe's substantive question — "can the 2D (P_+, P_−) operator carry rate-1/2 via the (1, 4)-eigendirection?" — is the **same open question R76 §11 itself ends with**. Phase 2 would be RECONSTRUCTING the (claimed) full operator T = T_diag + Off and trying to do what R76 §11's closing paragraph ("rigorous derivation of (1,4) eigenvalue = 1/2") explicitly says is open. That is not "verifying content"; that is "doing the open work R76 §11 declares open."

---

## (a) Exact text of R76 §11 — verbatim

R76 §11 is titled **"Class-resolved structural collapse (towards Result 77)"** in `result_76_conservation_law.md`. Three paragraphs carry the content:

**Para 1 (the structural collapse):**

> Computing the class-resolved bilinear moments P^{ab}(c) for (a,b) ∈ {+,−}², c ∈ {1,2}:
>
> > **For all n ≥ 2: P^{+−}(c) = 0 and P^{++}(1) = P^{++}(2), P^{−−}(1) = P^{−−}(2).**
>
> Cross-class moments vanish; class-c-symmetry holds exactly. Reduces 8-dim P-space to 2 free parameters (P_+, P_−) for n ≥ 2.

**Para 2 (the asymptotic targets and the (1, 4) direction):**

> **Asymptotic targets:** P_+ → 7/150, P_− → 14/75 = 28/150 (ratio 1:4 = (1/3)²:(2/3)² = squared class-mass ratio from R64.B).
>
> **Deviation direction:** (P_+ − 7/150, P_− − 14/75) is exactly proportional to **(1, 4)** at all observed levels, i.e., the deviation lives on a 1D subspace within the 2D (P_+, P_−) plane. The (1, 4) eigenvector preserves the squared class-mass ratio.

**Para 3 (the 1/2-eigenvalue claim and its caveat):**

> This 1D structural mode has eigenvalue **1/2** under the Tao recursion (the rate-½ identification). The orthogonal mode (breaking the 1:4 ratio) decays at faster rate (suppressed in observed data).
>
> **For full Result 77 closure:** derive the (1, 4) eigenvalue analytically from Tao's recursion combined with the class-mass conservation law (R66) and Plancherel structure. The key inputs are:
> - Asymptotic class fractions (1/3, 2/3) from R64.B
> - Mod-3 class transition rule from R66 (v even → class 1, v odd → class 2)
> - P(v even) = P(v odd) = 1/2 under Geom(2)
>
> The 1/2 rate emerges from P(v even) = 1/2, i.e., **the eigenvalue 1/2 of the Tao-recursion operator on the (1, 4)-eigendirection equals exactly the probability of v being even.**
>
> **Open:** rigorous derivation of leading coefficient 1/30 (numerical fit) and (1,4) eigenvalue = 1/2 (structural conjecture). Both reduce to algebraic identities from R66's chain dynamics.

That is the complete content of R76 §11.

---

## (b) Diagnosis — did the Candidate A agent overinterpret?

**Partially.** Here is what the Candidate A agent (CANDIDATE_A_DISPOSITION.md §"Alternative recommendation") actually wrote:

> A more substantive alternative: use **R76's class-resolved decomposition** (R76 §11) where π_n splits into π_{+, n} + π_{−, n} on the two mod-3 classes, and the (1, 4)-eigendirection is the rate-1/2 carrier. That framework is finite-dimensional (2D at each n) and might be a better testing ground for "where does rate-1/2 live spectrally" than the L²(Ẑ_3^×) framework. The R76 §11 framework was the original "operator at finite truncation" framing of R77.2; R77.4 ruled out the K_n form of that operator, but the **(P_+, P_−) 2D recursion remains an active anchor**.

The Candidate A agent's framing — "the (1, 4)-eigendirection is the rate-1/2 carrier" — is taken DIRECTLY from R76 §11's text. The agent did NOT invent the (1, 4)-eigendirection-with-eigenvalue-1/2 claim. R76 §11 names it explicitly.

What the Candidate A agent did NOT flag is that R76 §11 itself labels this claim as **conjectural and open**. The agent's "active anchor" phrasing reads as "this is an unsettled but live anchor", which is accurate to the source; but the parent task's pre-registration ("R76 §11 names ... as the rate carrier" / "the (1, 4)-eigendirection with eigenvalue 1/2 as R76 §11 claims") interprets the agent's phrasing as **content explicitly established**, not **content explicitly conjectured**.

So:
- R76 §11 has the **2D structure** explicitly (rigorous): the structural collapse P^{+−} = 0 + class-c-symmetry, reducing to 2 free parameters (P_+, P_−).
- R76 §11 has the **(1, 4)-eigendirection** explicitly (rigorous-as-observation through k=6, structural from R64.B + R66): the deviation lives on the (1, 4) subspace.
- R76 §11 has the **eigenvalue 1/2 on (1, 4)** claim explicitly but **conjecturally** ("structural conjecture", "open").
- R76 §11 does NOT have a rigorously-derived 2D operator with eigenvalue 1/2 in the (1, 4)-eigendirection.

---

## (c) Cross-reference: what the rigorous version of the 2D operator actually shows

The R77 sketch (`result_77_T_lead_spectrum.md`, also `result_77_T_diagonal.py` and `T_lead_2x2.py`) IS the project's attempt to rigorously construct the 2D operator R76 §11 conjectures about. That construction is incomplete in a specific, well-documented way:

> **Theorem 77.1 (rigorous):** The **diagonal-only** contribution of Tao's bilinear recursion to (P_+, P_−)_{n+1} is `T_diag = (1/5)·[[1, 1], [4, 4]]`.

Char poly of T_diag: `λ² − λ = 0`. Spectrum: **{0, 1}**, with:
- Eigenvector at λ = 1: **(1, 4)** — exactly the direction R76 §11 names
- Eigenvector at λ = 0: (1, −1)

So the eigendirection (1, 4) IS the slow-mode eigenvector of the rigorously-derived diagonal operator, but its eigenvalue is **1**, not 1/2.

The R77 sketch acknowledges this directly (verbatim, §2):

> T_diag alone gives S_{n+1} = S_n (eigenvalue 1 on (1, 4)). The actual S_n converges to 7/15 because of **off-diagonal corrections** ... where Off_n contains cross-frequency bilinear terms ...

And §3:

> **Conjecture 77.2:** The full operator T (T_diag + Off_n linearization) has subdominant eigenvalue λ_2 = 1/2 acting on the (1, 4) deviation subspace.

Section 6 of `result_77_T_lead_spectrum.md` ("Ledger of what's rigorous vs. empirical") puts the eigenvalue-1/2 claim on the empirical side of the ledger:

> ### Open (analytical work for fully rigorous closure)
> - ✗ Off-diagonal exact bilinear-sum analysis to confirm λ_2 = 1/2 from Tao's recursion
> - ✗ Rigorous derivation that 1/30 = S_∞/14 (combinatorial origin of 14)
> - ✗ Nisoli Theorem 2.15 application to certify lift from finite truncation T_N to limit T

This is the same "Open" R76 §11 declares.

So the rigorously-constructed 2D operator from §11's structure has eigenvalues {0, 1}. The 1/2 eigenvalue would belong to a different operator — T_diag + (some off-diagonal correction) — which is precisely what R77 sketch's §6 admits is unconstructed and §7 schedules as outstanding work.

**The 2D operator R76 §11 conjectures about exists; the operator that actually carries eigenvalue 1/2 on (1, 4) does not yet exist in the project.**

---

## (d) Can a 2D spectral object be constructed from what R76 §11 + R77 sketch provide?

A 2D operator can be constructed numerically. The R77 §6 / §7 path is well-defined: at each finite level n compute (P_+, P_−)_n and (P_+, P_−)_{n+1} exactly over Q, fit a 2x2 linear recursion connecting them, take its eigenvalues. `T_lead_2x2.py` already does this (numerically, fitting from k=2..5 deviations) and reports:

> Fitted T (from k=2→3, k=3→4):
>   eigenvalues to be computed; the numerical fit is consistent with one eigenvalue → 1/2 as more levels are used

This is the same numerical pattern the Off_n correction shows (§2 of result_77_T_lead_spectrum.md): ratios converging to ≈ 0.503 by k=5→6.

BUT: this is **fit-from-numerical-data**, not **derived-from-Tao-recursion**. It is the same operation as the empirical rate identification, with the same status: **conjectural rate-1/2 backed by ε_n through k=6, no rigorous operator-theoretic anchor**.

Constructing a 2D operator from what R76 §11 has, with rigorous eigenvalue 1/2, requires:
1. Rigorously deriving the off-diagonal correction Off_n as a bilinear operator on (P_+, P_−) (R77 sketch §5).
2. Showing the resulting full operator T = T_diag + Off_lin has eigenvalue exactly 1/2 on (1, 4).

Both are explicitly open in R77 sketch §10 and R76 §11's closing line. The project does not currently have them.

If the probe's Phase 2 PROCEEDED by fitting a 2x2 operator from finite-level numerical data (the only currently-tractable approach), it would be doing **reconstruction**, not **verification** — exactly the pattern the parent task warns against ("The R_K probe failed in part because R77.4 erratum §1's articulation was ambiguous and the agent had to reconstruct rather than verify; that ambiguity led to the R_K intractability. Avoid repeating that pattern.").

---

## (e) Phase 1 disposition

The pre-registration provided six dispositions. The cleanest match here is **INCONCLUSIVE** with the specific rationale:

> R76 §11 has the structural articulation but the eigenvalue 1/2 claim is articulated only as conjecture. The rigorously-derived version (T_diag) has eigenvalue 1, not 1/2, on the (1, 4) direction. The operator that would carry eigenvalue 1/2 (T_diag + Off_lin) is exactly what R76 §11 declares open. Proceeding to Phase 2 would mean reconstructing rather than verifying — the exact R_K probe failure mode the parent task warns against.

This is NOT H_R76_S11_DOESNT_HAVE_CLAIMED_CONTENT — the (1, 4) eigendirection is clearly named, and the eigenvalue 1/2 claim is articulated (just conjecturally).

This is NOT H_2D_DOESNT_CARRY_RATE_AT_ALL — the 2D structure may well carry rate-1/2 once Off is constructed; the issue is not "it doesn't carry rate" but "the operator that would carry it is itself the open derivation".

This is also NOT a confident H_2D_CARRIES_RATE_* — those require an explicit operator with eigenvalue 1/2 in hand, which R76 §11 does not provide and the project does not currently have.

The honest reading is **INCONCLUSIVE — R76 §11 articulates the conjecture but does not establish the operator; Phase 2 construction collapses to the very open derivation R76 §11 declares open.**

---

## (f) What would unblock Phase 2

Two articulations of R76 §11 would make Phase 2 tractable as a verification rather than a reconstruction:

1. **The R77 sketch §5 off-diagonal bilinear sum, written down as a closed-form 2x2 matrix Off_lin acting on (δ_+, δ_−)_n.** With T_diag + Off_lin as a concrete 2x2 matrix over Q, Phase 2 reduces to computing the spectrum and verifying eigenvalue 1/2 on (1, 4). This is the path R77 sketch §10 schedules ("1-2 hours of focused implementation").
2. **A separately-articulated 2x2 operator on (δ_+, δ_−) whose eigenvalue 1/2 on (1, 4) is derived (not conjectured), e.g., from R66's chain dynamics + R64.B's class mass + Geom(2)'s P(v even) = 1/2.** This is the path R76 §11's closing paragraph gestures at ("Both reduce to algebraic identities from R66's chain dynamics") but does not execute.

Neither is currently in the project. Constructing either is the next probe's task, not this one's — this probe's Phase 1 is content verification, and the verification result is "the content is articulated as conjecture; the operator is not constructed."

---

## (g) Routing recommendation

The natural next moves (NOT executed in this probe — surfacing for Nathan's decision):

- **Execute R77 sketch §7's implementation outline** (1-2 sessions): rigorously derive the off-diagonal bilinear correction Off_lin as a 2x2 matrix, compute spectrum of T_diag + Off_lin, verify eigenvalue 1/2 on (1, 4). If this succeeds, Phase 2 of this probe becomes well-defined and can be re-run.
- **Or route to Candidate B (Kozyrev wavelets)** as CANDIDATE_A_DISPOSITION.md's primary recommendation. This sidesteps the R76 §11 conjecture entirely.
- **Or recognize the structural pattern**: every spectral probe so far (T_3, K_n, R_k, W_k via φ_n, R76 §11 2D) has terminated at "the rigorous operator either doesn't exist yet or doesn't have the spectrum needed". The four-probe trajectory's lesson may be that **the project's framework is not currently equipped to deliver an operator with eigenvalue 1/2**, and the rate-1/2 phenomenon may genuinely live in the **branch-cut/density** structure of R77.6 (not in any spectrum).

---

## Adversarial check outcomes for Phase 1

**(A1) R76 §11 fidelity.** Quoted verbatim from `result_76_conservation_law.md` §11. The "structural conjecture" + "Open" labels are R76 §11's own framing, not this probe's. ✓

**(A2) The (1, 4)-eigendirection specifically.** R76 §11 articulates the basis explicitly: "(P_+ − 7/150, P_− − 14/75) is exactly proportional to **(1, 4)**". So the components are the (P_+ − target, P_− − target) deviation vector in the natural 2D basis. The basis is rigorous. The eigenvector identification is rigorous (empirical through k=6 + structural from R64.B class mass ratio). The associated **eigenvalue** is the conjectural part. ✓ — eigenvector basis articulated; eigenvalue conjecture flagged.

**(A3) Eigenvalue 1/2 over Q vs approximately.** R76 §11 claims 1/2 as a structural target; the rigorously-derived T_diag has eigenvalue **1 exactly over Q** on (1, 4) (computable from char poly λ² − λ = 0). The 1/2 eigenvalue is observed numerically in `T_lead_2x2.py` fits and in the off-diagonal ratio convergence. No rigorous operator over Q yields eigenvalue Fraction(1, 2) on (1, 4) at this time. ✓ — discrepancy flagged.

**(A4) Relationship between operator and moment functional.** R76 §11 + R77 sketch articulate this: (P_+, P_−) and S_n = 2(P_+ + P_−) and ε_n = 10 δ_+ where (δ_+, δ_−) is the deviation from (7/150, 14/75). The (1, 4)-mode IS the ε_n-carrying mode by Plancherel. This part is rigorous. ✓

**(A5) Conflict with R77.3, R77.4, R77.5.** No direct conflict. The 2D (P_+, P_−) operator R76 §11 articulates is a **different object** from K_n (the within-level Markov operator R77.4 ruled out as a rate carrier), T_3 (the order-3 companion R77.3 falsified as a clean three-mode geometric ansatz), and R_k (the inter-level residual vector R77.5 demonstrated to be a vector not an operator). Specifically: R76 §11's 2D operator is the **outer-product bilinear transfer operator on class-resolved moments**, which is yet another object. The non-conflict is welcome — the prior negatives don't rule it out. The conflict-shaped issue is that R76 §11's own derivation of T_diag (the rigorous component of its 2D operator) gives eigenvalue 1, not 1/2, on (1, 4). That's not conflict with prior R77.x results; that's a structural gap inside R76 §11's own framework. ✓ — no inter-result conflict; intra-§11 gap surfaced.

---

## Decision

**Phase 1 disposition: INCONCLUSIVE (with specific rationale, not vague "ambiguity").**

Proceeding to Phases 2–5 would mean re-doing the open work R76 §11 itself names. That is reconstruction, not verification, and the parent task explicitly warns against it.

The next deliverable, R76_S11_DISPOSITION.md, will record this disposition and surface routing recommendations.

Phases 2, 3, 4 deliverables are NOT produced (the gating condition was not met). This is the pre-registered behavior: "Phase 1 is gating. Do not proceed to Phase 2 if R76 §11's content is unclear."
