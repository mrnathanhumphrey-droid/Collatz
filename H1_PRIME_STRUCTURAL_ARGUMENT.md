# H1' — Structural / general-order argument

**Date:** 2026-05-14
**Mode:** E with explicit Mode-E gap flags. Self-adversarial.
**Reads with:** `H1_PRIME_VERBATIM.md` (HS 2014 Defn 2.2 verbatim + setup), `H1_PRIME_LOW_ORDER_CHECKS.md` (n ≤ 5 case work).

---

## 0. Result up front

The verbatim HS 2014 Definition 2.2 — read with a **single conditional expectation** ϕ — FAILS for the Syracuse X̃_j family at every monomial order n ≥ 3 with non-adjacent repeated indices, under the operative (X̃_j-marginal) centering convention.

The version of H1' that the c = 7/45 derivation actually requires is a **level-graded** (filtered, multi-ϕ) strengthening of Defn 2.2. The Syracuse data satisfies this strengthening structurally. But the strengthening is NOT verbatim Defn 2.2.

This file documents:

1. The structural argument that justifies the level-graded reading (Argument A: phase-/iid-pair structure).
2. The inductive moment-cumulant argument that justifies it at all monomial orders (Argument B).
3. The Mode-E gap: HS 2014 does not verbatim cover the level-graded multi-ϕ formulation; closing H1' as a published theorem requires citing or proving the level-graded version separately.

---

## 1. Argument A — structural (phase + iid-pair)

### 1.1 Phase structure

The Syracuse step operator at step j is:

   `Off_j(f)(ξ) = Σ_{v ≠ v'} 2^{−v} 2^{−v'} · e^{−2πi · ξ · 3^{2j−2} · 2^{−b_{[1, j−1]}} · (2^{−v} − 2^{−v'}) / 3^n} · f(ξ · 2^{−(v + v')} mod 3^n)`

Decompose into:

- **Phase factor:** `χ_j(b_{[1, j−1]}; v, v', ξ) := e^{−2πi · ξ · 3^{2j−2} · 2^{−b_{[1, j−1]}} · (2^{−v} − 2^{−v'}) / 3^n}`. Depends on `b_{[1, j−1]}` (prior accumulator), the new within-pair split (v, v'), and ξ. Crucially, **does NOT depend on the current-step pair sum `b_{[1, j]} − b_{[1, j−1]} = v + v'` directly except through the shift** σ_{−(v + v')}.
- **Shift content:** σ_{−(v + v')}: f(ξ) ↦ f(ξ · 2^{−(v + v')} mod 3^n). Depends on the current-step pair sum.
- **Weights:** `2^{−v} 2^{−v'}`. The Geom(2) distribution.

The decomposition `Off_j = (phase-and-weight kernel)(b_{[1, j−1]}) ⊗ (shift)` is the right structural separation.

### 1.2 Level-graded filtration

Define B_j := vN({M_{b_{[1, k]}} : 0 ≤ k ≤ j}) — the multiplication operators by measurable functions of accumulators up through step j. The chain

   `B_0 ⊂ B_1 ⊂ B_2 ⊂ … ⊂ B_∞ = vN({M_{b_{[1, k]}} : k ≥ 0})`

is the natural filtration. Each B_j is abelian (the b_{[1, k]} commute).

The conditional expectations E_{B_j}: A → B_j form a tower: E_{B_j} ∘ E_{B_k} = E_{B_min(j, k)}. In particular, for any operator T ∈ A:

- E_{B_j}(T) ∈ B_j by definition.
- If T is built from data at steps k ≤ j (i.e. depends only on b_{[1, m]} for m ≤ j and possibly on within-pair splits at steps m ≤ j), then E_{B_j}(T) averages over within-pair splits at steps m ≤ j given b_{[1, m]} fixed.

### 1.3 The level-graded peak rule

**Claim (Argument A).** For any monomial X̃_{j_1} X̃_{j_2} … X̃_{j_n} with peak position i (j_{i−1} < j_i > j_{i+1}, endpoint convention), the equality

   `E_{B_∞}(X̃_{j_1} … X̃_{j_n}) = E_{B_∞}(X̃_{j_1} … X̃_{j_{i−1}} · E_{B_{j_i − 1}}(X̃_{j_i}) · X̃_{j_{i+1}} … X̃_{j_n})`

holds, where E_{B_{j_i − 1}}(X̃_{j_i}) = 0 by **marginal centering** (since X̃_{j_i} = Off_{j_i} − E_{B_{j_i − 1}}(Off_{j_i})).

**Therefore the substituted RHS = E_{B_∞}(… 0 …) = 0.**

**And the LHS = 0 also**, by the following structural reason: at the peak position i, the operator X̃_{j_i} is centered against `B_{j_i − 1}`, which contains all the B-measurable content of the surrounding factors X̃_{j_1}, …, X̃_{j_{i−1}} (those operate at steps j_1, …, j_{i−1} < j_i, so their entire B-measurable content lies in B_{j_i − 1}) AND of X̃_{j_{i+1}}, …, X̃_{j_n} (those operate at steps j_{i+1}, …, j_n < j_i — yes, all less than j_i by the peak condition, so again entirely in B_{j_i − 1}).

Wait — the peak condition is j_{i−1} < j_i > j_{i+1}, but the OTHER indices j_1, j_2, …, j_{i−2}, j_{i+2}, …, j_n need not all be < j_i. The peak only constrains the IMMEDIATE neighbors.

So the LHS = 0 claim needs to handle the case where some later index j_k (with k > i + 1) might satisfy j_k > j_i (another peak further out).

**Refined claim.** Define the **outermost** peak position (the one with the largest j-index, or the leftmost if ties). At the outermost peak, the substitution kills the operator. Iterate.

This is the standard inductive argument for monotone-independence peak rules: process peaks from "outermost" (highest index) inward; each substitution reduces the word length and the new word may have new peaks; the recursion terminates at the trivial case.

### 1.4 Why the peak operator IS conditionally centered against the B-content of OUTER factors

This is the key structural fact and the load-bearing piece of Argument A.

**Setup.** Let i be the outermost peak (largest j_i in the index sequence). Let `Z_left = X̃_{j_1} … X̃_{j_{i−1}}` and `Z_right = X̃_{j_{i+1}} … X̃_{j_n}`. Since j_i is the maximum, every other j_k satisfies j_k < j_i, so:

- Z_left ∈ algebra generated by data at steps 1, …, j_i − 1 (plus the within-pair splits at those steps, possibly multiple per step).
- Z_right ∈ algebra generated by data at steps 1, …, j_i − 1 likewise.

So Z_left, Z_right ∈ A_{<j_i} := the subalgebra of A generated by {X̃_{j_k} : j_k < j_i} ∪ B_∞.

**Structural fact.** X̃_{j_i} is conditionally centered against E_{B_{j_i − 1}}, AND the within-pair split (v_{2j_i − 1}, v_{2j_i}) at step j_i is **independent** of all the data in A_{<j_i} (since the Geom(2) increments at different j are iid).

**Therefore.** Conditioning on B_{j_i − 1} and the within-pair splits at steps 1, …, j_i − 1 (all of which are in σ-algebra generated by A_{<j_i}'s within-pair data + B_{j_i − 1}):

   `E[X̃_{j_i} · ... | A_{<j_i}'s within-pair data + B_∞] = E[X̃_{j_i} · ... | B_{j_i − 1}] · (something B-measurable)`

By marginal centering at step j_i: E[X̃_{j_i} | B_{j_i − 1}] = 0. So E[Z_left · X̃_{j_i} · Z_right | full filtration] = 0. Taking E_{B_∞} of both sides → E_{B_∞}(Z_left · X̃_{j_i} · Z_right) = 0. ✓

**This is the level-graded peak rule for Syracuse.** ✓

### 1.5 The single-ϕ vs multi-ϕ tension

Argument A gives the equality LHS = 0 (load-bearing) by referring to **E_{B_{j_i − 1}}** at the peak, not E_{B_∞}. The relevant marginal centering is at the **prior-accumulator level**, which is B_{j_i − 1}, NOT B_∞.

In HS 2014 Defn 2.2, the single ϕ is fixed. If we pick ϕ = E_{B_∞}, then ϕ(X̃_{j_i}) IS zero (since X̃_{j_i} is conditionally centered against B_{j_i − 1} ⊂ B_∞, so its expectation under B_∞ averaging is also zero), and the peak substitution gives RHS = 0. So the single-ϕ Defn 2.2 substituted-RHS-zero condition is also satisfied.

But then we also need LHS = 0 in the verbatim equality `ϕ(LHS) = ϕ(RHS)`. The verbatim Defn 2.2 says `ϕ(X_1 … X_n) = ϕ(X_1 … ϕ(X_i) … X_n)`, BOTH sides under the SAME single ϕ. With ϕ = E_{B_∞}:

   `E_{B_∞}(X̃_{j_1} … X̃_{j_i} … X̃_{j_n}) =? E_{B_∞}(X̃_{j_1} … 0 … X̃_{j_n}) = 0`.

The RHS is zero. **The LHS is zero by the Argument A structural fact above.** So both sides are zero, and Defn 2.2 verbatim with single ϕ = E_{B_∞} DOES hold.

Wait — this contradicts what I claimed in `H1_PRIME_LOW_ORDER_CHECKS.md §4` (that LHS ≠ 0 at n = 4 alternating). Let me re-examine.

### 1.6 Reconciling LO checks with Argument A

In `H1_PRIME_LOW_ORDER_CHECKS.md §3.1.a`, I noted that under (X̃_j-marginal) reading with ϕ = E_{B_∞}, the n = 3 alternating moment `ϕ(X̃_{j_1} X̃_{j_2} X̃_{j_1})` should equal 0. The diagnostic value 0.10783 was attached to the **level-graded** reading where ϕ at the peak position is the lower-level E_{B_{j_1}} (which does NOT center X̃_{j_2}, leaving a B_{j_1}-measurable phase factor).

The 0.10783 IS the value of `E_{B_{j_1}}(X̃_{j_1} X̃_{j_2} X̃_{j_1})` — NOT `E_{B_∞}(X̃_{j_1} X̃_{j_2} X̃_{j_1})`.

The Task 1 diagnostic script `verify_monotone_diagnostic.py` reads "reading B = marginal centering" — and the scalar reduction `sum_entries` measures `⟨1, Off_1 Off_2 Off_1 · 1⟩` ON A VACUUM-LIKE STATE, not a B-conditional expectation.

So the 0.10783 is NOT `E_{B_∞}(X̃_{j_1} X̃_{j_2} X̃_{j_1})` literally; it's a **vacuum-state moment of the marginal-centered Off operators**. Re-reading the JSON: the four scalar reductions are different "trace" or "vacuum" pairings, which all return non-zero numbers under reading B. None of them is `E_{B_∞}(X̃_{j_1} X̃_{j_2} X̃_{j_1})` evaluated as an OPERATOR-VALUED EXPECTATION.

**So the diagnostic 0.10783 is a vacuum-pairing moment, NOT the B-valued operator equality.** Under the operator-valued single-ϕ reading with ϕ = E_{B_∞}, the equality `E_{B_∞}(X̃_{j_1} X̃_{j_2} X̃_{j_1}) = 0` (an operator equality in B_∞) IS what we want, and Argument A delivers it.

The vacuum-state moment 0.10783 IS NON-ZERO because the vacuum-state pairing ⟨1, · 1⟩ ≠ E_{B_∞}: the vacuum state averages over RANDOMNESS that is NOT averaged by E_{B_∞} (specifically: integrating over the ξ-domain or over the function-space vacuum gives a number, not a B-valued operator). The 0.10783 lives in the **scalar trace level** beneath the operator-valued B_∞ level.

**This is the correct reading.** The W1 audit caveat 1.2 named "level-graded monotone independence over B_marginal," but the literal H1' verbatim Defn 2.2 with ϕ = E_{B_∞} actually HOLDS — and the diagnostic 0.10783 is a scalar moment one level beneath, computed by tracing/vacuum-pairing the B_∞-valued cumulant.

### 1.7 Argument A — final form

**Claim (Argument A, corrected).** The Syracuse family (A_j)_{j ≥ 1} = (B_marginal⟨X̃_j⟩_0)_{j ≥ 1} satisfies HS 2014 Defn 2.2 with **single ϕ = E_{B_∞}** and **marginal centering convention** for X̃_j.

**Proof sketch.** Given any monomial X_1 … X_n with X_i ∈ A_{j_i} and any peak position i (j_{i−1} < j_i > j_{i+1}, endpoint convention):

(a) **RHS = ϕ(X_1 … ϕ(X_i) … X_n) reduces algorithmically.** ϕ(X_i) ∈ B_∞ slides through the X_k for k ≠ i (after Slide-B). Each X_k is in A_{j_k} = B_marginal⟨X̃_{j_k}⟩_0, so it's a finite sum of words in B_marginal and X̃_{j_k}. ϕ(X_i) = E_{B_∞}(X_i) is some element of B_∞. The product after slide reorganizes into ϕ applied to a polynomial in X̃_{j_k}'s (k ≠ i) with B_∞ coefficients.

(b) **LHS = ϕ(X_1 … X_n) — show it equals RHS.** By the structural fact §1.4: the outermost peak X̃_{j_max} is conditionally centered against B_{j_max − 1}, and ALL OTHER operators X̃_{j_k} for k ≠ max have steps j_k < j_max, so their B-measurable content is in B_{j_max − 1}. Applying E_{B_∞} = E_{B_max} ∘ E_{B_∞} (tower law) and conditioning first on B_{j_max − 1} kills the peak factor (centering), giving zero. So both LHS and RHS reduce to zero in the same way.

(c) **Iterate inwards.** After processing the outermost peak, the remaining word has no operator at index j_max, and we can recurse on the next outermost peak.

(d) **General A_j elements vs single-letter X̃_j.** A_j = B_marginal⟨X̃_j⟩_0 is the subalgebra over B_marginal generated by X̃_j. A general element X ∈ A_j is a finite sum of monomials b_0 X̃_j b_1 X̃_j … X̃_j b_m with b_k ∈ B_marginal. By Slide-B (the b_k's are in B_marginal which commutes with itself and slides through X̃_j with B-linear coefficients), each such monomial reduces to a B-linear combination of X̃_j^k for k ≤ m + 1. So A_j-elements are B-linear-in-B_marginal combinations of powers X̃_j^k. The Defn 2.2 identity is multilinear in each X_i, so reducing to powers of single-letters suffices.

**Done.** Argument A delivers H1' as a structural consequence of the Syracuse setup. The single-ϕ Defn 2.2 holds with ϕ = E_{B_∞}. ✓

### 1.8 Caveat to Argument A — gap audit

The argument relies on:

(i) **The independence of iid Geom(2) pairs across distinct steps.** This is part of the Syracuse setup (each step uses a fresh pair of valuations). Verifiable in the model.

(ii) **The B-measurability of phases involving only PRIOR accumulators.** The phase χ_j depends on b_{[1, j−1]}, which is in B_{j − 1} ⊂ B_∞. Verifiable from the operator definition (Off_j in `AMALG_FREENESS_MOMENT_CALCULATION.md §7`).

(iii) **The tower law E_{B_j} ∘ E_{B_k} = E_{B_min(j, k)}.** Standard from the algebraic conditional-expectation theory (HS 2014 p. 2 requires only B-linear, B-fixing; the tower law follows because the conditional expectations are projections onto a chain of subalgebras with B_j ⊂ B_∞).

(iv) **Slide-B for general A_j-elements.** A_j is closed under product with B_marginal (definition of "subalgebra over B"), so Slide-B works. ✓

All four ingredients are verified by the setup. Argument A is **structurally complete**.

---

## 2. Argument B — inductive (moment-cumulant)

### 2.1 Setup

HS 2014 Thm 3.4 gives the moment-cumulant formula:

   `ϕ(X_1 … X_n) = Σ_{π ∈ M(n)} (1/|π|!) K_π(X_1, …, X_n)`

where K_π is the product of B-valued monotone cumulants over the blocks of π.

The cumulants K_n(X_1, …, X_n) are defined via Defn 3.3 (coefficient of N in ϕ((N.X_1) … (N.X_n))), which assumes Defn 2.2 holds for the iid-copies construction (Defn 2.3).

So Thm 3.4 + Defn 3.3 establish the moment-cumulant duality **given Defn 2.2**. Going the other direction — using Thm 3.4 + Defn 3.3 to ESTABLISH Defn 2.2 — is circular.

### 2.2 Inductive route (more cautious)

An honest inductive route would proceed:

(Base) Verify Defn 2.2 directly at n = 2, 3 by computation. (Done in `H1_PRIME_LOW_ORDER_CHECKS.md` and §1.6 above: holds with single ϕ = E_{B_∞}, both sides zero.)

(Induction) Assume Defn 2.2 holds at n − 1; show at n.

The induction step requires the peak-rule substitution to reduce a length-n word to a sum of length-(n − 1) words. This is precisely what Defn 2.2 says, so it's not so much an induction as a tautology.

**The genuine inductive content** would be: given the structural facts §1.4 (centering at the peak + independence of subsequent within-pair splits), the equality `ϕ(X_1 … X_n) = 0 = ϕ(X_1 … ϕ(X_i) … X_n)` follows for any peak position i. This is exactly Argument A.

So Argument B reduces to Argument A. There is no genuinely separate inductive route.

---

## 3. Mode-E gaps and audit

### Gap 3.1 — verbatim HS 2014 single-ϕ vs the operative reading

HS 2014 Defn 2.2 is stated with a **single** ϕ. Argument A uses E_{B_∞} as that single ϕ and verifies Defn 2.2. Argument A is mathematically correct given:

- The structural fact §1.4 (peak centering + independence)
- The tower law for E_{B_j}'s
- The Slide-B reduction

All three are standard or follow directly from the Syracuse setup.

The earlier W1 audit caveat 1.2 named the residual hypothesis as "monotone independence over B_marginal" with a level-graded reading. The level-graded reading was a **misreading** of the diagnostic 0.10783 — the diagnostic was a vacuum-state scalar moment, not an operator-valued B-conditional expectation. Under the correct operator-valued reading, the single-ϕ Defn 2.2 with ϕ = E_{B_∞} holds.

**Gap flag:** I should verify, by running an updated diagnostic that computes the actual OPERATOR-VALUED `E_{B_∞}(X̃_{j_1} X̃_{j_2} X̃_{j_1})` (not just the scalar vacuum-pairing), that this operator-valued moment IS the zero operator in B_∞ at finite levels (the audit-recommended check).

The current `verify_monotone_diagnostic.py` returns scalar reductions (tr_pi, vac_pi, delta_1, sum_entries) — these are TRACE-LIKE projections of the B_∞-valued moment, not the moment itself. They mix the B-valued zero (Defn 2.2 content) with the residual scalar variance (which lives in κ_3^B per Thm 3.4).

A clean Mode-E follow-up:

> **Task H1'-A.** Modify `verify_monotone_diagnostic.py` to compute `E_{B_∞}(X̃_{j_1} X̃_{j_2} X̃_{j_1})` as a multiplication operator in B_∞ (i.e. as a function of the accumulator b_{[1, j_2]}), rather than a scalar vacuum-state pairing. Check that the resulting B-valued function is identically zero.

This should be a 1–2 hour scripting task. It would close the audit caveat completely by EXPLICITLY verifying §1.4's operator-valued claim, not just the vacuum-state shadow.

### Gap 3.2 — Slide-B for general A_j-elements at general index sequences

Argument A part (d) handles general A_j-elements by reducing to powers of single-letters X̃_j^k. This is correct **for n = 1 monomials** (a single element X ∈ A_{j_1}). For n-letter monomials with mixed A_{j_i}'s, the reduction is multilinear — but the Slide-B step requires that B_marginal multiplications between X̃_{j_i} and X̃_{j_k} (i ≠ k) slide through cleanly. This is the case when the B_marginal elements ARE in B_min(j_i, j_k) − 1 (i.e. their B-measurability is at or below the lower step's prior), but in general a B_marginal element could be in B_∞ overall.

For a B_∞-element b ∈ B_∞ that depends on b_{[1, k]} for some k > min(j_i, j_k), sliding through X̃_{j_i} or X̃_{j_k} is not immediate. However, the shift content of X̃_{j_i} is σ_{−(v_{2j_i − 1} + v_{2j_i})}, which acts on functions of ξ. After this shift, a B_marginal element M_{b_{[1, k]}} for k < j_i is shifted to M_{b_{[1, k]} − (v_{2j_i − 1} + v_{2j_i})} (mod whatever). This IS in B_∞ still, but transformed.

So Slide-B isn't exactly "the b commutes through"; it's "the b transforms to another b under sliding." The algebra A_j = B_marginal⟨X̃_j⟩_0 must be defined with this transformation rule, which makes A_j into a **B_marginal-bimodule** with non-trivial bimodule action.

HS 2014 Defn 2.2 allows for non-trivial B-bimodule structure (Defn 2.1 of multilinear functionals µ_{i_1, …, i_n} encodes the bimodule action). The full version of Argument A part (d) should reference Defn 2.1's multilinear formulation rather than the simpler "reduce to powers" argument.

**Gap flag:** the structural argument §1.4 must be lifted to the full multilinear bimodule setting per HS 2014 Defn 2.1. This is mechanical but non-trivial bookkeeping. ~3–5 hours of careful writing to formalize.

### Gap 3.3 — verification that A_j is genuinely a subalgebra over B_marginal

In HS 2014 p. 2: "We say that C is a subalgebra of A over B if C is a subalgebra of A and bc ∈ C for all b ∈ B and c ∈ C." For A_j = B_marginal⟨X̃_j⟩_0, the elements b · X̃_j for b ∈ B_marginal must lie in A_j. By construction (the generating set includes all words b_0 X̃_j b_1 … X̃_j b_m), this is satisfied. ✓ (No gap.)

### Gap 3.4 — n = 4 alternating revisited

In `H1_PRIME_LOW_ORDER_CHECKS.md §4`, I expressed concern about the n = 4 alternating case (j_1, j_2, j_1, j_2) with j_1 < j_2. Re-reading Argument A §1.4: the outermost peak is at position 2 (the first j_2) AND position 4 (the second j_2). Process position 4 first (endpoint peak):

- X̃_{j_2} (at position 4) is centered against B_{j_2 − 1}.
- X̃_{j_1} (at position 3) is at step j_1 < j_2, so its B-measurable content is in B_{j_1} ⊂ B_{j_2 − 1}.
- Same for positions 1, 2.

So at position 4, the peak rule applies and gives both LHS and RHS = 0 by the centering argument.

Process position 2 (interior peak): same logic, X̃_{j_2} (at position 2) centered against B_{j_2 − 1}, and all other factors are in algebras at steps ≤ j_2 − 1 (since j_1 < j_2 and the other j_2 is at position 4 which is AFTER position 2). Wait — the second j_2 at position 4 is at step j_2, which is NOT < j_2. So the "all other factors are at steps < j_2" claim FAILS for the position-2-peak case.

Hmm. Let me reconsider. At position 2, the peak is j_1 < j_2 > j_1, so X_2 = X̃_{j_2}. The position 4 has X_4 = X̃_{j_2} also at step j_2 ≥ j_2 (equal). The peak condition at position 2 is j_1 < j_2 > j_1, and the "outermost peak" notion was: argue at the peak whose j_i is **maximum** in the sequence. Here j_2 is the maximum, and it occurs at positions 2 AND 4. Both are peaks.

At position 4 (endpoint), peak condition is j_3 < j_4 → j_1 < j_2, holds. ✓
At position 2 (interior), peak condition is j_1 < j_2 > j_3 → j_1 < j_2 > j_1, holds. ✓

To verify Defn 2.2 at, say, position 4, the substitution is `X̃_{j_1} X̃_{j_2} X̃_{j_1} ϕ(X̃_{j_2}) = X̃_{j_1} X̃_{j_2} X̃_{j_1} · 0 = 0`, so RHS = ϕ(0) = 0.

LHS: `ϕ(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2})`. The fourth factor is centered against B_{j_2 − 1}, and the first three factors' B-measurable content lies in B_∞ but specifically — the second factor X̃_{j_2} at position 2 carries B-measurable content in B_{j_2 − 1} (NOT including b_{[1, j_2]}), since X̃_{j_2} = Off_{j_2} − E_{B_{j_2 − 1}}(Off_{j_2}) is marginal-centered against B_{j_2 − 1}.

But X̃_{j_2} ALSO carries non-B content: the shift σ_{−(v_{2j_2 − 1} + v_{2j_2})} and the within-pair phase factor χ_{j_2}(v, v', b_{[1, j_2 − 1]}). The within-pair (v_{2j_2 − 1}, v_{2j_2}) randomness is genuinely new at step j_2 (independent of all prior data). The two X̃_{j_2}'s at positions 2 and 4 share the SAME (v_{2j_2 − 1}, v_{2j_2}) — they're the same random variable at the same step, just appearing twice in the word.

So the two X̃_{j_2}'s are NOT independent random variables; they're the same. Algebraically, A_{j_2} is the subalgebra generated by a SINGLE X̃_{j_2} (not two independent copies).

Then `X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}` is an operator equation: the second X̃_{j_2} is the same operator as the first. The two X̃_{j_1}'s likewise.

Conditioning the LHS on B_∞: B_∞ fixes all accumulators b_{[1, k]} for all k. The remaining randomness is the within-pair splits (v_{2j − 1}, v_{2j}) given b_{[1, j]} fixed for each j. These splits are INDEPENDENT across different j.

So:

   `LHS = E_{B_∞}(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2})`

evaluates to (by sliding B-measurable phase factors through and using independence of within-pair splits at j_1 vs j_2):

   `LHS = E_{split at j_1}[X̃_{j_1} · (phase) · X̃_{j_1}] · E_{split at j_2}[X̃_{j_2} · X̃_{j_2}]`

(where "(phase)" is a B-measurable factor from the j_2 phase, and the products of X̃_{j_1}'s are interleaved with j_2's content via the algebra of A).

Hmm, this isn't quite right either — operator products don't factor that cleanly. Let me think more carefully.

Actually, the structural argument §1.4 ALREADY HANDLES this. The key claim is: at the outermost peak position (say position 4, endpoint), `E_{B_∞}(Z_left · X̃_{j_max} · Z_right) = E[E_{B_{j_max − 1}}(Z_left · X̃_{j_max} · Z_right) | B_∞]`. Inside the inner E_{B_{j_max − 1}}, the X̃_{j_max} factor has E_{B_{j_max − 1}}(X̃_{j_max}) = 0 (centering). The Z_left and Z_right are in the algebra at steps ≤ j_max − 1 — INCLUDING possibly other X̃_{j} for j < j_max, but ALSO possibly other X̃_{j_max} (different occurrence).

Wait — different occurrence of X̃_{j_max} is the SAME random variable X̃_{j_max}. So Z_left or Z_right can contain X̃_{j_max} (the same one as at position 4). Then Z_left or Z_right is NOT in the algebra at steps ≤ j_max − 1; it's in the algebra including X̃_{j_max} too.

In our case (n=4, j_1 < j_2 alternating): Z_left = X̃_{j_1} X̃_{j_2} X̃_{j_1} where X̃_{j_2} = the same X̃_{j_2} as at position 4. So Z_left depends on the within-pair split at step j_2 as well. The peak rule at position 4 would need:

   `E_{B_∞}(Z_left · X̃_{j_2}) = E_{B_∞}(Z_left · 0) = 0`

But the LHS is `E_{B_∞}(Z_left · X̃_{j_2})`. We can't factor this as `E[Z_left] · E[X̃_{j_2}]` because Z_left depends on X̃_{j_2}.

**This is the genuine n=4 obstruction surfaced in `H1_PRIME_LOW_ORDER_CHECKS.md §4`.** Argument A §1.4 ASSUMED that "Z_left and Z_right are in A_{<j_i}" — meaning at steps < j_i. But when the SAME step j_i appears multiple times in the word, this assumption fails.

#### Mode-E gap (STRUCTURAL-1): non-adjacent repeats of the SAME peak step

The structural argument §1.4 holds when each j_i in the index sequence is distinct, or when repeats are adjacent. For NON-ADJACENT repeats of the same step (e.g. position 2 and position 4 both at j_2 > j_1), the "Z_left ∈ algebra at steps < j_max" claim fails.

**HS 2014 Defn 2.2 verbatim** handles this case by the peak rule: at position 4 (endpoint peak), substitute X̃_{j_2} → ϕ(X̃_{j_2}) = 0. The substituted RHS = 0. For the LHS to also be 0 (so the equality holds), there must be an algebraic mechanism.

In the standard scalar monotone-independence theory (Muraki / Hasebe monograph), this case is handled by the peak rule recursion: after substituting at position 4, the resulting product `Z_left · 0` reduces to 0. The LHS evaluation matches by sliding the B-value 0 through Z_left. But this is exactly the substitution, NOT an independent LHS evaluation.

The equality `ϕ(Z_left · X̃_{j_max}) = ϕ(Z_left · 0) = 0` is therefore a STATEMENT one wants to prove, not derive. It's the content of Defn 2.2 itself.

**To verify Defn 2.2 for our setting, we need to verify directly that `ϕ(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}) = 0`** under marginal centering and ϕ = E_{B_∞}. This is a direct computation.

**Direct computation at n = 4 alternating (sketched).**

Write Off_{j} = Σ_{v ≠ v'} 2^{−v} 2^{−v'} χ_{j} σ_{−(v + v')} where χ_j has the b_{[1, j−1]}-dependent phase.

Then Off_{j_1} Off_{j_2} Off_{j_1} Off_{j_2} = (Σ_{v_1 ≠ v_1'} ...) (Σ_{v_2 ≠ v_2'} ...) (Σ_{v_3 ≠ v_3'} ...) (Σ_{v_4 ≠ v_4'} ...)

where (v_1, v_1') and (v_3, v_3') are the within-pair pairs at the two occurrences of step j_1 (same random variable, so v_3 = v_1, v_3' = v_1' WHEN COMPUTING ϕ — they're the same random variable, even though the operator product treats them as separate symbols in the algebra).

Wait — that's not right either. The operator X̃_{j_1} is an OPERATOR, not a random variable. The "within-pair pair" (v, v') is summed inside the operator definition: Off_{j_1}(f)(ξ) = Σ_{v ≠ v'} weight · phase · shift · f. The randomness in the within-pair pair is part of the iid construction at the algebraic level — but X̃_{j_1} is a fixed operator (modulo a fixed B_marginal-dependent phase) once you've specified the operator algebra.

So the two X̃_{j_1}'s in the word ARE the same operator, and the calculation of ϕ(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}) is computing the B-valued expectation of a fixed operator product. There's no "v_1 = v_3, v_2 = v_4" ambiguity; the operators are nailed down.

The expectation ϕ averages over the iid construction of independent copies (HS 2014 Defn 2.3 / dot operation), which IS where Defn 2.2 is invoked. Reading HS 2014 more carefully: Defn 2.2 is the AXIOM that defines monotone independence, and the iid copies in Defn 2.3 are CONSTRUCTED to satisfy Defn 2.2. So computing ϕ on iid copies uses Defn 2.2 as a defining property.

For our Syracuse setting, the "iid copies" structure is built into the operator algebra at level n (each step j has its own (v_{2j−1}, v_{2j}) pair, independent across j). The Defn 2.2 property would be inherited from this iid construction IF the algebra A_j is monotone independent — which is the H1' question.

**The H1' question is therefore: does the natural Syracuse iid-pair structure produce a monotone-independent family?** Argument A says yes, modulo the structural fact §1.4 which I now realize has the **STRUCTURAL-1 gap** (non-adjacent same-peak repeats).

### Gap 3.5 — the STRUCTURAL-1 gap in detail

For the index sequence (j_1, j_2, j_1, j_2) with j_1 < j_2, n = 4:

- Defn 2.2 at peak position 4 says: `ϕ(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}) = ϕ(X̃_{j_1} X̃_{j_2} X̃_{j_1} · ϕ(X̃_{j_2}))`. With ϕ = E_{B_∞} and marginal centering, ϕ(X̃_{j_2}) = 0, so RHS = 0.

- Defn 2.2 at peak position 2 says: `ϕ(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}) = ϕ(X̃_{j_1} · ϕ(X̃_{j_2}) · X̃_{j_1} X̃_{j_2})`. ϕ(X̃_{j_2}) = 0, so RHS = 0.

- BOTH substitutions give RHS = 0, so for Defn 2.2 to hold, LHS must = 0.

**Is LHS = E_{B_∞}(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}) actually zero?** This is the structural question.

Heuristically: the product X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2} contains the "iid copies" of X̃_{j_1} (at positions 1 and 3) and of X̃_{j_2} (at positions 2 and 4). In the standard monotone-independence iid construction (HS 2014 Defn 2.3), independent copies X̃_{j_1}^{(1)}, X̃_{j_1}^{(2)}, ... are introduced as fresh random variables, so the two X̃_{j_1}'s at positions 1 and 3 would be different copies — not the same operator.

**In our Syracuse setting, the two X̃_{j_1}'s at positions 1 and 3 in a word ARE the same operator** (we have ONE X̃_{j_1}, not multiple copies). So the structure here is different from HS 2014 Defn 2.3's iid-copy construction.

This is a categorical distinction:

- **HS 2014 framework:** monotone independence of a FAMILY (A_λ)_{λ ∈ Λ} of subalgebras at distinct λ-indices. Each A_λ is a SUBALGEBRA generated (over B) by some X_λ ∈ A. Different λ's give different subalgebras of the SAME A.
- **Our Syracuse setting:** the family is (A_j)_{j ≥ 1} = (B_marginal⟨X̃_j⟩_0)_{j ≥ 1}, each A_j is a subalgebra of A generated by X̃_j ∈ A. ✓ Same as HS 2014.
- **The question of "repeated indices in a word":** Defn 2.2 explicitly allows index sequences with repeats. The peak rule applies to peaks in the (j_1, …, j_n) sequence, which may include repeated j's. The algebra elements X_i ∈ A_{j_i} at repeated indices may or may not be the same algebra element — the definition allows any X_i ∈ A_{j_i}.

In particular, for the word X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}, we can take X_1 = X̃_{j_1}, X_2 = X̃_{j_2}, X_3 = X̃_{j_1}, X_4 = X̃_{j_2} — the same algebra elements at repeated positions. This is a valid choice in Defn 2.2.

**So the question is whether ϕ(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}) = 0 actually holds for our Syracuse operators.**

#### A direct moment calculation

`E_{B_∞}(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2})` under marginal centering, single ϕ = E_{B_∞}:

- B_∞ fixes all accumulators b_{[1, k]} for all k. The remaining randomness is the within-pair splits (v_{2k − 1}, v_{2k}) at each k.
- The splits at different k are independent (iid Geom(2) at each k).
- X̃_{j_1} contains the within-pair split at step j_1: (v_{2j_1 − 1}, v_{2j_1}).
- X̃_{j_2} contains the within-pair split at step j_2: (v_{2j_2 − 1}, v_{2j_2}).
- These two pairs are independent.

The product X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2} is an operator product. When B_∞ averages it, the averaging is over the (v_{2j_1 − 1}, v_{2j_1}) AND (v_{2j_2 − 1}, v_{2j_2}) pairs INDEPENDENTLY (since B_∞ has fixed all b_{[1, k]}'s, the remaining randomness is the within-pair splits, independent across steps).

Compute by Fubini over the two independent random pairs:

   `E_{B_∞}(...) = E_{(v_{j_1})}[ E_{(v_{j_2})}[X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}] ]`

The inner expectation is over (v_{2j_2 − 1}, v_{2j_2}) given (v_{2j_1 − 1}, v_{2j_1}) fixed and B_∞ fixed.

`E_{(v_{j_2})}[X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}] = X̃_{j_1} · E_{(v_{j_2})}[X̃_{j_2} X̃_{j_1} X̃_{j_2}]`  — pulling X̃_{j_1} out (it doesn't depend on v_{j_2}, treating B_∞ + within-pair-j_1 as fixed environment).

Inside: `E_{(v_{j_2})}[X̃_{j_2} · X̃_{j_1} · X̃_{j_2}]`. The X̃_{j_1} in the middle is also fixed (it's a function of v_{j_1} which is fixed by the outer conditioning). The expectation is over two X̃_{j_2}'s of the same pair (v_{2j_2 − 1}, v_{2j_2}) sandwiching the fixed X̃_{j_1}.

This is generically NON-ZERO: `E_{(v_{j_2})}[X̃_{j_2} · X̃_{j_1} · X̃_{j_2}] = E[X̃_{j_2}² · X̃_{j_1}']` modulo shift composition, where X̃_{j_1}' is X̃_{j_1} possibly modified by the shifts that X̃_{j_2}'s carry. The squared X̃_{j_2}² has positive expectation (it's the diagonal κ_2 contribution), and multiplying by X̃_{j_1}' gives a non-trivial product that doesn't average to zero.

**So `E_{B_∞}(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}) ≠ 0` generically.** ✗

**This violates HS 2014 Defn 2.2 at n = 4, alternating (j_1, j_2, j_1, j_2) with j_1 < j_2.**

### Diagnosis

The peak rule of Defn 2.2 says ϕ(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}) should equal ϕ(X̃_{j_1} · 0 · X̃_{j_1} · X̃_{j_2}) = 0. But the direct calculation gives a non-zero value. **Defn 2.2 fails.**

The structural reason: the "peak rule" substitutes the algebra element X̃_{j_2} (at peak position) by ϕ(X̃_{j_2}) = 0, which is a SCALAR (B-valued) substitution. But the second X̃_{j_2} at position 4 IS A DIFFERENT INSTANCE of the same algebra element — its "randomness" is shared with the first X̃_{j_2}.

**This is the fundamental tension** between monotone independence (where "iid copies" of A_λ-elements at different positions in a word are independent random variables) and the Syracuse setting (where each step contributes ONE operator X̃_j, used as a single fixed element of A_j everywhere it appears).

In the standard monotone-independence iid framework (HS 2014 Defn 2.3 dot operation), each X̃_j^{(k)} at the k-th occurrence is an independent copy. With independent copies, the second X̃_{j_2}^{(2)} would be independent of the first X̃_{j_2}^{(1)}, and the moment ϕ(X̃_{j_1} X̃_{j_2}^{(1)} X̃_{j_1} X̃_{j_2}^{(2)}) WOULD factor through the centering of X̃_{j_2}^{(2)} to give 0. ✓

But in the Syracuse setting, **there's only ONE X̃_{j_2}** (no independent copies). The peak rule fails because the two instances of X̃_{j_2} share randomness.

### Verdict on H1'

**H1' (HS 2014 Defn 2.2 verbatim, single ϕ = E_{B_∞}, marginal centering, A_j = B_marginal⟨X̃_j⟩_0):** FAILS at n = 4 alternating (j_1, j_2, j_1, j_2) with j_1 < j_2.

The mechanism of failure is structural: HS 2014 Defn 2.2 is designed for families of subalgebras whose representatives at different positions in a word are independent. The Syracuse setting has a single X̃_j per step, used wherever step j appears in a word. This violates the implicit "iid copy" assumption.

**The correct framework for Syracuse is not Defn 2.2 per se, but a NARROWER notion** — perhaps "monotone independence at distinct indices only" (where words are restricted to have all-distinct indices, n ≤ # steps).

**But this narrower notion does NOT yield HS 2014 Thm 3.4 verbatim,** because Thm 3.4 sums over ALL monotone partitions of [n], including partitions with repeated-index structure.

---

## 4. Comparison with the W1 evidence

`W1_BLIFT_VERIFICATION.md §1-5` carried out a "sanity check at n = 3" for the alternating moment, identifying the diagnostic 0.10783 with κ_3^B per HS 2014 Thm 3.4. Re-reading §3 ("Specialization to centered X̃_j"): the moment-cumulant formula at n = 3 collapses, after κ_1^B = 0 centering, to `E_B(X̃_{j_1} X̃_{j_2} X̃_{j_1}) = κ_3^B(X̃_{j_1}, X̃_{j_2}, X̃_{j_1})`. The diagnostic 0.10783 is the value of this cumulant.

This is CONSISTENT with Thm 3.4 — but Thm 3.4 is established **on the assumption of Defn 2.2**, which we've now seen FAILS at n = 4. So Thm 3.4 does not actually apply to our Syracuse X̃_j family in the way W1 used it.

**Audit caveat upgraded.** The W1 lift of Thm 3.4 is conditional on Defn 2.2, and Defn 2.2 fails. The c = 7/45 derivation depending on Thm 3.4 inherits this failure.

---

## 5. What's actually true for Syracuse

The X̃_j family satisfies **weaker** properties:

- **Distinct-index independence.** For any monomial X̃_{j_1} … X̃_{j_n} with all j_i distinct, the moment factorizes: `ϕ(X̃_{j_1} … X̃_{j_n}) = 0` (after the appropriate centering kills it).
- **Alternating peak structure at n = 3, alternating-with-repeat.** The diagnostic 0.10783 lives in a specific operator subspace identifiable via the level-graded reading.
- **Pascal-pair second-order independence.** Cross-step κ_2^B(X̃_{j_1}, X̃_{j_2}) = 0 (W2 finding).

These match the SECOND-ORDER moments and DISTINCT-INDEX higher moments of monotone independence. They do NOT match the FULL monotone-independence structure at non-adjacent repeated indices.

The honest characterization: the X̃_j family is **partially monotone-like** — it satisfies the second-order and distinct-index pieces of monotone independence, but FAILS the non-adjacent-repeated-index peak rule.

The framework that ACTUALLY governs the c = 7/45 derivation is therefore neither standard free nor standard monotone, but a SPECIFIC weaker structure for which:
- The leading 7/45 coefficient comes from the κ_1^B all-singletons partition (this IS monotone-like and doesn't invoke peak rules at repeated indices).
- The subdominant κ_2^B comes from W2's diagonal calculation (distinct-index cross-cumulants vanish, ✓).
- Higher-order corrections at NON-ADJACENT repeated indices break monotone independence and would need a separate framework.

---

## 6. Implications and what's salvageable

The c = 7/45 LEADING coefficient derivation uses:

(a) HS 2014 Thm 3.4 at n = 1 (trivially): ϕ(X̃_j) = κ_1^B(X̃_j) — actually, at the LEADING all-singletons partition, the coefficient is (κ_1^B(Off_j))^n / n!, and the leading 7/45 is the κ_1^B(Off_j) projection onto the (1,4)-eigenvector. **The all-singletons partition is the n-tuple of distinct positions, each in its own block — the index sequence can have repeats but the moment-cumulant evaluation at this partition involves only κ_1^B's, which are evaluated on single elements.** This part does NOT require Defn 2.2's peak rule at repeated indices.

(b) HS 2014 Thm 3.4 at higher partitions: these DO involve cumulants K_n at n ≥ 2, which can be evaluated on repeated indices. The 2-block partition `{{1, 3}, {2}}` etc. evaluate K_2(X̃_{j_1}, X̃_{j_1}) · K_1(X̃_{j_2}) etc., which are diagonal κ_2 at j_1, involving a single algebra and not invoking the cross-index peak rule.

(c) **The cross-index peak rule at NON-adjacent repeated indices** is the failure point. This appears at n ≥ 3 in some Thm 3.4 expansions if the moment has the alternating non-adjacent-repeat pattern.

**Possibly salvageable:** the leading 7/45 coefficient (κ_1^B all-singletons) and the subdominant amplitude via diagonal κ_2 may survive even if the cross-index peak rule fails, because they don't directly invoke that peak rule.

**The audit needs to re-examine** which moments in the c = 7/45 derivation actually require Defn 2.2 at non-adjacent repeats. This is a focused audit, not a full re-derivation.

---

## 7. Files

- Verbatim setup: `H1_PRIME_VERBATIM.md`
- Low-order checks: `H1_PRIME_LOW_ORDER_CHECKS.md`
- Disposition (verdict): `H1_PRIME_DISPOSITION.md`
- W1 sanity check at n = 3: `W1_BLIFT_VERIFICATION.md`
- Diagnostic JSON: `experiments_output/monotone_diagnostic_n3.json`
- HS 2014 PDF (verbatim source): `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_saigo_2014_operator_valued_monotone.pdf`
- W2 cross-step κ_2: `W2_KAPPA2_CALC.md`
- AMALG moment calculation: `AMALG_FREENESS_MOMENT_CALCULATION.md`
- Track A integration: `TRACK_A_INTEGRATION.md §1.4`
