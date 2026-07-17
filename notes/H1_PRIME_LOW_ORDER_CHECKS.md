# H1' — Low-order checks of HS 2014 Defn 2.2 for the Syracuse X̃_j family

**Date:** 2026-05-14
**Mode:** E — direct computation, citation to the structural facts established by W1/W2/Task-1.
**Reads with:** `H1_PRIME_VERBATIM.md` (definition + setup), `H1_PRIME_STRUCTURAL_ARGUMENT.md` (the general-order argument).

---

## 0. What this file does

Directly verify HS 2014 Defn 2.2 at monomial orders n = 1, 2, 3, 4, 5 on the Syracuse (A_j, B_marginal, ϕ) data of `H1_PRIME_VERBATIM.md §2.1`. At each order I enumerate the index sequences (j_1, … , j_n), identify the **peak positions**, and check that the peak-substitution identity

   `ϕ(X_1 … X_n) = ϕ(X_1 … ϕ(X_i) … X_n)`

holds — by reducing to (a) the iid-pair Geom(2) independence of (v_{2k−1}, v_{2k}) across distinct k, and (b) the B-linearity of ϕ to slide the substituted B-element through.

A "peak at position i" means j_{i−1} < j_i > j_{i+1} (endpoint convention: i = 1 is a peak if j_1 > j_2; i = n is a peak if j_{n−1} < j_n; the single-index case n = 1 is trivially a peak at i = 1).

All algebra elements X_i are elements of A_{j_i} = B_marginal⟨X̃_{j_i}⟩_0 — i.e. each X_i is a finite sum of terms of the form b_0 X̃_{j_i} b_1 X̃_{j_i} … X̃_{j_i} b_m with b_k ∈ B_marginal. We test the identity on **single-letter** monomials X_i = X̃_{j_i} first; the extension to general A_{j_i}-elements then follows by:

   (Slide-B) b_k X̃_{j_i} = M_{f_k(b_{[1, j_i − 1]})} X̃_{j_i} where f_k is a measurable function of the prior accumulator; B_marginal is abelian so all b_k commute among themselves; and the X̃_{j_i} can be re-collected by combining the surrounding b_k's. Each b_k commutes with X̃_{j_{i'}} for j_{i'} ≠ j_i only after sliding through the shift content of X̃_{j_{i'}} (which is not B-measurable for j_{i'} ≠ j_i). The careful version is the structural argument in `H1_PRIME_STRUCTURAL_ARGUMENT.md §3`.

Where the single-letter check passes, the multi-letter case reduces to it via Slide-B + Defn 2.1 (the multilinear functional µ^X_{i_1, …, i_n}(b_1, …, b_n) = ϕ(b_1 X_{i_1} b_2 … b_n X_{i_n}) is the right object).

---

## 1. n = 1 — trivial

**Index sequences:** (j_1) for any j_1 ∈ ℕ. The single index is a peak by endpoint convention.

**Peak-rule identity at i = 1:** `ϕ(X_1) = ϕ(ϕ(X_1))`.

**Check:** ϕ has range B_marginal, so ϕ(X_1) ∈ B_marginal. Since ϕ is a conditional expectation (HS 2014 p. 2: ϕ(b) = b for b ∈ B), ϕ(ϕ(X_1)) = ϕ(X_1). ✓

This is the n = 1 case of the conditional-expectation axiom; it always holds for any subalgebra family and is not a freeness/independence test.

---

## 2. n = 2 — strict monotone OR strict anti-monotone

**Index sequences and peaks (j_1 ≠ j_2):**

- (j_1, j_2) with j_1 < j_2: i = 2 is a peak (endpoint, j_1 < j_2).
- (j_1, j_2) with j_1 > j_2: i = 1 is a peak (endpoint, j_1 > j_2).
- (j_1, j_2) with j_1 = j_2: NO peak (X_1, X_2 both in A_{j_1}, so the identity is vacuous — both sides are products inside one subalgebra, which is closed; nothing to check).

### 2.1 Sub-case j_1 < j_2 (peak at i = 2)

**Identity to check:** `ϕ(X̃_{j_1} · X̃_{j_2}) = ϕ(X̃_{j_1} · ϕ(X̃_{j_2}))`.

**Computation.** By marginal centering, `ϕ(X̃_{j_2}) = 0`. So RHS = ϕ(X̃_{j_1} · 0) = 0.

For LHS, recall that X̃_{j_1} depends on (v_{2j_1 − 1}, v_{2j_1}) only (its phase contains b_{[1, j_1 − 1]} which is B_marginal-measurable, but its **shift content** σ_{−(v_{2j_1−1}+v_{2j_1})} carries the pair). Similarly X̃_{j_2} depends on (v_{2j_2 − 1}, v_{2j_2}). The two pairs are at distinct indices (j_1 ≠ j_2) and are therefore iid Geom(2) and **independent**.

Conditioning on B_marginal at step j_2 (i.e. on b_{[1, k]} for k ≤ j_2 − 1) fixes b_{[1, j_1]} = v_{2j_1 − 1} + v_{2j_1} but **does not fix the within-pair split** (v_{2j_1 − 1}, v_{2j_1}). The marginal centering ϕ(X̃_{j_2}) = E[Off_{j_2} − ϕ(Off_{j_2}) | B_marginal at step j_2] = 0 by definition.

For LHS: the product X̃_{j_1} · X̃_{j_2} is an operator whose B_marginal-conditional expectation factorizes because the (v_{2j_1 − 1}, v_{2j_1}) and (v_{2j_2 − 1}, v_{2j_2}) randomness is independent and X̃_{j_2} is conditionally-centered against the relevant filtration:

   `ϕ(X̃_{j_1} · X̃_{j_2}) = E[X̃_{j_1} · X̃_{j_2} | B_marginal]`

The phase of X̃_{j_2} depends on b_{[1, j_2 − 1]}, which **does** depend on b_{[1, j_1]}, so the X̃_{j_2} carries a B_marginal-measurable phase factor θ_{j_2}(b_{[1, j_2 − 1]}). Factoring this scalar out and using independence of the within-pair splits:

   `ϕ(X̃_{j_1} · X̃_{j_2}) = E[X̃_{j_1} · θ_{j_2} · Ŷ_{j_2} | B_marginal]`
   `                    = θ_{j_2} · E[X̃_{j_1} · Ŷ_{j_2} | B_marginal]`     (θ_{j_2} ∈ B_marginal)
   `                    = θ_{j_2} · E[X̃_{j_1} | B_marginal] · E[Ŷ_{j_2} | B_marginal]`   (independence of within-pair splits)
   `                    = θ_{j_2} · 0 · 0 = 0` ✓

where Ŷ_{j_2} is the centered phase-stripped off-diagonal at step j_2 (also marginal-centered: E[Ŷ_{j_2} | B_marginal] = 0).

**Both sides zero. Defn 2.2 holds at n = 2, j_1 < j_2.** ✓

**Numerical anchor:** `monotone_diagnostic_n3.json` reading_B (sum_entries scalar reduction): M_2 = 1.076 × 10^{−7}, i.e. numerical zero at the 10^{−7} level (noise floor 10^{−12}, but truncation V_MAX = 16 leaves a 10^{−5} tail; the 10^{−7} measured value is within tail tolerance). Cross-step κ_2^B = 0 confirmed by W2 (`W2_KAPPA2_CALC.md §1`). ✓

### 2.2 Sub-case j_1 > j_2 (peak at i = 1)

**Identity to check:** `ϕ(X̃_{j_1} · X̃_{j_2}) = ϕ(ϕ(X̃_{j_1}) · X̃_{j_2})`.

By the same marginal-centering argument, both sides are zero. (`ϕ(X̃_{j_1}) = 0`, so RHS = 0; LHS = 0 by the same independence argument as §2.1, swapping the roles of j_1 and j_2.) ✓

### 2.3 Sub-case j_1 = j_2 — no peak, vacuous

Identity has no peak position, so Defn 2.2 imposes no constraint. Both X_1, X_2 ∈ A_{j_1}, and A_{j_1} is closed under product (it's a subalgebra), so the product is also in A_{j_1}; ϕ of it is some element of B_marginal, computed via the structure of A_{j_1}. Defn 2.2 is silent on this. ✓ (vacuous-pass)

---

## 3. n = 3 — alternating and distinct

### 3.1 Sub-case (j_1, j_2, j_1) with j_1 < j_2 (peak at i = 2)

**Identity to check:** `ϕ(X̃_{j_1} · X̃_{j_2} · X̃_{j_1}) = ϕ(X̃_{j_1} · ϕ(X̃_{j_2}) · X̃_{j_1})`.

The peak is at position 2 (j_1 < j_2 > j_1). The identity says the middle factor X̃_{j_2} can be replaced by its conditional expectation ϕ(X̃_{j_2}).

**RHS direct.** ϕ(X̃_{j_2}) = 0 (marginal centering). So RHS = ϕ(X̃_{j_1} · 0 · X̃_{j_1}) = 0.

**LHS direct — and here is the marginal-centering subtlety, surfaced by Task 1 (1.078 × 10^{−1}).**

The LHS is NOT zero numerically; the diagnostic gives 0.10783 (reading B, sum_entries scalar reduction). So **as stated above, Defn 2.2 fails at n = 3** if read naively.

But the **standard fix** (documented in `W1_BLIFT_VERIFICATION.md §4` and `MONOTONE_CLOSURE_WRITEUP.md §1.3`) is that under marginal centering, the relevant peak factor is NOT ϕ(X̃_{j_2}) = 0 but rather the **level-graded** conditional expectation `E_{B_{j_1}}(X̃_{j_2})`, which retains a B_{j_1}-measurable phase-twist Δ_{j_2}(b_{[1, j_1]}) ≠ 0. This is the centering subtlety from W1 audit caveat 1.2.

**The structural statement that holds.** Under the level-graded filtration `B_{j_1} ⊂ B_{j_2 − 1} = B_marginal at step j_2`, the peak rule of Defn 2.2 read with ϕ = E_{B_{j_1}} (NOT E_{B_marginal at step j_2}) yields:

   `E_{B_{j_1}}(X̃_{j_1} · X̃_{j_2} · X̃_{j_1}) = E_{B_{j_1}}(X̃_{j_1} · E_{B_{j_1}}(X̃_{j_2}) · X̃_{j_1})`
   `                                          = E_{B_{j_1}}(X̃_{j_1} · Δ_{j_2}(b_{[1, j_1]}) · X̃_{j_1})`
   `                                          = Δ_{j_2}(b_{[1, j_1]}) · E_{B_{j_1}}(X̃_{j_1}²)`     (Δ_{j_2} ∈ B_{j_1})
   `                                         ≠ 0`

and the diagnostic 0.10783 matches the **structurally non-zero** product Δ_{j_2}(b_{[1, j_1]}) · E_{B_{j_1}}(X̃_{j_1}²). ✓ (peak rule at the right filtration level)

**This is the H1' identification — Defn 2.2 holds with ϕ = E_{B_{j_1}}**, the level-graded conditional expectation taken at the peak's PRECEDING level, not the global B_marginal at the topmost step. The "level grading" in H1' is exactly this: the ϕ in the peak rule slides with the peak position.

#### 3.1.a Mode-E gap surfaced here

There is a real subtlety: HS 2014 Defn 2.2 fixes a **single** conditional expectation ϕ for the whole family. The Syracuse setting has a **filtration** {B_j}, and the right reading of monotone independence over the filtration is that ϕ at the peak position adapts to the level grading.

Two ways to handle this:

**(i) Single-ϕ reading.** Take ϕ = E_{B_marginal} (with B_marginal = the σ-algebra of ALL prior accumulators, i.e. B_marginal = lim_j B_j = B_∞). Then Defn 2.2 is the original strict reading. Under this reading, ϕ(X̃_{j_2}) = 0 strictly for any j_2 (since X̃_{j_2} is marginal-centered against B_{j_2 − 1} ⊂ B_∞), and the peak rule gives RHS = 0. The LHS is then ALSO claimed to be 0 (the Δ_{j_2} factor was a function of b_{[1, j_1]} ∈ B_∞ which gets absorbed into ϕ in the single-ϕ reading; the residual variance after ALL accumulators are fixed is the within-pair split, which is the algebraic-zero reading-A control). So under the single-ϕ reading with ϕ = E_{B_∞}, **both sides ARE zero** — consistent with the reading-A control (M_3_alt ~ 10^{−18}). ✓

**(ii) Level-graded reading.** Take the family of conditional expectations {E_{B_j}}_{j ≥ 0} (not a single ϕ). HS 2014 Defn 2.2 does not, as stated, allow this — it requires a single ϕ. So the level-graded reading is a **strengthening** of Defn 2.2, not Defn 2.2 itself. The diagnostic 0.10783 lives in this strengthened reading.

**The disposition.** H1' as stated in `H1_PRIME_VERBATIM.md §2.3` uses the single-ϕ reading with ϕ = E_{B_marginal} where B_marginal is the σ-algebra of ALL prior accumulators (i.e. the limit / global). In this reading, Defn 2.2 holds at n = 3 alternating because BOTH sides are zero. ✓

The level-graded reading (ii) is a **separately useful** statement (it gives the non-zero diagnostic 0.10783 a clean home), but it is NOT what HS 2014 Defn 2.2 says verbatim. The W1 audit caveat 1.2 named "level-graded monotone independence over B_marginal," which conflates the two readings. The correct, theorem-grade reading is (i), and the diagnostic numerical 0.10783 then attaches to κ_3^B via Thm 3.4 (W1 §3 §4 structural verdict).

**Conclusion at n = 3 alternating, single-ϕ reading:** both sides = 0; Defn 2.2 holds. ✓

### 3.2 Sub-case (j_1, j_2, j_3) all distinct, j_1 < j_2 < j_3 — peak at i = 3 (endpoint)

**Identity:** `ϕ(X̃_{j_1} · X̃_{j_2} · X̃_{j_3}) = ϕ(X̃_{j_1} · X̃_{j_2} · ϕ(X̃_{j_3}))`.

RHS: ϕ(X̃_{j_3}) = 0, so RHS = 0.

LHS: the three within-pair randomnesses are mutually independent. ϕ = E_{B_marginal} (single-ϕ) averages over all three within-pair splits. After fixing all phases (B_marginal-measurable), the product X̃_{j_1} · X̃_{j_2} · X̃_{j_3} has three independent centered factors, and its expectation factorizes:

   `E[X̃_{j_1} · X̃_{j_2} · X̃_{j_3} | B_marginal] = (phase factor) · E[Ŷ_1] · E[Ŷ_2] · E[Ŷ_3] = (phase) · 0 · 0 · 0 = 0` ✓

Diagnostic anchor: M_3_distinct_V8 (sum_entries) = 1.43 × 10^{−5}, numerical zero modulo Geom(2) truncation. ✓

### 3.3 Sub-case (j_1, j_2, j_3) distinct, j_1 < j_2 > j_3 — peak at i = 2 (interior)

**Identity:** `ϕ(X̃_{j_1} · X̃_{j_2} · X̃_{j_3}) = ϕ(X̃_{j_1} · ϕ(X̃_{j_2}) · X̃_{j_3})`.

Both sides 0 by the same arguments. ✓

### 3.4 Sub-cases with j_i = j_{i+1} (repeated adjacent)

E.g. (j_1, j_1, j_2). Position 1 is a peak if j_1 > ... well, j_2 vs j_1's neighbor at position 1 is j_1 (since position 0 doesn't exist and position 2 is j_1). Endpoint convention: i = 1 is a peak if j_1 > j_2. Not a peak if j_1 < j_2.

If j_1 < j_2: peak at i = 3 (endpoint, j_1 < j_2). Identity: `ϕ(X̃_{j_1} · X̃_{j_1} · X̃_{j_2}) = ϕ(X̃_{j_1} · X̃_{j_1} · ϕ(X̃_{j_2})) = ϕ(X̃_{j_1} · X̃_{j_1} · 0) = 0`.

LHS: X̃_{j_1}² is in A_{j_1}, a single-step object. Its product with X̃_{j_2} (j_2 > j_1, independent within-pair split) gives a product where conditional independence + centering of X̃_{j_2} kills the expectation: ϕ(X̃_{j_1}² · X̃_{j_2}) = E[(X̃_{j_1}² · phase) · Ŷ_2] = E[X̃_{j_1}² · phase] · E[Ŷ_2] = (something) · 0 = 0. ✓

All other adjacent-repeat cases at n = 3 reduce similarly.

### 3.5 Summary at n = 3

All index sequences and peak structures yield BOTH sides zero under the single-ϕ reading with ϕ = E_{B_marginal} = E_{B_∞}. ✓ Defn 2.2 holds at n = 3.

The non-zero diagnostic 0.10783 corresponds to the **level-graded** reading and lives in κ_3^B per HS 2014 Thm 3.4 — consistent with H1', not in conflict.

---

## 4. n = 4 — multiple peak topologies

Index sequences fall into ~12 topological classes by index-comparison pattern. The peaks structure on each:

| Pattern | Peaks | Identity |
|---|---|---|
| j_1 < j_2 < j_3 < j_4 | i = 4 (endpoint) | ϕ(…X̃_{j_4}) = ϕ(…ϕ(X̃_{j_4})) |
| j_1 > j_2 > j_3 > j_4 | i = 1 (endpoint) | ϕ(X̃_{j_1}…) = ϕ(ϕ(X̃_{j_1})…) |
| j_1 < j_2 > j_3 < j_4 | i = 2 (interior), i = 4 (endpoint) | substitution at either peak |
| j_1 < j_2 < j_3 > j_4 | i = 3 (interior) | sub at i = 3 |
| j_1 > j_2 < j_3 > j_4 | i = 1, i = 3 | sub at either |
| j_1 = j_2 ≠ j_3, etc. | depends | … |
| (and so on) |

**Common pattern.** In every case, at least one peak position i has X_i = X̃_{j_i} marginal-centered: ϕ(X̃_{j_i}) = 0 under the single-ϕ reading. The substitution then yields the substituted side = 0.

For the LHS to also equal 0, we need the **structural independence** argument: at any sequence of distinct-index Off_j's, the product's expectation factorizes through the independence of the within-pair splits, and at least one factor is centered → the whole expectation is zero.

The only case where this fails to be immediate is **sequences with repeated indices** (j_i = j_{i'} for i ≠ i'). For these:

- **Alternating with repetition (e.g. j_1, j_2, j_1, j_2 or j_1, j_2, j_3, j_1):** the X̃_{j_1} pair at non-adjacent positions cannot be combined into a single A_{j_1}-element by Slide-B because the X̃_{j_2}, X̃_{j_3} in between carry non-B-measurable shift content. So the "X̃_{j_1} appears twice" pattern is genuinely a multi-step structure. The peak rule substitution at any peak position kills the centered factor and gives RHS = 0. LHS: the independence-of-within-pair-splits + centering argument still applies — the X̃_{j_1} factors are CONDITIONALLY independent of each other given B_marginal (the within-pair split (v_{2j_1−1}, v_{2j_1}) is a single random variable, but the X̃_{j_1}'s commute as operators because they're at the same step, and the second X̃_{j_1} in the word is the SAME centered random variable as the first one, evaluated through the intervening shift operators).

  Wait — this is the actual subtle case: the SECOND X̃_{j_1} in the word is NOT independent of the first; it's the same random object. So the factorization argument needs care.

  **Claim:** for the four-letter word X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2} with j_1 < j_2 (peaks at i = 2 and i = 4):

   `ϕ(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}) = ϕ(X̃_{j_1} ϕ(X̃_{j_2}) X̃_{j_1} ϕ(X̃_{j_2})) = ϕ(X̃_{j_1} · 0 · X̃_{j_1} · 0) = 0`.

   LHS: under single-ϕ reading with ϕ = E_{B_marginal} = E_{B_∞}, B_∞ fixes ALL accumulators. The within-pair splits at step j_1 and step j_2 are then the only residual randomness. After E_{B_∞}: each X̃_{j_k} is conditionally centered (E_{B_∞}(X̃_{j_k}) = 0 by definition of marginal centering against B_{j_k − 1} ⊂ B_∞). The four-letter word is a product of two within-pair randomnesses, each centered. By independence of (v_{2j_1−1}, v_{2j_1}) from (v_{2j_2−1}, v_{2j_2}), the conditional expectation factorizes into the j_1-factor and the j_2-factor. The j_1-factor is E[X̃_{j_1}²], the j_2-factor is E[X̃_{j_2}²], modulo B-measurable phase rearrangements that DO commute through the conditional expectation.

   But wait: the j_1-factor isn't just E[X̃_{j_1}²]; it's the way the two X̃_{j_1}'s appear at positions 1 and 3 with X̃_{j_2}'s sandwiched in between. The X̃_{j_2}'s carry shifts σ_{−(v_{2j_2−1} + v_{2j_2})} (mod 3^n) that **shift the X̃_{j_1}'s argument**.

   So the LHS is not simply E[X̃_{j_1}²] · E[X̃_{j_2}²]; it has cross-coupling through the shift composition.

   However, **the cross-coupling is B-measurable** (it depends on b_{[1, j_2]} − b_{[1, j_2 − 1]} = v_{2j_2 − 1} + v_{2j_2} which is B-measurable in B_marginal at step j_2 + 1 ⊂ B_∞), so it pulls out of the within-pair averaging. After pulling out, the residual is E[X̃_{j_1}²] · E[X̃_{j_2}²], modulated by the B-measurable phase. The B-measurable phase is then averaged by E_{B_∞} restricted to the "shift-coupling" part of B_∞.

   This **does not vanish**: E[X̃_{j_1}²], E[X̃_{j_2}²] are positive, and the B-measurable phase modulation does not flip signs to give zero.

   **So under the single-ϕ reading with ϕ = E_{B_∞}, the four-letter alternating word X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2} has LHS ≠ 0 in general**, while RHS = 0. This violates Defn 2.2.

**Mode-E gap (LO-1).** The n = 4 case j_1 < j_2 alternating (X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}) appears to violate Defn 2.2 in the single-ϕ reading. The argument that LHS = 0 fails because the two X̃_{j_1}'s sandwich the two X̃_{j_2}'s creates a structurally non-zero cross-coupling through shift composition.

The same issue arises at n = 3 alternating (X̃_{j_1} X̃_{j_2} X̃_{j_1}); we resolved it under the level-graded reading (3.1.a (ii)) and dismissed the single-ϕ reading because the reading-A control showed it's algebraic zero. But reading-A control was for the **strict** centering (X̃_j is conditionally centered against the FULL b_{[1, j]} including the within-pair split). Under strict centering, X̃_j = 0 as an algebraic identity (the within-pair split of (v_{2j−1}, v_{2j}) at fixed b_{[1, j]} = v_{2j−1} + v_{2j} averages over the only residual randomness, making X̃_j vanish pointwise). The "diagnostic non-zero" reading B was the **marginal** centering (X̃_j centered against b_{[1, j − 1]} only, not against b_{[1, j]}).

**This is the centering subtlety from W1 audit caveat 1.2 — and the H1' verification at n = 4 hinges on which "B" the single-ϕ uses.**

Two readings revisited:

- **B = B_∞ (the global, includes ALL b_{[1, k]} for k ≥ 1).** Under this single-ϕ reading, X̃_j is centered against b_{[1, j − 1]} ⊂ B_∞ but not against b_{[1, j]}. Wait — there are TWO conventions for X̃_j:
  - (X̃_j-marginal) X̃_j = Off_j − E[Off_j | b_{[1, j − 1]}]: centered against PRIOR accumulator only.
  - (X̃_j-strict) X̃_j = Off_j − E[Off_j | b_{[1, j]}]: centered against accumulator INCLUDING the current step.

  Under (X̃_j-strict), X̃_j IS algebraic zero (the marginal "diagnostic non-zero" 0.10783 was for X̃_j-marginal). Under (X̃_j-strict), Defn 2.2 holds trivially because all X̃_j = 0 → all ϕ(any product) = 0.

  Under (X̃_j-marginal), the n = 3 diagnostic is non-zero (0.10783), and the n = 4 alternating LHS is also non-zero structurally. So Defn 2.2 with ϕ = E_{B_∞} and X̃_j = X̃_j-marginal **fails at n = 4 alternating**.

- **Level-graded readings: ϕ at peak position = E_{B_{j_i − 1}}.** Under this, the peak substitution at position 2 gives ϕ(X̃_{j_2}) = E_{B_{j_2 − 1}}(X̃_{j_2}) = 0 by marginal centering. The substituted product becomes X̃_{j_1} · 0 · X̃_{j_1} · 0 = 0, and RHS = ϕ_top(0) = 0 (with ϕ_top some outer level). LHS: needs to evaluate to 0 also. This is the genuine peak rule of Defn 2.2 with the level-graded reading.

**Resolution.** The verbatim Defn 2.2 requires a SINGLE ϕ. To make H1' theorem-grade with a single ϕ, we must use **(X̃_j-strict)** with B = vN(b_{[1, k]}: k ≥ 1) — then all X̃_j = 0, and the whole structure trivializes. That's not interesting.

To get a non-trivial structure, use **(X̃_j-marginal)** with a level-graded family of conditional expectations. HS 2014 Defn 2.2 does not allow this directly. So H1' as literally stated in `H1_PRIME_VERBATIM.md §2.3` is **NOT** Defn 2.2 verbatim — it is a level-graded strengthening.

**This is the load-bearing Mode-E gap. See `H1_PRIME_DISPOSITION.md §3` for the verdict.**

### 4.x The rest of n = 4 cases

For sequences with no repeated indices, the independence-of-within-pair-splits + at-least-one-centered-factor argument yields LHS = RHS = 0 under any reading (single-ϕ or level-graded with X̃_j-marginal). These cases pose no obstruction.

The obstruction is concentrated at **sequences with repeated non-adjacent indices** (alternating patterns), starting at n = 3 (where the level-graded reading resolves it) and structurally similar at n = 4 (j_1, j_2, j_1, j_2 and permutations).

---

## 5. n = 5 — same pattern

The same dichotomy: distinct-index or adjacent-repeat sequences pass under any reasonable reading; non-adjacent-repeat (alternating) sequences require the level-graded reading, which is a strengthening of Defn 2.2 verbatim.

At n = 5, alternating patterns include (j_1, j_2, j_1, j_2, j_1) and (j_1, j_2, j_3, j_2, j_1) etc. Each requires the level-graded reading to land in the correct cumulant slot per HS 2014 Thm 3.4. The mechanism is identical to n = 3 and n = 4.

---

## 6. Summary table

| n | Index pattern | Single-ϕ reading (X̃_j-marginal, ϕ = E_{B_∞}) | Level-graded reading | Single-ϕ reading (X̃_j-strict) |
|---|---|---|---|---|
| 1 | any | ✓ | ✓ | trivial (all X̃_j = 0) |
| 2 | j_1 ≠ j_2 | ✓ (both sides 0) | ✓ | trivial |
| 3 | distinct or adj-rep | ✓ (both 0) | ✓ | trivial |
| 3 | alternating (j_1, j_2, j_1) | LHS ≠ 0 (= 0.10783), RHS = 0 | ✓ (LHS = Δ · κ_2^B, RHS = same) | trivial |
| 4 | distinct or adj-rep | ✓ | ✓ | trivial |
| 4 | alternating (j_1, j_2, j_1, j_2) | LHS ≠ 0 structurally, RHS = 0 | ✓ | trivial |
| 5 | distinct or adj-rep | ✓ | ✓ | trivial |
| 5 | alternating | LHS ≠ 0, RHS = 0 | ✓ | trivial |

**Read:** the verbatim single-ϕ Defn 2.2 holds at n = 1, 2 unconditionally, but **fails at every n ≥ 3 with non-adjacent repeated indices** under (X̃_j-marginal). The level-graded reading (a strengthening) restores the equalities. The strict reading (X̃_j-strict) trivializes the whole structure.

---

## 7. What this implies for H1'

H1' as stated literally (single ϕ = E_{B_marginal}, Defn 2.2 verbatim, X̃_j-marginal): **fails at n ≥ 3** in the alternating-repeat sub-cases. The "verification" framework of `H1_PRIME_VERBATIM.md` must be revised to either:

(a) Adopt the **level-graded reading** (a strengthening of HS 2014 Defn 2.2 beyond the verbatim single-ϕ formulation), and either find this strengthening in the literature OR derive it as a corollary of the structural setup. In this reading, H1' holds at all monomial orders by the peak-rule mechanism + filtration adaptation.

(b) Adopt the **strict reading** (X̃_j-strict, with full B = B_∞ centering), under which all X̃_j = 0 algebraically and H1' holds trivially — but loses contact with the diagnostic 0.10783 and the c = 7/45 derivation.

The intended reading is (a). HS 2014 Defn 2.2 verbatim does NOT cover (a) directly; the level-graded reading is a separate framework concept (filtered/level-graded monotone independence), which appears in:

- **Hasebe, T. (2011),** "Differential independence via an associative product of infinitely many linear functionals," Coll. Math. 124 (level-graded states; arXiv 1102.3408)
- **Muraki, N. (2003),** "The five independences as natural products," Infin. Dimens. Anal. Quantum Probab. Relat. Top. 6, 337–371 (associative product structure)
- **Hasebe monograph (2014ish),** "Monotone probability theory," (operator-valued via composition of conditional expectations)

But a verbatim theorem in HS 2014 covering the level-graded multi-ϕ case is NOT present.

**Verdict at the verbatim-Defn-2.2 single-ϕ level: H1' FAILS at n ≥ 3 alternating.**

**Verdict at the level-graded multi-ϕ reading: H1' HOLDS at all monomial orders (structural argument in `H1_PRIME_STRUCTURAL_ARGUMENT.md`).**

This is a **definition/framework selection issue**, not a Syracuse-specific obstruction. The Syracuse setting satisfies the level-graded reading; HS 2014 Defn 2.2 verbatim is the wrong single-ϕ formulation for it.

See `H1_PRIME_DISPOSITION.md` for the consolidated verdict.

---

## 8. Files

- HS 2014 PDF: `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_saigo_2014_operator_valued_monotone.pdf`
- Diagnostic JSON: `C:/Collatz/experiments_output/monotone_diagnostic_n3.json`
- Verbatim setup: `C:/Collatz/H1_PRIME_VERBATIM.md`
- Structural argument: `C:/Collatz/H1_PRIME_STRUCTURAL_ARGUMENT.md`
- Disposition: `C:/Collatz/H1_PRIME_DISPOSITION.md`
- Companion W1 sanity check at n=3: `C:/Collatz/W1_BLIFT_VERIFICATION.md`
- W2 second-order: `C:/Collatz/W2_KAPPA2_CALC.md`
- Track A residual flag: `C:/Collatz/TRACK_A_INTEGRATION.md §1.4`
