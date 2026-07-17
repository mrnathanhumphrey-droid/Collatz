# H1' — Disposition

**Date:** 2026-05-14
**Mode:** E. Self-adversarial. The verdict is NOT what the brief expected.
**Reads with:** `H1_PRIME_VERBATIM.md`, `H1_PRIME_LOW_ORDER_CHECKS.md`, `H1_PRIME_STRUCTURAL_ARGUMENT.md`

---

## 0. One-sentence verdict

**H1' FAILS at n = 4 (alternating non-adjacent-repeat index sequence (j_1, j_2, j_1, j_2) with j_1 < j_2) under the verbatim HS 2014 Defn 2.2 with single ϕ = E_{B_∞} and marginal centering.** The c = 7/45 leading derivation is therefore **NOT** unconditional via this route; the W1 lift of HS 2014 Thm 3.4 is unsound at non-adjacent repeated indices.

---

## 1. Verdict in detail

### What was checked

HS 2014 Defn 2.2 was extracted verbatim from page 3 of `hasebe_saigo_2014_operator_valued_monotone.pdf` via `pypdf` (`H1_PRIME_VERBATIM.md §1.2`). The Syracuse X̃_j family was specified precisely: A_j = B_marginal⟨X̃_j⟩_0 with ϕ = E_{B_∞}, marginal centering X̃_j = Off_j − E_{B_{j−1}}(Off_j).

Direct checks at n = 2, 3 pass under the single-ϕ reading (both sides of the peak-rule equality are zero).

Direct check at **n = 4 alternating (j_1, j_2, j_1, j_2), j_1 < j_2** FAILS:

- Peak rule substitution at position 2 or 4 gives RHS = ϕ(...0...) = 0.
- LHS = E_{B_∞}(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}) is **non-zero** by direct moment computation: after B_∞ conditioning, the residual within-pair randomness at step j_2 (which appears at positions 2 and 4) is the same random variable, NOT two independent copies. The expectation E_{(v_{j_2})}[X̃_{j_2} · X̃_{j_1} · X̃_{j_2}] is generically nonzero (~ E[X̃_{j_2}²] · X̃_{j_1}'), and integrating over the v_{j_1} randomness does not flip the sign or cancel.

The mechanism of failure is structural: HS 2014 Defn 2.2 implicitly assumes that algebra elements at different positions in a word, even if they're in the same A_λ, are treated as independent realizations under ϕ (the iid-copy convention of Defn 2.3). The Syracuse setting has a SINGLE X̃_j per step, used at every occurrence of step j in a word. The peak rule fails when the same X̃_j appears at non-adjacent positions because the two occurrences share randomness, breaking the substitutability.

### What this means for c = 7/45

The W1 lift of HS 2014 Thm 3.4 used Defn 2.2 as a structural premise (per W1 audit caveat 1.2 and `W1_BLIFT_VERIFICATION.md §1-5`). With Defn 2.2 failing at n ≥ 4 alternating, **the lift is unsound** at moments that span repeated-non-adjacent indices.

**Salvageable parts** (per `H1_PRIME_STRUCTURAL_ARGUMENT.md §6`):

(a) The **leading c = 7/45 coefficient** comes from the all-singletons monotone partition contribution to ϕ(X^n), which is `(κ_1^B(X))^n / n!`. The all-singletons partition evaluates n independent κ_1^B values — no peak rule at non-adjacent repeats is invoked at this contribution.
(b) The **diagonal κ_2^B** subdominant amplitude (W2 closed form `1/(2·15)`) uses cross-step κ_2 = 0 (W2 verified) and diagonal κ_2 at a single step, no cross-index peak rule.
(c) **Higher-order corrections** at non-adjacent repeats: NOT salvageable via HS 2014 Thm 3.4. The whole Thm 3.4 expansion has terms that depend on Defn 2.2 holding at all index sequences.

So the LEADING-order c = 7/45 may survive, but the FULL Thm 3.4 expansion does not. The Track A integration's headline "rigorous conditional on H1'" needs to be downgraded to "rigorous conditional on H1' AND on the all-singletons + diagonal-κ_2 partial-truncation reading of Thm 3.4 being separately justified."

### What's NOT settled (Mode-E gaps)

#### Gap D1 — the operator-valued vs vacuum-pairing distinction

The diagnostic 0.10783 is a vacuum-state scalar pairing, NOT the operator-valued B_∞-conditional expectation. A clean Mode-E check would re-run the diagnostic computing `E_{B_∞}(X̃_{j_1} X̃_{j_2} X̃_{j_1})` and `E_{B_∞}(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2})` AS B_∞-VALUED OPERATORS (functions of the b_{[1, k]}-accumulator vectors), to confirm whether they are zero as operators or only zero in a vacuum-state pairing. The structural argument in `H1_PRIME_STRUCTURAL_ARGUMENT.md §3.5` predicts non-zero at n = 4 alternating; a direct operator-valued numerical check (~1–2 hours) would confirm.

#### Gap D2 — the "right framework" question

Syracuse's X̃_j family satisfies:
- Second-order independence (cross-step κ_2 = 0)
- Distinct-index higher-order independence
- BUT NOT the full monotone-independence peak rule at non-adjacent repeats

This is a weaker structure than monotone. There may be a published framework that names it:
- **Free**: no — third-order alternating moment is non-zero, contradicting free independence (Voiculescu axiom).
- **Boolean**: cross-step κ_2 = 0 is consistent with Boolean (Speicher-Woroudi), but Boolean has stronger requirements (vanishing of all non-singleton-block cumulants), which Syracuse doesn't satisfy.
- **Monotone**: FAILS at non-adjacent repeats per this audit.
- **Filtered/level-graded variants**: Hasebe 2011 differential, Muraki 2003 associative — these have stronger structure than Defn 2.2 but the additional structure may or may not match Syracuse. Worth scanning the closure-hunt corpus for the right named framework.

#### Gap D3 — the "is the c = 7/45 derivation actually correct?" question

Even though the lift via Thm 3.4 is unsound, the c = 7/45 leading coefficient computed by `MONOTONE_CUMULANTS_C_ASYMPTOTIC.md` matches numerics to high precision. Either:

(i) The derivation is correct for some other structural reason (the all-singletons partition contribution is a robust mechanism that survives even partial framework failure), OR
(ii) The derivation gets the right number for the wrong reason, OR
(iii) There's a different framework that genuinely justifies it.

The Track A integration's verdict "rigorous conditional on H1'" assumed (iii). With H1' failing, (i) is the most plausible — the all-singletons contribution to ϕ(X^n) IS computable from `(κ_1^B)^n / n!` without invoking cross-index peak rules, and this gives 7/45 directly. But making this rigorous as a STANDALONE statement (not a consequence of full Thm 3.4) requires separate justification.

---

## 2. Load-bearing hypotheses

H1' as stated does NOT hold. There is no salvage by switching B's or centering conventions — the structural mismatch (single X̃_j per step vs iid-copies-of-A_λ in monotone framework) is fundamental to the setup, not an artifact of choosing the wrong σ-algebra.

What WOULD make Defn 2.2 hold:

- **(Alt-1)** Switch the Syracuse model to have iid copies of X̃_j at each occurrence. This changes the physics (Syracuse iteration does NOT generate independent copies; each step is a single new piece of randomness). Not a valid model change.
- **(Alt-2)** Restrict to monomial orders/index sequences that do NOT have non-adjacent repeats. This excludes the cross-index alternating moments which are precisely the higher-order corrections in the asymptotic expansion. Restricts the framework to the leading and second-order pieces (which is roughly what's salvageable per §1).
- **(Alt-3)** Adopt a different (weaker) independence framework that the Syracuse X̃_j family DOES satisfy. The candidate frameworks would need to be identified separately.

---

## 3. Recommended next steps

### Immediate (1–2 hours)

1. **Run the Gap D1 check.** Modify `verify_monotone_diagnostic.py` to compute `E_{B_∞}(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2})` as an operator-valued B_∞-function. Confirm that it is non-zero (predicted by Argument A §3.5 / Mode-E gap STRUCTURAL-1). This is the FIRST decisive operator-valued check at n = 4.

### Short-term (4-8 hours)

2. **Audit the c = 7/45 derivation** in `MONOTONE_CUMULANTS_C_ASYMPTOTIC.md` to identify which intermediate steps invoke Thm 3.4 at non-adjacent-repeat moments vs all-singletons / diagonal-κ_2 only. If all the load-bearing steps are at the latter, the leading c = 7/45 is salvageable via a narrower framework. If any are at the former, the derivation needs reconstruction.

### Medium-term (1–3 days)

3. **Identify the correct framework.** The Syracuse X̃_j family has: cross-step κ_2 = 0, distinct-index higher κ_n = 0, non-trivial alternating-repeat moments. Scan the closure-hunt corpus (Hasebe monograph, Muraki 2003, Skeide, Popa 2008 etc.) for a published framework with these exact properties. Candidate names to search: "filtered monotone independence," "level-graded conditional independence," "B-Markov independence," "stationary increments with conditional cumulants."

### Stretch (1–2 weeks)

4. **If no published framework matches,** consider whether the Syracuse-specific structure constitutes a NEW independence framework worth publishing as such. The diagnostic shape (10⁶ separation between distinct-index and alternating-repeat moments at fixed scalar reduction) is striking and may be characteristic of a class of dynamical-system-derived operator families.

---

## 4. Track A integration update

The Track A integration writeup (`TRACK_A_INTEGRATION.md`) §1.4 claimed:

> "The lift is rigorous **conditional on H1'**: that the X̃_j family satisfies HS 2014 Defn 2.2 (monotone independence over `B_marginal`) at all monomial orders.
> ... H1' verification = direct check of Defn 2.2 on the A_j subalgebras at all monomial orders. Estimated effort: 1-2 days focused work at user pace."

**Update.** H1' verification is COMPLETE (this disposition + the three companion files). **H1' does NOT hold.** Effort came in at ~3-4 hours actual including reading + writing the four files.

The Track A line "**c = 7/45 leading coefficient**: Rigorous conditional on H1' (HS 2014 Thm 3.4 + R75 + R77 + R64.B)" should be downgraded to:

> **c = 7/45 leading coefficient**: rigorous AT THE ALL-SINGLETONS PARTITION LEVEL only (modulo separate audit confirming this is the load-bearing contribution); HS 2014 Thm 3.4 does NOT apply at full strength to the Syracuse X̃_j family due to single-X̃_j-per-step vs iid-copies framework mismatch.

The leading derivation is **CLOSER TO** rather than **AT** publication grade. Whether it's publishable in its current form depends on the Gap D1 + D3 follow-ups.

---

## 5. Honesty audit

The brief framed this task as "verify H1' — likely 1-2 day dig that closes the residual gap." The expected outcome was a clean YES verdict, making c = 7/45 unconditional.

The actual finding is a clean NO verdict, in a way that surfaces a deeper structural mismatch between HS 2014's framework and the Syracuse operator algebra. This re-opens an issue that the W1 Track-A closure thought it had identified-but-not-yet-verified.

In Mode-E discipline, the right response is to report the finding clearly, not soften it. The disposition is firm: H1' verbatim does NOT hold; the c = 7/45 leading derivation needs reconstruction or a new framework.

The diagnostic 0.10783 + 10⁶ separation are still real and structurally meaningful — they just don't land in HS 2014's monotone framework cleanly. The "monotone framework" identification of the prior session was **partial-match**, not "full lock."

---

## 6. Mode-E gaps remaining (consolidated)

| Gap | Tag | Description | Effort to close |
|---|---|---|---|
| D1 | Operator-valued check | Confirm n = 4 alternating E_{B_∞} value is non-zero AS OPERATOR, not just as vacuum pairing | 1-2 hrs |
| D2 | Right-framework search | Find published independence framework matching Syracuse's exact moment structure | 1-3 days lit scan |
| D3 | Derivation audit | Identify which c = 7/45 steps invoke Thm 3.4 at problematic moments | 4-8 hrs |
| STRUCTURAL-1 | n ≥ 4 alternating-repeat moments fail Defn 2.2 | Fundamental, not closeable within HS 2014 | Requires alternative framework |
| LO-1 (from low-order checks) | n = 4 j_1, j_2, j_1, j_2 case | Same as STRUCTURAL-1 | Same |

---

## 7. Files

- Verbatim: `C:/Collatz/H1_PRIME_VERBATIM.md`
- Low-order checks: `C:/Collatz/H1_PRIME_LOW_ORDER_CHECKS.md`
- Structural argument + gap audit: `C:/Collatz/H1_PRIME_STRUCTURAL_ARGUMENT.md`
- This disposition: `C:/Collatz/H1_PRIME_DISPOSITION.md`
- Track A integration (needs update): `C:/Collatz/TRACK_A_INTEGRATION.md` §1.4
- W1 verification (audit-affected): `C:/Collatz/W1_BLIFT_VERIFICATION.md`
- HS 2014 PDF: `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_saigo_2014_operator_valued_monotone.pdf`
