# Track A Adversarial Audit: W1 + W2 closure claims

**Date:** 2026-05-14
**Auditor mode:** adversarial. PDF text extracted via `pypdf`, verbatim quotes only.
**Sources read directly:**
- `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_saigo_2014_operator_valued_monotone.pdf` (13 pages, full text)
- `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_monotone_probability_theory_monograph.pdf` (131 pages, section 3.5 extracted)
- `C:/Collatz/c_seven_forty_fifth.md` (R75)
- `C:/Collatz/result_76_conservation_law.md` (R76)
- `C:/Collatz/result_77_T_lead_spectrum.md` (R77)
- `C:/Collatz/experiments_output/monotone_diagnostic_n3.json` (W1 numerical anchor)
- `C:/Collatz/experiments_output/result_77_7_eps_exact_through_k8_v2_vec_pool.json` (W2 ε_k anchor)

---

## Closure 1 (W1) audit

### Verbatim check: HS 2014 Thm 3.4

HS 2014 Thm 3.4 (p. 7 of `hasebe_saigo_2014_operator_valued_monotone.pdf`, verbatim from `pypdf` extraction lines 277-282):

> "Theorem 3.4. The following moment-cumulant formula holds:
> ϕ(X_1 ⋯ X_n) = Σ_{π∈M(n)} (1/|π|!) K_π(X_1, ..., X_n)."

**The agent's quotation in W1_BLIFT_THEOREM.md §3 is faithful.** The formula holds in `(A, B, ϕ)` where A is any unital algebra, B any unital algebra ⊆ A, ϕ a B-linear conditional expectation. **No commutativity, faithfulness, or normality of ϕ is required**; HS 2014 §2.1 explicitly says "involutions on algebras are not essential in the scope of this paper, so we do not consider them below." So the abelian-B Syracuse setting is a strict specialization. **Verbatim claim SURVIVES.**

### Verbatim check: HS 2014 Prop 3.5

HS 2014 Prop 3.5 (p. 7, lines 285-288 of extraction):

> "Proposition 3.5. For X = (X_1, ..., X_n) ∈ A^n,
> K_n^{N.X} = N K_n^X."

**Verbatim faithful.** B-additivity over monotone-iid copies. **Verbatim claim SURVIVES.**

### Verbatim check: Hasebe monograph Defn 3.23 / Prop 3.25 (cardinality)

The "(n+1)!/2" cardinality formula is NOT in HS 2014. The agent attributes it to "Hasebe monograph Defn 3.23 / Prop 3.25". I verified by extracting `hasebe_monotone_probability_theory_monograph.pdf` directly:

**Defn 3.23 (verbatim from monograph p. 38, my extraction):**
> "Definition 3.23. Let T be a totally ordered finite set. An ordered set partition π = (ρ, ≤) of T is called a **monotone set partition** if
> • ρ is a noncrossing set partition,
> • if B, B′ ∈ ρ satisfies B ⪯ B' then B ≤ B'.
> The set of monotone set partitions of T is denoted by M(T)."

**Covering relation defined just before (verbatim):**
> "For nonempty subsets B_1, B_2 ⊆ T, we say B_1 covers B_2, denoted as B_1 ⪯ B_2, if min B_1 ≤ i ≤ max B_1 for all i ∈ B_2."

So `B_1 ⪯ B_2` means **B_1 is the OUTER block whose interval contains B_2**. Combined with Defn 3.23: OUTER (B) has SMALLER linear-order index (B ≤ B'). This matches HS 2014 p. 4 convention and matches the W1 / W2 interpretation.

**Prop 3.25 (verbatim from monograph p. 39):**
> "Proposition 3.25. Let T be a totally ordered finite set. The cardinality of M(T) is (|T|+1)!/2."

**Both citations verbatim faithful.**

### Sanity check at n=3 reproduction

The agent claimed `|M(3)| = (3+1)!/2 = 12` and enumerated via 5 underlying NC partitions:
- σ_1 = {{1,2,3}}: 1 ordering
- σ_2 = {{1},{2,3}}: 2 orderings (no nesting)
- σ_3 = {{1,2},{3}}: 2 orderings (no nesting)
- σ_4 = {{1,3},{2}}: 1 ordering ({1,3} outer covers {2}, so {1,3} ≤ {2})
- σ_5 = {{1},{2},{3}}: 6 orderings (no nesting)

Total: 1+2+2+1+6 = 12. ✓ **Reproduction is correct.**

The expansion of HS 2014 Thm 3.4 at n=3 with multi-index (j_1, j_2, j_1):
- σ_1 → κ_3^B(X̃_{j_1}, X̃_{j_2}, X̃_{j_1})
- σ_2 (2 orderings × 1/2!) → κ_1^B(X̃_{j_1}) · κ_2^B(X̃_{j_2}, X̃_{j_1})
- σ_3 (2 × 1/2!) → κ_2^B(X̃_{j_1}, X̃_{j_2}) · κ_1^B(X̃_{j_1})
- σ_4 (1 × 1/2!) → (1/2) κ_2^B(X̃_{j_1}, X̃_{j_1}) · κ_1^B(X̃_{j_2})
- σ_5 (6 × 1/3!) → κ_1^B(j_1)·κ_1^B(j_2)·κ_1^B(j_1)

After applying κ_1^B = 0 (centering), all terms with a κ_1^B factor vanish, leaving:
E_B(X̃_{j_1} X̃_{j_2} X̃_{j_1}) = κ_3^B(X̃_{j_1}, X̃_{j_2}, X̃_{j_1}).

This is the clean consistency identification. **The algebra survives.** Minor wrinkle: B-non-commutativity could affect block-cumulant orderings within a partition, but HS 2014 p. 5 says "we define A_π in the above way, neglecting the order structure of π" for the partition-cumulant recursion. The W1 enumeration honors this.

**JSON empirical re-verification (`monotone_diagnostic_n3.json`):**

| Reduction | M_2.abs | M_3_alt.abs | Ratio M3/M2 |
|---|---|---|---|
| tr_pi | 1.07e-5 | 5.09e-4 | ~50× |
| vac_pi | 5.14e-5 | 2.48e-4 | ~5× |
| delta_1 | 8.57e-6 | 2.41e-3 | ~280× |
| **sum_entries** | **1.08e-7** | **0.10783** | **~10⁶×** |

The "10⁶ separation" is verified — **but only in the `sum_entries` reduction**. The other three reductions give 5-280× separations. W1 cites "10⁶ separation" without flagging this is the most dramatic of four scalar reductions. **Mild cherry-pick caveat.**

### H1 gap assessment

W1_BLIFT_DISPOSITION.md §"Mode-E gaps remaining" item 1 acknowledges H1 (monotone independence of (A_j) over B in the sense of HS 2014 Def 2.2) is "not theorem-grade" and is "supported by framework identification + 10⁶ separation + peak-rule match" — i.e., **the agent is explicit that H1 is the load-bearing project-internal input, not closed by Route 2**.

Adversarial probe: is H1 ACTUALLY load-bearing, or did the agent understate? HS 2014 Def 2.2 requires the factorization

`ϕ(X_1 ⋯ X_n) = ϕ(X_1 ⋯ ϕ(X_i) ⋯ X_n)` whenever λ_{i-1} < λ_i > λ_{i+1}

for ALL n, ALL multi-indices (j_1, ..., j_n), ALL X_i ∈ A_{λ_i}. The Task 1 probe only tests n=3 with one specific peak pattern (j_1, j_2, j_1). **H1 at n=3 is supported numerically; H1 at all higher n is extrapolated.** The agent's framing ("strong evidence but not a verbatim theorem") is honest, not understated.

Additional adversarial point: W1_BLIFT_VERIFICATION.md §4 surfaces a **centering subtlety** — under marginal centering (reading B), `E_B(X̃_{j_2}) ≠ 0` (it's a non-trivial B-measurable function), which means κ_1^B(X̃_{j_2}) under the operative reading is NOT zero. The agent introduces a "refined hypothesis H1'" (level-graded monotone independence over B_marginal) to resolve this. This introduces **a new project-internal hypothesis** that didn't exist in the original closure framing. The disposition acknowledges this in item 2 ("centering subtlety should be propagated to Deliverables B and C of the MONOTONE_CUMULANTS_* chain for full consistency"). **This is a real residual gap — the lift is rigorous conditional on a hypothesis that itself was REFINED during the audit and now requires propagating.**

### Verdict: SURVIVED_WITH_CAVEAT

**Survived:**
- HS 2014 Thm 3.4 and Prop 3.5 are verbatim correctly cited.
- The abelian-B Syracuse setting IS a strict specialization (no extra hypotheses needed).
- The n=3 enumeration of M(3) = 12 monotone partitions is correct.
- The collapse to κ_3^B under κ_1^B = 0 (strict-centering reading A) is algebraically correct.

**Caveats:**
1. **The cardinality "(n+1)!/2" is from the Hasebe MONOGRAPH (Defn 3.23 / Prop 3.25), not HS 2014.** The agent's section header in `W1_BLIFT_VERIFICATION.md` line 35 says "Hasebe monograph Prop 3.25 cardinality formula `(n+1)!/2 = 4!/2 = 12`" — this is correct; the monograph IS the right citation. But the user's prompt referenced "Prop 3.25" attribution; I verified this lands in the right document.
2. **H1 (monotone independence of (A_j) over B) is acknowledged not-theorem-grade.** This is explicit in the disposition.
3. **Centering subtlety surfaces a NEW hypothesis H1' (level-graded monotone independence over B_marginal).** The original closure was framed as "rigorous conditional on H1"; in fact it's "rigorous conditional on H1' under reading B, OR conditional on H1 under reading A where the diagnostic is algebraically zero (10⁻¹⁸)". The reading-B version is the load-bearing one; H1' is not in HS 2014.
4. **The "10⁶ separation" headline is the best of four scalar reductions** (the others give 5-280×).

---

## Closure 2 (W2) audit

### Monotone partition count verification (n=2, 3, 4 manual enumeration)

I directly enumerated monotone partitions with exactly one 2-block + (n−2) singletons via Python (`itertools.combinations` over 2-block positions, then valid linear orderings under Defn 3.23 covering constraint):

| n | manual count (raw # monotone partitions) | (n−1)!·[n·H_{n−1} − (n−1)] | n·H_{n−1} − (n−1) |
|---|---|---|---|
| 2 | 1 | 1 | 1 |
| 3 | 5 | 5 | 2.5 |
| 4 | 26 | 26 | 4.333 |
| 5 | 154 | 154 | 6.417 |

**Finding:** The raw COUNT of monotone partitions with one 2-block + (n−2) singletons is `(n−1)! · [n·H_{n−1} − (n−1)]`, NOT `n·H_{n−1} − (n−1)`.

W2_PARTITION_COUNT.md is **algebraically correct in its derivation** (Step 3 line 113-116 derives `(n−1)! · [n·H_{n−1} − (n−1)]`), but then conflates "count" with "coefficient after dividing by |π|! = (n−1)!" in §4 line 152-156:

> `= (1/(n−1)!) · [number of such monotone partitions] · κ_2^B · (κ_1^B)^{n−2}`
> `= (1/(n−1)!) · (n−1)! · [n · H_{n−1} − (n−1)] · κ_2^B · (κ_1^B)^{n−2}`
> `= [n · H_{n−1} − (n−1)] · κ_2^B · (κ_1^B)^{n−2}`

This is correct: the **moment-cumulant coefficient is `n·H_{n−1} − (n−1)`** because of the `1/|π|!` factor in Thm 3.4. But the section header §4 says "Count = (n−1)! / (k+1)" then sums to `(n−1)! · [n·H_{n−1} − (n−1)]` as the partition count, and §5 line 175 says "The correct count is `n · H_{n−1} − (n−1) ~ n · ln(n)`" — this **conflates count with coefficient**.

W2_DISPOSITION.md §1 line 27-28: "the actual count is `n · H_{n−1} − (n−1)`" — same conflation.

**Severity:** This is a **terminological imprecision**, not a math error. The load-bearing conclusion — that the moment-cumulant CONTRIBUTION grows like `n · ln(n) · κ_2^B · (κ_1^B)^{n−2}`, polynomial-times-log, NOT `(1/2)^n` — survives. The redirect-to-T_M conclusion is unaffected.

### Hasebe Defn 3.23 noncrossing interpretation

**The CRITICAL question** the user raised: does Defn 3.23 force the 2-block to be an interval block, or can it be a non-interval like {1,3}?

**Verbatim Defn 3.23 (from monograph):**
> "An ordered set partition π = (ρ, ≤) of T is called a monotone set partition if
> • ρ is a noncrossing set partition,
> • if B, B' ∈ ρ satisfies B ⪯ B' then B ≤ B'."

A 2-block `{i, j}` with `j > i+1` PLUS surrounding singletons is **noncrossing** — singletons cannot cross anything, and a single 2-block alone never has a crossing. So {1,3} with {2} as singleton IS a valid noncrossing partition. Under Defn 3.23 covering: {1,3} ⪯ {2} ({1,3}'s interval contains 2), so {1,3} (outer) must come BEFORE {2} (inner) in the linear order.

**W2.B's interpretation is correct; the prior writeup MONOTONE_CUMULANTS_C §4 was wrong** (it claimed the noncrossing constraint forces 2-blocks to be intervals — false, as singletons don't cross anything). W2's correction is **verbatim consistent with the monograph.**

### R76 + R75 + R77 source check for 1/30 = 1/(2·15)

**R76 Thm 76.3 verbatim (line 75 of `result_76_conservation_law.md`):**
> "Theorem 76.3 (Leading-mode Identity). For every n ≥ 1, S_{n+1} = M_{n+1}(1) = −2 · M_{n+1}(1 + 3^n) = −2 · M_{n+1}(1 + 2·3^n)."

The **factor of 2** in W2's `1/(2·15)` is the `−2` factor in R76 Thm 76.3. **Verbatim verified.**

**R77 T_diag verbatim (line 12 of `result_77_T_lead_spectrum.md`):**
> "(P_+, P_−)_{n+1, diag} = T_diag · (P_+, P_−)_n, where T_diag = (1/5)·[[1, 1], [4, 4]]."

The **factor of 5** in W2's `15 = 3·5` is the `1/5` prefactor of T_diag. **Verbatim verified.**

**R75 §1 verbatim (line 46 of `c_seven_forty_fifth.md`):**
> "So c = lim ‖d_{k+1}‖² · 3^k = lim S_{k+1}/3 = S_∞/3."

The **factor of 3** in W2's `15 = 3·5` is attributed to "R75 Plancherel global factor (each level-n+1 high-freq mode has 3 lifts to level-(n+1) cosets)". The "3 lifts" claim is verbatim from R75 §3 line 93 ("three lifts (ξ_0, ξ_0+3^n, ξ_0+2·3^n)"). But the `1/3` factor itself shows up in `c = S_∞/3`, not directly as a subdominant amplitude prefactor.

**Adversarial reading:** W2's attribution `1/30 = 1/(2·15) = 1/(2·3·5)` is a **structural decomposition** of the empirical constant into three factors, each individually rigorous, but the **assembly** of these into `1/30` is essentially a numerical coincidence-check — the three factors don't algebraically multiply through R76/R75/R77 to yield `1/30` automatically. The agent acknowledges this in W2_DISPOSITION.md §6: "**The amplitude `α = 1/60` of `R_n − R_∞` (equivalent to the eigenvector amplitude of T_M at eigenvalue 1/2) is not a Hasebe combinatorial output**." So the decomposition `1/30 = 1/(2·15)` is **numerically consistent with the three rigorous factors but is not a rigorous derivation of `1/30`**. The agent is candid about this.

### Rate redirect claim assessment

W2 claims the rate `(1/2)^n` does NOT come from monotone-partition counting (since that gives `n·ln(n)` polynomial growth) but must come from **T_M's λ_2 eigenvalue**.

**R77 line 53:** "Conjecture 77.2: The full operator T (T_diag + Off_n linearization) has subdominant eigenvalue λ_2 = 1/2 acting on the (1, 4) deviation subspace." This is a **conjecture**, not a theorem.

**R77 line 64:** "Working through this: the leading off-diagonal eigenvalue of the v ≠ v' contributions is ~ 2·(1/4)·(weight) = 1/2 when summed over the leading bilinear couplings." This is a **heuristic** for why λ_2 = 1/2 plausibly emerges.

**R77 line 122:** "Computing it exactly over Q... and applying Nisoli's Theorem 2.15 with explicit error bound from Tao's Prop 1.17, would close the rigor gap. **Estimated effort: another session.**"

**Verdict on the rate-redirect:** W2 correctly identifies that monotone-partition combinatorics cannot supply the `(1/2)^n` rate. The redirect-to-T_M is HONEST: T_M's `λ_2 = 1/2` is empirical/conjectural (R77 Conj 77.2), not theorem-grade. **W2 does NOT close the rate; it redirects to an outstanding separate problem.** This is exactly what the agent claims ("Wrinkle 2 redirects to T_M closure, a separately-open problem"). **Not a punt — it's a clear-eyed redirection.**

### Cross-step κ_2^B = 0 check

W2_KAPPA2_CALC.md §1 line 30: "Cross-step κ_2^B at j_1 < j_2 vanishes structurally under marginal centering (per Deliverable B §2.2): `κ_2^B(X̃_{j_1}, X̃_{j_2}) = 0`."

**Adversarial test against W1 sanity check:** at n=3, the W1 sanity check identifies the non-zero moment M_3_alt = 0.10783 with `κ_3^B(X̃_{j_1}, X̃_{j_2}, X̃_{j_1})` — a THIRD-order cumulant, not second-order. So W2's `κ_2^B(X̃_{j_1}, X̃_{j_2}) = 0` at distinct steps does NOT directly contradict W1: the n=3 non-zero is in κ_3, not in cross-step κ_2.

**However:** if `κ_2^B(X̃_{j_1}, X̃_{j_2}) = 0` for j_1 ≠ j_2 structurally (Deliverable B §2.2 assumes "Pascal-pair independence at distinct steps"), this is a **strong** structural assumption. The W1 expansion (★) at n=3 has terms `κ_2^B(X̃_{j_2}, X̃_{j_1})` (mixed indices) which would vanish under W2's assumption — but these terms vanish anyway under W1's centering (they have a κ_1^B factor). So the two are **consistent in the n=3 reduction** (both yield M_3_alt = κ_3^B in the leading non-vanishing contribution), but **W2's cross-step κ_2 = 0 is an additional structural assumption beyond W1's H1, not derived from HS 2014 alone**.

### Empirical-value re-verification from JSON

I extracted `result_77_7_eps_exact_through_k8_v2_vec_pool.json` and computed `|ε_k|·2^k`:

| k | ε_k (exact, decimal) | sign | \|ε_k\|·2^k |
|---|---|---|---|
| 1 | +1/5 = 0.2 | + | 0.4 (transient) |
| 2 | +1/105 = 9.52e-3 | + | **0.0381** |
| 3 | -5191/1019445 = -5.09e-3 | − | **0.0407** |
| 4 | -2.45e-3 | − | **0.0392** |
| 5 | -1.15e-3 | − | **0.0369** |
| 6 | -4.98e-4 | − | **0.0319** |
| 7 | -1.18e-3 | − | 0.1504 (multi-spectral) |
| 8 | -7.46e-4 | − | 0.1909 (multi-spectral) |

**Discrepancy with W2_DISPOSITION.md §5 line 63:**
> "`|ε_n|·2^n` plateau (n=2..6): 0.0381, 0.0407, 0.0392, 0.0369, **0.0349**"

The actual k=6 value is **0.0319**, not 0.0349. W2 reports 0.0349. **This is a small empirical bookkeeping error.** Magnitude: 0.0319 vs 0.0349 — a ~10% deviation in the reported decimal value. The structural conclusion (k=2..6 plateau approaches 1/30 from above; k=7,8 deviate due to multi-spectral transient) is **unaffected** — in fact 0.0319 is BELOW 1/30 = 0.0333, suggesting either non-monotone convergence or that the multi-spectral transient onsets at k=6 rather than k=7.

Also: W2_KAPPA2_CALC.md line 152 lists "n=2..6 plateau approaches `1/30 ≈ 0.0333`" — this statement is roughly true with the corrected values (0.038, 0.041, 0.039, 0.037, 0.032), though now hitting 0.032 (below 0.033) at k=6 is consistent with the conclusion but raises a question about exactly when multi-spectral transient begins.

### Verdict: SURVIVED_WITH_CAVEAT

**Survived:**
- The combinatorial claim that monotone partitions of `[n]` with one 2-block + (n−2) singletons grow like `n·ln(n)` (polynomial-times-log) in n, NOT `(1/2)^n` — **verified by manual enumeration at n=2,3,4,5**.
- The correction to MONOTONE_CUMULANTS_C §4 (that the 2-block need NOT be an interval block under noncrossing) is **verbatim correct against Hasebe monograph Defn 3.23**.
- R76 Thm 76.3 `S_n = −2·R_n` verbatim verified — factor of 2 in `1/(2·15)` is rigorous.
- R77 T_diag = (1/5)·[[1,1],[4,4]] verbatim verified — factor of 5 in `1/(2·15)` is rigorous.
- The rate-redirect to T_M (R77 Conj 77.2 is empirical/conjectural; closing T_M is a separate session) is **honest and correctly attributed**.

**Caveats:**
1. **Terminological imprecision in W2_PARTITION_COUNT.md §4-§5 and W2_DISPOSITION.md §1**: "count = `n·H_{n−1} − (n−1)`" is wrong; raw count is `(n−1)! · [n·H_{n−1} − (n−1)]`. The moment-cumulant COEFFICIENT (after the 1/|π|! factor) is `n·H_{n−1} − (n−1)`. The math is correct in the derivation; the wording later in the doc conflates count and coefficient.
2. **Empirical `|ε_6|·2^6` value**: W2 reports 0.0349; actual is 0.0319. Small bookkeeping error; doesn't change conclusions.
3. **Cross-step κ_2^B = 0 is an additional structural assumption** beyond HS 2014. Consistent with W1 at n=3 (both yield M_3_alt = κ_3^B), but not derived; carried from "Deliverable B §2.2" (Pascal-pair independence claim).
4. **The `15 = 3·5` decomposition is a structural pattern-match, not an algebraic derivation.** W2 acknowledges in §6: amplitude `1/60` is the T_M eigenvector amplitude, empirical, not a Hasebe combinatorial output. The `1/30 = 1/(2·15)` decomposition is rigorous-as-a-pattern, not as a derivation.

---

## Cross-claim consistency check

### Do W1 and W2 contradict each other?

**No direct contradiction found.**

- W1 uses HS 2014 Thm 3.4 + Prop 3.5; W2 uses Hasebe monograph Thm 3.26 (the scalar form, p. 39). HS 2014 Thm 3.4 is the operator-valued version of monograph Thm 3.26 — they agree on the formula `ϕ(X^n) = Σ_{π∈M(n)} (1/|π|!) κ_π`.
- The cardinality `(n+1)!/2` (Monograph Prop 3.25) is consistent between W1 (used at n=3 to get 12) and W2 (cited in §1).
- The "Defn 3.23 covering relation" interpretation (outer-first / inner-later in linear order) is consistent between W1 enumeration of M(3) and W2 enumeration of one-2-block monotone partitions.
- W1's `κ_3^B` non-zero (lands the diagnostic) is at n=3 third-order; W2's `κ_2^B(j_1, j_2) = 0` at distinct steps is at second-order — different orders, no direct conflict.

### Does W1's moment-cumulant formula agree with W2's partition count?

**Yes.** W1 applies Thm 3.4 at n=3 with multi-index (j_1, j_2, j_1), sums over all 12 partitions in M(3), gets the explicit expansion. W2 sums over only the one-2-block partitions across general n and gets `[n·H_{n−1} − (n−1)] · κ_2^B · (κ_1^B)^{n−2}`. **W2's count is a substring of the full Thm 3.4 sum, consistent with W1's framework.** At n=3, the one-2-block partitions in M(3) are σ_2 (2 orderings), σ_3 (2 orderings), σ_4 (1 ordering) — total 5 partitions. Formula check: `n·H_{n−1} − (n−1) = 3·1.5 − 2 = 2.5`, and `(n−1)! · 2.5 = 5`. **Consistent.**

---

## Overall verdict on Track A so far (W4 pending)

**Both closures SURVIVED with non-fatal caveats.**

- **W1:** SURVIVED_WITH_CAVEAT. HS 2014 Thm 3.4 + Prop 3.5 are correctly cited and apply verbatim. The B-amalgamated lift is rigorous CONDITIONAL on H1 (monotone independence of (A_j) over B). H1 is supported by 10⁶ numerical separation (in the best of four scalar reductions) + framework identification + peak-rule match. The audit surfaced a refinement: H1' (level-graded monotone independence over `B_marginal`) is what's actually load-bearing under reading B, not the original H1 — this is acknowledged in the agent's writeup but introduces an additional project-internal hypothesis that needs propagation.

- **W2:** SURVIVED_WITH_CAVEAT. The combinatorial finding (one-2-block monotone partitions grow `n·ln(n)`, NOT `(1/2)^n`) is **verified by manual enumeration**. The correction to MONOTONE_CUMULANTS_C §4 is **verbatim correct** against Hasebe Defn 3.23. The factor-of-2 (R76 Thm 76.3) and factor-of-5 (R77 T_diag) attributions are rigorous. The factor-of-3 attribution (R75 Plancherel) is loose — it's a global-decay-rate factor showing up in `c = S_∞/3`, not directly an amplitude prefactor. The rate-redirect to T_M is honest (T_M's λ_2 = 1/2 is empirical Conj 77.2, not theorem-grade) — it does NOT close the `1/30` derivation; it redirects to a separately-open problem. Empirical `|ε_6|·2^6` is reported as 0.0349 but actually 0.0319.

**No fatal errors. No contradictions between W1 and W2. Both closures are honest about their residual gaps, with minor terminological/empirical bookkeeping issues that don't affect the load-bearing conclusions.**

The track has two genuine open problems escalated by these audits:
1. **W1 escalation:** H1' (level-graded monotone independence over B_marginal under reading B) needs theorem-grade proof or a switch to reading A (where the diagnostic is algebraically zero, and a different mechanism is needed for the non-trivial 0.10783).
2. **W2 escalation:** T_M spectral analysis (R75 §8 / R76 §6 / R77 §6) — Conj 77.2 (λ_2 = 1/2) and the eigenvector amplitude (= 1/60) need rigorous derivation. Estimated "another session" per R77 line 122.

W4 (Faure √3 ↔ cumulant op) is pending and not audited here.
