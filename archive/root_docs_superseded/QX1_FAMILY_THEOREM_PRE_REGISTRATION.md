# Pre-Registration: qx+1 Family-Level Plancherel Saturation Theorem (Move 2)

**Pre-registered:** 2026-05-11T11:15 EDT, before any compute on this task.
**Author:** Claude.
**Hypothesis frame:** prompt provided by user, reconciled with framework definitions below.

---

## §0. Framing reconciliation (critical context — locked before compute)

The prompt describes "qx+1 Plancherel saturation theorem at family level" with "Rate-1/2 universal across q." After reading `q_sweep_results.md`, `result_q_sweep_test_1_rate.md`, `result_q_sweep_test_3_decomposition.md`, `result_78_FINAL.md`, `result_78_extended.md`, and the supporting code, the project contains **two distinct rate-1/2 results that must not be conflated**:

**(A) ε_n rate-½ envelope** — `|ε_n^{(q)}| · 2^n` flat at q=3 only.
The empirical envelope test `|ε_n^{(q)}| · 2^n ∈ [0.032, 0.041]` for q=3 (R75/R77.x) does NOT generalize. `result_q_sweep_test_1_rate.md` Verdict: "RATE-MIXED" — for q ≥ 5, `S_n^{(q)}` itself diverges geometrically (growth ratio q/3 > 1), so there is no finite `S_∞^{(q)}` and no rate-½ envelope to test. The q=3 case is structurally the boundary (q/3 = 1, sequences just converge).
**The "rate-½ universal across q" claim at the ε_n level is FALSE per existing project data.**

**(B) K_p √N saturation exponent** — empirically universal across q ∈ {3, 5, 7, 11, 13}.
The Kalafatelis-sum saturation exponent `β` with `|K_p(r, c, m)| ~ N^β` is empirically `β ∈ [0.483, 0.518]` at canonical (c=1, m=0), spanning 5 primes (Pattern β classification per `q_sweep_results.md`). This IS family-level rate-½. Prefactor `C_p` varies ~1.4× across q with no closed-form arithmetic invariant identified.

**These are different objects.** (A) is convergence-rate-across-levels of a Plancherel sum. (B) is character-sum-magnitude-at-fixed-level.

The Move 2 prompt's "rate-1/2 universal" must refer to (B), the K_p saturation, because (A) is q=3-specific by data. With that substitution, the prompt's strategy becomes coherent: prove K_p √N saturation rigorously at family level, specialize to q=3 to recover an explicit `|μ̂_n(ξ)|` bound, plug into Nisoli framework as ε_K replacing Tao Prop 1.17's qualitative C_A.

**Locked pre-reg interpretation:** "rate-½" = K_p √N saturation exponent (object (B)), NOT the ε_n envelope (object (A)). All three claims are evaluated against this interpretation.

This framing distinction is documented in §0 of the eventual results document and cannot be retroactively modified.

---

## §1. The theorem to attempt (precise statement)

**THEOREM (qx+1 family-level Plancherel saturation, TARGET):** For every prime `p ≥ 3` and every `r ≥ r_0(p)` (some explicit threshold), the Kalafatelis-sum
> `K_p(r, c, m) = Σ_{u=0}^{N_p − 1} e_M(c·(1+p)^u − p²·m·u),  M = p^{r+1},  N_p = p^{r-1}`

satisfies, for `(c, m) ∈ Z/M × Z/M`:
> `|K_p(r, c, m)| ≤ f(p) · √N_p`

with `f(p)` an explicit closed-form function of `p`, derivable by structural argument (not curve-fitting).

**Three load-bearing claims:**

- **CLAIM 1 (rate-½ universal, exact):** the saturation exponent of `|K_p|` in `N_p` is **exactly ½** for every prime p ≥ 3 — structurally, not just empirically.
- **CLAIM 2 (prefactor closed form):** `f(p)` has an explicit closed form. Empirically `C_p ∈ [0.83, 1.17]` across q ∈ {3, 5, 7, 11, 13} (canonical metric) ranges by ~1.4×; the closed form must predict these values within their observation precision.
- **CLAIM 3 (q=3 specialization gives Nisoli input):** at p=3 the family-level bound specializes to an explicit `|μ̂_n(ξ)| ≤ f(3) · 3^{-(n-1)/2}` (or whatever the Plancherel-normalized version is), which translates to explicit ε_K = ‖T − T_K‖ feeding Nisoli Lemma 2.9 with η < 1 in the verified contour range. This dissolves R77.2's conditional reliance on Tao Prop 1.17's qualitative C_A.

**Pre-registered expectation: NULL** — at least one claim fails to close rigorously. The default expectation is that the family-level generalization of R78.1–78.3's q=3 proof encounters one of three known obstructions:

- (Obs1) R78.3's "support {a ≡ 1 mod 3}" sparsity argument is q=3-specific (uses the 4 = 1+3 principal-unit structure and 3-adic Pontryagin duality). Generalizing the support set to `{a ≡ 1 mod p}` may not yield the same exact-magnitude formula at general p.
- (Obs2) The 1.4× prefactor variation across q has no arithmetic invariant explanation identified in `q_sweep_results.md` §5. Closed-form `f(p)` may need additional structural input.
- (Obs3) Even if K_p saturation generalizes, translating it to a `|μ̂_n(ξ)|` bound at q=3 may not be in the right form for Nisoli's ε_K because the Plancherel-truncation operator's exact structure is distinct from the Kalafatelis-sum geometry.

Override of NULL requires structural proof closing all three obstructions, not empirical extension to more primes.

---

## §2. Locked procedure

### Phase 1: Structural formulation (pre-compute analysis)

- Read `result_78_FINAL.md`, `result_78_extended.md`, `r79b_*` notes for the q=3 R78.3 proof apparatus.
- Define precisely the family-level Plancherel formula for K_p(r, c, m).
- Identify the family-level structural object whose spectrum gives the universal saturation exponent (analogous to the principal-unit cyclic group structure at q=3 with order 3^r in (Z/3^{r+1})*).
- **Deliverable:** explicit statement of (1) the family-level Plancherel identity, (2) the family-level support of F̂_p, (3) the structural saturation magnitude formula.

### Phase 2: Claim 1 attempt — rate-½ universal, exact

Candidate routes (must commit to one route or document why all fail):

- **Route (a):** Direct algebraic generalization of R78.1–78.3 to family level. Prove for general prime p:
  - (78.1_p) complete-sum vanishing: `Σ_{u=0}^{p^{r+1}-1} e_{p^{r+1}}(c·(1+p)^u − p²mu) = 0` (Cochrane Theorem 2-style argument, no q=3-specific dependence).
  - (78.2_p) Fourier sparsity: `supp(F̂_p) ⊆ p·Z/p^{r+1}` with explicit cardinality.
  - (78.3_p) magnitude equidistribution: `|F̂_p(ξ)| = M/√|supp|` on its support, by Plancherel + uniformity.
- **Route (b):** Structural invariant from sibling-prime / Pattern β. Use `result_q_sweep_test_3_decomposition.md`'s lift-residual orthogonal decomposition (q-universal, proved over Q at all tested q) to derive a bound on K_p.
- **Route (c):** Family-level operator construction with universal ½ in its spectrum (an analog of T_3 from R77.2 at family level).

**Threshold:** Claim 1 PASSES at family level if at least one route produces a structural proof that the saturation exponent is exactly ½ for all primes p ≥ 3, without relying on observed values at q ∈ {3, 5, 7, 11, 13}. PARTIAL pass if a route works at general p but with a not-yet-effective constant.

### Phase 3: Claim 2 attempt — prefactor closed form

- If Phase 2 Route (a) closes: derive `f(p)` from R78.3_p's magnitude formula directly. At q=3 the formula `|F̂(ξ)| = 3√q` from R78.3 suggests `f(p) = ?√p` with explicit numerator — check against empirical `C_p`.
- If Phase 2 Route (b) or (c) closes: derive `f(p)` from the operator/orthogonal structure used in that route.
- Verify against empirical `C_p ∈ [0.83, 1.17]` at q ∈ {3, 5, 7, 11, 13}. Predictions must match within the observation precision (`σ_β ∈ [0.012, 0.021]`, comparable σ on C).
- **Threshold:** Claim 2 PASSES with closed-form `f(p)` matching empirical at all 5 tested primes within their CI. PARTIAL pass with closed form matching at q=3 only.

### Phase 4: Claim 3 attempt — Nisoli specialization at q=3

- If Claims 1 + 2 close, specialize the family-level bound to q=3 and verify it gives an explicit `|μ̂_n(ξ)| ≤ const · 3^{-(n-1)/2}` (or appropriately Plancherel-normalized).
- Translate to operator-norm bound `‖T − T_K‖ ≤ ε_K` via R77.2's framework.
- Verify `η = ε_K · M_3 < 1` for K in R77.2's verified range (K ≥ some K_0).
- Plug into Nisoli Lemma 2.9, derive explicit `|λ_2 − 1/2|` bound, conclude rate-½ for c=7/45 rigorously.

### Adversarial safeguard A2: verify at q ∈ {17, 19, 23}

If Claim 1 closes structurally with a closed-form prediction, run the empirical K_p sweep at p ∈ {17, 19, 23} at a depth comparable to existing p ∈ {11, 13} (r ≤ 7). Verify predicted β = ½ + ε holds within observation precision. **If predictions miss for any of these three primes, Claim 1 is falsified at family level.**

If Claim 1 does not close structurally, skipping A2 is permitted (no structural claim to verify against).

### Adversarial safeguards A1, A3, A4

- **A1:** distinguish structural proof from empirical extension. Each phase report names explicitly what is rigorous vs what is empirical-with-anchor vs what is conjecture.
- **A3:** route disagreement is a flag. If Phase 2 routes (a), (b), (c) disagree, the strongest doesn't override the weakest unless the weaker routes' failures are explained.
- **A4:** any deviation from this pre-registration is logged in writing with explicit reason before continuing.

### Decision rules (locked, null-favored)

- **THEOREM_PROVEN:** Claims 1, 2, 3 all close rigorously. A2 predictions hold at q ∈ {17, 19, 23}. Adversarial safeguards clean. Family-level theorem real; R77.2 conditional dissolves; c=7/45 rigorously closed via family-level mechanism instead of Tao §7 extraction.
- **CLAIM_1_ONLY:** rate-½ universal proven structurally but `f(q)` doesn't close. Standalone qx+1 saturation result (rate-½ universal) publishable; c=7/45 conditional on different input still.
- **CLAIMS_1_AND_2:** both close; Claim 3 has translation gap (Nisoli ε_K input form mismatched). Family theorem real; c=7/45 closure requires additional work.
- **EMPIRICAL_UNIVERSAL_NOT_STRUCTURAL:** rate-½ empirically universal across observed q but no structural proof at family level. Status unchanged from `q_sweep_results.md`; documents the obstruction.
- **STRUCTURAL_OBSTRUCTION_FOUND:** one of the claims falsifies under attempt or a specific known obstruction (Obs1 / Obs2 / Obs3) is precisely identified and shown to block the route. Family-level theorem as stated doesn't exist; the obstruction informs the next move.

---

## §3. Deliverables

- `QX1_FAMILY_THEOREM_PRE_REGISTRATION.md` — this document, committed before compute.
- `QX1_FAMILY_THEOREM_ATTEMPT.md` — disposition at top, all four phases, A1–A4 records.
- Per-phase artifacts (csvs / source) as appropriate.
- If A2 fires: q ∈ {17, 19, 23} verification data.
- If STRUCTURAL_OBSTRUCTION_FOUND: precise statement of the obstruction with mechanism.

## §4. Honest scope-of-attempt statement

This is a proof-attempt task with a serious mathematical target. The pre-registered NULL expectation is informed by:
1. The Move 2 prompt's framing required correction (§0 above) before the test could even begin.
2. R78.3's exact magnitude proof at q=3 uses q=3-specific Pontryagin/Plancherel structure that does not obviously generalize.
3. The 1.4× prefactor variation across q already has no arithmetic invariant explanation per `q_sweep_results.md` §5.

If a specific obstruction is identified precisely (mechanism named, not just "this didn't work"), that is itself a real result. The honest disposition spectrum is THEOREM_PROVEN → CLAIM_1_ONLY → CLAIMS_1_AND_2 → EMPIRICAL_UNIVERSAL_NOT_STRUCTURAL → STRUCTURAL_OBSTRUCTION_FOUND. Pre-registration locked.
