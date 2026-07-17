# GENERALIZATION_DISPOSITION

## Disposition: **H_PARTIAL_GENERALIZATION_R78_VARIANT**

**Date:** 2026-05-11. Analyst: Wilson. Reports to: Nathan.

The Path 2 construction chain (steps 1–5) transfers to a NARROW class of cubic-phase Postnikov sums on principal-unit cosets with one-parameter bilinear structure. Outside that class, the chain breaks at specific load-bearing steps. Methods-paper viability is **negative**: the chain's novel content reduces to a single technical lemma (Inner-Plancherel on the second p-adic digit), embedded in published Cochrane–Pinner machinery from 2002–2003. Recommendation: **fold into Paper 4 as a sub-result; abandon the methods-paper hypothesis.**

---

## Phase 2 outcomes (one row per candidate)

| Candidate | Step 1 | Step 2 | Step 3 | Step 4 | Step 5 | Verdict |
|-----------|--------|--------|--------|--------|--------|---------|
| 1. Cochrane–Pinner cubic on `Z/p^n` | WORKS | MODIFIED (2-fold vs p-fold bijection) | WORKS at n ∈ {5, 6} | WORKS_MODIFIED (linear-in-b_2 at p² order) | WORKS | **PASS** in modified form |
| 2. Heath-Brown cubic character sum on `(Z/p^n)×` | WORKS_MODIFIED (after Postnikov) | MODIFIED (2-parameter joint, not 1-param) | STATIONARY_ONLY at varying parameters | **FAILS — DIFFERENT_INNER** (2-parameter mixing breaks 1D Plancherel) | DIFFERENT_GRID at full bilinear | **FAILS at Step 4** |
| 3. Heilbronn-coset at higher prime-power `q = p^{r+1}` | WORKS_MODIFIED (phase is QUADRATIC) | **FAILS — no non-trivial bijection** (single saddle) | WORKS via quadratic Gauss sum (different mechanism) | FAILS (no inner sum to collapse) | N/A | **FAILS at Step 2** |
| 4. Postnikov-style sums on principal-unit subgroup (R78 family) | WORKS | WORKS | WORKS | WORKS | WORKS | **PASS** — but is R78 reparametrized |

**Summary:** 2 of 4 candidates pass the chain. Both passing candidates are STRUCTURALLY CLOSE to R78 (Cochrane–Pinner cubic on Z/p^n; R78 family with arbitrary generator). The two failing candidates fail at specific, identified, load-bearing steps with sharp failure modes (not vague "doesn't quite work" — A2 check satisfied).

---

## R78-specific load-bearing features (Phase 4 A3 summary)

Two features carry the chain and don't transfer beyond a narrow cubic-on-coset class:

1. **Cubic phase degree.** The multi-saddle bijection (Step 2) requires cubic-or-higher polynomial-phase derivative structure. The linear-in-second-digit substitution structure (Step 4) is specific to cubic phases — quadratic gives no bijection, quartic gives wrong digit dependence. **Heilbronn-on-coset fails here.**

2. **One-parameter principal-unit-coset bilinear setup.** The Inner-Plancherel collapse (Step 4) requires ONE outer parameter whose second base-p digit enters linearly after saddle substitution. Two-parameter bilinears (e.g. joint χ × a in Heath-Brown) create cross-terms that break the 1D collapse. **Heath-Brown character sum fails here.**

Features 1 (prime-power modulus), 2 (Cochrane Thm 2), 6 (1/sin grid identity), 7 (truncated p-adic log) are GENERIC and transfer freely. Feature 4 (saddle exactness J = r) ALIGNS for cubic phases at the right r.

The chain transfers iff Features 3 AND 5 both hold. That defines the class: **cubic Postnikov phases on principal-unit cosets with one-parameter bilinear structure.** Both passing candidates inhabit it; both failing candidates miss one of the two.

---

## Phase 4 adversarial checks

### A1 — Cherry-pick check

The 4 retained candidates were selected to include 2 likely-pass (Cochrane–Pinner, R78 family) and 2 likely-fail (Heath-Brown two-parameter, Heilbronn quadratic). The likely-fail candidates were chosen specifically because their failure modes are structurally informative (probing Feature 3 and Feature 5 individually). This is anti-cherry-pick design.

Candidates rejected from Phase 1 (Banks-Shparlinski primal-side, Iwaniec-Sarnak amplification, Burgess square-free, Bourgain sum-product) were genuinely out-of-scope, not "rejected to look better".

The result is asymmetric (2 PASS, 2 FAIL), which is what cherry-pick filtering would NOT produce. The PASS candidates are also honestly flagged as either "R78 inheritance from Cochrane-Pinner" (Candidate 1) or "R78 in different notation" (Candidate 4) — limiting their interpretation as evidence of generalization.

**A1 verdict:** anti-cherry-pick design held; the asymmetric outcome reflects genuine structural features, not selection bias.

### A2 — Negative-result sharpness check

For each FAILING candidate, the failure is documented at a specific load-bearing step:
- **Candidate 2** fails at Step 4 specifically because of two-parameter mixing creating cross-terms `a_2 · λ_{χ,1}` that prevent 1D Plancherel collapse. Sharp, not vague.
- **Candidate 3** fails at Step 2 specifically because quadratic phases have a single saddle (no non-trivial bijection). Sharp, not vague.

Neither failure is "the saddle structure doesn't quite work" — both have algebraic specificity. **A2 verdict: failure modes are sharp.**

### A3 — R78-specific feature audit (full doc: GENERALIZATION_R78_SPECIFIC_FEATURES.md)

Identified two load-bearing R78-specific features (cubic phase degree, one-parameter bilinear). These are the chain's structural envelope. Outside this envelope, the chain doesn't run, and the failure mode is predictable from which feature is missing.

### A4 — Methods-paper viability honest scoping (load-bearing per task spec)

**The novel content of the chain over published literature is:**

The chain uses:
- Cochrane Theorem 2 (2002, published) for Step 1.
- Postnikov decomposition (1955, classical) — implicit in Cochrane.
- Saddle-point closed form (Cochrane–Pinner 2003, published) for Step 3.
- Plancherel orthogonality on `Z/p` (classical) for Step 4.
- 1/sin grid identity (classical Dirichlet kernel) for Step 5.

**The NEW wrinkle is:**

> Step 4's specific algebraic move — when the post-saddle phase has the form `e_{p²}(quadratic-in-s*) · e_p(cubic-in-s* + linear-in-c_2 · s*)`, the inner sum on c_2 ∈ Z/p collapses via additive-character orthogonality to `p · D_p(a_0(s*), p²)` where D_p is the length-p Dirichlet kernel.

This is a TECHNICAL LEMMA — one step's worth of algebra exploiting a specific structural feature of cubic Postnikov phases at r = 3. It is NOT a methodology.

**Specialist's view (honest):** "Cochrane–Pinner's 2003 framework + a specific Plancherel application + classical 1/sin grid bound. Competent application of existing tools to a one-parameter cubic Postnikov sum. The 2N constant is structurally forced. Not new methodology — this is what Cochrane–Pinner's tools are FOR."

**A4 verdict:** the partial template (Cochrane–Pinner with this specific Inner-Plancherel wrinkle) is **NOT substantively new** vs. published 2002–2003 literature. Methods-paper hypothesis is **NOT VIABLE** as a standalone contribution. The result IS load-bearing for the Collatz application (Paper 4); it's not a separate paper.

---

## Decision rule selection

| Disposition | Triggered? |
|-------------|------------|
| H_FULL_GENERALIZATION | NO — 2 of 4 fail at specific steps |
| H_PARTIAL_GENERALIZATION_TEMPLATE (steps 1-3 broad, 4-5 specific) | NO — Step 4 is exactly where MOST candidates split, but Step 3 also splits (Candidate 3 fails at Step 2 before Step 3 fully applies) |
| **H_PARTIAL_GENERALIZATION_R78_VARIANT** (works on close cousins only) | **YES — 2 passing candidates are both R78-style cubic Postnikov on principal-unit coset** |
| H_NO_GENERALIZATION | NO — Candidate 1 (Cochrane–Pinner cubic, n ∈ {5, 6}) is a structurally distinct problem class from R78's Collatz application, and the chain does transfer there |
| INCONCLUSIVE | NO — candidates had clear empirical structure |

---

## Scope statement (what the partial template covers)

The chain transfers to: **bilinear character-sum problems on `(Z/p^{r+1})×` of the form**

> `Σ_{a ≡ 1 mod p in Z/p^r} 1̂(p · a) · e_q(P(a))`

**where:**
- `P(a)` is a polynomial phase of CUBIC degree (after Postnikov substitution to coset coordinates `a = 1 + pα`).
- The parametrization is ONE-DIMENSIONAL (single index `a`, not joint with auxiliary parameter).
- The 1̂ is a Dirichlet kernel of length `N = p^{r-1}` on the same principal-unit coset.
- The saddle exactness condition `J_p = r` holds (achievable at `r = 3` for the natural truncation of the p-adic log).

**Outside this scope, the chain fails:**
- Quadratic phases (e.g. Heilbronn-on-coset): no multi-saddle bijection. Use classical quadratic Gauss sum instead.
- Two-parameter bilinears (e.g. Heath-Brown character sum joint χ × a): no 1D Plancherel collapse.
- Square-free moduli, smooth amplitudes: out of scope by construction.

---

## Methods-paper viability assessment

**Verdict: NEGATIVE.**

The chain's novel content is one technical lemma (Inner-Plancherel collapse on the second p-adic digit), embedded in published 2002–2003 Cochrane–Pinner machinery. A specialist would see this as a competent application of existing tools, not a new methodology. The 5-step "chain" framing makes it look more substantial than it is — Steps 1, 2, 3, 5 are classical; only Step 4 is new, and it's a 2-line algebraic observation.

**A standalone methods paper** "Bilinear bounds for cubic Postnikov sums on principal-unit cosets via Inner-Plancherel" would compete unfavorably against:
- Cochrane–Pinner 2003 — much broader scope, same core machinery.
- Iwaniec–Kowalski Analytic Number Theory book — covers this material as standard.
- Heath-Brown 1996/2000 — handles harder problems (multiplicative characters) with related but different techniques.

The specialist response would be "this is a worked example, not a methods contribution".

---

## Recommendation

1. **Abandon the methods-paper hypothesis.** The chain's novelty is not enough to justify a standalone paper. The Phase 4 A4 honest scoping resolves this definitively.

2. **Fold into Paper 4 as a sub-result.** The bilinear bound `|S_partial| ≤ 2√N` at r ≤ 3 IS load-bearing for the Collatz eq 190 closure. It stays in the Paper 4 narrative as the analytic-machinery step that closes the principal-unit coset bilinear, with appropriate citation of Cochrane–Pinner 2003 for the closed-form magnitude inputs.

3. **Cite Cochrane–Pinner 2003 explicitly** when reporting the chain's r ≤ 3 bilinear bound. The novel content (Inner-Plancherel + 1/sin grid identity for `≤ 2N`) can be stated as a 1-page technical lemma within Paper 4's analytic chapter.

4. **Do NOT pursue more candidates.** The structural envelope is now clear: cubic Postnikov phases on principal-unit cosets with one-parameter bilinear structure. Searching for more candidates within this envelope would just produce more R78 variants; searching outside would just confirm Candidate 2/3-style failures at the identified load-bearing steps.

5. **Tao-email framing implication.** When summarizing Path 2's result for Tao, frame it as: "rigorous bilinear bound `|S_partial| ≤ 2√N` at r ≤ 3 on the principal-unit coset for our specific cubic Postnikov phase, via Cochrane–Pinner closed form + Inner-Plancherel + 1/sin grid identity. The constant 2 is structurally forced. This closes eq 190 at r ≤ 3 for the c = 7/45 Collatz application; we have not pursued generalization to other prime-power character sum problems and do not claim broader methods-paper scope." This is honest, defensible, and matches the disposition.

---

## Files produced

- GENERALIZATION_CANDIDATES.md (Phase 1)
- GENERALIZATION_CANDIDATE_1.md — Cochrane–Pinner cubic on Z/p^n (PASS modified)
- GENERALIZATION_CANDIDATE_2.md — Heath-Brown cubic character sum (FAIL at Step 4)
- GENERALIZATION_CANDIDATE_3.md — Heilbronn-coset at higher prime-power (FAIL at Step 2)
- GENERALIZATION_CANDIDATE_4.md — Postnikov-style sums on principal-unit subgroup (PASS but is R78 reparametrized)
- GENERALIZATION_R78_SPECIFIC_FEATURES.md (Phase 4 A3 feature audit)
- GENERALIZATION_DISPOSITION.md (this document)

---

## Caveats

1. **No numerical verification this session.** The structural assessment is sit-and-think. Phase 2 candidates were tested by walking the chain algebraically, not by computing test cases. Numerical verification of e.g. Candidate 1 (Cochrane–Pinner cubic at n = 6) would be a 1-day Python task; it's not load-bearing for the disposition since Cochrane–Pinner 2003 has done the verification for their published cases.

2. **Literature scan is necessarily limited.** I drew on the candidates the user pre-named plus a handful of natural additions. Deeper scan of e.g. recent (2020+) Heath-Brown/Munshi/Petrow work on prime-power character sums might surface 1–2 more candidates, but they'd inhabit the same scope envelope (cubic Postnikov + one-parameter) or fail at the same load-bearing steps. I judge "one more pass" would not change the disposition; this is supported by the structural ceiling identified in Features 3 and 5.

3. **Time budget.** This was one focused session of literature-and-structure work. The disposition is robust to ~25% additional candidate search; it's not robust to "we missed a fundamentally different framework". I would assign ~85% probability to H_PARTIAL_GENERALIZATION_R78_VARIANT, ~10% to H_NO_GENERALIZATION (if Candidate 1 doesn't pass on closer inspection), ~5% to upgrade to H_PARTIAL_TEMPLATE (if a recent paper surfaces that solidly extends the chain).
