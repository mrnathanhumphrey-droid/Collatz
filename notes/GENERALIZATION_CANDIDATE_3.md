# GENERALIZATION_CANDIDATE_3 — Heilbronn-coset sum at higher prime-power modulus

**Candidate object:** `H(a, q) := Σ_{n: n ≡ 1 mod p, n ∈ Z/p^r} e_q(a · n^p)` where `q = p^{r+1}`, restricted to the principal-unit coset (n ≡ 1 mod p).

**Background:** The classical Heilbronn sum is `Σ_{n=1}^p e_{p²}(a·n^p)` (Heath-Brown 1996). Heath-Brown showed `≪ p^{11/12}`. The natural higher-power-modulus extension `q = p^{r+1}` summed over the principal-unit coset is the "Heilbronn-on-coset" object.

**Empirical signal:** On principal-unit coset, n^p has special structure: `(1+pα)^p = 1 + p²α + ... mod p^3` — the p-th power map collapses to a "linear-shift-by-p" operator. This makes the phase nearly trivial at r = 2 and progressively more structured at r ≥ 3.

**Pre-assessment:** STEP1 likely pass (Postnikov-like collapse); STEP3 fragile (Heilbronn phase has different saddle structure than R78 cubic).

---

## STEP 1 — closed-form magnitude theorem

**R78 analog:** F̂_p closed form via Cochrane Theorem 2.

**Candidate test:** On principal-unit coset, substitute `n = 1 + pα` with `α ∈ Z/p^{r−1}`:
> `n^p = (1+pα)^p = Σ_{k=0}^p C(p, k) · (pα)^k`
> ≡ `1 + p · pα + C(p,2) · p²α² + ... mod p^{r+1}`
> ≡ `1 + p²α + p²·(p−1)/2 · α² · p + ... = 1 + p²α + p³·(p-1)α²/2 + ...`

Phase becomes:
> `e_q(a · n^p) = e_{p^{r+1}}(a · (1 + p²α + p³·(p-1)α²/2 + p⁴·... ))`

For `r = 3` (q = p⁴): the phase is `e_{p⁴}(a) · e_{p²}(a·α) · e_p(a·(p-1)α²/2)`.

This is a polynomial in `α` of degree 2 (after the leading constant `e_{p⁴}(a)` factors out). **Quadratic, not cubic** in the principal-unit-coset coordinate.

Cochrane Theorem 2 applies to polynomial phases. Closed-form magnitude: `|H(a, q)| = p^{(r-1)/2} · m` where m counts the saddle points of the quadratic.

For a quadratic phase, saddle is unique (the single critical point of a parabola). So `m = 1`, giving `|H(a, q)| = √N` magnitude on the principal-unit coset — saturated √N from a length-N coset.

**Verdict:** **STEP1_WORKS_MODIFIED** — phase is QUADRATIC not CUBIC on principal-unit coset. Cochrane closed form directly applies.

But: the quadratic structure changes the downstream story significantly. R78's cubic-with-r=3-saddle-exactness was the special alignment; for quadratic phase, saddle exactness happens at every r ≥ 2 (trivially, since the polynomial degree is 2 and ⌈r/2⌉ ≥ 1 for r ≥ 2).

---

## STEP 2 — bijection between parameter and saddle index

**R78 analog:** a ↔ C_a bijection.

**Candidate test:** Quadratic phase in α: `a · ((p-1)/2 · α² + α + ...) · scaling`. Saddle is `α* = −1/((p-1) · scaling stuff)` — single value, independent of `a` to leading order (since the quadratic completion eliminates the `a` dependence of the saddle position).

**Verdict:** **STEP2_FAILS** — there is NO non-trivial bijection between the parameter `a` and a saddle position, because the quadratic phase has a single fixed saddle (independent of `a` to leading order). The "p-fold bijection on the principal-unit coset" central to R78 doesn't have an analog here.

The chain stalls at Step 2 for quadratic phases. The bijection in R78 is fundamentally CUBIC-PHASE-SPECIFIC: it works because the cubic phase's `dP_a/ds = 0` saddle has multiple roots parametrized by the leading p-adic digit of the parameter.

---

## STEP 3 — saddle exactness

**R78 analog:** Saddle exact at r = 3.

**Candidate test:** Quadratic phase: saddle exactness is trivial — quadratic completion gives `Σ e_q(aα² + bα + c) = e_q(c − b²/(4a)) · √q · ε(a)` (classical quadratic Gauss sum). Exact.

**Verdict:** **STEP3_WORKS** trivially for quadratic phase — quadratic Gauss sum closed form.

But: this is a "different mechanism" than R78's saddle-point closed form. The quadratic Gauss sum is classical; no Cochrane-Pinner machinery needed. The R78 chain's Step 3 is OVER-POWERED for this case.

---

## STEP 4 — Inner-Plancherel reduction

**R78 analog:** Linear-in-c_2 phase enables Plancherel collapse.

**Candidate test:** Since Step 2 fails (no bijection to saddle parameter), there's nothing to Plancherel-collapse against. The quadratic Gauss sum is closed-form; no inner sum remains.

**Verdict:** **STEP4_FAILS** — the inner sum is not present for quadratic phase; the chain doesn't run.

---

## STEP 5 — 1/sin grid identity

**R78 analog:** Σ csc grid bound.

**Candidate test:** Not applicable — outer sum is already closed-form from quadratic Gauss sum.

**Verdict:** N/A (inherited from Step 4 failure).

---

## Five-step summary for Candidate 3

| Step | R78 (reference) | Candidate 3 | Outcome |
|------|-----------------|-------------|---------|
| 1 | F̂_p closed form | Cochrane Theorem 2 on quadratic phase | **WORKS_MODIFIED** (phase is quadratic) |
| 2 | p-fold bijection | No non-trivial bijection (single saddle) | **FAILS** |
| 3 | Saddle exact via Cochrane | Quadratic Gauss sum (different mechanism) | **WORKS** by classical machinery, not R78 chain |
| 4 | Inner-Plancherel | No inner sum to collapse | **FAILS** |
| 5 | 1/sin grid identity | N/A | N/A |

**Overall:** The chain FAILS at Step 2 for this candidate. The Heilbronn-on-coset object's phase, after Postnikov substitution, becomes QUADRATIC not CUBIC, and the multi-saddle bijection that drives R78 doesn't exist.

**However:** the candidate ADMITS a different, simpler closed form (classical quadratic Gauss sum) that gives saturation `√N`. So the *result* (√N saturation) generalizes by a different route. The *chain* doesn't generalize.

**Verdict:** STEP2 FAILS at a specific load-bearing point: the bijection requires the saddle to be parametrized by the leading digit of `a`, which requires cubic-or-higher phase degree. Quadratic phases have a single saddle and the chain doesn't run.

This is a clean negative result. The R78 chain depends on **cubic-or-higher polynomial phase structure**, which is one specific feature among many possible phase shapes in the analytic-number-theory literature. Heilbronn-on-coset has the WRONG phase shape for the chain.

**Anti-cherry-pick note:** I included this candidate specifically because its phase shape is QUADRATIC, anticipating it would fail. The failure mode is informative: it tells us exactly which structural feature of R78 (cubic phase degree) is load-bearing.
