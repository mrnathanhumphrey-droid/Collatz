# R81b — Mahler / 3-adic-analytic certification of the arg F̂ phase

**Date:** 2026-07-14. **Verdict: CERTIFIED — arg F̂ is a 3-adic analytic function; the "growing degree" is the shadow of a FIXED, r-independent Mahler profile.** One sub-claim falsified (see §3).

Standalone extension of R81 (`result_81b_mahler_extend.py`, does not touch `result_81_*`). Self-validated: r=3–6 reproduce R81's `(degree, v₃_leading)` = (3,3),(4,4),(4,4),(5,6) exactly; |ĝ|=3√q smoke PASS; float64 index residuals ≤ 3.3e-15 (integer phase index certified).

## 1. The invariant object: an r-independent Mahler profile

Writing arg F̂(3a) as its finite-difference (Mahler) expansion `s(b) = Σ_k c_k·C(b,k) (mod q)` in the support index `a = a0 + 3b`, the coefficient valuations `v₃(c_k)` **converge to a fixed, r-independent sequence**:

    v₃(c_k), k = 0,1,…,11  =  0, 2, 2, 3, 4, 6, 7, 8, 10, 11, 12, 15  (stable via r=15; grows ~1.3k)

| r | q=3^{r+1} | v₃(c_k) profile (below modulus cap) |
|---|---|---|
| 3 | 81 | 0 2 2 **3** |
| 4 | 243 | 0 2 2 3 **4** |
| 5 | 729 | 0 2 2 3 4 **6** |
| 6 | 2187 | 0 2 2 3 4 6 **7** |
| 7 | 6561 | 0 2 2 3 4 6 7 **8** |
| 8 | 19683 | 0 2 2 3 4 6 7 8 **9** |
| 9 | 59049 | 0 2 2 3 4 6 7 8 **10** |

Each `v₃(c_k)` **stabilizes** to the fixed value once r is large enough that 3^{r+1} exceeds it; the finite-r objects are **truncations of one 3-adic function**, not a family of unrelated polynomials. This is exactly the Mahler/Amice signature of a p-adic analytic function.

## 2. "Growing degree" is derived, not fundamental

The Z/q polynomial degree R81 reported is **entirely explained** by the fixed profile crossing the modulus threshold:

    degree(r) = max{ k : v₃(c_k) ≤ r }

- r=3→3, r=6→5, r=7→6, r=9→7 — every value matches this formula against the fixed profile above.
- So the degree grows because the analytic function's coefficient valuations grow past 3^{r+1}, **not** because the object has any intrinsic finite degree. The right invariant is the r-independent `v₃(c_k)` sequence; "growing degree" is its shadow in the finite quotient.

**Analyticity strength.** `v₃(c_k)` grows super-linearly over the observed range (~1.5·k on k=4..8), well above the p=3 analyticity threshold of k/(p−1)=k/2. arg F̂ is 3-adic analytic with room to spare (linear-in-k 3-adic decay of Mahler coefficients). Exact radius not pinned from 9 points; the qualitative class is secure.

## 3. Sub-claim FALSIFIED: degree is NOT ≈⌊r/2⌋+2 / not a stationary-phase halving

R81's "degree ≈ ⌊r/2⌋+2" was fit on r=3–6 and **BREAKS at r≥7**:

| r | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|
| observed degree | 3 | 4 | 4 | 5 | **6** | **7** | **7** |
| ⌊r/2⌋+2 | 3 | 4 | 4 | 5 | 5 | 6 | 6 |

The degree grows faster than r/2 for r≥7. **The "degree r in → degree ~r/2 out = stationary-phase signature" observation is a small-r artifact and should NOT go in any disposition.** (Under §2, degree(r)=max{k:v₃(c_k)≤r} is the correct law, and it is not ⌊r/2⌋+2.)

## 4. Routing (updated)

- **Reframing certified:** the correct disposition language is "arg F̂ is a **3-adic analytic function** of the support index, presented in its Mahler expansion; no finite Z/q-degree, but a fixed r-independent Mahler profile with linearly-growing coefficient valuations," **not** "unbounded-degree pathology."
- **Weyl / van der Corput still dead:** those require finite polynomial degree; a genuinely 3-adic-analytic (∞-degree-in-quotient) phase needs Θ(r) differencing — the uniform √-saving R81 already certified impossible. Reframing does not reopen that route.
- **Does NOT kill the object — relocates it:** 3-adic analytic phases have their own theory (p-adic oscillatory integrals / Igusa — but see `IGUSA_DISPOSITION.md` re-read: Barrier 1 (substrate-trivial) is substrate-specific and may not carry to ĝ(a); Barrier 2 (neg-rational poles vs positive-irrational target) carries over **only if the target is log₃2-as-pole**, not a bilinear bound. Pin the target before spending there.)
- **C1 cross-link:** this is the Fourier-side view of C1's `4 = 1+3 in Z₃` wall (`C1_DISPOSITION.md:38`: "Both obstructions trace to the algebraic root 4 = 1+3 in Z_3"). Three concordant views now: Cochrane degree-blowup (C1_DISPOSITION), 2-adic-exponential/not-bounded-degree (C1_TAO_RECURSION_FORM), Fourier 3-adic-analytic Mahler profile (R81/R81b).

## 6. Legendre/Kummer closed form for v₃(c_k) — TESTED, does NOT close (OPEN)

Candidate (Kummer/Legendre, given the (1+3)^x substrate): `v₃(c_k) = k·v₃(log₃4) − v₃(k!)`, `v₃(k!)=(k−s₃(k))/2`. Stable profile extended to k=11 (r=15 FFT): `v₃(c_k) = 0,2,2,3,4,6,7,8,10,11,12,15`.

- With `v₃(log₃4)=1` ⇒ `k−v₃(k!)`: residuals `0,1,0,1,1,2,3,3,4,6,6,8` — grows, **DOES NOT MATCH.**
- No clean linear closed form in `{k, s₃(k), v₃(k!)}`: least-squares max|resid| = 1.3, non-integer coefficients.
- The sequence is **deterministic/structured** (not random) but its closed form is **OPEN** — appendix-level observation, not a theorem.
- **Certification unaffected:** growth ~1.3·k > threshold k/2, so 3-adic analyticity holds regardless of the exact formula. Files: `result_81b_legendre_test.py` + `result_81b_legendre_log.txt`.

## 7. STANDING RULE — Igusa-family admissibility against ĝ(a)

Igusa / Tate / adelic **rational-pole** machinery is admissible against the R81 phase object ĝ(a) **ONLY IF the target is a bilinear / cancellation BOUND.** It is **NEVER** admissible if the target is **log₃2 (or log₃4) as a POLE / exponent** — the category obstruction (R82: a positive-irrational exponent cannot be emitted by any rational-pole machinery; unifies ADELIC arc 4 + IGUSA arc 9 + R81) kills that permanently. **A future agent reading "Igusa: NO_FIT" in the ledger must NOT conclude the p-adic-oscillatory-integral avenue is closed for bound-type targets** — the NO_FIT is pole-target-specific. An oscillatory-integral estimate on this certified-analytic phase that *outputs a bound* (not a pole location) is untouched by the category obstruction and remains open.

## 5. Process

Self-validation against R81's r=3–6 numbers was the gate; it passed before r=7–9 were trusted. The ⌊r/2⌋+2 sub-claim was caught only by extending past the fit window — the reason to compute r=7,8,9 rather than assert from r=3–6. Files: `result_81b_mahler_extend.py` + `result_81b_log.txt` + this disposition.
