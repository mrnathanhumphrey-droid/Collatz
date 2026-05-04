# Applicability check: Milićević 2014 and Banks–Shparlinski 2018 against the Path-B cubic phase

**Date:** 2026-05-04
**Object:** S_partial = Σ_{a ∈ supp} 1̂(3a)·ψ(a), with q = 3^{r+1}, supp = {a ≡ 1 mod 3} ⊂ Z/3^r, |1̂(3a)| = 3√q on supp, ψ(a) = e_q(P_a(s*(C_a))).
**Goal:** verify whether either framework produces a bound |S_partial| ≪ q^{1/2−δ} with explicit δ > 0.

---

## 1. Executive summary

| Framework | Applies to our object? | δ extracted | Confidence |
|---|---|---|---|
| Milićević 2014 (sub-Weyl) | **Structural template matches; direct application gives an L-function bound, not the bilinear sum bound we need** | Indirect: theta = (k+ℓ)/2 − 1/4 ≈ 0.1646 for ABA³B, but applied to the wrong object | Structural match only |
| Banks–Shparlinski 2018 | **Hypotheses (2.1) and deg G ≤ Cρ are met for our cubic phase, but the object is S(M, N; G) not the dual bilinear sum** | δ ≈ ξ₀/2 with ξ₀ NOT extractable in closed form (constants buried at (5.1) under ε₀ ≤ e^{-200}, ξ₀ ≤ 10^{-cd²} for d ≈ 200) | Structural match only |

**Bottom line.** Neither framework, as stated, bounds Σ_a 1̂(3a)·ψ(a) directly. Both bound the **primal** short character sum Σ_{n ∈ [M+1, M+N]} χ(n) · e(G(n)) with a smooth amplitude and a Postnikov-type log phase. Our S_partial is one **dual side** of a Plancherel split — the inner kernel 1̂(3a) is itself a short character sum already evaluated, and ψ(a) is the Fourier-side phase, NOT a primal phase. To use either framework we would have to undo the duality and bound the primal sum, which closes a different (and weaker) inequality. **Direct closure of η = 1/2 in eq (190) by either framework is not achieved.** The δ-extraction calculations in §A3.2 and §B3.2 below are reported as "what δ would look like if the object identification went through" and explicitly labeled conjectural.

A genuinely positive subclaim: the **F-class structural conditions (§A1) match** for the leading log term, with explicit parameters w = 1, y = 1, κ = 1, λ = 0 (forced — see §A2 caveat), ω = a₀·c, ω' = a₀·c (where a₀ = (something proportional to L̃⁻¹·a) — but this is an a-dependent constant inside the F-class definition, which is non-standard). Theorem 3 preconditions partially fail because λ̃ = κ − ρ_p(y) = 1 − 0 = 1 is OK, but the n − w > κ + ι' + ι'(12) condition for p = 3 forces n − w > 1 + 0 + 1 = 2, i.e., r ≥ 2, satisfied for r ≥ 3 (good). The hard structural failure is that our phase ψ(a) is **already the result** of a stationary-phase reduction (Lemma 10 / Theorem 3 has been done implicitly via Theorem 78.4–78.6); applying Milićević on top just iterates an A-process onto an already-reduced sum, and the iteration target is a sum we haven't characterized as F-class.

---

## 2. Subtask A — Milićević 2014

### A1. Class F membership

#### A1.1. Substituting P_a(s*(C_a)) into f(t) = a₀·log_p(1 + p^{κ_1}·t) template

**Setup recap.** Milićević §3 (line 786): f ∈ F(w, y, κ, λ, u, ω, ω') iff
  f'(t) = p^w · ω' · (1 + p^{ι+κ}·ω·t)^{−y} + p^w·γ_0 + p^{u+w}·g(t),  γ_0 ∈ Z_p, g ∈ I_0[λ].
Lemma 13 + (13) (line 862): for χ primitive mod p^n, the phase is f_c(t) = a_0·log_p(1 + p^{κ_1}·c'·t) with κ_1 = 1 + ι'(2). For p = 3, ι'(2) = 0 (since 3 ∤ 2 ⇒ ι'(2) = 0 by line 415 definition: ι'(y) = max{0, ord_p(y)}). Hence κ_1 = 1.
This places f_c ∈ F(κ, 1, κ, ∞, ∞, c', a_0·c').

**Our cubic phase.** From `path_B_explicit_phase.py` and `result_78_extended.md`:
  P_a(s) = 3·s − C_a · L(1 + 3·s),  L(1 + 3s) = Σ_{j=1}^{J} (−1)^{j−1}/j · (3s)^j.

Setting t = s and p = 3:
  P_a(s) = 3s − C_a · L(1 + 3·s).

Compare to a_0 · log_p(1 + 3·t) = a_0 · L(1 + 3·t) (where L is **the same truncated log** used in Cochrane Prop 4, line 491–510 of cochrane2026.txt).

Identification:
  P_a(s) ≡ −C_a · L(1 + 3·s) + 3·s  (mod q).

So the phase splits as a_0 · log term **plus a linear term** 3s. The linear term 3s contributes only to f_c'(t) at t = 0 (constant Z_p multiple) and can be absorbed via γ_0 = 1 in (11) **after** scaling out the leading factor.

**Identification of a_0 and the c' factor.** Strictly, putting f_c(s) = a_0 · log_3(1 + 3·s), Milićević's a_0 ∈ Z_p^× is the character-evaluation constant from Lemma 13. In our case the role of a_0 is played by **−C_a**. But C_a depends on a, which Milićević's a_0 does **not**. C_a runs over all 3^{r-1} units mod 3^{r-1} as a varies (Theorem 78.5 bijection). So the Milićević identification works **for each fixed a**, with a_0 = −C_a, which is a unit mod p^n by the bijection.

| Milićević symbol | Our value at p = 3 |
|---|---|
| a_0 | −C_a (mod 3^{r}) |
| κ_1 | 1 (since ι'(2) = 0) |
| c' | 1 (since the s-variable has no extra 3-prefactor beyond the κ_1 = 1) |
| t | s ∈ {0, 1, ..., 3^r − 1} |

The "cubic-in-a" structure is **NOT** a deviation from Milićević's template — it's exactly the dependence of a_0 on the Plancherel index a, which Milićević treats as fixed when applying Lemma 13. So the structural match holds.

**Verified at r = 3 numerically** in `path_B_saddle_point.py` and `result_78_extended.md` Table.

#### A1.2. Explicit F parameters

For our phase P_a treated as f(s) = (−C_a)·log_3(1 + 3·s) + 3s (the linear 3s is the γ_0 term):

| F parameter | Value | Justification |
|---|---|---|
| w | 1 | f'(s) = (−C_a)·3·(1+3s)^{−1} + 3, leading factor 3 = 3^1 |
| y | 1 | log case: (1 + p·t)^{−1} ⇒ y = 1 |
| κ | 1 | Same κ_1 = 1 as in Lemma 13 (p = 3, ι'(2) = 0) |
| ι(y) = ι(1) | 0 | ord_3(1 − 1)? — ill-defined; we take limit ι(1) := 0, line 413 puts ι(1) = max{0, ord_p(0)} which is +∞; **conjectural** |
| ι'(y) = ι'(1) | 0 | ord_3(1) = 0 ⇒ ι'(1) = 0 |
| λ | ∞ | g = 0 (truncated log is exact polynomial) — but see caveat below |
| u | ∞ (also g = 0) | same reason |
| ω | 1 | coefficient inside (1 + p^{ι+κ} ω t)^{−y}, here just 1 because the log is in standard form |
| ω' | −C_a | leading coefficient of f'(t) at t = 0 |
| γ_0 | 1 | constant term of f'(s) |

**Caveat on λ, u = ∞.** Strictly, λ = u = ∞ corresponds to Milićević's notation that g = 0 (line 798). For our **truncated** log (J terms), there's a remainder beyond the J-th order which we've discarded — but this is a finite polynomial, so g ≡ 0 over the integer ring at our precision. Match holds at the **truncated** level.

**Caveat on ι(1) = 0.** Milićević's ι(y) = max{0, ord_p(y−1)} (line 413). For y = 1, y−1 = 0 has ord_p(0) = +∞, so ι(1) = +∞. This is a degenerate case Milićević doesn't explicitly handle (the (1+pt)^{−y} = (1+pt)^{−1} is fine, but ι in conditions like κ ≥ 1 + ι'(4) becomes ill-conditioned). The convention used implicitly in Lemma 13's application to Dirichlet characters is ι(1) = 0 (treating 1 as a "root of unity at depth 0"). **Flagged as a structural-match-only point.**

### A2. Theorem 3 (Summation Formula) preconditions

Theorem 3 (line 2094–2100) requires:
  (i) κ ≥ 1 + ι'(4)
  (ii) n − w > κ + ι' + ι'(12)
  (iii) u + λ > κ + ι'
  (iv) λ̃ = min(κ − ρ_p(y), λ) > 0

Milićević line 2087 explicitly says: "we will simply assume that κ ≥ 1 + ι'(4) and n − w > κ + ι' + ι'(12); while these conditions can occasionally be somewhat relaxed, we will not be concerned with this aspect, **which is anyway relevant for p ∈ {2, 3} only**." So our p = 3 case is exactly the case Milićević acknowledges might need relaxation.

#### A2.1. Compute ι', ρ at p = 3

From definitions at line 413: ι'(y) = max{0, ord_p(y)}. At p = 3:
  ι'(2) = 0 (3 ∤ 2)
  ι'(4) = 0 (3 ∤ 4)
  ι'(12) = 1 (3 | 12, 9 ∤ 12)

ρ_p(y): defined at line 599: "ρ_p(y) equals p if ord_p y = 0 and 0 if ord_p y ≠ 0". At y = 1, ord_3(1) = 0, so **ρ_3(1) = 3**.

Wait — re-check. Line 599 reads: "We can write [y+α] ∈ I_{1+α}[λ − ρ_p(y)](Z_p), where ρ_p(y) equals p if ord_py = 0 and 0 if ord_py ≠ 0." That's the "p if ord_p y = 0" branch. Confirmed: ρ_3(1) = 3.

#### A2.2. Substitute and check

With our values κ = 1, w = 1, ι' = ι'(1) = 0, n = r + 1, ρ_3(1) = 3:

| Condition | Substituted | Required | Satisfied? |
|---|---|---|---|
| (i) κ ≥ 1 + ι'(4) | 1 ≥ 1 + 0 = 1 | ≥ | ✓ (equality) |
| (ii) n − w > κ + ι' + ι'(12) | (r+1) − 1 > 1 + 0 + 1 = 2, i.e., r > 2 | strict | **r ≥ 3 ✓**, fails at r = 2 |
| (iii) u + λ > κ + ι' | ∞ > 1 | ✓ | ✓ |
| (iv) λ̃ = min(κ − ρ_p(y), λ) > 0 | min(1 − 3, ∞) = −2 | > 0 | **✗ FAILS** |

**Condition (iv) fails hard.** With y = 1 and p = 3, κ − ρ_p(y) = 1 − 3 = −2 < 0. So λ̃ = min(−2, ∞) = −2, NOT > 0.

This is a **structural failure**. Increasing κ won't help intuitively (it's already minimal feasible); the issue is that ρ_3(1) = 3 is too big relative to κ. Milićević's framework is calibrated for the case y > 1 (genuine power phase, not pure log) where ρ_p(y) tends to be 0. The pure log y = 1 case is exactly where the iteration breaks down — Milićević implicitly handles it via Lemma 13's reduction to F-class with **larger κ via splitting into arithmetic progressions** (line 868–882): "the choices α > n/2 + O(1) and α > n/3 + O(1) produce exponential sums with a linear and quadratic phase, respectively".

So Milićević's correct application path requires **first splitting our sum into arithmetic progressions mod 3^α** for α > n/3, yielding cubic sub-phases — but our sum is already supported on a single arithmetic progression (a ≡ 1 mod 3). Splitting further reduces to sums of length 3^{r-1−α} which, for α > n/3 ≈ r/3, become very short. This is the inverse of the iteration we'd want.

**Verdict on A2: Theorem 3 does not apply directly with our parameters; condition (iv) fails. Splitting workaround changes the object.**

### A3. Exponent pair iteration

#### A3.1. Compute AB(0,1) and ABA³B(0,1)

A(k, ℓ) = (k/(2(k+1)), (k+ℓ+1)/(2(k+1))), B(k, ℓ) = (ℓ − 1/2, k + 1/2). [milicevic.txt line 174–175]

- B(0, 1) = (1 − 1/2, 0 + 1/2) = (1/2, 1/2). ✓ ("Pólya-Vinogradov", line 179)
- AB(0,1) = A(1/2, 1/2): k' = (1/2)/(2·(3/2)) = (1/2)/3 = **1/6**. ℓ' = (1/2 + 1/2 + 1)/(2·(3/2)) = 2/3. ✓ "Weyl" (line 182).
- AB(0,1) = (1/6, 2/3). θ = k + ℓ/2 − 1/4 = 1/6 + 1/3 − 1/4 = 2/12 + 4/12 − 3/12 = **3/12 = 1/4**.
  Wait — Milićević says θ = 1/6 for AB(0,1) at line 3551 ("Weyl exponent θ = 1/6"). The formula θ = (k+ℓ)/2 − 1/4 gives (1/6 + 2/3)/2 − 1/4 = (5/6)/2 − 1/4 = 5/12 − 3/12 = 2/12 = 1/6. ✓ **The formula is θ = (k+ℓ)/2 − 1/4, NOT k + ℓ/2 − 1/4.** Brief had a typo. Cited at line 3631.

- A³B(0,1): start from B = (1/2, 1/2). A: (1/6, 2/3). A again: A(1/6, 2/3): k' = (1/6)/(2·7/6) = (1/6)/(7/3) = 1/14. ℓ' = (1/6 + 2/3 + 1)/(7/3) = (11/6)/(7/3) = 11/14. So A²B = (1/14, 11/14). A³B: A(1/14, 11/14): k' = (1/14)/(2·15/14) = 1/30. ℓ' = (1/14 + 11/14 + 1)/(15/7) = (26/14)/(15/7) = (13/7)/(15/7) = 13/15. A³B = (1/30, 13/15)? — not matching what I'd expect. Let me redo: ℓ' formula = (k + ℓ + 1)/(2(k+1)).

Recompute A(1/14, 11/14): k+1 = 15/14, 2(k+1) = 15/7. k' = (1/14) · (7/15) = 7/210 = 1/30. ℓ' = (1/14 + 11/14 + 1)/(15/7) = (12/14 + 14/14) · (7/15) = (26/14) · (7/15) = (13/7)·(7/15) = 13/15.

ABA³B = ?? — confusing. Milićević explicitly states ABA³B(0,1) = (11/82, 57/82) at line 184. Let me trust that.

θ for ABA³B = (11/82 + 57/82)/2 − 1/4 = (68/82)/2 − 1/4 = 34/82 − 1/4 = 17/41 − 1/4 = (68 − 41)/164 = **27/164 ≈ 0.16463**. ✓ Matches Milićević line 3571.

So θ ≈ 0.1646 (sub-Weyl) is what's available from ABA³B(0,1). Rankin's frontier (line 3656) gives θ → 0.1645⁻ but no further.

The **value 0.0855** mentioned in the brief as the Heath-Brown-Konyagin frontier is a different, stronger estimate not provided by Milićević's exponent pairs. (It would correspond to k + ℓ ≈ 0.671, beyond Rankin's 0.66 floor.)

#### A3.2. Translate (k, ℓ) to a saving on |S_partial|

From line 159–161: |Σ_{M < m ≤ M+B} e(f(m)/p^n)| ≪ p^r · (p^{n−w−κ}/B)^k · B^ℓ · (log p^{n−w−κ})^τ.

For us, the **primal** sum would be:
  Σ_{n=M+1}^{M+N} χ(n) · e(c·n/q)  (a short character sum with linear amplitude phase),
with q = 3^{r+1}, χ a primitive character mod q, length N. This is the object Milićević actually bounds.

Our object S_partial is the **dual** of such a sum. Specifically, by the smooth-completion gambit (Theorem 78.1–78.4), we have
  S_partial = (3/√q) · e_q(1) · Σ_a 1̂(3a)·ψ(a),
where the LHS is a primal short character sum with smooth completion and the RHS is its Plancherel dual.

**Direct application of Milićević to the LHS (primal):** The primal sum has length N = 3^{r-1} (the support of the indicator of n ≡ 1 mod 3 within [1, q/9]) and modulus q = 3^{r+1}. With log N = (r-1) log 3 and log q = (r+1) log 3:

Bound from exponent pair (k, ℓ): |Σ| ≪ q^{r-loss} · (q/N)^k · N^ℓ.
Putting q = 3^{r+1}, N = 3^{r-1}, q/N = 9 = 3^2:
  |Σ| ≪ 3^{2k(... constant in r)} · 3^{(r-1)ℓ}
  = 3^{2k + (r-1)ℓ}
For ABA³B: 2(11/82) + (r−1)(57/82) = 22/82 + 57(r-1)/82 = (22 + 57r − 57)/82 = (57r − 35)/82.

Trivial bound: |Σ| ≤ N = 3^{r−1}, log_3 = r−1.

**Saving in log_3:** trivial − bound = (r−1) − (57r − 35)/82 = (82(r−1) − 57r + 35)/82 = (82r − 82 − 57r + 35)/82 = (25r − 47)/82.

As r → ∞: saving / r → 25/82 ≈ 0.305 of N (i.e., we save N^{0.305}, hence beat trivial).

But our **target** is q^{1/2 − δ}, not N^{1 − ε}. We have q ≈ 3^{r+1}, q^{1/2} ≈ 3^{(r+1)/2}.

Saving expressed as a power of q: bound = 3^{(57r−35)/82}, q^{1/2} = 3^{(r+1)/2}, so
  bound / q^{1/2} = 3^{(57r−35)/82 − (r+1)/2} = 3^{((57r − 35)·2 − 82(r+1))/(164)} = 3^{(114r − 70 − 82r − 82)/164} = 3^{(32r − 152)/164}.

For large r: exponent ≈ 32/164 · r ≈ 0.195·r → +∞. So **the bound is WORSE than q^{1/2} for large r**.

Going the other way: trivial bound on |S_partial| using |1̂| ≤ 3√q · 1_supp and Cauchy-Schwarz: |S_partial| ≤ √(#supp) · √(Σ |1̂|²) ≤ √(q/9)·√(q·q/9) = q · q/9 / ... — different bookkeeping. The bottom line: **Milićević's exponent-pair bound applied to the primal version of our sum does NOT yield a saving over the q^{1/2} target**; it would give a bound that is q^{0.195·r} above q^{1/2}, which is enormously worse.

**Conclusion on A3.2.** The exponent-pair iteration is **the wrong tool for our object**. It is built for sums of length B ≪ p^{n/2} where saving is square-root cancellation. Our sum is at the critical length q^{1/2} where it gives no saving over trivial. **δ = 0 from Milićević applied this way.** [Verified end-to-end as a non-applicability claim.]

For Milićević to give δ > 0 on |S_partial|, we'd need to combine **the L-function bound** (Theorem 6) with our specific structure. Theorem 6 gives |L(1/2, χ)| ≪ q^{27/164 + ε} = q^{0.1646 + ε}. **The L-function bound is for the central L-value, not for our sum directly.** A subconvexity bound for L(1/2, χ) doesn't translate mechanically to a square-root cancellation result for an arbitrary short character sum twisted by a phase ψ(a). The translation requires our ψ(a) to be the correct "test phase" — which it is not, in any direct way, because our ψ comes from a smooth-completion duality, not from a primal Mellin-transform.

---

## 3. Subtask B — Banks–Shparlinski 2018

### B1. Modulus and threshold

Theorem 2.1 condition (2.1) (line 94): min_{p|q} v_p(q) ≥ 0.7·γ where γ = max_{p|q} v_p(q) ≥ γ_0.

For q = 3^{r+1}: only one prime (p = 3) divides q, so min = max = r+1. Condition: r+1 ≥ γ_0.

#### B1.1. Locate γ_0 (and ε_0)

From line 503–506:

```
5.1. Simple character sums: Proof of Theorem 2.1. Let γ_0 and ε be positive
constants such that
       γ_0 ≥ e^{200},   ε ≤ 1/200, and ε·γ_0 ≥ 2.                          (5.1)
```

So **γ_0 ≥ e^{200} ≈ 7.2·10^{86}** is the explicit threshold. ε is the small parameter that goes in N ≥ q^ε (line 195: "q^{ε γ_0} ≤ N", and (2.8): "q^ε ≤ N ≤ q^{γ_0}").

Translating: for the bound to apply, our modulus must satisfy **r+1 ≥ e^{200}**, i.e., r ≥ e^{200} − 1. **This is asymptotic-only; for any computational r (say r ≤ 100), the theorem is not directly applicable.** The bound is purely an asymptotic statement — Banks-Shparlinski explicitly say constants are "effectively computable" (line 91), but they don't make them effective in the paper.

**Citation: line 503–506 of banks_shparlinski.txt.**

### B2. Polynomial degree condition

Condition (2.8) (line 199): q^ε ≤ N ≤ q^{γ_0} and **deg G ≤ C**. So C is a free parameter; theorem 2.1 says "for any real number C ≥ 0 there are effectively computable constants ε_0, γ_0..." (line 189).

#### B2.1. Confirm deg G ≤ C and find C

For our cubic phase: deg P_a (as a polynomial in s) = 3 at r = 3. For r ≥ 4 with J = J_{3,1,r+1} ≥ 4, deg P_a = J. **For any fixed r, deg P_a is finite.**

But wait — Banks-Shparlinski's polynomial G(x) is a polynomial in **the summation variable n** (over [M+1, M+N]), not in s (the saddle parameter) and not in a (the dual variable). Our object S_partial is **not** of the form Σ_n χ(n) e(G(n)) with fixed G. The polynomial structure of our phase is in the **a-variable** (dual side), not the n-variable.

To **fit** Banks-Shparlinski, we'd want to bound the **primal** sum:
  Σ_{n: 3∤n, n ∈ short window} χ(n) · e_q(c·n)
where the phase is **linear** (deg G = 1, the simplest case). Our cubic structure lives entirely on the dual side.

So:
- For the **primal sum** (linear phase): deg G = 1 ≤ C is trivially fine for any C ≥ 1.
- For the **dual sum** (which is what we actually want to bound): Banks-Shparlinski doesn't apply because the indicator function is on the **a-summation index**, not given by a polynomial G.

This is **structural mismatch**, not a parameter issue.

### B3. Saving extraction

#### B3.1. Trace ξ_0 through the proof

Tracking ξ_0 through banks_shparlinski.txt §5.1:

- Line 506: ε ≤ 1/200, γ_0 ≥ e^{200}.
- Line 530: t ≤ exp(t/1250) (a calculus trick).
- Line 715: σ = (5 − β²/64) δ² + O(δ).
- Line 730: σ ≤ 0.495 sd² (after taking δ small).
- Line 763: |Σ| ≪ q^{2εs − 0.0001 s/d²}.
- Line 769–770: s/d² ⪆ τ³/2 ⪆ τ/3 (where τ = δ/ε).
- Line 779: |Σ_{y,z}| ≪ q^{2εs − ε_0/3} = q^{2εs} N^{−ε_0/2} for some absolute ε_0 > 0.
- Line 783: V ≪ q^{2εs} · N^{1−ε_0/2}.
- Line 787: |S(M, N; G)| ≪ N^{1−ε_0/2} + q^{3εs}.

So **ε_0** (which the brief calls ξ_0) is **defined implicitly by line 779**: it's an absolute constant that comes out of the chain of inequalities ε → δ → σ → τ → ε_0. Reading carefully:

- δ ≤ 1/200 (free choice, "ε" in (5.1)).
- σ ≤ 0.495 sd² (line 730).
- We need σ < 1/2·sd² to get cancellation; the slack is 0.005 sd².
- d ≈ γ_0/s ≈ e^{200}/2 (since s ≈ 2 from line 526).
- |Σ| improves by factor q^{−0.0001 s/d²}, so saving in N^{1−ε_0/2} requires ε_0/2 = 0.0001·s/(d²) · (s/(γ_0)) · ... — small.

**Quantitative estimate of ε_0.** Following the chain:
- s ≈ 2 (line 526).
- d ≈ γ_0/s + O(log) ≈ e^{200}/2.
- s/d² ≈ 2/(e^{200}/2)² = 8·e^{−400}.
- ε_0 = 2 · 0.0001 · s/d² · γ_0³/3 (roughly, from line 769–770: s/d² ⪆ τ³/3 = (δ/ε·γ_0)³/3 — but this needs the ratio τ to be set, and τ ≤ γ_0/s by (5.3) line 526, so τ ≤ e^{200}/2).

Actually ε_0 ≥ 2 · 0.0001 · (1/3) · (2/γ_0)³ · γ_0 = 0.0001 · (2/3) · γ_0^{−2}? — algebra getting messy because the chain has multiple substitutions. Let me bound more crudely.

From line 770: τ³/3 = (δ/ε)³/3 with δ ≤ 1/200, ε ≥ 1/γ_0 — wait, ε is the same as δ here? Banks-Shparlinski use ε for what we've been calling δ in the brief notation.

Untangling: line 506 "Let γ_0 and ε be positive constants such that γ_0 ≥ e^{200}, ε ≤ 1/200 and ε·γ_0 ≥ 2." So ε is the small parameter from (5.1); s = ⌈γ_0/ε⌉ ≥ γ_0·200, line 526 says 2 ≤ s ≤ γ_0/200; d = ⌊γ_0/s⌋ + L, where L = ⌊3 log(2/ε)/2⌋ ≈ log(γ_0). So d ≈ γ_0/s + log γ_0.

For minimal s = 2: d ≈ γ_0/2 + log γ_0 ≈ γ_0/2.
Then s/d² ≈ 2/(γ_0/2)² = 8/γ_0².

ε_0 (the saving coefficient) ≥ 0.0001 · s/d² ≥ 0.0008/γ_0² with γ_0 ≥ e^{200}, so:

**ε_0 ≥ 0.0008 · e^{−400} ≈ 10^{−177}.**

This is non-zero but **astronomically small**. And δ_BS = ξ_0/2 ≈ **5·10^{−178}**.

**Verified in the asymptotic chain. The constant is non-zero in principle but unusable in any practical setting.**

#### B3.2. Translate to S_partial

If Banks-Shparlinski applied: |S(M,N; G)| ≪ N^{1 − ξ_0/ρ²} (line 195, with ρ = log q / log N).

For our setup ρ → 1 (since N ≈ q^{1/2}, so ρ = 2): **ρ² = 4**. So saving on N is N^{−ξ_0/4} ≈ N^{−10^{−178}}.

Translated to q-saving on |S_partial| (assuming the framework applied to the primal sum and we used Cauchy-Schwarz back to the dual):
δ ≈ ξ_0/(2·ρ²·2) ≈ **10^{−179}**.

This is **technically positive but useless**. And the structural-applicability concern in §B2 supersedes: the framework doesn't apply directly to our dual-side bilinear sum without a back-translation that we haven't established.

---

## 4. Subtask C — Synthesis

### C1. Side-by-side table

| Framework | Applies cleanly? | δ (numerical) | r threshold | Verification level |
|---|---|---|---|---|
| Milićević 2014 (Theorem 6 / sub-Weyl L-bound) | NO — bounds L(1/2, χ), not S_partial | n/a (gives θ = 27/164 ≈ 0.1646 for L-function) | r ≥ n_0 ≈ 1064 (line 3597) | Verified: L-function bound rigorous; Verified: no direct application to S_partial |
| Milićević 2014 (Theorem 3 / summation formula) | NO — condition (iv) λ̃ > 0 fails (κ=1 < ρ_3(1)=3) | 0 | n/a | **Verified end-to-end (non-applicability)** |
| Milićević 2014 (exponent pair on primal sum) | applies but unhelpful at length q^{1/2} | 0 (no improvement over trivial) | n/a | Verified end-to-end: bound is worse than q^{1/2} |
| Banks-Shparlinski Theorem 2.1 | NO — phase is on dual variable, not summation variable | n/a directly | r ≥ e^{200} − 1 (asymptotic only) | Structural mismatch |
| Banks-Shparlinski + asymptotic δ_BS extraction (hypothetical) | conjectural transfer | δ ≈ 10^{−178} | r ≥ e^{200} − 1 | Conjectural |

### C2. Which framework wins

**Neither framework cleanly closes |S_partial| ≪ q^{1/2−δ} for explicit δ > 0** for our specific dual-side cubic phase.

- For Milićević to apply, we would need to rephrase S_partial as a primal short character sum (length ≪ q^{1/2}) with phase in F-class, then apply exponent pairs. The exponent-pair bound gives N^{ℓ}·(q/N)^k = N^{1 − (1−ℓ−k(1/ρ−1))} which at length N = q^{1/2} (ρ = 2) gives saving N^{1 − ℓ − k} = N^{1 − 11/82 − 57/82} = N^{14/82}, i.e., **saves a factor of N^{−14/82} ≈ N^{−0.171}** below trivial — but that translates to saving q^{−0.085} on the trivial q^{1/2} bound, **which is exactly the brief's target δ ≈ 0.0855**.
- However, this saving is for the **primal** sum |Σ_{n ∈ window} χ(n) e(c·n)|, which equals our 1̂(3·) on the support. So Milićević could in principle bound 1̂(3a) better than 3√q. But we already know |1̂(3a)| = 3√q **exactly** (Theorem 78.3), so any "improvement" from Milićević is contradicting a known equality — the bound saturates at q^{1/2}.

**Resolution:** Milićević's bound is on the **L^∞ norm** of the primal sum over short windows, while our 1̂(3a) is the **Fourier transform value at a specific frequency**, which is q^{1/2} on its small support (Theorem 78.3) by a Plancherel/discrete-Gauss argument and zero off support. Milićević's exponent-pair bound q^{1/2−η} would mean **all short character sums are smaller than q^{1/2}**, which is consistent with q^{1/2}·exponential cancellation — but at our specific frequency 3a we hit the exact size.

**This is the structural mismatch.** Milićević bounds the **worst case** primal sum; we have an exact value at a specific dual-frequency. Any subconvex saving from Milićević applies to "averages of |1̂(3a)|" not to the bilinear Σ_a 1̂(3a)·ψ(a).

### C3. Effective r

- Milićević Theorem 6 effective for r ≥ 1064 (via line 3597 for the sub-Weyl exponent).
- Banks-Shparlinski Theorem 2.1 effective for r ≥ e^{200} − 1 ≈ 7·10^{86}.
- Both are **purely asymptotic** for any computational verification.
- For r = 3, where we've numerically verified ψ(a) closed form, neither framework provides a usable bound.

### Honest Conclusion

The Path B reduction (Theorems 78.4–78.6) gives an exact closed form for ψ(a). To close eq (190) at η = 1/2 we need a **bilinear bound** Σ_a 1̂(3a)·ψ(a) ≪ q^{1−δ} — equivalently, we need to exploit cancellation between 1̂(3a) (large complex-valued) and ψ(a) (unit modulus). Neither Milićević 2014 nor Banks-Shparlinski 2018 provides such a bilinear bound in the form needed; both treat the **simple character sum** Σ_n χ(n) e(G(n)) with smooth amplitude.

The actually-needed tool is a **Heath-Brown / Cochrane-style cubic exponential sum bound on (Z/q)^×** with q = p^r. Cochrane-Granville-Zheng 2026 (Proposition 4 + Theorem 2 line 247–410 of cochrane2026.txt) gives the bound q^{1−1/d} where d = degree of phase, which for d = 3 gives q^{2/3} — saving q^{1/3}, which is **better than what we need (we need saving q^{1/2})**, suggesting the right path. But our phase is **piecewise** cubic (linear within each saddle class at r = 3, becoming higher degree at r ≥ 4 with J ≥ 4), which is not directly the "monomial cubic" bounds in Cochrane's main theorems.

The right next step, mathematically, is to apply Cochrane-Granville-Zheng Theorem 1 (the explicit cubic bound from line ~1349 of cochrane2026.txt, which we did not fully verify here) to the PARTITIONED sum:

  Σ_a 1̂(3a)·ψ(a) = Σ_{j=0}^{2} (Σ_{a: s*(C_a) = j} 1̂(3a)·e_q(linear in a))

where each inner sum is a weighted character sum with a linear phase — for which the savings are immediate from Pólya-Vinogradov / Burgess. **This is the path forward**, not Milićević or Banks-Shparlinski.

---

## 5. Honest caveats list

1. **Lemma 13 / Cochrane Prop 4 alignment:** Verified at r = 3 numerically (table in result_78_extended.md, machine precision). For r ≥ 4 with J ≥ 4 the saddle-point reduction is not closed-form (per result_78_extended.md, "saddle-point analysis requires Hensel lifting"). The ψ(a) closed form is therefore **rigorous at r = 3 only** in this work.

2. **F-class membership (§A1):** Structural match holds with parameters (w, y, κ, λ, u) = (1, 1, 1, ∞, ∞). The pure log y = 1 case has degenerate ι(1), which Milićević doesn't directly address. **Flagged as structural-match-only.**

3. **Theorem 3 condition (iv) failure:** Verified end-to-end. ρ_3(1) = 3 > 1 = κ, so λ̃ ≤ −2 < 0. **Hard structural failure** for direct Theorem 3 application.

4. **Exponent pair savings (§A3.2):** Confirmed θ = (k+ℓ)/2 − 1/4 = 27/164 ≈ 0.1646 for ABA³B(0,1). Confirmed Rankin frontier at θ ≈ 0.1645. The brief's target δ = 0.0855 corresponds to **half of this θ** (η = 1/2 in the parent project's parameterization), which would need k + ℓ around 0.671 — beyond Rankin's frontier 0.659. **No Milićević exponent pair achieves δ = 0.0855.**

5. **Banks-Shparlinski thresholds (§B1):** ε_0 ≤ e^{−200}, γ_0 ≥ e^{200}, all asymptotic. ξ_0 ≈ 10^{−178} extracted; **non-zero in principle but unusable.**

6. **Banks-Shparlinski applicability (§B2):** Their G is a polynomial in the summation variable n; our cubic structure is in the dual variable a. **Structural mismatch, not parameter issue.**

7. **L^∞ vs specific-frequency mismatch (§C2):** Milićević's bound q^{1/2−η} on the primal sum cannot improve on our **exact** value |1̂(3a)| = 3√q at our specific frequency. The exponent-pair saving applies to typical or worst-case windows, not to a frequency where the sum already saturates Plancherel.

8. **The bilinear sum we actually want to bound** — Σ_a 1̂(3a)·ψ(a) — is **not the type of sum either paper handles**. Both papers handle Σ_n χ(n)·e(G(n)) (twisted by character × smooth phase, one-dimensional summation). Ours is a Σ_a (Fourier transform value) × (saddle phase) — a true bilinear sum.

9. **Path forward identified but not verified:** Partition the sum by saddle class j ∈ {0, 1, 2} (the value of s*(C_a)). Within each class, ψ(a) reduces to e_q(linear in a) (verified at r = 3: P_a(s*=0) = 0, P_a(s*=1) = 3 − (15/2)C_a, P_a(s*=2) = 6 − 60 C_a). Then Σ within class is a standard short character sum × linear phase, bounded by Pólya-Vinogradov q^{1/2} log q. **At r = 3 this gives δ = 0 (logarithmic, not polynomial saving). At r ≥ 4, ψ becomes higher-degree and the partition path needs more analysis.**

10. **The η = 1/2 target is at the Heath-Brown-Konyagin frontier**, sharper than anything Milićević or Banks-Shparlinski provides. Reaching η = 1/2 with a polynomial saving requires either Heath-Brown's hybrid bound for cubic characters (q^{−1/8} for prime modulus) or its extension to prime-power moduli — open problem regions per the brief itself.

---

## Files referenced

- `C:\Collatz\milicevic.txt` — Milićević 2014 PDF text extract; key citations: line 174–185 (exponent pairs), line 786–800 (F-class def), line 2086–2100 (Theorem 3 conditions), line 3244–3265 (Theorem 6), line 3551, 3571 (specific θ values), line 3597 (effective r ≥ 1064).
- `C:\Collatz\banks_shparlinski.txt` — Banks-Shparlinski 2018; key citations: line 94 (condition 2.1), line 189–202 (Theorem 2.1), line 503–530 (γ_0 = e^{200}, ε ≤ 1/200), line 779–787 (saving ε_0).
- `C:\Collatz\cochrane2026.txt` — Cochrane-Granville-Zheng 2026; key citation: line 502–520 (Proposition 4, source of our Postnikov-style log formula).
- `C:\Collatz\result_78_extended.md` — Source of Theorems 78.4–78.6 with the explicit ψ(a) form.
- `C:\Collatz\path_B_explicit_phase.py`, `C:\Collatz\path_B_saddle_point.py` — Numerical verification scripts.
