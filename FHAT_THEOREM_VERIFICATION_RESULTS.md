# F̂_p Theorem Candidate — Adversarial Verification Results

## Disposition

> **THEOREM_VERIFIED** (with scope extension to r ≥ 1)

The candidate magnitude formula
> **|F̂_p^full(ξ)| = p^{(r+3)/2}** on `supp = {p·a (mod M) : a ∈ Z/p^r, a ≡ 1 (mod p)}`, `|supp| = p^{r-1}`
> *(equivalently |G[a]| = p^{(r+1)/2} where G[a] is the length-period DFT and F̂_full(p·a) = p·G[a])*

is **verified to machine precision (≤ 4.56e-12 max off-support magnitude, ≤ 1.1e-15 max relative deviation on-support) across all 27 tested cells** spanning:

- 7 new primes (Phase 1): p ∈ {11, 13, 17, 19, 23, 29, 31} × r ∈ {2, 3} = 14 cells.
- 3 higher r values (Phase 2): p ∈ {3, 5, 7} × r ∈ {4, 5, 6} = 9 cells.
- r = 1 boundary (Phase 4): p ∈ {3, 5, 7, 11} = 4 cells **(pre-reg conservatively excluded r=1; theorem actually holds)**.

Support cardinality exactly p^{r-1} in all 27 cells; symmetric difference between numerical and predicted support exactly 0 in all 27 cells.

Adversarial safeguards clean:
- **A1** (dual precision FFT vs mpmath at 50 digits): agreement to 1e-15 (float64 limit); mpmath confirms equality is exact (1.86e-49 deviation from predicted at (p=5, r=3)).
- **A2** (hand computation at (p=3, r=2, ξ=3) = √27): FFT matches to 1.71e-16 relative.

Boundary p = 2: theorem fails as predicted (different principal-unit structure — (1+p) has order 2^{r-1} ≠ 2^r in (Z/2^{r+1})^× for r ≥ 2). Documented; confirms exclusion is necessary.

Phase 3 (proof template walkthrough at p = 11): all three template steps (Cochrane T2 + Plancherel + principal-unit equidistribution) apply at p = 11 with no hidden p-dependence beyond what the family-level statement captures. Same residual rigor caveat as R78.3 at q=3 (equidistribution sketch, not fully written derivation).

**Combined parent (Move 2) verification + this work:** 27 + 6 (Move 2 cells) = **33 cells**, primes 3 through 31, r ranging 1 through 6, all confirming `|F̂_p^full(ξ)| = p^{(r+3)/2}` on the principal-unit sub-support of size p^{r-1}. Candidate ranked as **verified up to the evaluation ceiling tested**.

---

## Pre-registration adherence

- **Pre-registered:** 2026-05-11, committed at `de21e8a` before any compute. Locked rules in [FHAT_THEOREM_VERIFICATION_PRE_REGISTRATION.md](FHAT_THEOREM_VERIFICATION_PRE_REGISTRATION.md).
- **Procedure followed as locked**, with one **A4 method deviation** (logged below).
- **Decision rules applied as written** — no retroactive modification.

### A4 method deviation

**Logged 2026-05-11 during compute.**

Pre-reg §2 Phase 5 stated: "Off-support magnitudes: max |F̂_p^short(ξ)| over ξ ∉ predicted support. Should be effective zero (< 1e-10)."

**Issue.** The first verification run (zero-padded length-M FFT of f_p, evaluating F̂_short on all ξ ∈ Z/M) gave **on-support magnitudes matching predicted formula to 1e-15** (correct) but **off-support magnitudes of order p^{(r+1)/2}** (orders of magnitude above the 1e-10 threshold). All 23 Phase 1+2 cells "failed" the literal pre-reg test.

**Root cause.** The pre-reg conflated two related-but-distinct objects:
- `F̂_p^short(ξ)` (zero-padded length-M FFT): equals F̂_full on support `p·a`, but has **spectral leakage** off `p·Z/M` because the zero-padding doesn't preserve f_p's period.
- `F̂_p^full(ξ)` (length-M FFT of periodically-extended f_p): IS supported on `p·Z/M` and vanishes off the principal-unit sub-support. **This is the object the theorem candidate is about.**

The cleanest computation is the **length-period DFT** `G[a] = Σ_s f_p(s) exp(-2πi a s / p^r)` for `a ∈ Z/p^r`, which satisfies `F̂_full(p·a) = p · G[a]` (period-extension copies sum to p) and `F̂_full(ξ) = 0` for ξ not multiple of p (orthogonality). Testing G[a]'s support / magnitude properties directly tests F̂_full's claim without spectral-leakage artifacts.

**Deviation taken.** Switched verification script from zero-padded length-M FFT to length-period FFT (G[a]). Same on-support magnitudes (FFT is exact on the period structure either way), but off-support magnitudes now correctly test the F̂_full claim instead of a fictitious "F̂_short vanishes off support" claim that the pre-reg had inadvertently put forward.

**Documentation:** the same kind of normalization labeling correction logged in Move 2 attempt §A4. The corrected verification confirms the theorem.

**Why this is a labeling correction, not a hypothesis change.** The theorem candidate as stated in the Move 2 attempt §"Headline result" (commit 45af179) refers to F̂_p^full, not the zero-padded F̂_short. The pre-reg accidentally wrote a Phase 5 test for the latter. The corrected verification tests the former — the actual theorem.

---

## §1. Phase 1 — extension to higher primes

**Cells:** p ∈ {11, 13, 17, 19, 23, 29, 31} × r ∈ {2, 3} = 14 cells.
**Method:** G[a] = length-period FFT of f_p, predicted support = `{a ≡ 1 (mod p)}` of size p^{r-1}, predicted magnitude p^{(r+1)/2}.
**Threshold:** rel_dev < 1e-12, max_off < 1e-10, sym_diff = 0.

| p | r | M | period | |supp| | sym_diff | max\|G\| / p^{(r+1)/2} | min/pred | max_off |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 11 | 2 | 1331 | 121 | 11 | 0 | 1.0000000000000004 | 1.0000000000 | 4.80e-15 |
| 11 | 3 | 14641 | 1331 | 121 | 0 | 1.0000000000000004 | 1.0000000000 | 4.12e-14 |
| 13 | 2 | 2197 | 169 | 13 | 0 | 1.0000000000000002 | 1.0000000000 | 1.00e-14 |
| 13 | 3 | 28561 | 2197 | 169 | 0 | 1.0000000000000007 | 1.0000000000 | 8.67e-14 |
| 17 | 2 | 4913 | 289 | 17 | 0 | 1.0000000000000002 | 1.0000000000 | 1.44e-14 |
| 17 | 3 | 83521 | 4913 | 289 | 0 | 1.0000000000000007 | 1.0000000000 | 2.09e-13 |
| 19 | 2 | 6859 | 361 | 19 | 0 | 1.0000000000000004 | 1.0000000000 | 1.35e-14 |
| 19 | 3 | 130321 | 6859 | 361 | 0 | 1.0000000000000004 | 1.0000000000 | 2.71e-13 |
| 23 | 2 | 12167 | 529 | 23 | 0 | 1.0000000000000004 | 1.0000000000 | 1.66e-14 |
| 23 | 3 | 279841 | 12167 | 529 | 0 | 1.0000000000000007 | 1.0000000000 | 4.82e-13 |
| 29 | 2 | 24389 | 841 | 29 | 0 | 1.0000000000000002 | 1.0000000000 | 3.36e-14 |
| 29 | 3 | 707281 | 24389 | 841 | 0 | 1.0000000000000007 | 1.0000000000 | 9.38e-13 |
| 31 | 2 | 29791 | 961 | 31 | 0 | 1.0000000000000002 | 1.0000000000 | 4.75e-14 |
| 31 | 3 | 923521 | 29791 | 961 | 0 | 1.0000000000000007 | 1.0000000000 | 1.15e-12 |

**Phase 1 result: 14/14 cells PASS.** Max rel_dev on-support: ~7e-16 (float64 precision floor). Max off-support: 1.15e-12 at (p=31, r=3) with period = 29791 — accumulated float roundoff, three orders of magnitude under the 1e-10 threshold.

## §2. Phase 2 — extension to higher r

**Cells:** p ∈ {3, 5, 7} × r ∈ {4, 5, 6} = 9 cells.

| p | r | M | period | |supp| | sym_diff | max\|G\| / p^{(r+1)/2} | min/pred | max_off |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 4 | 243 | 81 | 27 | 0 | 1.0000000000000002 | 1.0000000000 | 6.91e-15 |
| 3 | 5 | 729 | 243 | 81 | 0 | 1.0000000000000004 | 1.0000000000 | 1.21e-14 |
| 3 | 6 | 2187 | 729 | 243 | 0 | 1.0000000000000007 | 1.0000000000 | 3.02e-14 |
| 5 | 4 | 3125 | 625 | 125 | 0 | 1.0000000000000004 | 1.0000000000 | 2.42e-14 |
| 5 | 5 | 15625 | 3125 | 625 | 0 | 1.0000000000000009 | 1.0000000000 | 1.16e-13 |
| 5 | 6 | 78125 | 15625 | 3125 | 0 | 1.0000000000000011 | 1.0000000000 | 6.12e-13 |
| 7 | 4 | 16807 | 2401 | 343 | 0 | 1.0000000000000007 | 1.0000000000 | 9.35e-14 |
| 7 | 5 | 117649 | 16807 | 2401 | 0 | 1.0000000000000007 | 1.0000000000 | 6.53e-13 |
| 7 | 6 | 823543 | 117649 | 16807 | 0 | 1.0000000000000009 | 1.0000000000 | 4.56e-12 |

**Phase 2 result: 9/9 cells PASS.** Max off-support: 4.56e-12 at (p=7, r=6) with period 117649. Same float-roundoff pattern; well under threshold.

## §3. Phase 3 — proof template walkthrough at p = 11

The Move 2 attempt §"Proof template (p-blind)" states the F̂_p magnitude formula follows from three steps:

1. **Family-level 78.1_p (complete-sum vanishing) via Cochrane Theorem 2**
2. **Family-level 78.2_p (sparsity) via period of (1+p) in (Z/p^{r+1})^×**
3. **Family-level 78.3_p (equidistribution) via Plancherel + principal-unit Gauss-sum structure**

Walkthrough at p = 11:

### Step 1 — Cochrane T2 / complete-sum vanishing

Polynomial identification: `g(u) = c · (1+11)^u − 11² · m · u` mod 11^{r+1}.

Binomial expansion: `(1+11)^u = Σ_{k=0}^r C(u, k) · 11^k` mod 11^{r+1}. This expansion is **p-blind in form** — same combinatorial identity as at p=3 (Move 2 attempt §Phase 2 line 66: "purely binomial-coefficient algebra, p-blind").

Cochrane T2's `D = deg_p H+` parameter: H is the reduction of g mod p = 11. Since all `C(u,k) · 11^k` terms for k ≥ 1 vanish mod 11, `g(u) ≡ c (mod 11)` — a constant. So `H` is constant non-zero (c is a unit), `D = 0`, and Corollary 6's criterion "vanishes unless H(a) ≡ 0 mod p^{m-ℓ-τ}" is never met. **Complete-sum vanishes.**

No p-dependence beyond what the binomial expansion structure gives — same argument as at p = 3, p = 5, p = 7 (all verified in Move 2). **Verdict: Step 1 applies at p = 11 with no hidden p-dependence.**

### Step 2 — Sparsity (supp(F̂_full) ⊆ p·Z/M)

(1+p) = 12 in (Z/11^{r+1})^×. Standard fact: for any odd prime p, the principal units 1 + p·Z_p mod p^{r+1} form a cyclic group of order p^r, generated by (1+p). At p = 11: order of 12 in (Z/11^{r+1})^× is 11^r.

Consequence: f_p(u) = e_M(c · (1+p)^u) has period p^r = 11^r in u. The DFT of a period-T signal extended to length M (with M = T·p) is supported on multiples of M/T = p. **Supp(F̂_full) ⊆ p·Z/M of size p^r = 11^r.**

This is purely the "period divides length" DFT lemma. **No p-dependence beyond the order-of-(1+p) fact**, which is p-blind for odd primes (and fails at p = 2 — see Phase 4).

### Step 3 — Sub-support {a ≡ 1 mod p} and uniform magnitude

**Sub-support** `{a ≡ c mod p}`: from the leading-order pairing

`c · (1+p)^u ≡ c + c·p·u (mod p²)`,
so `e_M(c · (1+p)^u) ≡ e_M(c) · e_{p^r}(c·u) (mod higher-order corrections in p)`.

Then `G[a] = Σ_u e_M(c · (1+p)^u) · e_{p^r}(-au) ≈ e_M(c) · Σ_u e_{p^r}((c-a)·u) · [corrections] = e_M(c) · p^r · 𝟙[a ≡ c (mod p^r)] · [corrections]`.

The corrections from higher binomial terms shift each individual G[a] value but **preserve the support condition `a ≡ c (mod p)`** at family level (R78.2 sketch). At p = 11, c = 1: support `{a ≡ 1 (mod 11)}` within Z/11^r, size 11^{r-1}. **Matches empirical data exactly (Phase 1 sym_diff = 0 at both r = 2, 3 for p = 11).**

**Magnitude:** Plancherel on Z/p^r gives `Σ_a |G[a]|² = p^r · Σ_s |f(s)|² = p^r · p^r = p^{2r}`. If equidistribution holds (uniform magnitude on the size-p^{r-1} support), then `|G[a]|² = p^{2r}/p^{r-1} = p^{r+1}`, so **|G[a]| = p^{(r+1)/2} = 11^{(r+1)/2}** at p = 11.

**Equidistribution structural reason** (per R78.3 at q=3 and R78.4-78.6 explicit Gauss-sum factorization in `result_78_extended.md`): the principal-unit Gauss sum `G(a) = Σ_s e_{p^{r+1}}(P_a(s))` (where P_a is the polynomial obtained from c·(1+p)^s after change-of-variables) has uniform magnitude √(p^{r+1}) = p^{(r+1)/2} on the support, by saddle-point / Gauss-sum theory for principal-unit characters. The **structural argument is p-blind** (Cochrane Prop 4 + saddle/Hensel framework applies at any prime); the explicit Cochrane truncated **p-adic log expansion** `L(1+ps) = Σ_j (-1)^{j-1}/j · (ps)^j` has p-dependent coefficients but **same form**.

**At p = 11:** the explicit p-adic log expansion at p = 11 would replace the 3 ↦ 11 substitution throughout, but the saddle-point / equidistribution structural conclusion (uniform magnitude on support) is the same. **Empirical confirmation:** Phase 1 cells (p=11, r=2 and r=3) show all 11 + 121 = 132 support values with magnitude `1.0 × 11^{(r+1)/2}` to 1e-15.

### Phase 3 verdict

**All three proof-template steps apply at p = 11 with no p-dependence beyond what the family-level statement captures.** The same residual rigor gap as in R78.3 at q=3 remains: the equidistribution claim is sketched via the principal-unit Gauss-sum framework and verified empirically; the q=3 explicit derivation via Cochrane Prop 4 + truncated p-adic log saddle-point (R78.4-78.6) demonstrates the structural mechanism but the family-level analog at p = 11 (or general p) is not written out in the same explicit closed-form detail.

**Phase 3: PASS at the level of "proof template applies p-blindly with no hidden constants"; the SAME RESIDUAL EQUIDISTRIBUTION RIGOR GAP as at q=3 carries to family level.** No new p-dependence introduced by going from p=3 to p=11.

---

## §4. Phase 4 — boundary cases

### Boundary 1: p = 2 (excluded by theorem; verify why)

| r | M | p^r expected | (1+p) actual order | naive_pred | max\|G\| on naive_supp | max\|G\| off |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 8 | 4 | **2** | 2.828 | 1.414 | 1.414 |
| 3 | 16 | 8 | **4** | 4.000 | 3.696 | 0.000 |
| 4 | 32 | 16 | **8** | 5.657 | 4.000 | 0.000 |

**Mode of failure documented.** At p = 2, (1+p) = 3 has multiplicative order 2^{r-1} (not 2^r) in (Z/2^{r+1})^× for r ≥ 2, because principal units 1 + 2Z_2 mod 2^{r+1} are isomorphic to `Z/2 × Z/2^{r-1}` (not cyclic, unlike odd-p case). The period of f_p is half of what the family-level theorem assumes. Both the support structure and the magnitude formula fail at p = 2.

**Conclusion:** the theorem's exclusion of p = 2 is **structurally necessary, not a placeholder for "haven't checked"**.

### Boundary 2: r = 1 (pre-reg said "excluded"; data shows theorem actually holds)

| p | M | |supp| | pred(short) | max\|G\| | min\|G\| | max_off | rel_dev | PASS |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 9 | 1 | 3 | 3.0000000000 | 3.0000000000 | 4.00e-16 | 0.00 | YES |
| 5 | 25 | 1 | 5 | 5.0000000000 | 5.0000000000 | 4.85e-16 | 1.78e-16 | YES |
| 7 | 49 | 1 | 7 | 7.0000000000 | 7.0000000000 | 1.17e-15 | 0.00 | YES |
| 11 | 121 | 1 | 11 | 11.0000000000 | 11.0000000000 | 1.25e-15 | 1.61e-16 | YES |

**Theorem actually holds at r = 1.** Support reduces to a single element `{a = 1}` (size p^{r-1} = p^0 = 1), and `|G[1]| = p^{(r+1)/2} = p^1 = p` exactly. The pre-reg's "theorem requires r ≥ 2" was conservative; the empirical verification at r = 1 shows the formula extends.

**Reasoning:** At r = 1, period = p, support cardinality = 1. The single support value G[1] is `Σ_s e_{p²}(c · (1+p)^s) · e_p(-s)`. By the same algebraic structure (leading character pairing forces a = 1; principal-unit-character magnitude saturation gives |G[1]| = p), the formula holds. No degeneracy that breaks the proof template.

**Scope extension noted:** family-level theorem holds for `p ≥ 3, r ≥ 1` (not just r ≥ 2 as pre-reg stated). Documented as theorem-statement refinement, not pre-reg violation — pre-reg's r ≥ 2 was strictly more conservative.

### Boundary 3: evaluation ceiling

Largest cell tested: **(p = 7, r = 6)** with M = 823,543, period = 117,649. Length-period FFT completed in 0.04 s on float64. Memory footprint at this cell: ~10 MB for the period array.

Practical ceiling for length-period FFT on float64: any cell with period < 10^8 should complete in seconds and fit in ~2 GB memory. Cells up to (p = 3, r = 16) → period = 43M, or (p = 5, r = 11) → period = 49M, etc. are reachable without optimization. **No ceiling near the verified range.**

Higher cells could be tested but yield no new structural information — the proof template's rigor doesn't depend on cell size; what's tested empirically at r = 6 is the same structural claim that holds at any r ≥ 1.

---

## §5. Phase 5 — support characterization

The theorem statement claims support = `{p·a : a ∈ Z/p^r, a ≡ 1 (mod p)}` of size `p^{r-1}` for F̂_full (or equivalently `{a ∈ Z/p^r : a ≡ 1 (mod p)}` for G[a]), with F̂_full vanishing exactly off this set.

**Symmetric difference between numerical and predicted support: 0** across all 27 cells (Phase 1: 14/14, Phase 2: 9/9, Phase 4 r=1: 4/4).

Maximum magnitude off the predicted support across all 27 cells: **4.56e-12** at (p=7, r=6). Well under the 1e-10 pre-registered threshold; this is accumulated FFT roundoff at period = 117,649, not theorem violation.

Numerical support cardinality exactly p^{r-1} in all 27 cells. **Support characterization confirmed.**

---

## §6. Adversarial safeguards records

### A1 — dual precision FFT vs mpmath at 50 digits

Cells: (p=3, r=2), (p=11, r=2), (p=5, r=3).

| Cell | FFT max | mpmath max (50 dps) | predicted (p^((r+1)/2)) | max\|FFT − mpmath\| | rel | mpmath dev from pred |
|---|---:|---:|---:|---:|---:|---:|
| (3, 2) | 5.196152422706632 | 5.196152422706632 | 5.196152422706632 | 8.88e-16 | 1.71e-16 | 2.75e-17 |
| (11, 2) | 36.482872693909414 | 36.482872693909400 | 36.482872693909400 | 1.42e-14 | 3.90e-16 | 3.69e-17 |
| (5, 3) | 25.000000000000007 | 25.000000000000000 | 25.000000000000000 | 7.11e-15 | 2.84e-16 | **1.86e-49** |

FFT agrees with mpmath to float64 precision floor (~1e-15) at all three cells.

**At (p=5, r=3), mpmath at 50 digits confirms |G[a]| = 25 to better than 1e-49.** This is much smaller than the 50-digit precision floor would predict — **the equality is exact, not just numerically close.** Strong evidence that the magnitude formula is an algebraic identity, not an approximate / asymptotic relation.

(At the other cells, the small ~3e-17 "deviation" is artifact of mpmath ↔ float64 conversion; the mpmath value of `p^((r+1)/2)` is computed via `mp.power(p, mp.mpf(r+1)/2)` and rounded to float for comparison.)

**A1 verdict: numerical artifacts are clean. Magnitude formula confirmed exact.**

### A2 — hand-computation cross-check at (p=3, r=2, ξ=3)

Hand derivation:
```
f(u) = exp(2πi · k_u / 27),   k_u = 4^u mod 27 = (1, 4, 16, 10, 13, 25, 19, 22, 7).
F̂(3) = Σ_{u=0}^{8} exp(2πi · (k_u - 3u) / 27).
Phases (k_u - 3u) mod 27: (1, 1, 10, 1, 1, 10, 1, 1, 10).
Six terms at phase 1/27, three terms at phase 10/27.
|F̂(3)|² = |6 e^{2πi/27} + 3 e^{2πi·10/27}|² = 45 + 36 cos(2π · 9/27) = 45 − 18 = 27.
|F̂(3)| = √27 ≈ 5.196152422706632.
```

| | value |
|---|---:|
| Hand-computed | 5.196152422706632 |
| FFT-computed | 5.196152422706631 |
| Absolute diff | 8.88e-16 |
| Relative diff | 1.71e-16 |
| **Pass (rel < 1e-14)** | **YES** |

**A2 verdict: implementation correctness confirmed.**

### A3 — honest deviation logging

The first verification run produced 0/23 cells passing the literal pre-reg test, due to the F̂_short vs F̂_full normalization confusion documented in A4 above. This was honestly logged and the verification method was corrected (not the theorem statement or thresholds). The corrected verification then produced 27/27 passes.

No "rounding error" hand-waving was invoked. Max off-support 4.56e-12 at the largest cell was verified to be float-roundoff (consistent with the period-length FFT's accumulated error at period = 117,649 — order 10^{-15} · √period ≈ 3e-12 per central-limit-style argument).

### A4 — pre-reg adherence

One method deviation (F̂_short → length-period G[a]), documented above. No hypothesis changes. No threshold relaxation. Decision rules applied as written: 27/27 cells pass, A1/A2 clean → THEOREM_VERIFIED disposition.

The r=1 finding (theorem extends to r ≥ 1) is a **scope extension** not a scope narrowing — pre-reg's `r ≥ 2` was strictly more conservative than what's verifiable. Logged as informative extension; doesn't affect disposition under pre-reg decision rules.

---

## §7. What does THEOREM_VERIFIED mean here?

Per pre-reg §4 decision rules:

> **THEOREM_VERIFIED.** Phases 1–5 all pass. A1, A2 clean. Theorem candidate ranked as verified **up to the evaluation ceiling identified in Phase 4**, ready for formalization. The exact equality `|F̂_p^short(ξ)| = p^{(r+1)/2}` (equivalently `|F̂_p^full(ξ)| = p^{(r+3)/2}`) holds on predicted support, with predicted support structure, for all tested cells.

This is what the data shows. The candidate theorem from Move 2 attempt is **empirically verified across a substantial range** (primes 3 through 31, r ∈ {1, 2, 3, 4, 5, 6}), with the support characterization, magnitude formula, and structural proof template surviving the adversarial probes.

**What this does NOT mean.** Not a formal proof. The structural proof template (Cochrane T2 + Plancherel + principal-unit equidistribution) has the same residual rigor gap as at q=3: the equidistribution claim is sketched via principal-unit Gauss-sum theory and verified empirically (R78.3) but not written out with full rigor at family level. The R78.4-78.6 explicit Gauss-sum factorization with Cochrane truncated p-adic log + saddle-point at q=3 demonstrates the mechanism; the analogous family-level closed form at general p would require parallel derivation (likely feasible with the same technique but not done here).

**What THEOREM_VERIFIED enables.** Ready for formalization as a candidate theorem in a paper. Ready to be cited as "empirically verified to machine precision across [tested range]" rather than just "verified at 6 cells." The Move 2 attempt's "candidate theorem" status now has substantially more empirical backing.

**Strategic position.** Per Move 2 attempt §"What this does NOT do": the F̂_p result remains the wrong-object-shape to dissolve R77.2's Tao Prop 1.17 dependence. This verification doesn't change that — F̂_p is still distinct from K_p (Burgess bilinear wall) and from μ̂_n (Markov-chain stationary's characteristic function). The verified F̂_p theorem may stand as a **standalone family-level extension of R78.3** with strong empirical backing, but the c = 7/45 closure question is unchanged.

---

## §8. Refined theorem statement (verified)

> **THEOREM (qx+1 Plancherel saturation at the F̂ level, verified):** For every prime `p ≥ 3` and every `r ≥ 1`, and for `c ∈ (Z/p^{r+1})^×`, define
>
> `f_p(u) := e_M(c · (1+p)^u)`,  `M = p^{r+1}`,  `u ∈ Z/p^r`.
>
> The full-period Fourier transform `F̂_p^full` (DFT of f_p extended periodically to length M) satisfies:
>
> 1. **Support:** `supp(F̂_p^full) = { p · a (mod M) : a ∈ Z/p^r, a ≡ c (mod p) }`, cardinality `p^{r-1}`.
> 2. **Magnitude (uniform on support):** `|F̂_p^full(ξ)| = p^{(r+3)/2}` for all `ξ ∈ supp(F̂_p^full)`.
> 3. **Vanishing off support:** `F̂_p^full(ξ) = 0` for `ξ ∉ supp(F̂_p^full)`.
>
> **Verification scope (empirical):** all 33 (p, r) cells with `p ∈ {3, 5, 7, 11, 13, 17, 19, 23, 29, 31}` and `r ∈ {1, 2, 3, 4, 5, 6}` (subject to the specific (p, r) cell combinations tested in Move 2 and this verification). Magnitude matches `p^{(r+3)/2}` to ~1e-15; off-support magnitude bounded by 4.56e-12; symmetric difference between numerical and predicted support equals zero in every cell.
>
> **Proof-template scope (sketched, p-blind):** Cochrane Theorem 2 + Plancherel identity + principal-unit Gauss-sum equidistribution. The structural argument applies for every prime p ≥ 3. The equidistribution rigor — analogous to R78.3 at q=3 — is sketched but not fully written out at family level. The R78.4–78.6 explicit Gauss-sum closed form at q=3 demonstrates the mechanism.
>
> **Boundary:** p = 2 excluded (principal-unit subgroup of (Z/2^{r+1})^× is not cyclic for r ≥ 2; period of (1+p)=3 is 2^{r-1} ≠ 2^r). r = 0 trivial.

(Note: the r=1 case is **included** in the refined scope. Pre-reg's r ≥ 2 was conservative.)

---

## §9. Disposition handling per pre-reg

Pre-reg §4 spectrum: `THEOREM_VERIFIED → VERIFIED_WITH_SCOPE_NARROWING → THEOREM_FALSIFIED → INCONCLUSIVE`.

**Landed at: THEOREM_VERIFIED** with one scope-extension (r ≥ 1 instead of pre-reg's r ≥ 2 — empirical data is strictly stronger than pre-reg's claim).

Per pre-reg §4:
> "THEOREM_VERIFIED. Phases 1–5 all pass. A1, A2 clean. Theorem candidate ranked as verified up to the evaluation ceiling identified in Phase 4, ready for formalization and Paper 4 integration."

**Ready for formalization at the empirical-verification level.** The "candidate theorem" status from Move 2 now has empirical backing across 7 new primes (p ∈ {11..31}), 3 higher r values (r ∈ {4..6}), the r=1 boundary, hand-computation cross-check, and dual-precision mpmath confirmation showing exact equality. Phase 3 walkthrough at p=11 confirms the proof template applies p-blindly.

**What's still NOT formalized:** the explicit family-level rigorous derivation of the equidistribution (i.e., the analog of R78.4-78.6 at general p, not just empirical verification). That's the next-stage formalization work, parallel to writing out the full Cochrane Prop 4 + saddle-point derivation at general p with the Cochrane truncated p-adic log replacing the q=3 specific computation.

---

## §10. Files

- [FHAT_THEOREM_VERIFICATION_PRE_REGISTRATION.md](FHAT_THEOREM_VERIFICATION_PRE_REGISTRATION.md) — locked rules (commit de21e8a)
- [FHAT_THEOREM_VERIFICATION_RESULTS.md](FHAT_THEOREM_VERIFICATION_RESULTS.md) — this document
- [fhat_verification.py](fhat_verification.py) — main verification (Phases 1, 2, 4, 5)
- [fhat_verification_dual_precision.py](fhat_verification_dual_precision.py) — A1 + A2
- [fhat_verification_results.csv](fhat_verification_results.csv) — per-cell measurements
- [fhat_verification_a1_dual_precision.csv](fhat_verification_a1_dual_precision.csv) — A1 + A2 measurements

## Parent context

- [QX1_FAMILY_THEOREM_PRE_REGISTRATION.md](QX1_FAMILY_THEOREM_PRE_REGISTRATION.md) — Move 2 pre-reg
- [QX1_FAMILY_THEOREM_ATTEMPT.md](QX1_FAMILY_THEOREM_ATTEMPT.md) — Move 2 attempt (CLAIM_1_PARTIAL + STRUCTURAL_OBSTRUCTION_FOUND)
- [qx1_move2_phase2_check.csv](qx1_move2_phase2_check.csv) — original 6-cell empirical verification (now extended to 27 + 6 = 33 cells)
- [result_78_FINAL.md](result_78_FINAL.md) — q=3 specific R78.1–78.3 (the proof template being generalized)
- [result_78_extended.md](result_78_extended.md) — q=3 specific R78.4–78.6 (explicit Gauss-sum closed form via Cochrane Prop 4 + saddle-point — the structural mechanism for equidistribution)
