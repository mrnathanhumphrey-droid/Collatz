# Phase 2 — Candidate bridges Φ : F̂_p → μ̂_n

Per pre-reg, five categories enumerated. For each candidate: precise math statement, structural prediction, empirical test design, falsifiable criterion. Phase 3 will run the tests where they can be run; structural falsifications are noted here.

## Category (A): Spectral conjugacy via R77.2's T_diag / T_3

### Candidate A1: F̂_p eigenvalue in T_diag's spectrum

**Statement:** the diagonal operator T_diag = (1/5)·[[1,1],[4,4]] (R77.1 Theorem 77.1) has spec {0, 1}. The conjecture is that |F̂_p(ξ)|² normalised by some factor lies on a (1, 4)-type eigenvector of an extended diagonal operator that subsumes both F̂_p and the bilinear-pair moment M_n(η).

**Structural check.** T_diag is a 2×2 operator on (P_+, P_-)_n where P_±(c) = Σ_{ξ ≡ c} |μ̂_n^±(ξ)|² is a bilinear (squared) moment on the **class-resolved μ̂_n**. F̂_p is **not bilinear in any μ̂**; F̂_p is the linear Fourier transform of f_p, not the squared bilinear average of a μ̂. They live at different levels of the algebraic hierarchy (linear vs bilinear). For F̂_p to plug into T_diag, one would need a "diagonal P_+(c) for f_p" — but f_p is deterministic and class-resolution doesn't apply (there's no random class for a deterministic phase function).

**Prediction:** A1 is structurally falsified by level-of-algebra mismatch. F̂_p is a linear functional of the deterministic f_p; T_diag is a bilinear functional on the probabilistic μ̂_n^±. No match.

**Empirical test:** none meaningful — the objects are incommensurable. **Falsified at Phase 2.**

### Candidate A2: F̂_p eigenvalue in T_3 (order-3 companion)

**Statement:** T_3 from R77.2 §2.3 has companion-matrix spec {1/2, 1/4, 1/8}. Is the F̂_p magnitude formula p^{(r+3)/2} = p^{r/2}·p^{3/2} related to (1/2)^r somehow? Specifically, |F̂_p(ξ)|² = p^{r+3} grows like (1/2)^{-r} would predict if we invert.

**Structural check.** The spec {1/2, 1/4, 1/8} comes from R76 §10's least-squares fit ε_n ≈ A(1/2)^n + B(1/4)^n + C(1/8)^n on n=3,4,5. (1/2)^n is the leading rate of the μ̂_n bilinear-pair-moment system; (1/4) and (1/8) are sub-leading corrections.

The F̂_p magnitude p^{(r+3)/2} grows with r (it's not a decay rate; it's a saturation magnitude). The natural quantity to compare to (1/2)^n would be a *normalised* F̂_p that DECAYS in r. Candidate: |F̂_p(ξ)|²/M² = p^{r+3}/p^{2(r+1)} = p/p^{r-1} = p^{2-r}. At p = 3 this is 3^{2-r} — a decay rate of 1/3 per r, not 1/2.

So even the natural "normalisation by total Plancherel mass" doesn't produce the 1/2 rate. The candidate match would require:
> 3^{2-r}  ↔  (1/2)^n  with some r ↔ n correspondence.

These are not equal for any r ↔ n relationship; the bases (3 vs 2) differ.

**Prediction:** A2's exponents don't match; the F̂_p saturation decays at rate p (under natural normalisation) while μ̂_n's S_n - 7/15 decays at rate 1/2. **Falsified by rate mismatch at structural level.**

**Empirical test:** could compute |F̂_p|²/M² at p=3 and (S_n - 7/15) at matched levels and confirm the rate mismatch. Not informative beyond the structural answer.

## Category (B): Averaging/iteration

### Candidate B1: μ̂_n as average of F̂_p-type objects over (c, ξ)

**Statement:** maybe μ̂_n(ξ) = E_c[F̂_p(ξ, c)] for some average over c — i.e., the expectation in μ̂_n is over c (= initial state of the Markov chain) where F̂_p depends on c parametrically.

**Structural check.** F̂_p depends on c only as a phase shift (the central computation: G[a; c] · e_M(c) has |G[a; c]| = p^{(r+1)/2} independent of c on the c-support). So averaging F̂_p(ξ, c) over c gives:
- If we average **magnitudes:** E_c |F̂_p(ξ, c)| = p^{(r+3)/2} (a constant; no n^{-A} decay).
- If we average **complex values:** E_c F̂_p(ξ, c) is a sum of unit-modulus complex numbers Σ_c α(c) with possibly some cancellation, but the magnitude of the average is bounded by p^{(r+3)/2} from above with no obvious mechanism to push it to 0 faster than 1/√(number of terms).

In neither case do we get the Tao Prop 1.17 n^{-A} super-polynomial decay. Averaging F̂_p over c at fixed (p, r) cannot produce n-dependent super-polynomial decay because the underlying object doesn't have n.

**Prediction:** B1 cannot reproduce n^{-A} decay. **Falsified at Phase 2.**

### Candidate B2: μ̂_n as weighted sum over multiple resolutions

**Statement:** μ̂_n(ξ) = Σ_{r=1..n} w_r(ξ) · F̂_3^{(r)}(ξ) for some weights w_r and depth-r F̂_3 objects.

**Structural check.** Tao (1.26) gives Syrac(Z/3^n) ≡ Σ_{i=1..n} 3^{i-1}·2^{-a_{[1,i]}} (mod 3^n), so
> μ̂_n(ξ) = E[e^{-2πi ξ/3^n · Σ_i 3^{i-1}·2^{-a_{[1,i]}}}]
>         = E[ Π_i e^{-2πi ξ/3^n · 3^{i-1} · 2^{-a_{[1,i]}}} ]
>         = E[ Π_i e^{-2πi ξ·3^{i-1}/3^n · 2^{-a_{[1,i]}}} ].

This is **NOT** a sum of F̂_p-type objects. It's a product of n random-phase factors, where each factor depends on the partial sum a_{[1,i]} (not on a single Geom). The product structure means the n factors are NOT independent (a_{[1,i+1]} = a_{[1,i]} + a_{i+1}), which is precisely what makes the analysis hard (Tao §7).

A weighted-sum decomposition does NOT exist for μ̂_n; the natural decomposition is multiplicative-over-products-with-correlated-factors. So B2's premise (sum decomposition) doesn't match the object.

**Modified B2':** could one bound μ̂_n by a product of F̂_3-type magnitudes? Naive product bound:
> |μ̂_n(ξ)| ≤ Π_i |E_a[e^{-2πi ξ·3^{i-1}·2^{-a}/3^n}]|.

But the a_i are NOT independent across i (they enter via the partial sums a_{[1,i]}), so the product-of-expectations is NOT an upper bound. **Independence fails.**

If one IGNORES the partial-sum structure and PRETENDS the n contributions are independent — replace a_{[1,i]} by an independent Geom(2)-distributed variable b_i — then one gets the *toy* characteristic function
> μ̃_n(ξ) = Π_i E_{b_i ~ Geom(2)}[e^{-2πi ξ·3^{i-1}·2^{-b_i}/3^n}].

This toy is NOT μ̂_n. Tao 2022 explicitly notes this:
> (Tao paper, near (1.26)) "the expression (1.26) does not obviously resolve into such a sum of independent random variables, unfortunately." Tao §7 is precisely the analysis that overcomes this non-independence.

**Prediction:** B2 fails by the iid-vs-correlated obstruction. The replacement μ̃_n is not equal to μ̂_n; bounds on μ̃_n don't transfer to μ̂_n without controlling the dependence (which is Tao §7's job).

**Empirical test (informative even if falsifying):** compute μ̃_n empirically by sampling iid (b_1,...,b_n) ~ Geom(2)^n, and μ̂_n empirically by sampling correlated (a_1,...,a_n) per (1.26). Compare. If they agree to MC precision, that's surprising (and worth investigating); if they disagree, the gap quantifies the iid-violation contribution.

**Status:** **structurally pre-falsified by Tao's own remark**; would be Phase 3 quantification of the gap.

## Category (C): Coupling via p-adic lift

### Candidate C1: μ̂_n value = restriction of F̂_p to a coset

**Statement:** for p = 3 and n = r+1, can μ̂_n(ξ) be expressed as a specific value or coset-restriction of F̂_3 on Z/3^n?

**Structural check.** F̂_3 is supported on {3·a : a ≡ c (mod 3)} ⊂ Z/3^{r+1} = Z/3^n. μ̂_n's Tao Prop 1.17 is the decay at ξ with 3 ∤ ξ — the COMPLEMENT of F̂_3's support. So at the level-aligned scale, F̂_3(ξ) ≡ 0 on the set where μ̂_n's interesting behaviour lives. **They are concentrated on disjoint frequencies.**

To restrict F̂_3 to a coset where μ̂_n lives, one would need either:
- A different definition of F̂_3 (different f_p, not (1+p)^u) — abandons the verified theorem.
- A natural map from {ξ: 3∤ξ} into the F̂_p support {3·a: a ≡ c mod 3} — no such canonical map exists (the support shift ξ → ξ' = 3·ξ' mod 3^n is forced and reduces the problem mod 3^{n-1}, not solving anything).

**Prediction:** C1 falsifies by support disjointness. **Falsified at Phase 2.**

### Candidate C2: parity-sequence parameterisation

**Statement:** the Tao 2-adic valuation sequence (a_1, ..., a_n) for the Markov chain Syrac is related to the (1+p)^u orbit visited in F̂_p (which is a deterministic cyclic walk of period p^r).

**Structural check.** Syrac's chain step is x ↦ (3x + 1)/2^{a_1} where a_1 = v_2(3x+1) ∈ Geom(2). F̂_p's underlying function is f_p(u) = e_M(c·(1+p)^u) — a *deterministic* cyclic walk on Z/M generated by multiplication by (1+p), period p^r.

The (1+p) generator in F̂_p has nothing to do with the 2^{-a_i} multipliers in the Tao chain. (1+p) is on the principal-unit side of (Z/M)^×, while 2^{-1} is on a different cyclic subgroup of (Z/3^n)^× (the subgroup ⟨2⟩, which has order 2·3^{n-1}/gcd(...) — actually 2 has multiplicative order = 2·3^{n-1} mod 3^n for n ≥ 1).

The two walks (1+p)^u and 2^{-a_1-...-a_i} are unrelated cyclic walks on different subgroups of (Z/3^n)^×.

**Prediction:** C2 fails — the two parameterisations are on different multiplicative subgroups with no obvious morphism. **Falsified by group-theoretic mismatch.**

### Candidate C3: Tao p-adic perspective

**Statement:** Tao Remark 1.13 viewed Syrac(Z/3^n) as the projection of a 3-adic random variable Syrac(Z_3). F̂_p uses (1+p) which generates principal units in Z_p^×. Some 3-adic constant comparing the principal units in Z_3^× to Syrac(Z_3)?

**Structural check.** Syrac(Z_3) is a probability distribution on Z_3 (the closure under the level-n projections). F̂_p's (1+p)^u takes values in the principal-unit subgroup 1 + 3Z_3 ⊂ Z_3^×. These are subsets of different ambient sets (Z_3 vs Z_3^×) with non-overlapping support structures.

Even if one defines a 3-adic Fourier transform Z_3 → C and applies it to Syrac(Z_3)'s distribution, the answer is μ̂_n (in its 3-adic-continuation form). The 3-adic Fourier transform of (1+3)^u as a deterministic function would just be F̂_3 in 3-adic-continuation. They are still distinct objects in different categories (probability measure vs deterministic function).

**Prediction:** C3 is the closest "common ambient space" framing, but the objects within it remain distinct (random measure vs deterministic phase function). **Falsified by category-of-object mismatch (random vs deterministic).**

## Category (D): Parameter specialization (p = 3)

### Candidate D1: F̂_3 (p=3, r=n-1) ↔ μ̂_n

**Statement:** at p = 3, r = n - 1 to match scale 3^n, can the verified F̂_3 magnitude formula 3^{(r+3)/2} = 3^{(n+2)/2} translate into a bound on μ̂_n?

**Structural check.** F̂_3 at (p=3, r=n-1) has magnitude 3^{(n+2)/2} on its sparse support (size 3^{n-2}). Total Plancherel mass = (3^{(n+2)/2})² · 3^{n-2} = 3^{n+2} · 3^{n-2} = 3^{2n}, which equals M² for M = 3^n. So F̂_3's saturation accounts for ALL the Plancherel mass on Z/3^n.

μ̂_n has Plancherel mass S_n = Σ_{ξ:3∤ξ} |μ̂_n(ξ)|² → 7/15 (R75/R76). This is O(1), not 3^{2n}.

The two objects live on the same Z/3^n but with completely different total Plancherel masses (3^{2n} vs O(1)). They are not the same kind of object.

To go from F̂_3 to μ̂_n one would need a normalisation: μ̂_n is a *characteristic function* of a probability measure π_n; F̂_3 is the *Fourier transform of a unit-modulus function*. The natural normalised F̂_3 (dividing by M) gives values of magnitude 3^{(n+2)/2}/3^n = 3^{(2-n)/2} which decays *exponentially* in n — much faster than the n^{-A} decay of μ̂_n.

So if one tries D1's natural normalised candidate Φ(F̂_3) = F̂_3/M with prediction |Φ| ≤ 3^{(2-n)/2} on support, this is a STRONGER bound than Tao 1.17's n^{-A}. **But it's a bound on F̂_3/M, not on μ̂_n. There is no identity F̂_3/M = μ̂_n.**

The empirical question: at p=3, what is the relationship between F̂_3(ξ)/M at ξ ∈ supp and μ̂_n(ξ)? Both are complex numbers on Z/3^n. We can compute both at small (n, ξ) and compare.

**Empirical test (Phase 3):** at small n ∈ {2, 3, 4, 5}, compute F̂_3(ξ)/M and μ̂_n(ξ) at the same ξ. Either:
- F̂_3(ξ)/M ≠ μ̂_n(ξ) — D1 falsified.
- F̂_3(ξ)/M = μ̂_n(ξ) — would be remarkable; A1 safeguard flag (numerical match without derivation).

**Prediction (pre-empirical):** they will differ. F̂_3/M is supported only on ξ ∈ 3·Z/3^n; μ̂_n at ξ ∉ 3·Z/3^n is nonzero (per Tao 1.17 it's small but nonzero). So at ξ ∉ 3·Z/3^n, F̂_3(ξ)/M = 0 ≠ μ̂_n(ξ). **Almost certainly falsified.**

But at ξ ∈ 3·Z/3^n the comparison is meaningful: F̂_3 supports it with magnitude 3^{(2-n)/2}/M·M = 3^{(2-n)/2} (after the /M normalisation), and μ̂_n's value there equals (by Plancherel + Tao 1.22) μ̂_{n-1}(ξ/3) — i.e., μ̂_n at a multiple-of-3 frequency equals μ̂_{n-1} at the reduced frequency. This is a level-reduction identity, NOT a connection to F̂_3 saturation.

**Status:** D1 testable, expected to falsify on the support-disjointness front, but useful to confirm empirically. Run in Phase 3.

## Category (E): No-bridge baseline

### Candidate E1: F̂_p and μ̂_n are simply different objects

**Statement:** the structural analyses in (A)-(D) all fail because the two objects are *genuinely* different — deterministic cyclic-group character sum on Z/p^{r+1} vs Markov-chain stationary characteristic function on Z/3^n. The only connections are vocabulary (Plancherel, characters) and ambient group (some Z/p^k).

**Implication:** the bound |F̂_p(ξ)| = p^{(r+3)/2} is informationally orthogonal to |μ̂_n(ξ)| ≪_A n^{-A}. R77.2's Nisoli closure cannot be fed from F̂_p; it must come from either:
- Tao §7.2-7.4 bookkeeping (Route 1; concurrent agent's task).
- Burgess-type bilinear character-sum bound (Route 2; R78 wall).
- A truly novel third route not articulated in project documents (Route 3, the hypothetical this hunt tests).

If Phase 2 enumeration exhausts plausible Route-3 candidates and they all falsify, then E1 is the disposition: **NO_BRIDGE_FOUND**, and c=7/45 closure relies on Routes 1 and 2 (or new techniques).

## Falsification summary at Phase 2

| Candidate | Disposition pre-Phase-3 | Reason |
|---|---|---|
| A1 (T_diag eigenvector) | Falsified | Level-of-algebra mismatch (linear vs bilinear) |
| A2 (T_3 spectrum match) | Falsified | Rate mismatch (3-base vs 2-base exponentials) |
| B1 (μ̂_n = E_c[F̂_p]) | Falsified | No mechanism for n^{-A} decay |
| B2 (μ̂_n as product of F̂_p) | Falsified | Tao explicitly notes iid factorisation fails |
| C1 (coset restriction) | Falsified | Support disjointness at level alignment |
| C2 (parity-sequence) | Falsified | (1+p) vs ⟨2⟩ are different cyclic subgroups |
| C3 (3-adic ambient) | Falsified | Random vs deterministic category mismatch |
| D1 (F̂_3 ↔ μ̂_n at matched scale) | **Test in Phase 3** | Likely falsified by support disjointness, but worth confirming |
| E1 (no-bridge baseline) | Default null | Consistent with all of the above |

Phase 3 will: (i) compute F̂_3 and μ̂_n at matched (n, ξ) to confirm D1's falsification, and (ii) do the "toy independent product" test of B2' to quantify the iid-violation gap. (i) is the substantive test; (ii) is informative for understanding why the bridge fails.

## What would have to be true for BRIDGE_FOUND_*

For one of (A)-(D) to advance to Phase 4:
- The empirical test must show agreement at ALL frequencies ξ (both inside and outside F̂_p support), at multiple n, at multiple bases.
- A structural derivation of Φ must exist that doesn't import Tao §7.
- An out-of-sample cell (held out from any tuning) must agree.

Given the Phase 2 structural falsifications, **the pre-registered null (H_NULL: F̂_p and μ̂_n are structurally distinct, no bridge) is strongly favored.**
