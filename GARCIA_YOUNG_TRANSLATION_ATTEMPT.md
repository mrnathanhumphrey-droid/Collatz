# Garcia–Young 2023: direct translation attempt + adversarial check

**Date:** 2026-05-11. Companion to `GARCIA_YOUNG_TRANSLATION.md` (notation) and `GARCIA_YOUNG_MECHANISM.md` (where √q arises in GY's proof).

This file attempts the translation of our bilinear sum into GY's framework, applies their Lemma 2.15 (or analogous lemma), and reports the disposition with adversarial checks A1–A4 from the brief.

## 1. Setup

Our bilinear sum (PRECISE_ASK §3 Form B):
> **Σ := Σ_{a ≡ 1 (mod p) in Z/p^r} 1̂(p·a) · F̂(p·a)**.

Plancherel identity (`result_78_FINAL.md` §"What this gives for Kalafatelis's eq 190"):
> `S_partial(r, c=1, m=0) = (1/q) · Σ`,
i.e., our bilinear equals `q · S_partial`.

Equivalent Form A: `|S_partial| ≤ C · √N` with explicit `C`.

## 2. Reframe `F̂(p·a)` using GY's Lemma 2.15 (attempt)

`F̂(p·a)` is a length-`q` DFT of the period-`p^r` function `f(u) = e_q(c · (1+p)^u)` extended periodically to length `q = p^{r+1}`. By period-extension, `F̂(p·a) = p · G[a]` where
> `G[a] = Σ_{u (mod p^r)} e_q(c · (1+p)^u) · e_{p^r}(-au)`.

After the parametrization change `u ↦ φ(u)` discussed in `GARCIA_YOUNG_TRANSLATION.md` §3, this **could in principle** become a sum of GY's `S_{q,d}`-type form, but with two important deviations:

**Deviation 1: nonlinear coordinate.** The substitution `(1+p)^u = 1 + p·φ(u)` is **nonlinear in `u`** (since `φ(u) ≡ u + p·u(u-1)/2 + ... mod p^r`). Plugging into `G[a]`:
> `G[a] = Σ_{u} ψ(1+p·φ(u)) · e_{p^r}(-au)`,
where `ψ(1+p·x) = e_q(c · L_q(1+p·x))` is the Postnikov-rewrite of the character (with `a_ψ = c`).

To use GY's Lemma 2.15, we would need this sum to be in the form `Σ_u ψ(1+du) · e_q(dku)` for some `d, k`. With `d = p`, the prefactor matches (`ψ(1+p·φ(u))` is in GY's form **if and only if** we substitute `v := φ(u)`). Doing so:
> `G[a] = Σ_{v (mod p^r)} ψ(1+pv) · e_{p^r}(-a · φ⁻¹(v))`,
where `φ⁻¹` is the inverse bijection. **`φ⁻¹` is itself a Cochrane-style power series**, so `e_{p^r}(-a·φ⁻¹(v))` is **not** a linear character in `v`. It is a higher-degree exponential in `v`.

Concretely, `φ⁻¹(v) = v · (1 + O(p))`, so `e_{p^r}(-a·φ⁻¹(v)) = e_{p^r}(-a·v) · e_{p^r}(-a · (φ⁻¹(v) - v))`, and the correction `φ⁻¹(v) - v` is of order `pv²` (leading term). The "extra" phase is `e_{p^r}(-a · p · v · O(v))` — but since `p · v = O(p^r/p · p) = O(p^r)` mod `p^r` is a unit-magnitude phase: this is **not** a small perturbation in the `e_{p^r}` argument unless we keep careful track of it.

**Conclusion of Deviation 1:** GY's Lemma 2.15 in its literal form does **not** evaluate our `G[a]`. The `e_q(dku)` linear-twist hypothesis of Lemma 2.15 is broken by the parametrization mismatch.

**Deviation 2: hypothesis mismatch.** Lemma 2.15 requires `(q,3) = 1` and `d² | q | d³`. Our setting: `d = p, q = p^{r+1}`, so `d² = p²` and `d³ = p³`. The hypothesis `d² | q` requires `p² | p^{r+1}`, i.e., `r ≥ 1` (OK). The hypothesis `q | d³` requires `p^{r+1} | p³`, i.e., `r ≤ 2`. **Lemma 2.15 applies (literally) only at `r ≤ 2`.**

At `r = 2`, `q = p^3` and `d = p`; we have `d² = p² | q = p^3` ✓ and `q = p^3 | d³ = p^3` ✓. So Lemma 2.15 **does** apply at `r = 2`, modulo Deviation 1 (the parametrization issue).

At `r ≥ 3`, the Postnikov log `L_q(1+du)` has **cubic and higher terms in `u`** mod `q` (Lemma 2.13 only truncates to quadratic when `q | d³`). Lemma 2.15's evaluation relies on the quadratic truncation Lemma 2.13(2). **For `r ≥ 3`, GY's Lemma 2.15 does not apply.**

Lemma 2.16 (the weaker analogue for the regime `q | d²`) applies to `r = 1` (where `q = p², d = p`, so `q = p² = d²` ✓), but its output is simply `S_{q,d}(ψ, k) = q/d` on the supporting class — a trivial evaluation of magnitude `q/d = p`, not the Gauss-sum √q (since at this depth the Postnikov log truncates to **linear**, eq 2.13(1)). At `r = 1`, this matches our F̂ theorem (`FHAT_THEOREM_VERIFICATION_RESULTS.md` §4 Boundary 2): support size 1, magnitude `p` — exactly `q/d`.

## 3. The translation attempt, in detail

Let's force the question: at `r = 2` (where Lemma 2.15 hypotheses do apply), can we close the loop?

**At `r = 2`**:
- `q = p^3`, `d = p`, `N = p` (period = p^2, N = p^{r-1} = p).
- Our support has size `p^{r-1} = p`.
- Lemma 2.15 (if Deviation 1 is somehow absorbed) would give, for each support point `a`:
  `G[a] ≈ ε_q · √q · e_q(quadratic in a) · Jacobi`
  with `√q = p^{3/2}`.
- Compare empirical F̂ theorem: `|G[a]| = p^{(r+1)/2} = p^{3/2}` ✓ **(magnitudes match at r=2)**.

So at `r = 2`, even allowing for Deviation 1 absorbing the parametrization mismatch into the Jacobi-symbol / phase factor (which **would** require new work, since GY don't do this substitution), the *magnitudes* of GY's Lemma 2.15 evaluation and our F̂ theorem coincide.

**But the bilinear sum bound does not follow.** Even with `G[a]` evaluated in closed form, the bilinear is:
> `Σ = Σ_{a ≡ 1 (mod p)} 1̂(p·a) · p · G[a] = p · Σ_a 1̂(p·a) · e_q(quadratic in a) · (unit phase)`.

This is a **new** bilinear sum: `Σ_a 1̂(p·a) · e_q(Q(a))` where `Q(a)` is quadratic. To bound it by `≪ √(p^{r-1}) · √q = √N · √q` requires either:
- (i) a stationary-phase argument on the quadratic `Q(a)` — but `1̂(p·a)` is **not** a smooth amplitude, it is a Dirichlet kernel, so standard stationary phase doesn't apply directly,
- (ii) a Plancherel / Parseval re-expansion, but we've already used Plancherel once (that's how we got here), and re-applying gives back the primal sum,
- (iii) a Burgess-type estimate for the quadratic-phase character sum on the principal-unit coset — which **is the original ask**.

In other words: even at `r = 2` where Lemma 2.15 magnitudes match ours, Lemma 2.15 only gives a **closed-form expression for `F̂(p·a)`**, not a bound on the bilinear sum. We already have an equivalent closed form (Theorem 78.6 at q=3, generalized to F̂ theorem at family level). **Lemma 2.15 doesn't add new bilinear cancellation machinery.**

## 4. Adversarial checks

### A1 — Magnitude verify

Empirical: `|K|/√N ≈ 2.0` (constant across r=8..20, `r79b_S_partial_empirical.md`).

What would Lemma 2.15 give for `|K|/√N` if directly applied? **Lemma 2.15 gives no bilinear sum bound directly.** It gives `|F̂(p·a)| = p · √q = p^{(r+3)/2}` on support of size `p^{r-1}`, which we already have.

If we naively triangle-inequality the bilinear, plugging in `|F̂(p·a)| = p^{(r+3)/2}`:
> `|Σ| ≤ Σ_{a ≡ 1 (mod p)} |1̂(p·a)| · p^{(r+3)/2} ≤ p^{(r+3)/2} · Σ_a |1̂(p·a)|`.

Using Cauchy-Schwarz, `Σ_a |1̂(p·a)| ≤ √(p^{r-1}) · √(Σ_a |1̂(p·a)|²)`. By Plancherel, `Σ_a |1̂(p·a)|² ≤ q · N/q · (some factor) = N` or so; but actually the principal-unit-coset slice is `1/p` of `Z/p^r`, so `Σ_{a ≡ 1 mod p} |1̂(p·a)|² ≤ N · q/p · 1/p · p^r/q = N · p^{r-2}` ... this trivial routing gives a worse bound than the empirical √N.

The point: **plugging Lemma 2.15's evaluation into the bilinear gives nothing new beyond what we already have from the F̂ theorem.** Magnitude mismatch is not the issue; there is simply no additional structural information.

Verdict A1: **Lemma 2.15 magnitude doesn't conflict with our empirical** (since they're both `≈ √q` on each support point), **but it doesn't supply a bound on the bilinear sum either.** Inconclusive: doesn't fail magnitude verify, but doesn't deliver C either.

### A2 — Scope verify

Brief's check: Theorem 1.1 requires `d ≺ q ≼ d²`. Our: `d = p, q = p^{r+1}`. We need `q ≼ d²`, i.e., `p^{r+1} | p² · (unit power)`. Since `d` and `q` share prime factors (both are powers of `p`), `q ≼ d²` means `ν_p(q) ≤ 2·ν_p(d)`, i.e., `r+1 ≤ 2`, i.e., **`r ≤ 1`**.

So **Theorem 1.1 (Lemma 2.16 regime) applies only at `r = 1`**.

Theorem 1.2 requires `d² ≼ q ≼ d³`: `2 ≤ r+1 ≤ 3`, so **`r ∈ {1, 2}`**.

**Our empirical range is `r = 8..20`.** GY's framework, at `d = p`, is **out of scope for our depth.** This is a structural mismatch flagged in the brief's pre-reg.

What if we increase `d`? GY's `d` is a divisor of `q` with the depth-half condition. The depth-half condition `ν_p(d) ≥ ν_p(q)/2` for `q = p^{r+1}` gives `ν_p(d) ≥ (r+1)/2`, so `d = p^s` with `s ≥ (r+1)/2`. To map to our setting, what would `d = p^s` for `s ≥ (r+1)/2` correspond to? Our `1̂` lives at "depth-r" frequency: `1̂(p·a)` with `a ∈ Z/p^r`. So the "thin coset" in our setting would naturally be `d = p^r` (the depth where 1̂ is supported), not `d = p`.

**Re-translation with `d = p^r`:** Now `q = p^{r+1}, d = p^r, q/d = p, q/d² = p^{1-r}` — but `q/d² = p^{1-r}` is **not an integer for r ≥ 2**, violating Lemma 2.15's hypothesis `d² | q`. So Theorem 1.2's hypothesis fails. Theorem 1.1's hypothesis is `q ≼ d²`, i.e., `r+1 ≤ 2r`, i.e., `r ≥ 1` ✓. So Theorem 1.1 (Lemma 2.16 regime) applies with `d = p^r, q = p^{r+1}` for any `r ≥ 1`.

At `d = p^r`: Lemma 2.16 evaluates `S_{q,d}(ψ, k) = Σ_{u (mod q/d)} ψ(1+du) e_q(dku) = Σ_{u (mod p)} ψ(1+p^r · u) e_q(p^r · ku)`. This is a sum of length `q/d = p`, an extremely short sum.

The output (Lemma 2.16): `S_{q,d}(ψ, k) = q/d = p` if `k ≡ -a_ψ (mod q/d) = -a_ψ (mod p)`, else `0`. Magnitude `p`, not `√q`.

**This is much weaker than the empirical `√q` saturation.** Theorem 1.1's regime, at the natural `d = p^r` for our setting, gives Lemma 2.16 with a trivial magnitude `q/d = p`, not the Gauss-sum √q of Lemma 2.15. The √q in Theorem 1.1's `A` (eq 1.8: `A = ϕ(d)/d · √q · σ_0(|a_ψ|)/√|a_ψ|`) is the **net size after diagonalizing**, not the magnitude of `S_{q,d}` itself.

Verdict A2: **Hypothesis fit fails.** Lemma 2.15's regime (`d² | q | d³`) maps to `r ∈ {1, 2}` at `d = p`. Lemma 2.16's regime (`d | q | d²`) maps to `r = 1` at `d = p` or all `r ≥ 1` at `d = p^r`, but with trivial magnitude `q/d`. **Our `r ≥ 3` empirical range is not covered by GY's framework at any natural `d`.**

### A3 — Hidden infrastructure check

What does GY use that we **don't** have?

- **Approximate Functional Equation (Lemma 2.17, Iwaniec-Kowalski Theorem 5.3).** Converts `|L(1/2,χψ)|²` to a sum over `(m,n)` with weight `V(mn/q)`. We have no L-function setup, no AFE.
- **Orthogonality of characters mod `d`.** Used to detect `m ≡ ±n (mod d)`. We have **direct** access to the principal-unit coset (the F̂ support); no orthogonality detection needed.
- **Heath-Brown's hybrid bilinear bound (Lemma 2.2).** Used in §3.4 to handle the **balanced** range. Lives on a different object (`S(q; χ, h, n) = Σ_m χ(m+h)χ̄(m) e(mn/q)`, a Kloosterman-shape sum), and gives a sub-Weyl `√q · q_0^{-1/2}` saving in the bilinear, not a √N rate.
- **Postnikov formula (Lemma 2.14).** We have this implicitly via the Cochrane truncated p-adic log (`result_78_extended.md` Theorem 78.4).
- **Quadratic Gauss-sum evaluation (Lemma 2.4).** Standard, we can replicate.

**Is the AFE essential?** Yes — it's the entry point bringing the L-function moment to a sum `Σ ψ(m)ψ̄(n)/√(mn) · V(mn/q)`. The `1/√(mn)` amplitude is what makes the partial-summation argument in (1.9) yield `√q`. **Without AFE-style amplitudes, the (1.9) heuristic doesn't run.**

In our setting, `1̂(p·a)` has *no* amplitude `1/√...`; it's a Dirichlet kernel with values bounded by `min(N, |sin(πa/p^r)|^{-1})`. The (1.9) mechanism does not transfer to our bilinear.

Verdict A3: **AFE machinery is essential to GY's √q mechanism.** Removing the AFE (i.e., the `1/√(mn)` weight) removes the (1.9) partial-summation argument that produces the √q. Our setup has no AFE, so the mechanism does not transfer.

### A4 — Honest scope on extraction difficulty

The brief says: don't claim H_EXTRACTION_WORKS if it would require new lemma proofs.

What new work would be needed to **even attempt** to use GY's Lemma 2.15 for our bilinear?
1. **Parametrization conversion.** Re-derive Lemma 2.15 in the multiplicative-coordinate `(1+p)^u` rather than additive `1+pu`. This is more than relabeling: the linear twist `e_q(dku)` becomes a nonlinear-in-u phase. Re-doing the change of variables `u ↦ u + q/d²` and the quadratic Gauss-sum evaluation in the multiplicative coordinate is a non-trivial computation.
2. **Hypothesis extension to `r ≥ 3`.** Lemma 2.15's truncation Lemma 2.13(2) relies on `q | d³`. At `r ≥ 3`, the Postnikov log has cubic-and-higher terms. We have Theorem 78.4-78.6 (`result_78_extended.md`) that handles the cubic-and-higher case at p=3 specifically, via saddle-point + Hensel correction. Extending Lemma 2.15 to `r ≥ 3` would essentially mean **re-deriving Theorem 78.4-78.6 at family level for general p** — which is open.
3. **Bilinear bound itself.** Even if we had the closed-form evaluation of `G[a]` at family level (steps 1+2), we'd still need to bound the bilinear `Σ_a 1̂(p·a) · G[a] · e_q(quadratic_phase)`. This is the original ask — Lemma 2.15 doesn't bypass it.

**Verdict A4:** Each of steps 1, 2, 3 is non-trivial new mathematics. **Lemma 2.15 in the GY paper does not extract a bilinear bound for our setting.** What it does is **provide the structural template for a closed-form evaluation of F̂(p·a)** — and we already have that template at q=3 (Theorem 78.4-78.6) with a documented residual rigor gap at the family level (the equidistribution claim).

## 5. Disposition

**EXTRACTION_FAILS_STRUCTURAL.**

GY's secondary main term `A` of size `≈ √q` does not come from a bilinear-sum bound that translates to our setting. It comes from two combining sources:
1. The Gauss-sum magnitude `√q` of Lemma 2.15 (one length-`q/d` complete sum, quadratic Postnikov phase).
2. The Plancherel-duality `1/√q` prefactor from Poisson summation.

These compose with the V-weight Fourier-transform rapid decay and the divisor restriction `jn ≡ -a_ψ (mod d)` to give a deterministic main term, **not** a cancellation-based bilinear bound.

The mechanism diverges from our bilinear at the step where the (j,n)-summation is handled: GY rely on **Ŵ rapid decay** (a property of the smooth test function from AFE) to restrict the bilinear range to a divisor sum, then sum the residue contributions directly. We have **no smooth test function** — our `1̂(p·a)` is the Dirichlet kernel, which has slow decay (only `1/|a|`).

**Specific equation citation for divergence point:** Eq (3.11) in GY:
> `|Ŵ^±_{n,d}(x)| ≪_C (M/d) (1 + |x|M/d)^{-C}` for any `C > 0`.

This rapid-decay property of `Ŵ` (the smooth-amplitude Fourier transform) is what allows GY to truncate the `j`-sum at `|j| ≪ q^{1+ε}/M` with negligible error. **Our 1̂(p·a) does not have analogous rapid decay** — its Fourier transform is a Dirichlet kernel-of-Dirichlet-kernel which has only polynomial decay. The mechanism doesn't run.

## 6. What the structural parallel *does* tell us

Wilson's intuition that GY's `q^{1/2}` is structurally aligned with our √N is correct in **size**, but the parallel is between:
- **GY's `A` size**: a determined main term in a moment expansion, scaling as `ϕ(d)·√q/(d·√|a_ψ|)`, emerging from AFE + Poisson + Postnikov-Gauss saddle.
- **Our `|K|/√N`**: an empirical saturation constant, scaling as a **second moment** of `1̂·F̂`, emerging from Plancherel + Dirichlet-kernel concentration on principal-unit-coset frequencies.

Both involve √q at the "right" depth (`q = p^{r+1}` in our case). Both involve the **same Postnikov-Gauss-sum building block** (Lemma 2.15 ↔ Theorem 78.5/6 / F̂ theorem). But the **role** of √q is different:
- GY: √q is the **magnitude of a single Gauss sum** that becomes the **main term** after the bilinear range restricts to a divisor count.
- Us: √q is the **uniform magnitude of `F̂(p·a)` on support**, and the bilinear sum over a principal-unit-coset of size `p^{r-1}` saturates at `√N · √q` empirically — meaning `Σ_a |1̂(p·a)|² ≈ N` rather than `N²/3` (trivial).

These coincide in *size* (`√q · √N` vs `√q · √N`) because of the principal-unit-coset cardinality `p^{r-1} = N`. The **saving** in our bilinear is the analogue of *not* `A` but of a **cancellation among the `(a)` terms** — which is **what GY don't prove** (their `A` is the leading non-cancellation contribution, summed directly).

## 7. What would resolve this

Three potential paths, all open:
- **Path A.** Extend Theorem 78.4-78.6 to family level (general p) — handles our **closed-form** side of the equation. Doesn't directly close the bilinear bound but tightens our infrastructure.
- **Path B.** A dual Plancherel-side analogue of Milićević's F-class sub-Weyl bound — the gap flagged in `BURGESS_LITERATURE_FINDINGS.md` §5 row 4. GY doesn't supply this.
- **Path C.** A decoupling-style estimate for the principal-unit-coset DFT — `Σ_{a ∈ coset} |1̂(p·a)|² ≪ N` (instead of the trivial `N²/3`). This is **the** bilinear bound. None of GY/Milićević/Banks-Shparlinski/KMS/FKMS gives this directly.

GY 2023 is the **closest** in the literature (per `BURGESS_LITERATURE_FINDINGS.md`), but the closeness is at the level of *size of secondary term* and *Postnikov building block*, not at the level of bilinear-sum machinery. **The structural parallel is suggestive but not extractable.**
