# Pascadi 2025 — Direct translation attempt + adversarial checks

**Date:** 2026-05-11. Reader: extraction agent. Paper: arXiv:2511.08445.

## 1. The translation attempt

### 1.1 What we have

> `Σ_{a ≡ 1 mod p in Z/p^r}  1̂(p·a) · F̂(p·a)`,  `1̂(ξ) = Σ_{u=0}^{N-1} e_q(ξ·u)`,  `N = p^{r-1}`,  `q = p^{r+1}`,  `|F̂(p·a)| = p · √q` uniform on support.

### 1.2 What Pascadi assumes

> `Σ_{m ∈ I, n ∈ J, (mn, c) = 1}  α_m · β_n · S(am, n; c)`,  `|I|, |J| ≪ c^{1/2+o(1)}`, `α, β` arbitrary complex sequences, `S(m, n; c)` classical Kloosterman.

### 1.3 The attempted re-write

To even attempt to plug our object into Pascadi's framework, we need to identify three things:

- A modulus `c` ↔ our `q = p^{r+1}` — TENTATIVE.
- Two interval variables `m ∈ I, n ∈ J` ↔ ??? — **WE HAVE ONE COSET VARIABLE `a`**.
- A Kloosterman kernel `S(am, n; c)` ↔ ??? — **WE HAVE NO KLOOSTERMAN SUM ANYWHERE**.

There is no way to perform the identification. The closest possible attempt is:

**Attempt 1.3a.** Set `m = a`, `n = 0` (degenerate), `α_a = 1̂(p·a)`, β_0 = F̂(p·0) = 0 (since `0 ∉ supp(F̂)`). The Kloosterman factor `S(a·a', 0; c) = S(a·a', 0; c) = c · 1̂[(a·a',c) = 1] - 1` (Ramanujan sum at n=0), which has no relation to our F̂. Fails immediately.

**Attempt 1.3b.** Try to absorb `F̂(p·a) = p·√q·ψ(a)` into the weight: set `α_a := 1̂(p·a) · √q · p`, `β_n := ψ(n) · 1[n = a]` (constraint that `n` follows `a`). But this is a single-variable sum disguised as a two-variable sum with a delta-function constraint — Pascadi's mechanism cannot extract any saving from such a constraint because the §2.2 Fourier transform in `(m, n)` would just undo the disguise.

**Attempt 1.3c.** Insert a Kloosterman factor ourselves via Plancherel: write `1̂(p·a) = Σ_y (something involving S(p·a, y; q))`. The Kloosterman sum has its own Fourier-inversion representation, but `1̂` is the **additive Dirichlet kernel**, not a multiplicative character. There is no Plancherel identity that turns `1̂` into a Kloosterman sum — and even if there were a heavy machinery to do so, we'd then need to apply Pascadi's machine to a sum that exists in an entirely re-encoded form, and the re-encoding would not preserve the original cardinality structure.

**Attempt 1.3d.** Use Pascadi's amplifier `Γ_q(d)` on the abelian subgroup (Z/p^r)^× viewed inside SL₂(Z/qZ). The principal-unit coset `{a ≡ 1 (mod p)}` of (Z/p^r)^× sits inside the diagonal subgroup `D = {((a, 0),(0, a^{-1}))} ⊂ SL₂(Z/qZ)`. The non-trivial part of Pascadi's amplifier lives in the OFF-DIAGONAL elements of SL₂(Z/qZ), specifically in the upper-triangular generator `T = ((1,1),(0,1))`. Our sum has no off-diagonal content. The amplifier acts trivially on our sum.

### 1.4 Verdict on direct translation

**The translation does not even reach a starting line.** Pascadi requires:
- a two-variable bilinear sum over integer intervals (we have one coset variable),
- a Kloosterman kernel (we have no Kloosterman sum),
- a range `|I|, |J| ≪ c^{1/2 + o(1)}` (the analog of `N ≪ √q` fails for r ≥ 4),
- arithmetic embedded in the full SL₂(Z/cZ) group (our object lives in the diagonal abelian subgroup).

No notation correspondence consistent with the mechanism can be set up.

## 2. Adversarial checks (A1-A6 per pre-reg)

### A1. Smooth-weight obstruction (vs GY's Ŵ-rapid-decay)

**Test.** Does Pascadi's machinery require smooth weights with rapid Fourier decay?

**Result.** Pascadi explicitly allows arbitrary complex weights `(α_m), (β_n)` (§1, statement of Theorem 1.2: "complex entries"). The proof uses `‖α‖₂, ‖β‖₂` only via Cauchy-Schwarz against the Kloosterman matrix's spectral decomposition. **No smoothness required.**

**Comparison to GY.** GY's mechanism required Ŵ-rapid-decay to truncate the j-sum; Pascadi does not need this. **Pascadi RESOLVES the smooth-weight obstruction that killed GY.**

However: this is moot for our problem, because the obstruction in A4/A5 (below) blocks the translation before the weight question matters. We don't have weights `(α_m), (β_n)` to attach.

### A2. Scope obstruction (depth r ≥ 3)

**Test.** Does Pascadi's hypothesis range exclude `r = 8..20`?

**Result.** Pascadi Theorem 1.2 allows any composite c with factorization `c = d·d'·e` (d' | d, (d, e) = 1). For `c = p^r`, can take `d = p^{⌊r/2⌋}, d' = p^{r - ⌊r/2⌋}, e = 1`. **The theorem statement admits arbitrary prime powers.**

But: the saving for `c = p^r` is `(f / min(c, d²))^{1/6}` with `f ≍ p^{⌊r/2⌋}` (largest integer with `f² | c·d = p^{r + ⌊r/2⌋}`), and `min(c, d²) = p^{min(r, 2⌊r/2⌋)}`. For r ≥ 3, the saving exponent is roughly `-r/24` in c-units — non-trivial but **NOT √N saturation**.

**Comparison to GY.** GY's `d² | q | d³` excluded our `r ≥ 3` setting. Pascadi's hypothesis does not exclude r ≥ 3. **Pascadi RESOLVES the scope obstruction that killed GY.**

Again moot — A4/A5 kill the translation before scope matters.

### A3. AFE / hidden infrastructure obstruction

**Test.** Does Pascadi's saving secretly rest on amplitudes coming from elsewhere (L-function AFE, smooth amplitude, sheaf-theoretic input)?

**Result.** Pascadi's saving is purely **representation-theoretic + combinatorial**:
- §5: amplifier from normal-subgroup structure of `Γ_c(d) ⊲ SL₂(Z/cZ)`;
- §6: counting solutions to word equations `T^{h_1}S···T^{h_6}S ≡ ±I (mod p^k)`;
- §3.3 / Lemmas 5.2, 5.4: character bounds via Clifford theory.

**No AFE. No L-function. No sheaf-theoretic Weil/Deligne input.** §2 explicitly notes the method does NOT rely on Helfgott expansion, spectral gaps, ℓ-adic cohomology, or Ramanujan-Petersson. The saving is "purely combinatorial and character-theoretic, localized to properties of SL₂(Z/cZ) and the structure of its congruence kernels."

**Comparison to GY.** GY required the AFE to supply amplitudes `1/√(mn)` whose stationary-phase yielded `√q`. Pascadi does not. **Pascadi RESOLVES the AFE obstruction that killed GY.**

But the next obstruction kills the whole thing:

### A4. Cardinality vs amplitude × phase decomposition (the LOAD-BEARING check)

**Per GY agent's structural insight**: our `√N` is **a cardinality factor** from the principal-unit coset size `p^{r-1}`, combined with F̂'s magnitude `p·√q` at each support point. The decomposition is `cardinality × magnitude`. GY's `√q` was `amplitude × phase` (Gauss-sum magnitude × Plancherel prefactor).

**Test.** Does Pascadi's saving come from cardinality (matching ours) or from amplitude × phase (matching GY)?

**Result.** Pascadi's saving comes from **NEITHER**. It is:
- **Combinatorial counting** in SL₂(Z/p^k Z): the number of solutions to a word equation in the matrix group. The saving `p^{-1/6}` per prime is the ratio `|solutions mod p²| / (|solutions mod p²| expected from trivial counting)`.
- **Character orthogonality** at the level of irreducible representations of SL₂(Z/cZ).

Neither of these resembles a cardinality-of-coset saving. The amplifier `Σ_{ℓ ∈ Γ_c(d)} ρ'(ℓ)* ⊗ ρ(ℓ)` uses the NORMAL subgroup `Γ_c(d) ⊲ SL₂(Z/cZ)` non-trivially — and the abelian image of (Z/p^r)^× in SL₂(Z/cZ) (the diagonal subgroup) is **NOT a normal subgroup**. So Pascadi's amplifier cannot be re-deployed on our setting even at the level of "build an amplifier on (Z/p^r)^×."

**Crucial.** The abelian analog of Pascadi's amplifier is — per Pascadi's §2.3 explicit comment — **trivial**: "For abelian (Dirichlet) characters, the analogous amplifier weights `χ̄'(ℓ)χ(ℓ)`, which is trivial when only `{I}` is available." Our problem is abelian (the ambient group of summation is the abelian (Z/p^r)^×). The Pascadi amplifier has no non-trivial abelian analog. **A4 fails decisively.**

### A5. Magnitude verify

If we forced a translation (which A4 says we cannot), Pascadi gives a saving of `c^{1 - 1/12+o(1)} = q^{1 - 1/12}` for `c = q = p^{r+1}` (with the best choice of d). Our target is `N · √q = p^{r - 1 + (r+1)/2} = p^{(3r - 1)/2}`. In `c = q = p^{r+1}` units, target is `c^{(3r-1)/(2(r+1))}`. For r = 8: target `c^{23/18} ≈ c^{1.28}`. Pascadi gives `c^{1 - 1/12} ≈ c^{0.917}`. **Pascadi's bound is BELOW our target** — i.e., even if the translation worked, Pascadi's bound would be **WEAKER** than what we already have empirically `|K| ≲ √N ≪ N^{1/2}`.

Wait — re-check. Our target is `|Σ_a 1̂(p·a)·F̂(p·a)| ≤ C · N · √q`. The trivial bound (Cauchy-Schwarz + uniform F̂ magnitude) is `N · |F̂| · max|1̂| = N · p·√q · N = N² · p · √q = c · N · p / √q · √q = c · N` roughly. Pascadi's `c^{1 - 1/12}` gives a saving over `c` only; but the trivial bound for our sum is much larger than `c`. So Pascadi-style δ-saving over trivial doesn't reach our needed `√N · √q` rate.

**Order-of-magnitude verdict.** Even setting aside the structural mismatch, Pascadi's saving (δ-improvement over trivial) is the **wrong order of magnitude** for our needed rate (√N · √q saturation). The targets are not comparable.

### A6. Extraction difficulty / honest scope

**Test.** Even if we performed a non-trivial reformulation, would translation work?

**Result.** A "non-trivial reformulation" would require:
1. Embedding the cardinality-on-coset sum into a Kloosterman-matrix sum (no known mechanism — would require a new representation of our object that introduces SL₂(Z/cZ) content where there is none);
2. Generalizing Pascadi's non-abelian amplifier to a normal subgroup of an abelian group (impossible — every subgroup of an abelian group is normal, but Pascadi's amplifier becomes trivial in the abelian case per §2.3 explicit comment);
3. Or: proving a parallel saving via a completely different mechanism (a new theorem, not a translation).

**This is not "translation work." This is "do new mathematics in a different direction."** The literature scan's "method-shape match" is at the level of "both papers contain the words 'bilinear', 'amplification', 'composite moduli'." The methods themselves don't translate.

## 3. Decision per pre-reg

**EXTRACTION_FAILS_STRUCTURAL.**

Mapping to the three obstructions from GY:

| Obstruction | GY status | Pascadi status |
|---|---|---|
| (A1) Smooth-weight | KILLED GY | **Resolves** (Pascadi handles arbitrary weights) |
| (A2) Scope (r ≥ 3) | KILLED GY | **Resolves** (Pascadi admits all prime powers) |
| (A3) AFE / L-function | KILLED GY | **Resolves** (Pascadi has no L-function input) |
| (A4) Cardinality vs amplitude×phase | (load-bearing structural insight) | **NEW OBSTRUCTION**: Pascadi's saving is combinatorial-counting on SL₂(Z/cZ), neither cardinality nor amplitude×phase |
| (A5) Magnitude | n/a | Saving wrong order of magnitude even if translation worked |
| (Object shape) | Cardinality vs amplitude | **MORE FUNDAMENTAL**: Pascadi requires Kloosterman-matrix on two free variables; we have one coset variable with no Kloosterman sum |

Pascadi resolves **all three GY obstructions** (A1, A2, A3) but introduces a **new fundamental obstruction**: the object Pascadi bounds is a Kloosterman bilinear over two interval variables, and the saving is harvested from non-abelian structure of SL₂(Z/cZ). Our object is a one-variable coset sum with no Kloosterman content and lives in the abelian quotient (Z/p^r)^×.

**The lit-scan "method-shape match" claim is false** at the level of mechanism. The match is at the level of keywords ("bilinear," "amplification," "composite moduli," "Fourier"); the mechanism does not transfer.

## 4. What this rules out / leaves open

**Rules out.** Direct or nearly-direct translation of Pascadi's method to our problem.

**Does NOT rule out** (and is consistent with prior findings):
- Some entirely different abelian-amplifier construction tailored to (Z/p^r)^× and the principal-unit coset. The DFI / Heath-Brown abelian amplification scheme is the natural pre-cursor, but already known not to give √N saturation in the depth-r regime.
- A new bilinear bound proven by entirely different machinery (e.g., a direct evaluation of `Σ_a 1̂(p·a)·G(a)` using the explicit Cochrane-log structure of `G(a) = √q · e_q(P_a(s*))` per Theorem 78.6 — this is a "do the computation" path, not a "apply external theorem" path).
- Future work generalizing Pascadi's non-abelian amplifier to non-normal subgroups or to non-Kloosterman kernels — speculative and not visible in the published method.

The published Pascadi machine, as written, does not extract.
