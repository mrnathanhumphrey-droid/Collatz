# Pascadi 2025 — Mechanism trace

**Date:** 2026-05-11. Reader: extraction agent. Paper: arXiv:2511.08445.

## 1. Where does the saving come from?

The `c^{-1/12}` saving in Pascadi's Theorem 1.1/1.2 emerges from a **single combinatorial-counting bound** on solutions to a word equation in SL₂(Z/p^k Z). The chain is:

1. **§2.2 — Fourier-transform the Kloosterman matrix.**
   Pascadi treats `K = (S(m, n; c))_{m,n}` as a c×c matrix (or as a kernel on (Z/cZ)^×) and applies a unitary Fourier transform in `m, n`. The result re-expresses `K` in terms of matrices
   > `ρ_c(T^{h_1} S T^{h_2})`
   where `T = ((1,1),(0,1))`, `S = ((0,-1),(1,0))` are the standard generators of SL₂(Z/cZ), and `ρ_c` is a representation. The Kloosterman sum then appears as a **matrix coefficient on SL₂(Z/cZ)**.

2. **§2.3 — Build the non-abelian amplifier.**
   Pick a normal subgroup `Γ_c(d) = ker(SL₂(Z/cZ) → SL₂(Z/dZ))` and form the amplifier
   > `A(χ') := ‖ Σ_{ℓ ∈ Γ_c(d)} ρ'(ℓ)* ⊗ ρ(ℓ) ‖²_{S₂}`
   for each irreducible representation `ρ'` of SL₂(Z/cZ). Pascadi states this is "to the best of our knowledge the first instance of such a construction" in the non-abelian setting.

3. **§3.3 / Lemma 3.8 — Character orthogonality at the normal subgroup.**
   Use orthogonality of characters of the irreducible representations of SL₂(Z/cZ) restricted to `Γ_c(d)`, combined with Clifford theory (Lemma 5.2 / 5.4), to convert `A(χ')` into a count of group-theoretic solutions.

4. **§2.4 + §6 — The combinatorial counting bound.**
   The saving is reduced to counting solutions to the word equation
   > `T^{h_1} S T^{h_2} S T^{h_3} S T^{h_4} S T^{h_5} S T^{h_6} S ≡ ±I (mod p^k)`
   with `|h_i| ≤ H ≈ √c`. Elementary counting shows this is `O(p²)` mod `p²` and `O(p³)` mod `p` (for `c = p²`). The ratio `O(p²) / O(p^{6·k - ?})` extracts the **`p^{-1/6}` per-prime saving**, which after raising to the 6th-moment threshold becomes `c^{-1/12}` for `c ∈ {p², pq}`.

5. **§2.5 / Lemma 5.5 / 5.6 — Recombination.**
   The final bound is assembled from `A(χ')`-amplified character bounds via a large-sieve / Cauchy-Schwarz step, producing
   > `Σ_m Σ_n α_m β_n S(am, n; c) ≪ ‖α‖₂ · ‖β‖₂ · c^{1+o(1)} · (f / min(c, d²))^{1/6}`.

## 2. What input does the mechanism require?

The mechanism takes as input:

**I1.** The Kloosterman matrix `K = (S(am, n; c))` indexed by `m, n ∈ Z/cZ` (or by intervals after restriction). It is **essential** that the object is the c×c Kloosterman matrix, because the §2.2 Fourier transform IS the SL₂(Z/cZ) harmonic analysis on this matrix.

**I2.** **Two free interval variables** `m, n` over which to apply the unitary Fourier transform. The amplifier of §2.3 builds matrix coefficients `Σ_{ℓ ∈ Γ_c(d)} ρ'(ℓ)* ⊗ ρ(ℓ)`, which is a tensor product matching the two-variable bilinear shape `α_m β_n`. Without two free variables, there is no tensor product to build.

**I3.** **Range hypothesis** `|I|, |J| ≪ c^{1/2 + o(1)}` (the Polya-Vinogradov threshold). The Fourier transform in `m, n` reorganizes the bilinear sum into a sum over `(h_1, h_2)` of size `H ≈ √c`, after which the §6 combinatorial counting works. Outside this range, §2.5 notes "one must use a smaller value of q [the moment order]; one would then need to solve a counting problem with few variables."

**I4.** **A coprime factorization c = d·d'·e** with `d' | d`, `(d, e) = 1`. Pure prime-power `c = p^r` is technically admissible (with `e = 1`), but the saving exponent then degrades: `(f / d²)^{1/6}` with `f ≍ d` for prime powers, giving `d^{-1/6}` per choice of `d`.

## 3. What does the mechanism PRODUCE?

The output is

> `|B(α, β; a, c)| ≪ ‖α‖₂ · ‖β‖₂ · c^{1+o(1)} · (f / min(c, d²))^{1/6}`.

This is a **δ-saving over the trivial bound** `‖α‖₂·‖β‖₂·c` (which is what you'd get from `|S(am, n; c)| ≤ √c` Weil + Cauchy-Schwarz). It is **not** a √N saturation bound. It is a Burgess-style sub-convex improvement over the trivial bound, exactly in the spirit of the abstract: "beyond Pólya-Vinogradov."

For `c = p²`, `d = p`, this gives `c^{1 - 1/12}`. For `c = p^r` with `r ≥ 3` and `d = p^{⌊r/2⌋}`:
- `f ≍ p^{⌊r/2⌋}` (since `f² | c·d = p^{r + ⌊r/2⌋}`, so `f ≤ p^{(r + ⌊r/2⌋)/2}`; for r even `f = p^{3r/4}`, for r odd `f ≈ p^{(3r-1)/4}`).
- `min(c, d²) = min(p^r, p^{2⌊r/2⌋}) ≈ p^r` for r even, `p^{r-1}` for r odd.
- Saving exponent in c-units: roughly `c^{1 - r/24}` (worse as r grows; for r = 8..20 this is `c^{1 - r/24}`, marginally beating trivial).

## 4. What does the mechanism NOT do?

**N1.** It does **not** bound bilinear sums of the form `Σ_a f(a) · g(a)` with **one free variable**. The whole architecture is built on a c×c Kloosterman matrix indexed by two variables.

**N2.** It does **not** produce √N-saturation with explicit constant. It produces a δ-improvement over the trivial `M·N` bound — a different shape of bound entirely.

**N3.** It does **not** handle sums where the variable is restricted to a small coset of (Z/p^r)^×. The Fourier-transform step in §2.2 sums over all `m ∈ I, n ∈ J` (intervals in Z), and the amplifier in §2.3 is built on the full ambient group SL₂(Z/cZ).

**N4.** It does **not** handle inputs where one of the two "variables" has been already evaluated to a closed-form value (like our `F̂(p·a)`). The whole mechanism is about producing cancellation between two arbitrary weight sequences via Kloosterman matrix structure — if one factor is already a deterministic magnitude `p·√q`, there is no "non-abelian Fourier transform on this factor" to perform; the factor lives in the abelian group (Z/p^r)^× and the SL₂(Z/cZ) Fourier transform has nothing to act on.

## 5. The gap

Pascadi's mechanism is a method for breaking the Polya-Vinogradov barrier on bilinear sums whose **arithmetic core is the classical Kloosterman sum on a composite modulus**, where the saving is harvested from the **non-abelian structure of SL₂(Z/cZ)** acting on the Kloosterman kernel.

Our problem has:
- **No Kloosterman sum** in the object;
- **No two-variable structure** to Fourier-transform;
- **A different ambient group** (the abelian quotient (Z/p^r)^×, not SL₂(Z/cZ));
- **One factor (F̂) already evaluated to closed form** with no internal arithmetic to amplify;
- **The cancellation needed is among phases on a coset**, not among matrix coefficients of a representation.

This is the "method-shape" mismatch in plain terms: Pascadi's amplifier acts on **objects in SL₂(Z/cZ)**, ours lives in **(Z/p^r)^×**, and the non-abelian piece is exactly what makes Pascadi's saving work (per §2.3 explicitly: "the analogous amplifier weights `χ̄'(ℓ)χ(ℓ)`, which is trivial when only `{I}` is available" — i.e. the abelian case is empty).

The combinatorial-counting heart of the proof (§6: counting `T^{h_1}S···T^{h_6}S ≡ ±I (mod p^k)`) is a count **in SL₂(Z/p^k Z)**. The only SL₂-flavor present in our setting is the trivial fact that (Z/p^r)^× embeds in SL₂(Z/p^r Z) as the diagonal subgroup `{((a, 0),(0, a^{-1}))}` — but our object lives entirely on this abelian subgroup, never visits the off-diagonal entries, and so does not couple to the non-abelian saving Pascadi extracts.
