# Pascadi 2025 (arXiv:2511.08445) — Translation table

**Date:** 2026-05-11. Reader: extraction agent. Paper: Pascadi, "Non-abelian amplification and bilinear forms with Kloosterman sums," arXiv:2511.08445 (Nov 2025).

## 1. The Pascadi setting (their notation)

### 1.1 The bilinear sum

Theorem 1.1 / 1.2 of Pascadi bounds

> `B(α, β; a, c) := Σ_{m ∈ I, n ∈ J, (mn, c) = 1} α_m · β_n · S(am, n; c)`

where:
- `c ∈ Z₊` is the modulus (composite, with factorization c = d·d'·e, d' | d, (d, e) = 1)
- `a ∈ (Z/cZ)^×` is a fixed twist
- `I, J ⊂ Z` are integer intervals with `|I|, |J| ≪ c^{1/2+o(1)}`
- `(α_m), (β_n)` are arbitrary complex sequences (no smoothness)
- `S(m, n; c)` is the **classical Kloosterman sum**
> `S(m, n; c) := Σ_{x ∈ (Z/cZ)^×} e_c(m·x + n·x̄)`,  `x · x̄ ≡ 1 (mod c)`

### 1.2 The main bound (Theorem 1.2)

> `|B(α, β; a, c)| ≪ ‖α‖₂ · ‖β‖₂ · c^{1+o(1)} · (f / min(c, d²))^{1/6}`

where `f` = largest integer with `f² | c·d`.

### 1.3 The amplification mechanism (§5)

For each irreducible representation ρ' of SL₂(Z/cZ), the **amplifier**

> `A(χ') := ‖ Σ_{ℓ ∈ Γ_c(d)} ρ'(ℓ)* ⊗ ρ(ℓ) ‖²_{S₂}`

is built over the congruence kernel `Γ_c(d) := ker(SL₂(Z/cZ) → SL₂(Z/dZ))`, a **normal subgroup**. Saving comes from combinatorial counting of solutions to `T^{h_1} S T^{h_2} S ⋯ T^{h_6} S ≡ ±I (mod p^k)` in SL₂(Z/dZ).

## 2. Our setting

### 2.1 The bilinear sum (Form B from PRECISE_ASK.md)

> `K(r, c=1, m=0) = (1/q) · Σ_{ξ ∈ supp(F̂)} 1̂(-ξ) · F̂(ξ)`,
> `supp(F̂) = { p·a : a ∈ Z/p^r, a ≡ 1 (mod p) }`,  `|supp| = p^{r-1} = N`.

After Theorem 78.4 closed form:

> `K(r, 1, 0) = (3 · e_q(1) / √q) · Σ_{a ≡ 1 mod p in Z/p^r} 1̂(-p·a) · ψ(a)`

with `q = p^{r+1}`, `ψ(a) := G(a)/√q` unit-magnitude, and

> `1̂(ξ) := Σ_{u=0}^{N-1} e_q(ξ u)`  (length-N Dirichlet kernel, ξ ∈ Z/q).

### 2.2 The bound needed

> `|Σ_{a ≡ 1 mod p in Z/p^r} 1̂(p·a) · F̂(p·a)| ≤ C · N · √q`, explicit C.

## 3. Notation correspondence — attempted side-by-side

| Pascadi | Our setting | Match? |
|---|---|---|
| Modulus `c` (composite, c = d·d'·e) | `q = p^{r+1}` (prime power) | **Partial**: c = p^r included in Pascadi only via d ≍ p^{⌊r/2⌋}; saving degrades as r grows |
| Variables `m ∈ I, n ∈ J` (intervals in Z) | `a ∈ Z/p^r` subject to `a ≡ 1 (mod p)` (coset in Z/p^r) | **MISMATCH (load-bearing)**: intervals in Z vs. principal-unit coset in Z/p^r |
| Bilinear summand `S(am, n; c)` (Kloosterman) | `F̂(p·a)` (single F̂ value, **uniform magnitude p·√q on support**, verified theorem) | **MISMATCH**: Kloosterman across two free integer variables vs. one closed-form Gauss-sum-magnitude value |
| Weights `α_m, β_n` (arbitrary, two variables) | `1̂(p·a)` (Dirichlet kernel — **already a closed-form evaluation**, single variable) | **MISMATCH**: two free weight sequences vs. one explicit Dirichlet kernel |
| Saving via amplification over `Γ_c(d) ⊲ SL₂(Z/cZ)` (non-abelian normal subgroup) | Saving needed for cardinality factor `N = p^{r-1}` (size of principal-unit coset in (Z/q)^×) | **STRUCTURAL MISMATCH**: their saving lives in SL₂(Z/cZ); ours lives in the abelian quotient (Z/p^r)^× |
| Hypothesis `|I|, |J| ≪ c^{1/2+o(1)}` | `N = p^{r-1} ≈ q^{(r-1)/(r+1)} = q^{1 - 2/(r+1)}` — for r = 8..20 this is `q^{7/9} ... q^{19/21}`, which is FAR ABOVE `q^{1/2}` | **MISMATCH**: our range exceeds Pascadi's regime by a power |

## 4. Side-by-side of "bilinear sum shapes"

### 4.1 Pascadi's primal shape

```
Σ_{m ∈ [1, M]}  Σ_{n ∈ [1, N]}  α_m · β_n · S(am, n; c)
   |__________________|              |________________|
   product of two interval sums       Kloosterman across BOTH variables m, n
```

Structure: **two free interval variables, summand depends on both via the Kloosterman sum S(am, n; c)**. The Kloosterman sum is the "non-trivial" object; the weights α, β are arbitrary.

### 4.2 Our shape

```
Σ_{a ∈ {a ≡ 1 mod p in Z/p^r}}  1̂(p·a) · F̂(p·a)
       |______________________|     |_____| |_____|
       SINGLE variable a              Dirichlet  closed-form
       on a COSET of (Z/p^r)*         kernel    Gauss-sum-magnitude
                                      (length N)  ψ(a) · √q · p
```

Structure: **single coset variable a, summand is a product of two closed-form objects** (1̂ explicitly an oscillating Dirichlet kernel; F̂ already evaluated to be p·√q · ψ(a) with ψ unit-magnitude). The cancellation we need is among the **phases ψ(a) of the F̂ side weighted by 1̂**, summed over a.

## 5. Where the parallel breaks at level 0 (before any proof reading)

The Pascadi machine takes **two arbitrary weight sequences (α_m), (β_n) attached to two free interval variables m, n** and bounds their bilinear interaction through the Kloosterman sum `S(am, n; c)`. The non-abelian amplifier `Σ_{ℓ ∈ Γ_c(d)} ρ'(ℓ)* ⊗ ρ(ℓ)` extracts a saving from the GROUP STRUCTURE of SL₂(Z/cZ) on the Kloosterman matrix `K = (S(m,n;c))_{m,n}`.

Our object is a **single-variable sum on a coset** where one factor (F̂) is already evaluated to a unit-modulus phase times a deterministic magnitude. We have:
- **NO Kloosterman sum** anywhere in our object;
- **NO matrix structure** of the form `(K_{m,n}) = (S(am, n; c))` whose group-theoretic Fourier transform Pascadi exploits;
- **NO two free interval variables** to amplify against each other;
- **The cancellation we need is among phases ψ(a) on a coset** of size N in the abelian group (Z/p^r)^×, NOT among Kloosterman sums on a matrix indexed by SL₂(Z/cZ).

## 6. Scope-of-validity check

Pascadi's Theorem 1.1 requires:
- `M, N ≪ c^{1/2 + o(1)}`. **Our range exceeds this.** With `c ↔ q = p^{r+1}` and "cardinality of our sum" `N = p^{r-1}` (the coset size), the analog is asking whether `N ≪ q^{1/2}`, i.e. `p^{r-1} ≪ p^{(r+1)/2}`, i.e. `r - 1 ≤ (r+1)/2`, i.e. **r ≤ 3**.
- For our empirical range `r = 8..20`, the analog "interval length" `N` lives far above Pascadi's allowed regime by a power of p.
- Section 2.5: extension beyond `c^{3/8+o(1)}` requires a harder counting problem and the saving degrades.

Pascadi's Theorem 1.2 with c = p^r (r ≥ 3): optimal choice `d ≍ p^{⌊r/2⌋}` gives `f ≍ p^{⌊r/2⌋}`. Saving `(f / d²)^{1/6} = p^{-⌊r/2⌋/6}` — non-trivial **provided c can be factored as `c = d·d'·e` with `(d, e) = 1`**. For `c = p^r` (a pure prime power) the only available coprime decomposition is `d = p^a`, `e = 1` (and `d' = p^{r-a}` divides `d` only if `r - a ≤ a`, i.e. `a ≥ ⌈r/2⌉`). So `e = 1` always — but Pascadi's Theorem 1.2 hypothesis `(d, e) = 1` is trivially satisfied with `e = 1`. So in principle the theorem applies. But the **saving is then `p^{-⌊r/2⌋/6}`, NOT `√N`**.

Even if all other matches were perfect, the saving Pascadi delivers is δ-improvement over `c^{1+o(1)} · ‖α‖₂·‖β‖₂`, which is the **trivial bound** `M·N`. So Pascadi outputs `c · p^{-r/12}`. Our target is `N · √q = c^{1 - 2/(r+1)} · q^{1/2}`. For r = 8..20 the Pascadi exponent `1 - r/12` is FAR larger than our target exponent (1/2 + (r-1)/(r+1) in c-units) — Pascadi gives the trivial-side bound; we need the cardinality-side bound.

## 7. Summary

Surface match: the abstract mentions "bilinear forms," "Kloosterman sums," "composite moduli," "non-abelian amplification." All four are real, all four point at the right neighborhood. But at the level of the **object being bounded**:

- **Pascadi's object**: `Σ_m Σ_n α_m β_n · S(am, n; c)` — two free interval variables, Kloosterman sum across both, arbitrary weights.
- **Our object**: `Σ_a 1̂(p·a) · F̂(p·a)` on a coset — one free coset variable, no Kloosterman sum, one factor closed-form Gauss-sum-magnitude, the other a length-N Dirichlet kernel.

These are different bilinear-sum shapes. The translation question is whether Pascadi's amplifier mechanism can be re-deployed on our shape — examined in `PASCADI_MECHANISM.md` and tested in `PASCADI_TRANSLATION_ATTEMPT.md`.
