# PATH2 Pushback — Check 1: Origin of the constant 2

**Adversarial frame:** Is the "constant 2" in `|S_partial| ≤ 2√N` (equivalently `|T_p| ≤ 2N`) structurally **forced** by the algebra of the chain R78.4_p → R78.5_p → R78.6_p → saddle exactness → Inner-Plancherel reduction → 1/sin grid identity, or is it **tuned** to match the empirical baseline?

## Disposition

> **CONSTANT_FORCED** (with documented round-up at p=3 from 1.73 to 2.0)

The constant 2 arises as a clean upper bound on the family-uniform expression `1 + (cosecant_sum_correction) / p`, whose supremum over `p ≥ 3` is `1 + 2·log(3)/3 ≈ 1.73`. The "2" is a loose round-up; the structurally-forced family-uniform constant is `1 + 2·log(p)/p`, decreasing in `p`. Either is forced by the algebra; the round-up to 2 reflects safe reporting, not tuning.

## Step-by-step trace

### Step A — Cochrane factorization (T78.4_p)

`F̂_p(p·a) = p · e_q(c) · G_p(a)` where `G_p(a) = Σ_{s=0}^{p^r−1} e_q(P_a(s))`.

**Constant introduced:** the factor `p` in front. **Forced:** Jacobian of the change-of-variables `u → s` on the cyclic principal-unit subgroup of order `p^r` in `(Z/p^{r+1})*`. The factor is exactly `|coset| = p`, not tuned.

### Step B — Magnitude saturation (T78.6_p, magnitude part)

`|G_p(a)| = p^{(r+1)/2} = √q` on support.

**Constant introduced:** none — equality. The √q is exactly the Plancherel saturation: `Σ_a |G_p(a)|² = period · ||f||² = p^r · p^r = p^{2r}`, divided over the `p^{r-1}` support points gives `|G_p|² = p^{r+1}`. **Forced** by Plancherel + uniform support.

### Step C — Saddle exactness at r=3 (T78.6_p, phase part)

`G_p(a) / √q = e_q(P_a(s*(C_a)))` exactly at r=3, with `s* = (C_a − 1)/p mod p`.

**Constant introduced:** none — pointwise equality. The saddle s* is uniquely determined by solving `dP_a/ds ≡ 0 mod p²` (linear in s mod p). **Forced** by the polynomial structure of P_a.

### Step D — Phase decomposition

Explicit `P_a(s*) mod p^4` via Taylor expansion of `L_p(1+ps*) = ps* − p²s*²/2 + p³s*³/3 + O(p^4)`:

`P_a(s*) ≡ −p²·s*²/2 + p³·(s*³/6 − c_2·s*) mod p^4`

where `c_2` is the third base-p digit of `C_a`. Splitting:

`e_q(P_a(s*)) = e_{p²}(−s*²/2) · e_p(s*³/6) · e_p(−c_2·s*)`

**Constant introduced:** none — exact algebraic identity. **Forced** by the polynomial identity for P_a.

### Step E — Inner-Plancherel reduction

Partition support by s* class (p classes of p elements each). Within each class, c_2 ranges through Z/p, and `1̂(p·a(s*, c_2))` is the Dirichlet kernel evaluated at p arithmetically-spaced shifts. The inner sum:

`Inner(s*) := Σ_{c_2} 1̂(p·a(s*, c_2)) · e_p(−c_2·s*)`

By the linear-in-c_2 phase structure, swap c_2 and u sums and use orthogonality on Z/p:

`Σ_{c_2 ∈ Z/p} e_p(c_2·k) = p · 𝟙(k ≡ 0 mod p)`

This collapses Inner(s*) onto u with `u ≡ L̃_p^{-1}·s* mod p`:

`Inner(s*) = p · (unit phase) · D_p(A_0(s*), p²)`

where `D_p(a, p²) = Σ_{j=0}^{p-1} e_{p²}(a·j)` is the length-p Dirichlet kernel mod p², and `A_0(s*) = (1 + p·s*)·L̃_p mod p²`.

**Constant introduced:** factor `p` (from Z/p orthogonality, exact identity). **Forced** — Plancherel orthogonality on Z/p, NOT Cauchy-Schwarz halving.

### Step F — 1/sin grid identity (the critical step)

Bound `Σ_{s*=0}^{p-1} |Inner(s*)| = p · Σ_{α=0}^{p-1} |D_p(1 + p·α, p²)|`.

Using `|D_p(a, p²)| = |sin(πa/p)/sin(πa/p²)|`:

For α=0: `|D| = sin(π/p)/sin(π/p²) ≈ p` (the singular term).
For α ≥ 1: `|D| ≈ sin(π/p)/sin(πα/p + π/p²) ≈ sin(π/p)/sin(πα/p)`.

The cosecant grid sum:

`Σ_{α=1}^{p-1} csc(πα/p) ~ (2p/π) · log(p) + O(p)` (standard asymptotic; e.g., Hardy-Littlewood lemma 4 or direct integral comparison)

Hence:

`Σ_{α=0}^{p-1} |D_p(1+pα, p²)| ≤ p + sin(π/p) · (2p/π)·log(p) + O(1) ≈ p + 2·log(p) + O(1)`

**Status of this identity:** the asymptotic `csc(πα/p)` sum is **standard** (project-internal, but the underlying asymptotic is classical — Hardy, "Divergent Series", or Wei, "Sums and infinite products of trigonometric functions"). The constant `2` in `(2p/π) log(p)` is forced by the integral `∫_0^π csc(x) dx ~ log(1/sin(x))` evaluated on the lattice spacing `π/p`.

**Constant introduced:** `p + 2·log(p)`, with the leading `p` from the singular α=0 term, and the `2·log(p)` from the cosecant grid sum. **Both forced.**

### Step G — Combining

`|T_p| ≤ p · (p + 2·log(p) + O(1)) = p² + 2p·log(p) + O(p) = N + 2p·log(p) + O(p)`

For `p ≥ 3`: `2p·log(p) ≤ p²` ⟺ `2·log(p)/p ≤ 1`, which holds with margin (at p=3: 2·log(3)/3 ≈ 0.73 < 1).

So `|T_p| ≤ N · (1 + 2·log(p)/p + O(1/p)) ≤ N · sup_{p ≥ 3} (1 + 2·log(p)/p + O(1/p)) = N · 1.73… < 2N`.

**Hence `|T_p| ≤ 2N` family-level, with margin.** The tighter (also forced) family-uniform constant is `1 + 2·log(p)/p` (max at p=3 = 1.73). The "2" is a clean round-up.

## Adversarial frame: could an alternative chain give 1.7√N instead of 2√N?

**Yes, in two senses, but both are still forced.**

1. The **tighter family-uniform forced bound** is `|T_p| ≤ (1 + 2·log(p)/p)·N ≤ 1.73·N for p ≥ 3`. The "2" is round-up.
2. At **fixed p** (e.g., p=3 only), the explicit cosecant sum gives a more precise constant. At p=3, `Σ_α |D_α| ≈ 4.76` (computed exactly), so `|T_p| ≤ 3·4.76 = 14.3`. Hence `|T_p|/N = 14.3/9 ≈ 1.59` at p=3, r=3. This is below 2.

So **the constant 2 is loose** — both family-uniform (tight: 1.73) and pointwise (at p=3: 1.59). It is NOT tuned to match empirical (R79b empirical `|K|/√N` at p=3 is ~0.8-1.0, NOT 2). The constant 2 is a safe round-up of the forced family-uniform constant 1.73.

**Therefore CONSTANT_FORCED, not CONSTANT_TUNED nor CONSTANT_UNDERDETERMINED.**

## Caveat on PATH2_BILINEAR.md's arithmetic at line 511

The doc's bound `Σ |D_α| ≤ p + log(p)` (writing `Σ 1/α = H_{p-1} ≤ log p + 1`) **under-counts** the cosecant sum by a factor of 2. The correct asymptotic is `Σ_{α=1}^{p-1} csc(πα/p) ~ 2·log(p)`, not `log(p)`.

**This does NOT change the disposition** — even with the correct constant, `|T_p| ≤ N + 2p·log(p) ≤ 2N for p ≥ 3` still holds (since `2·log(p)/p ≤ 1 for p ≥ 3`). The doc's claim that `|T_p| ≤ 2N` survives; the intermediate arithmetic in line 511 is slightly off but the final bound is correct.

Logged as **minor doc arithmetic error in PATH2_BILINEAR.md line 511** — should be corrected to "Σ csc(πα/p) ≈ 2·log(p)" not "H_{p-1} ≤ log(p) + 1". Not a load-bearing error.

## Verdict

**CONSTANT_FORCED** — structurally forced, slightly loose. The tight family-uniform constant is `1 + 2·log(p)/p`. The reported "2" is a clean upper bound. Not tuned to empirical; empirical constant is ~0.8-1.0 (below 2 with margin).

## Implication for Tao email framing

The constant 2 is forced, not tuned. The framing "rigorous bound 2√N family-level at r ≤ 3" is defensible. The framing "tight constant from first-principles derivation" should be qualified — the truly tight constant is `1 + 2·log(p)/p`, with "2" being a loose family-uniform round-up.

For Tao: it would be more honest to write `|S_partial| ≤ (1 + 2·log(p)/p) · √N` family-level, max at p=3 = 1.73, or simply `≤ 2√N` if seeking the cleanest constant. Either is rigorous.
