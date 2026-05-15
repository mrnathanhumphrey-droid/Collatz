# ADELIC_I — Kontorovich 2014 (levels of distribution for affine sieve)

**Source:** C:/tmp/adelic/Kontorovich_2014_Levels_Distribution_Affine_Sieve.txt. Survey paper.

## Verbatim definition (Levels of Distribution, p. ~7)

> "Definition: Level of Distribution (for primes in progressions). We will say that the primes have a level of distribution Q if, for all A < ∞,
>    Σ_{q < Q} max_{(a,q)=1} |π(x; a, q) − Li(x)/φ(q)| = O_A(x / log^A x).
> When Q can be taken as large as x^{θ − ε} for some θ > 0, we call θ an exponent of distribution for the primes."

## Verbatim Bombieri-Vinogradov (p. ~7)

> "Bombieri-Vinogradov Theorem (1965). The primes have exponent of distribution θ = 1/2."

## Verbatim Affine Sieve setup (p. ~14)

> "Suppose we have an infinite set of natural numbers S ⊂ ℕ and wish to prove the existence and abundance of primes or R-almost-primes in S. Roughly speaking, all that is needed is an appropriate analogue [of Bombieri-Vinogradov on S]. In particular, suppose that S is fairly well distributed on average among multiples of q, in the sense that
>    #{n ∈ S ∩ [1, x] : n ≡ 0 (q)} = (1/q) #{S ∩ [1, x]} + r_q,
> (or perhaps with 1/q in (3.3) replaced by some analytically similar function like 1/φ(q)), where the errors are controlled by Σ_{q < Q} |r_q| = o(#{S ∩ [1, x]}), for some Q.
> Then the sieve technology (again very roughly) tells us that if Q can be taken as large as a power of x, say Q = x^{θ − ε} for some exponent of distribution θ > 0, then S contains R-almost-primes, with R = ⌈1/θ + ε⌉."

## Hypotheses isolated

- **h1 (SET S):** S ⊂ ℕ a thin orbit of a finitely generated semigroup Γ ⊂ GL_n(ℚ) acting on a base point.
- **h2 (DISTRIBUTION IN PROGRESSIONS):** S equidistributes on average in residue classes mod q with explicit error control.
- **h3 (EXPONENT θ > 0):** Such control holds up to Q = x^{θ − ε}.
- **CONCLUSION:** S contains R-almost-primes with R = ⌈1/θ + ε⌉.

## Hypothesis × input check

| Hyp | Syracuse |
|---|---|
| h1 (thin orbit of Γ ⊂ GL_n(ℚ)) | Could view (Z/3^n)* trajectories under Tao as orbits, but Tao is not a group action — it's an iterated *map* on (Z/3^n)*, not a *linear action* of a group. The map r ↦ (3r+1)/2^v is not the action of a fixed element of GL_n(ℚ); it depends on v. |
| h2 (distribution in progressions) | Tao trajectories DO equidistribute on (Z/3^n)* (per Tao's main theorem), but the "level of distribution Q" question is about how *uniformly* the equidistribution holds when restricting to congruence conditions. Kontorovich's setup is about *primes appearing in the orbit*, not about *the orbit's equidistribution rate*. Categorical mismatch. |
| h3 (θ > 0) | Open in the affine-sieve setting. Tao's μ_n equidistribution rate is what R75 controls; this rate isn't a "level of distribution" in Kontorovich's sense. |

## Disposition for I

**NO_FIT (categorical).**

Kontorovich's framework targets *primes or almost-primes in thin orbits* of subgroups of GL_n(ℚ). Tao recursion is *not a group action* (it's a stochastic iteration with v ~ Geom(1/2)). The Syracuse Plancherel mass S_n → 7/15 is not a primes-in-orbit statement.

The "levels of distribution" *language* (averaged equidistribution error in q) is structurally related to R75's Plancherel decomposition (averaged Fourier mass over ξ mod 3^n), but the formal correspondence requires:
- The set S = primes-in-some-thin-set;
- Distribution counted in arithmetic progressions mod q;
- Error control summed over q.

Tao gives instead:
- A measure μ_n on (Z/3^n)*;
- Plancherel decomposition over frequencies ξ;
- Convergence S_n → 7/15.

Different objects.

**Adelic factorization tag:** Same as BGS (candidate H): **GLOBAL_BUT_PLACE_BLIND** — Kontorovich uses archimedean ordering ∥γ∥ ≤ T to count orbit points, but conclusion is about prime/almost-prime structure in ℕ, not place-by-place factorization.

## Mode H circular fingerprint

None: hypotheses fail at h1 categorically.
