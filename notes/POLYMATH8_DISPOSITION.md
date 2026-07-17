# POLYMATH8 DISPOSITION

## Disposition: **POLYMATH8_NO_SUBTHREAD**

The Polymath8 chain confirms the broader obstruction pattern that has stacked across the five prior lit-scan candidates. No variant in the chain treats rough Fourier-completed amplitudes, prime-power moduli, or cardinality-side √-savings. The framework's three structural commitments — smooth amplitudes on Fourier-completed variables, square-free moduli, Weil-bound √-savings — appear unbypassable within the chain.

---

## One-paragraph rationale

Across all 28 blog posts, 3 arXiv papers (Polymath8a 1402.0811, Polymath8b 1407.4897, retrospective 1409.8361), the Type I/II/III framework holds three load-bearing conditions invariant: **(i) the modulus class `S_I` is defined as square-free numbers** (Polymath8a paper abstract verbatim: "large smooth squarefree moduli"; Type I/II post Definition 2 verbatim: "S_I denotes the square-free numbers whose prime factors lie in I"); **(ii) every Fourier-completed amplitude must be smooth** (Definition 1(iv); Type III amplitudes ψ₁, ψ₂, ψ₃ explicit smoothness requirement; the "rough α" allowance in Type I/II applies only to the outer Cauchy-Schwarzed variable, not the Fourier-completed side); **(iii) all power savings come from Weil/Deligne** (Type I/II post: "the Weil conjectures are the primary source of power savings"; Type III post: Bombieri-Birch for surfaces from Deligne's Weil II). Polymath8b never re-opens the equidistribution side — it consumes Polymath8a's estimates as a black box and innovates only on the multidimensional Selberg sieve (Maynard). Dense divisibility (singly, k-tuply) does *generalize* Zhang's y-smoothness, but stays within square-free and stays inside the Weil-bound framework — it is a relaxation of *how* you factor `q`, not of *what kind* of `q` is allowed. The chain ends with the September 2014 retrospective; no "Polymath8c" exists. Three independent structural obstructions block extension to `p^r` moduli (Weil-bound ramification at the special fibre, breakdown of CRT decomposition, degenerate dense-divisibility for prime-powers), and Ramanujan-sum cancellation — the one non-Weil saving source mentioned anywhere in the chain — itself depends on square-free moduli.

---

## Essential conditions confirmed unbypassable across the chain

| Condition | Polymath8a | Polymath8b | Subsequent retrospective | Verdict |
|-----------|-----------|------------|-------------------------|---------|
| Square-free modulus class `S_I` | Required (Def 2) | Inherited as black box | Restated | UNBYPASSED in chain |
| Smooth amplitudes on Fourier-completed side | Required (Def 1(iv); Type III ψᵢ) | Inherited | Restated | UNBYPASSED in chain |
| Siegel-Walfisz for the "non-completed" sequence | Required for β | Inherited | Restated | Cannot substitute for Fourier-completion smoothness on a prime-power-modulus target |
| Weil / Deligne √-savings | Source of all power savings | Inherited | Restated | UNBYPASSED in chain |
| Cauchy-Schwarz halving every bilinear step | Required architecturally | Inherited | Restated | UNBYPASSED in chain |
| Dense divisibility on q | Required (replacing y-smoothness) | Inherited | Restated | Internal generalization only; still inside `S_I` |
| CRT-compatible singleton congruence system | Required (Def 2) | Inherited | Restated | Requires square-free; cannot accept `p^r` |

---

## What this means for the Collatz target

The target bound `|Σ_{a ≡ 1 mod p in ℤ/p^r} 1̂(p·a) · F̂(p·a)| ≤ C · N · √q` fails the chain's structural requirements on all three counts simultaneously:

1. **Modulus.** `q = p^{r+1}` for fixed prime p, varying r ≥ 7 — exactly the non-square-free regime the chain excludes by definition.
2. **Amplitude.** `1̂(ξ) = Σ_{u=0}^{N-1} e_q(ξ u)` is a *Dirichlet kernel* — the canonical non-smooth Fourier-completed object. It cannot inhabit the "rough outer α" position because it IS the Fourier completion.
3. **Saving.** The empirical `√q` is saturating at `C ≈ 2.0` from a coset-cardinality count (`|principal-unit coset| · √q` factored differently), not from Weil-bound oscillatory cancellation, and not from an amplifier on a non-abelian group as in Pascadi 2025.

This puts Polymath8 in the same disposition class as the other five candidates already ruled out: structural mismatch, not parameter mismatch. Adapting Polymath8 to this problem would require **building a new variant from scratch** with: (i) a prime-power-modulus analog of the trace-function / Weil-bound machinery, (ii) a Fourier-completion tool that works for Dirichlet-kernel amplitudes, (iii) a cardinality-side saving accounting that the chain never developed. None of these is an extension that exists or has been outlined in the chain.

---

## Honesty caveat (per A4 in the brief)

A theoretical extension of dense divisibility's iterated-factorization idea to `p^r` is *conceivable* (the factors `p^i · p^{r-i}` exist), but the chain explicitly never develops it, and there is good reason it doesn't: every chain factorization needs to *land in a different range* for van der Corput to give gain, which is precisely what prime-power factorizations can't do. This is "subthread not written and would require new work," not "subthread found." It is filed under POLYMATH8_NO_SUBTHREAD, not POLYMATH8_INCONCLUSIVE.

---

## Citation anchor list (for downstream use)

- Tao, "Estimation of the Type I and Type II sums," 2013-06-12, https://terrytao.wordpress.com/2013/06/12/estimation-of-the-type-i-and-type-ii-sums/ — Definitions 1, 2; Theorem 3.
- Tao, "Estimation of the Type III sums," 2013-06-14, https://terrytao.wordpress.com/2013/06/14/estimation-of-the-type-iii-sums/ — Theorem 2; Bombieri-Birch.
- Tao, "An improved Type I estimate," 2013-07-27, https://terrytao.wordpress.com/2013/07/27/an-improved-type-i-estimate/ — Theorem 4; Definition 1 (k-tuply dense divisibility); Proposition 11.
- Tao, "Distribution of primes in densely divisible moduli," 2013-06-23 — dense divisibility introduced; modulus class still `S_I`.
- Tao, "Distribution of primes in doubly densely divisible moduli," 2013-07-07 — `MPZ″` improvement; still `S_I`.
- Polymath (D.H.J.), "New equidistribution estimates of Zhang type," arXiv:1402.0811, Algebra & Number Theory 8-9 (2014) 2067–2199 — definitive Type I/II/III statements. **Abstract verbatim: "large smooth squarefree moduli."**
- Polymath (D.H.J.), "Variants of the Selberg sieve, and bounded intervals containing many primes," arXiv:1407.4897 — sieve side only; uses 1402.0811 as black box.
- Polymath (D.H.J.), "The 'bounded gaps between primes' Polymath project — a retrospective," arXiv:1409.8361 — confirms no further equidistribution work in Polymath8 proper.
