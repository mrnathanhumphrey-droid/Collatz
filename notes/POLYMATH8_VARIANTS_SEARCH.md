# POLYMATH8 VARIANTS SEARCH — (a) rough inputs, (b) prime-power moduli, (c) cardinality savings

**Purpose.** Search the Polymath8 chain for any variant of Type I/II/III that addresses (a) rough amplitudes (e.g., Dirichlet kernels, characteristic functions of APs without smooth cutoff), (b) prime-power moduli `p^r` (`r ≥ 2`), or (c) abelian cardinality-side savings (saving a `√q` from a coset cardinality rather than a Weil bound or amplifier).

**Scope.** Only the Polymath8 chain proper as mapped in POLYMATH8_CHAIN.md. Descendants and parallel literature explicitly out of scope.

**Adversarial discipline (per brief A1–A4).** Theorem-level citation required. Passing mentions do not count. Distinguish "chain proper" from "descended-but-different-method." Don't speculate.

---

## (a) Rough amplitudes — is there a Polymath8 variant that avoids smooth-amplitude requirements?

### Best candidate within the chain: the Siegel-Walfisz alternative to smoothness for the β-side of Type I/II

**Status.** Polymath8a Definition 1(iii)/(iv) is explicit: a coefficient sequence may satisfy EITHER the Siegel-Walfisz property OR be smooth. Theorem 3 (Type I/II) and Theorem 4 (Improved Type I) demand only Siegel-Walfisz on β, not smoothness on β.

**Why this isn't actually a hit for Collatz.**

- The Siegel-Walfisz property `|Δ(α 1_{(·,q)=1}; a(r))| ≪ τ(qr)^O(1) N log^{-A} x` is about equidistribution of α's restriction to APs modulo `r` — it is a property of arithmetic progressions, NOT a Fourier-decay or Schwartz-control statement about α(n) as a function of n.
- The α (the *un-completed*, "rough" side of the bilinear sum) in Type I/II requires nothing about its Fourier behavior — it gets squared out by Cauchy-Schwarz. So in a *narrow* sense, "rough α" IS allowed.
- BUT: the bilinear-Cauchy-Schwarz step then converts the rough α into a divisor-function-bounded *quadratic* sum over `α(n)α(n')`, which is then handled by van der Corput on q-factors with completion of *smooth* sums on the *q-residue* side via Poisson. The smoothness has migrated from the amplitude to the q-residue completion — it has not been eliminated.
- Therefore the chain does NOT actually treat rough amplitudes in the Fourier-on-(ℤ/qℤ) sense the user needs. The "rough α" allowance is for the OUTER bilinear sequence, which is precisely the part that gets cancelled by Cauchy-Schwarz before any Fourier analysis touches it. The Fourier-decaying side IS always smooth.

### Type III post (2013-06-14) explicit statement

> "Only two of three [of ψ₁, ψ₂, ψ₃] need to have any sort of smoothness (because completion of sums is applied in these variables); the third sequence is eventually eliminated by Cauchy-Schwarz." (paraphrased from the chain's own discussion)

**This is the closest the chain comes to a relaxation, and it confirms the pattern, not breaks it.** Smoothness is required precisely on the variables that get *Fourier-completed*. The Collatz target IS a Fourier-completed sum — it cannot survive in the role of an "outer" Cauchy-Schwarzed variable.

### Conclusion for (a)

**No variant in the Polymath8 chain handles a rough amplitude that needs to be Fourier-completed against a non-smooth modular structure.** The chain's "rough α" allowance is an artifact of Cauchy-Schwarz position, not a genuine smoothness relaxation. The Dirichlet kernel `1̂(ξ) = Σ_{u=0}^{N-1} e_q(ξ u)` cannot enter as outer-α because the user's bound IS a Fourier-completed bound (the Dirichlet kernel IS the completion); there is no further Cauchy-Schwarz to absorb it into.

---

## (b) Prime-power moduli `p^r` (`r ≥ 2`) — is there any Polymath8 variant?

### Search procedure

1. Definition of `S_I` in Type I/II post: **"the square-free numbers whose prime factors lie in I."** Verbatim. Square-free is part of the *definition* of the modulus class.
2. Improved Type I post: same modulus class, with added k-tuple dense divisibility constraint.
3. DDM post (2013-06-23): replaces y-smoothness with dense-divisibility; **still inside `S_I` (square-free).**
4. DDDM post (2013-07-07): k-tuply densely-divisible; **still square-free.**
5. Polymath8a paper (arXiv:1402.0811) abstract: **"distribution estimates for primes in arithmetic progressions to large smooth squarefree moduli."** Verbatim.
6. Polymath8b paper (arXiv:1407.4897): uses 8a estimates as a black box; no new modulus class.

### Why the chain cannot extend to prime-power moduli

Three structural reasons identified across the chain:

1. **Weil bound applicability.** Deligne's Weil II gives `√q`-savings for exponential sums modulo `q` only when the underlying algebraic-geometric setup is *smooth and proper*. Going from `q = p` to `q = p^r` (`r ≥ 2`) introduces ramification at the special fibre — the trace-function machinery used in Polymath8a (and the more refined Fouvry-Kowalski-Michel inputs it cites) is set up for `(ℤ/pℤ)` not `(ℤ/p^rℤ)` for high r. The chain never proves a `p^r`-modulus analog.

2. **Chinese Remainder factorization.** Definition 2's singleton-congruence-class system is multiplicative across coprime square-free pieces. Prime-power moduli `p^r` cannot be CRT'd into smaller pieces of the same modulus — the entire factorization-and-recombination strategy collapses.

3. **Dense divisibility is incompatible with prime-power.** A modulus `p^r` is `y`-densely-divisible for `y ≥ p`, but it has only `r + 1` divisors (`1, p, p², ..., p^r`). The factorization-at-every-scale property that Polymath8a's iterated van der Corput needs is degenerate — every "factorization" is a power of `p`. The whole point of dense divisibility was to allow factoring `q = q_1 q_2 q_3 ...` with q_i in non-trivially different ranges; for `p^r` this gives nothing the trivial `p^r = p^{r-i} · p^i` doesn't.

### Conclusion for (b)

**No.** The Polymath8 chain is structurally constrained to square-free moduli by THREE independent mechanisms (Weil bound geometry, CRT decomposition, dense divisibility factoring). None of these has a documented prime-power analog within the chain.

---

## (c) Abelian cardinality-side savings — is there any saving extracted from coset cardinality rather than amplifier / Weil?

### Search

1. Type I/II post explicitly: "**The Weil conjectures are the primary source of power savings (x^{-c} for some fixed c > 0) in the argument**, but they need to overcome power losses coming from completion of sums."
2. Type III post: explicit Bombieri-Birch (Deligne for surfaces).
3. Improved Type I post: Proposition 11 is van der Corput; the input bound for one-dim exponential sums is again Weil/Deligne.
4. DDM, DDDM posts: same Weil-bound machinery, no parallel construction.
5. Polymath8a paper abstract: "exponent of distribution 1/2 + 7/300" — the `7/300` is what one gets *after* Weil-bound √-savings interact with q-factor combinatorics. No non-Weil saving source is described.

### Ramanujan-sum aside

The Type III post notes: "Ramanujan sums also play a critical role, exhibiting better-than-square-root cancellation." This is a structural cancellation in `c_q(n) = Σ_{d | (q,n)} d μ(q/d)`, but:
- It is NOT a cardinality-side saving on a coset.
- It is a CONSEQUENCE of multiplicativity over square-free moduli (`c_q` is multiplicative when `q` is square-free).
- For `q = p^r` with `r ≥ 2`, `c_{p^r}(n)` has a degenerate form that doesn't give the same cancellation regime.

So Ramanujan sums are not a cardinality-side analog of Weil; they are a specific tool that itself depends on the square-free assumption.

### Conclusion for (c)

**No.** Every power saving in the Polymath8 chain originates from Weil / Deligne or from Ramanujan-sum cancellation, both of which are structurally tied to square-free moduli. There is no parallel cardinality-side construction in the chain.

---

## Final cross-reference table

| Variant the user needs | Closest item in chain | Verdict |
|------------------------|------------------------|---------|
| (a) Rough amplitude on the Fourier-completed side | Type I/II's "rough α" allowance | NOT A HIT — rough α is outer-Cauchy-Schwarz position only, not Fourier-completed |
| (a') Siegel-Walfisz substitute for smoothness | Definition 1(iii); used for β of Type I/II | NOT A HIT — SW is an AP-equidistribution property, not what a fixed-prime-power Fourier sum needs |
| (b) Prime-power moduli `p^r` | nothing in chain | NOT A HIT — three structural obstructions |
| (c) Cardinality-side √-savings | Ramanujan-sum cancellation in Type III | NOT A HIT — itself requires square-free moduli |
