# Phase 1 — Structural Map: F̂_p vs μ̂_n

## Object definitions (precise)

### Object A: F̂_p^full(ξ) — verified family-level theorem

From `FHAT_THEOREM_VERIFICATION_RESULTS.md` §8 (verified 33 cells, primes p ∈ {3..31}, r ∈ {1..6}):

> f_p(u) := e_M(c · (1+p)^u),  M := p^{r+1},  u ∈ Z/p^r,  c ∈ (Z/M)^×.
>
> F̂_p^full(ξ) := Σ_{u=0}^{M-1} f̃_p(u) · e_M(-ξu),   f̃_p the period-p^r extension of f_p to Z/M.
>
> **Theorem:** |F̂_p^full(ξ)| = p^{(r+3)/2} on supp = {p·a (mod M) : a ∈ Z/p^r, a ≡ c (mod p)},
> |supp| = p^{r-1}, F̂_p^full ≡ 0 off supp.

Structural ingredients:
- Z/M cyclic group of order p^{r+1}, M = p^{r+1}.
- Single deterministic function f_p (no randomness).
- (1+p) is a generator of principal units 1 + pZ_p mod p^{r+1}, cyclic of order p^r (for p ≥ 3).
- Plancherel on Z/M: Σ_ξ |F̂(ξ)|² = M · Σ_u |f(u)|² = M·M = M².
- Cochrane Theorem 2 / saddle-point / principal-unit Gauss sum equidistribution.

### Object B: μ̂_n(ξ) — Tao Prop 1.17 object

From Tao 2022 (1.21)-(1.26), `R76 §1`, `R77.2 §3.4`:

> Syrac(Z/3^n) := stationary distribution of the Tao Markov chain K_n on (Z/3^n)\*.
>
> Decomposition (Tao 1.26):
>     Syrac(Z/3^n) ≡ 2^{-a_1} + 3·2^{-a_{[1,2]}} + 9·2^{-a_{[1,3]}} + ... + 3^{n-1}·2^{-a_{[1,n]}}  (mod 3^n)
>
> where a_i are iid Geom(2) random variables on N_+ (P(a=k) = 2^{-k} for k ≥ 1), and
>     a_{[1,j]} := a_1 + a_2 + ... + a_j.
>
> μ̂_n(ξ) := E[exp(-2πi · ξ · Syrac(Z/3^n) / 3^n)],   ξ ∈ Z/3^n.
>
> **Tao Prop 1.17 (qualitative):** for ξ with 3 ∤ ξ, |μ̂_n(ξ)| ≪_A n^{-A} for any A > 0.

Structural ingredients:
- Z/3^n cyclic group (specialized to p = 3, single prime).
- Probabilistic expectation over the iid Geom(2) sequence (a_1, ..., a_n).
- Decomposition (1.26) is a Z/3^n random variable built from n nested partial sums of independent geometrics.
- μ̂_n is the characteristic function of a Markov-chain stationary; the chain K_{n+1} maps Syrac(Z/3^n) → 3·Syrac/2^a + 1 (mod 3^{n+1}), so μ̂_{n+1} satisfies Tao's recursion (1.5)/(1.7):
>     μ̂_{n+1}(ξ) = Σ_{v ≥ 1} 2^{-v} · A_v(ξ) · μ̂_n(ξ · 2^{-v} mod 3^n)
- Plancherel on Z/3^n: S_n := Σ_{ξ: 3∤ξ} |μ̂_n(ξ)|² (R75 setup).

## Shared arithmetic substrate (commonalities)

Within technical limits, the two objects share these features:

**(a) Cyclic-group character framework.** Both live on a Z/p^k for p a prime; both are coefficients of characters e_{p^k}(-ξ·u) on those groups. (F̂_p on Z/p^{r+1}, μ̂_n on Z/3^n.)

**(b) Plancherel/Parseval normalisation.** Both satisfy a Plancherel-type identity:
- F̂_p: Σ_ξ |F̂_p|² = M² (exact, deterministic from |f_p(u)| = 1).
- μ̂_n: Σ_ξ |μ̂_n|² = 3^n · ‖π_n‖²_{ℓ²}, a different normalisation since π_n is a probability distribution not a unit-modulus function.

**(c) Principal-unit / multiplicative structure.** F̂_p's support uses (1+p) which generates the principal units. μ̂_n's recursion uses 2^{-v} mod 3^n, which traces out the cyclic subgroup ⟨2⟩ in (Z/3^n)^×. Both involve multiplicative cyclic subgroups of (Z/p^k)^×.

**(d) Conservation law.** R76 §2 (Theorem 76.1) is a Plancherel-style identity: Σ_j M_{n+1}(η_0 + j·3^n) = 0. F̂_p satisfies an analogous coset identity through its periodicity-in-u: F̂_p^short repeats p times around Z/M, summing to F̂_p^full = p·F̂_p^short on the principal coset.

**(e) Saturation at √(group order).** F̂_p saturates the Plancherel bound: |F̂_p|² = p^{r+3} per support point, |supp| = p^{r-1}, gives total Plancherel mass p^{r+3} · p^{r-1} = p^{2r+2} = M². μ̂_n cannot saturate Plancherel in the same way; it has Σ |μ̂|² = O(1) decay (since π_n has ‖π_n‖²_{ℓ²} → 0 like 1/3^n) — but Tao Prop 1.17 expresses a different (qualitative) decay regime.

## Divergences (the gap)

**(I) Deterministic character sum vs probabilistic expectation.**

F̂_p^full(ξ) is a single Z-linear combination of p^{r+1} unit-modulus complex numbers. Magnitude p^{(r+3)/2} comes from Cauchy-Schwarz + equidistribution. There is NO expectation operator anywhere in its definition.

μ̂_n(ξ) is the value of E[exp(-2πi·ξ·X/3^n)] where X = Syrac(Z/3^n) is a probabilistic mixture. Tao Prop 1.17's decay n^{-A} is a probabilistic phase-cancellation result on the (a_1, ..., a_n) law; Tao §7 proves it via white-points / black-region renewal-process arguments on the 2D Pascal random walk, not by Cauchy-Schwarz on a fixed sum.

This is the central structural divergence: the F̂_p saturation uses the Gauss-sum mechanism (deterministic principal-unit characters reach √M magnitude); the μ̂_n decay uses a fundamentally different mechanism (independent Geom(2) summands creating phase cancellation after averaging).

**(II) Group scale mismatch.**

- F̂_p lives on Z/M, M = p^{r+1}. Parameter r is the depth; p is the prime.
- μ̂_n lives on Z/3^n. p = 3 fixed; n is the depth.

If one attempts a parameter specialization p = 3, the F̂_p support is at *level r* with M = 3^{r+1}, indices ξ ∈ {3·a : a ≡ c (mod 3)} ⊂ Z/3^{r+1}. The μ̂_n object at *level n* lives on Z/3^n with index range ξ ∈ Z/3^n. To compare at a single (n, ξ) one would need to align levels (n ↔ r+1?) and pick which (M, ξ) pair to compare.

The natural alignment n = r+1 makes M = 3^n. But F̂_p's support is sparse: only at multiples of 3 inside Z/3^n, i.e. at ξ = 3a for a ∈ Z/3^{n-1}, a ≡ c (mod 3). μ̂_n by Tao Prop 1.17 is bounded at every ξ with 3 ∤ ξ — i.e. precisely the complement of F̂_p's support. **The two objects are concentrated on disjoint frequency sets when level-aligned.**

This is the most concrete structural obstruction. F̂_p has all its mass on ξ ∈ 3·Z/3^n; μ̂_n's Tao-Prop-1.17 decay is at ξ ∉ 3·Z/3^n. They literally don't see the same frequencies.

**(III) Time-scale role of n.**

In F̂_p, r is an algebraic parameter — the depth of the principal-unit tower; the theorem holds for every r ≥ 1 with explicit magnitude p^{(r+3)/2}.

In μ̂_n, n is also a depth parameter but in a different sense: it's the number of iid summands in the (1.26) decomposition. Increasing n adds new randomness; the n^{-A} decay is a long-time iteration/CLT-like rate-of-convergence on the renewal process, not an algebraic scaling.

**(IV) Cochrane vs Tao §7.**

F̂_p's proof template (Move 2 §"Proof template" + `FHAT_THEOREM_VERIFICATION_RESULTS.md` §3):
1. Cochrane Theorem 2 (complete-sum vanishing for the polynomial g(u) = c·(1+p)^u - p²·m·u).
2. Period structure (1+p) order p^r → supp ⊂ p·Z/M.
3. Principal-unit Gauss-sum equidistribution (R78.4-78.6 explicit at q=3).

μ̂_n's proof of Prop 1.17 (Tao §7):
1. Lemma 7.2 white-points (probabilistic phase cancellation on 2-adic words).
2. Lemma 7.4 deterministic triangles (combinatorial structure of black regions).
3. §7.3 Pascal random walk renewal-process analysis.
4. §7.4 α-parameter optimization.

These proof templates do not share machinery beyond "characters on Z/p^k." The arguments are not cousins.

## Cross-reference table

| Aspect | F̂_p^full | μ̂_n |
|---|---|---|
| Group | Z/p^{r+1} (p ≥ 3 free) | Z/3^n (p = 3 fixed) |
| Object type | Deterministic char-sum | E over iid Geom(2)^n |
| Magnitude/bound | = p^{(r+3)/2} exact (theorem) | ≪_A n^{-A} qualitative |
| Support / decay set | supp = {3·a : a ≡ c mod p}, sparse | ξ with 3 ∤ ξ, dense |
| Plancherel sum | M² (saturated) | S_n ∼ 7/15 (R75) |
| Proof tool | Cochrane T2 + Plancherel + Gauss equidist | White-points + black-region + renewal |
| Project source | `FHAT_THEOREM_VERIFICATION_RESULTS.md` | Tao 2022 §7, R75-R78 |

## R76 / R77.2 specific structural objects (potential bridge anchors)

Within the μ̂_n world, R76 distilled the rate problem to:
- **Bilinear pair moment:** M_n(η) = Σ_{ξ:3∤ξ} μ̂_n(ξ) · μ̂_n*(ξ·η).
- **Leading-mode identity (Th. 76.3):** S_{n+1} = -2·M_{n+1}(1+3^n).
- **Class-resolved 2×2 deviation operator on (P_+, P_-)** with rigorously-derived diagonal **T_diag = (1/5)·[[1,1],[4,4]]** (R77.1 Theorem 77.1) with spec {0, 1} and (1,4) the asymptotic-mass eigenvector.
- **Order-3 companion T_3** (R77.2 §2.3) with spec {1/2, 1/4, 1/8} — but spec(T_3) = {1/2, 1/4, 1/8} is **conditional** on the 3-mode fit ε_n = A(1/2)^n + B(1/4)^n + C(1/8)^n being structurally exact, not just an n ≤ 6 empirical regression.
- **R77.2 §3.4: ε_K from Tao Prop 1.17.** Nisoli needs ε_K = ‖T - T_K‖_op ≤ C_A²·N^{-2A}. Tao 1.17 is the input for this ε_K. **C_A is the qualitative constant.**

Within the F̂_p world, the parallel objects are:
- F̂_p^full(ξ) for ξ ∈ supp.
- K_p (short-window) = (1/M)·Σ 1̂_N(ξ)·F̂_p(ξ) — the q-sweep object.
- The Pólya-Vinogradov decomposition K_p ↔ F̂_p needs a bilinear / Burgess bound.

The **R76/R77.2 family** lives entirely on the μ̂_n side. The **F̂_p family** lives on the deterministic-character-sum side. **The two families do not share a single computational object** in the project — they share *language* (Plancherel, characters on Z/3^n) but not *machinery*.

## Tao §7's role in any bridge

Per Move 2 §Phase 4 explicitly:
> "The R78.3 magnitude formula sits at the deterministic-character-sum level. Even at q=3 it cannot substitute for Tao Prop 1.17's `|μ̂_n(ξ)| ≪_A n^{-A}` without an intermediate step that translates the deterministic-cyclic-group bound into a Markov-chain-stationary bound. That intermediate step is not in the project documents."

Any candidate Φ in Phase 2 that bridges F̂_p → μ̂_n must produce or replace that intermediate step. The intermediate step, in any honest derivation, is exactly Tao §7's white-points/black-region/renewal machinery — because that's the machinery that turns the (Geom(2))^n iid summand structure (1.26) into a character-sum decay. **A2 (anti-tautology):** any candidate Φ whose derivation imports Tao §7.2-7.4 reproduces the bookkeeping result; it is not an independent bridge.

## Summary of Phase 1

The two objects share **cyclic-group character framework + Plancherel notation + multiplicative subgroup language**. They diverge on **deterministic vs probabilistic, scale-of-group parameter, support-set, and proof-machinery**. Crucially, at the natural level alignment n = r+1, the F̂_p support and the μ̂_n decay-relevant frequencies are **disjoint** (Divergence II above).

This Phase 1 finding constrains Phase 2: any candidate bridge must explain how an object supported on multiples-of-p in Z/p^{r+1} (F̂_p) influences an object that decays on the complement (μ̂_n at coprime-to-p frequencies). This is a strong structural constraint and most natural candidates will fail it on inspection.

**Phase 1 does NOT terminate early.** The constraint above narrows candidates but doesn't trivially rule them all out — one can imagine constructions in category (B) "averaging" or (C) "coupling" that legitimately span the support divide. Proceed to Phase 2.
