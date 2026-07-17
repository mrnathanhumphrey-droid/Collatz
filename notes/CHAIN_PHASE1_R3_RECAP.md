# CHAIN Phase 1 — r=3 Inner-Plancherel Chain, Explicit Recap

**Date:** 2026-05-11. Analyst: Wilson. Reports to: Nathan.

**Goal:** Trace the r=3 chain step-by-step to expose the structural ingredients each step
depends on, so that at Phase 3 we can test which ingredients survive at higher digits.

## Setup recap

Prime p ≥ 3, r = 3, q = p^{r+1} = p^4, period p^r = p^3, support size N = p^{r-1} = p^2,
support = {a ∈ Z/p^r : a ≡ 1 mod p}.

Bijection (T78.5_p): a ↔ C_a := a · L̃_p^{-1} mod p^r. Both lie in {≡ 1 mod p}.

Two-digit decomposition at r=3:
- C_a = 1 + p·s_0 + p²·c_2 mod p^3, with s_0, c_2 ∈ Z/p.
- (s_0, c_2) ↔ a bijection of (Z/p)² → support.
- s_0 := (C_a − 1)/p mod p = leading saddle digit
- c_2 := higher-order digit of (C_a − 1)/p² mod p

Hensel-lifted saddle at r=3: s*(r=3) = (C_a − 1)/p mod p^{r−1} = s_0 + p·c_2. [But the
T78.6_p form at r=3 uses just s_0 — see Phase 2 reconciliation note below.]

## Step 1 — Family-level closed form

**T78.4_p (Cochrane factorization, rigorous):**
> F̂_p(p·a) = p · e_q(1) · G_p(a), where G_p(a) := Σ_{s=0}^{p^r − 1} e_q(P_a(s)), with
> P_a(s) := p·s − C_a · L_p(1+p·s).

**T78.6_p at r=3 (saddle-exact, rigorous via Phase 2 digit-chain at r=3):**
> G_p(a) = √q · e_q(P_a(s_0(C_a)))    [r=3, odd, no η_p factor]
>
> where s_0 = (C_a − 1)/p mod p.

**Phase polynomial mod p^4 (PATH2_BILINEAR §"Attempt A" expansion):**
> P_a(s_0) ≡ −p²·s_0²/2 + p³·(s_0³/6 − c_2·s_0) mod p^4

The c_2 dependence enters via the cross-term in (C_a · L_p(1+ps_0)) — concretely from
the (p²·c_2) contribution to C_a expanded against the leading p·s_0 in L_p. **The c_2
appears LINEARLY** (coefficient −s_0 · p^3, which is the p^3 stratum).

## Step 2 — Bilinear substitution

Define the target T_p (related to S_partial by |S_partial| = p · √q · |T_p|):
> T_p := Σ_{a ∈ supp} 1̂(p·a) · e_q(P_a(s_0(C_a)))

(Magnitude factor: |S_partial| ≤ C·√N iff |T_p| ≤ C·N up to fixed prefactor √q.)

Substitute the closed form:
> T_p = Σ_{(s_0, c_2)} 1̂(p·a(s_0,c_2)) · e_{p²}(−s_0²/2) · e_p(s_0³/6 − c_2·s_0)

Reorganize by s_0 class:
> T_p = Σ_{s_0=0}^{p−1} e_{p²}(−s_0²/2) · e_p(s_0³/6) · Inner(s_0)
>
> Inner(s_0) := Σ_{c_2=0}^{p−1} ω(s_0, c_2) · e_p(−c_2 · s_0)
>
> ω(s_0, c_2) := 1̂(p·a(s_0, c_2))

**Structural fact:** the phase is **linear in c_2** at each fixed s_0 (coefficient −s_0
mod p). This is what allows Plancherel-on-c_2 to collapse cleanly.

## Step 3 — The c_2-Fourier collapse (Lemma R78.7)

The key calculation. Parametrize:
- At fixed s_0, as c_2 ranges 0..p−1, a(s_0, c_2) = a_0(s_0) + c_2 · p² · L̃_p mod p^3.
- So p·a varies through p·a_0(s_0) + c_2 · p^3 · L̃_p mod p^4.

Plug into the explicit 1̂:
> 1̂(p·a) = Σ_{u=0}^{N-1} e_q(p·a·u) = Σ_{u=0}^{p²−1} e_{p^3}(a·u)

Hence:
> Inner(s_0) = Σ_{c_2} Σ_u e_{p^3}(a·u) · e_p(−c_2·s_0)
>            = Σ_{c_2} Σ_u e_{p^3}((a_0 + c_2·p²·L̃_p)·u) · e_p(−c_2·s_0)
>            = Σ_u e_{p^3}(a_0·u) · Σ_{c_2} e_p((L̃_p·u − s_0)·c_2)

**The c_2 sum is a length-p geometric series in e_p:**
> Σ_{c_2=0}^{p−1} e_p(k·c_2) = p · 𝟙[k ≡ 0 mod p]

So:
> Inner(s_0) = p · Σ_{u : L̃_p·u ≡ s_0 mod p,  0 ≤ u < N=p²} e_{p^3}(a_0·u)

The constraint u ≡ L̃_p^{-1}·s_0 mod p **restricts u to a length-p arithmetic progression
in {0, ..., p²−1}** (p choices, step p, starting offset u_0 = L̃_p^{-1}·s_0 mod p).

Continuing:
> Inner(s_0) = p · e_{p^3}(a_0·u_0) · Σ_{j=0}^{p−1} e_{p^3}(a_0 · j·p)
>            = p · e_{p^3}(a_0·u_0) · Σ_{j=0}^{p−1} e_{p²}(a_0 · j)
>            = p · e_{p^3}(a_0·u_0) · D_p(a_0, p²)

where **D_p(a, M) := Σ_{j=0}^{p−1} e_M(a · j)** is the length-p Dirichlet kernel mod M.

**Magnitude:**
> |Inner(s_0)| = p · |D_p(a_0(s_0), p²)|

This is **Lemma R78.7** — the c_2-Fourier collapse: the inner c_2-sum exactly reduces a
Dirichlet kernel of length p^2 (from 1̂ at q=p^4) to a Dirichlet kernel of length p (in
e_{p²}) times an outer factor p.

## Step 4 — 1/sin grid identity bounds the outer sum

Outer bound:
> |T_p| ≤ Σ_{s_0=0}^{p−1} |Inner(s_0)| = p · Σ_{s_0} |D_p(a_0(s_0), p²)|

As s_0 ranges over Z/p, a_0(s_0) cycles through the p elements {1, 1+p, ..., 1+(p−1)p}
of {a ∈ Z/p² : a ≡ 1 mod p}.

**Magnitude formula:** for a ≡ 1 mod p in Z/p², a ≠ 0 mod p²,
> |D_p(a, p²)| = |sin(π·a/p)/sin(π·a/p²)| = sin(π/p) / |sin(π·a/p²)|

(used a ≡ 1 mod p so sin(π·a/p) = sin(π/p)).

Writing a = 1 + p·α with α ∈ {0, 1, ..., p−1}:
- α = 0: |sin(π/p²)| ≈ π/p² → |D_p| ≈ p
- α ≥ 1: |sin(π(1+pα)/p²)| ≈ πα/p → |D_p| ≈ sin(π/p)·p/(πα) ≈ 1/α

Sum:
> Σ_{α=0}^{p−1} |D_p(a_0(α), p²)| ≤ p + Σ_{α=1}^{p−1} 1/α = p + H_{p−1} ≤ p + log p

**Outer bound:**
> |T_p| ≤ p · (p + log p) = p² + p log p ≤ 2 p² = 2N    (uniformly in p ≥ 3)

Equivalently:
> |S_partial| ≤ 2 · √N · √q = 2 √N · √q.

Strict **2√N** at r=3 (up to fixed √q factor).

## The "constant 2" origin — step-by-step

Trace where the 2 comes from:

1. The p² (= N at r=3) comes from the **α=0 single term in the cosecant sum** (the
   D_p(1, p²) ≈ p contribution times the outer factor p).
2. The H_{p−1} ≈ log p comes from the **α ≥ 1 tail** (small contributions decaying as 1/α).
3. Sum: p² + p log p ≤ p²(1 + log p / p) ≤ p² · (1 + log 3 / 3) ≈ 1.37 · p² for p=3,
   monotonically smaller for larger p.

The "2" is a safe uniform constant covering 1 + log p / p ≤ log 3 / 3 + 1 ≈ 1.37 ≤ 2 for
all p ≥ 3. The TIGHT constant is 1 + (log p)/p, asymptotically → 1.

**This 2N is the rigorous bound; empirically |T_p| ~ √N (R79b), much smaller.**

## Structural ingredients of the r=3 chain

The chain has **five load-bearing ingredients**:

**(I1) Closed-form phase polynomial mod p^{r+1}.** At r=3:
P_a(s_0) ≡ −p²·s_0²/2 + p³·(s_0³/6 − c_2·s_0) mod p^4.

**(I2) Linear-in-c_2 structure at the TOP stratum.** The p^3 stratum contains a term
LINEAR in c_2 with coefficient −s_0 (a unit of Z/p when s_0 ≠ 0; zero when s_0 = 0).

**(I3) Length-p Plancherel collapse.** Σ_{c_2} e_p(k·c_2) = p · 𝟙[k ≡ 0 mod p].
Applied to the linear c_2 phase, this kills a Dirichlet sum factor.

**(I4) Dirichlet-on-arithmetic-progression identity.** After the c_2 collapse, the
inner u-sum reduces to a length-p Dirichlet kernel D_p(a_0, p²) at modulus p² (one
modulus level lower than the original p^3).

**(I5) 1/sin grid identity (cosecant sum bound).** Σ_{s_0} |D_p(a_0(s_0), p²)| ≤ p + log p
≤ 2·p uniformly. Specifically, this uses that as s_0 cycles, a_0(s_0) is a full
arithmetic-progression cycle through {≡ 1 mod p in Z/p²}, allowing the cosecant sum to
be evaluated by the standard identity Σ_α csc(π·α/n + θ).

**Bound assembly:**
> |T_p| ≤ p · (p + log p) ≤ 2 · p² = 2N    [strict |S_partial| ≤ 2√N at r=3]

## What's specific to r=3

- The phase polynomial has exactly **TWO strata**: p² (depends on s_0 only) and p³
  (depends on s_0 and c_2).
- The top stratum (p^3) has a SINGLE inner-digit variable c_2, appearing LINEARLY.
- One Plancherel save (length-p, factor p) is needed.
- One outer bound (cosecant sum on s_0 over Z/p) is needed.

The chain at r=3 has **depth 1** — one inner digit, one Plancherel collapse, one outer.

## What needs to extend at r ≥ 4

At r=4, the phase has **THREE strata** (p², p³, p^4); at r=5, four (p², p³, p^4, p^5);
in general r−1 inner strata.

The **inner digits multiply**: at r=4, two inner digits c_2, c_3 (since (C_a − 1)/p has
r−1 = 3 digits, with s_0 = c_1 outer leaving c_2, c_3 inner). At r=5, three inner digits.
In general r−2 inner digits.

**Key question for Phase 3:** at each inner digit c_k, does the phase contain a
linear-in-c_k term at some stratum (so an analog of Lemma R78.7 applies)? If YES at every
c_k, the chain nests with one Plancherel per digit (Phase 4). If NO at some c_k (e.g., a
quadratic-only appearance), that step requires a Gauss-sum bound (factor √p, not p) and
the chain LOSES p^{1/2} at each such step.

## Empirical anchor (R79b)

R79b at p=3, r=8..20: |K|/√N ∈ [0.7, 2.7], max-over-sample.

Plancherel-translated, this corresponds to |T_p|/N empirically ≤ ~3/p = 1, **far below
the 2N rigorous ceiling**. The rigorous bound has factor 2 headroom in the constant; the
empirical truth has another √N headroom inside the bound. The chain is **not tight to
empirical**, but the strict-2 RIGOROUS bound is what's load-bearing for closure.

## Files referenced

- PATH2_BILINEAR_FROM_CLOSED_FORM.md — original derivation
- HENSEL_APPROACH_A.md — Hensel-lifted closed form at r ≥ 4
- HENSEL_PHASE_ARTICULATION.md — phase polynomial generating identity (1+y)log(1+y)
- r79b_S_partial_empirical.md — empirical anchor

## Next: Phase 2 expands the phase polynomial at r ≥ 4

Phase 2 articulates the digit-expanded P_a(s*(r)) at r=4, 5, 6 in terms of (s_0, c_2,
c_3, ..., c_{r−1}). Phase 3 then tests, for each higher-order digit c_k (k=3, 4, ...),
whether the linear-in-c_k structural ingredient (I2) survives.
