# CHAIN Phase 3 — Per-Digit Fourier-Collapse Test at r ≥ 4

**Date:** 2026-05-11. Analyst: Wilson. Reports to: Nathan.

**Goal:** For each higher digit c_k (k = 2, 3, ..., r−1), test whether an analogous
Fourier-collapse lemma applies. The chain nesting order is fixed: peel from the TOP
(largest k) DOWN.

## Setup

From Phase 1's bilinear substitution + ω parametrization:
> T_p = Σ_{c_1=0}^{p−1} e_{p²}(−c_1²/2) · e_p(...) · Inner(c_1)

where Inner(c_1) := Σ_{(c_2, ..., c_{r−1}) ∈ (Z/p)^{r−2}} ω(c_1, c_2, ..., c_{r−1}) ·
e_q(Q_3·p³ + Q_4·p^4 + ... + Q_r·p^r).

In e_q form (with q = p^{r+1}): each stratum m contributes e_{p^{r+1−m}}(Q_m). The
relevant strata for the inner digits are m = 3, 4, ..., r (= top stratum).

The plan: peel digits c_{r−1}, c_{r−2}, ..., c_2 in turn, evaluating one Plancherel per
digit at the stratum where it first appears LINEARLY.

## Test for c_{r−1} (the top digit) — COLLAPSE WORKS

**Phase 2 fact:** c_{r−1} first appears at stratum m=r, LINEARLY, with coefficient −c_1.

At stratum m=r, the phase factor is e_{p^{r+1−r}}(Q_r) = e_p(Q_r). With Q_r =
−c_1·c_{r−1} + (lower-digit-only stuff), we have:

> Sum over c_{r−1}:
>   Σ_{c_{r−1}=0}^{p−1} ω(c_1, c_2, ..., c_{r−2}, c_{r−1}) · e_p(−c_1·c_{r−1})

The ω(c_1, ..., c_{r−1}) is 1̂(p·a(c_1, c_2, ..., c_{r−1})). At fixed (c_1, ..., c_{r−2}),
as c_{r−1} varies, a shifts by c_{r−1} · p^{r−1} · L̃_p mod p^r, so p·a shifts by
c_{r−1} · p^r · L̃_p mod p^{r+1}.

Plug into 1̂ definition (length N=p^{r-1}):
> 1̂(p·a) = Σ_u e_q(p·a·u) = Σ_u e_{p^r}(a·u)

so
> 1̂(p·(a_base + c_{r−1}·p^{r−1}·L̃_p)) = Σ_u e_q(p·a_base·u + c_{r−1}·p^r·L̃_p·u)
>                                       = Σ_u e_q(p·a_base·u) · e_p(L̃_p·c_{r−1}·u)

Combining with the e_p(−c_1·c_{r−1}) phase factor (the linear-in-c_{r−1} from Q_r):

> Σ_{c_{r−1}} 1̂(p·a) · e_p(−c_1·c_{r−1})
>   = Σ_u e_q(p·a_base·u) · Σ_{c_{r−1}} e_p((L̃_p·u − c_1)·c_{r−1})
>   = p · Σ_{u : L̃_p·u ≡ c_1 mod p} e_q(p·a_base·u)
>   = p · Σ_{u : u ≡ L̃_p^{-1}·c_1 mod p, 0 ≤ u < p^{r-1}} e_q(p·a_base·u)

**This is the analog of Lemma R78.7 at the TOP digit.** It saves a factor p (from the
Plancherel) and restricts u to a length-p^{r−2} arithmetic progression (instead of full
length p^{r−1}).

**c_{r−1} COLLAPSE: WORKS, length-p, factor p save.**

Define the post-c_{r−1}-collapse object:
> T_p^{(1)}(c_1, c_2, ..., c_{r−2}) := Σ_{u ≡ u_0(c_1) mod p, 0 ≤ u < p^{r-1}} e_q(p·a_base·u)
>                                       · (residual c_2, ..., c_{r−2} phase from strata 3, ..., r)

where the residual phase at this point includes:
- Stratum r: e_p(Q_r − (c_{r−1}-linear part)) = e_p(Q_r |_{c_{r−1}=0})
- Strata 3, ..., r−1: e_{p^{r+1−m}}(Q_m), unchanged (didn't involve c_{r−1})

## Test for c_{r−2} (second-top digit)

**Phase 2 fact:** c_{r−2} appears at strata m = r−1 and m = r.
- At m = r−1: LINEAR in c_{r−2} (coefficient −c_1).
- At m = r: LINEAR in c_{r−2} (coefficient −c_2, from j=2 (2, r−2) composition; plus
  contribution from j=3 (1, 1, r−2) giving c_1²-coefficient, etc.).

The c_{r−2} phase at stratum r−1 is e_{p²}(Q_{r−1}). Stratum m = r−1 contributes
e_{p^{r+1−(r−1)}}(Q_{r−1}) = e_{p²}(Q_{r−1}).

Q_{r−1} = (j=2) −c_1·c_{r−2} − ... (other lower-digit terms not involving c_{r−1} since
c_{r−1} first appears at m=r).

Specifically (re-derived from compositions):
> Q_{r−1} = −c_1·c_{r−2} − (1/2)·Σ_{i+j=r−1, 1 ≤ i,j ≤ r−2, i ≠ j} c_i·c_j − (1/4)·Σ_{2i=r−1} c_i²·(missing terms from j≥3)

Wait — let me just use the r=4 stratum m=3 case as a concrete example:
> Q_3 (r=4) = −c_1·c_2 + c_1³/6.

So at r=4, the m=3 stratum has Q_3 = −c_1·c_2 + c_1³/6. The c_2 term is linear with
coefficient −c_1. The c_1³/6 is c_2-free.

The phase contribution at m=3 stratum: e_{p^{r+1−3}}(Q_3) = e_{p^{r−2}}(Q_3) =
e_{p²}(Q_3) at r=4.

> e_{p²}(−c_1·c_2 + c_1³/6)

This is a phase modulo p^2 in c_1 and c_2. The c_2 dependence: e_{p²}(−c_1·c_2).

**Note:** the modulus is **p², NOT p**. This is because c_{r−2} appears at the SECOND-TO-TOP
stratum (m=r−1), where the phase modulus is e_{p^{r+1−(r−1)}} = e_{p^2}.

Sum over c_{r−2}:
> Σ_{c_{r−2}=0}^{p−1} (post-c_{r−1}-collapse stuff) · e_{p²}(−c_1·c_{r−2})

The "post-c_{r−1}-collapse stuff" is T_p^{(1)}, which depends on c_{r−2} via the u-restriction
(from the c_{r−1} step) and via the higher-stratum residual phases.

**Now here's the subtle complication:** the phase e_{p²}(−c_1·c_{r−2}) lives at modulus p²,
**NOT modulus p**. To Plancherel-collapse this over c_{r−2} ∈ Z/p (a single base-p digit),
we'd need the phase to be linear at MODULUS p (i.e., the e_{p²} dependence on c_{r−2}
viewed mod p).

Let's check: e_{p²}(−c_1·c_{r−2}) = e^{2πi·(−c_1·c_{r−2})/p²}. For c_{r−2} ∈ {0,...,p−1},
the phase argument is −c_1·c_{r−2}/p² ∈ (−c_1/p, 0]. The factor of 1/p² (vs 1/p) means
this phase is "slow" in c_{r−2} — varying by ≈ 2π·c_1/p² per unit c_{r−2}.

**Cannot apply Σ_{c_{r−2}=0}^{p−1} e_p(...) = p·δ directly** — the phase is at e_{p²},
not e_p.

There are several options:

### Option A: Promote c_{r−2} to a Z/p² variable

If c_{r−2} were over Z/p², the sum Σ_{c_{r−2}=0}^{p²−1} e_{p²}(−c_1·c_{r−2}) = p²·δ(c_1
≡ 0 mod p²). But c_{r−2} is only over Z/p (a single base-p digit) — we can't extend.

### Option B: Recognize that the c_{r−2} variable combines with the post-c_{r−1} u-restriction

The u-restriction after c_{r−1} collapse: u ≡ L̃_p^{-1}·c_1 mod p, u ∈ {0, ..., p^{r-1}−1}.
Parametrize u = u_0 + p·v, with u_0 = L̃_p^{-1}·c_1 mod p fixed and v ∈ {0, ..., p^{r-2}−1}.

The remaining Inner-Plancherel structure at the c_{r−2} level needs to be derived from
the Inner-after-c_{r−1} sum:

> T_p^{(1)} = p · Σ_v e_q(p · a_base · (u_0 + p·v))
>           = p · e_q(p·a_base·u_0) · Σ_v e_q(p²·a_base·v)
>           = p · e_q(p·a_base·u_0) · Σ_v e_{p^{r−1}}(a_base · v)

(used e_q(p²·x) = e_{p^{r−1}}(x).)

So T_p^{(1)} = p · e_{prefix}(...) · 1̂_{(r−1)}(a_base) where 1̂_{(r−1)} is the
length-p^{r−2} indicator-FT at modulus p^{r−1}. **This is itself a Dirichlet kernel at one
modulus level lower.**

This means the post-c_{r−1}-collapse sum has the SAME STRUCTURE as the original Inner(s_0)
at the r=3 case, but with:
- Modulus reduced by p (from p^{r+1} = q to p^{r−1})
- Length reduced by p (from p^{r-1} to p^{r-2})
- One fewer inner digit (c_{r−1} already summed out)

**This is the recursion! The chain naturally nests.**

### Option C: Reading the recursion carefully

After the c_{r−1} collapse:
> T_p (at level r) = (outer factor) · p · Σ_{remaining stuff}

The "remaining stuff" includes an outer phase from the lower-stratum Q_3, ..., Q_{r−1}
phases (now constants with respect to c_{r−1}=0 substitution), plus a length-p^{r−2}
Dirichlet structure on v, plus inner sums over c_2, ..., c_{r−2}.

**The "remaining stuff" is structurally identical to T_p at level r−1, but with:**
- N → p^{r-2} = N/p
- Inner digits: c_2, ..., c_{r−2} (count r−3, one less than the original r−2)
- Phase: P_a(s*(r−1)) with similar polynomial structure

So the chain is a **proper recursion** with depth r−2:
- Level r: peel off c_{r−1}, get factor p, reduce to level r−1
- Level r−1: peel off c_{r−2}, factor p, reduce to level r−2
- ...
- Level 3: peel off c_2 (the r=3 base case), factor p, reduce to length-p²/p = p outer
  sum
- Outer sum on c_1: bounded by p + log p (the cosecant identity)

**Total factor accumulated:** p^{r−2} (from the r−2 peelings) · (p + log p) (outer) = p^{r−1}·
(1 + log p / p) = N · (1 + log p / p) ≤ 2N.

**This gives |T_p| ≤ 2N strict at r ≥ 4 IF the recursion proceeds at each level with
the c_{r−1}-collapse mechanism.**

## The CRITICAL test: does the recursion actually work clean?

The c_{r−1} step gave T_p^{(1)} = (factor) · Σ_v e_{p^{r−1}}(a_base · v) · (residual c_2,
..., c_{r−2} phases at lower strata).

But the residual phases at lower strata (m = 3, ..., r−1) STILL CONTAIN POWERS OF c_{r−2},
c_{r−3}, ..., c_2, including:
- Linear-in-c_{r−2} terms at multiple strata (m = r−1, r — but m=r had c_{r−1}, which is gone)
- Cross-terms c_2·c_{r−2}, etc.
- Higher-degree terms.

**For the recursion to mimic the r=3 base case at the next level, we need:**

(R1) The post-c_{r−1} reduction yields a Dirichlet sum at modulus p^{r−1} (✓ shown above).
(R2) The residual phase factor at the top stratum (now p^{r−1} after the modulus reduction)
     is LINEAR in c_{r−2} with non-zero coefficient.
(R3) The lower-stratum phases at lower depths play the role of the lower strata at level
     r−1 (modulus-shifted by p).

(R1) ✓ by direct calculation.
(R2) Let me check carefully.

### Phase at the new top stratum (m = r−1, originally; now in modulus p^{r−1})

The phase factor at stratum m=r−1 in the original modulus p^{r+1} is e_{p²}(Q_{r−1}).
After the c_{r−1} collapse and modulus reduction p^{r+1} → p^{r−1}, this phase factor
INSIDE the post-collapse Dirichlet sum is... actually it's outside the Dirichlet sum.
Let me re-examine.

Going back to the c_{r−1} collapse calculation:
> Σ_{c_{r−1}} 1̂(p·a) · e_q(P_a) = (factor p) · Σ_{u ≡ u_0 mod p} e_q(p·a_base·u) · e_q(P_a |_{c_{r−1}=0})

The e_q(P_a |_{c_{r−1}=0}) is the phase at c_{r−1}=0 — the lower-stratum residual that
survives the c_{r−1}-zeroing.

But wait — the e_q(P_a) before the c_{r−1} sum already had ALL strata combined:
e_q(P_a) = e_{p²}(Q_3·p / 1) wait let me redo the modulus.

e_q(P_a) with P_a = Σ_m Q_m · p^m, q = p^{r+1}:
> e_q(P_a) = ∏_m e_{p^{r+1}}(p^m · Q_m) = ∏_m e_{p^{r+1−m}}(Q_m)

So at stratum m, the phase is e_{p^{r+1−m}}(Q_m). At m=2: e_{p^{r−1}}(Q_2). At m=3:
e_{p^{r−2}}(Q_3). ... At m=r: e_p(Q_r). At m=r+1: e_1 = trivial.

The c_{r−1} appears in Q_r linearly (coefficient −c_1), so the c_{r−1} dependence is via
e_p(−c_1·c_{r−1}).

After the c_{r−1} sum, this δ-collapses giving u ≡ u_0(c_1) mod p restriction. The OTHER
phase factors e_{p^{r−1}}(Q_2), e_{p^{r−2}}(Q_3), ..., e_{p²}(Q_{r−1}) are independent of
c_{r−1}, so they pull out of the c_{r−1} sum as a multiplicative factor.

What's left to sum over c_{r−2}? At this point, the post-c_{r−1} sum looks like:
> p · Σ_v e_{p^{r−1}}(a_base · v) · e_{p^{r−1}}(Q_2) · e_{p^{r−2}}(Q_3) · ... · e_{p²}(Q_{r−1})

We've also substituted u = u_0 + p·v, so the "Dirichlet" Σ_v e_{p^{r−1}}(a_base·v) lives
at modulus p^{r−1}. This is the analog of the r=3 case's Σ_u e_{p^3}(a·u) but with modulus
ONE p-power smaller.

**Now we want to do c_{r−2} sum.**

c_{r−2} appears in:
- Q_{r−1} at stratum m=r−1: linear in c_{r−2} with coefficient −c_1 mod p (modulus p² since
  phase is e_{p²}(Q_{r−1})).
- Q_r at stratum m=r: c_{r−2} appears linear in c_{r−2} with coefficient (−c_2 + c_1²/2)
  mod p (from j=2 (2, r-2) and j=3 (1,1, r-2) compositions). BUT this is at the c_{r−1}=0
  evaluation, since c_{r−1} got collapsed. So Q_r |_{c_{r−1}=0} still has c_{r−2} terms.

Hmm — actually wait. After c_{r−1} sum forces u_0 = L̃_p^{-1}·c_1 mod p, the residual
phase from e_p(Q_r |_{c_{r−1}=0}) gets pulled out, including the c_{r−2} contributions
at stratum m=r.

But Q_r |_{c_{r−1}=0} includes c_{r−2}·(−c_2 + c_1²/2) (assuming r ≥ 4). At r=4 this is
Q_4 |_{c_3=0} = −c_2²/2 + c_1²·c_2/2 − c_1^4/12. At r=4, c_{r−2}=c_2 is the
"second-top digit", and at stratum m=4 it appears with degree 2 (the c_2²/2 term)! That's
QUADRATIC, not linear.

Hmm. This is the c_2 at stratum m=4 case for r=4. Let me check r ≥ 5.

For r=5, c_{r−2}=c_3 at stratum m=r=5 with c_{r−1}=c_4 collapsed. Q_5 |_{c_4=0}:
> Q_5 (r=5) |_{c_4=0} = 0 + (terms not involving c_4 from the j=2 (1,4) composition that
> vanish) + c_2·c_3 contributions ... 

Re-reading Phase 2 Q_5 at r=5:
> Q_5 (r=5) = −c_1·c_4 − c_2·c_3 + c_1²·c_3/2 + c_1·c_2²/2 − c_1³·c_2/3 + c_1^5/20

Set c_4 = 0:
> Q_5 (r=5)|_{c_4=0} = −c_2·c_3 + c_1²·c_3/2 + c_1·c_2²/2 − c_1³·c_2/3 + c_1^5/20

In this, c_{r−2} = c_3 appears:
- −c_2·c_3 (linear in c_3, coefficient −c_2)
- c_1²·c_3/2 (linear in c_3, coefficient c_1²/2)

Total coefficient of c_3: (−c_2 + c_1²/2) mod p. LINEAR in c_3 at stratum m=r=5.

For r=4, c_{r−2}=c_2 in Q_4 |_{c_3=0}:
> Q_4 (r=4) = −c_1·c_3 − c_2²/2 + c_1²·c_2/2 − c_1^4/12
> Q_4 |_{c_3=0} = −c_2²/2 + c_1²·c_2/2 − c_1^4/12

c_2 appears with degree 2 (the −c_2²/2 term). **QUADRATIC, not linear.**

This is a r=4-specific behavior because at r=4, c_{r−2} = c_2 and the (k, k) cross at
m = 2k = 4 = r is HIT (whereas for r ≥ 5, c_{r−2} = c_3 and 2·(r−2) = 2(r−2) = 2r−4 ≠ r
in general — equality only at r=4).

**Let me check r=5, c_{r−2}=c_3, at stratum m=r=5: does c_3² appear?** c_3² requires
composition with 2·index_3 = 6, so c_3² at m=6, NOT at m=5. So at m=5, c_3 is LINEAR. ✓

**For r=6, c_{r−2}=c_4, at stratum m=r=6: does c_4² appear?** c_4² at m=8, > 6, so NO. c_4
at m=6 is LINEAR. ✓

**For r=7, c_{r−2}=c_5, at stratum m=r=7: c_5² at m=10 > 7, so NO. LINEAR.** ✓

**General r: at stratum m=r, c_{r−2} appears LINEARLY (since 2·(r−2) = 2r−4 > r iff r > 4).
The r=4 case is exceptional: c_2² hits stratum m=4 because 2·2 = 4 = r.**

So the r=4 case has a quadratic-in-c_2 obstruction at stratum m=r=4. For r ≥ 5, the
linear pattern holds at stratum m=r for c_{r−2}.

### What about stratum m=r−1 for c_{r−2} at r ≥ 5?

c_{r−2} at stratum m=r−1: appears via j=2 (1, r−2) → −c_1·c_{r−2} LINEAR.

At stratum m=r−1, modulus p² (since e_{p^{r+1−(r−1)}} = e_{p²}).

Coefficient: −c_1 mod p². But c_1 ∈ {0, ..., p−1}, so c_1 mod p² IS just c_1 (a digit).
The coefficient −c_1 mod p² has magnitude 0 or p−1 — could be 0 if c_1 = 0, non-zero
otherwise.

**c_{r−2} appears at TWO strata (m = r−1 and m = r), both LINEARLY, with different
moduli and different coefficients.**

To sum over c_{r−2} ∈ Z/p, the combined phase in c_{r−2} is:
> e_{p²}(−c_1·c_{r−2}) · e_p(coefficient(c_2, c_1²)·c_{r−2})

Combine into a single modulus-p² phase: e_{p²}(−c_1·c_{r−2} + p·(coefficient·c_{r−2}))
= e_{p²}((−c_1 + p·coefficient)·c_{r−2}).

Sum: Σ_{c_{r−2}=0}^{p−1} e_{p²}((−c_1 + p·coefficient)·c_{r−2}) 

This is a length-p sum of a phase at modulus p². As c_{r−2} ranges 0..p−1:

> Σ = (e_{p²}(p·(−c_1 + p·coefficient)) − 1) / (e_{p²}(−c_1 + p·coefficient) − 1)
>   = (e_p(−c_1 + p·coefficient) − 1) / (e_{p²}(−c_1 + p·coefficient) − 1)

If −c_1 + p·coefficient ≡ 0 mod p (i.e., c_1 ≡ 0 mod p, since p·coefficient ≡ 0 mod p),
the numerator e_p(0) − 1 = 0, but the denominator is also small. By the same calculation
as r=2's |1̂(p·a)| analysis:

> |Σ_{c_{r−2}} e_{p²}((−c_1 + p·...)·c_{r−2})| = sin(π/p) / sin(π·(−c_1 + p·...)/p²)
>                                              ≈ p (when c_1 ≡ 0 mod p)
>                                              ≈ p/c_1 (when c_1 ≠ 0 mod p)

**This is a Dirichlet kernel — not a delta!** The c_{r−2} sum at level r−1's "top stratum"
gives a Dirichlet sum, not a clean Plancherel collapse.

**This is the key structural difference from the r=3 base case.**

In the r=3 base case: the c_2 sum at stratum m=r=3 was at modulus e_p(−c_1·c_2) — phase
modulus p, length p, exact δ via Σ_{c_2} e_p(−c_1·c_2) = p·δ(c_1 ≡ 0 mod p). Clean.

At r ≥ 4, when we get to the c_{r−2} level (the second peel), the relevant stratum is
m = r−1 (modulus p²) AND m = r (modulus p, surviving as residual phase). The combined
c_{r−2} phase is at modulus p² (the higher one), so the length-p c_{r−2} sum DOESN'T
fully resolve.

### What this means for the recursion

The clean recursion fails at the c_{r−2} step (for r ≥ 4) because the phase modulus is
larger than the digit range.

**However**, the Dirichlet sum can still be bounded — it gives |Σ| ≤ p (or |Σ| ≈ 1 for
generic c_1) and after summing over c_1 (the outer variable, eventually), the total bound
becomes a multi-Dirichlet sum that can be controlled by the cosecant identity ITERATED.

**This is the analog of the r=3 case but pushed up one level:** instead of one cosecant
identity at the outer, we get a NESTED cosecant identity at each level of the chain.

### Closed-form examination of the r=4 c_2 sum (special)

At r=4, after c_3 collapse, the c_2 sum at stratum m=4 has the −c_2²/2 quadratic. The
combined c_2 phase:
- At m=3 (modulus p²): e_{p²}(−c_1·c_2) linear in c_2.
- At m=4 (modulus p): e_p(−c_2²/2 + c_1²·c_2/2) — QUADRATIC and linear pieces in c_2.

Combined modulus is p²; rewriting at modulus p²:
> e_{p²}(−c_1·c_2 + p·(−c_2²/2 + c_1²·c_2/2))
> = e_{p²}((−c_1 + p·c_1²/2)·c_2 + p·(−c_2²/2))
> = e_{p²}((−c_1 + p·c_1²/2)·c_2) · e_p(−c_2²/2)

The c_2² piece is at modulus p (it doesn't get a 1/p reduction inside e_{p²}). Together:
> Σ_{c_2=0}^{p−1} e_p(−c_2²/2) · e_{p²}((−c_1 + p·c_1²/2)·c_2)

This is a length-p quadratic-phase × linear-phase sum. It's NOT a clean Gauss sum
(modulus mismatch), NOR a clean Dirichlet sum (the quadratic e_p(−c_2²/2) couples to the
sum). It needs separate analysis.

**For r=4 specifically, the c_2 step is structurally MORE COMPLEX than the analogous
r=3 c_2 step.** This breaks the clean "1 Plancherel per digit" pattern.

### General-r pattern recognition

The structural issue at the c_{r−2} step:
- Stratum m = r−1 contributes e_{p²}(−c_1·c_{r−2}), forcing the sum to a length-p
  geometric series at modulus p² (Dirichlet kernel, factor ≤ p).
- Stratum m = r contributes (e_p of polynomial in c_{r−2}) at modulus p.
- These combine to a phase at modulus p² with non-trivial structure.

**The chain DOES extend, but each peeling step at level k (for k < r−1) generates a
Dirichlet-kernel-like factor instead of a delta-collapse.**

Counting factors:
- Top peel (c_{r−1}): clean delta, factor p save. Total: p.
- Second peel (c_{r−2}): Dirichlet sum at modulus p², factor ≤ p (NOT delta). Total: p.
- Third peel (c_{r−3}): even worse — phase at modulus p^3 (stratum m=r−2 has modulus p^3).
- ...

Hmm — wait. Each successive peel sees a HIGHER modulus phase from its "earliest" stratum.
The c_k phase first appears at stratum m=k+1 with modulus e_{p^{r+1−(k+1)}} = e_{p^{r−k}}.
So:
- c_{r−1}: first appears at m=r, modulus e_p (modulus p). LINEAR in c_{r−1} at this modulus.
  Length-p sum at modulus p → clean delta.
- c_{r−2}: first appears at m=r−1, modulus e_{p²}. LINEAR at this modulus. Length-p sum
  at modulus p² → Dirichlet kernel (factor ≤ p, NOT delta).
- c_{r−3}: first appears at m=r−2, modulus e_{p^3}. Length-p sum at modulus p^3.
- ...
- c_2: first appears at m=3, modulus e_{p^{r−2}}. Length-p sum at modulus p^{r−2}.

**This is a CLEAN STRUCTURAL HIERARCHY:** the further into the chain we peel, the higher
the modulus of the dominant phase relative to the digit length. The "delta collapse" only
works at the TOP digit; subsequent digits give Dirichlet kernels of progressively higher
modulus.

## Reconciliation with HENSEL_APPROACH_A's claim

HENSEL_APPROACH_A.md claimed (with caveats) that the chain saves factor p at each digit
and produces strict 2√N at r ≥ 4.

After this careful re-derivation, that claim was OPTIMISTIC. The clean δ-collapse only
works at the top digit c_{r−1}. The second digit c_{r−2} already gives a Dirichlet
kernel at modulus p², not a delta.

**However:** the Dirichlet kernel has magnitude ≤ p (or smaller for non-zero c_1). The
chain still saves "approximately p per digit" in the AVERAGE sense, but not in the
WORST-CASE bound — the worst-case Dirichlet kernel value at c_1 ≡ 0 mod p is p, the same
as the delta save.

So in BOUND terms, the chain still saves factor p per peel (because Dirichlet kernel ≤ p),
but we need to bound the SUM over the outer variables more carefully, treating each
inner peel as a Dirichlet sum rather than a delta.

### Refining the bound

After c_{r−1} peel: T_p^{(1)} = p · Σ_v e_{p^{r−1}}(a_base·v) · (residual phases).

After c_{r−2} peel: T_p^{(2)} ≈ p · (Dirichlet at modulus p² from c_{r−2} sum) · Σ_v ... · ...

The Dirichlet sum from c_{r−2} sum is ≤ p in magnitude, but it doesn't restrict the v
sum (the way the δ at c_{r−1} restricted u → u_0 + p·v). Instead, it produces a phase
factor and a magnitude factor.

**The c_{r−2}, c_{r−3}, ..., c_2 sums each give a Dirichlet kernel factor in their
respective moduli, but they accumulate.**

The recursion structure becomes:
> |T_p| ≤ |outer over c_1| · ∏_{k=2}^{r−1} |Dirichlet kernel from c_k peel at modulus p^{r+1−(k+1)}|
> ≈ |outer| · ∏ p (worst case)
> = (p + log p) · p^{r-2}
> = N · (1 + log p / p)
> ≤ 2N

**Wait — this gives 2N back!** Each digit's Dirichlet kernel contributes ≤ p (the maximum
of a length-p Dirichlet kernel is p, achieved when the modulus phase argument ≡ 0 mod the
denominator).

So the WORST-CASE bound IS 2N. The chain extends. But the constants are not all the same
as r=3:

- At r=3: one Plancherel save (clean δ × cosecant outer). Constant ≈ 1 + log p / p.
- At r ≥ 4: r−2 Dirichlet-kernel saves (each ≤ p) × cosecant outer. Total constant ≈ 1 +
  log p / p (assuming each Dirichlet kernel attains its max at the same outer-c_1 value
  — which is the worst case).

**The WORST-CASE bound at all r is 2N.** ✓

**But the AVERAGE constant might be larger or smaller depending on how the Dirichlet
kernels interact.** This is where the rigorous-vs-near-rigorous distinction lives.

## Critical assessment

After this careful analysis, here's where I land:

**The chain DOES extend.** Each higher digit c_k has a linear-leading-coefficient structural
feature (I2 from Phase 1) at stratum m=k+1. Plancherel-like collapses apply at each level,
saving factor p (or ≤ p, in worst case) per digit.

**BUT the moduli differ at each level.** Only the TOP digit c_{r−1} gives a clean delta
(modulus p, matching the length-p digit range). All deeper digits give Dirichlet kernel
sums at progressively higher moduli (p², p³, ..., p^{r−2}).

**The bound assembly is:**
- Per-digit save (worst-case): factor p (= length-p Dirichlet kernel max).
- Outer cosecant on c_1: factor (p + log p) ≤ 2p (uniform in p).
- TOTAL: |T_p| ≤ p^{r-2} · (p + log p) = N · (1 + log p / p) ≤ 2N strict at r ≥ 4.

**This matches the r=3 bound shape exactly** — strict 2N (hence strict 2√N for |S_partial|).

The DIFFERENCE from r=3 is structural (Dirichlet kernels at progressively higher moduli),
but the BOUND is the same constant.

## Adversarial check: is the "Dirichlet kernel ≤ p" claim sharp?

A length-p sum Σ_{c=0}^{p−1} e_{p^L}(α·c) for arbitrary α has magnitude:
- = p when α ≡ 0 mod p^{L−1} (i.e., α = p^{L−1}·β with β integer, then the sum is
  Σ e_p(β·c) = p·δ(β ≡ 0 mod p)).
- = sin(π·α/p^{L−1}) / sin(π·α/p^L) ≤ csc(π·α/p^L) (Dirichlet kernel max), bounded by
  p^L/|α| / something ≈ p/|α/p^{L−1}|.

For generic α not divisible by high powers of p, the magnitude is O(1). For α ≡ 0 mod
p^{L−1}, magnitude is p.

In our case, α = −c_1 + p · (... lower-order corrections ...). The dominant behavior
depends on c_1. For c_1 ≢ 0 mod p, α is a unit times 1, so |α/p^{L−1}| ≈ 1/p^{L−1}, and
magnitude ≈ p^L/1 = p^L (wait this would be > p which is wrong).

Let me redo for L=2 (the c_{r−2} case): α = −c_1 + p·correction. 
- If c_1 ≢ 0 mod p: |α| ~ 1 (unit), |α/p²| ~ 1/p². sin(π·α/p²) ≈ π/p². sin(π·α/p) ≈ π·α/p ≈ π·c_1/p (a unit/p). Ratio: (π·c_1/p) / (π/p²) = c_1·p. So |Σ| ≈ c_1·p. 

That's bigger than p for c_1 > 1!

Hmm wait — Dirichlet kernel formula: Σ_{c=0}^{p−1} e_M(α·c) = (e_M(α·p) − 1)/(e_M(α) − 1)
= (e_{M/p}(α) − 1)/(e_M(α) − 1) when p divides M.

For M = p², |Σ| = sin(π·α/p)/sin(π·α/p²).

For α = c_1 with c_1 ∈ {1, ..., p−1}: sin(π·c_1/p)/sin(π·c_1/p²) ≈ (π·c_1/p)/(π·c_1/p²)
= p. **So magnitude is exactly p, not c_1·p.**

OK I had a calc error. Let me redo: sin(π·c_1/p)/sin(π·c_1/p²) ≈ (π·c_1/p)/(π·c_1/p²) = p^2/p = p. Yes, p.

So |Σ_{c_{r−2}} e_{p²}(c_1 · c_{r−2})| ≈ p for c_1 ∈ {1,...,p−1}.

For c_1 = 0: α = p · correction ~ p. sin(π·p/p)/sin(π·p/p²) = sin(π)/sin(π/p) = 0 / (π/p)
= 0. So |Σ| = 0 for c_1 = 0 (NO contribution, not p).

Wait that's surprising. Let me check: if α = 0 (c_1 = 0 and correction = 0), the sum is
Σ_{c=0}^{p−1} e_{p²}(0) = p. So magnitude p, not 0. Conflict with the above.

The issue: when α = 0 exactly, the sum is p. When α = p · (non-zero unit) (i.e., α is a
non-zero multiple of p but not p²), the sum is 0.

For our case (c_1 = 0), α = p·correction. If correction ≠ 0 mod p, then α = p·unit, so
|Σ| = 0. If correction = 0 mod p, then α ≡ 0 mod p², so |Σ| = p.

So the c_{r−2} sum magnitude is:
- c_1 ≠ 0 mod p (typical): |Σ| ≈ p (each c_1 gives the same large value).
- c_1 = 0 mod p: |Σ| = p (if correction ≡ 0 mod p) or 0 (if correction ≠ 0 mod p).

For c_1 ∈ {0, 1, ..., p−1} ranging over Z/p:
- p−1 values with c_1 ≠ 0: each gives |Σ| ≈ p.
- 1 value with c_1 = 0: |Σ| is either p or 0 (depending on correction).

The TOTAL outer sum over c_1 of |Σ| accumulates to ≤ (p−1)·p + p = p² = N.

Hmm but we still have OTHER inner digits to peel. Let me think about this more carefully
in Phase 4.

## Output of Phase 3

**Each higher digit c_k (k = 2, 3, ..., r−1) appears LINEARLY at its first stratum
m = k+1.** The structural ingredient (I2 from Phase 1) extends.

**The Fourier-collapse at each c_k saves at most factor p (the length-p Dirichlet kernel
max), BUT the modulus is p^{r−k} at the c_k step, NOT p as in the r=3 base case.**

Consequence: the "clean δ-collapse" only happens at the TOP digit c_{r−1}. The deeper
digits give Dirichlet-kernel reductions instead of pure deltas.

**Bound impact:** the worst-case bound is still |T_p| ≤ 2N strict at r ≥ 4 (each Dirichlet
kernel ≤ p), but the chain structure is more complex than at r=3 (nested Dirichlet kernels
at progressively higher moduli rather than a single δ).

**Caveats:**
- At r=4 specifically, c_{r−2} = c_2 has a QUADRATIC term at stratum m=r=4 (the c_2²/2),
  which couples to the c_2 Plancherel and produces a length-p quadratic-Gauss-mixed sum,
  NOT a clean Dirichlet kernel. This is an r=4-specific complication.
- At r ≥ 5, the c_{r−2} sum is a "clean" length-p Dirichlet kernel at modulus p² (no
  quadratic-in-c_{r−2} obstruction).
- At deeper digits (c_{r−3}, etc.), the moduli get higher (p^3, p^4, ...), and the
  Dirichlet kernels accumulate. Worst case still ≤ p per peel, but the iteration needs care.

## Status of the three pre-registered hypotheses

**H_CHAIN_EXTENDS_STRICT:** Provisionally supported. The Dirichlet kernel max at each peel
is ≤ p, summed against the outer cosecant the total is ≤ 2N. This is the same constant as
r=3. **However, the chain mechanism is different — Dirichlet kernels rather than deltas.**
The "strict 2√N" lands by the worst-case Dirichlet kernel attaining its max ≤ p, not by
clean deltas.

**H_CHAIN_EXTENDS_LOOSER:** Possible if the Dirichlet kernel doesn't attain p uniformly,
making the bound 2√p·√N or similar. Need careful r=4 analysis (quadratic-in-c_2 stratum-4
complication) and deeper-digit nested cosecant identity.

**H_CHAIN_FAILS:** Not triggered. Each c_k has linear leading appearance; collapses
exist.

**Currently leaning toward H_CHAIN_EXTENDS_LOOSER for r=4 (due to quadratic c_2 at m=r),
and H_CHAIN_EXTENDS_STRICT for r ≥ 5 (clean linear pattern).** The Phase 4 nesting works
out the precise constant.

## Files

- CHAIN_PHASE1_R3_RECAP.md — r=3 chain (delta collapse)
- CHAIN_PHASE2_DIGIT_EXPANSION.md — stratum tables
- This document — per-digit collapse test
- HENSEL_APPROACH_A.md — Wilson's prior claim of strict 2√N at r ≥ 4

## Next: Phase 4 assembles the nested chain bound

The actual nested chain bound depends on:
- Worst-case Dirichlet kernel magnitudes at each level (≤ p per).
- The r=4 quadratic-c_2 complication at stratum m=r=4.
- The cosecant identity on the OUTER c_1 sum (the r=3 base case's identity).

Phase 4 carries out this calculation explicitly and produces a constant.
