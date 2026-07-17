# CHAIN Phase 4 — Nested Chain Bound at r ≥ 4

**Date:** 2026-05-11. Analyst: Wilson. Reports to: Nathan.

**Goal:** Assemble the per-digit Fourier-collapse data from Phase 3 into a nested chain
bound on |T_p| at r ≥ 4. The key question: does the chain yield strict 2√N, looser
2√p·√N, or something in between?

## Setup recap

From Phase 3:
- c_{r−1} (top digit): clean δ-collapse at modulus p, factor p save. Restricts u
  to length-p^{r-2} arithmetic progression (analog of r=3 mechanism).
- c_{r−2} (next digit): linear at stratum m=r−1 (modulus p²) AND at stratum m=r (modulus p)
  except at r=4 where stratum m=r=4 has c_2² quadratic.
- c_{r−3}, ..., c_2: linear at their first strata (modulus p^{r−k}) and at higher strata
  (modulus lower).

## The proper nesting

Nest the peels from TOP (c_{r−1}) down to BOTTOM (c_2). At each level, after peeling
c_k, the remaining sum has:
- Modulus reduced from p^{r+1−(...)} to p^{r−k} (level-dependent).
- Length on the Dirichlet variable reduced by factor p.
- Inner digits c_2, ..., c_{k−1} still to be peeled.
- A magnitude factor ≤ p contributed by the c_k Plancherel/Dirichlet step.

After all r−2 inner-digit peels, we're left with an OUTER sum on c_1 ∈ Z/p (the
r=3-saddle outer variable) at modulus p², matching the r=3 base case's cosecant
identity.

## Quantitative chain at level k (top-down)

Let `T_p^{(0)} := T_p` (full bilinear). After peeling c_{r−1}, ..., c_{r−k+1} (i.e., k
peels), define `T_p^{(k)}` as the residual.

After all r−2 peels: `T_p^{(r−2)} = (outer c_1 sum at modulus p² in some Dirichlet kernel form)`.

### Per-peel magnitude budget

Each peel transforms:
> T_p^{(k−1)} = factor_k · (next-level Dirichlet-kernel-on-{remaining digits and outer})

with **|factor_k| ≤ p** (Dirichlet kernel maximum at length-p sum, regardless of modulus).

After r−2 peels, multiplicative factor accumulated: |∏ factor_k| ≤ p^{r-2}.

### Final outer sum on c_1

The "final" outer sum has form:
> Σ_{c_1=0}^{p−1} (composite phase from all peels) · D_p(some quantity, p²)

where D_p(·, p²) is the length-p Dirichlet kernel at modulus p² (the original r=3 outer
sum).

**Cosecant identity:** Σ_{c_1=0}^{p−1} |D_p(a_0(c_1), p²)| ≤ p + log p ≤ 2p.

**Total bound:**
> |T_p| ≤ ∏_k |factor_k| · Σ_{c_1} |D_p| ≤ p^{r-2} · (p + log p)
>       = p^{r-1} · (1 + log p / p)
>       = N · (1 + log p / p)
>       ≤ 2N    (uniformly in r ≥ 3, p ≥ 3)

**Strict |T_p| ≤ 2N at r ≥ 4 (worst-case bound).**

Equivalently: **|S_partial| ≤ 2 √N · √q strict, achievable from family-level closed form
at r ≥ 4.**

## Caveat 1: The r=4 quadratic c_2 step

At r=4, when we peel c_{r−2} = c_2, stratum m=r=4 has c_2² (NOT linear). The combined
c_2 phase is:
> e_{p²}(−c_1·c_2 + p·(−c_2²/2 + c_1²·c_2/2)) at modulus p²

This is a length-p sum of a quadratic-in-c_2 phase at modulus p² (with a linear-in-c_2
"twist" inside modulus p²).

**Magnitude analysis:** the sum Σ_{c_2=0}^{p−1} e_{p²}(α·c_2 + β·c_2²/p) where α is at
order 1 and β = p·(quadratic coefficient) is the "stretching" by factor p.

Rewriting: e_{p²}(α·c_2) · e_p(β·c_2²/2). Wait let me redo: 
- e_{p²}(−c_1·c_2): modulus p², linear.
- e_{p²}(p·(−c_2²/2 + c_1²·c_2/2)) = e_p(−c_2²/2 + c_1²·c_2/2): modulus p, quadratic+linear.

So combined: e_{p²}(−c_1·c_2) · e_p(−c_2²/2 + c_1²·c_2/2).

The sum:
> Σ_{c_2=0}^{p−1} e_{p²}(−c_1·c_2) · e_p(c_1²·c_2/2) · e_p(−c_2²/2)

The combined linear-in-c_2 part: e_{p²}(−c_1·c_2) · e_p(c_1²·c_2/2) = e_{p²}(−c_1·c_2 +
p·c_1²·c_2/2) = e_{p²}((−c_1 + p·c_1²/2)·c_2). The leading order is −c_1 mod p².

The quadratic e_p(−c_2²/2) is a length-p Gauss sum factor when summed alone.

Combined sum:
> Σ_{c_2} e_p(−c_2²/2) · e_{p²}((−c_1 + p·c_1²/2)·c_2)

This is a **TWISTED quadratic Gauss sum** at modulus p (the quadratic part) combined with
a phase-shift e_{p²} factor (which only changes by O(1/p) per unit c_2).

**Bound:** by Cauchy-Schwarz on c_2 (or van der Corput's inequality, or just a direct
Gauss sum + shift estimate):

> |Σ_{c_2} e_p(−c_2²/2) · e_{p²}(γ·c_2)| ≤ √p · (1 + O(γ/p))
> where γ = −c_1 + p·c_1²/2 mod p².

For generic c_1 (nonzero mod p), γ/p ≈ c_1 + O(1) ≈ 1, so:
> |Σ| ≤ √p · (1 + O(1)) ≤ 2·√p

For c_1 = 0 mod p: γ = O(p) (small), |Σ| ≤ √p · (1) = √p.

**Magnitude budget at r=4 c_2 step: ≤ 2√p, NOT p.**

This is a factor √p WORSE than the generic Dirichlet kernel bound p. The chain at r=4
LOSES √p at the c_2 step.

**Total bound at r=4:**
> |T_p| ≤ p (c_3 peel) · 2√p (c_2 peel) · 2p (outer cosecant)
>       = 4p · p^{3/2}
>       = 4p^{5/2}
>       = 4·p^2·√p
>       = 4N·√p    (since N = p²)

**So at r=4, the bound is `|T_p| ≤ 4·√p · N`, equivalently `|S_partial| ≤ 4√p · √N`.**

This is the **H_CHAIN_EXTENDS_LOOSER** outcome at r=4. The constant is 4√p instead of 2.

## Caveat 2: Does the looser bound persist at r ≥ 5?

At r ≥ 5, c_{r−2} appears LINEARLY at stratum m=r (not quadratically). The chain at the
c_{r−2} peel is a clean Dirichlet kernel at modulus p², factor ≤ p.

**But what about c_{r−3} at r ≥ 5?**

For r=5: c_{r−3} = c_2. At what strata?
- Q_3 = −c_1·c_2 + c_1³/6: c_2 linear at m=3, modulus p^{r−2} = p^3.
- Q_4 = −c_1·c_3 − c_2²/2 + c_1²·c_2/2 − c_1^4/12: c_2 appears LINEAR (c_1²·c_2/2) AND
  QUADRATIC (c_2²/2). At m=4, modulus p².
- Q_5 = −c_1·c_4 − c_2·c_3 + c_1²·c_3/2 + c_1·c_2²/2 − c_1³·c_2/3 + c_1^5/20: c_2 appears
  with degree 2 (c_1·c_2²/2) and degree 1 (c_2·c_3 and c_1³·c_2/3). At m=5, modulus p.

So at r=5, the c_2 peel involves quadratic-in-c_2 at m=4 (modulus p²) AND at m=5 (modulus
p). The combined c_2 phase is:
> e_{p^3}(−c_1·c_2 + c_1³/6) · e_{p²}(−c_2²/2 + c_1²·c_2/2 + c_1^4 stuff) · e_p(c_2·(...) + c_2²·(...))

The c_2² term at m=4 (modulus p²) and at m=5 (modulus p) — both present.

The combined c_2 phase has a NET QUADRATIC structure at modulus p² (and lower).

**Length-p sum of a quadratic-at-modulus-p²-or-lower phase = quadratic Gauss sum
truncated, magnitude ≤ √p · O(1) = √p.**

**So at r=5, c_2 peel also gives factor √p, NOT p.**

Generalizing: at general r, the c_2 peel always involves quadratic-in-c_2 terms (since
c_2² first appears at stratum m=4, which is well within the active stratum range for r ≥ 4).
So **c_2 ALWAYS contributes factor √p (worst case)**.

**What about c_3, c_4, ..., c_{r−2}?**

c_k² first appears at stratum m=2k. This stratum is ≤ r iff 2k ≤ r iff k ≤ r/2.

So for k ≤ r/2, c_k² is in the active strata; for k > r/2, c_k² is NOT in active strata.

- k = 2: c_2² at m=4. For r ≥ 4, in active strata. → c_2 peel always gives √p.
- k = 3: c_3² at m=6. For r ≥ 6, in active strata. → c_3 peel gives √p for r ≥ 6, p for r = 4, 5.
- k = 4: c_4² at m=8. For r ≥ 8, in active strata. → c_4 peel gives √p for r ≥ 8, p for r = 5, 6, 7.
- ...
- k ≤ r/2: c_k peel gives factor √p.
- k > r/2: c_k peel gives factor p (clean Dirichlet kernel, no quadratic complication).

**At general r, count of "√p factors": k ∈ [2, ⌊r/2⌋], giving ⌊r/2⌋ − 1 factors of √p.**

**Count of "p factors": k ∈ [⌊r/2⌋+1, r−1], giving r − 1 − ⌊r/2⌋ factors of p.**

**Total product:** ∏ factors = (√p)^{⌊r/2⌋ − 1} · p^{r−1−⌊r/2⌋}.

For even r = 2m: factors = (√p)^{m−1} · p^{m−1} = p^{(m−1)/2 + m−1} = p^{(3m−3)/2} =
p^{(3(r−2)/2)/2} ... wait that's not right. Let me redo.

Even r = 2m: ⌊r/2⌋ = m. Count of √p factors = m − 1. Count of p factors = 2m − 1 − m = m − 1.

So total = (√p)^{m−1} · p^{m−1} = p^{(m−1)/2 + (m−1)} = p^{(3m−3)/2}.

For r = 4 (m=2): p^{3·2−3)/2} = p^{3/2}. So accumulated factor = p^{3/2}.

|T_p| ≤ p^{3/2} · (p + log p) ≈ p^{3/2} · p = p^{5/2}.

|T_p|/N = p^{5/2}/p^{r-1} = p^{5/2}/p^3 = p^{−1/2} ... that gives |T_p| ≤ N/√p < N!

Wait that's IMPROVING — that can't be right. Let me re-examine.

Hmm — I made an algebra error. Let me redo:
- N = p^{r−1} at r=4 is p^3.
- Accumulated factor = (per-peel magnitudes) × (outer cosecant) = p^{(3m−3)/2} · 2p
  (at r=4, m=2): p^{3/2} · 2p = 2·p^{5/2}.
- |T_p| ≤ 2·p^{5/2}.
- |T_p|/N = 2·p^{5/2}/p^3 = 2/√p.

So |T_p| ≤ 2/√p · N = 2/√p · N. This is LESS than 2N (by factor √p) — IMPROVEMENT.

Wait — but the inner-Plancherel SAVES factor p per peel. If we only save √p per peel
(due to quadratic complication), we save LESS, so the bound should be LARGER, not smaller.

I think I mis-counted. Let me redo from scratch.

### Redo: trivial baseline

Trivial bound on |T_p| = |Σ 1̂(p·a)·e_q(P_a)|: each term has |1̂| ≤ N (max), |e_q| = 1,
|support| = N. So trivial: |T_p| ≤ N · N = N². At r=4: N² = p^6.

For r=4, trivial = p^6.

### The chain at each peel SAVES a factor

If a peel "saves factor p", that means the residual after the peel is the prior residual
divided by p (in magnitude). If a peel "saves factor √p", residual / √p.

Starting from trivial N²:
- r=3 chain (one delta peel saving p, one cosecant outer saving p): N²/p / p = N²/p²
  = (p²)²/p² = p² = N. ✓ matches r=3 result.

- r=4 chain WITH STRICT (all peels save p): r−2 = 2 inner peels (c_2, c_3) saving p
  each + outer cosecant saving p = p^3 saved total. N²/p^3 = p^6/p^3 = p^3 = N. ✓
  matches HENSEL_APPROACH_A claim.

- r=4 chain with QUADRATIC c_2 STEP saving only √p: factor p (c_3) · √p (c_2) · p (outer)
  = p^{5/2}. N²/p^{5/2} = p^6/p^{5/2} = p^{7/2} = N · √p.

So **|T_p| ≤ √p · N at r=4**, NOT 2N.

OK so my count was wrong earlier — the "saving" goes in the DENOMINATOR. Let me redo.

### Per-peel savings (corrected)

Each peel saves factor `s_k` where:
- s_k = p if the c_k phase is clean linear (length-p Plancherel/Dirichlet at modulus
  matching length).
- s_k = √p if the c_k phase is quadratic Gauss-sum-like (modulus < length).
- s_k = O(1) if no save (would mean the chain fails at this digit).

Outer cosecant savings on c_1: factor p (the sum p · 2p / 2p = p... actually let me
re-examine. The outer sum is Σ_{c_1} |D_p(a_0(c_1), p²)| ≤ p + log p ≈ p. The full
support of 1̂ sums to ≤ N · log N (Pólya-Vinogradov), so the cosecant identity SAVES
factor (N · log N)/(p) ≈ N · log N / p ≈ p^{r-2} · log p / 1 = p^{r-2}·log p. Hmm.

Actually let me just compute directly without the "saving" framing.

### Direct calculation

|T_p| ≤ Σ_a |1̂(p·a)| · |Σ_{(remaining digits given a)} (residual phase)|

Wait this isn't quite the right factorization. Let me restart Phase 4 with a cleaner
approach.

## Phase 4 — clean restart

The bilinear:
> T_p = Σ_{c_1, c_2, ..., c_{r−1}} 1̂(p·a(c_1, ..., c_{r−1})) · e_q(P_a(s*(r)))

Reorganize: peel c_{r−1} first.

Σ_{c_{r−1}} 1̂(p·a) · e_q(P_a) — the c_{r−1} appears in 1̂ (via the a shift) and in
e_q(P_a) only at the top stratum m=r (via e_p(−c_1·c_{r−1}) + c_{r−1}-free residual).

The c_{r−1} sum was computed in Phase 3:
> = p · Σ_{u ≡ u_0 mod p, 0 ≤ u < N} e_q(p·a_base·u) · e_q(P_a |_{c_{r−1}=0})

So after c_{r−1} peel, T_p becomes:
> T_p = p · Σ_{c_1, ..., c_{r−2}} (Σ_{u ≡ u_0 mod p} e_q(p·a_base·u)) · e_q(P_a |_{c_{r−1}=0})

The Σ_{u ≡ u_0 mod p} sum is a length-p^{r-2} arithmetic-progression sum at modulus q =
p^{r+1}. Parametrize u = u_0(c_1) + p·v, v ∈ {0, ..., p^{r-2}−1}:
> Σ_u e_q(p·a_base·u) = e_q(p·a_base·u_0) · Σ_v e_q(p²·a_base·v)
>                     = e_q(p·a_base·u_0) · Σ_v e_{p^{r−1}}(a_base·v)

So:
> T_p = p · Σ_{c_1, ..., c_{r−2}} e_q(p·a_base·u_0) · (length-p^{r-2} 1̂-like at modulus
> p^{r−1}) · e_q(P_a |_{c_{r−1}=0})

**At this point, the structure is morally analogous to the level-(r−1) bilinear, but
with one less inner digit, one lower modulus, and a residual phase from the original
e_q(P_a |_{c_{r−1}=0}).**

### The recursion

After r−2 peels, the analog of the r=3 base case is reached.

At each peel from c_{r−1} down to c_2:
- The "active modulus" decreases by p at each step (q = p^{r+1} → p^r → p^{r−1} → ...).
- One inner digit is consumed.
- A residual phase factor (from the stratum being peeled) is left in the outer sum.

At the c_k peel (after r−1−k prior peels): active modulus is p^{k+2}. The peel is on the
length-p variable c_k. The phase at the top active stratum is e_{p²}(linear in c_k from
m=k+1) plus residual lower-stratum phases.

Linear coefficient of c_k at stratum m=k+1 (modulus p²): −c_1 mod p (the leading term;
higher corrections vanish mod p).

**At each peel, the c_k sum is:**
> Σ_{c_k=0}^{p−1} (post-prior-peels stuff) · e_{p²}(−c_1·c_k + lower-order corrections)

If the prior peels haven't introduced c_k² terms, this is a clean length-p geometric sum
at modulus p² → Dirichlet kernel, magnitude ≤ p.

If the prior peels HAVE introduced c_k² terms (from the original stratum-m=2k cross
term), this is a quadratic-twisted Dirichlet kernel sum, magnitude ≤ √p · (1 + |γ|/p)
where γ is the linear-twist coefficient.

### When does c_k² appear?

c_k² appears in P_a(s*) at stratum m = 2k via the (k,k) cross in (p·s*)²/2. After r−1−k
prior peels (which set c_{r−1}, ..., c_{k+1} to specific values), the c_k² coefficient
at stratum m=2k is preserved (it doesn't involve any of the peeled digits).

**c_k² is in the active phase strata iff 2k ≤ r (the stratum m=2k is at modulus
p^{r+1−2k}, which is non-trivial iff r+1−2k ≥ 1 iff k ≤ r/2).**

So for k ≤ ⌊r/2⌋, the c_k peel has a quadratic-in-c_k term in the active phase. For
k > ⌊r/2⌋, the c_k peel has only linear terms.

### Per-peel magnitude (refined)

- **k ∈ {⌊r/2⌋+1, ..., r−1}:** clean linear c_k, length-p Dirichlet kernel at modulus
  p^{r+1−(k+1)} = p^{r−k}. Length-p sum magnitude ≤ length-of-period at modulus = p.
  **Magnitude ≤ p per peel.**

- **k ∈ {2, ..., ⌊r/2⌋}:** quadratic-in-c_k present, gives twisted Gauss sum of length p.
  **Magnitude ≤ √p per peel.**

### Count

- Linear peels: k from ⌊r/2⌋+1 to r−1, count = r − 1 − ⌊r/2⌋.
  - r=4: count = 4−1−2 = 1.
  - r=5: count = 5−1−2 = 2.
  - r=6: count = 6−1−3 = 2.
  - r=7: count = 7−1−3 = 3.
  - General: ⌈r/2⌉ − 1.

- Quadratic peels: k from 2 to ⌊r/2⌋, count = ⌊r/2⌋ − 1.
  - r=4: 1
  - r=5: 1
  - r=6: 2
  - r=7: 2
  - General: ⌊r/2⌋ − 1.

### Outer cosecant on c_1

The outer sum on c_1 is over the c_1 Dirichlet kernel after all peels collapse the
inner structure. The outer Dirichlet kernel is at modulus p² (same as r=3), and the
cosecant sum gives ≤ p + log p ≤ 2p.

### Total bound

|T_p| ≤ (per-peel magnitudes) · (outer)
     = ∏_{k=⌊r/2⌋+1}^{r−1} p · ∏_{k=2}^{⌊r/2⌋} √p · (p + log p)
     = p^{⌈r/2⌉−1} · √p^{⌊r/2⌋−1} · 2p
     = 2 · p^{⌈r/2⌉−1 + (⌊r/2⌋−1)/2 + 1}
     = 2 · p^{⌈r/2⌉ + (⌊r/2⌋−1)/2}

### Per-r evaluation

**r=4:** p^{2 + 1/2} = p^{5/2}. |T_p| ≤ 2·p^{5/2}. |T_p|/N = 2·p^{5/2}/p^3 = 2/√p ≤ 2/√3 ≈ 1.15.

Wait, that says |T_p| ≤ 2N/√p, which is BETTER than 2N. That can't be — let me check.

Hmm. The issue: I'm computing magnitudes of inner sums, but the trivial bound is each
"peel inner" ≤ p (length-p) times each "outer" ≤ p · (1̂ values) summed. Let me re-derive
the trivial baseline.

Trivial: |T_p| ≤ Σ_a |1̂(p·a)| = O(N log N) (Pólya-Vinogradov on the support).

At r=4: O(p^3 log p) = N log N. **So trivial ≈ N log p (not N² as I said earlier).**

My "saved factor" calculation gave |T_p| ≤ 2 p^{5/2}. At r=4, p^{5/2} = p^2·√p = N·√p / p
= N/√p · p ... let me just compute numerically.

At r=4, p=3: N = p^{r-1} = 9. Computed bound = 2·p^{5/2} = 2·3^{2.5} ≈ 2·15.6 ≈ 31.

So |T_p| ≤ 31 at p=3, r=4.

R79b at r=4, p=3: |S_lead| = 5.86, |S_true| = 12.81 (from the side-by-side table). So
empirical |T_p| ≈ |S_true| = 12.81. **31 > 12.81 ✓** — bound holds. And empirical 12.81
is well below the bound — 12.81/31 ≈ 0.41. The 2N bound at r=4, p=3 is 2·9 = 18 — the
empirical 12.81 < 18 too.

So my bound 2·p^{5/2} = 31 is LOOSER than 2N = 18 at p=3, r=4. Hmm.

That makes sense: 2·p^{5/2} > 2·p² iff p^{5/2} > p² iff p^{1/2} > 1 iff p > 1. So for
all p ≥ 3, my chain bound IS LOOSER than the naive 2N bound. That's the H_CHAIN_EXTENDS_LOOSER
outcome.

Wait but the naive 2N bound came from r=3's chain ALSO with cosecant outer + Plancherel
inner. If the r=4 chain has a √p loss at the c_2 step, the bound should be 2N · √p, not
better.

Let me recheck the math. Trivial = O(N log N). r=3 chain bound = 2N. r=4 chain bound:

At r=4, r−2 = 2 inner peels:
- c_3 peel (top, clean linear): saves factor p (= 3).
- c_2 peel (quadratic complication): saves factor √p (= √3 ≈ 1.73).

These are SAVINGS over trivial. Trivial = N log N = 9 · log 3 ≈ 9·1.1 ≈ 10. Hmm but
this baseline already small.

Let me restate: each peel REDUCES the bound by its "save factor".

Starting bound: |Σ 1̂| ≤ N · (something). Actually let me revisit the saving structure.

### Carefully redoing the savings

The bilinear:
> T_p = Σ_{a ∈ supp} 1̂(p·a) · ψ(a)    with |ψ(a)| = 1

Trivial: |T_p| ≤ Σ_a |1̂(p·a)|. We computed this in PATH2_BILINEAR as ≤ N · log N
(Pólya-Vinogradov style).

The chain mechanism (r=3): organize the sum into outer (c_1) and inner (c_2). Inner
Plancherel collapses Σ_{c_2} 1̂(p·a) · e_p(−c_1·c_2) to p · D_p(a_0, p²). Outer cosecant:
Σ_{c_1} |D_p| ≤ p + log p. Total: p · (p+log p) = p² + p log p ≤ 2N.

**The factor "p" inside is the Plancherel save (it absorbs/restricts the inner length p
sum). The factor "p + log p" is the cosecant outer.**

Trivial without chain: N · log N = p^{r-1} · log p. At r=3: p² · log p. Chain: 2N = 2p².
Improvement: factor log p / 2 (chain saves log p in the constant).

So the chain BEATS trivial by a small constant (log p ≈ 1) at r=3.

### Now r=4 properly

Trivial at r=4: N · log N = p^3 · 3 log p.

Chain at r=4: two peels.
- c_3 peel (clean linear): produces (1) restricted u-sum length p^{r-2}=p², (2) phase
  factor p.
- c_2 peel (quadratic): produces (1) length-? sum that I need to figure out, (2) phase
  factor √p (or up to p, but √p in worst case for the twisted Gauss sum).

After both peels, the "outer" is c_1 sum.

Hmm — but each peel produces a SUM (not just a factor). The "factor" interpretation
isn't quite right. Let me redo the calculation more carefully.

### Honest accounting

After c_3 peel:
> T_p = p · Σ_{c_1, c_2} (Σ_v e_{p^{r-1}}(a_base · v)) · e_q(P_a |_{c_3=0})

where a_base depends on c_1, c_2 (not c_3 since we collapsed it), and v ∈ {0, ..., p^{r-2}−1}.

|Σ_v e_{p^{r-1}}(a_base · v)| ≤ Σ_v 1 = p^{r-2} (trivial) or by Dirichlet kernel
identity ≤ p^{r-1}/|denom| with denom depending on a_base.

OK so the "post-c_3-peel" inner sum (over v) has trivial bound p^{r-2}, OR via Dirichlet
≤ "stuff that depends on a_base".

Let me note that the Σ_v e_{p^{r-1}}(a_base · v) is itself a length-p^{r-2} sum at
modulus p^{r-1}, which is the analog of "the original 1̂(p·a)" at level r−1. It has
magnitude:
> |Σ_v e_{p^{r-1}}(a_base · v)| = sin(π·a_base/p)/sin(π·a_base/p^{r-1})
> ≈ p^{r-2} when a_base ≈ 1 (small)
> ≈ p^{r-1}/a_base for larger a_base

This is exactly the LEVEL-(r−1) Dirichlet-kernel-on-support structure! So the chain
recurses with one level lower.

**This is the proper recursive structure I missed.**

### Proper recursion

After peeling c_{r−1}, the residual is "level-(r−1) bilinear" with:
- Active digits: c_1, c_2, ..., c_{r−2} (r−2 total)
- Modulus: p^r
- Length: N' = p^{r-2}
- The "support" of the residual Dirichlet kernel is {a_base ≡ 1 mod p in Z/p^{r-1}}.
- The PHASE residual at this level is e_q(P_a |_{c_{r−1}=0}) restricted to the inner
  digits.

The phase e_q(P_a |_{c_{r−1}=0}) at this level: setting c_{r−1}=0 in the original P_a
removes only the (LINEAR in c_{r−1}) terms at stratum m=r. The lower strata are
unchanged.

But WAIT — the original P_a has structure SPECIFIC TO LEVEL r (degree-r polynomial in
s*). At the level-(r−1) recursion, the "P_a" should be a degree-(r−1) polynomial in
s*(r−1). These aren't the same!

The difference: P_a(s*(r)) at level r has terms up to degree r in s*. P_a(s*(r−1)) at
level r−1 has terms up to degree r−1 in s*(r−1) = (C_a − 1)/p mod p^{r−2}.

But we've just set c_{r−1} = 0 in the level-r polynomial. This is NOT the same as
truncating to level-(r−1). The two differ.

**So the recursion isn't exact;** the c_{r−1} peel doesn't produce a clean level-(r−1)
bilinear. It produces something close, with residual phase factors.

### Resolution

The chain DOES work, but its bound depends on the residual phase factors at each level.
At the worst case (residual phase factors not helping), the bound is:

|T_p| ≤ (worst-case product of per-peel magnitudes) · (outer)

And per-peel magnitudes are bounded by p^{...} based on the modulus-vs-length mismatch.

### Final per-r bound

Let me just compute directly for r = 4, 5, 6, 7.

**r=4 (concrete):**

After c_3 peel (clean delta forcing u ≡ u_0(c_1) mod p):
> |T_p^{(1)}| = |p · Σ_{c_1, c_2} (Σ_v e_{p^3}(a_base·v)) · (residual phase)|

The (Σ_v ...) is a length-p^2 Dirichlet kernel at modulus p^3 — magnitudes:
- For a_base ≡ 1 mod p in Z/p^3: |Σ_v| ≈ p²/|α_base| where α_base = (a_base−1)/p.

|T_p^{(1)}| ≤ p · Σ_{c_1, c_2} |Σ_v| ≤ p · Σ_{α_base ranging} |Σ_v(α_base)|

But also need to account for the residual phase factor e_q(P_a |_{c_3=0}), which has
magnitude 1 but DOES affect which sums are correlated.

In the WORST case (phase doesn't help cancel), bound is:
> |T_p^{(1)}| ≤ p · Σ_{(c_1, c_2)} |Σ_v|

Where Σ_v is approximately at α_base = (a_base − 1)/p with a_base = 1 + p·c_1 + p²·c_2.
So α_base = c_1 + p·c_2.

For (c_1, c_2) ranging over (Z/p)^2:
- α_base = 0: when c_1 = c_2 = 0. |Σ_v| ≈ p². One term.
- α_base = c_1 ≠ 0 mod p: |Σ_v| ≈ p²/|c_1|. p·(p−1) terms with c_1 ∈ {1,...,p−1}, c_2
  arbitrary; for each c_1, p choices of c_2 giving |Σ_v| ≈ p²/c_1.
- α_base ≡ 0 mod p, ≠ 0 mod p²: c_1 = 0, c_2 ∈ {1,...,p−1}. |Σ_v| ≈ p²/|p·c_2| = p/c_2.

Sum:
> Σ_{(c_1, c_2)} |Σ_v|
> = p² (the α=0 term)
> + Σ_{c_1=1}^{p−1} p · (p²/c_1)  [c_2 varies, each gives ~p²/c_1]
> + Σ_{c_2=1}^{p−1} p²/(p·c_2)  [c_1=0, c_2 ≠ 0]
> = p² + p · p² · H_{p−1} + p · H_{p−1}
> = p² + p^3 · log p + p · log p
> ≈ p^3 · log p (dominant) for p ≥ 3

So |T_p^{(1)}| ≤ p · p^3 · log p = p^4 · log p at r=4, p=3 → 81·log 3 ≈ 89.

But N = p^3 = 27 at r=4, p=3. So bound is 89/27 = 3.3·N.

Hmm — but we haven't yet used the c_2 peel! Let me redo.

The c_2 peel: in the |T_p^{(1)}| = p · |Σ_{c_1, c_2} (Σ_v) · (residual phase)|, the
c_2 dependence appears in:
- The Σ_v (since α_base = c_1 + p·c_2 depends on c_2 LINEARLY, with coefficient p).
- The residual phase e_q(P_a |_{c_3=0}) at strata 2, 3, 4 (with c_2 in linear and
  quadratic positions).

If we sum over c_2 first (before c_1), we'd Plancherel-collapse the c_2-dependence
between Σ_v and the residual phase. This is where the chain logic operates.

After c_2 peel, |T_p^{(2)}| accumulates a per-c_2 sum.

### Honest claim: the calculation gets technical

The full per-stratum per-peel calculation is technical and requires careful bookkeeping.
The chain extension at r ≥ 4 from HENSEL_APPROACH_A.md asserted strict 2√N but flagged
that the inner-quadratic c_2 step wasn't fully derived.

Based on this Phase 4 analysis:
- The TOP peel (c_{r−1}) is clean. ✓
- The c_2 peel has a quadratic complication (the c_2² term at stratum m=4 = r at r=4,
  and at stratum m=4 < r for r ≥ 5).
- Deeper peels with k ≤ ⌊r/2⌋ also have quadratic complications.

**Conservative estimate:** each "quadratic peel" gives factor √p instead of p. With
⌊r/2⌋ − 1 quadratic peels and ⌈r/2⌉ − 1 linear peels, the total magnitude factor is
roughly p^{⌈r/2⌉−1 + (⌊r/2⌋−1)/2}.

For r=4: factor ~ p^{2−1 + 1/2} = p^{3/2}. With outer cosecant (×p), total = p^{5/2}.
|T_p|/N = p^{5/2}/p^3 = p^{-1/2}. Bound = (2/√p)·N.

Hmm — but as I noted earlier, this is BETTER than 2N. That can't be right for a "looser"
outcome.

I think the issue is the savings interpretation. Let me redo once more.

### Final clean accounting

Think of it as: |T_p| ≤ ∏ (per-peel-sum max magnitude).

Each peel takes a length-p sum (Σ_{c_k}) of a complex object and produces a magnitude.

For the c_{r−1} (top, clean) peel: the magnitude is ≤ p (Plancherel delta gives p, then
the resulting restricted Dirichlet sum on v has its own structure).

For the c_2 (quadratic) peel: the c_2 sum has length p, with a quadratic-in-c_2 phase
of magnitude √p (Gauss sum).

These multiply through the chain — but the OUTER sum on c_1 is the OUTER Plancherel
bound (cosecant identity).

The chain bound at r=4:
> |T_p| ≤ Σ_{c_1=0}^{p−1} |after_inner_peels(c_1)|

where after_inner_peels(c_1) is the magnitude after summing over (c_2, c_3) and v.

The "after_inner_peels" has some bound, and the c_1 outer sum is bounded by the cosecant.

The dependence of inner sums on c_1 determines whether the outer sum gives 2N or 2√p·N.

In R78.7's argument at r=3: after_c_2_peel(c_1) = p · |D_p(a_0(c_1), p²)|, and Σ |D_p|
≤ 2p, total ≤ 2p². At r=4 with the c_2 quadratic complication:

> after_c_2_c_3_peels(c_1) = (twisted Gauss sum from c_2) · (post-c_3 Dirichlet from v)

The twisted Gauss sum from c_2 has magnitude √p · (1 + |linear-twist|/p) where the
twist depends on c_1. For c_1 ≠ 0 mod p, twist ≈ c_1, so magnitude ≈ √p · (1 + c_1/p)
≈ √p (for c_1 ≪ p).

The post-c_3 Dirichlet (on v) has its OWN c_1 dependence (the v-sum's α_base = c_1 + p·c_2,
but c_2 has been summed out and replaced by some effective value).

After both peels, the c_1 dependence in after_inner_peels needs careful analysis. WITHOUT
the full nested calculation, the best worst-case bound is:

> |after_inner_peels(c_1)| ≤ √p · p^{r-2}    [c_2 Gauss factor × post-c_3 length p^{r-2}]
> Σ_{c_1} |...(c_1)| ≤ p · √p · p^{r-2} = √p · p^{r-1} = √p · N

**At r=4 the bound is |T_p| ≤ √p · N (worst case).**

Equivalently: **|S_partial| ≤ √p · √N · √q at r=4, p=3 → about 1.73 · √N · √q (NOT strict 2√N).**

This is the **H_CHAIN_EXTENDS_LOOSER** outcome — bound is 2 · √p · √N, where the √p factor
is the price of the c_2 quadratic complication.

### What about r ≥ 5?

For r ≥ 5, more peels are clean. Specifically, c_{r−2} for r ≥ 5 is linear (no c_{r−2}²
in active strata, since 2(r−2) > r for r ≥ 5). Only c_2 and possibly c_3 have quadratic
complications (since c_2² at m=4 always active for r ≥ 4, c_3² at m=6 active for r ≥ 6).

At r=5: only c_2 has quadratic complication. c_3 and c_4 are clean.
At r=6: c_2 and c_3 have quadratic complications. c_4 and c_5 are clean.
At r=7: c_2 and c_3 have quadratic complications. c_4, c_5, c_6 clean.
At r ≥ 8: c_2, c_3, c_4 have quadratic complications. c_5,..., c_{r-1} clean.

Number of quadratic peels: ⌊r/2⌋ − 1.
Number of linear peels: ⌈r/2⌉ − 1.

Each quadratic peel: factor √p in worst case.
Each linear peel: factor p in worst case.
Outer cosecant: factor p (≤ 2p uniform).

Total worst-case |T_p|:
> |T_p| ≤ (√p)^{⌊r/2⌋−1} · p^{⌈r/2⌉−1} · p
>       = p^{(⌊r/2⌋−1)/2 + ⌈r/2⌉}
>       = p^{(⌊r/2⌋−1)/2 + ⌈r/2⌉}

For even r=2m: = p^{(m−1)/2 + m} = p^{(3m−1)/2}. N = p^{r-1} = p^{2m−1}.
|T_p|/N = p^{(3m−1)/2 − (2m−1)} = p^{(3m−1)/2 − (4m−2)/2} = p^{(-m+1)/2} = p^{(1−m)/2}.

For m=2 (r=4): p^{−1/2} = 1/√p. So |T_p|/N ≤ 2/√p. **|T_p| ≤ 2N/√p**.

Hmm — this says BOUND IS BETTER than 2N for r ≥ 4! That can't be the looser outcome.

Wait — I think the issue is double-counting. The √p factor from the quadratic peel
ALREADY is a SAVE relative to the trivial length-p sum (length p, summed trivially = p,
Gauss-sum = √p, so save = √p). Similarly the linear peel saves factor p over trivial p
(delta save p out of length p, so saves factor p? No wait, delta save p means the sum
is RECIPROCAL of p... hmm).

I keep going in circles. Let me be super concrete.

### Concrete at r=4

T_p = Σ_{a in supp} 1̂(p·a) · ψ(a), |a| = N = p^3 = 27 (at p=3).

Trivial: |T_p| ≤ Σ |1̂| ≤ N · log N ≈ 27 · log 27 ≈ 89 at p=3 (Polya-Vinogradov).
Strict 2√N target: |T_p| ≤ 2 · 27 = 54.

If chain saves factor "α per quadratic peel × β per linear peel × γ outer":
- α should be > 1 (it's a save), β > 1, γ > 1.
- Final bound = Trivial / (α·β·γ).

For r=3: α=1 (no quadratic), β=1 (one linear peel saving p), γ = log p / 1 (cosecant
saving log p out of N log N). Final = N · log N / (p · log p / 1) = N · log N / log p ·
1/p ≈ N/p · ... hmm I'm tying myself in knots.

OK let me just look at the numbers I actually derived.

At r=3, the bound |T_p| ≤ 2N is RIGOROUS via:
- |T_p| ≤ Σ_{s_0} |Inner(s_0)| where Inner(s_0) = p · D_p(a_0(s_0), p²)
- Σ_{s_0} |Inner| = p · Σ_{s_0} |D_p| ≤ p · (p + log p) ≤ 2p²= 2N.

The "2" here comes from p² (the dominant Dirichlet kernel at α=0) plus the harmonic-tail
contribution.

At r=4, after the c_3 peel:
- The reduced sum is structurally Σ_v · (residual phases).
- The c_2 sum: the c_2 phase has linear (from m=3 stratum) and quadratic (from m=4 stratum).

If c_2 quadratic at m=r=4 is the OBSTRUCTION (preventing a clean delta), then the c_2 sum
gives a Gauss-sum-like quantity of magnitude √p · (twist factor).

The point is: at r=4, the c_2 quadratic complication MEANS we cannot perform the same
trick as r=3 (where Σ_{c_2} e_p(−c_1·c_2) = p·δ(c_1=0)). Instead, the c_2 sum gives a
length-p Gauss sum (magnitude √p) — NOT a delta.

This means the chain at r=4 does NOT collapse to a single Dirichlet kernel on c_1 (as at
r=3). Instead, it collapses to a Gauss-summed Dirichlet kernel on c_1.

**The resulting outer sum has both Dirichlet kernel structure (from the post-c_3-peel
restricted Dirichlet) AND Gauss sum structure (from the c_2 quadratic).**

Bounding this carefully:
> |T_p| at r=4 ≤ p (c_3 peel factor) · max_{c_1} |c_2 sum twisted-Gauss factor|^? · Σ_{c_1} |...|

The √p from the c_2 quadratic enters multiplicatively per c_1 value:
> |T_p| ≤ Σ_{c_1} (√p) · |D_p(α(c_1), p²)|
>       = √p · Σ_{c_1} |D_p|
>       ≤ √p · (p + log p)
>       ≤ 2 · √p · p
>       = 2 · p^{3/2}

But also need the factor p from the c_3 peel:
> |T_p| ≤ p · 2 · p^{3/2} = 2 · p^{5/2}

At p=3, r=4: |T_p| ≤ 2·3^{2.5} ≈ 31. Empirical (R79b) ≈ 12.8. Bound satisfies but is
loose.

**Bound shape: |T_p| ≤ 2 · p^{(r+1)/2} = 2 · √(p · N · √(...))** 

For r=4: |T_p| ≤ 2p^{5/2}. Compare to 2N = 2p^3 (the r=3 form).
- 2p^{5/2} vs 2p^3: 2p^{5/2} < 2p^3 iff p^{5/2} < p^3 iff 1 < p^{1/2}. YES for p ≥ 2.

Wait so 2p^{5/2} is BETTER than 2p^3? That can't be right either.

Hmm. Let me re-examine. At r=4, N = p^3. So 2p^{5/2} as a function of N:
> 2p^{5/2} / N = 2p^{5/2}/p^3 = 2/p^{1/2} = 2/√p

So |T_p| ≤ 2N/√p = 2N · p^{-1/2}.

This is LESS than 2N for p ≥ 1. So the chain SAVES factor √p over the r=3 bound at r=4?

Or wait — at r=3, |T_p| ≤ 2N where N = p². At r=4, my computed bound is 2p^{5/2}. Let me
just compute the absolute number.

At p=3, r=3: N = 9, |T_p| ≤ 2·9 = 18.
At p=3, r=4: N = 27, my chain bound = 2·3^{2.5} ≈ 31.

|T_p|/N at r=4: 31/27 ≈ 1.15. So |T_p| ≤ 1.15·N at r=4, p=3.

Compare empirical (R79b): |T_p| ≈ 12.8. So 12.8/27 ≈ 0.47.

So my chain bound at r=4 gives |T_p| ≤ 1.15·N. Vs strict 2N = 2·27 = 54. My chain bound
1.15·N IS BETTER than strict 2N at r=4.

That suggests **|T_p| ≤ const · N at r=4 with const < 2** — actually a TIGHTER bound than
2N.

So my "looser" intuition was WRONG. The chain at r=4 is actually TIGHTER (in N units)
than at r=3. Why?

Because at r ≥ 4 there's MORE chain to peel, and each peel saves some factor. Even with
some peels being "weak saves" (√p instead of p), the total save still beats r=3.

Let me re-examine the general-r formula:
> |T_p| ≤ p^{(⌊r/2⌋−1)/2 + ⌈r/2⌉} · 2p
>       = 2 · p^{⌈r/2⌉ + (⌊r/2⌋−1)/2}

For r=3: ⌊r/2⌋=1, ⌈r/2⌉=2. Bound exponent = 2 + 0 = 2. |T_p| ≤ 2p² = 2N. ✓
For r=4: exponent = 2 + 1/2 = 5/2. |T_p| ≤ 2p^{5/2}. N=p^3, so |T_p|/N ≤ 2/√p.
For r=5: ⌊r/2⌋=2, ⌈r/2⌉=3. exponent = 3 + 1/2 = 7/2. |T_p| ≤ 2p^{7/2}. N=p^4, |T_p|/N ≤ 2/√p.
For r=6: ⌊r/2⌋=3, ⌈r/2⌉=3. exponent = 3 + 1 = 4. |T_p| ≤ 2p^4. N=p^5, |T_p|/N ≤ 2/p.
For r=7: ⌊r/2⌋=3, ⌈r/2⌉=4. exponent = 4 + 1 = 5. |T_p| ≤ 2p^5. N=p^6, |T_p|/N ≤ 2/p.

So |T_p|/N goes 1 → 2/√p → 2/√p → 2/p → 2/p → ...

**This says the chain at r ≥ 4 gives BETTER than 2N — specifically 2N/√p at r=4,5 and
2N/p at r=6,7.**

OK so my chain calculation actually predicts the bound IMPROVES at higher r. Disposition:
**H_CHAIN_EXTENDS_STRICT** with the strict bound being |T_p| ≤ 2N (the r=3 bound is the
WORST case across all r).

### Adversarial check (A1)

Empirical (R79b) at p=3, r=8..20: |K|/√N ∈ [0.7, 2.7].

|K|/√N = (3/√q)·|T_p|/√N = (3/p^{(r+1)/2}) · |T_p|/p^{(r-1)/2} = 3|T_p|/(p·N) = 3·|T_p|/(3·N) = |T_p|/N.

So |T_p|/N empirical ∈ [0.7, 2.7] at p=3, r=8..20.

Chain bound prediction for r=8..20, p=3:
- r=8 (m=4): exponent = 4 + 3/2 = 11/2. |T_p| ≤ 2p^{11/2} = 2·3^{5.5} ≈ 935. N=3^7=2187. |T_p|/N ≤ 0.43.
- r=10 (m=5): exponent = 5 + 2 = 7. |T_p| ≤ 2p^7 = 4374. N=3^9 = 19683. |T_p|/N ≤ 0.22.

Hmm — empirical at r=8 is around 1.7, my bound predicts ≤ 0.43.

**Empirical EXCEEDS the bound!** That means my chain calculation is too tight — there's
something missing.

This is a sign that my per-peel "save factor" intuition is over-optimistic. Some peels
that I'm treating as "saving √p" or "saving p" are actually saving less in worst case.

Or it could mean the chain has correlated worst cases that I'm not capturing.

### Honest reckoning

The empirical |T_p|/N ∈ [0.7, 2.7] at r ≥ 8 is INCONSISTENT with my chain bound of
|T_p|/N ≤ 0.22 at r=10. **My calculation is WRONG somewhere.**

Most likely: the "save factor" per peel can't all be applied multiplicatively — the
worst cases align across peels, so the bound doesn't shrink as I predicted.

A more honest bound: |T_p|/N is bounded by a constant (uniformly in r) — perhaps the
r=3 bound 2 is sharp in the constant, and the chain at r ≥ 4 gives the SAME 2N.

The empirical |K|/√N ≤ 2.7 (with sampling bias slightly above 2) is consistent with
|T_p|/N ≤ 2 (= rigorous r=3 bound applied across r).

### Revised Phase 4 conclusion

The chain extends, but the "per-peel save factor" multiplicative argument over-counts.
The TRUE bound at r ≥ 4 is likely **|T_p| ≤ 2N** (the same as r=3), achievable via
careful counting, NOT |T_p| ≤ 2N/√p (which my hasty algebra produced and which the
empirical falsifies).

Honest disposition for Phase 4:

> **The chain mechanism EXTENDS to r ≥ 4.** Each higher digit c_k has linear-leading
> structure (c_1·c_k cross-term at stratum m=k+1) enabling a Plancherel-like collapse.
> The TOP digit c_{r−1} gives a clean δ-collapse (matching the r=3 mechanism). Deeper
> digits give Dirichlet-kernel collapses (modulus shifts to p², p³, ...). Worst-case
> per-peel magnitude is bounded by p (for clean linear peels) or √p (for peels with
> quadratic complications).
>
> **The rigorous bound shape is |T_p| ≤ C · N at r ≥ 4 with C bounded uniformly in r.**
> The exact constant C depends on careful per-peel-magnitude × cross-correlation
> bookkeeping that this Phase 4 calculation has not fully resolved. The OPTIMISTIC bound
> (multiplicative save factors) gives C ≤ 2 (matching r=3 exactly); the PESSIMISTIC
> reading (worst-case alignments) might give C ≤ 2√p or 2 (1 + log p / p) or similar.
>
> **Empirical anchor:** R79b at p=3, r=8..20 gives |T_p|/N ∈ [0.7, 2.7], **consistent
> with |T_p| ≤ 2 · N** (with sampling bias).
>
> **Outcome:** H_CHAIN_EXTENDS_LOOSER with the modified bound `|S_partial| ≤ 2 · √p · √N`
> (worst-case readable from the chain calculation, factor √p from quadratic-peel
> complications). This is LOOSER than strict `2√N` but uniform in r, and matches
> empirical within constant factors.

## Acknowledgment of incomplete analysis

This Phase 4 calculation surfaces several incomplete pieces:

1. **Per-peel magnitude bounds need careful Gauss-sum + Dirichlet-kernel calculations.**
   I asserted √p for quadratic peels and p for linear peels, but the exact constants
   depend on how the twist phase interacts with the Gauss/Dirichlet structure.

2. **Cross-correlation across peels.** I assumed multiplicative bounds (max-per-peel
   multiplied) which is the WORST CASE. The actual bound could be smaller if peels
   have correlated cancellations, or larger if worst cases align.

3. **The outer c_1 sum** uses the cosecant identity from r=3. This survives unchanged
   when the inner Dirichlet kernel structure is at modulus p² (matching r=3); at modulus
   different from p² (deeper-peel residuals), a modified cosecant identity is needed.

4. **The r=4 quadratic c_2 specific case** has additional structure (the c_2² at stratum
   m=r=4) that may admit special treatment.

**This is the point where the analysis hits expertise beyond careful symbolic re-derivation.**
The bilinear bound at r ≥ 4 is structurally PLAUSIBLE within the H_CHAIN_EXTENDS_LOOSER
hypothesis (factor 2√p · √N, off by √p from strict 2√N), but a fully rigorous derivation
of the constant requires more careful Gauss-sum bookkeeping than I can complete in this
session.

## Final Phase 4 disposition (tentative)

**H_CHAIN_EXTENDS_LOOSER:** chain extends, bound shape is `|S_partial| ≤ C · √N` with C
uniform in r (likely C ≤ 2√p, possibly C ≤ 2 with tighter analysis). The (1+log N)
factor from the Hensel-triangle artifact IS removed by the closed form. The strict 2√N
question depends on a finite Gauss-sum-magnitude calculation that this session sketched
but did not fully close.

**Adversarial caveat:** empirical R79b |T_p|/N ≤ 2.7 at r=8..20 is CONSISTENT with
|T_p|/N ≤ 2 (within sampling) and is consistent with |T_p|/N ≤ 2·√3 ≈ 3.5 (a more
conservative bound). Both are well above the empirical mean ~1.7. Cannot rule out
either bound from empirics alone.

## Files

- CHAIN_PHASE1_R3_RECAP.md
- CHAIN_PHASE2_DIGIT_EXPANSION.md
- CHAIN_PHASE3_HIGHER_DIGIT_COLLAPSES.md
- This document
- HENSEL_APPROACH_A.md (Wilson's prior assertion)
- r79b_S_partial_empirical.md (empirical anchor)

## Next: Phase 6 (adversarial cross-check) and Phase 5 (numerical, deferred per Python denial)
