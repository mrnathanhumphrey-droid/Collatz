# CHAIN Phase 2 — Digit Expansion of P_a(s*(r)) at r ≥ 4

**Date:** 2026-05-11. Analyst: Wilson. Reports to: Nathan.

**Goal:** Expand the Hensel-lifted phase polynomial P_a(s*(r)) by p-power stratum in the
inner digits c_2, c_3, ..., c_{r−1}, exposing which strata contain which digits and at
what polynomial degree. This is the data Phase 3 tests for Fourier-collapse.

## Notation

- q = p^{r+1}, period = p^r, support size N = p^{r-1}
- C_a ∈ Z/p^r digit decomposition: **C_a = 1 + p·c_1 + p²·c_2 + p³·c_3 + ... + p^{r−1}·c_{r−1} mod p^r**
- s*(r) := (C_a − 1)/p mod p^{r−1} = c_1 + p·c_2 + p²·c_3 + ... + p^{r−2}·c_{r−1}
- s_0 := c_1 (outer / saddle / Plancherel-outer digit). Use s_0 and c_1 interchangeably.
- The "inner digits" we Fourier-collapse over are c_2, c_3, ..., c_{r−1} (count: r−2).

## Closed-form polynomial (Hensel)

For p > r (clean range, no small-prime denominator merges):

> **P_a(s*(r)) ≡ Σ_{j=2}^{r} (−1)^{j−1} · (p·s*)^j / (j(j−1)) mod p^{r+1}**

Coefficients: c_j := (−1)^{j−1} / (j(j−1)) for j ≥ 2.
- c_2 = −1/2,  c_3 = +1/6,  c_4 = −1/12,  c_5 = +1/20,  c_6 = −1/30, c_7 = +1/42.

Expanding p·s* = p·c_1 + p²·c_2 + p³·c_3 + ... + p^{r−1}·c_{r−1}, we have

> (p·s*)^j = Σ_{(i_1,...,i_j) : 1 ≤ i_k ≤ r−1, ordered}  p^{i_1+...+i_j} · c_{i_1}·...·c_{i_j}

Total p-power of a monomial c_{i_1}·...·c_{i_j} is **i_1 + i_2 + ... + i_j** (the sum of
digit-positions).

The full P_a phase is the sum over j ∈ {2,...,r} and over all index multisets. We
collect contributions by **p-power stratum m** mod p^{r+1}:

> P_a(s*(r)) = Σ_{m=2}^{r+1} p^m · Q_m(c_1, c_2, ..., c_{r−1}) mod p^{r+1}

where Q_m is a polynomial in the digits whose coefficients are determined by the
generating identity.

(Note: m starts at 2 because the j=2 term contributes at minimum p^2 (i_1=i_2=1 case).
m ≤ r+1 is the modulus boundary.)

## Stratum-by-stratum: total p-power m = i_1+...+i_j

For each j ≥ 2 and each multiset {i_1, ..., i_j} of digit-positions, the contribution is
at stratum m = Σ i_k.

The monomial c_{i_1}·c_{i_2}·...·c_{i_j} (with all i_k ≥ 1) has total degree j and total
weight m. Constraints: m ≥ j (each i_k ≥ 1), m ≤ j·(r−1) (each i_k ≤ r−1).

**Stratum m has contributions from all j with j ≤ m ≤ j·(r−1), i.e., j ≥ m/(r−1) and
j ≤ m.**

Group by the digit-index sequence (compositions of m into j parts, each part ≥ 1, each
part ≤ r−1):

> Q_m = Σ_{j=2}^{min(m, r)} c_j · M_{j,m}

where M_{j,m} := Σ_{(i_1,...,i_j) composition of m, parts in [1,r−1]} c_{i_1}·...·c_{i_j}
(symmetric in the indices since ordering doesn't matter for the product; expand
multinomially).

Equivalently, M_{j,m} = j!/(α_1!·α_2!·...) · (monomial in c_i^{α_i} with Σ α_i = j, Σ
i·α_i = m), summed over valid (α_1, α_2, ...).

## Per-r explicit strata

### r = 4 (q = p^5)

Digits available: c_1, c_2, c_3. p-power stratum m runs over 2, 3, 4, 5 (mod p^5).

**Stratum m=2:** j=2, only composition is (1, 1). Q_2 = c_2 · c_1² = (−1/2) · c_1².
> Q_2 = −c_1² / 2

**Stratum m=3:** 
- j=2, composition (1,2) and (2,1): M_{2,3} = 2·c_1·c_2. Contribution: c_2 · 2·c_1·c_2 = (−1/2) · 2 c_1·c_2 = −c_1·c_2.

  Wait — c_2 is doing double duty as the coefficient (−1/2) and as the digit. Let me
  rename coefficients to avoid collision.

  **Renaming:** B_j := (−1)^{j−1}/(j(j−1)) are the polynomial coefficients (call them
  Bernoulli-like). The digit name remains c_2, c_3, etc.

  B_2 = −1/2, B_3 = +1/6, B_4 = −1/12, B_5 = +1/20, B_6 = −1/30, B_7 = +1/42.

- j=2, composition (1,2) and (2,1): M_{2,3} = 2·c_1·c_2. Contribution to Q_3: B_2 · 2·c_1·c_2 = −c_1·c_2.
- j=3, composition (1,1,1): M_{3,3} = c_1³. Contribution: B_3 · c_1³ = c_1³/6.
> **Q_3 = −c_1·c_2 + c_1³ / 6**

**Stratum m=4:**
- j=2, compositions (1,3), (3,1), (2,2): M_{2,4} = 2·c_1·c_3 + c_2². Contribution: B_2·(2·c_1·c_3 + c_2²) = −c_1·c_3 − c_2²/2.
- j=3, compositions: (1,1,2), (1,2,1), (2,1,1): M_{3,4} = 3·c_1²·c_2. Contribution: B_3 · 3·c_1²·c_2 = c_1²·c_2/2.
- j=4, composition (1,1,1,1): M_{4,4} = c_1^4. Contribution: B_4 · c_1^4 = −c_1^4 / 12.
> **Q_4 = −c_1·c_3 − c_2²/2 + c_1²·c_2/2 − c_1^4/12**

**Stratum m=5:**
- j=2: compositions (2,3), (3,2): but also wait (1,4)? max index is r−1=3, so (1,4) forbidden. Only (2,3) and (3,2): M_{2,5} = 2·c_2·c_3.
  Contribution: −c_2·c_3.
- j=3: compositions of 5 into 3 parts each ≥1 ≤3: (1,1,3) perms = 3; (1,2,2) perms = 3; (1,3,1), (3,1,1) absorbed in perms; (2,1,2), (2,2,1) absorbed. So M_{3,5} = 3·c_1²·c_3 + 3·c_1·c_2².
  Contribution: B_3 · (3·c_1²·c_3 + 3·c_1·c_2²) = c_1²·c_3/2 + c_1·c_2²/2.
- j=4: compositions of 5 into 4 parts each in [1,3]: (1,1,1,2) perms = 4. M_{4,5} = 4·c_1³·c_2.
  Contribution: B_4 · 4·c_1³·c_2 = −c_1³·c_2/3.
- j=5: (1,1,1,1,1) only. M_{5,5} = c_1^5. Contribution: B_5 · c_1^5 = c_1^5 / 20.
> **Q_5 = −c_2·c_3 + c_1²·c_3/2 + c_1·c_2²/2 − c_1³·c_2/3 + c_1^5/20**

### Summary table at r=4

| stratum m | polynomial in (c_1, c_2, c_3) | degree in c_3 (top digit) | degree in c_2 |
|---|---|---|---|
| 2 | −c_1²/2 | 0 | 0 |
| 3 | −c_1·c_2 + c_1³/6 | 0 | 1 |
| 4 | −c_1·c_3 − c_2²/2 + c_1²·c_2/2 − c_1^4/12 | **1** | 2 |
| 5 | −c_2·c_3 + c_1²·c_3/2 + c_1·c_2²/2 − c_1³·c_2/3 + c_1^5/20 | **1** | 2 |

(Note stratum 5 = r+1 is the modulus boundary; mod p^{r+1} this stratum contributes via
e_{p^{r+1−m}} = e_1 = trivial. Stratum 5 thus DOES NOT CONTRIBUTE to e_q(P_a). The
relevant strata are m = 2, 3, 4.)

**Active strata at r=4: m ∈ {2, 3, 4}.**

### r = 5 (q = p^6)

Digits c_1, c_2, c_3, c_4. Active strata m ∈ {2, 3, 4, 5} (m=6 mod p^6 trivial).

Q_2 = −c_1²/2 (same as r=4)
Q_3 = −c_1·c_2 + c_1³/6 (same as r=4, since c_3, c_4 don't reach m=3 with j=2)
Q_4 = −c_1·c_3 − c_2²/2 + c_1²·c_2/2 − c_1^4/12 (same as r=4)

Q_5 (now with c_4 available at index 4):
- j=2 compositions of 5 with parts in [1,4]: (1,4), (4,1), (2,3), (3,2). M_{2,5} = 2·c_1·c_4 + 2·c_2·c_3.
  Contribution: −c_1·c_4 − c_2·c_3.
- j=3 compositions of 5 in [1,4]: (1,1,3) → 3; (1,2,2) → 3; (1,3,1), (3,1,1) absorbed. Plus (1,1,3) already and now (1,1,3) is fine because max=4. M_{3,5} = 3·c_1²·c_3 + 3·c_1·c_2².
  Contribution: c_1²·c_3/2 + c_1·c_2²/2.
- j=4: (1,1,1,2) perms = 4. M_{4,5} = 4·c_1³·c_2. Contribution: −c_1³·c_2/3.
- j=5: (1,1,1,1,1). M_{5,5} = c_1^5. Contribution: c_1^5/20.
> **Q_5 (r=5) = −c_1·c_4 − c_2·c_3 + c_1²·c_3/2 + c_1·c_2²/2 − c_1³·c_2/3 + c_1^5/20**

### Summary table at r=5

| stratum m | polynomial | top-digit (c_4) degree | c_3 degree | c_2 degree |
|---|---|---|---|---|
| 2 | −c_1²/2 | 0 | 0 | 0 |
| 3 | −c_1·c_2 + c_1³/6 | 0 | 0 | 1 |
| 4 | −c_1·c_3 − c_2²/2 + c_1²·c_2/2 − c_1^4/12 | 0 | **1** | 2 |
| 5 | −c_1·c_4 − c_2·c_3 + c_1²·c_3/2 + c_1·c_2²/2 − c_1³·c_2/3 + c_1^5/20 | **1** | 1 | 2 |

### r = 6 (q = p^7)

Digits c_1, c_2, c_3, c_4, c_5. Active strata m ∈ {2, 3, 4, 5, 6}.

Q_2, Q_3, Q_4 same as r=4, 5.

Q_5 with c_5 available (index 5 now allowed): but index 5 > r−1 = 5 wait r=6 so r−1 = 5, so c_5 IS available.

Actually wait — for r=6, s*(r) = (C_a−1)/p mod p^{r−1} has digits c_1, ..., c_{r-1} = c_1, ..., c_5. So c_5 is the top digit.

Q_5 (r=6):
- j=2 compositions of 5 in [1,5]: (1,4), (4,1), (2,3), (3,2). Same as r=5 (since c_5 not in m=5 with j=2 needing pair-sum=5). M_{2,5} = 2·c_1·c_4 + 2·c_2·c_3.

  Wait, can we have (5,0)? No, since each part ≥ 1. Can we have (5)? Only if j=1, but j ≥ 2.

  So M_{2,5} = 2·c_1·c_4 + 2·c_2·c_3 (unchanged).
- j=3, j=4, j=5 same as r=5.
> Q_5 (r=6) = Q_5 (r=5) ... unchanged at m=5 since c_5 first appears at m=5 only via the linear term (j=1) which isn't in our sum.

Hmm — let me reconsider. The j=1 term (linear in c_5 via p·c_5) is the "constant" linear
term of the original L_p expansion that gets cancelled by the p·s* in P_a = p·s − C_a·L.
So c_5 LINEAR (single appearance) at p^5 stratum is NOT in our polynomial.

But via j=2: c_5 would appear at m = i_1+i_2 with one of i_1, i_2 = 5. Smallest is
(1,5), giving m=6, contributing 2·c_1·c_5. So c_5 FIRST APPEARS at stratum m=6 in
P_a(s*(r=6)).

Q_6 (r=6, NEW stratum for r=6):
- j=2 compositions of 6 in [1,5]: (1,5), (5,1), (2,4), (4,2), (3,3). M_{2,6} = 2·c_1·c_5 + 2·c_2·c_4 + c_3².
  Contribution: B_2·M_{2,6} = −c_1·c_5 − c_2·c_4 − c_3²/2.
- j=3 compositions of 6 in [1,5]: (1,1,4)→3; (1,2,3)→6; (2,2,2)→1; (1,3,2) covered; (1,4,1) covered; (3,1,2) covered etc. M_{3,6} = 3·c_1²·c_4 + 6·c_1·c_2·c_3 + c_2³.
  Contribution: B_3·M_{3,6} = c_1²·c_4/2 + c_1·c_2·c_3 + c_2³/6.
- j=4 compositions of 6 in [1,5]: (1,1,1,3)→4; (1,1,2,2)→6. M_{4,6} = 4·c_1³·c_3 + 6·c_1²·c_2².
  Contribution: B_4·M_{4,6} = −c_1³·c_3/3 − c_1²·c_2²/2.
- j=5 compositions of 6 in [1,5]: (1,1,1,1,2)→5. M_{5,6} = 5·c_1^4·c_2.
  Contribution: B_5·M_{5,6} = c_1^4·c_2/4.
- j=6 compositions of 6 in [1,5]: (1,1,1,1,1,1). M_{6,6} = c_1^6.
  Contribution: B_6·c_1^6 = −c_1^6/30.

> **Q_6 (r=6) = −c_1·c_5 − c_2·c_4 − c_3²/2 + c_1²·c_4/2 + c_1·c_2·c_3 + c_2³/6 − c_1³·c_3/3 − c_1²·c_2²/2 + c_1^4·c_2/4 − c_1^6/30**

### Summary table at r=6

| stratum m | top-digit (c_5) degree | c_4 degree | c_3 degree | c_2 degree |
|---|---|---|---|---|
| 2 | 0 | 0 | 0 | 0 |
| 3 | 0 | 0 | 0 | 1 |
| 4 | 0 | 0 | 1 | 2 |
| 5 | 0 | **1** | 1 | 2 |
| 6 | **1** | 1 | **2** | 3 |

## Pattern recognition: where each digit c_k first appears

Empirically from the tables:

| digit | first stratum | degree at first appearance | structure of leading term |
|---|---|---|---|
| c_1 (=s_0) | m=2 | 2 (quadratic) | −c_1²/2 |
| c_2 | m=3 | 1 (linear) | −c_1·c_2 |
| c_3 | m=4 | 1 (linear) | −c_1·c_3 |
| c_4 | m=5 | 1 (linear) | −c_1·c_4 |
| c_5 | m=6 | 1 (linear) | −c_1·c_5 |
| c_k | m=k+1 | 1 (linear) | −c_1·c_k (from j=2 cross-term) |

**General pattern (PROVEN by counting compositions):**

> The digit c_k FIRST APPEARS at stratum m = k+1 (via j=2 cross-term `2·c_1·c_k·p^{k+1}`
> from (p·s*)²/2 expanded), where it appears **LINEARLY** with coefficient −c_1 mod p
> (since B_2 = −1/2 and the cross-term has multinomial coefficient 2 cancelling the 1/2).
>
> At strata m > k+1, c_k can also appear via higher-degree terms or other compositions.
> But at the FIRST appearance m=k+1, the term is **uniquely linear in c_k** (the
> composition (1, k) of m=k+1 into j=2 parts is the ONLY composition involving c_k at
> j=2 stratum m=k+1, and higher j requires m ≥ j ≥ 3 > j=2 minimal, but k ≥ 2 and
> stratum m=k+1 with j=3 requires three indices summing to k+1, each ≥ 1 — possible
> only if k+1 ≥ 3, i.e., k ≥ 2; e.g., k=2, m=3, j=3 composition (1,1,1) gives c_1³.
> That doesn't involve c_2. So c_2 appears at m=3 only via j=2 (1,2) composition).

**Refined statement:** c_k at stratum m=k+1 appears ONLY through compositions of m=k+1
where one index equals k and the rest equal 1. For j=2: (1, k). For j=3: (1, 1, k−1)
needs index k−1 not k, so doesn't involve c_k. Etc. So at stratum m=k+1, c_k appears
ONLY via j=2 composition (1, k), uniquely.

**The linear-in-c_k coefficient at stratum m=k+1 is −c_1.** Specifically, the j=2 term
M_{2,k+1} contributes 2·c_1·c_k (multinomial factor 2 for the two orderings of (1,k))
times B_2 = −1/2, giving −c_1·c_k.

## The Top Stratum of the Modulus

The phase modulus is p^{r+1}, so the topmost meaningful stratum is m = r (giving e_p
factor). Stratum m = r+1 wraps to e_1 = trivial.

At m = r:
- Via j=2 composition (1, r−1): term −c_1·c_{r−1} (LINEAR in the top digit c_{r−1}).
- Via j=2 compositions (i, r−i) for 1 < i ≤ r−1 with both i, r−i ≤ r−1: these introduce
  cross-terms of c_i·c_{r−i}.
- Via j=3, 4, ..., r: various polynomial terms in c_1, c_2, ..., c_{r−2}.

**So the top stratum m=r contains the linear-in-c_{r−1} term with coefficient −c_1, plus
"lower-digit-only" polynomial terms.**

## Per-stratum digit profile at general r

| stratum m | active digits | leading top-digit term | top-digit-degree |
|---|---|---|---|
| m=2 | c_1 | −c_1²/2 | (c_1 quadratic) |
| m=3 | c_1, c_2 | −c_1·c_2 | (c_2 linear) |
| m=4 | c_1, c_2, c_3 | −c_1·c_3 | (c_3 linear) |
| ... | ... | ... | ... |
| m=k+1 (k ≥ 1) | c_1, ..., c_k | −c_1·c_k | (c_k linear) |
| ... | ... | ... | ... |
| m=r | c_1, ..., c_{r−1} | −c_1·c_{r−1} | (c_{r−1} linear) |

## Key observation for Phase 3

**Every higher digit c_k (k = 2, ..., r−1) first appears LINEARLY at stratum m=k+1 with
coefficient −c_1.** This is exactly the structural feature (I2) from Phase 1's r=3 chain.

This linearity is a direct consequence of:
1. The j=2 term in the (1+y)·log(1+y) generating series has coefficient −1/2.
2. The cross-term 2·c_1·c_k in (c_1·p + c_2·p² + ... + c_k·p^k + ...)² appears at
   p^{k+1} stratum with multinomial coefficient 2.
3. These combine to give exactly −c_1·c_k at stratum m=k+1.

**The structure is UNIVERSAL across r and across k — the (1+y)·log(1+y) generating identity
is the load-bearing input.**

## Caveat: "leading top-digit term" is not the FULL story per stratum

At each stratum m=k+1, beyond the linear-in-c_k cross-term, other polynomial terms
involving (c_1, c_2, ..., c_{k−1}) also live at the same stratum. These "lower-digit-only"
terms are CONSTANTS with respect to c_k (don't depend on c_k). When we sum over c_k,
they factor out of the c_k Plancherel collapse — they live in the "outer phase" relative
to c_k.

So at each peeling step (summing over c_k), the inner-Plancherel sees a phase of the
form e_p(−c_1·c_k + lower-digit-only stuff). The c_k sum collapses by Σ_{c_k}
e_p(−c_1·c_k) = p·𝟙[c_1 ≡ 0 mod p] regardless of the c_k-free parts (which just become
an outer phase factor).

**This is the load-bearing feature: at each stratum m=k+1, the c_k variable appears
LINEARLY with coefficient −c_1 (mod p), and is decoupled at the c_k-level from all the
lower digits at this stratum.**

## Issue: but each c_k ALSO appears at HIGHER strata m > k+1

At m=k+2: c_k can appear via:
- j=2 composition (2, k): term 2·c_2·c_k (linear in c_k, coefficient depending on c_2).
- j=3 composition (1, 1, k): term 3·c_1²·c_k (linear in c_k, coefficient depending on c_1²).

So c_k appears at strata m = k+1, k+2, ..., up to k + (something). The dependence is
still linear (or up to quadratic via (k,k) cross at m=2k), but at MULTIPLE strata.

At stratum m=k+1, c_k coefficient is −c_1 mod p (from B_2 · 2·c_1·c_k).
At stratum m=k+2, c_k coefficient picks up contributions from j=2 (2,k) and j=3 (1,1,k):
- j=2 (2,k): B_2 · 2·c_2·c_k = −c_2·c_k (linear in c_k, coefficient −c_2 mod p).
- j=3 (1,1,k): B_3 · 3·c_1²·c_k = c_1²·c_k / 2 (linear in c_k, coefficient c_1²/2 mod p).
- TOTAL at stratum m=k+2: c_k coefficient = (−c_2 + c_1²/2) mod p, still LINEAR in c_k.

Pattern: **c_k appears at strata k+1, k+2, ..., 2k (where 2k = the (k,k) cross-term at j=2 makes c_k² appear)**.

Actually let me check: c_k² appears at m = 2k via j=2 (k,k) composition, with multinomial
coefficient 1, contribution B_2·c_k² = −c_k²/2. So at m=2k, c_k appears with degree 2.

For m < 2k: c_k appears LINEARLY (linear in c_k, with the coefficient being a polynomial
in lower digits c_1, ..., c_{k−1}).

For m = 2k: c_k appears QUADRATICALLY.

For m > 2k up to some upper bound: c_k can appear with higher degrees (cubic at m = 3k,
quartic at m = 4k via j ≥ 3 compositions).

**Bound on top stratum reached:** c_k can appear at stratum up to j·k, but j ≤ r and
stratum m ≤ r. So c_k appears at strata m=k+1 up to m = r (or earlier if 2k > r).

## Linear range of c_k

**c_k is LINEAR in c_k at strata m ∈ [k+1, 2k−1] (when 2k ≤ r+1, else up to m=r).**
**c_k is QUADRATIC at stratum m = 2k (if 2k ≤ r).**
**c_k is cubic at m = 3k (if 3k ≤ r).**

For the "outermost / top stratum" question relevant for Plancherel collapses: at
stratum m=r, the c_{r−1} (top digit) appears LINEARLY only — since c_{r−1}² requires
m ≥ 2(r−1) > r for r ≥ 3.

So **c_{r−1} appears ONLY at stratum m=r, LINEARLY.**

**c_{r−2} appears at strata k+1 = r−1 and k+2 = r (if 2(r−2) > r, i.e., r > 4; for r=4
have 2·(r−2) = 4 = r, so c_2 appears at m=3 linearly AND at m=4 with c_2² term).**

This is the key structural data for Phase 3.

## Output table — TOP-DIGIT linearity at each stratum for r = 4, 5, 6, 7

| r | strata | top digit at stratum r | c_{r−1} appearances | c_{r−1} polynomial structure |
|---|---|---|---|---|
| 3 | 2, 3 | c_2 (at m=3) | m=3 only | linear |
| 4 | 2, 3, 4 | c_3 (at m=4) | m=4 only | linear |
| 5 | 2, 3, 4, 5 | c_4 (at m=5) | m=5 only | linear |
| 6 | 2, 3, 4, 5, 6 | c_5 (at m=6) | m=6 only | linear |
| 7 | 2, 3, 4, 5, 6, 7 | c_6 (at m=7) | m=7 only | linear |

**The TOP DIGIT c_{r−1} is always LINEAR at the unique stratum where it appears (m=r).**

This is the structural ingredient (I2) at the top digit. The c_{r−1}-Fourier collapse
should work analogously to the c_2-collapse at r=3.

## What about the deeper digits c_2, c_3, ..., c_{r−2}?

These appear at MULTIPLE strata (linear at some, possibly quadratic at others, ...). The
Plancherel structure at the deeper digits is more complex.

**Pre-output of Phase 3:** the c_{r−1} collapse works cleanly (analog of r=3's c_2
collapse). The deeper-digit collapses need analysis. The KEY QUESTION is whether the
nested chain collapses properly, peeling off c_{r−1}, then c_{r−2}, ..., down to c_2.

## Numerical check (one cell): r=4, p=7

At r=4, p=7, take a = 1 (so C_a = L̃_7^{-1} mod 7^4). We need to verify the polynomial
expansion matches by hand.

(Hand check deferred to Phase 5 if Python is available.)

## Files

- HENSEL_APPROACH_A.md, HENSEL_PHASE_ARTICULATION.md — closed form and stratum analysis
- CHAIN_PHASE1_R3_RECAP.md — r=3 chain
- This document — digit expansion at r ≥ 4

## Next: Phase 3 tests each c_k for Fourier-collapse

For each higher digit c_k (k = 2, 3, ..., r−1), test:
(a) The LEADING APPEARANCE at stratum m = k+1 has c_k linear with coefficient −c_1.
    Plancherel-on-c_k at this stratum gives a clean length-p collapse → δ(c_1 ≡ 0 mod p).
(b) Subsequent appearances at strata m = k+2, k+3, ... (linear or quadratic in c_k).
    These are HIGHER-STRATUM contributions. The Plancherel at THE c_k stratum (m=k+1)
    is the dominant one; higher-stratum contributions appear in the OUTER PHASE after
    c_k collapse.
(c) Whether after the c_k collapse forces c_1 ≡ 0 (or another linear condition), the
    remaining phase admits the next-c_{k−1} collapse.
