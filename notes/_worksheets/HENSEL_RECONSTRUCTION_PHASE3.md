# HENSEL_RECONSTRUCTION_PHASE3 — Phase polynomial P_a(s*(r))

**Date:** 2026-05-11. Independent re-derivation; NO use of original Hensel-lift files.

## Setup (carried from Phase 2)

Saddle: `s*(r) = (C_a − 1)/p mod p^{r-1}`. Formal-p-adic identity: ps* = C_a − 1, hence 1 + ps* = C_a.

Phase polynomial:
> P_a(s*) = ps* − C_a · L_p(1 + p·s*)    mod q = p^{r+1}

Substituting ps* = C_a − 1 in the formal-p-adic ring:
> P_a(s*) = (C_a − 1) − C_a · L_p(C_a)

But L_p is the **truncated** Cochrane log. To J_p terms:
> L_p(C_a) = Σ_{j=1}^{J_p} (−1)^{j-1}/j · (C_a − 1)^j = Σ_{j=1}^{J_p} (−1)^{j-1}/j · (p·s*)^j

In the formal-Z_p ring, L_p(C_a) = log(C_a) (the convergent p-adic log on principal units), with truncation introducing O(p^{J_p+1}). Since J_p ≥ r, this is O(p^{r+1}) ≡ 0 mod q. So we can compute formally and reduce.

## Phase polynomial in closed form

Define
> **M(y) := y − (1+y)·log(1+y)** (formal Z_p series, valid for y ∈ p·Z_p)

Then P_a(s*) = (C_a − 1) − C_a · log(C_a) = M(C_a − 1) = **M(p·s*)** mod q.

## Series coefficients of M

Computed in Phase 2:
> **M(y) = Σ_{j ≥ 2} (−1)^{j−1}/(j·(j−1)) · y^j**

Derivation (redo for cleanliness):

log(1+y) = Σ_{j ≥ 1} (−1)^{j-1}/j · y^j.

(1+y)·log(1+y) = log(1+y) + y·log(1+y)
              = Σ_{j ≥ 1} (−1)^{j-1}/j · y^j + Σ_{j ≥ 1} (−1)^{j-1}/j · y^{j+1}

Coefficient of y^j in (1+y)·log(1+y):
- For j = 1: (−1)^0/1 = 1.
- For j ≥ 2: (−1)^{j-1}/j + (−1)^{j-2}/(j-1) = (−1)^{j-1} · [1/j − 1/(j-1)] = (−1)^{j-1} · [(j-1 − j)/(j(j-1))] = (−1)^{j-1} · [−1/(j(j-1))] = (−1)^j / (j(j-1)).

So (1+y)·log(1+y) = y + Σ_{j ≥ 2} (−1)^j/(j(j-1)) · y^j.

M(y) = y − (1+y)·log(1+y) = − Σ_{j ≥ 2} (−1)^j/(j(j-1)) · y^j = Σ_{j ≥ 2} (−1)^{j-1}/(j(j-1)) · y^j.

**Coefficient pattern:**
- j=2: (−1)^1/(2·1) = −1/2 → coefficient is +(−1)^{j-1}/(j(j-1)) = +(−1)/(2) = −1/2. So −1/2 at j=2.
- j=3: (−1)^2/(3·2) = 1/6 → coefficient +(−1)^{j-1}/(j(j-1)) = +(1)/(6) = 1/6. So +1/6 at j=3.
- j=4: (−1)^3/(4·3) = −1/12 → coefficient −1/12.
- j=5: 1/20.
- j=6: −1/30.

So:
> **M(y) = −y²/2 + y³/6 − y^4/12 + y^5/20 − y^6/30 + ...**

The factor (j(j-1)) in the denominator is the "double factorial-like" pattern arising naturally from differentiating (1+y)·log(1+y).

## Substitute y = p·s*(r)

> **P_a(s*) = M(p·s*) = Σ_{j=2}^{?} (−1)^{j-1}/(j(j-1)) · (p·s*)^j  mod p^{r+1}**

Where the truncation `?` is determined by v_p((p·s*)^j / (j(j-1))) ≥ r+1, which (for s* a unit-level quantity) means j + v_p((j(j-1))^{-1}) ≥ r+1.

For p ≥ 5: v_p(j(j-1)) = 0 for j ≤ p. So truncation at j = r+1 (since (p·s*)^{r+1} has v_p ≥ r+1).

For p = 3: 1/(j(j-1)) can have v_3 = −1 (when 3 | j(j-1)), so terms at j = r+1 might still contribute. The truncation is j where j − v_3(j(j-1)) ≥ r+1.

In either case the truncation is "exactly enough digits of s* to cover the precision r+1", consistent with s* being known mod p^{r-1} (which gives p·s* known mod p^r — falling short by one digit for p^{r+1} precision, but the (p·s*)^j powers absorb this).

## Verification at r=3

At r=3, we want P_a(s*) mod p^4.

M(p·s*) = −(ps*)²/2 + (ps*)³/6 + O((ps*)^4)

(ps*)² = p²·s*². s* = (C_a − 1)/p mod p² = t_1 + p·t_2 (mod p²). p·s* = p·t_1 + p²·t_2 mod p³. (p·s*)² = p²·t_1² + 2·p³·t_1·t_2 mod p^4.

(ps*)³ = p³·t_1³ + O(p^4) mod p^4.

So:
M(p·s*) mod p^4 = −(p²·t_1² + 2p³·t_1·t_2)/2 + p³·t_1³/6 + O(p^4)
              = −p²·t_1²/2 − p³·t_1·t_2 + p³·t_1³/6
              = −p²·t_1²/2 + p³·(t_1³/6 − t_1·t_2)

**Match with Phase 1 Step 3:** Phase 1 gave P_a(s*) ≡ −p²·t_1²/2 + p³·(t_1³/6 − t_1·t_2) mod p^4. ✓ **Identical.**

## Verification at r=4

At r=4, want P_a(s*) mod p^5. s* = t_1 + p·t_2 + p²·t_3 mod p^3 (3 digits). p·s* = p·t_1 + p²·t_2 + p³·t_3 mod p^4.

(p·s*)² = (p·t_1 + p²·t_2 + p³·t_3)² mod p^5
        = p²·t_1² + 2·p³·t_1·t_2 + p^4·(t_2² + 2·t_1·t_3) + O(p^5)

(p·s*)³ = p³·t_1³ + 3·p^4·t_1²·t_2 + O(p^5)

(p·s*)^4 = p^4·t_1^4 + O(p^5)

(p·s*)^5 = O(p^5) — vanishes mod p^5.

M(p·s*) = −(p·s*)²/2 + (p·s*)³/6 − (p·s*)^4/12 + O(p^5)
       = −[p²·t_1² + 2p³·t_1·t_2 + p^4·(t_2² + 2t_1·t_3)]/2
         + [p³·t_1³ + 3p^4·t_1²·t_2]/6
         − [p^4·t_1^4]/12
       = −p²·t_1²/2 − p³·t_1·t_2 − p^4·(t_2² + 2t_1·t_3)/2
         + p³·t_1³/6 + p^4·t_1²·t_2/2
         − p^4·t_1^4/12

Collect by level:
- p²: −t_1²/2 ✓
- p³: −t_1·t_2 + t_1³/6 ✓ (matches Phase 2 derivation)
- p^4: −(t_2² + 2 t_1·t_3)/2 + t_1²·t_2/2 − t_1^4/12 = −t_2²/2 − t_1·t_3 + t_1²·t_2/2 − t_1^4/12 ✓ (matches Phase 2)

**Perfect match with the direct expansion in Phase 2.** The closed form M(p·s*) = Σ_{j ≥ 2} (−1)^{j-1}/(j(j-1)) · (p·s*)^j reproduces the phase polynomial exactly.

## Connection to "(1+y)·log(1+y)" pattern

The function M(y) = y − (1+y)·log(1+y) has its negative as the integrand of the dilogarithm:
> M'(y) = 1 − log(1+y) − 1 = −log(1+y)

So M is related to the antiderivative of −log(1+y). Specifically:
> M(y) = ∫_0^y −log(1+t) dt − ∫_0^y 0 dt   …  wait

Actually:
> d/dy [y − (1+y)·log(1+y)] = 1 − log(1+y) − (1+y)·(1/(1+y)) = 1 − log(1+y) − 1 = −log(1+y).

So M(y) = ∫_0^y M'(t) dt = ∫_0^y −log(1+t) dt.

Integration: ∫_0^y log(1+t) dt = (1+y)·log(1+y) − y. So M(y) = y − (1+y)·log(1+y) = −∫_0^y log(1+t) dt. ✓

**M(y) = −∫_0^y log(1+t) dt** is an elegant closed-form interpretation. The coefficient (−1)^{j-1}/(j(j-1)) of y^j arises from integrating term-by-term: ∫ (−1)^{j-1}/j · t^j dt evaluated at y gives (−1)^{j-1}·y^{j+1}/(j(j+1)). Shift index j+1 → j: (−1)^{j-2}·y^j/((j-1)·j) = −(−1)^{j-1}·y^j/(j(j-1)). With the leading sign flip from M = − ∫ log(1+t) dt: M(y) = Σ_j (−1)^{j-1}·y^j/(j(j-1)). ✓

## Final phase polynomial

> **P_a(s*(r)) ≡ Σ_{j ≥ 2}^{r+1} (−1)^{j-1}/(j(j-1)) · (p·s*(r))^j   mod p^{r+1}**

with `s*(r) = (C_a − 1)/p mod p^{r-1}`.

Or equivalently:
> P_a(s*(r)) ≡ M(p·s*(r))    mod p^{r+1}

where **M(y) = y − (1+y)·log(1+y)** is the "(1+y)·log(1+y)-deficit" series.

## Phase 3 conclusion

The phase polynomial at the saddle has a **clean closed form**:
- Coefficient pattern: `(−1)^{j-1} / (j·(j-1))` for j ≥ 2.
- Variable: `(p·s*(r))^j` where p·s* = C_a − 1.
- Equivalent to: `M(p·s*) = (C_a − 1) − C_a·log(C_a)`.

**This matches the structural form anticipated in the original Hensel-lift claim's (1+y)·log(1+y) coefficient pattern (per the system message's description).**

## Notes on truncation

The series M(p·s*) truncates implicitly at j such that v_p((p·s*)^j / (j(j-1))) ≥ r+1. For p ≥ 5 and small r, this is at j = r+1 (each (p·s*)^j contributes p^j times unit, divided by j(j-1) which is a unit for p ≥ 5 when j ≤ p). For p = 3 and growing r, the 1/3-denominators in j(j-1) can lower the v_p of the j-th term, extending the truncation.

The original Hensel claim (per system message) states truncation at j = r:
> P_a(s*(r)) = Σ_{j=2}^{r} (−1)^{j−1}·(p·s*(r))^j / (j·(j−1)) mod p^{r+1}

My derivation gives Σ_{j=2}^{r+1} (or higher for p=3 corner cases). The TRUNCATION DEPTH differs by 1 from the system-message description. Let me check if r+1 vs r matters.

At j = r+1 (with p ≥ 5): (p·s*)^{r+1} has v_p ≥ r+1 (since v_p(p·s*) ≥ 1). So contribution to P_a mod p^{r+1} is ≡ 0 — vanishes by reduction.

So **j=r+1 contributes nothing mod p^{r+1}**, the truncation at j = r is functionally equivalent to j = r+1 for p ≥ 5. ✓

For p = 3 at low r, the j=r+1 might contribute via the 1/3 denominator, but typically not at the lowest relevant r.

**Match with system-message form: yes, modulo the harmless choice of truncation endpoint (j = r or r+1 — both give the same result mod q for p ≥ 5).**

## Possible discrepancy — direction of comparison

System message: "P_a(s*(r)) = Σ_{j=2}^r (−1)^{j-1}·(p·s*(r))^j / (j·(j-1)) mod p^{r+1}"

My derivation: Σ_{j=2}^{r+1} same coefficients.

Note: at r=3, my Σ_{j=2}^{4} vs claim's Σ_{j=2}^{3}. The j=4 term: (−1)^3·(ps*)^4/12 = −(ps*)^4/12. At r=3, (ps*)^4 has v_p ≥ 4. Divided by 12 (a unit for p ≥ 5, p^{-1} for p=3): for p ≥ 5, v_p((ps*)^4/12) ≥ 4 = r+1, so vanishes mod p^4. For p=3: v_3((ps*)^4/12) = 4 − 1 = 3, which is mod p^4 NON-trivial.

Hmm so at p=3 r=3, my j=4 term contributes mod p^4, but the system-message form truncates at j=3 and misses this.

Let me check Phase 1 explicit derivation. Phase 1 Step 3 gave:
P_a(s*) at r=3 = −p²·t_1²/2 + p³·(t_1³/6 − t_1·t_2) mod p^4. 

This used J_p = 3 (truncated log to 3 terms). The j=4 term would come from L_p truncation level 4, but R78.6 says J=3 at q=3 r=3.

OK so at p=3 r=3, the truncated L_p has only j=1,2,3 terms. So my "j up to r+1" framing was overgeneral — the actual truncation comes from J_p of the Cochrane log, which is `max j with j − v_p(j) < r+1`.

For p=3 r=3: j=3: 3−1=2<4 ✓. j=4: 4<4 fails. So J_3=3.

So the j up to 3 only (not 4) at p=3 r=3. **My "M(p·s*)" expansion is actually capped at j = J_p, not r+1.**

For p ≥ 5 r=3: j=3: 3<4 ✓. j=4: 4<4 fails. So J_p=3. **At r=3, J_p=3 for all p ≥ 3.**

At r=4: J_p (max j with j − v_p(j) < 5):
- p=3: j=4: 4<5 ✓. j=5: 5<5 fails. j=6: 6−1=5<5 fails. So J_3=4.
- p≥5: j=4: 4<5 ✓. j=5: 5<5 fails. So J_p=4.

**At r ≥ 2, J_p = r for p ≥ 5 (since v_p(j) = 0 for j < p, j−v_p(j) = j; cutoff j = r+1 fails first, so J_p = r).**

**At r ≥ 2 and p = 3, J_p ≈ r for small r but eventually J_3 grows faster.**

So the truncation in M(p·s*) for the saddle phase is at j = J_p, which equals r for p ≥ 5 r ≤ p-1 (the generic regime).

**This matches the system-message form: Σ_{j=2}^{r} (with the convention J_p = r).**

## Final answer

> **P_a(s*(r)) ≡ Σ_{j=2}^{J_p} (−1)^{j-1}/(j(j-1)) · (p·s*(r))^j   mod p^{r+1}**

with J_p = r for p ≥ 5 (generic). Equivalently `P_a(s*) ≡ M(p·s*) mod p^{r+1}` with M(y) = y − (1+y)·log(1+y) truncated to J_p terms.

Matches the (1+y)·log(1+y) coefficient pattern claimed in the original.

## Files
- HENSEL_RECONSTRUCTION_PHASE3.md (this)
