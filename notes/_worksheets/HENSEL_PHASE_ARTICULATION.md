# HENSEL Phase 1 — Articulation of the Hensel-corrected phase at r ≥ 4

**Date:** 2026-05-11. Analyst: Wilson. Reports to: Nathan.

## Goal

Make the "Hensel-correction series" concrete enough at family level to (a) name what closure would look like and (b) feed Phase 2 Approach A. The deliverable here is symbolic, not numerical.

## Recap of the r=3 saddle (T78.6_p, c=1)

For prime p ≥ 3, r=3, q = p^{r+1} = p^4, period p^3, support {a ≡ 1 mod p in Z/p^r}:

- `C_a := a · L̃_p^{-1} mod p^r` (T78.5_p — bijection on support).
- `P_a(s) := p·s − C_a · L_p(1+ps) mod q`, with `L_p(1+ps) = Σ_{j=1}^{J_p} (-1)^{j-1}/j · (p·s)^j`.
- `s*(C_a) = (C_a − 1)/p mod p`.
- T78.6_p (r=3): `G_p(a) = p^{(r+1)/2} · e_q(P_a(s*(C_a)))`.

Explicit closed form (from PATH2_BILINEAR_FROM_CLOSED_FORM.md §"At r=3"):

> `P_a(s*) ≡ −p²·s*²/2 + p³·(s*³/6 − c_2·s*) mod p^4`

where `C_a = 1 + p·s* + p²·c_2 mod p^3` (the s* and c_2 are the base-p digits of (C_a−1)).

## Why the saddle works at r=3 (mechanistic recap)

Around the saddle `s*` (mod p), the Taylor expansion of P_a is:

> P_a(s* + h) = P_a(s*) + P'_a(s*)·h + (1/2)·P''(s*)·h² + (1/6)·P'''(s*)·h³ + (1/24)·P''''(s*)·h⁴ + ... mod q

Derivatives evaluated at the saddle s* (which satisfies 1+ps* ≡ C_a mod p²):

| Derivative | Formal series in 1/(1+ps) | v_p at saddle | Form at saddle, c=1 |
|---|---|---|---|
| P'_a(s) | `p·c − C_a · p / (1+ps)` | `≥ 2` mod q (≡ 0 mod p²) | small |
| P''_a(s) | `C_a · p² / (1+ps)²` | exactly 2 | `p² / C_a` |
| P'''_a(s) | `−2·C_a · p³ / (1+ps)³` | exactly 3 | `−2·p³ / C_a²` (note: extra factor `c=1` cancellation) |
| P''''_a(s) | `6·C_a · p⁴ / (1+ps)⁴` | exactly 4 | `6·p⁴ / C_a³` |

(The d^k P / ds^k coefficient at s* picks up the geometric-series factor `C_a · p^k · (sign·k!) / (1+ps)^k`, evaluated at 1+ps ≡ C_a mod p^k+1, so each derivative at s* drops down to `±k!·p^k/C_a^{k−1}`.)

**Gaussian saturation at r=3:** The quadratic term `(1/2)·P''(s*)·h² ≡ p²·h²/(2·C_a)` summed over h ∈ Z/p^3 is a length-p^3 sum of period-p^2 quadratic phases. By standard p-adic Gauss-sum theory, this sums to `p · √(p²) · (root of unity) = p² = p^{r−1}`. Combined with the saddle-fluctuation prefactor, the inner sum reaches magnitude `p^{(r+1)/2} = p²` — Plancherel-saturation.

**Cubic absorption at r=3:** The cubic term `(1/6)·P'''(s*)·h³ ≡ −p³·h³/(3·C_a²)` mod q=p^4 is `e_{p^4}(p³·X) = e_p(X)`. This depends on h mod p ONLY (since (h+p)³ ≡ h³ mod p). It splits as a class-constant phase per h mod p — exactly what produces the `e_p(s*³/6 − c_2·s*)` term in the r=3 closed form.

**Quartic and higher are absorbed by q=p^4 precision:** P''''(s*) has v_p = 4, so `(1/24)·P''''(s*)·h⁴ ≡ 0 mod p^4`. Similarly for higher orders — they drop out.

**So at r=3, the saddle integration TERMINATES exactly. T78.6_p holds rigorously.**

## What changes at r=4

At r=4, q = p^5, period = p^4, support size N = p^3. The saddle equation `c·(1+ps) ≡ C_a mod p^4` now has solution `s* mod p^3` (not mod p), with three p-adic digits.

Define the **Hensel-lifted saddle**:

> `s*(r=4)(C_a) := (C_a − 1)/p mod p^3 = s_0 + p·s_1 + p²·s_2`

where `s_0 = (C_a − 1)/p mod p` (the r=3 saddle), and `s_1, s_2 ∈ {0,...,p−1}` are the next base-p digits.

Equivalently, writing `C_a = 1 + p·c_1 + p²·c_2 + p³·c_3 mod p^4`:

> `s_0 = c_1,   s_1 = c_2,   s_2 = c_3`.

**The Hensel correction is just digit-extraction: there is no extra δ_1, δ_2 series — these "corrections" are simply higher base-p digits of (C_a − 1)/p.** This collapses the abstract `δ_k(C_a)`, `ε_k(C_a)` notation.

## Phase derivation at r=4

We compute `P_a(s*(r=4)) mod p^5` directly, expanding L_p(1+p·s*) with J_p ≥ 5.

Let `s* = s_0 + p·s_1 + p²·s_2 mod p^3`, so `p·s* = p·s_0 + p²·s_1 + p³·s_2 mod p^4`. Note `p·s* ∈ p·Z/p^4` (i.e., `v_p(ps*) = 1` generically when s_0 ≠ 0).

Each Taylor coefficient at s*:
- L_p(1+ps*) = ps* − (ps*)²/2 + (ps*)³/3 − (ps*)⁴/4 + (ps*)⁵/5 mod p^5.

Expand `(p·s*)^j` mod p^{5+1}: needed up to order 5.

Using `C_a = 1 + p·s* mod p^4` (T78.5_p, by definition of s* as (C_a−1)/p):

> `C_a · L_p(1+ps*) = (1+ps*) · L_p(1+ps*) mod p^5`

Let `y := p·s*`. Then `C_a · L_p(1+y) = (1+y) · (y − y²/2 + y³/3 − y⁴/4 + y⁵/5) mod p^5`.

Expand:
- `(1+y) · y = y + y²`
- `(1+y) · (−y²/2) = −y²/2 − y³/2`
- `(1+y) · (y³/3) = y³/3 + y⁴/3`
- `(1+y) · (−y⁴/4) = −y⁴/4 − y⁵/4`
- `(1+y) · (y⁵/5) = y⁵/5 + y⁶/5`  → `y⁶ ≡ 0 mod p^5` always when v_p(y) ≥ 1 (i.e., s_0 ≠ 0 generic). When s_0 = 0, y has v_p ≥ 2, so y⁶ has v_p ≥ 12, certainly 0 mod p^5.

So:
> `C_a · L_p(1+y) = y + y²/2 − y³/6 + y⁴·(1/3 − 1/4) − y⁵·(1/4 − 1/5) mod p^5`
>                `= y + y²/2 − y³/6 + y⁴/12 − y⁵/20 mod p^5`

(Hidden simplification: `(1+y)·log(1+y) = ...` has a nice form. The coefficients are `1, 1/2, −1/6, 1/12, −1/20, ...` — these are `1/(j(j−1))` for j ≥ 2 with appropriate signs. Specifically, `(1+y)·log(1+y) = y + Σ_{j≥2} (−1)^j · y^j / (j·(j−1))` — this is a known generating-function identity.)

Then:
> `P_a(s*) = p·s* − C_a · L_p(1+ps*) = y − [y + y²/2 − y³/6 + y⁴/12 − y⁵/20] mod p^5`
>        `= −y²/2 + y³/6 − y⁴/12 + y⁵/20 mod p^5`

**Substituting y = p·s*:**

> **`P_a(s*(r=4)) ≡ −p²·s*²/2 + p³·s*³/6 − p⁴·s*⁴/12 + p⁵·s*⁵/20 mod p^5`**

The last term is `≡ 0 mod p^5` since p⁵·s*⁵/20 has v_p ≥ 5 (denominator 20 = 4·5; for p ≥ 5, 1/20 might pick up −1 from 1/5, but with p⁵ in front the net v_p ≥ 4; close call — let me re-check at p=5).

For p=5: p⁵/20 = 5⁵/20 = 3125/20 = 156.25 — non-integer at p-adic level. The "issue" is denominator includes 5 = p. So `p⁵/20 = p⁴/4`, which has v_p = 4. Hence `p⁵·s*⁵/20 = p⁴·s*⁵/4` mod p^5. Non-zero! So the term DOES contribute.

For p ≥ 7: 1/20 has v_p = 0, so p⁵·s*⁵/20 ≡ 0 mod p^5. Drops.

For p=3: 1/20 has v_3 = 0 (20=4·5, neither is 3), so drops.

So the p=5 case has an extra term. **Note this for adversarial check (A2) — the Hensel correction is NOT purely uniform across p.**

Generic formula (for p ∉ {2, 5}):
> **`P_a(s*(r=4)) ≡ −p²·s*²/2 + p³·s*³/6 − p⁴·s*⁴/12 mod p^5`**

Compare to r=3 formula:
> `P_a(s*(r=3)) ≡ −p²·s*²/2 + p³·s*³/6 mod p^4` ... but wait — at r=3 there was the `−c_2·s*` term I haven't accounted for here.

**Crucial subtlety.** The r=3 formula in PATH2_BILINEAR §"Attempt A" derived `P_a(s*) ≡ −p²·s*²/2 + p³·(s*³/6 − c_2·s*) mod p^4`. The `−c_2·s*` term comes from expanding `C_a = 1 + p·s_0 + p²·c_2 mod p^3` and tracking carefully — `s* = s_0` mod p (only the leading digit), and the `c_2·s_0` cross term arises in `C_a · L_p(1+p·s_0)`.

The clean computation I just did above uses `s* = s_0 + p·s_1 + p²·s_2 mod p^3` and `C_a = 1 + p·s* mod p^4` *exactly*. Let me reconcile.

**Reconciliation:** At r=3, the saddle `s_0 = (C_a−1)/p mod p` satisfies `1 + p·s_0 ≡ C_a mod p²` only. The TRUE identity `1 + p·s* = C_a` only holds when `s* = (C_a−1)/p` (no mod). When we evaluate `P_a(s_0)` (using only the leading saddle digit) instead of `P_a(s*(r=4)) = P_a((C_a−1)/p mod p^3)`, we get residual cross-terms involving the higher digits c_2, c_3.

So the question becomes: **at r=4, does T78.6_p(r=4) hold with s* = (C_a−1)/p mod p^3 (the Hensel-lifted saddle)?**

If yes, then `P_a(s*(r=4))` is the closed-form phase, and the formula I derived above (`−p²s*²/2 + p³s*³/6 − p⁴s*⁴/12 mod p^5`) IS the explicit Hensel-corrected phase. No δ_k, ε_k series — just the Hensel-lifted saddle plugged into a polynomial of degree 4 in s*.

This is the candidate closed form for Approach A.

## Phase derivation at r=5, 6 (pattern)

At r=5: q=p^6, period p^5, s*(r=5) = (C_a−1)/p mod p^4. The same calculation gives:

> **`P_a(s*(r=5)) ≡ −p²·s*²/2 + p³·s*³/6 − p⁴·s*⁴/12 + p⁵·s*⁵/20 mod p^6`**

(with the p=5 small-prime caveat on the last term).

At r=6: q=p^7, period p^6, s*(r=6) = (C_a−1)/p mod p^5. Same pattern, one more term:

> **`P_a(s*(r=6)) ≡ −p²·s*²/2 + p³·s*³/6 − p⁴·s*⁴/12 + p⁵·s*⁵/20 − p⁶·s*⁶/30 mod p^7`**

**General pattern (conjectural):**

> **`P_a(s*(r)) ≡ Σ_{j=2}^{r} (−1)^{j−1} · p^j · s*^j / (j·(j−1)) mod p^{r+1}`**

valid when p > j for all terms (no small-prime denominator issues), i.e., for p > r.

Verify: coefficients `1/(j(j−1))` matches the `(1+y)·log(1+y)` expansion `Σ_{j≥2} (−1)^j · y^j / (j·(j−1))` (after sign-flip from p·s* − C_a · L; sign is `−(−1)^j = (−1)^{j+1}`).

Actually let me double-check: y² coefficient is `+1/(2·1) = 1/2`, my expansion has `−y²/2` after the minus from `P = p·s − C_a · L`. So sign is `−`, matching `(−1)^{j−1} · 1/(j(j−1))` for j=2 gives `(−1)^1 / 2 = −1/2`. ✓

**So at family level (p > r), the candidate Hensel-lifted closed form is:**

> **`G_p(a) (r) =? p^{(r+1)/2} · e_{p^{r+1}}(Σ_{j=2}^{r} (−1)^{j−1} · (p·s*)^j / (j·(j−1)))`**

where `s* = s*(r)(C_a) = (C_a − 1)/p mod p^{r−1}`.

This is the **HENSEL-LIFTED CLOSED-FORM CANDIDATE**. Phase 2 Approach A tests its exactness.

## Identifying cancellation vs accumulation terms

For Approach B's cancellation hunt, decompose the phase by p-adic depth.

In `P_a(s*(r)) = Σ_{j=2}^{r} (−1)^{j−1} · p^j · s*^j / (j·(j−1))`, plug `s* = s_0 + p·s_1 + p²·s_2 + ... mod p^{r-1}` and look at terms of each p-power mod p^{r+1}:

- **Order p² (`j=2`):** `−p²·s*²/2 = −p²·(s_0 + p·s_1 + p²·s_2 + ...)²/2`
  - `s_0²`: order p² (leading)
  - `2·s_0·s_1·p`: order p³
  - `2·s_0·s_2·p² + s_1²·p²`: order p^4
  - `2·s_0·s_3·p³ + 2·s_1·s_2·p³`: order p^5
  - etc.

- **Order p³ (`j=3`):** `+p³·s*³/6`
  - `s_0³`: order p³ (leading)
  - `3·s_0²·s_1·p`: order p^4
  - `3·s_0²·s_2·p² + 3·s_0·s_1²·p²`: order p^5
  - etc.

- **Order p⁴ (`j=4`):** `−p⁴·s*⁴/12`
  - `s_0⁴`: order p^4
  - `4·s_0³·s_1·p`: order p^5
  - etc.

- **Order p⁵ (`j=5`):** `+p⁵·s*⁵/20` (drops mod p^5 generically, contributes at higher r)

**Collecting by p-power in P_a(s*) mod p^{r+1}:**

The phase at order p^k mod p^{r+1} reads as e_{p^{r+1-k}}(coefficient). Higher-k contributions are "fine" phases on smaller mods.

Specifically:
- p² coefficient: −s_0²/2 → depends on s_0 ONLY → constant on s_0-class
- p³ coefficient: −s_0·s_1 + s_0³/6 → depends on (s_0, s_1) → "fine" within s_0-class
- p^4 coefficient: −(s_0·s_2 + s_1²/2) + (s_0²·s_1)/2 − s_0^4/12 → depends on (s_0, s_1, s_2)
- p^5 coefficient: similar pattern; depends on (s_0, ..., s_3)

The "linear-in-c_2 term" at r=3 was `−c_2·s_0` (with c_2 ≡ s_1 in the new notation). This is the p^3·(−s_0·s_1) coefficient × (−1). Same here. **The phase at r=4, 5, 6 is NOT purely cubic — it gains higher-order terms but their structure is class-by-class linear in the highest digit at each depth.**

### Cancellation structure (for Approach B)

Within each s_0-class, sum over (s_1, s_2, ..., s_{r-2}) = (the inner support). The phase at depth k=3,4,...,r+1 is **piecewise multilinear in the higher digits**.

Inner-Plancherel candidates:
- Depth 3 phase: `−s_0·s_1` is **linear in s_1** at fixed s_0 → Inner-Plancherel on s_1 (length p) collapses cleanly.
- Depth 4 phase: `−s_0·s_2 − s_1²/2 + (s_0²·s_1)/2 − s_0^4/12` — depends on (s_1, s_2), **linear in s_2** at fixed (s_0, s_1) → Inner-Plancherel on s_2 collapses; the s_1² term is quadratic, requires a Gauss-sum sub-bound.
- Depth 5: linear in s_3, quadratic in s_2, cubic in s_1.

**The structural pattern:** at depth k, the phase is **linear in the (k−2)-th digit** (the newest digit appearing). Earlier digits appear with degrees 2, 3, ... in a controlled way.

This is exactly the structure of a **Vinogradov mean-value problem at depth r−2**. The "decoupled" formulation matches the Bourgain-Demeter-Guth 2016 setup at degree (r−1).

### Verdict for cancellation structure

At depth k, the inner-Plancherel on digit s_{k−2} should give a length-p save (factor `p`), IF the phase is exactly linear in s_{k−2} at fixed earlier digits. This requires no quadratic term in s_{k−2}, which is what the explicit expansion gives (s_{k−2}'s leading appearance is `−s_0·s_{k−2}·p^k` from j=2 term — linear).

**Implication:** if T78.6_p generalizes to r ≥ 4 as conjectured (with s* = (C_a−1)/p mod p^{r−1}), then a NESTED inner-Plancherel applies, peeling off ONE digit at a time. Each peeling gives a factor p save. After r−1 peelings, the bound becomes:

> `|T_p| ≤ p^{r-1} · |T_p^{outer}|` = N · (bound on s_0-class outer sum)

where the outer sum (over s_0 ∈ Z/p, the leading saddle digit) is the same `Σ |D_p(a_0(s*), p²)| ≤ p + log p` from PATH2_BILINEAR Attempt G+.

Combined: `|T_p| ≤ N · 2 = 2N` family-level, NO log N factor. **This would close the bound to strict 2√N.**

But this REQUIRES T78.6_p at r ≥ 4 (the Hensel-lifted saddle exactness). That's the gate.

## Cancellation-vs-accumulation decomposition

| Term in P_a(s*) | Depth (p-power) | Structure | Inner-Plancherel save? |
|---|---|---|---|
| −s_0²·p²/2 | p^2 | constant per s_0-class | N/A (outer) |
| −s_0·s_1·p³ | p^3 | linear in s_1 | YES — factor p save |
| s_0³·p³/6 | p^3 | constant per s_0-class | N/A |
| −s_0·s_2·p^4 | p^4 | linear in s_2 (highest digit) | YES — factor p |
| −s_1²·p^4/2 | p^4 | quadratic in s_1 | Quadratic Gauss sum, factor √p (not p) |
| (s_0²·s_1)·p^4/2 | p^4 | linear in s_1 (already integrated) | redistributes within Inner |
| −s_0^4·p^4/12 | p^4 | constant per s_0-class | N/A |
| ... at higher depths | p^k | linear in s_{k−2} (top digit), polynomial in lower digits | linear top gives factor p; lower-digit polynomial may need finer analysis |

**Cancellation-bearing terms:** the linear-in-top-digit terms `−s_0·s_{k−2}·p^k` at each depth.

**Accumulation-bearing terms:** the multi-digit cross-terms at the SAME depth that appear in the inner-Plancherel residual (they don't kill the bound but they obstruct trivial linearity).

## Family-level p-blindness

All steps are p-blind EXCEPT:
- (Q4-redux from PATH2_FAMILY_EXTENSION): the `1/j(j−1)` denominators in the expansion of `(1+y)log(1+y)`. For j(j−1) to be coprime to p, need p > j. **At depth r+1 we need p > r.** For p=3, r ≥ 3 already triggers a small-prime issue at j=3 (1/3·2 = 1/6, v_3 = −1). The `(p·s*)^j / (j(j−1))` for j=3, p=3 is `p^3·s*^3 / 6 = 27·s*^3 / 6 = 9·s*^3/2`. v_3 = 2. So mod p^4 = 81, this is `9·s*^3/2 mod 81` — non-zero, well-defined since gcd(2, 81)=1.

So small-prime quirks DON'T break the formula; they just shift the v_p of each term by `−v_p(j(j−1))`. The `closed form` is still valid as a p-adic identity; the **p-power structure** of each digit-coefficient gets shifted, which changes which inner-Plancherel digit is which.

For p=3, r=4 specifically: the j=3 term `p³·s*³/6` has v_p = 2 (not 3), so it MERGES INTO the depth-2 stratum with `−p²·s*²/2`. The "saddle phase mod p²" gets the cubic contribution as well. This affects the inner-Plancherel layering but doesn't break it.

## What this articulation gives Approach A

The candidate closed form `G_p(a) (r=4) =? p^{(r+1)/2} · e_{p^5}(P_a(s*(r=4)))` is now explicit:

> `P_a(s*(r=4)) ≡ −p²·s*²/2 + p³·s*³/6 − p⁴·s*⁴/12 mod p^5` (for p ∉ {2,5}; p=5 picks up extra term)
>
> where `s* = (C_a − 1)/p mod p^3`.

This is a TESTABLE prediction at the cell level. Phase 2 Approach A computes G_p(a) directly and compares to this prediction at p ∈ {3,5,7,11} × r ∈ {4,5,6}.

If the prediction holds to machine precision, **Approach A produces a candidate closed form and Approach B/C aren't needed**.

## What this articulation gives Approach B

If A's candidate misses (T78.6_p saddle exactness fails at r=4), the closed form gains residual terms `ε_k(C_a)·p^{r+k}`. From the calculation above, the residual at r=4 would come from the saddle equation `c·(1+ps*) ≡ C_a mod p^4` being satisfied only mod p^3 (not mod p^4), introducing a leftover term `≈ P'_a(s*)·δ` where `δ = (C_a − (1+ps*))/p` and v_p(δ) ≥ 2 by saddle definition. Tracking this gives the recursive series for δ_k.

This is the technically demanding step; deferred unless A misses.

## What this articulation gives Approach C

If A and B miss, the bilinear |Σ 1̂(p·a) · F̂_p(p·a)|² as a 6th moment:

`|Σ 1̂(p·a) · F̂_p(p·a)|^2 = Σ_a Σ_b 1̂(p·a) · conj(1̂(p·b)) · F̂_p(p·a) · conj(F̂_p(p·b))`

Plugging the family-level G_p(a) = p^{(r+1)/2} · e_q(P_a(s*(C_a))) (assuming the Hensel-lifted form holds approximately), the F̂_p·conj(F̂_p) phase is `e_q(P_a(s*_a) − P_b(s*_b))`. The cubic-in-(s_a, s_b) part suggests a Vinogradov-mean-value structure on (s_0, s_1, ..., s_{r-2}) at degree r−1.

Bourgain-Demeter-Guth 2016 proved the conjectured VMV exponent for degree-3 polynomials (cubic case). For our degree r−1 ≥ 3 phase, BDG gives the right scaling — but **only "right exponent", not explicit constants**. The constant matters for the bilinear closure to strict 2√N.

Deferred unless A, B both miss.

## Adversarial flags pre-loaded for Phase 4

- **(A1) Triangle inequality re-examination:** If Approach A produces exact closed form at r=4, the inner-Plancherel argument peels off ONE digit (factor p save). Combining across all r−1 inner digits: `|T_p| ≤ |Inner_outer| ≤ N + p log p ≤ 2N` family-level uniformly. **Strict 2√N. The log N IS removable** IF saddle exactness holds.
- **(A2) Hensel correction order matching:** Each digit gets ONE inner-Plancherel peel. At r=4 we peel off s_1, s_2 (two inner digits); at r=5 we peel s_1, s_2, s_3 (three). The recursive structure means the same proof handles all r ≥ 4 IF the saddle is Hensel-exact.
- **(A3) Empirical anchor:** R79b shows |K|/√N ∈ [0.7, 2.7] at p=3, r=8..20. Prediction from closed form: |T_p|/N ≤ 2 uniformly, so |K|/√N ≤ const · 2 ≤ a few. **Match.**
- **(A4) VMV literature check:** BDG 2016 gives sharp VMV exponents for cubic phases on integer-cube intervals. Our setup is a finite-modulus principal-unit Gauss sum — different. BDG translates only structurally. If Approach C is needed, expect "structural-match-only" verdict similar to Milicevic-Banks.

## Summary

Phase 1 produced:
1. A specific **candidate closed form** for the Hensel-lifted phase at r ≥ 4: `P_a(s*(r)) ≡ Σ_{j=2}^{r} (−1)^{j−1} · (p·s*)^j / (j·(j−1)) mod p^{r+1}` with `s* = (C_a−1)/p mod p^{r−1}`. The series `(1+y)·log(1+y) = Σ y^j/(j(j−1))` underlies this.
2. A specific **inner-Plancherel layering**: each digit s_k of s* admits a length-p Plancherel save, peeling depth-by-depth from the outer s_0-saddle bound.
3. A specific **family-level closure prediction**: |T_p| ≤ 2N strict (no log N) at r ≥ 4, IF the Hensel-lifted closed form holds.
4. **Small-prime caveats:** For p ∈ {2, 5} at certain r, denominators 1/j(j−1) pick up small-prime factors that shift p-powers; the closed form survives but the digit-layering shifts.

Phase 2 Approach A now has a concrete prediction to test.
