# HENSEL_RECONSTRUCTION_PHASE2 — Saddle extension to r ≥ 4 (independent derivation)

**Date:** 2026-05-11. Independent re-derivation; NO use of original Hensel-lift files.

## Setup

From Phase 1 at general r ≥ 2:

- q = p^{r+1}, period = p^r
- P_a(s) = ps − C_a · L_p(1 + ps) mod q, with L_p the truncated p-adic log to J_p terms.
- J_p satisfies: J_p = max j with j − v_p(j) < r+1.
- The saddle equation `dP_a/ds ≡ 0` should now hold modulo a higher power of p than at r=3.

For G_p(a) to be evaluable in closed form at r ≥ 4 via stationary phase, the saddle congruence must be "deep enough" — i.e., dP_a/ds vanishes modulo p^? where the orthogonality argument of Phase 1 closes the non-saddle sums.

## The Hensel question

At r=3, dP_a/ds ≡ p²·(s − t_1) + O(p³) (Phase 1 Step 2). Setting ≡ 0 mod p³ gave s ≡ t_1 mod p. The orthogonality then vanished non-saddle classes mod p^{r-1} = p².

At r=4: the orthogonality argument needs the non-saddle phase to be a non-trivial complete sum on Z/p^{r-1} = Z/p³. So we need dP_a/ds ≡ 0 mod **p^{r}** for the saddle, not just mod p³ — wait, let me think again.

**Careful re-derivation of the orthogonality requirement.**

Going back to Phase 1 Step 5 logic at general r:
- Partition s ∈ Z/p^r by s mod p: s = s_0 + p·u, u ∈ Z/p^{r-1}.
- Fix s_0; compute P_a(s_0 + p·u) − P_a(s_0) mod q.
- For the sub-sum over u to be a complete sum of a non-trivial additive character on Z/p^{r-1} (and hence vanish by orthogonality unless trivial), need the leading u-dependence to be `(p / p^{r-1}-divisible)·u·(unit)`.

Specifically: P_a(s_0 + p·u) − P_a(s_0) = dP_a/ds|_{s_0} · p·u + (1/2)·d²P_a/ds²|_{s_0} · p²·u² + ... mod q.

We need this mod q to be a non-trivial linear function of u in Z/p^{r-1} (so that Σ_u of e_q vanishes by orthogonality). The leading linear-in-u coefficient is `p · dP_a/ds|_{s_0}`. For Σ_u e_q(...) to vanish, this coefficient (read mod q) should be a non-zero unit times p^{r}/p^{r-1} = p · (representative in Z/p^{r-1}). Multiplied by u and read mod q, the orthogonality gives:

Σ_{u=0}^{p^{r-1}-1} e_q(p · dP_a/ds|_{s_0} · u) = p^{r-1} if p · dP_a/ds|_{s_0} ≡ 0 mod q = p^{r+1}, else 0.

So **vanishing condition** for the non-saddle sub-sum (just from leading linear-in-u term):
> dP_a/ds|_{s_0} ≡ 0 mod p^r  ⟺  saddle.

In Phase 1 we computed dP_a/ds = p(1 − C_a) + C_a·p²·s − C_a·p³·s² + ... For C_a = 1 + p·t_1 + p²·t_2 + ... + p^{r-1}·t_{r-1} (digit expansion mod p^r):

1 − C_a = −p·t_1 − p²·t_2 − ... − p^{r-1}·t_{r-1} mod p^r

So:
> dP_a/ds = p(1 − C_a) + C_a·p²·s − C_a·p³·s² + C_a·p^4·s³ − ...

The contribution at each order in p:

- p^2 level: p(1 − C_a) ≡ −p²·t_1 mod p³. Plus C_a·p²·s ≡ p²·s mod p³.
  Combined at p² level: p²·(s − t_1).
- p^3 level: from p(1 − C_a): −p³·t_2. From C_a·p²·s = (1 + p·t_1)·p²·s + ... : extra p³·t_1·s. From −C_a·p³·s²: −p³·s² (leading).
  Combined at p^3 level: p³·(−t_2 + t_1·s − s²).
- p^4 level: from p(1 − C_a): −p^4·t_3. From C_a·p²·s = (... + p²·t_2)·p²·s: p^4·t_2·s. From −C_a·p³·s² = −(1+p·t_1+...)·p³·s² : extra −p^4·t_1·s². From +C_a·p^4·s³: p^4·s³ (leading).
  Combined at p^4 level: p^4·(−t_3 + t_2·s − t_1·s² + s³).

**Pattern observed:**

At p^{k+1} level (k ≥ 1):
> coefficient_at_p^{k+1} of dP_a/ds = −t_k + t_{k-1}·s − t_{k-2}·s² + ... + (−1)^{k-1}·t_1·s^{k-1} + (−1)^k·s^k
> = Σ_{j=0}^{k} (−t)^... wait let me redo this combinatorially.

Actually a cleaner way: the truncated p-adic log derivative satisfies
> dL_p(1+ps)/d(ps) = Σ_{j≥1} (−1)^{j-1} · (ps)^{j-1} = 1 / (1 + ps)   (formal p-adic series)

So `dL_p(1+ps)/ds = p · d/d(ps)·L_p = p/(1+ps)` (formally).

Then `dP_a/ds = p − C_a · p/(1+ps) = p · (1 − C_a/(1+ps)) = p · (1 + ps − C_a) / (1+ps)`.

This is the cleaner formulation. Setting dP_a/ds = 0:
> p · (1 + ps − C_a) / (1+ps) = 0 ⟺ 1 + ps − C_a = 0 ⟺ ps = C_a − 1 ⟺ s = (C_a − 1)/p

**The exact (formal-p-adic) saddle equation is `s = (C_a − 1)/p`!**

Since C_a ≡ 1 mod p (so C_a − 1 ∈ p·Z_p), the saddle `s = (C_a − 1)/p` is a well-defined element of Z_p. Reducing mod p^r:

> **s*(r) := (C_a − 1)/p mod p^{r-1} ∈ Z/p^{r-1}**

(Where C_a is taken in its canonical lift to Z/p^r; (C_a − 1)/p ∈ Z/p^{r-1} since C_a − 1 is a multiple of p.)

**Note:** at r=3, this gives s*(3) = (C_a − 1)/p mod p², not just mod p. But we showed in Phase 1 that the LEADING residue `(C_a − 1)/p mod p = t_1` is what controls the saddle vanishing at r=3 — the higher-precision lift to mod p² doesn't matter at r=3 because the next-order terms in dP_a/ds are O(p^4) which is ≡ 0 mod q.

For r ≥ 4, we need the saddle TO HIGHER PRECISION. Specifically, dP_a/ds must vanish mod p^r, not just mod p^3. This corresponds to determining s* mod p^{r-1}, not just mod p.

## The exact saddle is "clean," not a series

**Claim (independently derived):**
> The saddle of P_a at r ≥ 2 is **s*(r) = (C_a − 1)/p mod p^{r-1}**, a clean closed-form integer division.

**Proof sketch.** The formal p-adic saddle equation is `s = (C_a − 1)/p` (derived above from dL_p/ds = p/(1+ps) and dP_a/ds = 0). Reducing this formal saddle to Z/p^{r-1} gives the saddle representative for the truncated polynomial P_a(s) mod q.

**Verification that the truncation respects this:** The truncated log L_p(1+ps) differs from the true p-adic log by terms of order p^{J_p+1}, which contribute to dP_a/ds at order p^{J_p+1−1} = p^{J_p}. For the saddle condition dP_a/ds ≡ 0 mod p^r to be unaffected by truncation, need J_p ≥ r. 

For r ≥ 4:
- p = 3: J_p = max j with j − v_3(j) < r+1. At r=4: j − v_3(j) < 5. j=3: 3−1=2<5 ✓. j=6: 6−1=5 NOT<5. j=4,5: 4<5✓, 5−0=5 NOT<5. So J_3 at r=4 is 4. Hmm wait j=5: 5 − v_3(5) = 5, not < 5, so j=5 fails. j=4: 4−0=4 < 5 ✓. So J_3(r=4) = 4.
- For p ≥ 5: J_p = r (since v_p(j)=0 for j < p, and j=r+1 fails j<r+1 for the cutoff). So J_p(r) = r at p ≥ 5 r ≤ p-1.

In both cases, J_p ≥ r for r ≥ 2 (verified for small r at p=3, generally for p≥5 at r ≤ p-1).

**For p=3 r ≥ 4, careful:** at r=4, J_3=4 (just barely ≥ r = 4). At r=5, J_3 needs j with j − v_3(j) < 6. j=5: 5<6 ✓. j=6: 6−1=5<6 ✓. j=7: 7<6 fails. So J_3(r=5) = 6. Hmm but that's J=6, which is > 5 = r.

OK so J_p ≥ r at all relevant r and p ≥ 3. **Truncation doesn't kill the saddle approximation.**

## Refining: what does s*(r) = (C_a − 1)/p mod p^{r-1} look like?

This is **digit extraction** of C_a:
- C_a ∈ Z/p^r, C_a ≡ 1 mod p, so write C_a = 1 + p·X mod p^r where X = (C_a − 1)/p ∈ Z/p^{r-1}.
- s*(r) = X = (C_a − 1)/p mod p^{r-1} — literally the upper r-1 digits of C_a (after stripping the leading "1 + p·" structure).

**At r=3, s*(3) = (C_a − 1)/p mod p² ∈ Z/p².** Phase 1 showed only the **leading digit** s*(3) mod p = t_1 matters for the saddle at r=3. But for r ≥ 4, the **full** s*(r) ∈ Z/p^{r-1} matters.

## Verification: does the saddle at r ≥ 4 give a single coherent class?

At r=4, q = p^5, period = p^4. We want the saddle s* ∈ Z/p^4 such that:
- For s ≡ s* mod p^{r-1} = p^3: P_a(s) ≡ P_a(s*) mod q (so all p^{r-1} terms in the saddle class add coherently).
- For s ≢ s* mod p^{r-1} = p^3: the sub-sum vanishes by orthogonality.

The number of saddle-class elements in Z/p^r = Z/p^4 is p^r / p^{r-1} = p. So there are p representative saddle classes mod p^{r-1}, but only ONE of these (the one matching s* = (C_a − 1)/p) gives the coherent phase. The other p − 1 classes vanish.

Wait — I have a potential issue. At r=3, I partitioned s by s mod p (Phase 1 partition was mod p, with p^{r-1} = p² terms per class). At r=4, do I partition by s mod p^{r-1} = p³ (with p terms per class), or by s mod p (with p^{r-1} = p³ terms per class)?

The orthogonality argument's "shift step" was p, giving p^{r-1} terms per class. This was at r=3. The relevant fact is: the sub-sum over u ∈ Z/p^{r-1} of e_q(linear-in-u shift) vanishes when the shift has v_p < r+1 and is non-zero.

**Generalized argument for r ≥ 4:**

P_a(s + p·u) − P_a(s) = p · dP_a/ds|_s · u + (p²/2) · d²P_a/ds²|_s · u² + ... mod q

For Σ_u e_q(...) over u ∈ Z/p^{r-1} to vanish:
- need the polynomial in u (mod q) to NOT be the zero polynomial in Z/p^{r-1}.

If dP_a/ds|_s is a unit times p^k for some k < r, then the linear-in-u term `p·dP_a/ds|_s·u` is a unit times p^{k+1} times u, which mod q gives a non-trivial character on Z/p^{r-1} → orthogonal sum vanishes (when k+1 < r+1, i.e., k < r, which means dP_a/ds|_s ≢ 0 mod p^r).

For Σ_u to NOT vanish (coherent addition), need dP_a/ds|_s ≡ 0 mod p^r — the EXACT saddle condition.

**So the saddle condition at r ≥ 4 is dP_a/ds ≡ 0 mod p^r, and the solution is exactly s ≡ (C_a − 1)/p mod p^{r-1}** (via the formal-p-adic factorization 1 + ps − C_a = 0).

**Number of saddle solutions mod p^r:** p (the p lifts of (C_a − 1)/p mod p^{r-1} to mod p^r).

## So how many "saddle representatives" are there?

There are p representatives s = s*(r) + p^{r-1} · v for v ∈ Z/p. All these have dP_a/ds ≡ 0 mod p^r. **They all give the SAME value of P_a(s) mod q** (modulo possible quadratic corrections — let's check).

P_a(s* + p^{r-1}·v) − P_a(s*) = p·dP_a/ds|_{s*}·(p^{r-1}·v) + (p²/2)·d²P_a/ds²|_{s*}·(p^{r-1}·v)² + ...

= p^r · dP_a/ds|_{s*} · v + (p^{2r-1}/2) · d²P_a/ds²|_{s*} · v² + ...

At saddle, dP_a/ds|_{s*} ≡ 0 mod p^r, so p^r · dP_a/ds|_{s*} ≡ 0 mod p^{2r} ≡ 0 mod q (for r ≥ 2).

The next term: p^{2r-1}/2 · d²P_a/ds²|_{s*} · v². For r ≥ 2, 2r-1 ≥ r+1 = order of q, so this term is also ≡ 0 mod q ... unless r=2, where 2r-1 = 3 = r+1.

Hmm so at r=2, the **quadratic** correction kicks in at level p^{r+1} = q — non-trivial.

**At r ≥ 3, all p saddle representatives give P_a(s) ≡ P_a(s*(r)) mod q.**

**At r = 2, the p saddle reps have phase varying by a Gaussian factor (the d²P_a/ds² · v² term).**

This explains the **r=2 vs r ≥ 3 split** seen in the empirical claim. The Gauss-sum factor at r=2 is the quadratic Gaussian integration over the p representatives.

## Counting the magnitude

G_p(a) = Σ_{s ∈ Z/p^r} e_q(P_a(s))

By orthogonality, only the p^{r-1}-element saddle sub-set contributes coherently (well, actually it's the p saddle representatives times p^{r-1}/p? Let me recount).

I think I conflated two partitions. Let me redo cleanly at r ≥ 4.

**Setup:** s ∈ Z/p^r. Partition s by s mod p^{r-1} ∈ Z/p^{r-1}: each class has size p.

For each class indexed by ζ ∈ Z/p^{r-1}, the elements are s = ζ + p^{r-1}·v for v ∈ Z/p.

P_a(ζ + p^{r-1}·v) − P_a(ζ) computed above: at saddle ζ = s*(r) = (C_a − 1)/p mod p^{r-1}, this difference is ≡ 0 mod q (for r ≥ 3). At non-saddle ζ, dP_a/ds|_ζ ≢ 0 mod p^r, so the linear-in-v term is non-trivial mod q, and Σ_v e_q vanishes... wait, but v ranges over Z/p (size p), not Z/p^{r-1}. So I need:

Σ_{v ∈ Z/p} e_q(p^r · dP_a/ds|_ζ · v + ...) 

= Σ_v e_q(p^r · (unit · p^k) · v) where k = v_p(dP_a/ds|_ζ) for k < r.

For ζ at saddle: dP_a/ds|_ζ ≡ 0 mod p^r, so p^r · (this) is ≡ 0 mod p^{2r} which is ≡ 0 mod q iff 2r ≥ r+1, i.e., r ≥ 1. So at saddle, all p reps give identical phase mod q.

For ζ NOT at saddle: v_p(dP_a/ds|_ζ) ≤ r-1 (since ζ ≠ saddle to that precision). Then p^r · (unit · p^k) = unit · p^{r+k}. For this to give non-trivial e_q on Z/p, need r+k < r+1, i.e., k = 0. **But k could be > 0** for non-saddle ζ that are partial saddles to lower precision.

Hmm, this is more delicate. Let me reconsider.

**Restated partition:**

Partition s ∈ Z/p^r by s mod p (yielding classes of size p^{r-1}). For each class indexed by s_0 ∈ Z/p, elements are s = s_0 + p·u, u ∈ Z/p^{r-1}.

Compute P_a(s_0 + p·u) − P_a(s_0) mod q. The linear-in-u term: p · dP_a/ds|_{s_0} · u.

For Σ_u to be non-zero (saddle sub-sum): need dP_a/ds|_{s_0} ≡ 0 mod p^r.

For Σ_u to vanish: need dP_a/ds|_{s_0} ≢ 0 mod p^r. The orthogonality:
Σ_{u ∈ Z/p^{r-1}} e_q(p · K · u) = Σ_u e_{p^r}(K·u) — note q/p = p^r, so e_q(p·K·u) = e_{p^r}(K·u).

This sum is **p^{r-1}** if K ≡ 0 mod p^r, else 0... wait, the sum over u ∈ Z/p^{r-1} of e_{p^r}(K·u) is p^{r-1} if K ≡ 0 mod p^r/gcd(p^r, p^{r-1}) = p, else needs more care.

Standard identity: Σ_{u=0}^{n-1} e^{2πi·k·u/m} where n divides m. The sum is `(1 − ω^n)/(1 − ω)` where ω = e^{2πi·k/m}. For ω = 1 (i.e., k ≡ 0 mod m), sum = n. For ω^n = 1 but ω ≠ 1 (i.e., k ≡ 0 mod m/n, k ≢ 0 mod m), sum = 0. For ω^n ≠ 1, sum = (ω^n − 1)/(ω − 1) (non-zero, bounded).

In our case: n = p^{r-1}, m = p^r, k = K. ω = e^{2πi·K/p^r}. ω^n = e^{2πi·K/p}. 
- ω^n = 1 iff K ≡ 0 mod p.
- ω = 1 iff K ≡ 0 mod p^r.

So Σ_u e_{p^r}(K·u):
- = p^{r-1} if K ≡ 0 mod p^r (i.e., K = 0 in Z/p^r)
- = 0 if K ≡ 0 mod p but K ≢ 0 mod p^r (because then ω^n = 1, ω ≠ 1, geometric sum gives 0)
- = (ω^n − 1)/(ω − 1) ≠ 0 if K ≢ 0 mod p (small magnitude, bounded by p/|sin(π K/p^r)| etc.)

Hmm so if dP_a/ds|_{s_0} is divisible by p (so K = dP_a/ds|_{s_0} is divisible by p), the sub-sum vanishes. If dP_a/ds|_{s_0} is a UNIT mod p (not divisible by p), the sub-sum is non-zero but small.

So I need to identify ALL s_0 such that dP_a/ds|_{s_0} ≢ 0 mod p (these are the "non-saddle" residues that give small but non-zero sub-sums) versus those with dP_a/ds|_{s_0} ≡ 0 mod p (saddle to first order).

dP_a/ds = p(1 − C_a) + C_a·p²·s − ... mod p^r. Mod p (leading term mod p), dP_a/ds ≡ 0 mod p for ALL s (since dP_a/ds is divisible by p — every term has at least one factor of p). So dP_a/ds is ALWAYS divisible by p.

**Hmm that means K = dP_a/ds is always divisible by p, so the orthogonality sub-sum is either 0 or p^{r-1}.**

Specifically, write dP_a/ds = p · K' where K' = (1 − C_a) + C_a·p·s − C_a·p²·s² + ... .

The sub-sum over u of e_q(p · K · u) = e_q(p²·K'·u) = e_{p^{r-1}}(K'·u). This vanishes if K' ≢ 0 mod p^{r-1}, equals p^{r-1} if K' ≡ 0 mod p^{r-1}.

So the **first-order saddle condition** (at the p-adic level) is K' ≡ 0 mod p^{r-1}, i.e., dP_a/ds|_{s_0} ≡ 0 mod p^r.

K' = (1 − C_a) + C_a·p·s − C_a·p²·s² + C_a·p³·s³ − ...

Setting K' ≡ 0 mod p^{r-1}:
- Mod p: −p·t_1 + p·s ≡ 0 mod p, but everything is divisible by p, so mod p² (dividing by p): −t_1 + s ≡ 0 mod p ⟹ s ≡ t_1 mod p.

Wait, I think I'm confusing myself. Let me restart cleanly.

**Clean restart.** Write the formal-p-adic identity:
> dP_a/ds = p · (1 + ps − C_a) · (1 + ps)^{-1}    (formal in Z_p[[s]])

This is an EXACT identity in the formal p-adic ring, derived from dL_p(1+x)/dx = 1/(1+x).

For the truncated polynomial, the difference is O(p^{J_p+1}). At J_p ≥ r, the truncation matches modulo p^{r+1} = q.

**Saddle condition dP_a/ds ≡ 0 mod p^r (for orthogonality on Z/p^{r-1}-shift):**
> p · (1 + ps − C_a) · (1 + ps)^{-1} ≡ 0 mod p^r

Since (1 + ps)^{-1} is a unit in Z_p[[s]] (i.e., invertible), multiply through:
> p · (1 + ps − C_a) ≡ 0 mod p^r
> 1 + ps − C_a ≡ 0 mod p^{r-1}
> **ps ≡ C_a − 1 mod p^{r-1}**
> **s ≡ (C_a − 1)/p mod p^{r-2}**

Hmm — saddle determined mod p^{r-2}, not p^{r-1}. Let me double-check.

If 1 + ps − C_a ≡ 0 mod p^{r-1}, then ps ≡ C_a − 1 mod p^{r-1}. Dividing by p (valid since C_a − 1 ≡ 0 mod p):
s ≡ (C_a − 1)/p mod p^{r-1}/p = p^{r-2}.

So s* mod p^{r-2}. Number of saddle representatives in Z/p^r: p^r / p^{r-2} = p².

Hmm but at r=3, this gives s* mod p^1 = mod p, matching Phase 1. ✓

At r=4: s* mod p² (not mod p³).

## What about the "second saddle correction"?

We've found s ≡ s_0 mod p^{r-2}, with p² lifts in Z/p^r. Among these p² lifts, does the phase coherently add, or do they cancel?

Compute P_a(s_0 + p^{r-2}·w) − P_a(s_0) for w ∈ Z/p²:

P_a(s_0 + p^{r-2}·w) = P_a(s_0) + dP_a/ds|_{s_0}·(p^{r-2}·w) + (1/2)·d²P_a/ds²|_{s_0}·(p^{r-2}·w)² + ...

Term 1: dP_a/ds|_{s_0} · p^{r-2}·w. At saddle s_0 = s*, dP_a/ds|_{s*} ≡ 0 mod p^r (by saddle equation), so this term is ≡ 0 mod p^r · p^{r-2} = p^{2r-2}. For 2r-2 ≥ r+1 (i.e., r ≥ 3), this is ≡ 0 mod q. Term 1 vanishes for r ≥ 3.

Term 2: (1/2)·d²P_a/ds²|_{s*}·(p^{r-2})²·w² = (p^{2r-4}/2)·d²P_a/ds²|_{s*}·w². For r ≥ 3, 2r-4 ≥ 2. d²P_a/ds² is computed from the formal saddle: d/ds [p(1+ps−C_a)·(1+ps)^{-1}].

Compute d²P_a/ds² = d/ds [p − C_a · p/(1+ps)] = −C_a · d/ds[p/(1+ps)] = −C_a · p · (−p)/(1+ps)² = C_a · p² / (1+ps)².

At s = s*, 1 + ps* = C_a (exact in Z_p, by saddle equation), so (1+ps*)^{-2} = C_a^{-2}. So:

d²P_a/ds²|_{s*} = C_a · p² / C_a² = p² / C_a

**So d²P_a/ds²|_{s*} = p² / C_a in Z_p.**

Substituting:
Term 2 = (p^{2r-4}/2) · (p²/C_a) · w² = p^{2r-2}/(2·C_a) · w²

For r=3: p^4/(2 C_a) · w² mod q = p^4. So Term 2 = p^4·w²/(2 C_a) ≡ 0 mod p^4 = q. So vanishes at r=3.

For r=4: p^6/(2 C_a) · w². q = p^5. p^6/(2 C_a) · w² mod p^5 = p^6 mod p^5 = 0. So vanishes at r=4 too. ✓

In general: Term 2 = p^{2r-2}/(2 C_a) · w² mod p^{r+1}. 2r-2 ≥ r+1 iff r ≥ 3. So Term 2 vanishes at r ≥ 3.

**Higher order terms:** Term k = (1/k!) · d^k P_a/ds^k|_{s*} · (p^{r-2}·w)^k. Each d^k P_a/ds^k|_{s*} has p-adic valuation ≥ k (from C_a · p^k structure of the truncated log derivatives). So Term k ~ p^k · (p^{r-2})^k = p^{k(r-1)}. For k ≥ 1 and r ≥ 3, k(r-1) ≥ r-1; for k ≥ 2 and r ≥ 3, k(r-1) ≥ 2(r-1) = 2r-2 ≥ r+1 (since r ≥ 3). So all Term k for k ≥ 2 vanish mod q at r ≥ 3. Term 1 vanishes too (as computed).

**So all p² saddle representatives at r ≥ 3 give P_a(s) ≡ P_a(s*) mod q — coherent addition.**

Wait, this gives me MORE saddle representatives than I had at Phase 1 (where I claimed s ≡ s* mod p gave p^{r-1} reps with coherent phase). Let me reconcile.

At r=3: saddle mod p^{r-2} = p^1 = p. Number of reps in Z/p^r = Z/p^3: p². So my Phase 1 derivation matches.

At r=4: saddle mod p^{r-2} = p². Number of reps in Z/p^4: p². Same number of reps as at r=3.

Hmm but then the magnitude calculation:
|G_p(a)| = p² · |saddle phase coherent factor|

But empirical (FHAT) says |G_p(a)| = p^{(r+1)/2}. At r=4: |G_p| = p^{2.5} = p² · √p. So magnitude p² alone is OFF BY √p at r=4.

**This is where the Gauss-sum factor enters at r ≥ 4 (for even r in particular).** The "p² coherent saddle reps" calculation misses a √p factor.

Let me reconsider what could give √p.

## Where the missing √p comes from at even r

Going back to careful Hensel-saddle analysis at r=4 (q=p^5, period=p^4):

Partition s ∈ Z/p^4 by s mod p^2 (= s mod p^{r-2}), with p² classes, each of size p² (the p² lifts within Z/p^4).

Among the p² classes (indexed by ζ ∈ Z/p²):
- One class is the saddle, where K' = (1 − C_a + C_a·p·ζ − ...) ≡ 0 mod p^{r-1} = p³.
- The other p² − 1 classes are non-saddle.

But what about classes where the saddle condition is met to ORDER p² (i.e., K' ≡ 0 mod p²) but NOT to order p^3? These are "partial saddles" — first-order satisfied but second-order missed.

For non-saddle class ζ: write K'|_ζ = p^k · (unit) with k < r-1.
- The shift step is p^{r-2}, so the sub-sum over the p² lifts within the class is Σ_w e_q((p^{r-2})·(p·K')·w + quadratic) = Σ_w e_q(p^{r-1}·K'·w + ...) = Σ_w e_{p²}(K'·w + ...).

Hmm getting messy. Let me try a cleaner approach.

## Cleaner approach: full saddle analysis via Cochrane Proposition

**Approach (structural):** rely on Cochrane's exponential sum theorem for the polynomial P_a.

The standard p-adic stationary-phase result (Cochrane, also Igusa, also Denef-Loeser style):

For an additive character sum Σ_s e_q(f(s)) where f is a polynomial mod q = p^M, the value is determined by:
- Hensel-lifting the critical points (where df/ds ≡ 0 mod p^M up to appropriate height)
- Computing the leading Taylor expansion at the lifted critical point
- Multiplying by a Gauss-sum factor at each critical point depending on the multiplicity / order of vanishing of df/ds

In our case, f(s) = P_a(s), df/ds = p · (1 + ps − C_a)/(1 + ps). The critical point is at s* = (C_a − 1)/p, a SIMPLE critical point in the formal-p-adic sense (since d²f/ds²|_{s*} = p²/C_a, a unit times p², non-zero).

For a simple critical point, the standard Cochrane-style result is:

> **Σ_{s ∈ Z/p^M} e_{p^M}(f(s)) = p^{(M - v_p(f''(s*)))/2} · e_{p^M}(f(s*)) · η**

where η is a quadratic Gauss-sum factor depending on whether the Hessian f''(s*) (here just second derivative) is a quadratic residue mod p, and on the parity of M − v_p(f'') (Gauss sum of length p multiplied in for odd parity).

In our case:
- M = r+1
- f''(s*) = p²/C_a, v_p = 2
- M − v_p(f'') = r+1 − 2 = r − 1

So:
> G_p(a) = p^{(r-1)/2} · e_q(P_a(s*)) · η_p

where η_p depends on parity of r − 1 (equivalently, on parity of r):
- r odd: r−1 even, p^{(r-1)/2} is integer, η_p incorporates a quadratic-character factor (Legendre symbol).
- r even: r−1 odd, p^{(r-1)/2} is half-integer power, η_p incorporates a √p quadratic Gauss sum.

**Compare to empirical magnitude:** |G_p(a)| = p^{(r+1)/2} = √q.

The standard Cochrane formula gives |G_p(a)| = p^{(r-1)/2} · |η_p|.

For these to match: |η_p| = p^{(r+1)/2 − (r-1)/2} = p^1 = p.

Hmm, so |η_p| = p. Let me reconsider; maybe I mis-stated the Cochrane formula.

The standard p-adic stationary-phase formula (Igusa, Cochrane book): for f with simple critical point s* with v_p(f''(s*)) = k:
> Σ_s e_{p^M}(f(s)) = p^{(M+k)/2} · (Gauss-sum-factor) · e_{p^M}(f(s*))
   
   Hmm let me try harder to recall the correct exponent. Actually the cleanest statement: for f(s) = a + b(s−s*)² · (unit) + higher, mod p^M, the sum is p^{M/2} · (sign + character) · e_{p^M}(a). But with v_p(b) = k > 0, it's p^{(M+k)/2}.

Let me derive from first principles. Near s*:

f(s* + h) = f(s*) + (1/2) f''(s*) · h² + O(h³)
            = f(s*) + (p²/(2 C_a)) · h² + O(h³)

For h ∈ Z/p^r (since s ∈ Z/p^r), the sum:
Σ_h e_{p^{r+1}}(f(s* + h)) = e_{p^{r+1}}(f(s*)) · Σ_h e_{p^{r+1}}((p²/(2C_a)) · h²)

The h³ and higher terms: cubic term has v_p ≥ 3 in f''' (from log expansion d³L/ds³ at s* gives p³ · -2/(1+ps*)^3 = -2p³/C_a³). Times h³: v_p(p³·h³) ≥ 3 + 0 = 3 for unit h, but h ∈ Z/p^r so h³ can be anything in Z/p^r. The cubic term mod p^{r+1} is p³ · h³ · (unit)/6 mod p^{r+1}, which is non-trivial unless v_p ≥ r+1, requiring 3 + 3·v_p(h) ≥ r+1. For h a unit, 3 ≥ r+1 requires r ≤ 2.

So for r ≥ 3, the cubic+ terms DO contribute mod q at general h, complicating the "Gauss sum" reduction.

This complication is exactly what the Hensel-lifting machinery is supposed to handle. **The Hensel-corrected saddle s*(r) "absorbs" the contributions from higher-order terms.**

## The Hensel correction interpretation

Standard Hensel-Newton iteration: at the formal-p-adic saddle s* = (C_a − 1)/p, the saddle is exact in the FORMAL p-adic sense. Truncating P_a to a polynomial mod q means the saddle equation 1 + ps − C_a = 0 is satisfied EXACTLY in Z_p, and the truncated version is satisfied to high precision in Z/p^r.

**Claim:** the Hensel-lifted saddle is s*(r) = (C_a − 1)/p mod p^{r-1} (literally the integer-division of (C_a − 1) by p, with C_a's full r-digit representation).

**Why:** the formal p-adic root of `1 + ps − C_a = 0` is `s = (C_a − 1)/p`, and this is already in Z_p (well-defined since C_a − 1 ∈ p·Z_p by C_a ≡ 1 mod p). The truncated polynomial P_a's "true" saddle is exactly this, reduced mod p^{r-1}.

**Importantly: the saddle s*(r) is a CLEAN closed-form digit extraction, NOT a Hensel-iteration series.**

This is the structural surprise. The standard "Hensel-lift" framing would compute s*(r) = s*(3) + p²·δ_1 + p³·δ_2 + ... where δ_k solve the iteration at each order. But because P_a's saddle equation FACTORS so cleanly (`(1+ps − C_a)·(1+ps)^{-1}` has a linear factor in s), the iteration terminates trivially at the first step: s = (C_a − 1)/p is exact, no iteration needed.

**The "Hensel series" is identically s*(r) = (C_a − 1)/p mod p^{r-1}, no δ_k corrections.**

## Verification: at r=4, does the magnitude work?

Now at r=4 with s*(4) = (C_a − 1)/p mod p^3:
- We have p² saddle representatives in Z/p^4 (mod p^{r-2} = p²).
- These all add coherently (from the analysis above showing higher-order terms vanish mod q for r ≥ 3).

Wait, but I just realized I need to re-examine the saddle multiplicity argument. At r=4:
- Saddle determined mod p^{r-2} = p² (from saddle equation `1 + ps − C_a ≡ 0 mod p^{r-1} = p³`).
- Number of saddle representatives in Z/p^4 = p^4 / p² = p².
- Phase across p² reps: by the d²f/ds² Taylor analysis, the variation is (p²/(2C_a))·(p^{r-2}·w)² = (p²/(2C_a))·p^{2r-4}·w² = p^{2r-2}/(2C_a)·w². For r=4: p^6/(2C_a)·w² mod p^5 = 0 (since 6 > 5). So phase is exactly constant on the saddle class.

**So at r=4, p² coherent saddle reps give magnitude p². But empirical says √q = p^{5/2} = p²·√p. Missing factor √p!**

Let me reconsider. Maybe I'm missing non-saddle classes that contribute partially.

**Non-saddle classes:** ζ such that K'(ζ) = 1 + p·ζ − C_a is NOT ≡ 0 mod p^{r-1} = p³. But there are "partial saddles" where K'(ζ) ≡ 0 mod p² (not mod p³). Each ζ where v_p(K'(ζ)) = k gives a sub-sum with magnitude p^k.

Count of ζ ∈ Z/p² with v_p(K'(ζ)) = k:
- K'(ζ) = (1 − C_a) + p·C_a·ζ + O(p²) mod p² (treating ζ ∈ Z/p²; here K'(ζ) ∈ Z/p²).
- v_p(K'(ζ)) = k means K'(ζ) ≡ 0 mod p^k but not mod p^{k+1}.

For k = 0 (K' a unit): ζ such that 1 − C_a + p C_a ζ ≢ 0 mod p, i.e., 1 − 1 = 0 mod p plus p stuff which is 0 mod p — so K'(ζ) ≡ 0 mod p ALWAYS (for C_a ≡ 1 mod p). So k ≥ 1 always.

For k = 1: K'(ζ) ≡ 0 mod p (always), not ≡ 0 mod p². Mod p²: K'(ζ) = (1 − C_a) + p·C_a·ζ. Setting ≡ 0 mod p²: this is one congruence in ζ mod p (since dividing by p, get (1−C_a)/p + C_a·ζ ≡ 0 mod p, which is t_1 ≡ C_a · ζ mod p, i.e., ζ ≡ t_1/C_a ≡ t_1 mod p). So exactly ONE residue of ζ mod p satisfies K'(ζ) ≡ 0 mod p². 

Hmm but ζ ranges over Z/p², so given ζ mod p, there are p lifts. So among Z/p² ζ's, exactly p of them have K'(ζ) ≡ 0 mod p², and these are the saddles to higher precision.

Now within these p second-order saddles, do they all give K' ≡ 0 mod p³? Need to check the next-order:

K'(ζ) mod p³ for ζ ≡ t_1 mod p, ζ = t_1 + p·ζ_1: K'(ζ) = (1 − C_a) + p·C_a·(t_1 + p·ζ_1) + O(p²) where O(p²) comes from the quadratic term in K's full expansion (the −p²·s²·C_a term in dP_a/ds = p·K'·...).

Actually, recall K' = (1 + ps − C_a)/(1+ps) in formal Z_p. So K'(ζ) ≡ 0 mod p^k iff (1 + pζ − C_a) ≡ 0 mod p^k. Setting 1 + p·ζ ≡ C_a mod p^k means ζ ≡ (C_a − 1)/p mod p^{k-1}.

So:
- v_p(K'(ζ)) ≥ k ⟺ ζ ≡ (C_a − 1)/p mod p^{k-1}
- The saddle to order p^k corresponds to ζ matching s* to k-1 digits.

In Z/p^{r-1} (where we want K' ≡ 0 mod p^{r-1} for the FULL saddle):
- Setting v_p(K') ≥ r-1 requires ζ ≡ (C_a − 1)/p mod p^{r-2}. Number of such ζ in Z/p^{r-1}: p (one per lift). Hmm wait Z/p^{r-1} elements ζ; ζ ≡ s* mod p^{r-2} has p representatives in Z/p^{r-1}.

Hmm so my count of "p² coherent saddle reps" was for the partition by s mod p^{r-2}, giving p² classes each with p reps. Let me redo with the partition by s mod p (giving p classes each with p^{r-1} reps).

Hmm this is getting confused. Let me try a different tactic.

## Direct Plancherel argument

Empirical: |G_p(a)| = √q = p^{(r+1)/2}. This is FHAT-verified. The empirical fact is established.

From R78.3 / Theorem 78.3 (and family-level analog): |F̂(ξ)| = q/√(|supp|) on the support of size q/p². But |supp| corresponds to {a ≡ 1 mod p in Z/p^r}, of size p^{r-1}. The original F̂ is the period-p^r DFT-transformed: |F̂(ξ)|² summed over support = q² · ... 

OK let me re-examine R78.3 logic, applied family-level:
- f(u) = e_q(c·(1+p)^u) is p^r-periodic in u.
- F̂(ξ) = Σ_{u=0}^{q-1} e_q(c·(1+p)^u − ξu) — sum over q = p^{r+1} values of u, periodic with period p^r, so F̂ = p · G where G is the period sum.
- |F̂|² summed over ξ in Z/q = q · Σ_u |f|² = q² (Plancherel).
- F̂ supported on |supp(F̂)| = p^{r-1} elements (the principal-unit sub-support).
- So |F̂(ξ)|² · |supp| = q² ⟹ |F̂(ξ)| = q/√|supp| = q · p^{-(r-1)/2} = p^{r+1} · p^{-(r-1)/2} = p^{(r+3)/2}.
- |G_p(a)| = |F̂(p·a)|/p = p^{(r+3)/2}/p = p^{(r+1)/2} = √q. ✓

**This Plancherel + uniform-magnitude argument gives |G_p(a)| = √q independent of saddle analysis.** The empirical magnitude is structurally forced.

## So where does √q come from in the saddle picture?

The saddle gives |G_p(a)| = (number of saddle reps) · (variation factor). If saddle has p² coherent reps (size p^{r-2} mod p^r partition gives p² classes, each with p² lifts — wait I keep getting confused).

Let me carefully recount.

**At r=3:** Partition s ∈ Z/p^3 by s mod p (p classes, each of size p²). Phase 1 derivation: only ONE class (s ≡ t_1 mod p) gives coherent sum; other classes vanish. Coherent magnitude = p² (class size). Empirical |G_p| = p² = √q (since at r=3, √q = p^{(r+1)/2} = p²). **Match: saddle class size = √q exactly at r=3.**

**At r=4:** Empirical √q = p^{5/2} = p²·√p. The saddle analysis at r=4 must give a coherent contribution of magnitude p² · √p, not just p².

**Possibilities for where the √p comes from at r=4:**

(a) Multiple saddle classes — instead of ONE saddle class, there are √p saddle classes each contributing the same magnitude (then total = √p · p² = p² · √p). But this would require ~√p ≈ p^{0.5} classes, which is non-integer — doesn't quite make sense.

(b) A Gauss-sum factor of magnitude √p from incomplete saddle (quadratic completion of square at the saddle, mod p — classical Gauss sum has magnitude √p).

(c) Combination of saddle multiplicity (extra factor p) and √p Gauss sum (e.g., p · √p · √p = p · p = p², adjusting other factors).

The Gauss-sum interpretation (b) is the classical mechanism. Let me think about how it arises.

At r=4, the saddle equation is 1 + ps − C_a ≡ 0 mod p^{r-1} = p^3. The saddle "determines s mod p^{r-2} = p²", with p² lifts to Z/p^4 (p² coherent reps).

BUT the quadratic-correction term `(p²/(2C_a))·h²` where h = s − s* gives, at h = p^{r-2}·w = p²·w (w ∈ Z/p²), the quadratic phase `p²/(2C_a) · p^4 · w² = p^6/(2C_a) · w²`. Mod q = p^5: p^6 ≡ 0 mod p^5. So quadratic correction VANISHES exactly mod q at r=4.

But wait — maybe I should partition more finely. What if I partition by s mod p^{r-1} = p^3 instead? Each class has p elements. There are p^3 classes. Among these p^3, the saddle is at s ≡ s* mod p^{r-1} — but wait, s* is only defined mod p^{r-2} = p². So "saddle mod p^{r-1}" has p saddles (the lifts of s* mod p^{r-2} into Z/p^{r-1}).

For each saddle ζ ∈ Z/p^{r-1} = Z/p^3 (p saddles, each a lift of s*), partition lifts s = ζ + p^{r-1}·v for v ∈ Z/p (each class size p).

P_a(ζ + p^{r-1}·v) − P_a(ζ): leading term p · dP_a/ds|_ζ · p^{r-1}·v = p^r · dP_a/ds|_ζ · v. At ζ a saddle (K'(ζ) ≡ 0 mod p^{r-1}, so dP_a/ds = p·K' has v_p ≥ r), p^r · dP_a/ds|_ζ has v_p ≥ 2r ≥ r+1 (for r ≥ 1). So linear-in-v term vanishes mod q.

Quadratic in v: (1/2)·d²P_a/ds²|_ζ · (p^{r-1})² · v² = (p²/(2C_a)) · p^{2r-2} · v² = p^{2r}/(2C_a) · v². For r=4: p^8/(2C_a) · v² mod p^5 = 0. So quadratic also vanishes.

So at the p saddles in Z/p^{r-1}, the lifts in Z/p^r all coherently add. Total saddle contribution to G_p(a): p · p = p² (p saddles × p lifts each).

Same as before. Still missing √p.

**OK so the standard "Hensel-saddle p² contribution" gives magnitude p² at r=4, but empirical says p² · √p.** There must be something else.

## Reconsidering — non-saddle contributions

What about ζ NOT at saddle (where K'(ζ) has v_p = k < r-1)? The sub-sum over the p lifts of ζ in Z/p^r gives:

P_a(ζ + p^{r-1}·v) − P_a(ζ) = p^r · dP_a/ds|_ζ · v + O(higher).

Term 1 contribution to e_q: e_q(p^r · dP_a/ds|_ζ · v). With dP_a/ds = p·K' and v_p(K'(ζ)) = k, dP_a/ds|_ζ = p · p^k · (unit) = p^{k+1} · (unit). So p^r · dP_a/ds|_ζ = p^{r+k+1} · (unit). Mod q = p^{r+1}: this is p^{k} · (unit) (after factoring out p^{r+1}).

Wait, e_q(p^{r+k+1} · X) = e_{p^{r+1−(r+k+1)}}(X) = e_{p^{−k}}(X) — but that's a "negative-modulus" character, which makes no sense (or equivalently, the phase is trivial when r+k+1 ≥ r+1, i.e., k ≥ 0, always).

So actually for any ζ (saddle or not), the linear-in-v term is ≡ 0 mod q because p^r · (anything in Z_p) is at least p^r which gets absorbed.

Hmm, then by my analysis the **EVERY partition class** has all p lifts coherently adding, giving Σ_{ζ ∈ Z/p^{r-1}} e_q(P_a(ζ)) · p.

So G_p(a) = p · Σ_{ζ ∈ Z/p^{r-1}} e_q(P_a(ζ))? That's just splitting Z/p^r into Z/p^{r-1} (the "ζ" variable) × Z/p (the "v" variable), with v contributing trivially.

But Z/p^{r-1} has p^{r-1} elements, not p² (only at r=3 is p^{r-1} = p²). Hmm let me check.

Oh wait — at r=4, p^{r-1} = p^3, not p². So Z/p^{r-1} = Z/p^3, and "saddle class size" of p^{r-1} = p^3 corresponds to total coherent magnitude p^3, not p² as I miscomputed.

Hmm let me redo from scratch. The partition I want:
- Partition s ∈ Z/p^r by s mod p (p classes, each of size p^{r-1}).
- For each class (indexed by s_0 ∈ Z/p), sum over u ∈ Z/p^{r-1} of e_q(P_a(s_0 + p·u)).

The shift-step is p, so the orthogonality argument from Phase 1 should still apply.

**Linear-in-u term in P_a(s_0 + p·u) − P_a(s_0): p · dP_a/ds|_{s_0} · u = p² · K'(s_0) · u mod q.**

K'(s_0) ∈ Z_p has v_p(K'(s_0)) = v_p(1 + p·s_0 − C_a). For s_0 ∈ Z/p:
- 1 + p·s_0 = 1 + p·s_0 (in Z_p).
- C_a ∈ Z/p^r, lift to Z_p.
- 1 + p·s_0 − C_a: leading term 1 − C_a = −p·t_1 − p²·t_2 − ... mod p^r. So 1 + p·s_0 − C_a = p(s_0 − t_1) − p²·t_2 − ... mod p^r.

So K'(s_0) = p(s_0 − t_1) − p²·t_2 − p³·t_3 − ... mod p^r (treating s_0 ∈ Z/p so v_p(s_0) ≥ 0 typically).

v_p(K'(s_0)):
- If s_0 ≠ t_1 mod p (i.e., s_0 − t_1 is a unit): K'(s_0) has v_p = 1 (leading p·(s_0 − t_1)).
- If s_0 = t_1 mod p: K'(s_0) = −p²·t_2 − p³·t_3 − ... has v_p ≥ 2.

Substitute back into the sub-sum over u ∈ Z/p^{r-1}:

p² · K'(s_0) · u mod q. K'(s_0) = p^{v_p(K')} · (unit). So p² · K'(s_0) = p^{2 + v_p(K')} · (unit) · u mod q = p^{r+1}.

For s_0 ≠ t_1: 2 + 1 = 3, so p² K' = p^3 · (unit). For r ≥ 3, the linear-in-u term is non-trivial mod q (i.e., contributes a character of Z/p^{r-1} that's non-trivial). Σ_u vanishes by orthogonality.

For s_0 = t_1: 2 + 2 = 4. For r ≥ 4: p^4 · u mod p^{r+1}, with r+1 ≥ 5, so this is p^4·u·(unit) which is non-trivial mod p^5 (linear character of period p mod p^5). Σ_u over u ∈ Z/p^{r-1} = Z/p^3 of e_{p^5}(p^4·(unit)·u) = Σ_u e_p((unit)·u) summed over Z/p^3 = p^2 · Σ_{u in Z/p} e_p((unit)·u) = p² · 0 = 0.

**So at r=4, s_0 = t_1 sub-sum ALSO vanishes? That gives G_p(a) = 0, contradicting |G_p(a)| = √q.**

Something is wrong. Let me re-examine.

I missed the **quadratic-in-u** correction. At s_0 = t_1, the linear term vanishes by orthogonality, but the quadratic term can contribute.

P_a(s_0 + p·u) − P_a(s_0) = p · K'(s_0) · p · u + (1/2) · (p² K''(s_0)) · (p·u)² + ...

Wait, dP_a/ds = p · K'(s)·(1+ps)^{-1} from the formal identity. And K'(s) = 1 + ps − C_a (linear in s). And dK'/ds = p (constant). So d²P_a/ds² = d/ds [p · K'(s)/(1+ps)] = p · [p/(1+ps) − K'(s)·p/(1+ps)²] = p²/(1+ps) · [1 − K'(s)/(1+ps)] = p²/(1+ps) · [(1+ps − K'(s))/(1+ps)] = p²·(1+ps − (1+ps−C_a))/(1+ps)² = p²·C_a/(1+ps)².

At s_0 = t_1 (mod p, so 1 + p·t_1 = first-order-saddle), 1 + ps_0 = 1 + p·t_1. In Z_p, this is approximately C_a (close to C_a mod p², specifically C_a mod p² = 1 + p·t_1).

d²P_a/ds²|_{t_1} = p² · C_a / (1 + p·t_1)² ≡ p²·C_a/C_a² mod higher = p²/C_a mod higher.

So second-derivative at first-order saddle ≈ p²/C_a, v_p = 2 (unit divided into p²).

Quadratic term in P_a expansion: (1/2)·d²P/ds²|_{t_1}·(p·u)² = (1/2)·(p²/C_a)·p²·u² = p^4·u²/(2 C_a).

Mod q = p^5 at r=4: p^4·u²/(2 C_a) mod p^5 = p^4 · (u²·(2C_a)^{-1}) mod p^5. This is p^4 · (a unit) · u² mod p^5.

So the sub-sum becomes:
Σ_{u ∈ Z/p^3} e_{p^5}( ... + p^4 · (2C_a)^{-1} · u² )

Now P_a's linear-in-u contribution at s_0 = t_1 was p²·K'(t_1)·u with K'(t_1) = −p²·t_2 + O(p³), giving p²·(−p²·t_2 + ...) · u = −p^4·t_2·u + O(p^5·u). Mod p^5: −p^4·t_2·u.

So total: e_{p^5}(P_a(t_1) − p^4·t_2·u + p^4·(2C_a)^{-1}·u²)
       = e_{p^5}(P_a(t_1)) · e_{p^5}(p^4·(−t_2·u + (2C_a)^{-1}·u²))
       = e_{p^5}(P_a(t_1)) · e_p(−t_2·u + (2C_a)^{-1}·u²)

Sum over u ∈ Z/p^3:
Σ_u e_p(−t_2·u + (2C_a)^{-1}·u²) — but e_p depends on u mod p only, so this is p² · Σ_{u in Z/p} e_p(−t_2·u + (2C_a)^{-1}·u²) (each value of u mod p appears p² times in Z/p^3).

Σ_{u ∈ Z/p} e_p(−t_2·u + (2C_a)^{-1}·u²) = Σ_u e_p(quadratic in u).

This is a classical **quadratic Gauss sum** of length p! Its magnitude is √p (when the leading quadratic coefficient is a unit mod p, which (2C_a)^{-1} is for p ≥ 3).

So the sub-sum at s_0 = t_1 contributes:
p² · e_{p^5}(P_a(t_1)) · G_{quad}(p, ...) where G_{quad} has |G_{quad}| = √p.

Magnitude: p² · √p = p^{5/2} = √q at r=4. ✓ **Match.**

**So at r=4, the √p Gauss-sum factor emerges from the quadratic-correction term in the Taylor expansion at the first-order saddle.**

## Generalization

At general r, partition by s mod p:
- s_0 ≠ t_1: sub-sum vanishes (orthogonality on linear term, which has v_p < r+1).
- s_0 = t_1: sub-sum contributes (e_q(P_a(t_1))) · (quadratic-in-u Gauss sum).

The structure of the quadratic-in-u Gauss sum depends on r:

The full expansion of P_a(t_1 + p·u) in u has:
- Constant term: P_a(t_1) (which is the "phase part" we computed, modulo Hensel correction).
- Linear in u: K'(t_1)·p²·u — non-trivial for r ≥ 4 (gives the e_p(−t_2·u) etc., which combines with quadratic to give shifted Gauss sum).
- Quadratic in u: p^4·u²/(2 C_a).
- Higher terms: cubic and beyond, kicking in at higher r.

Looking at the pattern: at r=4, quadratic gives √p. At r=5: cubic gives ... ?

Actually wait, at r=5 the quadratic term `p^4·u²/(2 C_a)` is still at level p^4, but now mod q = p^6, so the quadratic term is at level p^{4}/p^6 = level −2 (well, it's level 4 below 6, so contributes non-trivially mod p^6 as p^4·(unit)·u² mod p^6 = p^4·(unit) mod p^6 · u² which is a length-p² character on u² ... getting complicated).

Let me try a different organization.

**Parity-dependent structure:**

The standard p-adic stationary-phase result says:
> Σ_{s ∈ Z/p^M} e_{p^M}(f(s)) = p^{(M+k)/2} · e_{p^M}(f(s*)) · κ

where k = v_p(f''(s*)) and κ is a Gauss-sum factor. The exponent (M+k)/2 must be an integer; if not, the Gauss sum κ contributes √p to absorb the half-integer.

In our case M = r+1, k = v_p(d²P_a/ds²|_{s*}) = 2.

(M+k)/2 = (r+1+2)/2 = (r+3)/2.

So expected magnitude p^{(r+3)/2}? But empirical |G_p(a)| = p^{(r+1)/2}. Off by p.

Hmm. Let me re-derive. The classical formula for the p-adic complete Gauss sum:

Σ_{s ∈ Z/p^M} e_{p^M}(f(s)) where f has a simple zero of df/ds at s*, with v_p(f''(s*)) = k:

For M > 2k+1: Σ = p^{(M+k)/2} · e_{p^M}(f(s*)) · ε_p,k,M

where ε_p,k,M is a unit Gauss-sum factor. (Reference: Igusa, "Forms of higher degree", or Cochrane's textbook.)

In our case: M = r+1, k = 2. Need r+1 > 5, i.e., r ≥ 5, for the formula. For r=4 boundary, separate analysis.

For r ≥ 5: |Σ| = p^{(r+3)/2} · |ε| = p^{(r+3)/2} (since |ε| = 1).

But empirically |G_p(a)| = p^{(r+1)/2}, not p^{(r+3)/2}. **Off by factor p.**

Hmm. Let me reconsider. Maybe d²f/ds²|_{s*} actually has v_p ≠ 2.

d²P_a/ds²|_{s*} computed: = p²/C_a in the formal-p-adic sense, with v_p = 2.

But wait — let me check the formula d²P_a/ds² = p²·C_a/(1+ps)². At s = s*, 1+ps* = C_a (exact), so d²P_a/ds²|_{s*} = p² · C_a / C_a² = p²/C_a. v_p = 2.

So in the Igusa formula, k = 2.

Hmm let me recheck the Igusa formula. Actually, in many references the formula reads:

For f(s) = b(s − s*)² + higher (with simple critical point, b = (1/2) f''(s*)):
Σ_{s ∈ Z/p^M} e_{p^M}(f(s)) = ε · p^{M/2 + v_p(b)/...}

Hmm I think I'm mis-remembering. Let me derive.

After the change s = s* + h:
f(s* + h) = f(s*) + b·h² + cubic + ...    where b = (1/2)f''(s*) = p²/(2C_a) at our case.

For h ∈ Z/p^r: Σ_h e_{p^{r+1}}(b·h² + cubic·h³ + ...).

For r large enough that cubic dominates and beyond are negligible, this is just Σ_h e_q(b·h²) = Σ_h e_q(p²·h²/(2C_a)).

Let h = p^j · h' with v_p(h) = j. As j increases, the contribution to e_q(p²·h²) becomes trivial. Specifically:
- For v_p(h) = 0 (h a unit): p²·h²/(2C_a) has v_p = 2, contributing e_{p^{r-1}}(h²/(2C_a)).
- Sum over h ∈ Z/p^r with h a unit: ... messy.

Better: factor h = p^j · h_0 with h_0 a unit, then sum over j ≥ 0 and h_0.

Actually the cleanest: change variables h → h such that h ranges over Z/p^r. Let h' = h mod p^? for some appropriate truncation. After absorbing the p² factor, the sum becomes:

Σ_{h ∈ Z/p^r} e_{p^{r+1}}(p²·h²·(2C_a)^{-1}) = Σ_h e_{p^{r-1}}(h²·(2C_a)^{-1})

This is the sum of e_{p^{r-1}} over h ∈ Z/p^r. But h ∈ Z/p^r vs the character on Z/p^{r-1}: each value of h mod p^{r-1} appears p times. So:

Σ_{h ∈ Z/p^r} e_{p^{r-1}}(h²/(2C_a)) = p · Σ_{h' ∈ Z/p^{r-1}} e_{p^{r-1}}((h')²/(2C_a))

(where the squaring is mod p^{r-1}; actually no, h² mod p^{r-1} for h ∈ Z/p^r — need to think...)

Hmm, h ∈ Z/p^r, h² ∈ Z/p^r (well, h² is determined by h mod p^? — specifically h² is determined by h mod p^{r-1} when r ≥ 1; since (h + p^{r-1}·k)² = h² + 2·h·p^{r-1}·k + p^{2r-2}·k², which mod p^{r-1} ≡ h² for r-1 ≤ r-1 ✓).

Wait, I want h² mod p^{r-1} for the inner exponential. If h ∈ Z/p^r, h² ∈ Z/p^r in general. We're using e_{p^{r-1}}(h²/(2C_a)), so only h² mod p^{r-1} matters. And h² mod p^{r-1} is determined by h mod p^{(r-1+1)/2}... hmm no.

Let h = h' + p^j · k with h' ∈ Z/p^j and k ∈ Z/p^{r-j}. h² = h'² + 2·h'·p^j·k + p^{2j}·k². Mod p^{r-1}: if 2j ≥ r-1 (i.e., j ≥ (r-1)/2), then h² mod p^{r-1} = h'² + 2h'·p^j·k mod p^{r-1}.

OK this is getting deep. Let me just claim the standard result and check magnitudes.

**Standard p-adic Gauss sum (for simple Gaussian):**

Σ_{h ∈ Z/p^N} e_{p^N}(α · h²)  =  ε(α, N, p) · p^{N/2}    (for α a unit mod p, p ≥ 3)

where ε is a unit (specifically a Legendre symbol or quartic root of unity depending on parity).

For our case: rewrite p²·h²/(2C_a) mod p^{r+1} as α·h² where α = p²/(2C_a). v_p(α) = 2.

Let h = p^j · h_0 with h_0 ranging over units (or h_0 ∈ Z/p^{r-j} more carefully). Actually, decompose:

Σ_{h ∈ Z/p^r} e_q(α h²) where α = p²·u with u a unit.

Substitute h: h = p^j · h_1 with j = v_p(h). Then α·h² = p²·u·p^{2j}·h_1² = p^{2+2j}·u·h_1². For this to contribute non-trivially mod q = p^{r+1}: need 2+2j < r+1, i.e., j < (r-1)/2.

For j ≥ (r-1)/2: α·h² ≡ 0 mod q, contributes 1 to the sum. Number of such h: ... let me just count for r=4 specifically.

r=4: q = p^5, α = p²·u. Σ_{h ∈ Z/p^4} e_{p^5}(α·h²):
- h = 0: contributes 1.
- h = p·h_1, h_1 ∈ (Z/p^3)*: α·h² = p² · u · p² · h_1² = p^4·u·h_1² mod p^5 → e_p(u·h_1²). Sum over h_1 ∈ (Z/p^3)* of e_p(u·h_1²): each value of h_1 mod p in (Z/p)* appears p²·(p−1)/(p−1)·... well there are (p-1)·p² units in Z/p^3 (units = coprime to p), so each h_1 mod p in (Z/p)* (i.e., (p-1) values) appears p² times. Sum = p² · Σ_{h_1 mod p ∈ (Z/p)*} e_p(u·h_1²). The last sum is the classical Gauss sum minus the h_1=0 contribution: Σ_{h_1=0}^{p-1} e_p(u·h_1²) − 1 = G(u, p) − 1 where G(u, p) = √p · (Legendre symbol stuff).

Wait this is getting messy. Let me just trust the structural result.

**Standard result (Igusa's local zeta function approach):**

For Σ_{s ∈ Z/p^M} e_{p^M}(f(s)) with f having a simple critical point s* (i.e., f'(s*) = 0, f''(s*) ≠ 0 mod p^... something), the sum equals:

p^{(M+k)/2} · sign factor · e_{p^M}(f(s*))    where k = v_p(f''(s*)/2) appropriately.

I'm not nailing down the exact formula but the structural shape is:
- exponent in p is (M + k)/2 (approximately)
- a Gauss-sum factor of magnitude 1 if (M+k) is even, magnitude √p if (M+k) is odd (this introduces the parity dependence).

For our case M = r+1, k = ?. If k = 1 (instead of 2 as I computed), then (M+k)/2 = (r+2)/2:
- r even: (r+2)/2 = r/2 + 1 — half-integer-issue? r even means r/2 integer, (r+2)/2 = r/2 + 1 integer.
- r odd: (r+2)/2 half-integer; Gauss sum gives √p.

Magnitude: p^{(M+k)/2 + (Gauss-correction)} = p^{(r+2)/2 + (0 or 1/2)}. Empirical p^{(r+1)/2}. Match if (r+2)/2 + correction = (r+1)/2 + integer.

(r+2)/2 − (r+1)/2 = 1/2. So if correction = −1/2, match. Gauss sum contributes √p factor in magnitude, equivalent to p^{+1/2}, not −1/2.

Hmm I'm not matching. Let me just abandon trying to recall the exact formula and reason from first principles.

## First-principles magnitude calculation

I'll directly compute |G_p(a)| at r=4 using the saddle structure and confirm it gives √q = p^{5/2}.

G_p(a) = Σ_{s ∈ Z/p^4} e_{p^5}(P_a(s)).

Partition s by s mod p:
- For s_0 ≠ t_1 (p − 1 classes): sub-sum vanishes (linear-in-u orthogonality as argued).
- For s_0 = t_1 (one class): sub-sum is Σ_{u ∈ Z/p^3} e_{p^5}(P_a(t_1 + p·u)).

Within this saddle class, expand P_a(t_1 + p·u) in powers of u:
- Constant: P_a(t_1)
- Linear in u: dP_a/ds|_{t_1} · p · u. Recall dP_a/ds = p · K'/(1+ps). At s = t_1, K'(t_1) = (1 + p·t_1) − C_a = 1 + p·t_1 − (1 + p·t_1 + p²·t_2 + p³·t_3) = −p²·t_2 − p³·t_3. So K'(t_1) = −p²·(t_2 + p·t_3) + O(p^4). And 1 + p·t_1 = C_a − p²·t_2 + O(p³), so (1 + p·t_1)^{-1} = C_a^{-1}·(1 − (−p²·t_2/C_a) + ...) ≈ C_a^{-1} mod p² (high precision).

dP_a/ds|_{t_1} = p · K'(t_1) / (1 + p·t_1) = p · (−p²·t_2 − p³·t_3) · (C_a^{-1} + O(p²))
              = −p³·t_2·C_a^{-1} − p^4·t_3·C_a^{-1} + O(p^5)

p · dP_a/ds|_{t_1} · u (the linear-in-u in P_a expansion) = (−p^4·t_2·C_a^{-1} − p^5·t_3·C_a^{-1} + O(p^6)) · u.

Mod q = p^5: −p^4·t_2·C_a^{-1}·u (the p^5 and beyond drop).

- Quadratic in u: (1/2)·d²P_a/ds²|_{t_1}·(p·u)². d²P_a/ds²|_{t_1} = p²·C_a/(1+p·t_1)² ≈ p²·C_a/C_a² = p²/C_a mod p². So quadratic-in-u = (1/2)·(p²/C_a)·(p²·u²) = p^4·u²/(2 C_a). Mod q: p^4·u²/(2 C_a).

- Cubic and higher: d³P_a/ds³ has v_p ≥ 3 (from log derivative), giving p^3·u³·(p)³? Actually d³P_a/ds³ = d/ds [p²·C_a/(1+ps)²] = p² · C_a · (−2p)/(1+ps)³ = −2p³·C_a/(1+ps)³, v_p = 3. So cubic-in-u in P_a is (1/6)·(−2p³C_a/C_a³)·(p·u)³ = (−p³/(3 C_a²))·p³·u³ = −p^6·u³/(3 C_a²). Mod q = p^5: ≡ 0 (since p^6 > p^5).

So at r=4, cubic and higher contributions vanish mod q. We have:

P_a(t_1 + p·u) ≡ P_a(t_1) − p^4·t_2·C_a^{-1}·u + p^4·u²/(2 C_a) mod p^5
             = P_a(t_1) + p^4·(−t_2·u + u²/2)/C_a   mod p^5

Sub-sum:
Σ_u e_{p^5}(P_a(t_1) + p^4 · (−t_2·u + u²/2)/C_a)
= e_{p^5}(P_a(t_1)) · Σ_{u ∈ Z/p^3} e_p((u²/2 − t_2·u)/C_a)
= e_{p^5}(P_a(t_1)) · p² · Σ_{u ∈ Z/p} e_p((u²/2 − t_2·u)/C_a)

(since e_p only depends on u mod p, each u mod p value appears p² times in Z/p^3.)

The remaining sum is a length-p Gauss sum:
Σ_{u ∈ Z/p} e_p((u²/2 − t_2·u)/C_a)

Complete the square: u²/2 − t_2·u = (1/2)·(u² − 2 t_2·u) = (1/2)·((u − t_2)² − t_2²).

So:
Σ_u e_p((u − t_2)²/(2 C_a) − t_2²/(2 C_a)) = e_p(−t_2²/(2 C_a)) · Σ_u e_p((u − t_2)²/(2 C_a))
                                          = e_p(−t_2²/(2 C_a)) · Σ_{w ∈ Z/p} e_p(w²/(2 C_a))      [w = u − t_2]
                                          = e_p(−t_2²/(2 C_a)) · G_{quad}(p; 1/(2 C_a))

where G_{quad}(p; α) := Σ_{w ∈ Z/p} e_p(α w²). For p ≥ 3 and α a unit mod p, |G_{quad}| = √p, and the precise value is √p · (Legendre symbol · root of unity).

**Combining:**
G_p(a) at r=4 = p² · e_{p^5}(P_a(t_1)) · e_p(−t_2²/(2 C_a)) · G_{quad}(p; 1/(2 C_a))

|G_p(a)| = p² · √p = p^{5/2} = √q. ✓ **Match.**

**The √p factor is exactly the quadratic Gauss sum G_{quad}.**

## Phase factor structure

> **G_p(a) at r=4 = p² · e_{p^5}(P_a(t_1) − p^4·t_2²/(2 C_a)) · G_{quad}**

The "phase" part e_{p^5}(P_a(t_1) − p^4·t_2²/(2 C_a)) absorbs the Gaussian-completion-of-square shift.

**Hmm, but this is the phase with the t_2² term absorbed. Is this expressible in terms of s*(r)?**

Recall s*(r) at r=4 should be (C_a − 1)/p mod p^{r-1} = mod p^3. C_a = 1 + p·t_1 + p²·t_2 + p³·t_3 mod p^4. So (C_a − 1)/p = t_1 + p·t_2 + p²·t_3 mod p^3.

So s*(r=4) = t_1 + p·t_2 + p²·t_3 mod p^3.

Now compute P_a(s*(r=4)) directly. P_a(s) = ps − C_a·L_p(1+ps) mod p^5.

s*(r=4) = t_1 + p·t_2 + p²·t_3, so ps* = p·t_1 + p²·t_2 + p³·t_3.

1 + ps* = 1 + p·t_1 + p²·t_2 + p³·t_3 = C_a − p^3·... + ... hmm wait, C_a = 1 + p·t_1 + p²·t_2 + p³·t_3 mod p^4, so 1 + ps* = C_a mod p^4. So 1 + ps* − C_a ≡ 0 mod p^4 — but C_a is only defined mod p^4 = p^r at r=4 anyway, so really 1 + ps* ≡ C_a mod p^r.

Hmm actually we want P_a(s) mod q = p^{r+1} = p^5. C_a is determined mod p^r = p^4 only. Different lifts of C_a to mod p^5 differ by p^4·k.

Let me lift C_a to a specific Z_p representative C_a' (the formal-Z_p version). Then 1 + ps* = C_a' mod p^? (it depends on the lift of s* too).

Actually the cleanest is: define s*(r=4) ∈ Z/p^3 by (C_a − 1)/p mod p^3 where C_a is treated as the unique element of Z/p^4. Then ps* ∈ Z/p^4 has ps* = C_a − 1 mod p^4 (i.e., exactly, by the definition of /p for an element divisible by p). So 1 + ps* ≡ C_a mod p^4, i.e., (1 + ps* − C_a) ≡ 0 mod p^4.

But we want P_a(s) mod p^5 = mod q. So lifts of s* and C_a from mod p^4 to mod p^5 matter.

Let s*' be a lift of s* to Z (or Z/p^4), and similarly C_a' a lift of C_a to Z/p^5. Then (1 + p·s*' − C_a') mod p^5 = ?

By construction, 1 + p·s*' − C_a' = (a multiple of p^4). Different lifts give different multiples.

This lift ambiguity is a real complication but doesn't affect the final form modulo physical interpretation. Let me compute P_a(s*) directly using the formal-saddle:

P_a(s*) ≡ p·s* − C_a · L_p(1 + p·s*) mod p^5.

L_p(1 + p·s*) = L_p(C_a + p^4·ε) for some ε ∈ Z/p (from lift ambiguity).

Using log(C_a + p^4·ε) = log(C_a) + p^4·ε/C_a + O(p^8) mod p^5 = log(C_a) + p^4·ε/C_a mod p^5.

So L_p(1+ps*) = log(C_a) + p^4·ε/C_a mod p^5.

Now log(C_a) where C_a ∈ 1 + p·Z_p: log(C_a) = Σ_{j≥1} (−1)^{j-1}/j · (C_a − 1)^j (formal series). C_a − 1 = ps* (by saddle equation). So log(C_a) = log(1 + ps*) = Σ_j (−1)^{j-1}/j · (ps*)^j.

Then C_a · L_p(1+ps*) = C_a · log(1+ps*) + p^4·ε mod p^5.

P_a(s*) = ps* − C_a · log(1+ps*) − p^4·ε mod p^5.

This still has the ε ambiguity. But the ε term doesn't affect e_q(P_a(s*)) because e_{p^5}(p^4·ε) = e_p(ε), which can take p values — these are real differences in the phase!

Wait, that can't be right — G_p(a) is a well-defined value, so P_a(s*) mod q can't depend on arbitrary lift choices.

OH I see my error. The original P_a polynomial is P_a(s) = ps − C_a · L_p(1+ps) where L_p is the truncated Cochrane log. The truncation level J_p makes L_p a polynomial in (ps), not a formal series. So my use of the "infinite" log(C_a) is wrong — I should use the truncated L_p.

Let me redo. L_p(1+ps*) = Σ_{j=1}^{J_p} (−1)^{j-1}/j · (ps*)^j.

For r=4, J_p (formula j − v_p(j) < r+1 = 5):
- p=3: j=3: 3-1=2<5✓. j=4: 4-0=4<5✓. j=5: 5-0=5, NOT<5. So J_3=4.
- p≥5: J_p = r = 4 (max j with j<5).

Both give J_p = 4 at r=4.

L_p(1+ps*) = ps* − (ps*)²/2 + (ps*)³/3 − (ps*)^4/4 mod p^5.

Hmm but for p=3, (ps*)^3/3 has the 1/3 issue. Let's see: (3·s*)^3 / 3 = 27·s*³/3 = 9·s*³ — clean integer. OK so the 1/3 is absorbed.

C_a · L_p(1+ps*) needs computation. With C_a = 1 + p·t_1 + p²·t_2 + p³·t_3 mod p^4 (lift to mod p^5 by adding p^4·t_4 for an arbitrary t_4 ∈ Z/p — different lifts):

(ps*)^j needs s* ∈ Z/p^3 (mod p^{r-1}=p^3). Let's lift s* to a specific representative s*'. Then ps*' ∈ Z/p^4. Wait, want ps* in Z/p^5 (i.e., mod q). p·s* with s* ∈ Z/p^3, p · s* ∈ Z/p^4. To get mod p^5, lift to mod p^5 (one more digit).

This lift ambiguity again... but the form is consistent: P_a(s*) mod p^5 is well-defined modulo a CONSISTENT choice of lifts of s* and C_a.

The clean way: work in Z_p. Pick canonical Z_p lifts (e.g., 0 ≤ digit < p for each base-p digit). Then everything is computed in Z_p and reduced mod p^5 at the end.

**This is the canonical way and gives a well-defined value of e_q(P_a(s*(r=4))) mod q.**

## Conjectured phase polynomial

Looking at the structure at r=3 (Phase 1): P_a(s*) = −p²·s*²/2 + O(p³).

At r=4 (after the saddle expansion at t_1 plus Gauss-completion):
P_a(t_1) − p^4·t_2²/(2 C_a) = ?

Hmm this is the Gauss-completed phase, with s = t_1 (just the first digit of s*).

Compare to P_a(s*(r=4)) where s*(r=4) = t_1 + p·t_2 + p²·t_3 (all three digits).

Let me see if these are the same. Compute P_a(s*(r=4)) directly:

s* = t_1 + p·t_2 + p²·t_3 in Z/p^3 representation.

p·s* = p·t_1 + p²·t_2 + p³·t_3.

(p·s*)² = (p·t_1)² + 2·(p·t_1)·(p²·t_2) + ... = p²·t_1² + 2·p^3·t_1·t_2 + p^4·(t_2² + 2·t_1·t_3) + O(p^5)

(p·s*)³ = (p·t_1)³ + 3·(p·t_1)²·(p²·t_2) + ... = p³·t_1³ + 3·p^4·t_1²·t_2 + O(p^5)

(p·s*)^4 = (p·t_1)^4 + ... = p^4·t_1^4 + O(p^5)

L_p(1+ps*) = ps* − (ps*)²/2 + (ps*)³/3 − (ps*)^4/4 mod p^5
           = [p·t_1 + p²·t_2 + p³·t_3] − (1/2)·[p²·t_1² + 2p³·t_1·t_2 + p^4·(t_2² + 2t_1·t_3)]
             + (1/3)·[p³·t_1³ + 3p^4·t_1²·t_2] − (1/4)·[p^4·t_1^4] + O(p^5)
           = p·t_1 + p²·(t_2 − t_1²/2) + p³·(t_3 − t_1·t_2 + t_1³/3) + p^4·(−t_2²/2 − t_1·t_3 + t_1²·t_2 − t_1^4/4) + O(p^5)

C_a · L_p(1+ps*) with C_a = 1 + p·t_1 + p²·t_2 + p³·t_3 mod p^4. Need C_a mod p^5 — but it's only defined mod p^4, so different lifts give phase differing by p^4·t_4 · L_p mod p^5 = p^4·t_4·(ps* + ...) mod p^5 = p^5·t_4·s* + ... ≡ 0 mod p^5. So lift of C_a doesn't matter (good).

Take C_a as 1 + p·t_1 + p²·t_2 + p³·t_3 (canonical representative). Multiply with L_p(1+ps*):

C_a · L_p(1+ps*) = (1 + p·t_1 + p²·t_2 + p³·t_3) · L_p(1+ps*)

To order p^5 (the mod we care about), need to track all products up to and including p^4 (since we're computing P_a(s*) mod p^5 and the leading p·s* is at p^1).

Let A_0 = p·t_1, A_1 = p²·(t_2 − t_1²/2), A_2 = p³·(t_3 − t_1·t_2 + t_1³/3), A_3 = p^4·(−t_2²/2 − t_1·t_3 + t_1²·t_2 − t_1^4/4). So L_p(1+ps*) = A_0 + A_1 + A_2 + A_3 + O(p^5).

C_a · L = (1 + p·t_1 + p²·t_2 + p³·t_3) · (A_0 + A_1 + A_2 + A_3)

Compute terms by p^k level:
- p^1: A_0 = p·t_1 → p·t_1
- p^2: 1·A_1 + p·t_1·A_0 = p²·(t_2 − t_1²/2) + p·t_1·(p·t_1) = p²·(t_2 − t_1²/2) + p²·t_1² = p²·(t_2 + t_1²/2)
- p^3: 1·A_2 + p·t_1·A_1 + p²·t_2·A_0 = p³·(t_3 − t_1·t_2 + t_1³/3) + p·t_1·p²·(t_2 − t_1²/2) + p²·t_2·p·t_1
       = p³·(t_3 − t_1·t_2 + t_1³/3) + p³·t_1·(t_2 − t_1²/2) + p³·t_1·t_2
       = p³·[(t_3 − t_1·t_2 + t_1³/3) + (t_1·t_2 − t_1³/2) + t_1·t_2]
       = p³·[t_3 + t_1·t_2 + t_1³·(1/3 − 1/2)]
       = p³·[t_3 + t_1·t_2 − t_1³/6]
- p^4: 1·A_3 + p·t_1·A_2 + p²·t_2·A_1 + p³·t_3·A_0
       = p^4·(−t_2²/2 − t_1·t_3 + t_1²·t_2 − t_1^4/4) 
         + p·t_1·p³·(t_3 − t_1·t_2 + t_1³/3)
         + p²·t_2·p²·(t_2 − t_1²/2)
         + p³·t_3·p·t_1
       = p^4·(−t_2²/2 − t_1·t_3 + t_1²·t_2 − t_1^4/4 + t_1·t_3 − t_1²·t_2 + t_1^4/3 + t_2² − t_1²·t_2/2 + t_1·t_3)
       
Collecting:
- t_1·t_3 terms: −t_1·t_3 + t_1·t_3 + t_1·t_3 = t_1·t_3
- t_1²·t_2 terms: t_1²·t_2 − t_1²·t_2 − t_1²·t_2/2 = −t_1²·t_2/2
- t_1^4 terms: −t_1^4/4 + t_1^4/3 = t_1^4·(−1/4 + 1/3) = t_1^4·(1/12)
- t_2² terms: −t_2²/2 + t_2² = t_2²/2

So at p^4: p^4·(t_2²/2 + t_1·t_3 − t_1²·t_2/2 + t_1^4/12).

P_a(s*) = p·s* − C_a · L_p(1+ps*) mod p^5

p·s* = p·t_1 + p²·t_2 + p³·t_3 + 0·p^4 (s* has only digits t_1, t_2, t_3).

P_a(s*) at p^k levels:
- p^1: p·t_1 − p·t_1 = 0 ✓
- p^2: p²·t_2 − p²·(t_2 + t_1²/2) = −p²·t_1²/2
- p^3: p³·t_3 − p³·(t_3 + t_1·t_2 − t_1³/6) = p³·(−t_1·t_2 + t_1³/6)
- p^4: 0 − p^4·(t_2²/2 + t_1·t_3 − t_1²·t_2/2 + t_1^4/12) = −p^4·(t_2²/2 + t_1·t_3 − t_1²·t_2/2 + t_1^4/12)

**So at r=4:**
> **P_a(s*) = −p²·t_1²/2 + p³·(t_1³/6 − t_1·t_2) + p^4·(−t_2²/2 − t_1·t_3 + t_1²·t_2/2 − t_1^4/12) mod p^5**

Compare to the saddle-Gauss-completion form at t_1: e_{p^5}(P_a(t_1) − p^4·t_2²/(2 C_a)). Compute P_a(t_1):

s = t_1, so ps = p·t_1.
L_p(1+p·t_1) = p·t_1 − (p·t_1)²/2 + (p·t_1)³/3 − (p·t_1)^4/4 mod p^5
            = p·t_1 − p²·t_1²/2 + p³·t_1³/3 − p^4·t_1^4/4.

C_a · L_p(1+p·t_1) = (1 + p·t_1 + p²·t_2 + p³·t_3) · (p·t_1 − p²·t_1²/2 + p³·t_1³/3 − p^4·t_1^4/4) mod p^5

Compute by level:
- p^1: p·t_1
- p^2: −p²·t_1²/2 + p·t_1·p·t_1 = −p²·t_1²/2 + p²·t_1² = p²·t_1²/2
- p^3: p³·t_1³/3 + p·t_1·(−p²·t_1²/2) + p²·t_2·p·t_1 = p³·t_1³/3 − p³·t_1³/2 + p³·t_1·t_2 = p³·(t_1³·(1/3 − 1/2) + t_1·t_2) = p³·(t_1·t_2 − t_1³/6)
- p^4: −p^4·t_1^4/4 + p·t_1·p³·t_1³/3 + p²·t_2·(−p²·t_1²/2) + p³·t_3·p·t_1
       = −p^4·t_1^4/4 + p^4·t_1^4/3 − p^4·t_1²·t_2/2 + p^4·t_1·t_3
       = p^4·(t_1·t_3 − t_1²·t_2/2 + t_1^4·(−1/4 + 1/3))
       = p^4·(t_1·t_3 − t_1²·t_2/2 + t_1^4/12)

P_a(t_1) = p·t_1 − C_a·L_p(1+p·t_1):
- p^1: 0
- p^2: −p²·t_1²/2
- p^3: −p³·(t_1·t_2 − t_1³/6) = p³·(t_1³/6 − t_1·t_2)
- p^4: −p^4·(t_1·t_3 − t_1²·t_2/2 + t_1^4/12)

So:
> P_a(t_1) = −p²·t_1²/2 + p³·(t_1³/6 − t_1·t_2) + p^4·(−t_1·t_3 + t_1²·t_2/2 − t_1^4/12) mod p^5

**Difference P_a(s*) − P_a(t_1):**
- p^2: same (both −p²·t_1²/2)
- p^3: same (both p³·(t_1³/6 − t_1·t_2))
- p^4: P_a(s*) has −p^4·(t_2²/2 + t_1·t_3 − t_1²·t_2/2 + t_1^4/12)
       P_a(t_1) has −p^4·(t_1·t_3 − t_1²·t_2/2 + t_1^4/12)
       Difference: P_a(s*) − P_a(t_1) = −p^4·t_2²/2 (since the t_1·t_3, t_1²·t_2, t_1^4 terms cancel).

So **P_a(s*) = P_a(t_1) − p^4·t_2²/2 mod p^5.**

And our Gauss-completion gave G_p(a) = p² · e_{p^5}(P_a(t_1)) · e_p(−t_2²/(2 C_a)) · G_{quad}.

Note e_{p^5}(−p^4·t_2²/2) = e_p(−t_2²/2).

Hmm wait: e_{p^5}(P_a(t_1) − p^4·t_2²/2) = e_{p^5}(P_a(t_1)) · e_{p^5}(−p^4·t_2²/2) = e_{p^5}(P_a(t_1)) · e_p(−t_2²/2).

But the Gauss-completion gave e_p(−t_2²/(2 C_a)), not e_p(−t_2²/2). These differ if C_a ≠ 1 mod p. But C_a ≡ 1 mod p, so 1/C_a ≡ 1 mod p, so 1/(2 C_a) ≡ 1/2 mod p, so e_p(−t_2²/(2 C_a)) = e_p(−t_2²/2) mod p. ✓ Match.

**So P_a(s*(r=4)) = P_a(t_1) − p^4·t_2²/2 mod p^5, and the e_{p^5} value at the "true" saddle s*(r=4) ALREADY ABSORBS the Gauss-completion phase shift.**

This means we can write:
> G_p(a) = p² · e_{p^5}(P_a(s*(r=4))) · G_{quad}    at r=4

with the phase polynomial evaluated at the TRUE saddle s*(r=4) = (C_a − 1)/p mod p^3 (clean digit extraction, no Hensel series).

## Phase polynomial structure — coefficient pattern

Looking at P_a(s*(r=4)) = −p²·t_1²/2 + p³·(t_1³/6 − t_1·t_2) + p^4·(−t_2²/2 − t_1·t_3 + t_1²·t_2/2 − t_1^4/12) mod p^5

The leading-order coefficient at each level: 
- p²: −t_1²/2 = −s*²/2 + O(p) (since s* = t_1 + O(p))
- p³: ... involves t_1³/6, t_1·t_2 which are higher-order in s*

Let me try to express this directly in terms of `ps* = p·t_1 + p²·t_2 + p³·t_3`:

(ps*)² / 2 = (p·t_1 + p²·t_2 + p³·t_3)² / 2 = p²·t_1²/2 + p³·t_1·t_2 + p^4·(t_2²/2 + t_1·t_3) + O(p^5)

Hmm let me compute P_a(s*) cleanly using a series formula.

**Key observation:** s*(r) = (C_a − 1)/p, so ps* = C_a − 1 (in Z_p exactly, for s* the Z_p representative).

Then 1 + ps* = C_a, so L_p(1+ps*) = log(C_a) in the formal-p-adic sense.

And C_a · L_p(1+ps*) = C_a · log(C_a).

So P_a(s*) = ps* − C_a · log(C_a) = (C_a − 1) − C_a · log(C_a).

**Definition:** for y ∈ pZ_p (so 1 + y ∈ principal units), define
> M(y) := y − (1+y) · log(1+y)

Then P_a(s*) = M(ps*) = M(C_a − 1).

**This is the clean closed form for the phase!**

Compute M(y):
M(y) = y − (1+y)·log(1+y)
     = y − (1+y) · (y − y²/2 + y³/3 − y^4/4 + ...)
     = y − [y − y²/2 + y³/3 − y^4/4 + ... + y² − y³/2 + y^4/3 − ... ]
     = y − y − (−y²/2 + y²) − (y³/3 − y³/2) − (−y^4/4 + y^4/3) − ...
     = −y²/2 − y³·(1/3 − 1/2) − y^4·(−1/4 + 1/3) − ...
     = −y²/2 − y³·(−1/6) − y^4·(1/12) − ...
     = −y²/2 + y³/6 − y^4/12 + y^5/20 − ...

Wait sign pattern. Let me recompute the j-th coefficient.

M(y) = y − (1+y)·log(1+y). Let log(1+y) = Σ_{j≥1} (−1)^{j-1}·y^j/j.
(1+y)·log(1+y) = log(1+y) + y·log(1+y) = Σ_j (−1)^{j-1}·y^j/j + Σ_j (−1)^{j-1}·y^{j+1}/j

= y^1/1 − y²/2 + y³/3 − y^4/4 + ... + y²/1 − y³/2 + y^4/3 − ...

= y + y²·(−1/2 + 1) + y³·(1/3 − 1/2) + y^4·(−1/4 + 1/3) + y^5·(1/5 − 1/4) + ...

= y + y²/2 + y³·(−1/6) + y^4·(1/12) + y^5·(−1/20) + ...

Coefficient of y^j for j ≥ 2: (−1)^{j-1}/j + (−1)^{j-2}/(j-1) = (−1)^{j-1}·[1/j − 1/(j-1)] = (−1)^{j-1}·[-1/(j(j-1))] = (−1)^j/(j(j-1)).

Hmm let me re-check for j=2: should be 1/2 from my calculation. Formula: (−1)^2/(2·1) = 1/2. ✓
For j=3: should be −1/6 from my calculation. Formula: (−1)^3/(3·2) = −1/6. ✓
For j=4: should be 1/12. Formula: (−1)^4/(4·3) = 1/12. ✓
For j=5: should be −1/20. Formula: (−1)^5/(5·4) = −1/20. ✓

So (1+y)·log(1+y) = y + Σ_{j≥2} [(−1)^j/(j(j-1))] · y^j.

Then M(y) = y − (1+y)·log(1+y) = − Σ_{j≥2} (−1)^j/(j(j-1)) · y^j = Σ_{j≥2} (−1)^{j-1}/(j(j-1)) · y^j

**So:**
> **M(y) = Σ_{j≥2} (−1)^{j-1}/(j(j-1)) · y^j**

This is the "(1+y)·log(1+y)" derivative-style series, with very clean coefficients!

## Closed form for the phase

**At general r ≥ 3 (saddle-exact regime):**
> **P_a(s*(r)) = M(p·s*(r)) = Σ_{j≥2} (−1)^{j-1}/(j(j-1)) · (p·s*(r))^j mod p^{r+1}**

with s*(r) = (C_a − 1)/p mod p^{r-1}.

The series M(p·s*) is truncated implicitly: terms with v_p ≥ r+1 vanish mod q. (p·s*)^j has v_p ≥ j (for s* a unit) up to v_p(j!) corrections from the 1/j(j-1) denominators. Effectively the series truncates at j = r+1 or so for p ≥ 5 (and slightly later for p=3 due to 1/3 etc.).

**This is the answer to Phase 3 (phase polynomial).** Let me hold the result for that phase.

## Phase 2 conclusion

**Saddle at r ≥ 4 (independently derived):**
> **s*(r) = (C_a − 1)/p mod p^{r-1}** (digit extraction, no Hensel iteration)

This is a **clean closed-form saddle**, not a perturbative series. The reason: the saddle equation `1 + ps − C_a = 0` (from dP_a/ds = p(1+ps − C_a)/(1+ps) = 0) is LINEAR in s, so its unique formal-p-adic root is `s = (C_a − 1)/p`, exact in Z_p. Truncation to Z/p^{r-1} gives the discrete saddle.

**Approach (b) "Structural" wins over Approach (a) "Perturbative":** the perturbative ansatz s* = s*(3) + p²·δ_1 + p³·δ_2 + ... has δ_k = (C_a's k+1-th digit), which is exactly what (C_a − 1)/p mod p^{r-1} encodes via base-p expansion. No δ_k is "computed" by a Hensel-Newton iteration; they're all read off C_a's digit sequence directly.

**Why is the closed form so clean?** The Cochrane polynomial P_a(s) = ps − C_a·L_p(1+ps) has saddle equation `dP/ds = p·(1+ps−C_a)/(1+ps) = 0`, which factors with a linear factor `(1+ps − C_a)`. This is because L_p is structurally a logarithm — and the saddle of `linear − C·log` is at the natural point where `1 + (linear/(...)) = C`. So the saddle is "intrinsically clean" due to the log structure.

## Cross-check vs PATH2_FAMILY_EXTENSION

PATH2 doc says (Phase 1 line 95-96):
> "Hensel correction at r ≥ 4 (R79b documents this at p=3). Family-level analog should hold structurally (Hensel lifting is p-blind). The specific class-correlated deviation (j=0 anomalous, j≥1 regular) needs empirical verification at general p..."

PATH2 doc treats Hensel correction as "open at family level" — but doesn't claim a specific structure. The Inner-Plancherel argument at r=3 in PATH2_BILINEAR didn't extend; it BUMPS into the Hensel problem at r ≥ 4.

**My independent derivation gives a CLEAN structure: s*(r) = (C_a − 1)/p mod p^{r-1}, phase P_a(s*) = M(ps*) = Σ (−1)^{j-1}/(j(j-1))·(ps*)^j.**

This is **structural**, not perturbative — the "Hensel correction" is just a DIGIT EXPANSION of C_a, not a series of corrections.

Whether this matches the original Hensel-lift agent's claim is the Phase 5 comparison.

## Files
- HENSEL_RECONSTRUCTION_PHASE2.md (this)
