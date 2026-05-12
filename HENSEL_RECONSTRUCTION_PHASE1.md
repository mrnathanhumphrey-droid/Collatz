# HENSEL_RECONSTRUCTION_PHASE1 — Re-derivation of r=3 saddle from R78.4-78.6 (family-level)

**Date:** 2026-05-11. Analyst: independent re-derivation agent (mirror of Path 2 Pushback Check 6).

**Constraint:** This document deliberately does NOT use the original Hensel-lift derivation (HENSEL_APPROACH_A.md, HENSEL_PHASE_ARTICULATION.md, HENSEL_DISPOSITION.md, HENSEL_NUMERICAL_VERIFICATION.md). It is derived from the foundational materials only: R78.4-78.6, PATH2_FAMILY_EXTENSION, PATH2_BILINEAR_FROM_CLOSED_FORM, FHAT verification, R79b.

## Goal of Phase 1

Reproduce, independently, the saddle-point closed form at r=3 stated in R78.6 (q=3) and PATH2_FAMILY_EXTENSION.md (family p ≥ 3):
- saddle: `s*(C_a) = (C_a − 1)/p mod p`
- phase: `G_p(a) = √q · e_q(P_a(s*(C_a)))` (exact at r=3 with J_p = 3)

## Starting line

From PATH2_FAMILY_EXTENSION T78.4_p, the inner Gauss sum at family p ≥ 3, fixed unit c (taking c=1 for clarity, matching R78 convention):

> G_p(a) = Σ_{s=0}^{p^r − 1} e_q(P_a(s))
>
> P_a(s) = ps − C_a · L_p(1+ps),    L_p(1+ps) = Σ_{j=1}^{J_p} (−1)^{j−1}/j · (ps)^j

with `q = p^{r+1}`, `period = p^r`, `C_a = a · L̃_p^{−1} mod p^r`, `L̃_p = L_p(1+p)/p`. By T78.5_p, `a ↔ C_a` is a bijection on {a ≡ 1 mod p in Z/p^r}, and `C_a ≡ 1 mod p` for a in this coset.

**For r=3, p ≥ 3:** the truncation `J_p` is the maximum j with `j − v_p(j) < r+1 = 4`.
- For p = 3: J = 3 (since 3 − v_3(3) = 3 − 1 = 2 < 4; 4 − v_3(4) = 4 < 4 fails).
   Wait: 4 − 0 = 4, not less than 4. So J = 3 for p=3 at r=3.
   Re-check: actually j=4: 4 − v_3(4) = 4 < 4? No, 4 < 4 is false. So J=3 yes.
- For p ≥ 5: v_p(j) = 0 for j < p, so j − v_p(j) = j. J = max j with j < 4 ⟹ J = 3.

**So at r=3, J_p = 3 for all p ≥ 3.** Consistent across the family.

## Step 1: write out P_a(s) mod p^4 at r=3

Expand:
- L_p(1+ps) = ps − (ps)²/2 + (ps)³/3 mod p^4  (since J=3)
- P_a(s) = ps − C_a · L_p(1+ps)
- = ps − C_a · [ps − p²s²/2 + p³s³/3]
- = ps − C_a·ps + C_a·p²s²/2 − C_a·p³s³/3
- = ps(1 − C_a) + (C_a·p²/2)·s² − (C_a·p³/3)·s³

The arithmetic of `1/2` and `1/3` mod p^4 needs p ≥ 5 to be entirely clean; for p=3 the 1/3 introduces v_3 = −1 which is absorbed by the p³ factor (giving net v_3 = 2). At family level for p ≥ 5, both 2 and 3 are units mod p^4, so the polynomial sits in Z/p^4 cleanly.

## Step 2: derive saddle by stationary-phase

A saddle point of P_a (the place where the inner sum concentrates) is where dP_a/ds = 0 mod some appropriate power of p. Differentiate:

> dP_a/ds = p − C_a · dL_p(1+ps)/ds = p − C_a · (p − p²s + p³s² − ...) = p(1 − C_a) + C_a·p²·s − C_a·p³·s² + O(p^4)

**Leading order analysis.** Mod p²:
- dP_a/ds ≡ p(1 − C_a) mod p²

Since C_a ≡ 1 mod p (by T78.5_p restriction to coset), write C_a = 1 + p·t_1 + p²·t_2 + ... where t_1, t_2 ∈ Z/p. Then 1 − C_a = −p·t_1 mod p², so dP_a/ds ≡ −p²·t_1 mod p², which is ≡ 0 mod p² automatically.

Mod p³:
- dP_a/ds ≡ p(1 − C_a) + C_a·p²·s
- ≡ −p²·t_1 − p³·t_2 + (1 + p·t_1)·p²·s mod p³
- ≡ −p²·t_1 + p²·s + p³·(...) mod p³
- ≡ p²·(s − t_1) mod p³

**Setting dP_a/ds ≡ 0 mod p³ gives `s ≡ t_1 mod p`.**

Recall `t_1 = (C_a − 1)/p mod p` — this is exactly the digit-extraction `s*(C_a) := (C_a − 1)/p mod p`.

**Independently derived: at r=3, s*(C_a) = (C_a − 1)/p mod p.** ✓ Matches PATH2_FAMILY_EXTENSION T78.6_p phase prediction.

## Step 3: evaluate P_a(s*) mod p^4

Substitute s = s* + p·s_1 + ... where the leading s_1 is yet-to-determine if needed for r=3 phase closure. (We'll show the leading representative s = s* is sufficient at r=3 because higher-order corrections drop out mod p^4.)

For r=3 mod p^4:

P_a(s) = p·(1 − C_a)·s + (C_a·p²/2)·s² − (C_a·p³/3)·s³ mod p^4

Substitute C_a = 1 + p·t_1 + p²·t_2 mod p^3 (recall C_a lives in Z/p^r = Z/p^3; the lift to Z/p^4 has ambiguity p^3·Z/p^4 which contributes p^3·s to P_a — but at our final exponential, e_q(p^4·k) = 1, so this lift ambiguity affects P_a(s) only mod p^4 once s is fixed; we'll see this carefully).

Actually, careful: in the sum G_p(a) = Σ_s e_q(P_a(s)), each summand uses C_a interpreted mod p^r. The polynomial P_a(s) — for the inner sum mod q=p^{r+1} — needs C_a known mod p^{r+1} effectively. The trick: any lift of C_a from Z/p^r to Z/p^{r+1} differs by p^r·k for k ∈ Z, and P_a depends on C_a via terms like C_a·ps (giving p^{r+1}·k·s ≡ 0 mod q), C_a·p²s²/2 (giving p^{r+2}·... ≡ 0), etc. So any lift suffices. We pick the "natural" lift via the L_p formula extended to L_p(1+p) (which is a Z_p element); the resulting C_a is well-defined mod p^r and lifts don't affect e_q(P_a).

OK with that resolved, substitute s = t_1 (the leading representative of s* mod p) into P_a:

P_a(t_1) = p(1 − C_a)·t_1 + (C_a·p²/2)·t_1² − (C_a·p³/3)·t_1³ mod p^4

Compute 1 − C_a = −(p·t_1 + p²·t_2) mod p^3, so:
- p(1 − C_a)·t_1 = −p²·t_1² − p³·t_1·t_2 mod p^4
- C_a·p²/2 · t_1² = (1 + p·t_1)·p²·t_1²/2 + O(p^4) = p²·t_1²/2 + p³·t_1³/2 mod p^4
- C_a·p³/3 · t_1³ = p³·t_1³/3 + O(p^4) mod p^4

Sum:
- P_a(t_1) = −p²·t_1² − p³·t_1·t_2 + p²·t_1²/2 + p³·t_1³/2 − p³·t_1³/3 mod p^4
- = p²·t_1²·(−1 + 1/2) + p³·(−t_1·t_2 + t_1³/2 − t_1³/3)
- = −p²·t_1²/2 + p³·(t_1³/6 − t_1·t_2) mod p^4

**Result:**

> **P_a(s* = t_1) ≡ −p²·t_1²/2 + p³·(t_1³/6 − t_1·t_2) mod p^4**

This matches the derivation in PATH2_BILINEAR_FROM_CLOSED_FORM.md lines 207-217 (with the identification t_1 = s*, t_2 = c_2).

## Step 4: account for the saddle multiplicity / Gaussian-integration factor

The saddle is `s ≡ t_1 mod p`. But the full sum G_p(a) ranges over s ∈ {0, 1, ..., p^r − 1} (i.e., p^r values of s in Z/p^r). The saddle congruence `s ≡ t_1 mod p` selects p^{r-1} values of s out of p^r. For each fixed s ≡ t_1 mod p, write s = t_1 + p·u for u ∈ Z/p^{r-1}.

**Now evaluate P_a(t_1 + p·u) mod p^4 (at r=3):**

P_a(s) = p(1−C_a)·s + (C_a·p²/2)·s² − (C_a·p³/3)·s³ mod p^4

Substitute s = t_1 + p·u:
- s² = t_1² + 2p·t_1·u + p²·u²
- s³ = t_1³ + 3p·t_1²·u + 3p²·t_1·u² + p³·u³

So:
- p(1−C_a)·s = p(1−C_a)·t_1 + p²(1−C_a)·u
- (C_a·p²/2)·s² = (C_a·p²/2)·(t_1² + 2p·t_1·u + p²·u²) = (C_a·p²/2)·t_1² + C_a·p³·t_1·u + O(p^4)
- (C_a·p³/3)·s³ = (C_a·p³/3)·t_1³ + O(p^4)

So:
- P_a(t_1 + p·u) = P_a(t_1) + [p²(1−C_a) + C_a·p³·t_1] · u + O(p^4)

Substitute 1 − C_a = −p·t_1 − p²·t_2 mod p^3:
- p²(1 − C_a) = −p³·t_1 − p^4·t_2 ≡ −p³·t_1 mod p^4
- C_a·p³·t_1 = (1 + p·t_1)·p³·t_1 + O(p^5) = p³·t_1 + O(p^4) mod p^4

Sum:
- [p²(1−C_a) + C_a·p³·t_1] · u = (−p³·t_1 + p³·t_1)·u + O(p^4) = O(p^4)

**So at r=3, P_a(t_1 + p·u) ≡ P_a(t_1) mod p^4 — independent of u!**

This is the crucial fact: at the saddle, the phase P_a is constant across all s ≡ t_1 mod p in Z/p^r. There are p^{r-1} = p² such s values, all contributing the SAME phase.

For s ≢ t_1 mod p: the linear-in-shift term `p²·(s − t_1) mod p³` from dP_a/ds is non-zero, so P_a(s) cycles through values that average to zero (a complete sum of e_q over a coset of p²·Z/q ≅ Z/p, which has cardinality p — vanishes by orthogonality unless s ≡ t_1 mod p).

## Step 5: combine — closed form at r=3

G_p(a) = Σ_{s=0}^{p^r − 1} e_q(P_a(s))

Partition s mod p:
- For s ≡ s_0 mod p with s_0 ≠ t_1: the sub-sum vanishes (orthogonality on the p^{r-1} values; phase varies linearly with quotient).
- For s ≡ t_1 mod p: the phase is constant = e_q(P_a(t_1)), and there are p^{r-1} such s values.

Wait — I need to be careful. The argument that "phase varies linearly with quotient" needs verification. Let me redo for s = s_0 + p·u with s_0 ≠ t_1:

dP_a/ds at s = s_0 ≡ p²·(s_0 − t_1) mod p³ (from Step 2's expansion).

P_a(s_0 + p·u) = P_a(s_0) + p²·(s_0 − t_1)·u + O(p^3 · u) + O(p^4)

Hmm, the dependence on u is at level p² for s_0 ≠ t_1. Then:
- e_q(P_a(s_0 + p·u)) = e_q(P_a(s_0)) · e_q(p²·(s_0 − t_1)·u + O(p^3·u))

Sum over u ∈ Z/p^{r-1}:
- Σ_u e_q(p²·(s_0 − t_1)·u + ...) = Σ_u e_{p^{r-1}}((s_0 − t_1)·u) · [correction]

Wait, e_q(p²·k·u) = e_{p^{r-1}}(k·u) when q = p^{r+1}. So this is a complete sum over u ∈ Z/p^{r-1} of e_{p^{r-1}}((s_0 − t_1)·u), which is **p^{r-1}** when (s_0 − t_1) ≡ 0 mod p^{r-1} and **0** otherwise. Since s_0 − t_1 ∈ {1, ..., p-1} (units mod p), we have (s_0 − t_1) ≢ 0 mod p^{r-1} for r ≥ 2 (and the deeper Hensel-correction terms could contribute but at leading order this is exact).

Actually wait, this needs more care for r ≥ 4 where the Hensel correction kicks in. At r=3 specifically:
- For r=3, the leading sub-sum is over u ∈ Z/p² (since r-1 = 2).
- The leading u-dependence in P_a(s_0 + p·u) is `p²·(s_0 − t_1)·u`, giving e_{p^{r-1}=p²}((s_0 − t_1)·u) when divided by p² inside e_q.
- This is a complete sum of e_{p²} over u ∈ Z/p² — vanishes unless s_0 − t_1 ≡ 0 mod p².
- Since s_0, t_1 ∈ {0, ..., p-1}, the difference can be 0 only when s_0 = t_1.

**So at r=3, the sub-sum vanishes for s_0 ≠ t_1.** ✓

For s_0 = t_1: the phase is exactly constant P_a(t_1) mod p^4 by the Step 4 computation, all p² = p^{r-1} terms add coherently.

**Combine:**

> **G_p(a) = p^{r-1} · e_q(P_a(t_1))** at r=3

with `t_1 = s*(C_a) = (C_a − 1)/p mod p`.

## Step 6: check the magnitude prediction

|G_p(a)| = p^{r-1} = p²  at r=3.

But the empirical (FHAT-verified) magnitude is `|G_p(a)| = √q = p^{(r+1)/2} = p²` at r=3.

These match! p^{r-1} = p^{(r+1)/2} requires 2(r-1) = r+1, i.e., r=3. ✓

**This is the saddle-exact case: r=3 with J_p=3 gives p^{r-1} = √q exactly.**

For r=2: p^{r-1} = p but √q = p^{(r+1)/2} = p^{3/2}, discrepancy factor √p. The saddle calculation at r=2 picks up an additional Gauss-sum factor of magnitude √p (this is Phase 4's job to derive).

For r ≥ 4: the Step 4 "phase constant on coset" fails because higher-order Hensel-correction terms make P_a(t_1 + p·u) depend on u at level p^{r-1} or below. The leading-order saddle t_1 needs Hensel correction. This is Phase 2's job.

## Step 7: closed form at r=3

> **G_p(a) = p^{r-1} · e_q(P_a(s*(C_a))),    s*(C_a) = (C_a − 1)/p mod p, r=3**

Equivalently:
> **G_p(a) / √q = e_q(P_a(s*(C_a)))** at r=3 (since p^{r-1} = √q at r=3).

This is the **R78.6 statement** verbatim, generalized to family-level p ≥ 3, derived independently.

## Cross-check with PATH2_FAMILY_EXTENSION T78.6_p

Quoted from foundational doc:
> Define s*(C_a) = (C_a − 1)/p mod p.
> At r=3: G_p(a) = p^{(r+1)/2} · e_q(P_a(s*(C_a))) exact (analog of R78.6 at q=3 r=3 with J=3).

**Independent derivation lands at the same place.** ✓

## Honesty check: what was used / not used

**Used (foundational):**
- T78.4_p factorization G_p(a) = Σ e_q(P_a(s)), P_a = ps − C_a·L_p(1+ps).
- C_a = a·L̃_p^{−1} mod p^r, C_a ≡ 1 mod p (T78.5_p).
- Cochrane truncated p-adic log expansion L_p(1+ps) = Σ (−1)^{j−1}/j (ps)^j.
- Orthogonality of additive characters: Σ_{u ∈ Z/n} e_n(k·u) = n·δ_{k≡0}.

**NOT used:**
- Any specific Hensel-lift construction (the original derivation files).
- Any non-Cochrane machinery.
- Any deep arithmetic-geometry input.

**The derivation is pure stationary-phase + character orthogonality on R78.4-78.6's polynomial.**

## Phase 1 conclusion

Saddle at r=3 (family-level, p ≥ 3): **s*(C_a) = (C_a − 1)/p mod p**, **G_p(a) = √q · e_q(P_a(s*))**. Independent derivation reproduces T78.6_p saddle prediction.

The derivation chain is Cochrane factorization → stationary-phase to first order → orthogonality vanishing of non-saddle classes → coherent addition on saddle class → closed form.

**This is exactly the standard p-adic stationary-phase argument applied to the Cochrane polynomial.** Nothing exotic; nothing q=3-specific.

Next: Phase 2 — extend to r ≥ 4 where the saddle picks up Hensel corrections.
