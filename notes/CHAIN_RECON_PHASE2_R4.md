# CHAIN_RECON Phase 2 — Independent re-derivation at r=4

**Date:** 2026-05-11. Analyst: Wilson (chain-reconstruction track).

This document extends Phase 1's r=3 c_2-collapse to r=4 INDEPENDENTLY from any CHAIN_*.md material. Foundational input: Hensel closed form (HENSEL_APPROACH_A.md, doubly confirmed).

---

## Setup at r=4

- q = p^5, period = p^4, N = p^{r-1} = p^3, |supp| = p^3.
- C_a = 1 + p·s_0 + p²·c_2 + p³·c_3 mod p^4, with (s_0, c_2, c_3) ∈ (Z/p)³ parametrizing supp bijectively.
- Hensel-lifted saddle: s*(4) = (C_a−1)/p mod p^3 = s_0 + p·c_2 + p²·c_3.
- η_p(4) factor present (a-independent quadratic Gauss coefficient from the digit-chain inner Gauss sum; absorbed into the η_p(4) prefactor, modulus 1).

Hensel polynomial:
> P_a(s*(4)) = Σ_{j=2}^{4} (−1)^{j−1}·(p·s*)^j/(j(j−1)) mod p^5
>            = −p²·s*²/2 + p³·s*³/6 − p^4·s*^4/12 mod p^5.

---

## Stratum decomposition of P_a(s*(4))

Expand each term in p-adic digits of s* = s_0 + p·c_2 + p²·c_3.

**−p²·s*²/2 mod p^5:**
- s*² = s_0² + 2p·s_0·c_2 + p²·(2·s_0·c_3 + c_2²) + 2p³·c_2·c_3 + p^4·c_3².
- p²·s*² mod p^5 = p²·s_0² + 2p³·s_0·c_2 + p^4·(2·s_0·c_3 + c_2²) + O(p^5).
- −p²·s*²/2 mod p^5 = −p²·s_0²/2 − p³·s_0·c_2 − p^4·(s_0·c_3 + c_2²/2).

**+p³·s*³/6 mod p^5:**
- s*³ mod p² = s_0³ + 3p·s_0²·c_2.
- p³·s*³ mod p^5 = p³·s_0³ + 3p^4·s_0²·c_2.
- p³·s*³/6 mod p^5 = p³·s_0³/6 + p^4·s_0²·c_2/2.

**−p^4·s*^4/12 mod p^5:**
- s*^4 mod p = s_0^4.
- p^4·s*^4 mod p^5 = p^4·s_0^4.
- −p^4·s*^4/12 mod p^5 = −p^4·s_0^4/12.

**Total P_a(s*(4)) mod p^5:**
> P_a(s*(4)) ≡ p²·(−s_0²/2) + p³·(−s_0·c_2 + s_0³/6) + p^4·(−s_0·c_3 − c_2²/2 + s_0²·c_2/2 − s_0^4/12) mod p^5.

**Stratum-wise role of each digit:**

| Digit | p² stratum | p³ stratum | p^4 stratum |
|---|---|---|---|
| s_0 (class) | −s_0²/2 | s_0³/6 | −s_0^4/12 |
| c_2 (middle) | 0 | −s_0·c_2 (LIN) | s_0²·c_2/2 (LIN) + **−c_2²/2 (QUAD)** |
| c_3 (top inner) | 0 | 0 | −s_0·c_3 (LIN) |

**Critical observation:** at the p^4 stratum, the digit c_2 has BOTH a linear contribution (with coefficient s_0²/2) AND a quadratic contribution (with coefficient −1/2). The c_3 digit appears ONLY LINEARLY at p^4.

**This c_2² term — absent at r=3 because there it lived at the modulus stratum p^4 = q which drops — is the source of the chain agent's √p loss.**

---

## Two-step inner Plancherel chain

T_p := Σ_{a ∈ supp} 1̂(p·a) · G_p(a)/√q = η_p(4) · Σ_{(s_0,c_2,c_3)} 1̂(p·a(s_0,c_2,c_3)) · e_q(P_a(s*(4)))

Pull out the s_0-only phase (the "class constant"):
> T_p = η_p(4) · Σ_{s_0} const(s_0) · Σ_{c_2, c_3} 1̂(p·a) · e_{p²}(−s_0·c_2 + s_0³/6) · e_p(−s_0·c_3 − c_2²/2 + s_0²·c_2/2 − s_0^4/12)

with const(s_0) := e_{p^3}(−s_0²/2)·e_{p^2}(s_0³/6)·e_p(−s_0^4/12) (modulus 1).

### Step 1: c_3 Plancherel-collapse (top inner digit)

The c_3 phase is purely LINEAR: e_p(−s_0·c_3). And a(s_0,c_2,c_3) = a_0(s_0,c_2) + c_3·p^3·L̃_p mod p^4. So:
> 1̂(p·a) = Σ_u e_q(p·a_0·u) · e_p(c_3·L̃_p·u).

Σ_{c_3 ∈ Z/p}: Σ_{c_3} 1̂(p·a) · e_p(−s_0·c_3) = p · Σ_{u : u ≡ s'_0 mod p, 0≤u<N} e_q(p·a_0·u), where s'_0 = L̃_p^{-1}·s_0 mod p.

Substituting u = s'_0 + p·v with v ∈ Z/p^{r-2}=Z/p² (one new digit-of-u peeled off):
> Σ_{c_3} ... = p · e_q(p·a_0·s'_0) · Σ_{v=0}^{p²−1} e_q(p²·a_0·v)
>            = p · e_q(p·a_0·s'_0) · D_{p²}(a_0; p^4)

where D_{p²}(a; p^4) := Σ_{v=0}^{p²−1} e_{p^4}(a·v) (length p² mod p^4).

**Result of Step 1:**
> Inner_{c_3}(s_0, c_2) := p · e_q(p·a_0(s_0,c_2)·s'_0) · D_{p²}(a_0(s_0,c_2); p^4).

Factor p ✓ — same as r=3's c_2-collapse. Clean linear collapse.

### Step 2: c_2 sum with quadratic phase

The c_2 sum:
> H(s_0) := Σ_{c_2 ∈ Z/p} Inner_{c_3}(s_0, c_2) · e_{p²}(−s_0·c_2) · e_p(−c_2²/2 + s_0²·c_2/2)

Expand D_{p²}(a_0(s_0,c_2); p^4) using a_0(s_0,c_2) = a_0(s_0,0) + c_2·p²·L̃_p:
> D_{p²}(a_0+c_2·p²·L̃_p; p^4) = Σ_{v=0}^{p²−1} e_{p^4}(a_0(s_0,0)·v) · e_{p²}(c_2·L̃_p·v).

Also: e_q(p·a_0(s_0,c_2)·s'_0) = e_q(p·a_0(s_0,0)·s'_0) · e_q(p·c_2·p²·L̃_p·s'_0)
                                = [const(s_0)] · e_{p²}(c_2·L̃_p·s'_0)
                                = [const(s_0)] · e_{p²}(c_2·s_0)    [using L̃_p·s'_0 = s_0]

The e_{p²}(c_2·s_0) factor exactly cancels with the e_{p²}(−s_0·c_2) from the original phase:
> e_{p²}(c_2·s_0 − s_0·c_2) = 1.

So the surviving c_2 phase (after the c_3-collapse cancellation):
> H(s_0) = p · [const(s_0)] · Σ_v e_{p^4}(a_0(s_0,0)·v) · G(v, s_0)

where
> G(v, s_0) := Σ_{c_2 ∈ Z/p} e_{p²}(c_2·L̃_p·v) · e_p(c_2·M(v, s_0) − c_2²/2),
> M(v, s_0) := s_0²/2 mod p (with v_1-dependent corrections at sub-leading order).

(Writing v = v_0 + p·v_1, the e_{p²}(c_2·L̃_p·v) factor decomposes as e_{p²}(c_2·v_0)·e_p(c_2·(v_1 + λ·v_0)) where L̃_p = 1+p·λ; so M effectively absorbs v_1 + λ·v_0.)

### Magnitude of G(v, s_0): the √p loss

**Case A: v_0 = 0.** Then e_{p²}(c_2·v_0) = 1; G reduces to a complete quadratic Gauss sum at modulus p:
> G(v_0=0, v_1, s_0) = Σ_{c_2=0}^{p−1} e_p(c_2·M_A − c_2²/2)

with M_A = v_1 + s_0²/2 mod p. **|G| = √p exactly** (standard quadratic Gauss).

**Case B: v_0 ≠ 0.** e_{p²}(c_2·v_0) is a genuine p²-th root; G is an incomplete quadratic Gauss sum at modulus p² over the length-p interval c_2 ∈ {0,...,p−1}. By standard Polya-Vinogradov for quadratic phases (completing the square, then Dirichlet kernel bound):
> |G(v_0 ≠ 0, v_1, s_0)| ≤ 2√p uniformly.

**Uniform bound: |G(v, s_0)| ≤ 2√p for all v ∈ Z/p², s_0 ∈ Z/p.**

This is the **single √p loss at the c_2-quadratic stratum** flagged by the chain agent. My independent derivation reproduces it AT r=4.

### Bounding |H(s_0)| at r=4

H(s_0) = p · [const(s_0)] · Σ_v e_{p^4}(a_0(s_0,0)·v) · G(v, s_0).

|H(s_0)| ≤ p · Σ_v |e_{p^4}|·|G(v, s_0)| ≤ p · p² · 2√p = 2 · p^{7/2}.

Tighter (Cauchy on v): |H(s_0)| ≤ p · √(p²) · √(Σ_v |G|²) ≤ p · p · √(p²·4p) = 2 · p^{7/2}. Same.

Using the Plancherel-on-v_1 structure (G depends on v through (v_0, v_1)):
- For v_0=0 (p values of v): G = √p · η · e_p(M_A²/2). The v_1 dependence is in M_A = v_1 + s_0²/2. Σ_{v_1} e_{p^3}(a_0·v_1) · √p · e_p((v_1+s_0²/2)²/2) is a length-p mixed-modulus Gauss sum, magnitude ≤ √p. So Σ_{v: v_0=0} |...| ≤ √p · √p = p.
- For v_0 ≠ 0: trivially Σ_{v_1} ≤ p · 2√p = 2 p^{3/2}. Σ over (p−1) nonzero v_0: ≤ (p−1)·2·p^{3/2} ≈ 2·p^{5/2}.

Sum: Σ_v |e_{p^4}(a_0·v) · G| ≈ p + 2·p^{5/2} ≈ 2·p^{5/2}.

|H(s_0)| ≤ p · 2·p^{5/2} = 2·p^{7/2}.

### Outer s_0 sum

|T_p| ≤ Σ_{s_0} |H(s_0)| ≤ p · 2·p^{7/2} = 2·p^{9/2}.

**At r=4, N=p^3, so p^{9/2} = N·p^{3/2}.** Hence

> **|T_p| ≤ 2·p^{3/2}·N at r=4** from direct triangulated bound.

This is **a factor √p WORSE than the chain agent's claim of 2√p·N**. The discrepancy comes from the v-sum triangle being LOOSE — the chain agent presumably extracts an additional √p of cancellation from the v-direction via Plancherel against the e_{p^4}(a_0·v) phase.

### Attempting to close the √p gap (Plancherel-on-v argument)

Try Σ_v e_{p^4}(a_0·v) · G(v, s_0) more carefully. The function G(v, s_0) has magnitude ≤ 2√p uniformly. Apply Plancherel on Z/p²:

> |Σ_v e_{p^4}(a_0·v) · G(v, s_0)|² ≤ ||e_{p^4}(a_0··)||²·||G||² (by Cauchy)
>                                     ≤ p² · Σ_v |G|²
>                                     ≤ p² · p² · (2√p)² = 4·p^5.

So |Σ_v| ≤ 2·p^{5/2}, **matching the Case B trivial bound**. Cauchy doesn't help here because both factors saturate at their bounds simultaneously.

**Alternative:** treat the v sum as a "diagonal" Plancherel of G against the Dirichlet kernel D_{p²}. Specifically:

Σ_v e_{p^4}(a_0·v) · G(v, s_0) ≈ p² · max-value of G, with the e_{p^4}(a_0·v) ≈ Dirichlet kernel "selecting" the v where G is largest. This doesn't help either.

**My best rigorous bound at r=4: |T_p| ≤ 2·p^{3/2}·N = 2·p·√p·N**, equivalent to |S_partial| ≤ 2·p^{3/2}·√N.

Looser than chain agent's claim by a factor of p (in my direct triangulation).

---

## Discussion: where the chain agent might get strict 2√p·N

Re-reading the chain agent's pre-reg disposition (without consulting CHAIN_*.md files): the claim invokes a "nested inner-Plancherel chain on (c_2, c_3, ..., c_{r-1})" and asserts each digit's Plancherel gives factor p (linear) or √p (quadratic Gauss).

At r=4: c_3 collapse gives p (verified above). c_2 sum gives √p (verified above). Net Plancherel save: p × √p = p^{3/2}. Trivial bound on the bilinear: |T_p| ≤ N² (each |1̂| ≤ N, |supp| = N, no cancellation). With save factor p^{3/2}: |T_p| ≤ N²/p^{3/2}. At r=4: N²/p^{3/2} = p^6/p^{3/2} = p^{9/2} = √p · p^4 = √p · p · N.

Hmm — the "trivial / save = p^{9/2}" matches the chain agent's claim of 2√p·N IF we identify N = p^4 (i.e., q rather than p^{r-1}). The chain agent might be parametrizing N differently.

**Alternative interpretation:** if "√N" in the chain agent's claim means "√(supp size) = √(p^{r-1}) = p^{(r-1)/2}", then 2√p·√N at r=4 = 2√p·p^{3/2} = 2p^2. And |T_p| ≤ 2·p^2 = 2·N^{2/3} would be much TIGHTER than my bound. That doesn't match either.

**Most likely interpretation:** the chain agent's "|S_partial| ≤ 2√p · √N" maps to |T_p| ≤ 2√p · N in our notation, equivalently |K|/√N ≤ 2√p (matching R79b's empirical max 2.7 ≈ 2·√3 at p=3, r large).

My bound 2·p^{3/2}·N would correspond to |K|/√N ≤ 2·p^{3/2}/p = 2√p · ... wait let me recompute.

|K| = (p/√q)·|T_p|. With |T_p| ≤ C·N: |K| ≤ C·p·N/√q = C·p·p^{r-1}/p^{(r+1)/2} = C·p^{(r-1)/2 + 1 − (r+1)/2 + 1} = C·p^{1/2 + 1 − 1/2} = C·p.

Hmm so |K| ≤ C·p with |T_p| ≤ C·N. Then |K|/√N ≤ C·p/√N = C·p/p^{(r-1)/2}. For r=4: C·p/p^{3/2} = C/p^{1/2}. At p=3: C/√3 ≈ C·0.58.

For C = 2 (strict): |K|/√N ≤ 1.15. R79b empirical max at p=3 large r: 2.7. **My bound predicts |K|/√N ≤ 1.15, but R79b shows up to 2.7 → over-prediction by my rigorous bound (rigorous bound is LARGER than empirical, which is fine).**

Wait, I confused myself. Let me re-derive. |T_p| ≤ C·N → |K| ≤ (p/√q)·C·N = C·p·N/√q. For r=4, p=3: C·3·27/√243 = 81·C/15.6 ≈ 5.2·C. And √N = √27 ≈ 5.2. So |K|/√N = 5.2·C/5.2 = C.

So |T_p| ≤ C·N corresponds to |K|/√N ≤ C. R79b max at p=3, r=8..20 is ~2.7. So C = 2 (strict) gives |K|/√N ≤ 2 — tight at p=3 (empirical max ~2.7 has sampling noise).

My bound |T_p| ≤ 2·p^{3/2}·N corresponds to |K|/√N ≤ 2·p^{3/2}. At p=3: 2·3·√3 ≈ 10.4. Way looser than R79b empirical max 2.7 — my bound is consistent with R79b but very loose.

The chain agent's 2√p·N corresponds to |K|/√N ≤ 2√p. At p=3: 2·√3 ≈ 3.46. Slightly above R79b max 2.7 — TIGHT.

**So both bounds are consistent with R79b (which only sets an EMPIRICAL lower bound on the rigorous bound), but the chain agent's claim is tighter than my direct triangulated bound by a factor of p.**

---

## Hypothesis on the gap

The chain agent's 2√p·N requires an additional √p of cancellation in the v-direction that I'm losing via triangle inequality. Specifically, the v-direction has a Plancherel-on-v_1 structure (length-p Dirichlet kernel inside a length-p² sum) that, combined with the quadratic Gauss G(v, s_0), produces a √p save.

**My honest assessment:** I see the √p loss from the c_2-quadratic stratum (matches chain agent). I don't see how to extract the second √p save from the v-direction without a more delicate argument. The chain agent may be correct OR may be assuming a cancellation that doesn't hold rigorously.

The HENSEL_APPROACH_A.md's analog argument (lines 254-358) sketches the same chain at r=4 but ALSO doesn't fully rigorously close to strict 2√N or 2√p·√N — it concludes with "2·p^{(r+3)/2} = 2·√p·N" at r=4 with caveat "sub-trivial but does NOT reach strict 2N". The doc notes that whether the deeper-digit Plancherel chain works to save the additional p depends on "whether the deeper-digit Plancherel chain introduces additional structural factors (like a log per digit) or stays at constant per digit".

**My independent derivation at r=4:** the c_2 quadratic stratum produces a √p loss. Whether further cancellation reduces the bound from `√p · p · N` to `√p · N` depends on the v-direction structure, which I cannot rigorously establish from the materials I've consulted.

**Provisional Phase 2 result:**

> **At r=4:** |T_p| ≤ 2·p^{3/2}·N rigorously (my triangulated bound). The chain agent's tighter claim 2·√p·N depends on an additional v-direction cancellation not rigorously established here. **At r=4 specifically, the chain agent's claim is PLAUSIBLE (single √p loss at c_2) but the cleaner extra-√p save in v-direction is unverified in my derivation.**

---

## Adversarial checks at r=4

(A1) **Empirical anchor.** R79b at p=3, large r: |K|/√N ≤ 2.7. My bound: |K|/√N ≤ 2·p^{3/2} = 10.4 at p=3. Chain agent: 2·√p = 3.46 at p=3. Both consistent with R79b (R79b is below both bounds). Chain agent is tighter.

(A2) **r=3 reduction.** At r=3, c_2² lives at p^4 stratum = q → drops mod q. So no √p loss at r=3. My Phase 1 derivation gives strict |T_p| ≤ 2N. ✓ Consistent.

(A3) **Top-digit cleanness.** c_3 (top inner digit at r=4) appears LINEARLY only. Plancherel-collapse gives clean factor p. ✓ Chain agent correct on this.

(A4) **Honesty check.** I did not consult any CHAIN_*.md or TIGHTEN_*.md. The HENSEL_APPROACH_A.md sketch of the r=4 chain at lines 254-358 is consistent with my independent derivation (both reach |T_p| ≤ √p·p·N triangulated; both flag the cleaner save as not rigorously established). ✓

---

## Summary at r=4

- The c_2-quadratic stratum p^4 produces a single √p loss in the c_2 Plancherel sum, **CONFIRMED** independently.
- The c_3 (top inner digit) is linear → clean factor p collapse, **CONFIRMED**.
- My direct rigorous bound: |T_p| ≤ 2·p^{3/2}·N at r=4 (one factor √p looser than chain agent).
- Chain agent's claim 2·√p·N requires additional v-direction cancellation that I cannot rigorously establish from the closed form alone.

**At r=4 the chain agent's claim is PLAUSIBLE but my independent derivation can only confirm |T_p| ≤ √p·p·N.** Phase 3 examines whether the situation worsens (more √p losses) at r ≥ 6.
