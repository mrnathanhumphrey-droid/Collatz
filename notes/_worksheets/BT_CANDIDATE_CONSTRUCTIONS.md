# BT_CANDIDATE_CONSTRUCTIONS — Phase 1 disposition for Candidates A / B / C

**Date:** 2026-05-12. Probe BT, Phase 1. Gate: at least one of A/B/C must show feasibility for Phase 2.

---

## Setup

Bruhat-Tits tree T_3 for PGL_2(ℚ_3):
- Vertices = equivalence classes of ℤ_3-lattices in ℚ_3² (Lubotzky 2013, Thm 1.2.1).
- Each vertex has degree q + 1 = 4.
- Boundary ∂T_3 = ℙ¹(ℚ_3); the standard vertex v_0 corresponds to the standard lattice [L_0].
- An odd integer r corresponds to the boundary point [r : 1] ∈ ℙ¹(ℚ_3).

Tao recursion: r odd, r → r' = (3r + 1) / 2^v where v = v_2(3r + 1).

---

## Candidate A — single SL_2(ℚ_3) lift via Möbius transformations

**Construction:** Find M_T ∈ SL_2(ℚ_3) (or PGL_2(ℚ_3)) such that Tao recursion is r → M_T · r in projective coordinates.

**Analysis:** Tao at fixed v corresponds to the linear-fractional map z → (3z + 1) / 2^v, with matrix

    M_v = [[3, 1], [0, 2^v]]   in GL_2(ℚ_3).

Projectivized, [M_v] ∈ PGL_2(ℚ_3) depends on v: M_1 = [[3,1],[0,2]] and M_2 = [[3,1],[0,4]] have different normalized forms (top rows [3/2, 1/2] vs [3/4, 1/4]) and represent different elements of PGL_2(ℚ_3).

Since 2 ∈ ℤ_3^* has infinite multiplicative order in (ℤ_3^*)_torsion-free, no two M_v with v ≠ v' are projectively equal.

**Pre-registered hypothesis H_BT_A_DIRECT:** FALSIFIED.

The Tao map does not lift to a single SL_2(ℚ_3) matrix. It requires the family {M_v : v ≥ 1} indexed by the 2-adic valuation of each step.

**Salvageable structural finding (not the hypothesis, but worth noting):**

Each M_v has discriminant (trace² − 4 det) = (3 + 2^v)² − 12·2^v = (2^v − 3)², an exact square in ℚ. So M_v is **hyperbolic** in PGL_2(ℚ_3) with eigenvalues {3, 2^v}, both rational integers. Translation length on T_3 equals |v_3(eigenvalue ratio)| = |v_3(3 / 2^v)| = 1.

So every M_v translates by length 1 along its own translation axis. The axes are different for each v (one endpoint always [1:0] = ∞, the other at [1 : 2^v − 3]); only M_2 has [1:1] (= the Collatz attractor) as a fixed point.

This is *suggestive* of T_3 structure but does not save Candidate A's hypothesis — the single-matrix lift fails.

**Disposition: FAIL.**

---

## Candidate B — random walk on T_3 with v-randomness

**Construction:** Generators {M_v : v = 1, 2, 3, …}. Walk at step n: choose v_n ~ Geom(1/2), apply M_{v_n}. Hypothesis: walk is transient with stationary harmonic measure concentrated at the "1-cusp" [1 : 1] ∈ ∂T_3.

**Analysis:** I tested the **3-adic convergence** of the Tao trajectory of r_0 = 27 to 1:

| step | r_n | r_n − 1 | v_3(r_n − 1) |
|------|-----|---------|--------------|
| 0 | 27 | 26 | 0 |
| 1 | 41 | 40 | 0 |
| 2 | 31 | 30 | 1 |
| 7 | 121 | 120 | 1 |
| 8 | 91 | 90 | 2 |
| 33 | 577 | 576 | 2 |
| 35 | 325 | 324 | 4 |
| 36 | 61 | 60 | 1 |
| 39 | 53 | 52 | 0 |

v_3(r_n − 1) bounces around between 0 and 4 with no monotone trend. The Tao trajectory does **NOT** converge to 1 in the 3-adic topology. It converges archimedean-ly only after the integer recursion terminates.

Equivalent diagnostic: r_n mod 9 for r_0 = 27 is

    [0, 5, 4, 2, 8, 8, 8, 4, 1, 2, 4, 2, 8, 4, 2, 8, 8, 4, 5, 8, …]

This is essentially a uniform walk on the unit residues mod 9, **not** convergent to 1 mod 9.

**Mechanism (verified algebraically):** r_{n+1} mod 3 = (3 r_n + 1)/2^v mod 3 = 1 / 2^v mod 3 = 2^{−v} mod 3 ∈ {1, 2} depending on parity of v. So r_{n+1} mod 3 depends only on v_n parity. **All higher 3-adic digits of r_n+1 are *independent* of r_n's higher 3-adic digits.** (Cross-reference: project_collatz_prefix_nonpropagation memory.)

**Pre-registered hypothesis H_BT_B_RANDOM:** FALSIFIED.

The induced random walk on T_3 does **not** have stationary harmonic measure concentrated at the 1-cusp. It mixes uniformly over the 3-adic residues coprime to 3, with no preferred boundary direction.

**Disposition: FAIL.**

---

## Candidate C — arithmetic quotient of T_3

**Construction:** Identify (ℤ/3^n ℤ)^* with depth-n vertices in some arithmetic quotient Γ\T_3. Tao recursion = depth-walk going one step deeper. Hypothesis: trajectories descend toward an apex (1-cusp).

**Test 1 — strong cardinality bijection:**

| n | (ℤ/3^n)^* | depth-n sphere in T_3 (= 4·3^{n−1}) |
|---|-----------|--------------------------------------|
| 1 | 2 | 4 |
| 2 | 6 | 12 |
| 3 | 18 | 36 |
| 4 | 54 | 108 |

The depth-n sphere of T_3 is exactly 2× the size of (ℤ/3^n)^*. **Cardinality bijection FAILS** at every n.

**Test 2 — weak encoding as subset:**

(ℤ/3^n)^* injects into the depth-n sphere of T_3 (≅ ℙ¹(ℤ/3^n)) via r → [r : 1]. Image has 2·3^{n−1} elements ⊂ 4·3^{n−1}-element sphere.

Tao recursion preserves this subset: given r ∈ (ℤ/3^n)^*, the image [3r+1 : 2^v] = [(3r+1)·2^{−v} : 1] lies in the unit subset because 3r+1 ≡ 1 (mod 3) and 2^v is a unit mod 3^n. **Encoding is well-defined.** ✓

**Test 3 — is it depth-walk?**

Hand-computed tree distance d_{T_3}([r:1], [r':1]) at n = 6 (so depth-6 sphere of T_3):

| r | Tao(r) = r' | v_3(r − r') | tree distance at n=6 |
|---|------|------|------|
| 5 | 1 | 0 | 12 |
| 7 | 11 | 0 | 12 |
| 11 | 17 | 1 | 10 |
| 17 | 13 | 0 | 12 |
| 27 | 41 | 0 | 12 |
| 35 | 53 | 2 | 8 |

The tree distance at depth-n is dominated by the maximal value 2n (typical case v_3(r − r') = 0). The dynamics is **NOT** a depth-walk — it's near-maximal lateral motion at depth n. No descent.

Even worse: the "depth" of r in T_3 is fixed by the choice of resolution n; Tao on (ℤ/3^n)^* maps depth-n → depth-n, not depth-n → depth-(n+1) or depth-(n−1).

**Bass-Serre gap (Phase 0 Gap A):** A *true* arithmetic quotient by Γ_0(3) or similar would require checking the stabilizer structure of Γ in PGL_2(ℚ_3). Γ_0(3) is not cocompact, so Lubotzky 2013 §1.2 does not apply directly. We cannot verify "the (ℤ/3^n)^*-bijection at depth n" rigorously without Bass-Serre machinery we don't have.

**Pre-registered hypothesis H_BT_C_QUOTIENT:** FALSIFIED (in its strong form), PARTIAL (in weak encoding form).

The strong claim "Tao recursion is depth-walk descending toward root" is FALSE: dynamics is lateral at fixed depth, not depth-monotone.

The weak claim "(ℤ/3^n)^* encodes into T_3 and Tao is conjugate to a well-defined map on the encoded subset" is TRUE but **vacuous as a structural insight** — it merely restates that Tao is a self-map of (ℤ/3^n)^*, which we already knew. No new T_3 geometry is gained.

**Disposition: FAIL (strong form) / PARTIAL (weak form, but vacuous).**

---

## Gate verdict for Phase 1

**All three candidates FAIL their pre-registered hypotheses.**

- A: FAIL (no single matrix lift; only a v-indexed family)
- B: FAIL (no concentrated harmonic measure on 1-cusp; trajectories mix uniformly on 3-adic residues)
- C: FAIL strong / PARTIAL vacuous (no descent dynamics; encoding tautologically re-describes Tao)

**Gate triggers H_BT_NONE_FIT.**

---

## What the Tao recursion actually looks like in PGL_2(ℚ_3)

Salvageable structural picture (not enough for the probe to proceed, but worth documenting):

- The family {M_v : v ≥ 1} ⊂ PGL_2(ℚ_3) consists of **hyperbolic isometries of T_3** with eigenvalues {3, 2^v}, translation length 1 each.
- All M_v share the boundary fixed point ∞ ∈ ∂T_3 (the eigenvector (1, 0)).
- The second fixed point of M_v is 1/(2^v − 3) ∈ ∂T_3, varying with v.
- The Collatz attractor "1" is the second fixed point of M_2 specifically (since 2² − 3 = 1).
- All M_v's translation axes pass through ∞; they form a *pencil* of axes through a common boundary point.

This pencil-of-axes structure is the **only** clean T_3 picture extractable from Tao. It does not yield a billiard, a Coxeter group, or a cusp-attractor monotonic dynamics. It is consistent with: Tao is a stochastic "step along axis selected by v_n" walk in PGL_2(ℚ_3) — but the harmonic-measure analysis (Candidate B) shows this walk does not have the 1-attractor property in the 3-adic topology.

The 3-adic / Bruhat-Tits geometry simply does not see the Collatz attractor.

---

## Cross-check: consistency with T_lead 43/45 algebraic input

T_lead at 43/45 on (1, 4) is a within-level Q-rational eigenvalue derived from cross-frequency W_+(g) weights. T_lead is a 2×2 operator on the class-resolved moment space (P_+, P_-), **not** a PGL_2(ℚ_3) element. There is no direct T_lead-to-PGL_2(ℚ_3) bridge in our setup, so consistency is vacuous. T_lead's algebra lives in a different category (a tensor of moment-class spaces × cross-frequency Fourier-on-(ℤ/3^n)), not in PGL_2(ℚ_3) acting on T_3.

This *reinforces* the H_BT_NONE_FIT verdict: Tao's natural algebraic home (cross-frequency Fourier on (ℤ/3^n)^*) is a different category from PGL_2(ℚ_3)'s natural home (3-adic places of a global number field acting on a totally-disconnected boundary).

**Same category-of-object barrier that closed Clusters 1, 2, and the 5-probe Fourier-decay arc.**
