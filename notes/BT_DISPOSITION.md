# BT_DISPOSITION — Bruhat-Tits tree / BKL billiards probe headline

**Date:** 2026-05-12. Probe BT, top-level disposition. Working dir: C:/Collatz/.

---

## Headline

**H_BT_NONE_FIT.** Probe closed at Phase 1.

The Tao recursion on (ℤ/3^n)^* does not lift to a natural action on the Bruhat-Tits tree T_3 of PGL_2(ℚ_3) in any of the three pre-registered senses (single SL_2(ℚ_3) matrix, geometric random walk with 1-cusp harmonic measure, depth-walk on arithmetic quotient). Phase 2 (BKL billiard structure) and Phase 3 (universal bounds, including the all-N candidate 3c) were not reached.

The Bruhat-Tits / BKL route is closed as a Collatz substrate. It joins the Cluster 1 (Cochrane / Burgess-style exp sums), Cluster 2 (BMP/PSF cut-and-project), and 5-probe Fourier-decay arc as a fourth category-of-object barrier on the modern-framework-transplant route.

---

## Pre-registered probability vs. realized outcome

| Hypothesis | Pre-Phase-0 | Post-Phase-0 (updated) | Realized |
|---|---|---|---|
| H_BT_NONE_FIT | 30% | 45% | **REALIZED** |
| H_BT_STRUCTURE_PARTIAL | 30% | 30% | NOT |
| H_BT_UNIVERSAL_FAILS | 25% | 15% | NOT (would have required Phase 2 PASS) |
| H_BT_UNIVERSAL_PARTIAL | 10% | 7% | NOT |
| H_BT_ALL_N_CANDIDATE | 5% | 3% | NOT |

Phase 0 prior update was based on the literature category-gap: BKL papers (DHN, HPS, DH) live in continuous-hyperbolic-Riemannian setting; Bruhat-Tits papers (Lubotzky 2013, Casadio-Tarabusi & Picardello) live in p-adic / tree / Hecke-spectrum setting. Zero overlap in the pulled corpus. Asking whether Tao fits both simultaneously is a third-category-of-object request — same structural pattern as the prior probe failures.

---

## Phase that closed the probe

**Phase 1.** All three candidate constructions (A, B, C) failed their pre-registered hypotheses on hand-computed n = 1..6 examples.

- **A (single SL_2(ℚ_3) lift):** FAILED. Tao requires a v-indexed family {M_v = [[3, 1], [0, 2^v]]} in PGL_2(ℚ_3), not a single matrix. Each M_v is hyperbolic with eigenvalues {3, 2^v} and translation length 1 in T_3, but the family is parameterized by 2-adic data that PGL_2(ℚ_3) cannot collapse.

- **B (random walk with 1-cusp harmonic measure):** FAILED. The trajectory r_0 = 27 has v_3(r_n − 1) bouncing in {0..4}; r_n mod 9 is essentially uniform on unit residues. The walk does NOT 3-adically converge to 1 — it mixes uniformly over the 3-adic boundary directions. (Direct connection to memory `project_collatz_prefix_nonpropagation`: prefix signatures don't propagate. The 3-adic Bruhat-Tits geometry literally does not see the Collatz attractor.)

- **C (arithmetic quotient depth-walk):** FAILED strong / PARTIAL vacuous. Cardinality bijection between (ℤ/3^n)^* (size 2·3^{n−1}) and depth-n sphere of T_3 (size 4·3^{n−1}) fails by factor of 2. The weak encoding r → [r:1] is well-defined and Tao preserves it, but the dynamics is **lateral motion at fixed depth-n** (tree distances 8–12 at n=6, vs depth 6 — i.e., near-maximal lateral). Not depth-walk. No descent toward root. Bass-Serre gap (Γ_0(3) not cocompact in PGL_2(ℚ_3)) precludes rigorous Lubotzky 2013 §1.2 verification, but hand-computation closes the question independently.

Gate "≥ 1 of A/B/C must show feasibility for Phase 2" triggers H_BT_NONE_FIT directly.

---

## What was learned (substantive content, not just negative)

1. **Tao's algebraic home is not Bruhat-Tits.** The natural T_3 = PGL_2(ℚ_3)-tree is parameterized by 3-adic data; Tao's recursion is driven by 2-adic data (v_2(3r + 1)). The two valuations interact only via the "3r + 1 ≡ 1 (mod 3)" identity that strips the 3-adic prefix at every step. This is the same prefix-non-propagation phenomenon already documented.

2. **Salvageable T_3 picture for Tao:** {M_v} forms a pencil of hyperbolic axes in T_3 sharing the boundary point ∞. M_2 alone fixes the Collatz attractor 1; M_v for v ≠ 2 has its own "fake attractor" at 1/(2^v − 3). The probabilities of each v firing (P(v = 2) = 1/4 under Geom(1/2)) determine how often the walk steps along the "real" axis vs. fake axes — but this never resolves into a concentrated harmonic measure.

3. **The BKL connection was always a stretch.** Phase 0 made this explicit: DHN/HPS/DH papers live in continuous β-space with Lorentz reflections off hyperplane walls; Lubotzky / Casadio-Tarabusi live in discrete tree with Hecke / adjacency operators. The Coxeter / Weyl-chamber structure that DHN derives does not manifest on a tree — trees are 1-complexes with valence q+1, not (≥2)-dimensional hyperbolic polyhedra. Even *if* the Tao lift had worked (it didn't), forcing BKL machinery onto it would have been a separate category-jump.

4. **Mode E discipline held.** Phase 0 flagged the inherited-claim trap explicitly; Phase 1 closed cleanly on hand-computed examples without trying to import any BKL theorem.

5. **The H_C2_ENCODING_PARTIAL pattern recurs.** Just as Cluster 2 found that "support layer fires, weight layer is circular," Candidate C here finds that "encoding-into-T_3 works, but the encoded dynamics is vacuous." This is now a *recognized pattern* across the framework-transplant probes.

---

## 3 follow-up questions (if it had been positive)

N/A — probe closed negative.

## 3 follow-up questions (negative case)

1. **Is there a *different* tree on which Tao acts naturally?** The natural candidate would be the **Hensel tree of (ℤ/3^n)^*** — a 3-regular tree (branching factor 3 from depth 1 onward) where Tao IS a self-map at each depth, but the tree is NOT the Bruhat-Tits tree of any p-adic group. It's an ad hoc combinatorial tree. Worth a brief look as a "wrong tree, right object" replacement, but doesn't connect to the Lubotzky / spherical-function spectral machinery.

2. **Does the pencil-of-hyperbolic-axes structure ({M_v} sharing fixed point ∞) connect to anything?** A pencil of hyperbolic isometries through a common boundary point generates a *parabolic-like* subgroup of PGL_2(ℚ_3). The stabilizer of ∞ in PGL_2(ℚ_3) is the upper-triangular Borel subgroup B. Our M_v all sit in B. So Tao's algebraic content is "a discrete random walk in the Borel B ⊂ PGL_2(ℚ_3)." This is a known object (it's just the affine group acting on ℚ_3), and its harmonic analysis is standard (it's not the spherical-function setting because there's no maximal compact stabilizer K — B is opposite to K). Worth flagging that this re-routes back to *affine-group / Heisenberg-style* analysis, which is in scope of the Bourgain-Konyagin discrete sum-product secondary route already on the docket (per POLYNOMIAL_IN_A_LANDSCAPE.md §Forward Approach (b)).

3. **Is the "1-attractor" categorically incompatible with any p-adic / Bruhat-Tits substrate?** Strong reading: yes. The 1-attractor is an *archimedean* phenomenon — r_n stops being > 1 archimedean-ly. The 3-adic place can't see archimedean convergence; the 2-adic place is consumed driving the recursion forward. So the attractor lives at the *archimedean* place, which the p-adic Bruhat-Tits machinery is built to ignore. Any substrate that can see the attractor must include the archimedean place — i.e., must work over **adelic** or **global** geometry, not p-adic-only. This is a sharper version of the "category-of-object barrier" message and might be the strongest single takeaway from this probe.

---

## Routing implications

**Stays open per POLYNOMIAL_IN_A_LANDSCAPE Forward Approach:**

- (a) **Tauberian arc** — primary. Lives on the generating series E(z) = Σ ε_n z^n, escapes the Bruhat-Tits category trap entirely (different object).
- (b) **Bourgain-Konyagin discrete sum-product on ℤ/3^n ℤ** — secondary. Categorically correct. The {M_v} pencil structure noted above actually *re-routes* to this approach: Borel subgroup B's harmonic analysis on ℚ_3 reduces to multiplicative-energy / sum-product questions on (ℤ_3)^*.
- (c) Genuinely new technique — tertiary. The category-of-object barrier pattern is now stable enough across probes that the "no existing framework fits" reading is the leading honest summary.

**Closed by this probe:**

- Any Collatz-substrate-in-Bruhat-Tits-tree route.
- Any BKL-cosmological-billiards-as-Collatz-dynamics route.

**Not opened, but flagged:**

- Adelic / global-geometry substrate. The 1-attractor's archimedean character suggests that any substrate that *can* see the attractor must include all places of ℚ — i.e., a fully adelic setting. This is a strategic note, not a concrete probe — adelic Collatz substrates are not in any of our pulled corpora.

---

## Deliverables

In C:/Collatz/:

- `BT_PHASE0_THEOREMS.md` — Phase 0 verbatim extraction (T1–T9, gaps A–B flagged).
- `BT_CANDIDATE_CONSTRUCTIONS.md` — Phase 1 disposition for Candidates A/B/C.
- `BT_BKL_STRUCTURE.md` — not produced (Phase 2 gated out).
- `BT_UNIVERSAL_BOUNDS.md` — not produced (Phase 3 gated out).
- `BT_DISPOSITION.md` (this file) — headline.

Per the brief: no `git commit` or `push`. Nathan commits manually.
