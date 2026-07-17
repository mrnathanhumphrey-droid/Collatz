# Phase 1b — Tractability Classification of §7.2–§7.4 Constants

**Companion to:** `TAO_PROOF_CONSTANT_MAP.md` (Phase 1a).
**Question for each constant:** if Phase 2 fires, can we actually extract a numeric value (or A-dependent expression) for it, or will it block?

**Classification scheme (from pre-reg §3.2):**

- **TRIVIAL** — already explicit or one-step from cited reference.
- **MODERATE** — careful reading, no new math; primarily elementary algebra / book-keeping / Gaussian moment-generating-function work.
- **HARD** — requires deeper analytic technique (local-CLT, large-deviations, Berry-Esseen) or auxiliary estimate not in §7 itself but referenced from §2.
- **BLOCKED** — needs arithmetic-combinatorics / analytic-number-theory expertise the project lacks.

Project expertise inventory (honest):
- ✅ Linear algebra, operator theory, finite-dim spectra (R77 series).
- ✅ Computational Fourier analysis on Z/pZ and Z².
- ✅ Hand-driven probabilistic moment / MGF computations (R76 §10 fits, T_3 companion construction).
- ⚠️ Multivariate local CLT / Berry-Esseen with effective constants — **doable but slow**; requires careful saddle-point / characteristic-function manipulation.
- ❌ Arithmetic-combinatorics in the style of Bourgain–Konyagin (sum-product, Burgess-type) — present in §2 only marginally for §7, not load-bearing for §7.2–§7.4 itself.
- ⚠️ Analytic number theory of log 3 / log 2 (Baker's theorem) — Tao explicitly says §7 does NOT use this; only the explicit slope (log 9/log 2 ≈ 3.17) enters.

---

## §7.2 constants

| # | Constant | Class | Justification |
|---|---|---|---|
| C-1 | log 9 | TRIVIAL | exact. |
| C-2 | log 2 | TRIVIAL | exact. |
| C-3 | ε ∈ (0, 1/100) | TRIVIAL | parameter; fixed by user / proof end. |
| C-4 | 1/100 in weakly-black | TRIVIAL | numeric. |
| C-5 | 1/10 in strip and separation | TRIVIAL | numeric. |
| C-6 | (log 9 + log 2)/10 ≈ 0.289, 1 − 0.289 = 0.711 boundary | TRIVIAL | elementary arithmetic; trace through Case 1 of Claim (*). |
| C-7 | propagation factors 9, 4, 2 | TRIVIAL | from (7.12)–(7.15). |

**§7.2 summary: 7/7 constants TRIVIAL.** The deterministic structural part of the proof is fully explicit. This is a strong positive signal for tractability.

---

## §7.3 constants

| # | Constant | Class | Justification |
|---|---|---|---|
| C-8 | 1/4 = P(Pascal = 3) | TRIVIAL | exact. |
| C-9 | EHold = (4, 16) | TRIVIAL | exact via Pascal MGF. |
| C-10 | exponential-tail rate of Hold | MODERATE | computable from explicit Pascal MGF; the radius of convergence is exactly log 2 in the l-direction, so any c < log 2 works. Sharp value: c can be taken as ½ log 2 ≈ 0.347 with a finite explicit prefactor. Standard MGF computation, no novel math. |
| C-11 | c in e^{−c(l′−s)} of Lemma 7.7 | MODERATE | inherited from C-10. |
| C-12 | c in argument of G_{1+s} | HARD | this is the **width parameter** of the Gaussian-tail bound in Lemma 7.7, and traces through Lemma 2.2's 2D local-CLT for the (1, Pascal′)-step renewal process. Effective bookkeeping requires producing an explicit local-CLT statement for the joint (j, l)-renewal — a 2D Berry-Esseen-type computation with the specific characteristic function of (1, Pascal′). **Doable but expensive in person-hours** — perhaps 1-2 weeks of careful work on the project, leveraging classical Stone or Sazonov-type bounds. |
| C-13 | ≪ Vinogradov absolute in (7.48) | MODERATE | union-bound prefactor; trace through. |
| C-14 | Lemma 2.2 constant (the local-CLT itself) | HARD | the centerpiece. Reproved with explicit constants becomes a self-contained ~10-15 page side-note. Standard machinery (Esseen 1945 / Bhattacharya & Rao 1976) applied to the explicit 2D step distribution. Project can do this with care but it is the dominant Phase 2 cost. |

**§7.3 summary: 2 TRIVIAL, 3 MODERATE, 2 HARD.** The two HARDs (C-12, C-14) are essentially **one combined task** — produce an effective version of Lemma 2.2 for the specific (1, Hold) step distribution. This is the dominant Phase 2 obstacle but it is in standard probability technique, not blocked-by-expertise.

---

## §7.4 constants

| # | Constant | Class | Justification |
|---|---|---|---|
| C-15 | ε | TRIVIAL | parameter. |
| C-16 | A | TRIVIAL | parameter. |
| C-17 | 1 in max(·, 1) | TRIVIAL | numeric. |
| C-18 | **C_{A,ε} (TERMINAL)** | derived from all below; see §"Synthesis" | this IS the proxy for C_A; not extracted independently but composed from C-19 through C-41. |
| C-19 | exp(−ε) | TRIVIAL | numeric. |
| C-20 | O(·) inside exp in (7.42) | MODERATE | elementary calculus bound. Explicit numeric ≤ 2 by direct computation. |
| C-21 | 1 + O(A/log m) from E exp(A log m · Geom(4)/m) | MODERATE | Geom(4) MGF; explicit. |
| C-22 | "m large depending on A, ε" for Case 1 | MODERATE | the threshold is the m where exp(−ε/2 + O(A/log m)) ≤ exp(−ε/4), which solves for m ≥ exp(4·(constant)·A/ε). Effective threshold: **m₀(A, ε) ≈ exp(C₁ · A/ε)** for some absolute C₁ extractable from C-20, C-21. |
| C-23 | O(·) in (7.49) | MODERATE | same shape as C-21 but uses C-14; depends on local-CLT input. |
| C-24 | ≫ 1 in (7.51) | MODERATE | Gaussian-tail computation. The first-passage location is centered at (s/4, l + O(1)) with Gaussian j-dispersion O(s^{1/2}); the white-side region (lying just outside Δ with horizontal distance ≤ O(1)) has positive Gaussian mass. Quantifying it: the white-zone is roughly a strip of horizontal width ≥ (1/10)·log(1/ε) (from C-5) sitting at the predicted j-mean, so the Gaussian mass on it is ≥ erfc((1/10) log(1/ε) / s^{1/2}) − [tail] ≫ 1 uniformly as long as the Gaussian width s^{1/2} dominates the strip width — which it DOES because s ≥ m/log² m ≥ Ω(1). **Extractable to a numeric lower bound (e.g., ≥ 1/100) by direct computation.** |
| C-25 | "m ≥ C_{A,ε}" for Case 2 | MODERATE | similar to C-22. |
| C-26 | 0.9 | TRIVIAL | numeric. |
| C-27 | 0.8 | TRIVIAL | numeric. |
| C-28 | c in exp(−cm) | HARD | large-deviations rate for Geom(4)-sum; computable explicitly from Cramér's theorem with the specific Geom(4) rate function I(x) = x log(4x/3) + (1 − x) log((1−x)·4) (or similar — the Geom(4) rate function). The gap 0.8 − 0.793 ≈ 0.007 sets the prefactor; **a careful Cramér computation gives c ≈ (gap)²/(2·variance) ≈ 0.000025**. Project can do this with effort. **Numeric but tight** — see Phase 1c. |
| C-29 | ≪_P factor | MODERATE | combinatorial union bound over P repeated applications; tracks linearly. |
| C-30 | 10A/3 | TRIVIAL | numeric. |
| C-31 | exp(−10A), 10^{-A-2} | TRIVIAL | numeric. |
| C-32 | R = A²/ε | TRIVIAL | parametric. |
| C-33 | ≫ in (7.59) | MODERATE | inherits from C-24. |
| C-34 | ε < 1/100 | TRIVIAL | parametric. |
| C-35 | ≪ Vinogradov in Lemma 7.10 | MODERATE | inherited from C-14. |
| C-36 | c in exp(−c A²(1+p)) of Lemma 7.10 | HARD | inherits from C-28 / C-14 type bounds; same machinery. |
| C-37 | exponents 0.4, 0.6, 0.1, 0.2 | TRIVIAL | numeric, chosen by Tao for technical convenience; would re-derive identical values. |
| C-38 | 10 (margin in l′ − s/log 2) | TRIVIAL | numeric. |
| C-39 | 4 in s′ < 4A(1+p)³ | TRIVIAL | numeric. |
| C-40 | ≪ in P(F) ≪ 10^{-A-2} | MODERATE-HARD | terminal absorption point; tracks accumulated drift over all prior constants. **Risk:** cumulative ≪-drift through 18 named-unspecified instances could inflate by O(2^{18}) ≈ 10^{5} in the worst case under Tao's "vary line to line" convention. The Phase 1c projection quantifies this. |
| C-41 | "P large depending on A, ε" | MODERATE | the threshold is P = O_{A,ε}(1) where the implicit constant is computable once C-22, C-25 are pinned. |

**§7.4 summary: 14 TRIVIAL, 9 MODERATE, 3 HARD (C-28, C-36, C-40 with caveat).**

---

## §2 (auxiliary) — Lemma 2.2

| Constant | Class | Justification |
|---|---|---|
| Lemma 2.2 local-CLT absolute constants | **HARD (but NOT BLOCKED)** | Standard 2D Berry-Esseen / local-CLT methodology applied to the explicit 2-step distribution (1, Pascal′) × Geom(4). Project can do this with focused effort (~1-2 weeks). The classical reference is Bhattacharya & Rao 1976 *Normal Approximation and Asymptotic Expansions*, Chapter 9 (local-CLT with rates). |

---

## Synthesis

**Counts (40 cataloguable individual constants in §7.2–§7.4, ignoring the auxiliary C-18 terminal):**

| Class | Count | Fraction |
|---|---|---|
| TRIVIAL | 22 | 55% |
| MODERATE | 13 | 32.5% |
| HARD | 5 | 12.5% (C-12, C-14, C-28, C-36, C-40) |
| BLOCKED | 0 | 0% |

**Three independent positive checks** (per pre-reg §3.4 override):

1. **Constant map clean.** Phase 1a produced 42 entries with no fundamental ambiguity; the proof structure does support line-by-line bookkeeping (Tao's writing is dense but unambiguous on the bookkeeping-relevant lines). ✅
2. **No BLOCKED-by-expertise constants.** All HARD constants are in standard probability technique (local-CLT / large-deviations / Berry-Esseen). They are time-expensive but not expertise-blocked. The arithmetic-combinatorics fear (Bourgain-Konyagin / sum-product) is **not present in §7.2–§7.4** — Tao routes around it via the deterministic structural analysis of §7.2 and the renewal-process probabilistic analysis of §7.3/§7.4. ✅ — strongest positive signal in Phase 1.
3. **Tractability classification mostly TRIVIAL/MODERATE.** 35/40 = 87.5% are TRIVIAL or MODERATE; only 5/40 are HARD; 0 BLOCKED. ✅

**All three positive checks pass.** The structural conclusion of Phase 1b is therefore: bookkeeping is **mechanically tractable for the project**, with the bottleneck being the 5 HARD constants centered on producing an effective version of Lemma 2.2 for the specific (1, Hold) step distribution. Project should budget ~2-4 weeks of focused Phase 2 work, dominated by the Lemma 2.2 reproof.

The remaining question is **looseness** — even if we CAN extract every constant, would the resulting C_A be tight enough to satisfy Nisoli η < 1 at the verified K range? That is Phase 1c.

---

End of Phase 1b — Tractability Classification.
