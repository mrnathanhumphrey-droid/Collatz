# READING_A_SCOPING_RECOMMENDATION — Phase 3 ranked ordering and rationale

**Date:** 2026-05-11. Wilson (analyst) reporting to Nathan. Phase-3 deliverable of the Reading A scoping probe.

---

## Top-line ranking

1. **Candidate A — Hilbert spaces of locally constant functions on Ẑ_3^×.** Most tractable entry point. Minimum-viable test specified (READING_A_SCOPING_MIN_VIABLE_TEST.md). External reading is light and standard.
2. **Candidate B — Wavelet-like frames (Kozyrev) on Ẑ_3^× adapted to {W_k}.** Becomes preferred *if Candidate A's test reveals the W_k basis is the wrong decomposition* (signature F2/F3 of Phase 2). Heavier external reading.
3. **Candidate C — Transfer-operator analysis on Syracuse coherent extension.** Conceptually the cleanest framework for rate-1/2 (= spectral gap of Φ_∞), but blocked at step (a) — the project has not specified which extension of the Syracuse map to Ẑ_3 to analyze. Highest external-reading cost; least focused literature. Last resort within the scoping ordering.

---

## Rationale (per criterion, condensed from Phase 1)

| Criterion | (A) Locally constant fns | (B) p-adic wavelets | (C) Transfer operator |
|---|---|---|---|
| Construction tractability | HIGH | MEDIUM | LOW |
| Basis tractability at k=2..5 | HIGH | MEDIUM | N/A |
| φ_n articulation | MEDIUM | LOW-MEDIUM (conditional) | POTENTIALLY HIGH (conditional on a) |
| Min-viable test definable | YES, one session | YES, heavier | NO until (a) |
| External reading | LOW (Folland / Tate / Vladimirov) | MEDIUM (Kozyrev '02, Khrennikov '09) | HIGH and scattered |
| Single Hilbert space (A2 audit) | ✓ | ✓ | ✓ if defined |
| R77.5 §7 fidelity (A1 audit) | direct | direct | direct, but drift risk |
| External-machinery honesty (A3) | clean | clean (specific papers) | clean (gap honestly named) |
| Min-viable test falsifiability (A4) | ✓ | ✓ | not yet articulable |

Candidate A wins on every axis except "φ_n articulation could be cleaner under (C) if step (a) were resolved." That's a real upside for C, but it's gated behind work the project hasn't done.

---

## Rank #1 — Candidate A: proposed full-construction project scope

**Full-construction scope** (after H_A_CONFIRMED from the minimum-viable test):

### Phase I — Standard-apparatus port (2-3 days)

- Establish notation: L²(Ẑ_3^×, μ) with μ the 3-adic Haar measure normalized to 1.
- Establish the filtration: V_k ⊂ V_{k+1} via the lift T (R77.5 already has this); W_k = V_{k+1} ⊖ T(V_k); ⨁_k W_k = L²(Ẑ_3^×) by profinite Fourier / Plancherel.
- Map R77.5's V_k / W_k to the standard L²(profinite) framework. Cite Tate's thesis / Folland for the standard apparatus.

### Phase II — Operator Φ_∞ definition (3-5 days)

- Define the projective-limit operator Φ_∞ on L²(Ẑ_3^×) such that, restricted to each V_k, it agrees with the level-k Markov transition K_k (R77.5 anchor).
- Verify Φ_∞ is well-defined on the dense subspace ⊕_k V_k and extends to a bounded operator on L²(Ẑ_3^×).
- Identify Φ_∞'s symmetry / non-symmetry structure (R77.6's branch-cut finding implies non-self-adjoint; this informs the spectral apparatus).

### Phase III — Spectral characterization of Φ_∞ (1-2 weeks)

- Compute / estimate the spectrum of Φ_∞ at small k truncations (using R77.5 / R77.4's existing K_k spectra as data).
- Identify whether rate-1/2 is encoded as: (a) an isolated eigenvalue of Φ_∞ (predicted unlikely given R77.6 branch cut); (b) a branch-cut feature in the resolvent (R77.6 evidence supports this); (c) a spectral-density feature of the absolutely continuous part of the spectrum.
- Build the Mellin / Laplace transform of ε_n's generating function in this framework.

### Phase IV — Closure attempt (1-2 weeks)

- With Φ_∞ characterized, attempt either: (a) a Nisoli-style closure on Φ_∞'s resolvent (replacing the failed M_3 calculation on T_3); or (b) a direct branch-cut / Tauberian extraction of rate-1/2 from ε_n's generating function — this latter route is conceptually cleaner if Φ_∞'s spectrum has a branch cut at z=2 (R77.6 prediction).
- Combine with the delivered bilinear bound (PATH2 + HENSEL) to attempt full c=7/45 closure.

**Total: 3-4 weeks of focused construction work.** Honest framing: this is the scope of a research subproject, not a session task.

### External reading required

- **Folland G.B., _A Course in Abstract Harmonic Analysis_, 2nd ed. (2016)**, Ch. 3-6: Haar measure, profinite groups, Plancherel theorem.
- **Tate J., "Fourier analysis in number fields and Hecke's zeta-functions"** (1950 thesis, reprinted in Cassels-Frohlich), §2: local L²-theory on Ẑ_p^×.
- **Vladimirov V.S., Volovich I.V., Zelenov E.I., _P-adic Analysis and Mathematical Physics_** (World Scientific 1994), Ch. 6: integral transforms on Ẑ_p.

Light, standard, citable. No "p-adic wavelets are well-known" framing.

### Expected deliverables

- A characterization of Φ_∞ as a bounded operator on L²(Ẑ_3^×).
- A spectral statement: rate-1/2 is encoded as [eigenvalue | branch cut | spectral density feature] of Φ_∞.
- Either: an M_3-equivalent resolvent bound for Φ_∞ that feeds into Nisoli closure; OR a direct Tauberian extraction of rate-1/2 from the generating function — replacing the Nisoli framework with a different rigor route.
- A potentially publishable result on the spectral theory of profinite Markov chains on locally compact totally disconnected groups, independent of the c=7/45 application.

---

## Rank #2 — Candidate B: when it becomes preferred

Candidate B is the **second probe** if Candidate A's minimum-viable test triggers **H_A_FALSIFIED_WRONG_BASIS** (the F2 / F3 signature in Phase 2: L²(Ẑ_3^×) is the right Hilbert space, but the W_k filtration is not the basis in which rate-1/2 is encoded).

The hand-off rationale: if per-k contributions all decay at the trivial 1/√3 cardinality rate (F2), or if k*(n) wanders erratically with n (F3), then the rate-1/2 feature is "delocalized in the W_k basis" — it might be localized in a different basis on the same Hilbert space. The Kozyrev wavelet basis is the natural candidate because:

- It diagonalizes the Vladimirov fractional-derivative operator D^α — *if* Φ_∞ has any algebraic relation to D^α (open question), the wavelet basis would be the natural diagonalizing basis.
- It's an explicit named orthonormal basis (no construction gap relative to Candidate A's W_k basis).
- It has a built-in multi-scale structure (`j` index for scale, `n` index for translation, `ε` index for direction) that's richer than R77.5's plain W_k filtration.

**If H_A_CONFIRMED on Candidate A's test, Candidate B becomes either (i) redundant or (ii) a parallel-derivation exercise.** No urgent need to probe.

### Cost relative to Candidate A

External reading: roughly 2x (Kozyrev's specific construction is ~50 pp; Khrennikov-Shelkovich-Skopina is ~60 pp; understanding the relationship between Kozyrev wavelets and the V_k / W_k filtration requires careful translation).

Compute: roughly 2-3x (basis construction needs roots-of-unity arithmetic — either SymPy or a `RootOfUnity(3^k)` class in `fractions`-like exact form).

Honest scope: ~5-6 weeks of construction work after a confirming Candidate-B scoping probe.

---

## Rank #3 — Candidate C: why it's last entry-point priority

Despite being conceptually the cleanest framework (rate-1/2 = spectral gap of a transfer operator is the textbook ergodic-theoretic framing), Candidate C is last in the scoping ordering because:

1. **Step (a) is unresolved.** The Syracuse map at integer level (`n ↦ (3n+1)/2^{v_2(3n+1)}`) does NOT have a canonical extension to Ẑ_3. Multiple non-equivalent extensions exist (Tao 2019 has one, Lagarias 1985 has a 2-adic one, ergodic-theory treatments have various shift-space extensions). The probe would have to *choose* a specific extension before constructing its transfer operator. Without that choice the test is not well-posed.

2. **R77.5 §7 names the candidate but does not specify the construction.** This is the most under-specified of the three §7 candidates, and the probe respects A1 (do not invent a fourth candidate to fill the gap).

3. **External reading is the heaviest and most scattered.** Baladi's transfer-operator book (~300 pp) plus an open-ended search for the right Syracuse extension. Reading scope ≥ 2-3 weeks just to enter the framework, before any project-specific work.

4. **Single-Hilbert-space framing OK, but the operator itself is ambiguous.** All three candidates respect A2 (avoid the multi-Hilbert-space pathology of R_K). But Candidate C's operator depends on the choice of extension, and different choices give different operators with potentially different spectra. The probe would be testing "this particular extension," not "the Syracuse transfer operator."

**When does Candidate C become preferred?** If Candidate A produces H_A_CONFIRMED with Pattern A2 (phase-cancellation envelope) but the cancellation structure suggests a clean transfer-operator interpretation — i.e., the W_k contributions look like they would simplify under a *dynamical* re-organization rather than a static-filtration one. In that case, Candidate C becomes the natural next probe.

But this is conditional and downstream. **Not the entry point.**

---

## What this probe does NOT recommend

- Do not pursue Reading A as a full construction without first running the Candidate-A minimum-viable test. The test is ~2-3 hours of compute; the full construction is ~3-4 weeks. The test gates the construction.
- Do not invent a fourth candidate framework (Hochschild cohomology of profinite Markov chains, non-commutative geometry of Ẑ_3, etc.). R77.5 §7 lists three; the probe respects the bound (A1).
- Do not abandon Candidate A in favor of "let's just try external machinery" (e.g., jumping to Kozyrev or to Ruelle transfer operators without first checking whether the W_k filtration on plain L²(Ẑ_3^×) carries the rate). Candidate A is the lowest-cost entry; Candidates B and C only become preferred under specific falsification signatures from Candidate A.

---

## Bottom line for Nathan

**Run Candidate A's minimum-viable test in a focused session (~2-3 hours).** It will either:

- **Confirm** the W_k filtration carries rate-1/2 (Pattern A1 or A2), in which case full Reading A construction becomes a well-defined 3-4 week subproject with named external reading and a clear deliverable;
- **Falsify** in the F2/F3 sense, in which case Candidate B (Kozyrev wavelets) becomes the next scoping probe;
- **Falsify** in the F1 sense (decomposition sum doesn't match ε_n), in which case R77.5's structural identity has a bug — but this is essentially impossible given R77.5's exact-Q verification of c_k = 0;
- **Inconclusive** at n=2..6 budget, in which case we know the test direction is right but need more level depth.

Any of these outcomes resolves the immediate question of where the c=7/45 spectral-completion work should go next. The probe's recommendation is to **gate the next move on Candidate A's test, not commit weeks to construction without it.**
