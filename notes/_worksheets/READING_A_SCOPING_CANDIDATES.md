# READING_A_SCOPING_CANDIDATES — Phase 1 candidate-by-candidate tractability assessment

**Date:** 2026-05-11. Wilson (analyst) reporting to Nathan. Scoping probe for R77.5 §7's three function-space candidates. NOT a full Reading A construction.

---

## Probe framing (fidelity to R77.5 §7)

R77.5 §7 names three candidates for the function-space framework that would carry the rate-1/2 phenomenon now that single-operator-at-finite-truncation framings (T_3 in R77.3, K_k within-level in R77.4, R_k inter-level in R_K probe) have all been shown structurally inadequate:

- **(A) Hilbert spaces of locally constant functions on Ẑ_3^×** — "the natural completion of ⊕_k V_k via the projective system."
- **(B) Wavelet-like frames on Ẑ_3^× adapted to the {W_k} multi-resolution decomposition.**
- **(C) Transfer-operator analysis on the action of the Syracuse map's coherent extension to Ẑ_3.**

The probe assesses tractability of each as an *entry point* — not whether each can be fully constructed.

Anchors used: `result_77_5_inter_level_residual.md` (§§1-3 and §10 are load-bearing), `result_77_5_compute_R_k.py` (constructive grounding of V_k / W_k / lift T as exact rationals), `result_77_6_generating_function.md` (independent evidence: branch cut at z=2, not pole → non-self-adjoint Φ_∞), `M3_DISPOSITION.md` and `R_K_DISPOSITION.md` (precedents — both failed by trying single operator on either finite Hilbert space or rectangular inter-level map).

---

## Candidate A — Hilbert spaces of locally constant functions on Ẑ_3^×

**R77.5 phrasing:** "the natural completion of ⊕_k V_k via the projective system."

### (a) Construction tractability — HIGH

The construction is *essentially already implicit* in R77.5. The chain

```
V_1 ⊂_T V_2 ⊂_T V_3 ⊂_T ... ⊂_T V_k ⊂_T ...
```

with isometric (up to √3) embeddings T_{k→k+1} : V_k → V_{k+1} (R77.5 §3.1) is exactly an inductive system of finite-dimensional Hilbert spaces. The colimit / direct limit ⋃_k T_k(V_k) sits inside `L²(Ẑ_3^×, μ)` where μ is the unique 3-adic Haar measure on Ẑ_3^× normalized to 1.

The orthogonal direct sum decomposition

```
L²(Ẑ_3^×, μ) = ⨁_{k≥0} W_k
```

(with W_0 := V_1 = "constants on Ẑ_3^× / nontrivial residue mod 3", treated as the "scale-0" subspace, and W_k as defined in R77.5 §3.1) is the standard scale-filtration of `L²(K)` for K a profinite abelian group. This is textbook profinite-Fourier-analysis territory: K. Iwasawa / Tate-style local-field analysis.

The functions in W_k are exactly the locally constant functions on Ẑ_3^× that are constant on cosets of 3^{k+1}Ẑ_3 ∩ Ẑ_3^× but have mean zero on cosets of 3^k Ẑ_3 ∩ Ẑ_3^× — i.e., R77.5's "mean zero on each 3-fiber" condition (§3.1) re-expressed as a 3-adic locally-constant condition.

**Verdict:** the construction is essentially done already. What's missing is *the explicit statement of the construction* — naming the Haar measure, naming the filtration, citing that this is the standard scale-filtration of L²(profinite). No new mathematics required.

### (b) Basis tractability — HIGH at small k

For each k ≥ 1, W_k is finite-dimensional with `dim W_k = N_{k+1} − N_k = 2·3^k − 2·3^{k−1} = 4·3^{k−1}` (R77.5 §4).

An explicit basis can be constructed entry-by-entry over Q via the existing R77.5 computation:
- enumerate coprime r' ∈ Z/3^{k+1};
- for each coprime r ∈ Z/3^k, the 3 lifts {r, r+3^k, r+2·3^k} span a 3-dimensional subspace of V_{k+1}; T(V_k) picks out the "all equal" direction; W_k contributes the 2 mean-zero directions.

So W_k has a natural basis of "fiber-wise mean-zero" indicator-difference vectors, indexed by (coprime r ∈ Z/3^k) × (2 mean-zero directions per fiber). All entries are rational over Q. Explicit at k = 1..5 already (uses the same machinery as `result_77_5_compute_R_k.py`).

**Verdict:** basis can be constructed at k=2..5 using project-internal tools, in fractions arithmetic, within a session. No external theory required.

### (c) φ_n articulation tractability — MEDIUM

The bilinear pair-form moment φ_n from R76 needs to be re-expressed as a locally constant function on Ẑ_3^×. R77.5 §5 sketches this:

```
ε_n = ⟨φ_n, π_n − π_∞⟩ = Σ_{k} ⟨φ_n, lift_n(R_k)⟩.
```

The functional φ_n lives at level n (it's a function on coprime states in Z/3^n, hence locally constant at scale 3^{-n} on Ẑ_3^×). Its projection onto each W_k subspace is computable over Q from R76's definition + R77.5's lift machinery.

**What's *not* yet done:** writing out the closed-form projection `proj_{W_k}(φ_n) = ?`. R77.5 §5 names this as the open question. It's a finite linear-algebra problem at each (n, k) — tractable at small k but the answer is not yet derived.

**Verdict:** medium tractability. Projecting φ_n empirically is a finite-rational computation at small (n,k). Deriving a closed form would require additional structural work (likely doable but not free).

### (d) Minimum-viable-test definability — YES

A minimum-viable test exists and is project-internal:

1. Compute basis of W_k for k=1..5 over Q (extension of `result_77_5_compute_R_k.py`).
2. Compute φ_n for n = 4, 5, 6 over Q from R76's bilinear pair-form definition (already exists in R76 anchor).
3. Compute `c_{n,k} := ⟨φ_n, lift_n(R_k)⟩` over Q for all (n, k) with k < n ≤ 6.
4. Check: does the dominant-W_k contribution `max_k |c_{n,k}|` decay like 2^{-n} as n grows? If yes, **W_k filtration carries the rate-1/2 envelope**, validating the Candidate-A framework direction. If no — if the contributions decay at different rates per k, or if the rate is not 1/2 — the W_k filtration is not the right carrier.

**Falsifiable signature:** rate-1/2 must show up either as (i) per-k geometric decay at rate 1/2, OR (ii) phase-cancellation across k's contributions producing rate-1/2 in the sum. Either pattern is detectable; absence of both falsifies.

**Verdict:** clean minimum-viable test, project-internal, one focused session.

### (e) External reading requirements — LOW

For *full construction* (Reading A in the R77.5 sense — the rigorous closure path):

- **Tate's thesis (1950) or Folland _A Course in Abstract Harmonic Analysis_ Ch. 4-6** — for the L²(profinite abelian group) framework, Haar measure, Fourier decomposition into characters.
- **Vladimirov, Volovich, Zelenov _P-adic Analysis and Mathematical Physics_** — for explicit p-adic L² basis constructions; this is the standard reference for our setting.

These are *named, citable, standard*. Not "p-adic wavelets are well-known"; specifically the L²(Ẑ_p^×) framework via Mellin transform / Tate. The project would need ~2-3 days of focused reading to lift the standard apparatus, then ~1-2 weeks to assemble it for the specific Syracuse problem.

**Verdict:** external reading is light, named, standard.

---

## Candidate B — Wavelet-like frames on Ẑ_3^× adapted to {W_k}

**R77.5 phrasing:** "Wavelet-like frames on Ẑ_3^× adapted to the {W_k} multi-resolution decomposition."

### (a) Construction tractability — MEDIUM

P-adic wavelets are a developed but specialized theory. The canonical construction is **Kozyrev (2002)** — "Wavelet analysis as a p-adic spectral analysis" (Izv. Math. 66.2), which constructs an orthonormal wavelet basis ψ_{j,n,ε} on L²(Q_p) as eigenfunctions of the Vladimirov fractional-derivative operator D^α.

For our setting on Ẑ_3^× (compact, not all of Q_3), the relevant adaptation is **Khrennikov, Shelkovich, Skopina** ("p-adic refinable functions and MRA-based wavelets", 2009) which gives p-adic multiresolution analysis (MRA) frameworks on compact subsets.

The *concept* maps cleanly to R77.5's W_k: a p-adic MRA produces scale subspaces W_j orthogonal to coarser scales, and the W_j wavelets diagonalize natural translation-and-dilation operators. R77.5's W_k *is* such an MRA — but the project hasn't certified that the abstract Kozyrev W_j on L²(Ẑ_3^×) coincides up to isomorphism with R77.5's W_k.

**Verdict:** medium tractability. Construction exists in literature; matching it to R77.5's specific W_k requires translation work the project hasn't done.

### (b) Basis tractability — HIGH (in principle), MEDIUM (in execution)

Kozyrev's wavelets ψ_{j,n,ε}(x) = p^{−j/2} χ(p^{-j} ε x) · Ω(|p^{-j} x|_p) (where χ is the standard additive character, Ω is the indicator of Z_p) form an explicit orthonormal basis. Restricted to Ẑ_3^× (kill the j=0 / r=0 cosets) and adapted to the project's normalization, this gives an explicit basis for each W_k.

**Verdict:** basis exists explicitly in literature; encoding it in fractions-arithmetic for project-internal computation at k=2..5 is doable but requires ~1-2 sessions of careful translation. NOT trivial drop-in.

### (c) φ_n articulation tractability — LOW-MEDIUM

φ_n decomposed in the Kozyrev wavelet basis is in principle computable, but the *advantage* is supposed to be that the Vladimirov operator D^α has a clean spectral decomposition in this basis — and we don't yet know whether the Syracuse-relevant transfer operator Φ_∞ has any relation to D^α. If yes (the Syracuse coherent extension respects 3-adic dilation), Φ_∞ might inherit a clean spectral form in the wavelet basis. If no, the wavelet basis is just one of many equivalent bases for L²(Ẑ_3^×) with no special advantage.

**Verdict:** advantage of Candidate B over Candidate A is *conditional* on the Vladimirov / Kozyrev structure being compatible with the Syracuse map. This compatibility is not established.

### (d) Minimum-viable-test definability — YES, but heavier

A minimum-viable test would:
1. Implement Kozyrev's ψ_{j,n,ε} on Ẑ_3^× over Q (need exact representatives of additive character values; doable but requires roots of unity in fractions arithmetic, which Python's `fractions` does not natively support — would need a `RootOfUnity(p^k)` class or use SymPy).
2. Project R_k onto the Kozyrev basis at small k.
3. Check whether the dominant components carry rate-1/2 scaling.

The test is well-defined but the basis-construction overhead is larger than Candidate A's.

**Verdict:** falsifiable, but a heavier first session.

### (e) External reading requirements — MEDIUM

- **Kozyrev S.V. (2002)** "Wavelet analysis as a p-adic spectral analysis" — Izv. Math. 66.2, arXiv:math-ph/0012019. Specific paper, specific construction.
- **Khrennikov A.Yu., Shelkovich V.M., Skopina M.A. (2009)** "p-adic refinable functions and MRA-based wavelets" — J. Approx. Theory.
- **Vladimirov V.S. (1988)** "Generalized functions over the field of p-adic numbers" — Russian Math. Surveys; standard reference for the Vladimirov fractional-derivative operator.

These are *specific named papers*, not "p-adic wavelets are well-known." Reading scope: ~1 week for Kozyrev + Khrennikov, plus ~1 week to translate to project notation. Heavier than Candidate A.

---

## Candidate C — Transfer-operator analysis on the Syracuse coherent extension to Ẑ_3

**R77.5 phrasing:** "Transfer-operator analysis on the action of the Syracuse map's coherent extension to Ẑ_3."

### (a) Construction tractability — LOW

The Syracuse map at integer level is `n ↦ (3n+1)/2^{v_2(3n+1)}` on odd integers. The "coherent extension to Ẑ_3" needs first an extension of the dynamics to a *measurable map* on a profinite or p-adic space, then a transfer operator built from that map.

Multiple non-equivalent extensions exist in the literature:
- **Tao (2019)** "Almost all Collatz orbits attain almost bounded values" — uses a Markov-process model on 2-adic / 3-adic integers, but the relevant operator is not characterized as a transfer operator on a single Hilbert space.
- **Lagarias (1985)** classical Syracuse 2-adic conjugation.
- **Various ergodic-theoretic** treatments (Allouche, Sander) on shift spaces.

None of these is "the" Syracuse coherent extension to Ẑ_3 in a form whose transfer operator on L²(Ẑ_3^×) is characterized. R77.5 names this as a candidate but does not name a specific construction.

**Verdict:** the *map* whose transfer operator we'd analyze isn't yet specified. This is conceptually the biggest gap of the three candidates.

### (b) Basis tractability — N/A until (a) resolved

Without a definite map on Ẑ_3 we have no transfer operator to find a basis for.

### (c) φ_n articulation tractability — POTENTIALLY HIGH (if (a) resolved)

If the Syracuse coherent extension Σ : Ẑ_3 → Ẑ_3 is defined and L_Σ : L²(Ẑ_3^×) → L²(Ẑ_3^×) is its transfer operator (L_Σ f(y) = Σ_{Σ(x) = y} f(x) / |Σ'(x)|_3 in suitable sense), then ε_n = ⟨φ_n, π_n − π_∞⟩ rewrites as a moment of L_Σ^n applied to some initial functional, and rate-1/2 would directly correspond to the *spectral gap of L_Σ being 1/2*.

This is the *cleanest* framework conceptually — it's the standard ergodic-theoretic / Ruelle-transfer-operator framing where rate-1/2 has an immediate operator-theoretic meaning. But it depends on (a).

### (d) Minimum-viable-test definability — UNCLEAR

Cannot articulate a minimum-viable test without first specifying the map Σ. Different candidate extensions of the Syracuse dynamics to Ẑ_3 yield different transfer operators with potentially different spectra. The probe would be testing *a particular choice of extension* — without an external anchor, the test isn't well-posed.

**Verdict:** not yet falsifiable at the scoping level. First need to fix a specific extension.

### (e) External reading requirements — HIGH and unfocused

- **Baladi V., _Positive Transfer Operators and Decay of Correlations_ (2000)** — standard reference for transfer-operator spectral theory. Long, dense.
- **Tao (2019)** — for the candidate measure-preserving Syracuse model, though not in the transfer-operator framing we'd need.
- **Hutchinson J. (1981)** "Fractals and self-similarity" — for iterated-function-system transfer operators (the IFS analog of our 3-pre-images).
- Some specific work on **p-adic dynamical systems' transfer operators** — references here are scattered; this is an active but unconsolidated literature.

Reading scope: ~2-3 weeks for the standard transfer-operator background plus an open-ended chase to find / construct the right Syracuse extension. The "named, citable" criterion is partially failed — there's no single reference that gives "the Syracuse transfer operator on L²(Ẑ_3^×)" off the shelf.

**Verdict:** highest external-reading cost; least focused.

---

## Summary table

| Criterion | (A) Locally constant fns | (B) p-adic wavelets | (C) Transfer operator |
|---|---|---|---|
| Construction tractability | HIGH (essentially done) | MEDIUM (literature exists, needs translation) | LOW (map not yet defined) |
| Basis tractability | HIGH (extension of R77.5 machinery) | MEDIUM (Kozyrev, but roots-of-unity arithmetic) | N/A until (a) resolved |
| φ_n articulation | MEDIUM (finite Q linalg) | LOW-MEDIUM (advantage conditional) | POTENTIALLY HIGH if (a) | resolved |
| Min-viable test | YES (one session) | YES, but heavier | UNCLEAR (no fixed Σ) |
| External reading | LOW (Folland / Tate / Vladimirov) | MEDIUM (Kozyrev 2002 + Khrennikov 2009) | HIGH and unfocused |
| Single-Hilbert-space framing (A2) | ✓ L²(Ẑ_3^×) is single space | ✓ L²(Ẑ_3^×) is single space | ✓ if defined, naturally single space |
| Risk of fourth-candidate drift (A1) | none — directly (A) | none — directly (B) | watch: "specific extension" choices could drift to non-R77.5 framework |

---

## Adversarial check report

**(A1) R77.5 §7 fidelity.** The three candidates above are exactly R77.5 §7's three listed candidates. The probe is NOT proposing a fourth.

**(A2) Single-Hilbert-space framing.** All three candidates respect the L²(Ẑ_3^×) framing — they're three different framings of operators on the *same* space, NOT a sequence of rectangular maps between different spaces. The R_K probe's failure mode is structurally avoided.

**(A3) External-machinery honesty.** Candidate A names Folland and Tate / Vladimirov as standard references. Candidate B names Kozyrev (2002) and Khrennikov et al. (2009) as specific citable papers. Candidate C is the most honest about its gap — the relevant Syracuse extension is *not* characterized in named literature in the form we'd need.

**(A4) Falsifiability of minimum-viable tests.** Candidate A's test has a clean signature: per-k decay rate of `⟨φ_n, lift_n(R_k)⟩` either does or does not match rate-1/2. Candidate B inherits the same falsifiability with heavier setup. Candidate C cannot yet articulate a minimum-viable test until (a) is resolved.

---

## Phase 1 conclusion

Candidate A is by every assessed dimension the most tractable entry point. Candidate B is feasible but heavier; it becomes preferred *if* Candidate A's test reveals that the W_k subspaces aren't carrying the rate, *and* the Kozyrev structure adds value over the plain W_k basis. Candidate C is conceptually the cleanest framework (rate-1/2 = spectral gap of transfer operator) but is blocked at step (a) — without a fixed Syracuse extension on Ẑ_3, we cannot articulate a falsifiable scoping test.

Phase 2 will articulate Candidate A's minimum-viable test in full detail. Phase 3 ranks all three.
