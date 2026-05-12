# Polynomial-in-A Fourier-Bound Research Arc — Landscape

**Last updated:** 2026-05-12 (post-5-probe consolidation)
**Status:** Five structural negatives across modern Fourier-decay frameworks. Routing forward to Tauberian arc + Bourgain-Konyagin.

---

## Context

After the seven-probe spectral trajectory + eighth/ninth probes (T_lead at 43/45 + Nisoli closure at corrected rate), the c=7/45 closure landscape collapsed from "three independent obstructions" to a **single isolated obstruction**: a novel polynomial-in-A Fourier bound on |μ̂_n(ξ)| outside Tao 1909.03562's iterated-cubic recursion (§7.4 Case 3). Tao's bookkeeping forces C_A ≥ A^{O(A)}, blowing the Nisoli closure budget at A=3 by 18×. Under optimistic C_A = 1, the inequality |K_bil|·K^{−A}·M_3'' < 1 fires at (r=3, A=3, K=6) with product 0.679.

A 27-PDF literature corpus targeting "polynomial Fourier decay on dynamically-defined measures outside Tao's iterated-cubic method" was assembled at `C:/Users/Nate/OneDrive/Documents/polynomial_in_a/` (INDEX flags the top six candidate frameworks). The research arc below probed the corpus's primary candidates structurally.

---

## The Five-Probe Map

| Probe | Framework | Source | Disposition | Load-bearing failure |
|---|---|---|---|---|
| 1 | L²-flattening | Baker-Khalil-Sahlsten 2407.16699 | H_L2_FLATTENING_FAILS | Step 3 separation: log\|D_v\| = −v log 2 on single AP; Plancherel collapse in discrete setting |
| 2 | SL_2 / Furstenberg | DKW 2108.06006 + Hochman-Solomyak 1610.02641 + Frostman 2601.14061 | H_SL2_EMBEDDING_DOESNT_EXIST | T_lead det=0 rank-1 (G1); even constructed 2-atom lift fails T1 transfer (first moment ≠ Furstenberg ν) |
| 3 | Cocycle Dolgopyat | Algom-Rodriguez Hertz-Wang 2306.01275 | H_COCYCLE_DOLGOPYAT_LINEAR_EXCLUSION | Φ_3 = {x→(x+a)/3} affine, c(I,x) = \|I\|·log 3 constant in x, UNI fails identically |
| R1 | ARHW + Syracuse-derived smoothing | (same paper, new strategy) | H_DELTA_EXISTS_TRANSFER_BROKEN | Five candidates; C5 (T_lead two-branch IFS) clears C² + UC + UNI; transfer pincer closes |
| R2 | Drift conditions | Glynn-Zeevi + arxiv:2005.08145 + Hairer notes | H_DRIFT_EXISTS_BUT_SPECTRAL_GAP_DOESNT_TRANSFER | Spectral gap bounds *transient* mixing P^k → π_n; three independent transfer failures (object / non-negativity / no uniform-in-n gap) |

Five pre-registered hypotheses, five honest pre-registered most-likely outcomes confirmed. No talked-up candidates, no rushed dispositions. Probes 3, R1, R2 had full verbatim-quotation discipline (pypdf 6.10.2 working method); Probes 1, 2 caveated at abstract+intro level pending re-runs.

Deliverables (all `C:/Collatz/`): [L2_FLATTENING_DISPOSITION.md](L2_FLATTENING_DISPOSITION.md), [SL2_EMBEDDING_DISPOSITION.md](SL2_EMBEDDING_DISPOSITION.md), [COCYCLE_DOLGOPYAT_DISPOSITION.md](COCYCLE_DOLGOPYAT_DISPOSITION.md), [SMOOTHING_SCOPING_DISPOSITION.md](SMOOTHING_SCOPING_DISPOSITION.md), [DRIFT_SCOPING_DISPOSITION.md](DRIFT_SCOPING_DISPOSITION.md).

---

## Three Killer Structural Findings (paper-worthy)

### 1. Plancherel collapse in discrete L²-flattening (Probe 1, §2.2 of L2_FLATTENING_VERIFICATION)

‖π_n * π_n‖₂² = (1/3^n)·Σ|π̂_n(ξ)|⁴ on Z/3^n Z. In BKS's continuous setting the averaging–flattening–separation triad is non-redundant because flattening is applied to the *derivative-cocycle pushforward*, a different object than μ. In the discrete setting without smooth cocycle, the three strategy-steps collapse and the would-be subordinate flattening estimate becomes the target Fourier bound. **L²-flattening as a strategy is structurally equivalent to the target in discrete settings.**

*Caveat:* the agent's specific Cauchy-Schwarz implication step is a heuristic equivalence-of-difficulty observation; tightening to a rigorous obstruction theorem requires another pass with full PDF discipline (Probe 1 ran before pypdf was triaged).

### 2. The ARHW pincer (Probe R1, Candidate C5)

C5 = Φ' = {φ_+, φ_−, φ_2} two-branch IFS derived from T_lead's (1,4) eigenvector. C² + uniform contraction + UNI all PASSED — the first Syracuse-encoded smoothing to clear ARHW's entry hypotheses non-vacuously.

Then the transfer gate revealed a structural pincer: transfer ν → μ_n requires a C² conjugacy h: Φ' → Φ_3 (base-3 affine) via base-3 expansion κ. But **if h exists, then Φ' is C²-conjugate-to-linear, which is exactly what ARHW Theorem 1.1's non-linearity hypothesis forbids.** The conjugacy that gives the transfer is the conjugacy that breaks the framework's applicability. Not ad-hoc-ness — structural impossibility.

### 3. The object-category mismatch (all five probes)

Each framework delivers polynomial Fourier decay for an object Syracuse μ_n is not:

- BKS → Patterson-Sullivan / Gibbs / smooth IFS measures
- Furstenberg → P¹-stationary of continuous SL_2(ℝ) random walks
- ARHW → self-conformal measures of C² IFS (smooth-dynamical)
- ARHW + smoothing → requires C² conjugacy that itself violates ARHW's hypothesis
- Drift conditions → time-to-stationarity rates of P^k → π_n (transient, not stationary)

Two especially load-bearing common failure modes recur across Probes 2, 3, R1: **real Fourier characters χ_q on ℝ don't restrict to 3-adic characters on Z_3 under base-3 expansion κ**. Even when a candidate lift gets polynomial decay on a real measure ν, the decay is of the wrong object. Smoothing the source doesn't bridge this — the gap is in the **Fourier-analytic category** of the target.

---

## The Unifying Meta-Pattern

**All five probes fail via category-of-object mismatch.** The polynomial-in-A bound on |μ̂_n(ξ)| is a *discrete-arithmetic Fourier-decay statement about the stationary distribution itself*. Modern Fourier-decay frameworks were built for objects π_n is structurally not:

- Continuous/smooth-dynamical target (Probes 1, 2, 3, R1)
- Transient-mixing target (Probe R2)

Discrete-arithmetic Markov-chain stationary measures do not fit either category. The polynomial-in-A unblocker therefore requires either:

(a) A **different object** — generating-series / Plancherel-trace Tauberian, escapes the target-category trap
(b) A **different category** — discrete-arithmetic sum-product, natively correct
(c) Genuinely **new technique** outside existing literature

This is itself a paper-worthy no-go boundary. The five-probe arc is a clean structural negative on transplanting modern Fourier-decay machinery to Syracuse stationary measures.

---

## Forward Approach

### (a) Tauberian arc — PRIMARY

Already opened per memory (`project_collatz_r78_bilinear_cracked.md`). Framework: Flajolet-Sedgewick Ch. VI singularity analysis + Chevalier 2507.15394 Thm 1.16 (meromorphic h with pole of order M at 0 → coefficient n^{M−3/2}).

Operates on generating series E(z) = Σ ε_n z^n / Plancherel-trace Σ_k|π̂_k(ξ)|², **not on the dynamical operator T_lead**. Different object → escapes the category-of-object trap that closed all five framework probes.

Literature bundle: `C:/Users/Nate/Documents/burgess/literature/tauberian/` (Hank-curated INDEX, 18 PDFs).

**Pending blocker:** ε_7 exact-rational compute extended to k=8. R77.7 v2 solver redesign needed (k=8 fired earlier and killed at 2.8 GB memory blowup with stdout-buffer deadlock; parallel-primes / sparse-mod-p / Wiedemann's algorithm are candidate redesigns). Single-theorem selection (Chevalier Thm 1.16 vs siblings) cannot complete without it.

### (b) Bourgain-Konyagin discrete sum-product on Z/3^n Z — SECONDARY

Bundle already at `C:/Collatz/Bourgain-Konyagin/`. **Categorically correct** — discrete-arithmetic question → discrete-arithmetic framework, no category mismatch. The 2-is-primitive-root-mod-3^n fact lives in (Z/3^n Z)* multiplicative structure, which is exactly Bourgain-Konyagin's setting. Technically heavier than the corpus's modern Fourier-decay frameworks but doesn't carry their target-category structural barrier.

Specific candidate: discrete sum-product on the orbit {2^v mod 3^n : v ≥ 0} gives multiplicative-energy bounds that may translate to polynomial-in-A character-sum control on Σ_x μ_n(x) χ(x).

### (c) Genuinely new technique — TERTIARY / strategic note

The five-probe arc demonstrates the modern Fourier-decay literature targets the wrong category for Syracuse. The polynomial-in-A unblocker may require a new discrete-arithmetic Fourier-decay technique that doesn't yet exist as a packaged theorem. This is the high-risk/high-reward branch and is the natural fallback if Tauberian + Bourgain-Konyagin both close negative.

---

## Big-Picture View

**Where this sits in c=7/45 closure landscape:**

| Component | 2026-05-11 | 2026-05-12 (post-5-probe) |
|---|---|---|
| Operator-theoretic obstruction (rate-specific) | Open | LIFTED — T_lead at 43/45, eigenvector (1,4), M_3'' = 24.426 exact |
| Bilinear bound \|K\| | Delivered (r ≤ 3: 2√N strict; r ≥ 4: 2√p·√N polylog-free via Hensel) | Same |
| Tao Prop 1.17 effective C_A | Open | Open, **single remaining obstruction** |
| Modern Fourier-decay literature transplant | Untested | NO-GO across 5 framework families (this document) |
| Tauberian arc | Opened | Live, gated on ε_7→ε_8 compute |
| Discrete sum-product (Bourgain-Konyagin) | Untouched | Pre-cleared categorically, technically open |

**What the polynomial-in-A unblocker still requires:**

A bound |μ̂_n(ξ)| ≤ C(n) · A^{−γ} for some γ > 0, with C(n) controlled, that:
1. Lives in the Fourier-analytic category of π_n on Z/3^n Z (not on real or smooth-IFS measures)
2. Has C_A polynomial-in-A dependence (not Tao's iterated-cubic A^{O(A)})
3. Fires at A ≥ 3 to close Nisoli at (r=3, K=6)

**Strategic state:** the five-probe arc consolidates a definite no-go boundary on one route while opening / pre-clearing two alternative routes. The polynomial-in-A unblocker is harder than literature-transplant but the search space has been mapped.

---

## Cross-references

- Memory: `project_collatz_r78_bilinear_cracked.md`, `project_collatz_prefix_nonpropagation.md`, `project_collatz_move2_qx1_family.md`
- Session disposition: `C:/Users/Nate/Documents/burgess/SESSION_DISPOSITIONS_2026_05_12.md`
- Tao 1909.03562 §7.4 (the Prop 1.17 obstruction): in corpus at `polynomial_in_a/pdfs/`
- ARHW 2306.01275 full extract: `C:/tmp/arhw_full.txt`
- Probe deliverables: `C:/Collatz/{L2_FLATTENING,SL2_EMBEDDING,SL2_FRAMEWORK,COCYCLE_DOLGOPYAT,SMOOTHING_SCOPING,DRIFT_SCOPING}_*.md`
- Bundle root: `C:/Users/Nate/OneDrive/Documents/polynomial_in_a/`
