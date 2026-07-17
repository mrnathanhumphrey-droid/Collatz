# Result 39 (qx+1 paper) — the q=3 critical point is a genuine ORDER-2 EXCEPTIONAL POINT (non-Hermitian defective coalescence). The "critical exponent" analog is the EP order, not a fit to r_q.

**Date:** 2026-07-16. **Verdict: ★ H_EP CONFIRMED — q=3 is a defective (order-2 exceptional-point) coalescence: eigenvectors align (cos→1), left/right self-orthogonalize (biov→0, κ→∞), δ→0. q=5,7 are ordinary gaps (cos≈0, κ≈1.3). H_DELTA: δ(L) reported, no exponent fit.**

**Headline: R32 called q=3's top-two-eigenvalue merge a "Jordan block in the limit." This confirms and sharpens it — q=3 is a non-Hermitian EXCEPTIONAL POINT of order 2. As L→∞ the two amplitude-carrying eigenvalues coalesce DEFECTIVELY: their right eigenvectors become parallel (cos: 0.9983→0.99999) and each eigenvalue self-orthogonalizes (biov=|⟨l|r⟩|: 3.8e-2→3.6e-3, condition κ: 26→274). At q=5,7 the top-two modes are ORTHOGONAL (cos≈0) and well-conditioned (κ≈1.3) — an ordinary spectral gap. The EP is specific to the boundary. Since q=3 is isolated in the primes (no continuous approach ⇒ no r_q exponent fit), the universal critical datum is the EP ORDER (=2), and it is confirmed STRUCTURALLY.**

Probe: `probe_39_exceptional_point.py`. Log: `result_39_exceptional_point_log.txt`. Runtime: ~2 min.

## The question (thread 1) and why the naive version is a trap

r_q closes the gap at q=3 (r_3=1). The stat-mech instinct: extract the critical exponent — how `(1−r_q)` scales approaching criticality. But **q=3 is isolated in the primes** — the nearest is q=5, then q=7 is *farther* — so "as q→3" is 2 unrefinable points (r_5≈0.62, r_7≈0.39) and fitting an exponent is exactly the trap `r2_cannot_discriminate_monotone_fits` warns against (this arc: quantitative priors 0-for-many). The answerable, universal question is the **nature of the degeneracy at criticality**: is the coalescence *defective* (an exceptional point — the two eigenvectors collapse to one) or *semisimple* (an ordinary degeneracy — two independent eigenvectors sharing an eigenvalue)? The EP order is the universal invariant, and it is what a universality class would share (thread 3 tie-in).

## Method

`build_M` (probe_25), the gate-validated cascade operator. For the top-2 eigenvalues by |λ| (the coalescing pair at q=3), matched right (`r`) and left (`l`) eigenvectors, measure:
- `δ = |λ₁−λ₂|` (eigenvalue splitting) → 0 at an EP.
- `cos(r₁,r₂) = |⟨r₁|r₂⟩|/(‖r₁‖‖r₂‖)` → **1** (parallel) at an EP; bounded < 1 for an ordinary gap.
- `biov = |⟨l₁|r₁⟩|/(‖l₁‖‖r₁‖)` → **0** (self-orthogonality) at an EP; `κ = 1/biov` (eigenvalue condition number) → ∞.

q=3 at L=2 (dense) and L=3 (sparse top-12); q=5,7 at L=2 (sparse) as ordinary-gap controls.

## Results

| q | L | dim | \|λ₁\| | \|λ₂\| | δ | cos(r₁,r₂) | biov | κ |
|---|---|---|---|---|---|---|---|---|
| **3** | 2 | 324 | 0.34683 | 0.34392 | 2.91e-3 | 0.998344 | 3.79e-2 | 26.4 |
| **3** | 3 | 8748 | 0.333336 | 0.333236 | 9.96e-5 | **0.999986** | **3.65e-3** | **274** |
| 5 | 2 | 10000 | 0.333334 | 0.326311 | 1.36e-1 | **0.000000** | 0.771 | 1.30 |
| 7 | 2 | 21609 | 0.333334 | 0.326941 | 1.30e-1 | 0.000086 | 0.762 | 1.31 |

**q=3 (L: 2→3):** δ falls 29×, cos rises 0.9983→0.99999 (→1), biov falls 10× (→0), κ rises 26→274 (→∞). Every EP signature, monotone in L. **The two eigenvalues merge into a single defective eigenvector — an order-2 exceptional point.**

**q=5,7 (controls):** cos ≈ 0 (the top-2 amplitude modes — Perron λ₁ and the gap mode r_q·λ₁ — have *orthogonal* eigenvectors), biov ≈ 0.77 (O(1)), κ ≈ 1.3 (well-conditioned). δ = 0.13 ≈ the real gap `(1−r_q)λ₁`. **An ordinary, non-defective spectral gap — not an EP.**

## Reading

- **H_EP CONFIRMED.** The q=3 phase boundary is a genuine non-Hermitian exceptional point (order 2). This upgrades R32's "Jordan block in the limit" from a spectral-value statement to a full eigenvector-geometry statement: the defect is the *coalescence of eigenvectors*, quantified by cos→1 and the diverging condition number.
- **The critical-exponent question, answered in the right form.** There is no closed form for r_q (R28) and no continuous approach to q=3, so a numerical exponent is not extractable. But the **EP order** is the invariant that plays the role of the critical exponent — it classifies the transition type — and it is **2**, established structurally (eigenvector coalescence), not fit.
- **δ(L) scaling:** reported (29× per L), not fit — 2 L-points and doubly-exponential tower ratios (`x_j = 2^{−2·3^{j−1}}`) make a clean power law unlikely and unjustified from the data. The EP *type* is the robust deliverable, not a δ-exponent.

## Composition with R40 (universality) — a unique, isolated EP

R40 refuted the universality hypothesis: `ord_q(w)=2` does **not** create a q=3-type critical point elsewhere ((q=5,w=4) gaps weakly, X_k converges). So q=3 is *uniquely* critical. R39 says that unique critical point is an order-2 EP. Together: **q=3 is an isolated non-Hermitian exceptional point — a single defective coalescence, with no universality class of siblings.** The two results are complementary: R39 gives the critical point's *type*, R40 shows it does not *generalize*.

## Not at stake
R10–R38. This characterizes the critical point's spectral nature; the r_q gaps (R27/R32), the d=2 boundary, and L3's structure stand. If anything, it reinforces L2's "eigenvalues collide only in the limit" (the EP is a limit phenomenon; at finite L the operator is diagonalizable).

_Reporting discipline: the semisimple-vs-defective falsifier was pre-registered (if q=3's cos stayed < 0.9 / biov stayed O(1) while δ→0, it would be an ordinary degeneracy, contradicting R32 — reported either way). It confirmed the EP. No δ-exponent fit committed (isolated critical point, 2 L-points). L=4 (q=3) was dropped: its operator has ~230M nonzeros (multi-GB) and was thrashing — L=2,3 give the two-point trend and the EP structure unambiguously; the drop is disclosed, not silent. Controls q=5,7 use sparse top-12 (only the top pair is needed for the geometry)._
