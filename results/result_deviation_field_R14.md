# Probe R14 — the deviation-field retarget (ψ-resolution) — **R13-C RESOLVED: no non-uniform ψ; 2 flags**

**Date:** 2026-07-21  Read-only verification of Wilson's ψ-resolution session (follows R13). Probe
`probes/probe_deviation_field_R14.py` (reuses R7/R9/R10). Gates the analytical resolution of R13-C's open
ψ-existence decider.

## The resolution (CONFIRMED): no non-uniform ψ exists
γ_n(τ_m) is a **partial sum of a frozen series**: γ_n(τ_m) = 1 + Σ_{r≤n} A_r(m), with
**A_r(m) := C_{r+1}(m)/3 = γ_r(τ_m) − γ_{r−1}(τ_m)** (the R13-C successive differences), and A_r(m)/S_r = the m-th
Fourier coefficient of the normalized layer-r profile.

**(1) Partial-sum structure — exact.** A_r(m) = C_{r+1}(m)/3 = γ_r(τ_m) − γ_{r−1}(τ_m) verified for m∈{1,2,3,9},
r=2…5 (all OK). The A_r are exactly R13-C's diffs.

**(2) Boundedness ⟹ no non-uniform ψ — confirmed.** R9-C gave γ_n(τ_m) bounded in n. Bounded partial sums ⟹
A_r(m) has **no nonzero limit** (else partial sums ~ nL → ∞), and Cesàro means γ_n/n → 0 (measured: m=1
0.667→0.102, m=9 1.667→0.299 over n=1…7). Therefore if the normalized profile converges at all, every Fourier
coefficient A_r(m)/S_r → 0 (S_r→7/15≠0), i.e. **it converges to UNIFORM. A non-uniform limiting angular profile ψ
is ruled out.** (A_r may still oscillate without converging — R13-C's non-shrinking diffs — so the alternative is
"→ uniform or does not converge"; either way, no non-uniform shape.) **R13-C is resolved: the limiting-shape
reading is dead**, and the object is the **deviation field δ_r** := (normalized layer-r profile) − uniform, not a
limiting ψ.

## The retarget
Because ⟨Re w⟩ = 0 (continuous mean, verified), the uniform part of the profile is annihilated by the weight, so
the theorem Σ_{r≥1} Λ_r = −1/10 becomes a statement about the deviation field alone. **Exact per-r identity
(verified r=2…7):**

> **Λ_r = S_r·⟨δ_r, Re w⟩ + Λ_r^unif**,  Λ_r^unif = (S_r/M)·Σ_{k prim} Re w(k/N) = R10-D's uniform baseline.

## ⚠️ Flag A — the uniform part is not *exactly* nothing at finite r (r=2 carries 52%)
Wilson's "Λ_r = S_r⟨δ_r, Re w⟩ exactly, every r" is a **tail statement**, exact only as r→∞. The discrete
primitive-angle mean of Re w is not zero (only the continuous mean is), so Λ_r^unif ≠ 0 at finite r:

| r | Λ_r | S_r⟨δ_r,Re w⟩ | Λ_r^unif | note |
|---|---|---|---|---|
| 2 | −0.0073079 | −0.0035313 | **−3.78e−3** | Λ^unif = 52% of Λ₂ — NOT negligible |
| 3 | +0.0013199 | +0.0013207 | −8.80e−7 | negligible |
| 4 | +0.0006503 | +0.0006503 | −1.96e−17 | dead |
| 5 | +0.0003269 | +0.0003269 | −5.7e−18 | dead |
| 6 | −0.0003387 | −0.0003387 | −7.0e−18 | dead |
| 7 | +0.0002148 | +0.0002148 | −6.8e−18 | dead |

Λ_r^unif is doubly-exponentially small for r≥3 (= R10-D's uniform baseline, which collapses as
1/(4^{3^r}−1)), but **at r=2 it is over half of Λ₂**. So the clean exact target keeps the Λ^unif term:
**Σ_r Λ_r = −1/10 = Σ_r S_r⟨δ_r,Re w⟩ + Σ_r Λ_r^unif**, with ΣΛ_r^unif ≈ −3.78e−3 (essentially all from r=2). The
deviation-field-only sum Σ_r S_r⟨δ_r,Re w⟩ ≈ −0.0962, not −1/10. The retarget is correct in spirit (the object IS
δ_r) but the r=2 uniform-baseline term is real and must be carried.

## ⚠️ Flag B — the sign story fails (it is a two-band competition, not the near-0 lobe)
Wilson's "Λ_r's sign is just the sign of δ_r's near-zero lobe … one localized feature drives the entire ledger" is
**not** supported. sign(Λ_r) = sign(near-x≈0 contribution) only at r=4,6 — **False at r=3,5,7**:

| r | Λ_r | near-x≈0 (\|x\|<0.1) contrib | Re w>0 band | Re w<0 band | sign match? |
|---|---|---|---|---|---|
| 3 | +1.32e−3 | **−2.36e−3** | −6.67e−4 | +1.99e−3 | ✗ |
| 4 | +6.50e−4 | +1.56e−3 | +2.94e−4 | +3.56e−4 | ✓ |
| 5 | +3.27e−4 | **−2.28e−3** | −9.49e−4 | +1.28e−3 | ✗ |
| 6 | −3.39e−4 | −4.70e−3 | −1.89e−3 | +1.55e−3 | ✓ |
| 7 | +2.15e−4 | **−6.66e−4** | −3.99e−4 | +6.14e−4 | ✗ |

The depletion **is** near x≈0 (R13-B stands), but Re w>0 there, so a depletion (δ<0) pushes Λ *negative* — while
Λ_r is often positive, its sign coming from the **Re w<0 band (x≈½)**. The sign is a competition between the two
weight-sign bands (x≈0/1 vs x≈½), each O(1e−3) and oscillating, with the near-0 term frequently *opposing* Λ_r.
"One near-zero feature drives the ledger" oversimplifies; the full δ_r profile across both bands sets the sign.

## Note — Wilson's mid-turn ± reframe (both flags via the conjugation-fixed points), gated
Wilson reframes both flags through the ± symmetry: Re w has exactly two extrema (+1/3 at x=0, −1/5 at x=½), which
are the two fixed points of x↦−x (the conjugation ± class, R12-D); he argues (A) δ_r must concentrate its unbalanced
mass at those pinned angles, and (B) Λ_r^unif is the {0,½} endpoint-atom imbalance of the uniform profile. **The
extrema observation is exactly true**, and it correctly names *why* Re w is sign-unbalanced (peak +1/3 vs trough
−1/5). But two discrete caveats keep it a heuristic, not the mechanism:
- **Odd lattice has no self-conjugate primitive angles.** N=3^r is odd, so x=0 is the trivial character (not
  primitive) and x=½ is not a lattice point at all; every primitive angle pairs (k, N−k) with **no fixed points**.
  So Λ_r^unif is not a literal endpoint-atom sum — it is the **primitive-angle sampling residual** of ⟨Re w⟩=0
  (= the R10-C trace 3^r/(4^{3^r}−1)−3^{r−1}/(4^{3^{r−1}}−1)), largest at r=2 (only 6 angles sample the circle) and
  doubly-exp small after. Direct r=2 check: the 6 primitive angles give mean Re w = −0.0077, ×S₂ = −3.66e−3 ≈ Λ₂^unif ✓.
- **Symmetry does not force concentration at {0,½}.** A real symmetric δ_r carries mass symmetrically at *every*
  pair (x,−x), not only at fixed points; and Flag B's two-band competition shows the value collects across both the
  x≈0 band and the x≈½ band (with intermediate bulk structure, R13-B x≈0.21 enhanced). The ± reframe captures the
  *source* of the imbalance, but "one localized {0,½} feature drives the ledger" remains unsupported by the data.

## Status
**R14: R13-C RESOLVED.** Wilson's central argument is verified — γ_n(τ_m) is a partial sum of the frozen A-series
(A_r = C_{r+1}(m)/3, exact), bounded ⟹ no nonzero A_r-limit ⟹ **no non-uniform ψ exists** (limiting-shape reading
dead; the object is the deviation field δ_r). The retarget identity Λ_r = S_r⟨δ_r,Re w⟩ + Λ_r^unif holds exactly.
**Two honest flags:** (A) the uniform-baseline term Λ_r^unif is not zero at r=2 (52% of Λ₂; negligible r≥3), so the
deviation-field sum is −1/10 minus ΣΛ^unif ≈ −3.78e−3, not exactly −1/10; (B) the sign of Λ_r is a two-band
(x≈0 vs x≈½) competition, not the near-zero lobe alone — Wilson's sign mechanism fails at r=3,5,7. **Still owed
(pen):** the rate/phase/localization of the deviation field δ_r — Σ_r S_r⟨δ_r,Re w⟩ (+ the r=2 baseline) = −1/10,
now the whole theorem, via the explicit transport operator (R13). No fitting; exact partial-sum gate, labeled
numeric for the identity/sign; both flags reported as flags, not smoothed.
