# Probe R15 — the endpoint split and the bulk correlation — **A gate PASS, B forced, C/D measured**

**Date:** 2026-07-21  Reuses R7/R9/R10; exact where marked. Probe `probes/probe_endpoint_R15.py`. Gates Wilson's
corrected retarget (R14): the constant splits into a bulk correlation of the deviation field plus a measure-free
endpoint/sampling term.

## R15-A — ENDPOINT IDENTITY (forced, exact): **GATE PASS**
**Λ_r = S_r·b_r + Λ_r^unif** exact, r = 2…7, with Λ_r^unif computed **independently** from the R10-C trace
(measure-free) and the bulk S_r·b_r = Λ_r − Λ_r^unif:

Λ_r^unif = S_r · [3^r/(4^{3^r}−1) − 3^{r−1}/(4^{3^{r−1}}−1)] / (2·3^{r−1}).

Anchors byte-exact: **Λ₂^unif = −110/29127**, **S₂·b₂ = −240/67963**, sum = **−1490/203889 = Λ₂** ✓. Sum ==
Λ_r verified all r=2…7. Walk-back #34 not incurred; the endpoint decomposition is exact.

## R15-B — FIXED-POINT CONCENTRATION (forced): **no self-conjugate primitive angle (odd lattice)**
Direct check, r=2…5: **neither x=0 nor x=½ is a primitive angle.** x=0 is the trivial character (order 1, excluded
from the order-3^r primitives); x=½ ↔ k=3^r/2 is **not an integer** (N=3^r odd). Every primitive angle pairs
(k, N−k) with **no fixed point**. Therefore:

> **Λ_r^unif is the primitive-SAMPLING RESIDUAL** — S_r·(mean of Re w over primitive angles) = the R10-C trace
> 3^r/(4^{3^r}−1) − 3^{r−1}/(4^{3^{r−1}}−1) — **not a literal endpoint atom.** The "both ends" picture requires the
> near-{0,½} *neighborhoods*, not atoms at the fixed points.

(Numeric cross-check Σ_{k prim}Re w vs the trace matches to float precision; at r≥4 both are below the ~1e-15 float
floor while the exact trace is ~1e-47 — the exact identity is the R15-A rational gate, not the float here.) This
confirms R14's caveat as a forced gate — Wilson anticipated the branch (R15-B: "if x=½ absent, needs the
neighborhood"), and x=½ is indeed absent.

## R15-C — BULK CORRELATION LEDGER (measurement): **the −1/10 is ~99% uniform baseline**
Layer 1 is entirely uniform: Λ₁^unif = −2/21 = Λ₁, so **bulk₁ = 0**. The bulk sequence S_r·b_r = Λ_r − Λ_r^unif:

| r | S_r·b_r | float | running Σ_{r≥2} S_r·b_r |
|---|---|---|---|
| 2 | −240/67963 | −3.531e−3 | −0.00353133 |
| 3 | (exact) | +1.321e−3 | −0.00221059 |
| 4 | (exact) | +6.503e−4 | −0.00156033 |
| 5 | (exact) | +3.269e−4 | −0.00123341 |
| 6 | (exact) | −3.387e−4 | −0.00157208 |
| 7 | (exact) | +2.148e−4 | −0.00135723 |

**Target: Σ_{r≥2} S_r·b_r = −1/10 − Σ_{r≥1}Λ_r^unif ≈ −9.845e−4** (exact; ΣΛ_r^unif = −2/21 − 110/29127 + ~−8.8e−7
= −0.0990155). Signs of S_r·b_r: −,+,+,+,−,+ (r=2…7), the same oscillation as Λ_r (bulk ≈ Λ for r≥3 where Λ^unif is
dead). The bulk sequence b_r, not Λ_r, is the clean object whose sum is the theorem.

**⚠️ CORRECTION (Wilson, post-R15) — the "99% uniform baseline" is NOT progress; it localizes the difficulty, it
does not reduce it.** The Σ Λ_r^unif ≈ −0.099 that makes up 99% of −1/10 is precisely the part the corpus **already
has derivations for** — S₁ = 2/3 (four routes) and S₂ = 10/21 — so nothing new is explained by it. And the residual
deviation-field sum is *worse* conditioned, not better:
- r=2 bulk alone: S₂·b₂ = −240/67963 = **−3.531e−3**
- total bulk: **−9.845e−4**  ⟹  **Σ_{r≥3} bulk = +2.547e−3**

The tail (r≥3) is **2.59× the answer in magnitude and opposite in sign, cancelling 72% of the r=2 term.** A small
sum from a large sign-alternating series is the **hardest** convergence to establish, not the easiest — no
truncation argument can work. The reframe **relocates the entire difficulty into 1% of the target**; it is not
"99% solved." This is the operative fact, not the baseline size.

## R15-D — TWO-BAND DECOMPOSITION (measurement, NO fit): **the clean 5/3 story does NOT hold**
Split b_r at the Re w zeros (cos2πx=1/4, x=arccos(1/4)/2π=0.2088): a_r = Re w>0 band (near x≈0/1), c_r = Re w<0
band magnitude (near x≈½):

| r | a_r (Re w>0) | c_r (Re w<0) | b_r=a_r−c_r | a_r/c_r | vs 5/3 |
|---|---|---|---|---|---|
| 2 | −1.270e−2 | −9.169e−3 | −3.531e−3 | 1.385 | < |
| 3 | −6.665e−4 | −1.987e−3 | +1.321e−3 | 0.335 | < |
| 4 | +2.938e−4 | −3.564e−4 | +6.503e−4 | −0.824 | < |
| 5 | −9.487e−4 | −1.276e−3 | +3.269e−4 | 0.744 | < |
| 6 | −1.889e−3 | −1.550e−3 | −3.387e−4 | 1.219 | < |

**The ratio a_r/c_r oscillates (1.39, 0.34, −0.82, 0.74, 1.22), never crosses 5/3, and goes negative** — because
c_r (and a_r) are not sign-definite: the band-integrated δ_r·Re w changes sign *within* each band. So Wilson's
predicted mechanism "b_r = (1/3)a_r − (1/5)c_r with a_r,c_r > 0 masses, sign flips at a_r/c_r = 5/3" does **not**
hold as stated — the two-band contributions are themselves oscillating and signed. The sign of b_r is not captured
by a clean ratio of two positive near-{0,½} masses. Reported verbatim, no law fitted (R32/#32/#33 rule in force).

## R15-E — CONVENTION + FEASIBILITY
(i) **Second gate on the anchor:** Λ₂ from the μ₂ table = −1490/203889, Λ₂^unif (trace) = −110/29127, bulk =
−240/67963 — reproduced independently of R14's arithmetic. (ii) **r=8 wall:** exact b₈ needs Λ₈ = (ε₉−ε₈)/2 but ε₉
is float-only (k≥9 wall), and the D-band at r=8 needs the μ₈ profile (supp 4374, ~19M-Fraction autocorr). **Exact
bulk-sequence ceiling = r=7.** Deferred.

## Status
**R15: A gate PASS** (endpoint identity Λ_r = S_r·b_r + Λ_r^unif exact r=2…7, anchors byte-exact; #34 not
incurred), **B forced** (no self-conjugate primitive angle — odd lattice — so Λ^unif is a sampling residual, not an
atom; the "both ends" picture needs neighborhoods, confirming R14's caveat), **C measured** (the −1/10 is ~99% the
r=1 uniform baseline −2/21; the deviation-field bulk sum is the small exact target ≈−9.845e−4, the clean object
b_r, signs −,+,+,+,−,+), **D measured** (the two-band 5/3 story does NOT hold — a_r/c_r oscillates, never crosses
5/3, bands sign-changing). **Still owed (pen):** the rate/phase of the bulk sequence b_r whose sum is
−1/10 − ΣΛ_r^unif ≈ −9.845e−4 — now the (small, cleaner) whole theorem, via the transport operator (R13). No
fitting; exact endpoint gate, labeled numeric bulk/bands; the D non-fit reported as a non-fit.
