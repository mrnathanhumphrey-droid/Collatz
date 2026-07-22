# Probe R16 — the transport recursion — **A gate PASS (operator route real); B/C/D delivered; two campaigns, one crux**

**Date:** 2026-07-21  Reuses R7/R10; exact where marked. Probe `probes/probe_transport_R16.py`. Follows R15's
walk-backs #34 (the 5/3 end-lobe threshold — a/c are signed, not masses) and #35 (the endpoint-atom/trapezoid
mechanism — no self-conjugate primitive angle exists). Both were pre-registered LOW with kill conditions that fired
as specced (the #33 standing rule working twice in one probe).

## R16-A — TRANSPORT RECURSION GATE (forced): **GATE PASS**
One transport step from μ_{r−1} + Geom(2) reproduces the layer-r frequency measure (dlog domain) **exactly**,
r = 2…6: `θ_r(t) = E_v E_{X~μ_{r−1}}[ 1(dlog₄(1+3·2^{−v}X) = t) ]` equals the frozen θ_r from μ_r, and the derived
Λ_r matches (measure == frozen ✓, Λ(transport) == frozen ✓, all r). Walk-back #36 not incurred.

**Consequence: b_r is not an unexplained sequence — it is the output of an operator chain whose every link is
PASS-gated:** β is the renewal (R13-D), U is certified (R11-A/R12-A), the weight D is explicit. And b₁ = 0 (R16-D),
so the deviation field has **no initial condition** — it is generated entirely by the renewal, from nothing.

**Technical crux, named honestly (Wilson):** the recursion maps between *growing* spaces — layer-r characters live
on μ_r, layer-(r−1) on μ_{r−1} — so it is **not a self-map on a fixed space**. The analytic content is therefore a
**uniform contraction estimate on the tower**, which is the *same species* as R5's open step in the qx+1 universal-
rate paper. **Two campaigns (the 7/15 constant and the qx+1 rate) converge on one obstruction** — flagged as a
possibility eight probes ago, now arrived from the other direction.

## R16-B — BULK SEQUENCE VERBATIM (measurement, NO fit): the 72% cancellation, term by term
| r | S_r·b_r | running sum | remaining to target | sign |
|---|---|---|---|---|
| 2 | −3.531e−3 | −0.00353133 | +2.547e−3 | − |
| 3 | +1.321e−3 | −0.00221059 | +1.226e−3 | + |
| 4 | +6.503e−4 | −0.00156033 | +5.759e−4 | + |
| 5 | +3.269e−4 | −0.00123341 | +2.490e−4 | + |
| 6 | −3.387e−4 | −0.00157208 | +5.876e−4 | − |
| 7 | +2.148e−4 | −0.00135723 | +3.728e−4 | + |

Target Σ_{r≥2} S_r·b_r = **−9.845e−4**. The r=2 term is −3.531e−3; Σ_{r≥3} = **+2.547e−3 = −0.72× the r=2 term**
(opposite sign, cancelling 72%). The running sum overshoots and oscillates toward the tiny target — **a small sum
from a large sign-alternating series, the hardest convergence, not the easiest.** (Confirms R15's corrected
framing: the reframe localizes the difficulty, it does not reduce it.) No rate, envelope, or period extracted.

## R16-C — q-SWEEP SIGN DISCRIMINATOR (structural only, NO period/rate)
q=3 build validated (S₁,S₂,S₃ reproduce the frozen shell exactly). Layer-mass increment signs, deepest exact r:

| q | S_k increment signs | S_∞ behaviour |
|---|---|---|
| **3** | **- - + + +** (r=1…5) | bounded, → 7/15 (oscillates around the fixed point) |
| 5 | + + + (r=1…3) | grows ~(5/3)^k (S at r=4 ≈ 3.76) |
| 7 | + + + (r=1…3) | grows ~(7/3)^k (S at r=3 ≈ 23.15) |

**Reading: the sign oscillation is a q=3 fixed-point/criticality feature.** Only q=3 has a bounded marginal fixed
point (S_k → const at the λ=½ EP), so only q=3 can *oscillate around* it; at q≠3, S_k ~ (q/3)^k grows monotonically
(all-+ increments) — there is no fixed point to deviate from, so no analogous oscillation. This is nearest the
"q=3 criticality" branch of the pre-registration (the pattern does **not** transfer to q≠3 as a halving-lattice
universal, nor shift systematically as a 3-adic mechanism). **Honest caveat — insufficient depth:** only 3 terms at
q=5,7 (the monotone growth phase); R1 noted deep band-ringing at q=5,7 by k≈8, so I cannot rule out sign structure
in the (q/3)^k-normalized deviations at depth I can't reach. The period is *not measured* here and is not claimed —
this is a structural sign statement only (per the standing steer: stop measuring the period in every
representation).

## R16-D — FORCING CHECK (forced): dim δ_r = 3^{r−1} − 1, δ₁ = 0
Layer r has 2·3^{r−1} primitive characters in 3^{r−1} conjugate pairs, so the deviation field has dimension
3^{r−1} − 1 (the pair-values minus the sum-zero constraint). Verified by count, r=1…5: dim = **0, 2, 8, 26, 80**. At
r=1: |ν̂(χ₁)|² = |ν̂(χ₂)|² exactly (conjugates), so **δ₁ = 0 — forced by conjugation, independent of the measure.**
Certifies "the deviation field has no initial condition" as computation, not argument.

## R16-E — R85 rung-1 feasibility (statement, NO run)
Unchanged from R13-E: R85 rung-1 (= this transport) n=8 extension is one dedicated Bluestein/support-pruned probe
(O(N log N) per block; the β/U tables to r=5 exist, support law prunes to ~1/6 density) — cheaper than July, not
free. Deferred.

## Status
**R16: A gate PASS** (transport recursion exact r=2…6 — b_r is the output of a fully PASS-gated operator chain, not
an unexplained sequence; #36 not incurred), **D forced** (dim δ_r = 3^{r−1}−1, δ₁ = 0 — the deviation field has no
initial condition, generated by the renewal), **B measured** (the bulk sum −9.845e−4 assembled by 72%
sign-cancellation, term by term), **C structural** (the sign oscillation is a q=3 criticality feature — q=5,7 grow
monotonically with no fixed point; insufficient depth to probe normalized q≠3 ringing; no period measured).
**The crux is now named and shared:** the transport is a tower map between growing spaces, so the theorem is a
**uniform contraction estimate on the tower — the same obstruction as R5's open qx+1 step.** Two campaigns, one
open problem. **Still owed (pen):** that contraction estimate (equivalently the rate/phase of the bulk b_r whose
sum is −9.845e−4). No fitting; exact transport/forcing gates, labeled numeric bulk, structural-only q-signs;
walk-backs #34/#35 logged, #36 not incurred.
