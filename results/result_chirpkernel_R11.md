# Probe R11 — the chirp kernel (the two-campaign merger) — **A gate PASS, C/D PASS, B honest partial, E measured**

**Date:** 2026-07-21  Numeric U-pipeline @1e-12 + exact Λ/weld (reuses R7 build_mu/mu1/cram, R10
dlog_table/autocorr_dlog/A_N/Lambda_r/layer_mass). Probe `probes/probe_chirpkernel_R11.py`. Gates Wilson's merger
session: the layer series Λ_r is **both** campaigns' endgame — **2Λ_r = S_{r+1}−S_r = ε_{r+1}−ε_r** (ε_k = S_k−7/15
= R5's deviation d_k), so the 7/15 value question (ΣΛ_r = −1/10) and the corpus's ε-rate question (ρ≈0.984,
period-9.2) are the sum and the tail of **one** sequence of frozen exact reals.

**The bridge.** β(z) := dlog₄(1+3z) is a tower-compatible bijection of ℤ/3^r; U(k,ξ) = 3^{−r}Σ_z e((kβ(z)−ξz)/3^r)
is the chirp-DFT unitary; Λ_r = ⟨μ̂, K_r μ̂⟩, K_r = U*·diag(w)·U on the primitive shell, w(k)=1/(4e(k/3^r)−1).

## R11-A — KERNEL GATE: **PASS**
| r | (i) max\|U μ̂ − ν̂(R10)\| | (ii) max\|U(k,ξ)\|, v₃(k)≠v₃(ξ) | verdict |
|---|---|---|---|
| 2 | 1.2e−15 | 6.5e−16 | ✅ transport + block-diag |
| 3 | 3.1e−15 | 2.0e−15 | ✅ |
| 4 | 7.4e−15 | 4.8e−15 | ✅ |
| 5 | 2.3e−14 | 1.7e−14 | ✅ |

(i) The chirp unitary transports the additive Fourier transform μ̂ to R10's multiplicative ν̂ exactly. (ii) U is
**shell block-diagonal** — U(k,ξ)=0 whenever v₃(k)≠v₃(ξ), machine-zero — so β's tower-compatibility is confirmed and
layer mass = shell mass is a *theorem of the bijection* (U is a shell isometry), not a numerical weld. Walk-back #31
**not** incurred; the bridge lives.

## R11-B — FLATNESS vs Th 78.3 (measurement, labeled): **honest PARTIAL match**
| r | measured \|U\| on primitive block | Th 78.3 dense-flat 3^{−r/2} | verdict |
|---|---|---|---|
| 2 | sparse: 50% support, nonzeros **all = 3^{−(r−1)/2} = 0.57735** | 0.33333 | flat-on-support ✅ / dense-flat ❌ |
| 3 | sparse: 50% support, nonzeros all = 0.33333 | 0.19245 | " |
| 4 | sparse: 50% support, nonzeros all = 0.19245 | 0.11111 | " |
| 5 | sparse: 50% support, nonzeros all = 0.11111 | 0.06415 | " |

**Reported without forcing the identification (R17):** U's nonzero entries on the primitive block ARE flat (all
equal magnitude 3^{−(r−1)/2}) — the "flat magnitude" reading holds *on the support*. But U is a **sparse** unitary
(exactly half the primitive block is machine-zero; each row has 3^{r−1} nonzeros of magnitude 3^{−(r−1)/2},
|row|²=1), **not** the dense maximally-spreading unitary the Th 78.3 reading suggested. Honest mismatch on the
"maximally-spreading" sub-claim; the true structure is *block-diagonal (by shell) × flat-on-a-half-support*. This
does not affect A/C/D (all exact); it corrects the geometric picture of U.

## R11-C — QUADRATIC-FORM PIPELINE: **PASS**
⟨μ̂, K_r μ̂⟩ = Σ_{k prim} w(k)|θ̂(k)|² = Λ_r to ≤1e-15, r=2…5 (r=2 −0.0073078979, r=3 +0.0013198640,
r=4 +0.0006502557, r=5 +0.0003269206, all matching R10's exact Λ_r). The kernel K_r = U*diag(w)U is certified as
the pipeline the endgame runs on — an **additive-side quadratic form with the chirp-conjugated weight kernel written
down explicitly.**

## R11-D — LEDGER EXTENSION + THE TWO-CAMPAIGN WELD: **exact weld PASS, extended to r=6**
2Λ_r (character side) = d_{r+1} − d_r (ε-increment) **exact, r = 1…5** — d_k = S_k − 7/15 = R5's deviation = the
corpus's ε-sequence. The value campaign and the rate campaign are welded: Λ_r **is** the ε-increment ledger.

**Exact extension (new):** built μ₇ from the renewal and computed S₇ exactly via R10's layer-mass identity (S₆
cross-check = frozen S₆ ✅), giving a **new exact ledger value Λ₆ = (S₇−S₆)/2 = −3.38666e−4.** The per-r ledger,
all exact:

| r | Λ_r | Λ_r^unif | sign |
|---|---|---|---|
| 1 | −9.52381e−02 | −9.52381e−02 | − |
| 2 | −7.30790e−03 | −3.77656e−03 | − |
| 3 | +1.31986e−03 | −8.80e−07 | + |
| 4 | +6.50256e−04 | −1.3e−17 | + |
| 5 | +3.26921e−04 | −4.0e−50 | + |
| **6** | **−3.38666e−04** | (dead) | **−** |

**Λ signs −,−,+,+,+,− :** Λ₆ flips back to negative — the first turnover after the +,+,+ run, i.e. the
period-9-in-k oscillation of the ε-sequence appearing directly in the layer ledger (2Λ₆ = d₇−d₆ < 0, the deviation
dips further below 7/15 at k=7). **The wall:** μ₇ (3⁷ states) is the feasible exact ceiling for the renewal build;
r≥7 character-side needs the cached π₁₃…₁₆ + FFT pipeline (a permutation-reindex of the ε-computation). The
historical float ε-table (k=7…16) is the *same* object in a prior convention — the exact weld through r=6 certifies
the identification; the float tail is not re-certified here.

## R11-E — ANGULAR PROFILE (measurement, NO fit, NO detector)
Accumulation Re[Σ_{k≤K} w(k)|θ̂(k)|²] as K sweeps the primitive angles (25/50/75/100%):

| r | Λ_r | K=25% | K=50% | K=75% | K=100% | top-decile \|θ̂\|² mass |
|---|---|---|---|---|---|---|
| 4 | +6.50e−4 | +1.85e−2 | +3.25e−4 | −1.69e−2 | +6.50e−4 | 19.6% |
| 5 | +3.27e−4 | +1.78e−2 | +1.64e−4 | −1.72e−2 | +3.27e−4 | 24.0% |
| 6 | −3.39e−4 | +1.72e−2 | −1.69e−4 | −1.74e−2 | −3.39e−4 | 29.4% |

**The value is carried entirely by within-layer cancellation.** The partial sum swings to ≈+0.017 over the first
quarter of angles, then to ≈−0.017 over the third quarter, the two nearly annihilating to leave the tiny Λ_r
(~1e−4) — Λ_r is the residue of a large near-cancellation between angle-halves, i.e. the **correlation between the
profile |θ̂|² and the weight phase w along the χ(4)-angles**, exactly Wilson's reading. The |θ̂|² mass itself
concentrates mildly (top-decile 20→29% as r grows), but the *value* is all cancellation, and its sign (the r=6
flip) lives in the phase alignment, not the mass. Raw material for locating period-9-in-r as an angular phenomenon;
pen adjudicates.

## Status
**R11: A gate PASS** (chirp unitary transports μ̂→ν̂ exactly and is shell block-diagonal — layer=shell mass is now
a theorem of β), **C/D PASS** (quadratic-form pipeline ⟨μ̂,K_rμ̂⟩=Λ_r certified; two-campaign weld 2Λ_r=d_{r+1}−d_r
exact, **extended to a new exact Λ₆=−3.387e−4** via μ₇), **B honest partial** (U flat-on-support but sparse — the
"maximally-spreading" reading corrected), **E measured** (Λ_r carried by angle-half cancellation; Λ₆ sign flip =
period-9 appearing). The merger stands: the 7/15 value (ΣΛ_r=−1/10) and the corpus ε-rate (ρ≈0.984, period-9.2)
are the sum and tail of one frozen sequence, now expressible as an additive quadratic form ⟨μ̂,K_rμ̂⟩ with a
chirp-conjugated kernel — the vocabulary where Tao's Prop-1.17 machinery lives (its H₂ sibling with the kernel
written down). **Still owed (pen):** Σ_{r≥1}Λ_r = −1/10 / the decorrelation rate of ⟨profile, weight-phase⟩ (=
the corpus's ρ) in closed form. No fitting; exact rationals for the ledger, labeled numeric for the pipeline/angular
measurements; the R11-B mismatch reported as a mismatch.
