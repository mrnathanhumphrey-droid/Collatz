# Result — PROBE E: the compressed (e_ρ,γ)-chain, frozen + L-series to L=4. Perron → 1/3 monotone (but q-generic); the compression error DOMINATES the true gap — the uniform-average convention is too crude to resolve the coalescence.

**Date:** 2026-07-16. Instrument work: freeze the compressed chain's averaging convention (E1), run the L-series L=2,3,4 (E2, the compression's payoff — L=4 reachable where the 236k operator walls), emit the entry decomposition (E3), q=7 control (E4). No proof, no rate fit. Probe `probes/probe_phase2b_E.py`, log `logs/probe_phase2b_E_log.txt`, dumps `outputs/compressed_q3_L{2,3}_frozen.tsv` + `_decomp.tsv`.

**Headline: the compressed Perron → 1/3 monotone (dist 1.06e-2 → 9.78e-4 → 1.18e-4 at L=2,3,4) — but this is q-GENERIC (it is c₀=Σw², which → 1/3 for every q). The q=3-special signal is the compressed GAP (partner vs c₀), which shrinks at q=3 (2.9e-3 → 1.18e-4) but stays OPEN at q=7 (0.174). ⚠️ CRITICAL CAVEAT: the compression error (1.7e-2 → 3.2e-3) is LARGER than the true gap (1e-4 at L=3), so the compressed partner lands on the WRONG SIDE of c₀ and the compressed gap does NOT resolve the true coalescence. The uniform-average convention is too crude for the two-limit program as posed.**

## E1 — FROZEN AVERAGING CONVENTION (the blocker, documented)
`Lmat[src_class, dst_class] = (1/|src_class|) · Σ_{src∈class} Σ_{dst∈class} M[dst, src]`, with `|class| = D`.
- **Source-side UNIFORM average** of outgoing weights: each source state in a class is weighted `1/D`.
- Equivalently the Galerkin projection `P Mᵀ P` with `P` = orthogonal projection onto (e_ρ,γ)-class indicator functions **under the uniform (counting) inner product** on each class. Left-action; eigenvalues approximate M-eigenvalues whose **left** eigenvectors are (e_ρ,γ)-class functions.
- **The single free choice, NAMED: uniform source weighting.** Not stationary-weighted, not right-eigenvector-weighted. (This choice is what the entry decomposition E3 derives against, and — see the caveat — is likely too crude: a stationary- or partner-weighted inner product is the natural next lever.)
- Class index `idx = e_ρ·q^L + γ`, `e_ρ ∈ ℤ/D`, `γ ∈ ℤ/q^L`. Verified: `Lmat[(0,0)→(0,0)] = 65/189 = Σw²` (matches Probe C exactly).

## E2 — L-series of the compressed chain (q=3)
| L | dim | Perron | dist to 1/3 | c₀ | compressed-partner | compressed gap | compression error (rel to true) |
|---|---|---|---|---|---|---|---|
| 2 | 54 | 0.343915 (=c₀) | 1.06e-2 | 0.343915 | 0.341016 | 2.90e-3 | 1.7e-2 |
| 3 | 486 | 0.334312 (=partner) | 9.78e-4 | 0.333336 | 0.334312 | 9.76e-4 | 3.2e-3 |
| 4 | 4374 | 0.333451 (=partner) | 1.18e-4 | 0.333333 | 0.333451 | 1.18e-4 | (true partner = G) |

- **(a) Perron → 1/3, monotone. ✅ (pre-registered direction).** BUT this is **q-generic**: the Perron is c₀ (L=2) then the compressed-partner (L≥3), and both → 1/3 because `c₀ = Σw² → 1/3` for every q (the "3" is q-blind, R5). **Perron→1/3 is NOT by itself the q=3 signal** — exactly the q-generic trap the pre-registration warned about.
- **(b) top-6, family-flagged:** only **c₀** appears as a family member in the compressed spectrum (`*` in the log) — expected, since only k=0 (ℓ₀, no twist) is gauge-invariant / in the (e_ρ,γ) subspace. The compressed-partner is non-family (c₀-masquerade criterion applied).
- **(c) compressed gap (partner vs c₀):** 2.90e-3 → 9.76e-4 → 1.18e-4 (the **derived-side** sequence; the true 2.9e-3, 1.0e-4 stays untouched/unfitted). It decreases — but see the caveat.
- **(d) compression error:** 1.7e-2 (L=2) → 3.2e-3 (L=3), decreasing (~0.19×). **L=4 true partner not computed (that's G).**

**⚠️ CRITICAL CAVEAT (deviation, reported).** At L=3 the compression error (3.2e-3) is **>> the true gap (1.0e-4)**. Concretely: the true partner is 0.333236 (**below** c₀); the compressed-partner is 0.334312 (**above** c₀) — wrong side, off by 3.2e-3. So the compressed gap (9.76e-4) is **dominated by compression error, not the true coalescence** — it does NOT resolve the true partner–c₀ approach. **The two-limit program as posed (compressed Perron→1/3 AND compression error→0) does not close with the uniform-average convention**: the error must vanish *faster* than the true gap, and it does not (3.2e-3 vs 1e-4 at L=3). The natural next lever is E1's named free choice — a stationary- or partner-weighted inner product — or the true partner directly (G).

## E3 — entry decomposition (the derivation's judge)
Dumped `outputs/compressed_q3_L{2,3}_decomp.tsv`: for each nonzero `Lmat[src,dst]`, the list of `(δa,δb)` move-pairs and their w-products. Example (the Δ self-loop, verifying the gate algebra):
```
(0,0)→(0,0)  = 0.34391534 = Σw²   via (1,1):0.258 (2,2):0.0645 (3,3):0.0161 (4,4):0.00403 (5,5):0.00101 (6,6):0.00025   [the equal-move (δ,δ) pairs]
(0,0)→(2,1)  = 0.0529 = 10/189·... via (1,5),(2,6),(3,1),(4,2),(5,3),(6,4)   [same-parity leaks]
```
So the self-loop is exactly the diagonal `Σ_δ w_δ²` (c₀) and the leaks are the same-parity move-pairs — the algebra Claude derives entry-by-entry. Common-denominator (189-family) exact-rational matrices: `compressed_q3_L{2,3}_frozen.tsv` (L=2 exact; L=3 rationalized).

## E4 — q=7 control: the gap stays OPEN (q=3 is special)
| q=7 L | dim | Perron | dist to 1/3 | compressed-partner | gap (partner, c₀) |
|---|---|---|---|---|---|
| 2 | 1029 | 0.333334 (=c₀) | 3.18e-7 | 0.15890 | **0.174** (rel 3.1e-3) |
| 3 | 50421 | — | — | — | **DENSE EIG WALLS (~20GB) — reported** |

- **✅ Pre-registration confirmed (L=2):** at q=7 the compressed **gap stays open (0.174 ≈ r₇)** — the partner (0.159) sits far below c₀ (0.333), no coalescence. **Contrast:** q=3 gap shrinks (2.9e-3 → 1.18e-4), q=7 gap open (0.174). **The q=3 approach IS special against the control** — the compressed gap, not the Perron, carries the specialness (the Perron → 1/3 is q-generic for both).
- q=7 L=3 (dim 50421) walls on dense eig; L=2 stands. No extrapolation.

## Adjudication
| item | verdict |
|---|---|
| E1 | convention frozen: source-side uniform average (the named free choice). |
| E2a | Perron → 1/3 monotone ✅ — but q-GENERIC (it is c₀→1/3), not the q=3 signal. |
| E2c/d | ⚠️ compression error (3.2e-3) DOMINATES the true gap (1e-4) at L=3 → the compressed gap does NOT resolve the true coalescence; two-limit program doesn't close with the uniform convention. |
| E3 | entry decomposition dumped; self-loop = Σw² verified; the algebra is localizable per move-pair. |
| E4 | q=7 gap stays OPEN (0.174) ⇒ q=3 gap-shrinkage is special ✅; q=7 L=3 walls. |

**Redirect:** the compressed chain reproduces c₀ and the Perron→1/3 cleanly, and confirms the q=3-vs-q=7 gap contrast — but the **uniform-average convention is too crude to resolve the true coalescence rate** (error ≫ true gap). The derivation should either (i) reweight the inner product (E1's free choice — stationary/partner-weighted) so the compressed partner lands on the correct side of c₀, or (ii) target the true partner via a direct L=4 solve (G). The entry decomposition (E3) and the frozen convention (E1) stand as the handoff for deriving the c₀-row and the leak structure; the *rate* needs a better-conditioned compression than uniform.

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, H, J, L-A, L-B, Real-T1, R1-STOP, Probe P (P1/P3/P4), Probe C. No `r_q` value changes; **no rate-law fit** (2.9e-3, 1.0e-4 untouched, and explicitly NOT identified with the compressed gap).

_Reporting discipline: E2a's "Perron→1/3" is reported as CONFIRMED but flagged q-generic (not the q=3 signal), per the pre-registration's own warning. The compression-error-dominates-gap caveat is reported as a deviation that BLOCKS the two-limit program under the uniform convention — not smoothed over. The c₀-masquerade criterion was applied at every extraction. q=7 L=3 is a reported wall, no extrapolation. Exact rationals provided. The compressed gap is explicitly NOT identified with the true rate sequence._
