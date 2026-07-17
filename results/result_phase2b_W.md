# Result — PROBE W: the QUASI-STATIONARY-weighted compression FIXES Probe E's blocker. The partner resolves to rel err ~1e-5, correct side of c₀. The two-limit program closes with this convention.

**Date:** 2026-07-16. Reweighted compression (Probe E's named next lever). Recon; no proof, no rate fit. Direct/dense at q=3. Probe `probes/probe_phase2b_W.py`, log `logs/probe_phase2b_W_log.txt`.

**Headline: replacing E1's UNIFORM source average with a QUASI-STATIONARY weight — μ(src) = |M's dominant right eigenvector| (the surviving-mass measure) — resolves the compressed partner to rel err 9.8e-7 (L=2) / 1.4e-5 (L=3), ON THE CORRECT SIDE of c₀. That error is now BELOW the true gap (1e-4 at L=3), so the compressed gap tracks the true coalescence — exactly what the uniform convention could not do. c₀-deflation FAILS (wrong magnitude). The free choice from E1 is REPLACED: quasi-stationary, not uniform.**

## Result
`Lmat_μ[c,c'] = (Σ_{src∈c} μ(src)·Σ_{dst∈c'} M[dst,src]) / (Σ_{src∈c} μ(src))`. Compressed-partner = top non-family eigenvalue distinct from c₀ (c₀-masquerade criterion). Compared to the true partner (subspace-extracted).

| L | true partner (side, gap to c₀) | U uniform | **Q quasi-stationary** | D c₀-deflate |
|---|---|---|---|---|
| 2 | 0.346827 (above, 2.9e-3) | 0.341016 — **below** ✗, rel 1.7e-2 | **0.346827 — above ✓, rel 9.8e-7** | 0.133 ✗, rel 0.62 |
| 3 | 0.333236 (below, 1.0e-4) | 0.334312 — **above** ✗, rel 3.2e-3 | **0.333231 — below ✓, rel 1.4e-5** | 0.272 ✗, rel 0.18 |

- **Q (quasi-stationary) WINS decisively.** μ = |dominant right eigenvector of M| (the quasi-stationary / surviving-mass distribution — the standard reweighting for a non-lumpable chain). Partner captured to **9.8e-7 (L=2), 1.4e-5 (L=3)**, both on the **correct side of c₀**. At L=3 the residual (1.4e-5) is **an order of magnitude below the true gap (1.0e-4)** ⇒ the compressed gap under Q now RESOLVES the true coalescence (the compressed partner sits at 0.333231, gap-to-c₀ ≈ 1.05e-4 ≈ the true 1e-4).
- **U (uniform, E's baseline) puts the partner on the WRONG side** (below c₀ at L=2, above at L=3) with error 1.7e-2 / 3.2e-3 ≫ the true gap — reproducing E's blocker.
- **D (c₀-deflation via Real-T1 ℓ₀) FAILS** (rel 0.62 / 0.18) — deflating the known c₀ mode does not isolate the partner cleanly here (the near-degenerate right eigenvector r₀ is EP-contaminated; the rank-1 subtraction distorts). Not the fix.

## What this resolves and what it opens
- **Resolves E's caveat:** the two-limit program (compressed Perron→1/3 AND compression error→0) DOES close — but with the **quasi-stationary-weighted** compression, not uniform. Under Q the compression error (~1e-5) vanishes faster than the true gap, so the compressed chain faithfully tracks the partner and its coalescence with c₀.
- **Reframes the Lambda question (open, next instrument step):** applying Q at L=4 needs μ = the dominant right eigenvector of the full 236k-state operator. That is the **LARGEST-modulus** mode, obtainable by **power iteration = pure sparse matvec** — NOT the banned interior-near-EP solve, and GPU-friendly. So a Q-weighted L=4 compressed chain is reachable via (i) power-iteration μ on the 236k operator (Lambda GPU, or big-RAM CPU), (ii) direct accumulation of Lmat_μ (dim 4374). This is a *cleaner* Lambda use than the walled interior LU — but only IF the L=4 point is wanted (rate-fitting stays banned; it would be a confirmatory cross-check, not load-bearing).

## Adjudication
| weight | verdict |
|---|---|
| U uniform (E1's frozen choice) | wrong side of c₀, error ≫ true gap — E's blocker. |
| **Q quasi-stationary (μ=\|dominant right eigvec\|)** | **partner to rel ~1e-5, correct side; closes the two-limit program.** New frozen convention. |
| D c₀-deflation | fails (EP-contaminated r₀ / rank-1 distortion). |

**⟹ The E1 convention is amended: use the quasi-stationary source weight.** The hand-derivation of the entry formula (E3) is unchanged (the entries are the same gate algebra); only the class inner product / averaging weight changes from uniform to quasi-stationary — and that is what makes the compressed Perron/partner faithful.

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, H, J, L-A, L-B, Real-T1, R1-STOP, Probe P/C/E (E's uniform-convention finding stands — it's why we reweighted). No `r_q` value changes; **no rate-law fit** (the true 2.9e-3, 1.0e-4 untouched; the Q-compressed gap now agrees with it but is NOT fitted to it).

_Reporting discipline: Q's success is reported with the correct-side check AND the rel-err vs true-gap comparison (1.4e-5 ≪ 1e-4), not just a distance. U reproduces E's blocker (consistency check). D is reported as a failed candidate, not omitted. The c₀-masquerade criterion was applied at every extraction. The L=4 route via power-iteration μ is flagged as OPEN (not run), with the honest note that it is confirmatory (rate-fitting banned)._
