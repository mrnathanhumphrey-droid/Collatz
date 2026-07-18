# Probe 2c6-L4 — the cap-rung L-flow test (Lambda A10, GPU SpMV)

**Date:** 2026-07-18  L=4 joint (g₀,g₁) corrector search, β*=3/5 frozen. Ran on Lambda gpu_1x_a10 (ML-CF lane):
tower 233,280 states, full op 229.6M nnz, build 247s (CPU), search 110s (CuPy SpMV), instance terminated. ≈$0.65.

## The L-sequence (the answer to Q3 — does the joint bracket breathe with L?)
| | L=2 | L=3 | L=4 |
|---|---|---|---|
| baseline (rung-1 9/49) | 0.183673 | 0.183673 | 0.183673 |
| v₀ alone | 0.160488 | 0.160585 | 0.160622 |
| **JOINT (g₀,g₁)** | **0.151896** | **0.160353** | **0.160365** |
| joint g₁ | (nonzero) | → 0 | → 0 |
| joint bracket | — | [0.269048, 0.429401] | [0.269040, 0.429405] |

## Verdict: the joint bracket is **L-INVARIANT once the structure is complete** — it does NOT breathe.
- The L=2→L=3 apparent motion (Δ=8.5e-3) was **entirely the D9-empty degeneracy at L=2** (the v₁ trit isn't fully present there). With the complete structure at both L=3 and L=4, the **joint width is flat to ~1e-5** (0.160353 → 0.160365) and the bracket **endpoints are stable to ~1e-5**. g₁→0 reproduces at L=4.
- The residual is off the trit ladder at L=4 too, and **more sharply**: bad-key count *grows* with γ-resolution (a9·γ9·e6 → 324; γ27 → 972; γ81 → 2880) — finer γ never resolves it; the coupling splits every class.

## What this means (refines the shell picture, honestly)
- **Below-cap rungs are L-locked — confirmed at three levels.** baseline and the trit dressings give the same width at L=2,3,4. Exactly the shell-picture prediction.
- **But the joint bracket plateaus at an L-invariant floor (~0.1604), set by the off-ladder v₀↔v₁ coupling residual.** The partner's real L-flow — the distinctive 2.9e-3 → 1.0e-4 descent of ρ(M_tower) to 1/3 — is **3+ orders of magnitude below this floor** and is **invisible to the trit-resolution C-W bracket.** The tiny width drifts that do exist (v₀-alone 0.16049→0.16059→0.16062) are generic c₀-type, decelerating, and go the **wrong way** (up) for a convergence-to-1/3.
- So the honest reading of "does anything finally move with L": **not in this instrument.** The cap-rung rate law is not disproved — it lives *below* the coupling floor, which this bracket cannot resolve. **The lever is not another rung on the ladder; it is killing the off-ladder coupling residual (~0.16) so the partner's 1e-4 flow becomes visible.**

## Status
Two-trit corrector chain characterized to L=4: closed-form trits, L-locked, spent; residual owned by the off-ladder v₀↔v₁ coupling; bracket floor L-invariant. The contraction/limit theorem does **not** close by adding ladder rungs — the next object is the **coupling itself** (the sub-trit, non-γ-graded mixing that dressing v₀ injects into v₁). That is where ρ(M_tower)→1/3 must ultimately be read, and where Wilson's pen goes next.

Probe `probes/l4_joint_gpu.py` (GPU); mirrors `probes/probe_phase2c6.py` (L=2,3).
