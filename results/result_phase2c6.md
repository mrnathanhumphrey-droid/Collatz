# Probe 2c6 — THE JOINT (g₀, g₁) CORRECTOR SEARCH (reveal)

**Date:** 2026-07-17  β*=3/5 frozen; corrector `h = h_β·(1+g₀[τ] on v₀)·(1+g₁[cls] on v≥1)`, both trits mean-zero.
Wilson's blind joint-optimum ruling was on the record first (decoupled = reference/direction; the two
mod-27 couplings named as the deviation's owners). This is the search.

## Q1 — the real joint shrink (modest; the trits are nearly spent)
| | L=2 | L=3 |
|---|---|---|
| baseline (rung-1, 9/49) | 0.183673 | 0.183673 |
| v₀ trit ALONE (2c4) | 0.160488 (1.144×) | 0.160585 (1.144×) |
| **JOINT (g₀,g₁)** | **0.151896 (1.209×)** | **0.160353 (1.145×)** |

At **L=3 the joint optimum sets g₁ → 0** (verified a true optimum: the best g₁ at the g₀-optimum buys only
0.0002, 0.15%). The class-3 v₁ dressing adds essentially **nothing** at L=3 over v₀-alone. (At L=2 g₁ helps —
but L=2 is degenerate: the **D9 class is empty**, so the v₁ trit is incomplete there.)

## Q2 — which structure carries the residual  →  the coupling, off the trit ladder
- Residual carrier (widest cell): **O,v1** at L=3 (odd, v≥1), width 0.1604 = the full bracket.
- The post-joint residual is **NOT a function of the trit coordinates at any depth**: on (a mod9, γ mod M, e mod6)
  the bad-key count is **identical (936) for M = 27, 81, 243** — deeper γ-resolution does not resolve it.
  So the binding spread lives **within** the trit classes, injected by the **v₀↔v₁ coupling** (dressing v₀
  perturbs the v₁ ratios by an amount that is *not* a trit-coordinate function). This is exactly the coupling
  Wilson flagged blind as the deviation's owner — the search confirms the trit dressings are spent and the
  couplings own what's left.

## Q3 — does anything finally MOVE with L?   (the one that outranks the rest)
- **Rungs 1 & 2 are L-locked**: baseline width L=2 = L=3 to 2e-16; v₀-alone drifts only 1e-4 (generic c₀).
- **The JOINT width is the FIRST corrector quantity that is NOT L-invariant**: 0.1519 (L2) → 0.1604 (L3),
  Δ = 8.5e-3 — ~100× the v₀-alone drift, and it lives in the **coupling residual**, precisely the address the
  shell picture fingered (shell meets fold at the cap rung; everything below-cap is L-blind, as observed).
- **HONEST CAVEAT (does not overclaim):** the L=2↔L=3 comparison is **confounded** — D9 is empty at L=2, so the
  v₁ trit + its coupling only *fully exist* at L=3, and the motion is *upward*, not a clean partner-descent to
  1/3. We have exactly **one** level (L=3) with the complete structure. The **decisive** L-flow test is **L=4**
  (the first pair of levels both carrying the complete v₁ trit + coupling) — which is the **cap-rung** the shell
  picture predicts should carry the rate law. **L=4 tower is heavy (L=3 was 8424 states / 892k nnz; L=4 ~10×) →
  needs a greenlight / Lambda, not fired here.**

## Verdict
The joint search independently **lands on Wilson's blind ruling**: the two trits (v₀ τ, v₁ b/η) are each
closed-form and L-locked but **nearly spent** together (joint shrink only 1.21×/1.145×); the residual is owned
by the v₀↔v₁ **coupling**, which is **off the trit ladder** (unresolved by any γ-depth); and the coupling
residual is the **sole L-carrier** in the whole chain — the only thing that moved with L. The shell picture's
central prediction (rate law lives at exactly one address, the cap rung, and nowhere else) is **consistent with
every measurement here**, but the actual **breathing with L is not yet confirmed** — that is the L=4 cap-rung test.

Probe `probes/probe_phase2c6.py`; log `logs/probe_phase2c6_log.txt`.
