# W4 falsification — final disposition

**Date:** 2026-05-14 late evening
**Supersedes:** W4_DISPOSITION.md "PARTIAL IDENTIFICATION" verdict
**Status:** FALSIFIED by M_4 + ‖d_k‖² diagnostic (auditor deadlock-breakers #1+#2)

---

## Summary

The W4 agent's claim that Faure 2009's √3 ≈ 1.732 corresponds to `1/r_s(T_dev)`, with T_dev the bilinear deviation propagator on `{M_n(η) : η ≠ 1}`, **does not survive** the M_4 + ‖d_k‖² diagnostic recommended by the W4 adversarial audit.

## Evidence

**Computed (level k=1..4) from existing `bilinear_pair_operator.py` infrastructure:**

| k | S_k | ‖d_k‖² | ‖d_k‖ | n_modes | normalized (‖d_k‖/√n_modes) |
|---|---|---|---|---|---|
| 1 | 0.66667 | 0.11111 | 0.33333 | 1 | 0.33333 |
| 2 | 0.47619 | 0.16780 | 0.40963 | 5 | 0.18319 |
| 3 | 0.46157 | 0.38310 | 0.61895 | 17 | 0.15012 |
| 4 | 0.46421 | 0.45432 | 0.67403 | 53 | 0.09259 |

**Level-to-level ratios:**

| k→k+1 | raw ‖d_{k+1}‖/‖d_k‖ | normalized ratio | n_{k+1}/n_k |
|---|---|---|---|
| 1→2 | 1.229 | **0.5496** ← W4 anchor | 5.00 |
| 2→3 | 1.511 | **0.8194** ← climbing | 3.40 |
| 3→4 | 1.089 | **0.6168** ← oscillating | 3.12 |

**Targets:**
- 1/√3 = 0.5774 (W4's claimed asymptotic limit)
- √3 = 1.7321 (pure counting-tautology raw rate)

## Falsification logic

W4's claim required the **normalized ratio** to converge to 1/√3 ≈ 0.577 as k → ∞.

- Level 1→2 gave 0.550 (close to 1/√3, the W4 agent's evidence)
- Level 2→3 gave 0.819 (climbing AWAY from 1/√3)
- Level 3→4 gave 0.617 (oscillating back down, not toward 1/√3)

**The normalized ratio is oscillating, not converging.** The W4 agent observed a single data point that happened to land near 1/√3 by coincidence and extrapolated without computing the next level.

## What also is NOT supported

- **Pure counting tautology** (W4 audit's alternative hypothesis): raw ratios 1.229, 1.511, 1.089 average ~1.28, NOT √3 ≈ 1.732. So this is not "‖d_k‖ ~ √n_modes" either.
- **Bounded spectral decay** (the "W4 lives" reading): raw ‖d_k‖ grows from 0.333 to 0.674 (~2×) over four levels, not bounded.

## What IS supported

- **S_k → 7/15 = 0.4667** is reproduced cleanly (0.667, 0.476, 0.462, 0.464). This confirms W1's leading-order 7/45 derivation independently. The diagnostic positively validates W1's headline.
- **‖d_k‖² growth is sub-linear in n_modes.** Per-mode squared deviation ‖d_k‖²/n_modes = 0.111, 0.034, 0.023, 0.0086 — decaying, but not at any obvious rate.

## Implication for Faure √3

**Faure √3 has no demonstrated identification with a Syracuse cumulant or bilinear operator scale.** W4 is removed from Track A's "delivered" list.

Open questions left unresolved by Track A:
- What IS Faure √3, if anything, in the Syracuse picture? (Open)
- Is √3 a numerical coincidence with no Syracuse meaning, or does it correspond to a more subtle operator scale not captured by T_dev / κ_k^B? (Open)

The PADE asymptotic at z ≈ 1.016 (slow-mode rate ρ ≈ 0.984) remains the project-internal asymptotic prediction; √3 was an intermediate spectral candidate that did not survive.

## Updated Wrinkle inventory (post-Track A)

| Wrinkle | Pre-Track-A status | Post-Track-A status |
|---|---|---|
| W1 (B-amalgamated lift of HS) | conjectural | **CLOSED + AUDITED** (HS 2014 Thm 3.4) conditional on H1' |
| W2 (−1/30 closed form) | open | **PARTIAL CLOSED + AUDITED**: 1/30 = 1/(2·15) rigorous; rate (1/2)^n redirects to T_M λ_2 |
| W3 (PADE complex pair period 9.2) | borderline-pivot | Unchanged — Track A did not address |
| W4 (Faure √3 ↔ cumulant operator) | open | **FALSIFIED** — no identification with κ_k^B or T_dev |

## Files

- `C:/Collatz/w4_diagnostic_M4.py` — diagnostic script (computes M_k for k=1..4)
- `C:/Collatz/experiments_output/w4_diagnostic_M4.json` — full numerical outputs
- `C:/Collatz/W4_ADVERSARIAL_AUDIT.md` — audit that surfaced the falsification path
- `C:/Collatz/W4_DISPOSITION.md` (now superseded) — original "PARTIAL IDENTIFICATION" claim
