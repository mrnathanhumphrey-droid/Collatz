# Probe ILEDGER2 — un-normalized interference ledger to j=12 — **the profile path validates (g_r exact to r=12, reproduces the Λ ladder), but the DIAGONAL KERNEL is under-specified: `|D|²=1/(5−4cosθ)` is the first-moment symbol (DC=1, mass-EXPANDING, ρ_prop≈1.17), NOT the mass-contracting v=v′ diagonal Σ_v P(v)²=1/3 (θ-independent constant). The fresh-source sign is kernel-dependent — both prior leans (0.477 and 7/15) are RETRACTED until the diagonal is pinned.**

**Date:** 2026-07-25. Probe `probes/probe_ILEDGER2.py`. Rebuilds the interference ledger with the intended
un-normalized (mass-contracting) diagonal and extends δ_j via `build_nu`→dlog→|FFT|² to j=12 (past the mu-ladder
wall at 8), for the late-window rate (G1–G4).

## G0 — profile path VALIDATED
`⟨δ_j, Re w⟩` (build_nu→dlog pushforward→|FFT|²) reproduces the known `g_j = (Λ_j−Λ^unif)/S_j` to <10⁻⁷ for
j=2..7. **The FFT profile path is correct, and the exact deviation coupling `g_r` is now available to r=12.**

## P4 / G2 — telescoping holds; g_r reproduces the exact Λ ladder
`Σ_j A[j,r−j] = g_r = ⟨δ_r, Re w⟩` for r=2..12 (telescopes for any adjoint pair). And `g_r/g_{r−1}` matches the
exact `Λ_r/Λ_{r−1}` to ~10⁻³ (r=4..12: 0.490/0.501/−1.03/−0.635/1.716/0.985/1.071/0.988/0.870). So the ledger's
overall object *is* the Λ chain — as it must be. **This is the useful, kernel-independent gain: g_r extended to
r=12 via a validated path.**

## ⚠️ The kernel is wrong — G4 is the tell
Row rate `ρ_prop = A[j,k+1]/A[j,k]` came out ≈ **+1.17 (>1, EXPANDING)**, not a sub-1 contraction — and
`Σ_{a} wD(children of DC) = 1 + 1/7 + 1/7 = 9/7`, not the mass-contracting **1/3** derived from `Σ_v P(v)²`. Root
cause is a genuine conflation:

- `wD = 1/(5−4cos2πθ) = |Σ_v P(v)e^{ivθ}|²` is the **first-moment** transport symbol: **DC value 1**, and it
  contains the **v≠v′ cross terms**.
- The **second-moment diagonal** (v=v′ part of `|A|²+|B|²`) is `Σ_v P(v)² = Σ_v 4^{−(v+1)} = 1/3` — **θ-independent
  (a constant), not `1/(5−4cos)`**. (This is R17-B/C's measured ⟨|D|²⟩≈1/3 — but that is the *average* of the
  symbol, i.e. the collision probability, not the pointwise kernel.)

So `wD(0)=1` and the diagonal-DC-weight `1/3` are **different objects**. The same-parity diagonal `|A|²+|B|²` is
`(constant 1/3 from v=v′) + (same-parity cross)`, and the opposite-parity source is **not** the full `1/(5−4cos)`.
**The |D|² symbol is not the diagonal kernel — the diagonal is under-specified.**

## Consequence — the fresh-source sign is kernel-dependent; no verdict
`⟨s_j, Re w⟩` (the fresh coupling) depends on which kernel defines `T̃_diag`:

| run | kernel | fresh-source sign, j=2… |
|---|---|---|
| ILEDGER (normalized avg) | stochastic (÷Z) | `− − − − − −` (all neg, j≤7) |
| ILEDGER2 (un-norm wD) | expanding (×wD) | `− + − − − + + + + + +` (pos for r≥7) |

Different kernels → different signs. **So neither lean survives:** the earlier "leans 0.477" (ILEDGER, wrong tail
factor + normalized kernel) and the corrected-rate "reverses to 7/15" (which assumed a mass-contracting kernel that
this build shows `wD` is not) are **both retracted.** The sign rides entirely on the diagonal kernel, which is not
yet pinned. G1 (source rate) and G3 (margin constant) are **not meaningful** under the wrong kernel and are not
reported as results.

## What is needed (pen)
The exact **same-parity diagonal branch weights** — the v=v′ constant `1/3` plus the same-parity (even-even,
odd-odd) cross terms — as a kernel on the Plancherel profile, distinct from the opposite-parity source
`2Re(A B̄)`. Only then is `T̃_diag` well-defined, `s_j` the true interference term, and the fresh-source sign
computable. The candidate that would be *mass-contracting* is the one whose 3-branch weights sum to ≤1 (the v=v′
diagonal sums to 1/3); `wD` summing to 9/7 disqualifies it.

## Status
**ILEDGER2:** (G0) profile path via build_nu→dlog→|FFT|² **validated** — `⟨δ_j,Re w⟩=g_j` to 10⁻⁷, and **g_r now
exact to r=12**, reproducing the Λ ladder (G2). (⚠️) The diagonal kernel is **mis-specified**: `wD=1/(5−4cosθ)` is
the first-moment symbol (DC=1, **mass-expanding**, ρ_prop≈1.17), while the mass-contracting v=v′ diagonal is the
**constant** `Σ_v P(v)²=1/3` — two different operators. **The fresh-source sign is kernel-dependent** (normalized:
all-neg; un-normalized: pos for r≥7), so **both the 0.477 lean and the 7/15 lean are RETRACTED** — the interference
ledger cannot deliver a sign until the same-parity diagonal is defined exactly (pen). Not at stake: R1–R30,
R80–R82; the g_r ladder to r=12 (validated). Still open and now correctly scoped: the diagonal-kernel definition,
which is the sole thing between here and a computable coupling sign.
