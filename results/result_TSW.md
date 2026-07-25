# Probe TSW — transport-adjoint direction convergence — **the slow mode is REAL (no complex pair), which CLOSES the "oscillate-forever → 0.477" route but does NOT by itself give 7/15: reality ⟹ finitely many sign changes, and the crux is now the SIGN of the slowest contributing mode's coupling — a real-mode question, not a complex one. The observed g_r (single-signed + since r=7, decaying at exactly the slow-mode rate) currently leans 0.477.**

**Date:** 2026-07-25. Probe `probes/probe_TSW.py`. **Exact-operator analysis** (not the S2 splitting estimate) —
its output may touch the exact ladder. Certified operators only: freq-domain identity from R14/R11 (Re w, profile,
g_r); the linear transport operator is the **R28/R29 gap matrix** `M` (the one certified operator advancing the
renewal r−1→r linearly on the fixed lattice {R(2m)}). No fresh construction.

## TSW-A — ground-truth gate: PASS
`g_r = ⟨δ_r, Re w⟩ = Λ_r/S_r` (up to `Λ^unif/S_r`, doubly-exp small) confirmed exact r=1..7 (all OK; `Re w(x) =
15/(2(17−8cos2πx)) − ½` verbatim R14). Sign sequence:

`sign(g_r), r=1..16:  − − + + + − + + + + + + + + + +`   (matches the pre-registered `−−+++−+…` on r=1..15).

**Last sign change at r=7** (Λ_6<0, Λ_7>0); single-signed **positive** through r=16 (and, per S2, through r≈35).

## TSW-B/C — the decider: the slow mode is REAL (robust across truncation)
Power-iterate `T* = M^T` on the mean-zero complement (leading/uniform mode Hotelling-deflated; `|proj on leading|`
held ≤10⁻¹⁵ every step — mean-zero guardrail satisfied). Seed = the certified gap-domain `W_gap[m]=4^{−|m|}`.

| D | direction overlap ⟨Ŵ_r,Ŵ_{r−1}⟩ (it 2→14) | Rayleigh → | #complex eigenvalues | slow eig |
|---|---|---|---|---|
| 8 | 0.971 → **0.99991** | +0.975 | **0** | +0.980 REAL |
| 10 | 0.968 → **0.99945** | +0.987 | **0** | +1.005 REAL |
| 12 | 0.968 → **0.99944** | +0.979 | **0** | +1.001 REAL |

**Direction converges monotonically to 1 (no rotation); the full spectrum has ZERO complex eigenvalues at every
reachable D.** The slow eigenvalue is **real ≈0.98–1.00** (leading ≈1.083). This is the robust real/complex read
(direction converges far faster than value), and it agrees with R29-A's direct diagonalization (all-real top-6).

**Consistency:** the observed local decay `Λ₁₆/Λ₁₅ = 0.894` equals the relative rate `λ_slow/λ_leading =
0.98/1.083 = 0.905` — so the observed g_r tail (positive, decaying) **is** the slow real mode, already dominant by
r≈12–16.

## What REAL actually means here — and why it does NOT hand us 7/15
The naive reading "real ⟹ 7/15" is **too quick**. Correctly:

- **Real spectrum ⟹ finitely many sign changes ⟹ g_r eventually single-signed** — Wilson's premise, confirmed.
- **The "g_r oscillates with lengthening period, never converges → 0.477" route is now CLOSED** — there is no
  complex eigenvalue to produce a genuine oscillation. The crossings (r≈2.65, 9.01) and any r≈36 rollover are
  **real-mode competition** (a finite number of sign changes as one real mode overtakes another), **not a rotating
  complex mode.** This directly corroborates the **lengthening-transient seam** (a finite transient, asymptotically
  frozen) and **refutes** the complex-pair / genuine-log-periodic reading (and, with the S2 waveform retraction,
  the ×3 story is doubly dead).
- **But single-signed does not mean single-signed *negative*.** 7/15 (ε_∞=0) needs g_r eventually **negative** (to
  pull ε back down from its current +9×10⁻³). A real **positive** slow mode with **positive** coupling gives g_r
  single-signed **positive** → ε rises to a positive limit → **S_∞≈0.477**. The decider is therefore the **SIGN of
  the slowest *contributing* mode's coupling**, not real-vs-complex.

## TSW-D — the sign is the crux, and it is NOT decided here (Wilson's derivation)
The coupling `⟨δ_0, φ_slow⟩` cannot be read cleanly from this probe: the eigenvector's sign convention is
arbitrary, so the data-projection sign is not meaningful without the correct δ_0 pairing. What the **data** says:
g_r has been single-signed **positive since r=7**, decaying at exactly the slow-mode rate, with **no complex mode**
available to flip it. On its face this points to **g_r stays positive → S_∞≈0.477 (NOT 7/15)** — the "eventually
single-signed" appears to have already happened at r=7.

For **7/15** to survive, there must be a **slower** contributing real mode (|λ| closer to the leading ≈1.083 than
the 0.98 currently seen) carrying a **negative** coupling, which overtakes the +0.98 mode and forces **one more
(+ → −) crossing** — the r≈36 rollover. Whether such a mode exists and carries a negative coupling is the exact
**positivity/sign question** — scalar, on the certified operator, and **Wilson's pen.** TSW has reduced 7/15 to
exactly this one inequality and removed the complex-oscillation escape hatch.

## Status
**TSW (exact-operator analysis):** (A) gate PASS — `g_r=⟨δ_r,W⟩=Λ_r/S_r`, signs `−−+++−+` then all + (last change
r=7). (B/C) **the slow mode is REAL** — direction overlap →0.9999, **zero complex eigenvalues** at D=8,10,12,
Rayleigh→0.98, and `λ_slow/λ_lead=0.905` matches the observed `Λ₁₆/Λ₁₅=0.894`. **This CLOSES the complex-
oscillation route to 0.477** and confirms **finitely many sign changes** (lengthening-transient seam corroborated,
complex/×3 reading refuted). **But reality ⟹ 7/15 is false as stated:** with the observed g_r single-signed
**positive** since r=7 and no complex mode to flip it, TSW currently **leans S_∞≈0.477** — the same direction as
the S2 deep run, now from an independent *exact-operator* argument. **7/15 survives only if a slower real mode with
negative coupling forces one final (+→−) crossing (the r≈36 rollover)** — the crux inequality, handed to the pen.
**Caveat (R29-B):** |λ_slow| value is not D-converged (0.98→1.00 across D=8..12); its *reality* is robust at
reachable D, but if λ_slow is a continuous-spectrum edge the value is not a discrete eigenvalue (the reality still
holds numerically). Not at stake: R1–R30, R80–R82 (exact identities, M-reality). Under pressure from two
independent methods (S2 splitting + TSW exact operator): the value **7/15**; the surviving path is a single
sign/positivity inequality on the slowest real mode.
