# RESULT — EPS-LATTICE (Option 1): the stopping-time renewal residue ε_S carries NO robust base-3 log-periodicity; the apparent signal was artifacts. Honest NULL (2026-07-28)

**Probe:** `probes/probe_eps_lattice.py`. **Question (Option 1 of the ladder↔spectral bridge):** does the stopping-time renewal residue `ε_S(k)` (`k=log₂N`) carry the **same log2/log3 incommensurate log-periodicity** the LATTICE result found in the Plancherel tower `Λ_i`? If yes, the stopping-time arc and the S_∞ arc share one oscillation mechanism. **Verdict: NO robust signal — the apparent 4–8σ was three stacked artifacts; clean regime = 1.6σ. NULL (not detected, not excluded).**

## The prediction (sharp, falsifiable)
LATTICE (`result_LATTICE.md`): the tower `Λ_i` (indexed by level `i`, modulus `3^i`) has the base-3 lattice mode **aliased to a constant at integer `i`**, and the base-2 (`÷2^v`) oscillation **visible** at period ≈9 = `2π/log2`; `log3/log2=1.585` irrational ⟹ quasi-periodic. `ε_S` is the **exact mirror** (indexed by `k=log₂N`): a log-periodic-base-`b` signal in `log N` has period `log(b)/log2` in `k` — base-2 → period **1.000** (aliases to DC at integer `k`, which is why the old `2^{32,34,36}` data looked like non-monotone noise), base-3 → **1.585**. Sampling `k` at fine/irrational spacing resolves both.

## What was run
Dense per-octave Monte-Carlo of `ε_S(k) = ⟨σ_S − logN/log(4/3)⟩` over starts in band `[2^{k-1}, 2^k)` (per-octave sampling removes the huge cross-octave `σ_S` variance ⟹ SE~0.008 at 8–10M orbits/band vs 50M for the `[3,N]` convention). Numba walk kernel matches certified `experiments/38_eps_S_log4_test.py` (`σ_S` = # Syracuse steps to 1). Progressively stricter analyses.

## The finding — the signal evaporates as artifacts are removed
```
   analysis                                    base-3 (P=1.585)   what it included
   loose (all k, linear trend)                     7.8σ           trend-leakage + low-k hump
   k≥16, quadratic trend, 0.25 grid                3.9σ           peak was actually P=1.0 (binning artifact, 6.4σ)
   irrational k-grid, k≥16                          4.8σ           binning killed; low-k hump still in
   k≥20 (past hump), irrational grid, SE-inflated   1.6σ           CLEAN regime; peak misplaced at P=1.405
```
Three artifacts stacked to fake the signal: **(a) trend-leakage** — a long-period cosine mopping up the strong `1/N` transient (giant false peaks at P=5.7/9.1); **(b) the low-k hump** (k=16–20, where `ε_S` rises to ~1.52 then settles — the *tail of the finite-N convergence*, not a persistent oscillation); **(c) an octave-binning artifact at P=1.000** — bins tied to powers of 2 make `frac(k)` period-exactly-1.0, so any binning systematic aliases to "base-2"; killed by the irrational grid (`√5−2` spacing), which dropped the P=1.0 tone `6.4σ→3.2σ→0.5σ`. In the clean regime (past the transient, artifact-immune grid, SE inflated ×1.39 for the mild ×1.9 excess scatter), the base-3 tone is **1.6σ and the periodogram peak sits at 1.405, not 1.585.**

## Honest scope — "not detected," NOT "excluded"
Detection floor here is amp ≈ 0.01 (N≤2³², 10M orbits/band). The tower's oscillation is a ~10–20% modulation of `log Λ`; the *same mechanism* in `ε_S` could be ≲0.005 (a <0.4% wiggle on 1.38) — **below the floor**. So the clean read is: `ε_S`'s finite-N wiggle is consistent with an ordinary `1/N` transient + MC noise, and any shared log-periodicity with the tower is too small to see at reachable N. Pushing it (larger N, 4× orbits) only halves the floor for a signal there's no positive reason to expect — low expected value.

## Net
- **The arcs are NOT visibly bridged through `ε_S`'s oscillation.** The "unexplained finite-N oscillation" that `closed_form_findings` flagged (3 mechanisms ruled out) looks, on dense sampling, like plain convergence transient — not a log2/log3 fingerprint.
- Validation the probe is sound: `ε_S(3,2)` trend asymptote and the on-record ~1.37–1.38 band reproduce; the null is about the *oscillation*, not the value.
- **The bridge, if real, is in the transcendence STRUCTURE, not the finite-N data** — which is exactly what Option 2 (the two-walls / one-idele verdict, `result_TWO_WALLS.md`) then established.
- **Not at stake:** S_∞≈0.475 (floor 0.473177), MAHLER, PHYDRA_FAMILY, MIRROR, GARSIA, DENOM, SOLSTICE, R1–R30. This is a clean negative on a bridge hypothesis, not a change to any value or structure.
