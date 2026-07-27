# RESULT — P6G: the kernel is three points and the formula is EXACT — Λ_i = telescoped [4R_e(2)−R_e(0)], no boundary (2026-07-26)

**Probe:** `probe_p6g.py` (+ inline high-precision checks). Wilson assembled the kernel K to a closed form in base-2 lag
Fourier and gave a falsifiable telescoped formula for Λ_i. This gates it. The kernel structure and positivity are
exactly right; one factor-4 normalization is corrected; the boundary is **identically zero**.

## Wilson's kernel algebra — EXACT (structure), off by a clean factor 4 (normalization)
`K̂ = ¼(25−16c²)·[(15/2)/(25−16c²) − ½] = 2c² − 5/4 = cos2θ − ¼` verified to 3e-16 on a θ-grid. But the certified
per-level pairing is **4×** that:
```
P̃_i := Σ_{k≥1} 4⁻ᵏ ⟨ρ_i, shiftₖ ρ_i⟩   (certified full base-4 autocorrelation pairing)
     = 4 R_e(2) − R_e(0)     EXACTLY (boundary 1e-18, ratio 4.0000 vs R_e(2)−¼R_e(0), all i=2..6)
```
So the true kernel is **K = {n=±2: 2, n=0: −1}, K̂ = 4cos2θ − 1** — Wilson's three-point support (±2 and 0) is dead
right, the slip is a single factor-4 in the composition (the `Re w` channel-weight normalization / a base-2↔base-4
conversion). The positivity condition is **unchanged** by the factor: `P̃_i > 0 ⟺ R_e(2) > ¼R_e(0)` — Chebyshev/covariance
family, the m=0 proof's family. Holds at every level i=1..6 (all `R_e(2) > ¼R_e(0)`).

## THE FORMULA — machine-exact, NO boundary
```
Λ_i = 3ⁱ[4 R_e⁽ⁱ⁾(2) − R_e⁽ⁱ⁾(0)] − 3ⁱ⁻¹[4 R_e⁽ⁱ⁻¹⁾(2) − R_e⁽ⁱ⁻¹⁾(0)]
```
where `R_e⁽ⁱ⁾` = autocorrelation of the even-branch sub-measure ν_e at level i, base-2 lags 0 and 2 (ratio-4
autocorrelation and self-collision). Gated vs certified Λ_i (shellA):

| i | Λ_i certified | recon (R_e only) | diff |
|---|---|---|---|
| 2 | −0.00730790 | −0.00730790 | −4e-16 |
| 3 | +0.00131986 | +0.00131986 | −1e-16 |
| 4 | +0.00065026 | +0.00065026 | +4e-16 |
| 5 | +0.00032692 | +0.00032692 | +6e-16 |
| 6 | −0.00033867 | −0.00033867 | −4e-16 |

Anchors: Λ_1 = −2/21 exact (convention pin); telescoping `3ⁱP̃_i − 3ⁱ⁻¹P̃_{i−1} = Λ_i` ratio 1.0000; 3ⁱ normalization
validated (`3ⁱ⟨ρ,shift₁⟩ → 0.733` = M₋ + cascade). **The boundary term is 1e-18 — identically zero, cleaner than the
derivation predicted.** Every arrow verified; the constant stands on two autocorrelation values per level.

## Where it leaves 7/15
`Σ_{i=2..6}Λ_i = −0.00535` (partial). Target `7/15 ⟺ Σ_{i≥2}Λ_i = −1/210 = −0.00476`. The partial sum is in range but
the tail (i>6) and its convergence are the open question — the per-level *formula* is closed and exact; whether the
telescoped `[4R_e(2)−R_e(0)]` difference sums to exactly `−1/210` is now a statement purely about the level-sequence of
`R_e(2), R_e(0)` (one sub-measure, two lags per level, no free maps, no boundary). **Pen next (Wilson):** reconcile the
factor 4 in the composition (structure/positivity unaffected), and the asymptotics of `3ⁱ[4R_e⁽ⁱ⁾(2)−R_e⁽ⁱ⁾(0)]`.
Not at stake: P6D/P6E/P6F chain, P6/P6B/P6C, P1LVL, BRIDGE2, CHANNEL_ID, dichotomy, R1–R30.
