# Probe FIBER — the fiber-mean / fiber-fluctuation split (replaces the retracted interference ledger) — **GATE 1 ✓ (M Re w = Re w with 4→64, machine-clean); Λ^unif closed form ✓ (its doubly-exp death IS the 4→4^{3^k} ladder). ⚠️ The GATE-3 sign claim below is CORRECTED by `result_modes.md`: the fiber-mean is NOT always subdominant (it DOMINATES at r=6, driving the g₆ crossing), so "no crossing in 10 levels / single channel / 0.477" was WRONG — g_r crosses at r=2 and r=6. See the correction note.**

> **⚠️ CORRECTION (2026-07-25, `result_modes.md`).** The GATE-3 table below started at r=7 and **omitted r=6**, where
> `⟨dbar,Re w⟩=−1.17e−3` DOMINATES `⟨dfluct,Re w⟩=+4.4e−4`, so **g₆ is negative — a real crossing** (matching the
> corrected TSW-A `0 − + + + − +`). The fiber-mean is subdominant only for r≥7; at r=6 it overpowers the fluctuation
> precisely because the fluctuation is near its minimum there. So "fluctuation positive ⟹ no crossing ⟹ 0.477" is
> retracted: the fluctuation coupling IS positive, but g_r crosses. The mode decomposition (`result_modes.md`) shows
> the real structure — the dominant mode Re δ̂(1) is single-signed positive, and the crossings are n=2 / n=3
> excursions (the lengthening-period mechanism). Also: the fiber ladder and the S2 splitting are **one method
> extended** (S2 continues this exact g_r ladder past r=16), not two independent methods — the honest gain here is
> that the exact side now reaches r=12.

**Date:** 2026-07-25. Probe `probes/probe_fiber.py`. The interference-ledger split (ILEDGER/2) is retracted: in the
dlog domain the branch `(−2)^{−v}` acts by translation (⟨−2⟩=⟨4⟩=G_r cyclic odd order 3^r; β(−2)=(3^r+1)/2
verified r=2,3,4), so `T̃_diag = (1/3)·Id` and the ledger telescoped to `g_r` with zero content. The real 3-fold
structure is the **level-lift fiber**, not the branch. δ_r via the validated build_nu→dlog→|FFT|² path (r=2..12).

## GATE 1 — fiber map correct (machine-clean)
Primitive `χ_k` (level r) → `ψ_m` (level r−1), `m=k mod 3^{r−1}`, 3-to-1, fiber `{m,m+3^{r−1},m+2·3^{r−1}}`.
Fiber-average of Re w: `M(Re w)(m) = ⅓Σ_a Re w((m+a·3^{r−1})/3^r)`. **`M(Re w) = Re w` with `4→4³=64` to
7.8×10⁻¹⁷** (r=2..7) — the 3-average kills modes not divisible by 3 and reindexes `4^{−|k|}→4^{−3|k|}`.

## GATE 2 — the Λ^unif closed form is validated (a bonus result); the channel-identification was imprecise
**Validated:** the **primitive mean of Re w** over order-3^r characters is
`(1/(2·3^{r−1}))[3^r/(4^{3^r}−1) − 3^{r−1}/(4^{3^{r−1}}−1)]`, numerically **−1/7 (r=1), −0.00793079 (r=2),
−1.907×10⁻⁶ (r=3)** — matching the banked `Λ_r^unif/S_r` (the R10-C −1/7 trace; the byte-anchor −0.00793078).
**So `Λ^unif`'s long-unexplained doubly-exponential death is the `4→4^{3^k}` ladder** — a previously-banked
observation now given a closed form.

**Correction to the gate as posed:** `⟨fiber-mean of δ_r, Re w⟩ ≠ Λ^unif/S_r`. Since δ_r is already mean-zero
(the uniform `1/M` subtracted), its fiber-mean is level-(r−1) *deviation* content, **O(10⁻³)** (r=3:
`⟨dbar,Re w⟩=−3.2×10⁻³`), not the doubly-exp `−1.9×10⁻⁶`. The `Λ^unif` closed form lives in the *uniform*
pairing, which was already removed from g_r. The closed form is right; it is one layer over from the fiber-mean of δ.

## GATE 3 — the sign collapses to one dominant channel, and it is POSITIVE
`g_r = ⟨fiber-mean, Re w⟩ + ⟨fiber-fluctuation, Re w⟩`:

| r | ⟨dbar,Re w⟩ | ⟨dfluct,Re w⟩ | g_r | \|dbar\|/\|fluct\| |
|---|---|---|---|---|
| 7 | −1.10e−4 | **+5.71e−4** | +4.62e−4 | 0.19 |
| 9 | −3.71e−5 | **+8.17e−4** | +7.80e−4 | 0.045 |
| 10 | −2.17e−5 | **+8.57e−4** | +8.36e−4 | 0.025 |
| 12 | −1.34e−5 | **+7.32e−4** | +7.18e−4 | 0.018 |

- **Fiber-mean is subdominant and asymptotically dead:** it decays ~0.5/level (faster than the fluctuation), so
  `|dbar|/|fluct|` falls 0.19 → **0.018** by r=12. Your structural claim holds (the "doubly-exp" reason was
  overstated — the real rate is ~0.5 — but it is still negligible against the fluctuation).
- **The fluctuation coupling `⟨dfluct, Re w⟩` is POSITIVE for every r=3..12**, and its late-window rate is **~0.9**
  (r=10→12: 1.049, 0.972, 0.878) — the ladder rate, exactly what a sole surviving channel must show.

**Therefore `sign(g_r) → sign⟨fluctuation, Re w⟩ = +`.** For the first time in the arc the sign is **one clean
quantity**, not a competition: a single positive channel, decaying at 0.9<1, with no competitor and no crossing in
10 consecutive levels. **This points to S_∞ ≈ 0.477, not 7/15.**

## Status
**FIBER:** the interference ledger is retracted (branch = translation ⟹ `T̃_diag=(1/3)Id`, telescoping null); the
real structure is the level-lift fiber. **GATE 1 ✓** (fiber map machine-clean, `M Re w = Re w|_{4→64}`). **Λ^unif
closed form ✓** — its doubly-exp death is the `4→4^{3^k}` ladder (bonus result; the "fiber-mean=Λ^unif" wording is
corrected — the fiber-mean of the mean-zero δ is O(10⁻³), the closed form is the uniform pairing one layer over).
**GATE 3:** sign(g_r) reduces to the **fiber-fluctuation** coupling (fiber-mean subdominant, ~0.5/level, 1.8% at
r=12), which is **POSITIVE for all r=3..12 at rate ~0.9** ⟹ **S_∞≈0.477**, a single clean channel with no crossing
in 10 levels. Caveat unchanged: r≤12 extrapolation; a positive channel at rate 0.9 stays positive, but r>12 is
unobserved — a very-late crossing has zero evidence and now no structural mechanism (no competing channel) to
produce it. Consistent with the S2 deep-splitting lean (0.477). Not at stake: R1–R30, R80–R82, the g_r ladder to
r=12, and the Λ^unif closed form. **The 7/15-vs-0.477 decider is now a single, computed, single-signed channel —
and it reads 0.477.** Pen: whether the fluctuation coupling's positivity is provable (a real single-signed mode) or
merely observed to r=12 is the last open question.
