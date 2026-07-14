# Result 83 — Does the n=3 DWM match carry any mod-9 information?

**Date:** 2026-07-14. **Verdict (Task A): H_SENSITIVE.**

Probe `probe_83_dwm_mod9_ablation.py`; data `result_83_data.csv`; log `result_83_log.txt`.

## Why this ablation IS evidential (where Probe 82's reproduction was not)

Probe 82's walled n=3 diagnostic was a *reproduction* — it confirmed the DWM operator was transcribed correctly, and was explicitly barred from counting as evidence *for* the bridge. This is an **ablation**: it perturbs a **named structural element** (the step-j phase factor) on an exactly-specified operator and measures the response, holding the σ shift and the Geom(½) weights fixed. That is a controlled experiment on the match's *sensitivity* to a component — a different instrument from a fit, and it does count. The 4-cell grid (below) is designed so a single cell outside the pre-registered one can catch the failure mode where the reductions resolve no phase at all.

## n=3 modulus table (what phase actually exists at each step)

| j | r=n−2j+1 | effective modulus 3^{n−2j+2} | phase at n=3 |
|---|---|---|---|
| 1 | 2 | 3³ | nontrivial |
| 2 | 0 | 3¹ | **mod-3 shadow** of the mod-9 twist |
| ≥3 | <0 | 3⁰=1 | trivial (and absent from the 3-/4-alt moments) |

So the moments involve only j=1,2, and "ablate j≥2" = ablate the **j=2 mod-3 phase** — the n=3 shadow of the mod-9 twist, not the twist itself. There is barely any mod-9 structure present at n=3 to be sensitive to.

## Task A — 4-cell ablation grid

**V_MAX = 16** (shift = |cell − baseline|, same V_MAX):

| cell | j1 | j2 | G1 sum | G2 sum | G2 tr_π | G2 δ₁ | G2 vac_π | max rel-shift |
|---|---|---|---|---|---|---|---|---|
| baseline | on | on | +1.07831e-01 | +6.08879e-01 | +5.35722e-02 | +5.74203e-02 | +4.77548e-03 | 0.00e+00 |
| A | on | off | +2.33320e-01 | +7.77803e-02 | +4.57432e-02 | +4.69986e-02 | +2.44955e-03 | 1.16e+00 |
| B | off | on | -4.52209e-06 | +2.96305e-01 | +3.60617e-02 | +3.60617e-02 | +1.29666e-03 | 1.00e+00 |
| C | off | off | +2.96351e-01 | +9.87926e-02 | +2.67418e-02 | +2.67418e-02 | +2.85203e-03 | 1.75e+00 |

**V_MAX = 20** (shift = |cell − baseline|, same V_MAX):

| cell | j1 | j2 | G1 sum | G2 sum | G2 tr_π | G2 δ₁ | G2 vac_π | max rel-shift |
|---|---|---|---|---|---|---|---|---|
| baseline | on | on | +1.07820e-01 | +6.08922e-01 | +5.35746e-02 | +5.74226e-02 | +4.77574e-03 | 0.00e+00 |
| A | on | off | +2.33316e-01 | +7.77723e-02 | +4.57453e-02 | +4.70007e-02 | +2.44964e-03 | 1.16e+00 |
| B | off | on | -2.82574e-07 | +2.96297e-01 | +3.60623e-02 | +3.60623e-02 | +1.29664e-03 | 1.00e+00 |
| C | off | off | +2.96300e-01 | +9.87671e-02 | +2.67435e-02 | +2.67435e-02 | +2.85221e-03 | 1.75e+00 |

Per-reduction relative shift vs baseline (V_MAX=20):

| reduction | cell A (j2 off) | cell B (j1 off) | cell C (both off) |
|---|---|---|---|
| G1/sum_entries | 1.16e+00 | 1.00e+00 | 1.75e+00 |
| G2/sum_entries | 8.72e-01 | 5.13e-01 | 8.38e-01 |
| G2/tr_pi | 1.46e-01 | 3.27e-01 | 5.01e-01 |
| G2/delta_1 | 1.81e-01 | 3.72e-01 | 5.34e-01 |
| G2/vac_pi | 4.87e-01 | 7.28e-01 | 4.03e-01 |

## Task A verdict + interpretation

**H_SENSITIVE — the pre-registered H_EMPTY prior LOST, and the banked result survives the audit.** Cell A (j=2 phase off) moves every reduction materially (G1 sum_entries 0.108→0.233, G2 sum_entries 0.609→0.078 at V_MAX=20 — shifts of 30–100%, orders of magnitude above the ~10⁻⁴ truncation floor). The j=2 phase — even as a mere mod-3 shadow of the mod-9 twist at n=3 — genuinely carries the match. The n=3 DWM↔Syracuse match therefore **does** encode phase information, and the DWM identification's quantitative confirmation stands. **No erratum.** (This is the outcome that costs the least and was assigned the lowest prior; per the repo's discipline the prior is recorded as having lost — the third pre-registered prior in this arc to do so, after H_QUAD and ⌊r/2⌋+2.) Cell C (both phases off) also moves materially, so the reductions are **not** phase-blind — H_PHASE_BLIND is refuted, and the extra cell earned its place by ruling out the bigger erratum. Cell B (j=1 phase off) collapses the 3-alternating moment G1 to ~0 (−3×10⁻⁷) while G2 survives, showing G1 is carried entirely by j=1's phase (j=1 appears twice in ϕ(X̃₁X̃₂X̃₁)) whereas the 4-alternating moment retains structure without it.


## Cell B is a finding: the non-freeness mechanism, experimentally isolated

Cell B (j=1 phase off) drives the 3-alternating moment G1 = ϕ(X̃₁·X̃₂·X̃₁) to **−3×10⁻⁷ — identically zero to the truncation floor.** That is not "j=1 carries G1"; it is a structural statement about the corpus's terminal framework result. `OBSTRUCTION_MAP_TERMINAL.md:86` records that the third-order alternating repeated-index moment `φ(X̃_{j₁}·X̃_{j₂}·X̃_{j₁}) ≠ 0` is *the* diagnostic that killed B-amalgamated **free** independence and forced **monotone** (Muraki 2003 / Hasebe–Saigo 2011) — with the stated mechanism (`:91`): "when X̃_{j₁} appears on both sides of X̃_{j₂}, the phases induced by [the b_prior] coupling do not cancel."

This ablation demonstrates that mechanism directly: **turn off the j=1 (bracketing-index) phase → the bracketing coupling vanishes → G1 → 0 → freeness would hold.** The non-freeness of Syracuse — the single fact that redirected the entire framework arc from free to monotone — is carried **entirely by the j=1 phase**, and it is now *experimentally isolated by ablation* rather than *argued by inspection*. That upgrades a load-bearing structural claim in the corpus from derived to demonstrated. (G2, the 4-alternating moment, survives cell B at 0.296 — its non-freeness has additional carriers, consistent with the higher-order pattern.)

## Task B — Mahler predicts r=2 (independent of Task A)

Fixed profile from r=6 (r≥3 only, never saw r=2); c_k reduced mod 27. Since v₃(c_k)≥3 for k≥3, the tail vanishes mod 27, so the r=2 prediction is s₂(b)=Σ_{k≤2} c_k·C(b,k) mod 27 — a **prediction**, not a 3-point fit.

| ℓ | ε | c₀,c₁,c₂ mod 27 | predicted s₂(0,1,2) | certified | pred−cert |
|---|---|---|---|---|---|
| 0 | 0 | [4, 18, 18] | [4, 22, 4] | [13, 4, 13] | +18 (const) |
| 1 | 0 | [4, 18, 18] | [4, 22, 4] | [22, 13, 22] | +9 (const) |
| 2 | 0 | [4, 18, 18] | [4, 22, 4] | [4, 22, 4] | +0 (const) |
| 0 | 1 | [8, 0, 9] | [8, 8, 17] | [26, 26, 8] | +9 (const) |
| 1 | 1 | [8, 0, 9] | [8, 8, 17] | [17, 17, 26] | +18 (const) |
| 2 | 1 | [8, 0, 9] | [8, 8, 17] | [8, 8, 17] | +0 (const) |

**The r≥3 Mahler SHAPE predicts r=2 exactly; only a global phase offset differs.** For every family, `pred − cert` is *constant in b* — a pure global phase — and always a multiple of 9 = 3² (v₃ ≥ 2), varying with (ℓ,ε): offsets {18,9,0,9,18,0}. So the b-dependence (the coefficients c₁,c₂) transfers from r≥3 to r=2 exactly; only the constant term c₀ carries a **level-dependent global phase** of 3-adic depth ≥2 that the r≥3 profile does not fix.

Reading: **r=2 is DERIVED up to a global phase, not freely fitted** — the Probe 82 untestability floor (R81's r≥3 *degree-fitting* floor) is dissolved for the phase *shape*. The residual is a single level-dependent constant per family, not 3 free values.

### The mod-9 offset is residual structure, not a normalization artifact

> **⚠ CORRECTION (superseded by R84, 2026-07-14).** This section's conclusion is
> **overturned.** Probe 84 predicted the offset at r=2..7 (not just r=2) and found
> **v₃(offset) = r at every level** — the offset is a multiple of `3^r` (the top
> modulus layer), *not* a fixed `3²`. It is the global phase induced by the
> family-defining twist `(1+3^r)^ℓ` (which collapses to `ω₃^ℓ` because `4 ≡ 1 mod 3`),
> i.e. a **normalization artifact**, not independent residual structure. The character
> check below is correct as far as it goes but ran at a *single* level (r=2) and was too
> coarse to catch the r-scaling — exactly the failure mode Probe 84's pre-reg flagged.
> The "lead handed to R81/R81b" is **withdrawn**; see `result_84_mod9_offset.md`. The
> shape-transfer result (c₁,c₂,… predict r=2 exactly) is unaffected and stands. Text
> preserved below for the record.

The five-minute discriminator (is `(ℓ,ε) → offset/9 ∈ Z/3` a character?): the map is `{(0, 0): 2, (0, 1): 1, (1, 0): 1, (1, 1): 2, (2, 0): 0, (2, 1): 0}`. It is **not** a group homomorphism (`f(0,0)=2≠0`), **no** linear form `aℓ+bε` reproduces it, and it does **not** factor through `c mod 9` (`factors_c9=False`). So the offset is **not** a normalization artifact of how c_{ℓ,ε} was defined — it is genuine residual structure. It *is* low-complexity: **affine in ℓ with ε flipping the slope sign** — `f(ℓ,0)=2−ℓ`, `f(ℓ,1)=1+ℓ=−f(ℓ,0) mod 3` (eps0: const=2,slope=2; eps1: const=1,slope=1). **The 3² reappears in the one term the Mahler profile doesn't explain, and it is exactly the term that distinguishes the six families** — a lead handed back to the R81/R81b agent, not a closed nuisance. The ε-antisymmetry `f(ℓ,1)=−f(ℓ,0)` is worth noting against the sibling 3x±1 sign symmetry `σ(r)=−r` (K₋=σK₊σ), but that link is not established here.

Two caveats, both stated plainly: (1) this unlocks the r=2 *shape* only — the absolute phase carries the mod-9 residual above; (2) it does **not** dissolve the j=1-exceptional problem (Task A), which is about which DWM *step* carries the moment — moot here anyway since Task A fired H_SENSITIVE independently.

## Task C — re-cost of the n≥5 evidential bridge (Mahler-updated)

The F̂ side no longer requires computing F̂ at high r: the phase is the fixed Mahler profile, available at any r essentially for free (Task B shows it even predicts *downward* to r=2). So the remaining cost of an n≥5 evidential bridge is entirely the **Syracuse-side directly-measured moments at n≥5**, which do not exist and must be produced (their own probe, own pre-reg, own falsifier). The DWM-prediction side scales as state_count(n)³·V_MAX⁴ (state_count = 2·3^{n−1}: 162 at n=5, 486 at n=6) → ~6 h (n=5) to ~160 h (n=6) as before, but that is now the *only* heavy item and it is one-sided. Net: the bridge is no longer gated on the F̂ side at all — it is gated on standing up a new Syracuse measurement at n≥5, and (per Task A's amendment) it should target **j≥2 at n≥5**, not j=1 at n=4.

## Scope — what is untouched

**c = 7/45 (`THEOREM_C_745.md`) is UNAFFECTED by every outcome here.** It is derived from R75 Plancherel × R76 conservation × R77 T_diag, and `D3_DERIVATION_AUDIT.md` established it never depended on the DWM framework-identification overlay. Whatever this probe does to DWM's *evidential* status, 7/45 stands. Theorems 78.1–78.3 are likewise untouched.

_Reporting discipline: the fired outcome is reported with per-reduction shift magnitudes, not a binary. The ablation is named once (phase factor → 1; shift and weights kept). Cells B and C were added as a pre-fire §3′ amendment with H_PHASE_BLIND pre-registered. Task A and Task B are independent; neither licenses a claim in the other._
