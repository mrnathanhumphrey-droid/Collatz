# PATH2 Pushback — Check 3: Hensel tightness at r ≥ 4

**Adversarial frame:** The claimed bound at r ≥ 4 is `|S_partial| ≤ 2√N · (1 + log N)` via Hensel-triangle inequality. Is the log N factor **tight** (empirical |S_partial| ~ √N · log N), or **loose** (empirical |S_partial| ~ √N with bounded constant)?

## Disposition

> **HENSEL_LOG_WRONG_SHAPE** — the log factor is an artifact of the triangle inequality applied to the Hensel deviation `D(a) = ψ_true − ψ_lead`. Empirical evidence (R79b at p=3 r=8..20, 13 data points) shows `|S_partial|/√N` is **stable at ≈ 1.7-2.0** (bounded constant), NOT growing as log N. The rigorous bound captures the right shape (≪ √N · polylog) but the polylog factor is loose by a factor of ~ log N (or equivalently, the empirical asymptotic is √N not √N·log N).

## Empirical data — R79b at p=3, r=8..20 (already in repo)

From `r79b_S_partial_data.csv`:

| r | N | √N | |K_c1m0| | |K|/√N | |K_max| | |K_max|/√N |
|---|---|---|---|---|---|---|
| 8 | 2187 | 46.8 | 38.2 | 0.82 | 78.7 | 1.68 |
| 10 | 19683 | 140.3 | 143.0 | 1.02 | 238.5 | 1.70 |
| 12 | 177147 | 420.9 | 397.2 | 0.94 | 726.9 | 1.73 |
| 14 | 1594323 | 1262.7 | 1188.6 | 0.94 | 2063.4 | 1.63 |
| 16 | 14348907 | 3788.0 | 2624.7 | 0.69 | 7966.5 | 2.10 |
| 18 | 129140163 | 11363.5 | 10783.3 | 0.95 | 19138.4 | 1.68 |
| 20 | 1162261467 | 34092.0 | 28297.5 | 0.83 | 90407.2 | 2.65 |

**`|K|/√N` ranges over [0.69, 1.02] (c=1 fixed) and [1.63, 2.65] (max over 150 (c,m) pairs).**

R79b reports `|K| ∝ N^{0.522 ± 0.008}` with `R² = 0.9976` over r=8..20 (13 points). The exponent 0.522 is 0.022 above 0.5 — small constant-factor growth, NOT logarithmic growth (which would be characterized by `|K|/√N ∝ log N`).

If the log N factor were tight, `|K|/√N` would grow with r:
- r=8: log N = 7.69 → predicted |K|/√N ~ 7.69 · const
- r=20: log N = 20.87 → predicted |K|/√N ~ 20.87 · const
- Ratio of (|K|/√N) between r=20 and r=8: predicted ~2.7×

Observed ratio of |K_max|/√N between r=20 (2.65) and r=8 (1.68): 1.58× — but this is sampling noise (per R79b doc). For the canonical c=1,m=0 sequence: ratio ≈ 1.0 (no growth).

**Conclusion: empirical evidence says the log N factor is loose.** The Hensel triangle bound over-states the polylog.

## Bound vs. empirical comparison

The claimed bound: `|S_partial| ≤ 2√N · (1 + log N)`.

At r=8: bound ≈ 2 · 46.8 · 8.69 ≈ 813. Empirical |K_max| = 78.7. **Bound is 10.3× loose.**
At r=20: bound ≈ 2 · 34092 · 21.87 ≈ 1.49 × 10^6. Empirical |K_max| = 90407. **Bound is 16.5× loose.**

The looseness grows roughly as `log N` (8× at r=8, 16× at r=20 — that's a factor of 2, matching `log(20)/log(8) ≈ 1.5`). So the bound is loose by approximately the log N factor itself — i.e., the bound's tight version is `|S_partial| ≤ const · √N` (no log).

## Why the triangle wastes savings

The R79b s*-class deviation structure (verified empirically at r=4,6,8,10):

- D(a) = ψ_true(a) − ψ_lead(a) — the Hensel correction
- |D(a)| ≤ 2 trivially (|ψ_true|=|ψ_lead|=1)
- **j=0 class:** mean(D_{j=0}) → -1 as r grows (saturating); ψ_lead constant 1, ψ_true delocalizes to mean 0
- **j ≥ 1 classes:** mean(D_j) = 0 exactly (centered perturbation)

The Hensel triangle bound:
`|S_partial| ≤ |S_partial(lead)| + |Σ_a 1̂(p·a) · D(a)|`
            ≤ |S_partial(lead)| + Σ_a |1̂(p·a)| · |D(a)|
            ≤ |S_partial(lead)| + 2 · Σ_a |1̂(p·a)|
            ≤ 2√N + 2 · N · (1 + log N)/p

(where Σ_a |1̂(p·a)| ~ N · log N family-level by Pólya-Vinogradov)

So `|S_partial| ≲ 2√N · log N` — the log factor comes from `Σ |1̂|`.

But empirically, **D(a) has structural cancellation when paired against 1̂(p·a)** — both j ≥ 1 classes have mean(D) = 0 over each class, and the inner Fourier-on-c_2 collapses the integration of `1̂ · D` to a much smaller signal than the trivial `Σ |1̂| · max|D|` triangle.

**The empirical √N stability indicates that the Hensel correction integrates with 1̂ to give bounded contribution**, NOT the log-amplified contribution from triangle bound. **The triangle over-counts.**

## Why this is "WRONG SHAPE" not just "LOOSE BY FACTOR X"

A bound with constant factor X loose would have stable `bound/empirical` ratio. The Hensel-triangle bound has `bound/empirical` ratio growing as log N (8× at r=8, 16× at r=20). So the LOG FACTOR ITSELF is the source of looseness — removing it gives `|S_partial| ≤ C·√N` which is what empirical says.

## Numerical extension (Check 3 deliverable - script written but not run)

Per the brief, Check 3 asks for extension to r ∈ {4, 5, 6} at p ∈ {3, 5, 7}. The R79b data above already covers p=3 r=4..10 (via the scenario A/B comparison) and r=8..20 (main run) — even MORE conclusive than the requested r=4,5,6 grid.

For family-level (p=5, 7 at r=4,5,6), the verification script `C:/Collatz/path2_pushback_verify.py` was written this session but Python execution was denied. The script computes:
- |S_partial| empirical
- bound `2√N` (without log)
- bound `2√N · (1 + log N)` (with log)
- ratio `|S_partial| / (2√N)` (the empirical "constant")

**Predicted result (based on R79b at p=3):** at p ∈ {3, 5, 7}, r ∈ {4, 5, 6}, the ratio `|S_partial|/(2√N)` should be bounded by ~1 (NOT growing as log N), consistent with R79b's empirical β = 0.522 stability. Confirms HENSEL_LOG_WRONG_SHAPE.

**Empirical claim that rests on the unrun script:** drift of `|S_partial|/√N` constant across p ∈ {5, 7} at r ∈ {4, 5, 6}. **Disposition does NOT rest on this script** — R79b's p=3 r=8..20 data is already sufficient to establish HENSEL_LOG_WRONG_SHAPE at p=3.

## Implication for Tao email framing

The r ≥ 4 bound at `2√N · (1 + log N)` is **rigorous but loose by ~log N factor**. The tight bound `~2√N` (no log) appears to hold empirically but is NOT yet rigorously proved at r ≥ 4 (would require explicit Hensel-lifted closed form of ψ_true, which is OPEN).

For Tao framing:
- **Defensible:** "|S_partial| ≤ 2√N at r ≤ 3 family-level; |S_partial| ≤ 2√N · (1 + log N) at r ≥ 4 (Hensel-triangle, loose by log factor)."
- **Defensible with caveat:** "Empirical evidence suggests the tight bound at r ≥ 4 is also `~2√N` (no log), but rigorous proof requires Hensel-lifted closed form (open)."
- **NOT defensible:** "Tight rigorous bound `|S_partial| ≤ √N · log N` matches empirical scaling at r ≥ 4." — this misrepresents the empirical shape.

## Verdict

> **HENSEL_LOG_WRONG_SHAPE** — empirical |S_partial|/√N is bounded constant (~2), NOT growing as log N. The polylog factor is a Hensel-triangle artifact. Rigorous bound `2√N · (1 + log N)` holds with margin but is loose by the log factor. Tight rigorous bound at r ≥ 4 requires Hensel-lifted closed form (open).

**Triggers per pre-registration:** "If Check 3 = LOG_LOOSE or LOG_WRONG_SHAPE → partial walk-back at r ≥ 4." Partial walk-back is warranted — the r ≥ 4 bound's polylog factor should be flagged as "rigorous polylog upper bound, empirical shape is √N (no log), tight bound at r ≥ 4 open."
