# C3: Band-l¹ analysis on the dangerous band — outcome (saturated, path closed)

**Date:** 2026-05-04. Companion to R79b empirical study and saddle_class_subsum_analysis.md. Closes the third sub-path of Path C ("direct band-l¹ analysis") for Kalafatelis 2026 eq 190.

## Verdict: **eq 190 cannot be closed via pointwise band-l¹ analysis at observed r.**

The metric `N_r^{-1/2} · ‖ĥ_{r,ℓ}‖_{ℓ¹(D_{r,t}(η))}` (which eq 190 requires to be `≪ η^{1/2+δ}`, i.e., uniformly bounded in r) is empirically observed to **grow exactly linearly in N_r** across every (ℓ, t, η) combination tested:

| η | mean β (slope of log metric vs log N_r) | std | n |
|---|---|---|---|
| 0.5 | **1.0003** | 0.0021 | 9 |
| 0.25 | **0.9996** | 0.0022 | 9 |
| 0.125 | **1.0008** | 0.0032 | 9 |
| 0.0625 | **1.0005** | 0.0060 | 9 |

R² ≈ 1.0000 in all 36 individual fits. Slope = 1.0 with virtually no variance across ℓ ∈ {0, 1, 2}, t ∈ {0, π/2, π}, η ∈ {1/16, 1/8, 1/4, 1/2}.

**The dangerous band exhibits no inter-m cancellation.** `‖ĥ‖_{ℓ¹(D)}` saturates at the trivial `|D| · sup|ĥ| ≈ √η · N_r^{3/2}` bound. Equivalently: `|ĥ(m)|` is approximately uniformly `√N_r` over m ∈ Z/N_r (Plancherel-saturated), including on the dangerous band — no concentration to exploit, no phase cancellation to harvest.

## Setup

Per Kalafatelis 2026 §4.4–4.5:
- q = 3^{r+1}, N_r = 2·3^{r-1}, ω_r := 2^{N_r} mod q  (primitive cube root of unity mod 3^{r+1})
- h_{r,ℓ}(j) := e_q(ω_r^ℓ · 2^j),  j ∈ Z/N_r,  ℓ ∈ {0, 1, 2}
- ĥ_{r,ℓ}(m) := Σ_{j=0}^{N_r-1} h_{r,ℓ}(j) · e^{−2πi mj/N_r}  (forward DFT)
- m_r,t(θ) := Σ_{b=1}^{N_r} 2^{−b} e^{ib(t−θ)}  (multiplier)
- D_{r,t}(η) := { m ∈ Z/N_r : |m_r,t(2πm/N_r)| > 1−η }  (dangerous band)
- |D_{r,t}(η)| ≪ √η · N_r + 1 (Prop 22)

Eq 190 conjecture (Kalafatelis Prop 24 input):
$$
\sup_{r \geq 1, \ell, J: |J| \ll \sqrt{\eta}\,N_r} N_r^{-1/2} \cdot \|\hat{h}_{r,\ell}\|_{\ell^1(J)} \ll \eta^{1/2+\delta}
$$

The LHS must be **bounded uniformly in r** at fixed η, decaying to 0 as η → 0. Empirical r-slope of the metric reveals which regime holds.

## Reference scalings

| Regime | metric scaling vs N_r | β (slope of log metric vs log N_r) | Eq 190? |
|---|---|---|---|
| Trivial / no cancellation | N_r¹ | 1.0 | NO |
| Cauchy-Schwarz on band (same as trivial up to η factor) | N_r¹ | 1.0 | NO |
| Square-root cancellation in m within band | N_r^{1/2} | 0.5 | INSUFFICIENT |
| Eq 190 holds | N_r⁰ (bounded) | 0.0 | YES |

Empirical: β = 1.0 across all configurations. Trivial/saturated.

## Per-(ℓ, t, η) breakdown

All 36 sweeps gave β within [0.989, 1.010], R² ≥ 0.9999. The full table is in [band_l1_log.txt](band_l1_log.txt) and [band_l1_data.csv](band_l1_data.csv).

Examples:

| ℓ | t | η | β | R² |
|---|---|---|---|---|
| 0 | 0 | 0.5 | 0.9982 | 1.0000 |
| 0 | π/2 | 0.25 | 1.0024 | 1.0000 |
| 0 | π | 0.125 | 0.9982 | 1.0000 |
| 1 | π/2 | 0.0625 | 1.0093 | 0.9999 |
| 2 | 0 | 0.5 | 1.0029 | 1.0000 |

The fit is **incredibly clean** — log metric is essentially linear in log N_r over r ∈ {6, 8, 10, 12, 14}, with slope exactly 1, intercept depending only on η.

## η-dependence at fixed r

Intercepts α(η) (where log metric = α + 1.0 · log N_r) follow approximately α ≈ 0.79 · log η + const, suggesting `metric ≈ η^{0.79} · N_r · const`. The η-power 0.79 is **between** the random model (η^{0.5} = √η, expected from |D| ≈ √η·N_r and uniform |ĥ|) and the saturated bound (η^1, expected if |ĥ| concentrated on the band peak).

This intermediate η-power 0.79 indicates ĥ has **slight concentration** on the dangerous band — the band peak attracts marginally more |ĥ| mass than uniform — but **not enough to save** the N_r^{1} growth.

## Why the saturation is decisive

Plancherel: `Σ_m |ĥ(m)|² = N_r · ‖h‖²_{ℓ²} = N_r · N_r = N_r²` (since |h(j)| = 1).

Average |ĥ(m)|² over m ∈ Z/N_r: **N_r exactly**. So `|ĥ_typical| ≈ √N_r`.

On a band of size |D| ≈ √η·N_r:
- `Σ_{m ∈ D} |ĥ(m)|² ≈ |D| · N_r ≈ √η · N_r²` (uniform mass on band)
- `Σ_{m ∈ D} |ĥ(m)| ≈ |D| · √N_r ≈ √η · N_r^{3/2}` (no inter-m cancellation)

Then `N_r^{-1/2} · Σ ≈ √η · N_r`, which grows like N_r^{1.0}. Matches our empirical β = 1.0.

To break this and get bounded-in-r: would need `|ĥ(m)|` distributionally NOT uniform over m, OR phase cancellation between m-values when summing in |·|. Empirically neither occurs at observed r.

## Implications

**For c=7/45 closure via Path A**: the saddle-class partition (saddle_class_subsum_analysis.md) and direct band-l¹ analysis (this) are both closed as closure paths. Combined with Cochrane (R78) and van der Corput (R79) being closed, this exhausts all "elementary" approaches to bounding |S_partial|.

**For Path C remaining sub-paths:**
- C1 (5x+1 sibling-attack reframing): not yet attempted — significant rebuild required
- C2 (Bourgain-Konyagin sum-product on ⟨4⟩): the BK technique gives bounds like `|Σ_{x ∈ H} ψ(x)| ≤ |H|^{1−δ}` for multiplicative subgroups H ⊂ (Z/p)^×. For our setup, H = ⟨4⟩ (principal-units mod 3^{r+1}, size 3^r) and the sum is over the support of F̂. **BK could give a true polynomial saving on the primal Kalafatelis sum directly** — but this is a different attack route than band-l¹ on ĥ. C3 closure does not block C2.

**For framework integrity**: the C3 saturation is consistent with R79b's β = 0.522 finding (square-root scaling of |K| against N) — both reflect that ĥ is "Plancherel-uniform" with no concentration. No surprise; just a clean confirmation from the dual side.

## Honest caveats

1. **Five r-points only**: r ∈ {6, 8, 10, 12, 14} gives 5 points per fit. Extending to r=16 is feasible (~30 sec) and would tighten R². But β = 1.000 ± 0.006 across 36 fits leaves no room for ambiguity.
2. **Sample of t ∈ {0, π/2, π}**: doesn't cover all dangerous-band locations. But the multiplier modulus `|m_r,t|` is symmetric under t ↔ -t and quasi-periodic, so {0, π/2, π} is a reasonable coverage of qualitatively distinct band positions.
3. **No alternative band definitions tested**: Kalafatelis's specific D_{r,t}(η) is the structural target. Other definitions (e.g., level sets of |ĥ| itself) might exhibit different behavior — but those wouldn't connect to eq 190's framework.
4. **Inter-m cancellation might appear at much larger r**: r → ∞ is an unproven asymptotic. Conceivable that some non-trivial structure emerges at r ≥ 20+ but not at r ≤ 14. Unlikely given the cleanness of β = 1.0 fit, but not formally excluded.

## Files

- `band_l1_analysis.py` — measurement script
- `band_l1_data.csv` — 540 rows = 5 r × 3 ℓ × 3 t × 4 η × 9 (per-fit data within sweep, actually 5×3×3×4=180 rows in the simple table)
- `band_l1_log.txt` — full stdout
- `band_l1_analysis.md` — this writeup

## Decision tree resolution (continuing R79b's)

Per the morning consolidation framing established in R79b:

- **Path A (closed)**: Cochrane (R78), van der Corput (R79), saddle-class partition (saddle_class_subsum_analysis.md), direct band-l¹ analysis (this).
- **Path B (fallback)**: consolidate-and-publish c=7/45 with empirical β=0.522, structural anchors T78.1–T78.6, and the now-comprehensive obstruction map.
- **Path C remaining**: C1 (5x+1 sibling) untested; C2 (Bourgain-Konyagin sum-product on the primal sum) unattempted.

C3 outcome adds a fourth "elementary closure path closed" entry to the obstruction map. The closure of c=7/45 either requires:
- Sum-product / additive-combinatorics machinery (BK level), OR
- A genuinely new approach not in the dual-Fourier or Cochrane-derived families.

The framework's structural anchors (T78.1–T78.6) remain intact; the empirical certification (β = 0.522, |ε_n|·2^n ≤ 0.04 through k=6) remains intact; the rate-1/2 rigorous proof remains the single open problem.
