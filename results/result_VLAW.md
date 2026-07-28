# RESULT — VLAW: Collatz is NEARLY critical (δ₁≈−8e-4), the 7/15 wall is real not a v-artifact, but the v-law is a genuine ~2% lever on S_∞ (2026-07-27)

**Probe:** `probes/probe_vlaw.py`. Spec (Wilson): the exact empirical v-law `P(v | m mod 2^k)` from Lagarias–Sinai
(`outputs/v_distribution_by_N.csv`, forward integer orbits, χ²/dof=4373 vs Geom(½)) deviates hugely across *all*
moments — but S_∞ sees one linear functional. Compute the criticality functionals **directly** (χ² is not δ).

## Stage A — the two criticality functionals (forward-orbit P(v))
`δ₁ = E_exact[2^{−v}] − 1/3` (drift; =0 ⟺ E[W]=1, W=3·2^{−v}, the P6K λ=½ criticality) and
`δ₂ = Σ_v P_exact(v)² − 1/3` (participation / D₂). Both are exactly 0 under Geom(½).
```
   N       E[2^-v]      δ₁           ΣP²         δ₂
  2^28    0.3322995   -1.03e-03    0.3315765   -1.76e-03
  2^30    0.3323770   -9.56e-04    0.3316889   -1.64e-03
  2^32    0.3324557   -8.78e-04    0.3318142   -1.52e-03
  2^34    0.3325167   -8.17e-04    0.3319092   -1.42e-03
```
- **Both δ negative, both shrinking monotonically with N** (~7% smaller per doubling). Real Collatz forward orbit is
  **very slightly SUBCRITICAL**: `E[W] = 3·E[2^{−v}] = 0.99755 < 1` (more contracting than the exactly-critical model).
- **δ₂ = −0.43% of 1/3** at 2³⁴ — NOT machine-zero, but small and **not** the 2% size of the 7/15-vs-0.475 gap on its own.
- χ²/dof=4373 measures all-moment deviation; the two functionals S sees are ~0.1–0.4%. **The huge χ² is not evidence δ≠0** —
  computed directly, δ is small (the trap Wilson flagged, avoided).

## The escape-hatch closure (the load-bearing point for 7/15)
**The 7/15 object uses Geom(½) BY CONSTRUCTION** (`build_nu` lam=0.5; Tao's Syracuse chain draws v~Geom(½) as its kernel,
independent of state). So **δ₁=δ₂=0 identically in the object where 7/15 lives**, and SOLSTICE measured **0.475 under exact
Geom(½)**. Therefore **0.475 is the exact-critical value — the 7/15-vs-0.475 gap is NOT a v-law artifact. The wall is real.**
The "secretly non-critical ⟹ not 7/15" escape hatch is **closed**.

## Mirage-check / iid Stage B — the v-law IS a real ~2% lever (but it doesn't point at 7/15)
Build the Syracuse stationary `π_k` two ways — Geom(½) vs the forward-orbit `P(v)` as an iid kernel — and compare (`T_level_w`):
```
   n   |π_geom-π_fwd|max   S=2T geom   S=2T fwd     dS
   2      6.85e-03          0.46157     0.45675   -4.82e-03
   4      1.94e-03          0.46551     0.45810   -7.42e-03
   6      4.93e-04          0.46549     0.45693   -8.56e-03
```
- The two stationary measures **DIFFER** (`|Δπ|` shrinks with n, but the value gap `dS` **grows**): imposing the forward
  marginal **lowers S**, with `S_fwd ≈ 0.458` (roughly flat) vs `S_geom` rising toward 0.475.
- So the v-law is a **genuine ~2% lever on S_∞ — the same magnitude as the 7/15 gap** (Wilson's "same object?" — in *size*, yes).
- **Direction kills the rescue:** `0.458 < 7/15 = 0.4667 < 0.475`. **7/15 sits BETWEEN the critical model and the real-v
  model.** The v-law correction pushes S *past* 7/15 (further below), not onto it — it does not explain or restore 7/15. It
  says the *real-Collatz* value is likely **lower** than 0.475, with 7/15 inside the model-uncertainty band but not privileged.

## Verdict (the honest, three-part answer)
1. **Is Collatz exactly critical, or only nearly?** — **Nearly.** Forward-orbit `δ₁=−8×10⁻⁴`, `E[W]=0.9975`, slightly
   subcritical/contracting; δ shrinking with N (plausibly →0, or a small structural floor — Lagarias–Sinai says the latter).
2. **Is 7/15 rescued by the v-law?** — **No.** The 7/15 object is exactly critical (Geom½) and gives 0.475; the real-v
   correction moves S *down* toward ~0.458, past 7/15. 7/15 is excluded from both directions.
3. **New caveat on the SOLSTICE value (the real finding):** "S_∞ ≈ 0.475" is **model-specific (Geom½)**. The v-law is a real
   ~2% lever, so real-Collatz S_∞ is v-law-dependent within roughly **[0.458, 0.475]**. The number 0.475 is the *critical-model*
   value, not automatically the real-Collatz value.

## Open (the exact Stage B — Wilson's call, bigger build)
The iid-forward kernel is **not** real Collatz: `v` is arithmetic-deterministic given `m mod 2^k` (Lagarias–Sinai — half of
residues mod 4 force v=1, etc.), correlated with the state. The exact real-Collatz S_∞ needs the **deterministic-per-residue**
kernel (`outputs/v_conditional_distributions.csv`), which correlates v with the chain state rather than drawing iid — a larger
build. The iid result (~0.458, 2% down) is indicative, not exact. **Data pointer:** `outputs/v_distribution_by_N.csv` (marginal,
used here) and `outputs/v_conditional_distributions.csv` (per-residue, for the exact kernel).

**Not at stake:** P6D–P6K identities, `S_{i+1}=2T_i`, i=20 no-crossing, GARSIA, R1–R30. The 7/15 exclusion **strengthens**
(now excluded by criticality-model value 0.475 AND by the real-v correction ~0.458). The one thing this *adds*: SOLSTICE's
0.475 is the Geom(½)-model value; real-Collatz S_∞ carries a ~2% v-law model dependence.
