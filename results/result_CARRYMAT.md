# RESULT — PROBE CARRYMAT: the per-cell angular matrix (Wilson's amended Option 1) (2026-07-26)

**Probe:** `probes/probe_carrymat.py`. Move the sign into the class that can see angle-sign — a quadratic form.
In the 2-D mean-zero plane of ℤ/3 (orthonormal basis), per cell c (rotation cm=c%3), form the μ-averaged
`S_c = sym E_μ[ P(T_c δπ_{4x}) ⊗ P(δπ_x) | c ]`; `U_c = p_c·tr(S_c)`, `Σ_c U_c = q_r(1)-1/3` (GATE).
Wilson's amendment: report eigenvalues+**eigenvectors**, the observed direction histogram, and the observed-mass
overlap with S_c's **negative eigendirection** — because we need positivity on the *observed cone*, not PSD.

**GATE PASS:** `Σ_c p_c tr(S_c) = +4.178936e-4` vs banked `q_16(1)-1/3 = +4.1789e-4`, rel `8.5e-6`. Reframe exact.
(Note: an earlier run's `cos_θ`/`D` column was corrupt from a 2-D rotation-sign shortcut; fixed to compute the
defect in 3-vector space = CARRYCOV's D. Numbers below are the clean rerun; the S_c matrices were correct throughout.)

## Two confirmations for Wilson (the κ-fix and the "one event")

**cos_θ_true rises — the κ=1 shortcut was the whole discrepancy.** With `κ²=VAR_L/VAR` and
`cos_θ_true=(1+κ²−D²)/(2κ)` (vs the geometric `1−D²/2` which sets κ=1):

| r | κ² | D | 1−D²/2 | cos_θ_true | q−1/3 | U0+U3 |
|---|-----|-----|--------|-----------|-------|-------|
| 8  | 0.9872 | 1.3873 | 0.0377 | 0.0315 | 1.011e-3 | −1.39e-3 |
| 11 | 1.0012 | 1.3904 | 0.0333 | 0.0339 | 6.93e-4 | −2.16e-4 |
| 12 | 1.0050 | 1.3902 | 0.0337 | 0.0361 | 6.40e-4 | +3.25e-5 |
| 16 | 1.0189 | 1.3921 | 0.0311 | 0.0401 | 4.18e-4 | +1.09e-4 |

`cos_θ_true` tail r=9→16: `0.031→0.040`, **monotone rising**, while `D` sits rock-stable 1.383–1.392 (~1.5% below √2,
angle ~88°). The geometric `1−D²/2` *falls* only because it drops the growing `(κ²−1)/2` term. **The right variable
has a stable/rising margin** — the conditioning gain is real (a constant lower bound is possible on ~0.04; impossible
on something decaying at 0.89). It remains a reparametrization, not new evidence — but a far better-conditioned one.

**κ² crosses 1 at r=11; U0+U3 crosses 0 at r=12 — same window.** Wilson's "might be one event" is **supported**:
the profile-spread asymmetry flipping (dpL more spread than dpx, κ²>1) coincides with the rotation-0 group turning
positive. A concrete single-mechanism lead for the pen.

## The matrix reframe: exact and correctly-typed, but it does NOT reduce at q=3

Per-cell S_c at r=16 (`iso = tr/2` = the isotropic drag; anisotropy = `(λ+−λ−)/2`):

| cell | rot | p_c | U_c | iso (drag) | λ+ @ang | λ− @ang | f_neg | obs axis |
|------|-----|-----|-----|------------|---------|---------|-------|----------|
| c=0 | 0 | 0.247 | **−6.4e-5** | −0.0001 (+1) | +0.0001 @156° | −0.0004 @66° | 0.503 | −155° |
| c=1 | 1 | 0.262 | +2.04e-4 | +0.0004 (−½) | +0.0030 @150° | −0.0022 @−120° | 0.473 | +109° |
| c=2 | 2 | 0.243 | +1.05e-4 | +0.0002 (−½) | +0.0003 @46° | +0.0001 @136° | 0.509 | +104° |
| c=3 | 0 | 0.248 | +1.73e-4 | +0.0003 (+1) | +0.0007 @149° | −0.0000 @−121° | 0.470 | +127° |

Three findings, none a kill, all honest negatives on the *mechanism*:

1. **PSD is false where Wilson expected — but the −½ drag it feared is ABSENT.** c=0 has net-negative trace, c=1 is
   indefinite (λ−=−0.0022). Yet `iso ≈ 0` in every cell (−0.0001…+0.0004), **not** the `cos(120·cm)=−½` of cells
   1,2. The isotropic drag only exists if `M ≈ I` (profile carried similarly) — and CARRYCOV already showed `M ≉ I`
   (`D≈√2`, near-independent). So **there is no −½ to beat**; Wilson's "V1 > 0 beats the −½ of mere similarity" is
   moot. V1>0 is carried by S_1's eigenvalue asymmetry (λ+ @150° exceeds |λ−| @−120°), not by defeating a drag.

2. **The observed direction distribution is ~ISOTROPIC — no cone-avoidance mechanism.** `f_neg` (observed direction
   energy on S_c's negative eigendirection) = 0.503 / 0.473 / 0.509 / 0.470 — all within 0.03 of the isotropic 0.5;
   the arg-histograms are spread, no sharp peak. By Wilson's own gate — "if observed mass avoids the negative
   direction, that avoidance IS the mechanism" — **the mass neither avoids nor sits on it; it is diffuse.** So the
   sign collapses to `sign(tr S_c)` = whether λ+ outweighs |λ−|, which is the original coherence scalar. **No
   reduction.** This is the *third* exact, correctly-typed reframe (D, covariance, now matrix) that repackages the
   signed coherence without lowering it: isotropic directions + indefinite matrix ⟹ sign = sign(trace).

3. **One geometric regularity:** the positive eigendirection λ+ clusters at ~150° across cells 0,1,3 (156/150/149°),
   with observed axes ~30–45° off (104–127°). Real structure, but with f_neg≈0.5 it does not translate to a clean
   avoidance the sign could ride.

## Verdict on the multiverse gate

Wilson set the criterion: *"If the alignment story holds at q=3, the multiverse spend is justified. If the alignment
isn't there, we'll have learned that before paying."* **The alignment/avoidance story does NOT hold at q=3** (f_neg≈0.5,
isotropic directions). So the matrix-cone-avoidance route does not survive as a mechanism to generalize — exactly the
finding his gate was built to catch before paying for base-q carry machinery. **The multiverse spend is not justified
on this basis.** (Alternative reading, not a recommendation: the isotropy of the direction distribution might itself
be a q=3 *degeneracy* — at gapped q≥5 the directions could be anisotropic and the cone real. That would re-open the
multiverse, but as a fresh conjecture, not off this result.)

**What survives, sharpened:** the sign continues to live in `sign(tr S_c)` = the signed coherence = CARRYLEMMA's
`T1+T2`; the reframes keep confirming it is exactly, and only, that. The genuine gains are (a) the well-conditioned
variable — `cos_θ_true` stable/rising ~0.04, angle ~88°, D stable below √2; (b) the κ²↔U0+U3 "one event" lead at
r≈11–12.

**Not at stake:** CHANNEL_ID (`d1=A_r(1)/S_r`) + CARRYLEMMA identities, R1–R30, R80–R82, all Thread-3 probes.
`build_nu(0.5,16)` (360s). `scratchpad/carrymat.json` saved.
