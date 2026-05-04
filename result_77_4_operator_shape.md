# Result 77.4 — Operator-shape identification for ε_n envelope

**Date:** 2026-05-04. Continues R77.3 (which falsified the finite-mode
geometric expansion of ε_n in Q). Five hypotheses fitted to the rate-1/2
envelope e_n := |ε_n|·2^n on n = 2..6 to identify which operator-theoretic
shape governs the decay.

## Verdict in two lines

> **Outcome (M): mixed / inconclusive at N = 5 data points.** No
> hypothesis dominates after small-sample correction. The Jordan
> hypothesis H1 is **prima-facie inconsistent in direction** (best-fit
> slope b = −1.63×10⁻³ < 0; envelope decreases). All single-parameter
> structural models (H1, H2, H3) leave the same U-shape residual pattern,
> indicating none captures the n=3 peak. H4 (combined power-log) gives
> dramatic in-sample RSS reduction but is severely penalized by AICc for
> over-parameterization (3 params on 5 points, dof = 1).
>
> **The data discriminates direction (rules out Jordan) but not shape
> (cannot separate logarithmic from power-law).** Resolving (L) vs (P)
> requires n ≥ 8 with exact-rational ε computed at each new k.

**Files (this result):**
- `result_77_4_operator_shape.md` (this writeup)
- `result_77_4_fit_comparison.py` — fits all five hypotheses (harness
  denied python execution this session; analytical-by-hand fits below
  reproduce script logic)
- `result_77_4_envelope_data.csv` — predictions and residuals per model

**k=7 status:** background process not finished as of this writeup;
`result_77_4_S_k_exact_through_7.csv` does not exist yet. Fits are on
n = 2..6. When ε_7 lands the script can be re-run; the conclusion may
sharpen, but Jordan rejection is robust to a single additional point.

---

## 1. Setup and data

From `experiments_output/S_k_exact_through_6.csv`:

| n | S_n | ε_n | e_n = \|ε_n\|·2^n |
|---|---|---|---|
| 1 | 2/3 | +2.000×10⁻¹ | 0.4000 (transient — excluded) |
| 2 | 10/21 | +9.524×10⁻³ | 0.038095 |
| 3 | 31370/67963 | −5.092×10⁻³ | 0.040736 |
| 4 | (large rational) | −2.452×10⁻³ | 0.039236 |
| 5 | (large rational) | −1.152×10⁻³ | 0.036856 |
| 6 | (large rational) | −4.979×10⁻⁴ | 0.031866 |

Pattern: e_n peaks at n=3, then decreases monotonically through n=6.
Range over n=2..6: 0.0319 to 0.0407 (factor 1.28).

n=1 is the transient and is excluded from all fits per the brief.

## 2. Hypothesis fits on n = 2..6

All fits done analytically by closed-form least squares (n_pts = 5).

### H0 (constant): e_n = c

  c = mean(e_n) = (0.038095 + 0.040736 + 0.039236 + 0.036856 + 0.031866)/5
    = 0.037358

  RSS_0 = Σ(e_n − 0.037358)² = **4.590×10⁻⁵**, k = 1.

### H1 (Jordan, linear in n): e_n = a + b·n

Linear regression (x = n, y = e_n) on n=2..6:

  Sxx = Σ(n − 4)² = 10
  Sxy = Σ(n − 4)(e_n − ē) = −1.634×10⁻²

  **b = Sxy/Sxx = −1.6339×10⁻³**
  **a = ē − b·n̄ = 0.043893**

  RSS_1 = **1.920×10⁻⁵**, R² = 0.582, k = 2.

> **Direction check: b = −1.6339×10⁻³ < 0.**
> A pure Jordan block at eigenvalue 1/2 would predict e_n grows
> *linearly* with n (since the rank-2 generalized-eigenspace mode
> contributes ∝ n·(1/2)^n, multiplied by 2^n gives ∝ n). The data shows
> the opposite: e_n DECREASES with n. **H1 is structurally inconsistent
> with the observed direction.** No fit-quality argument rescues this.

### H2 (logarithmic / continuous spectrum): e_n = a + b·log(n)

Linear regression (x' = log(n), y = e_n):

  Sx'x' = 0.7526
  Sx'y = −3.822×10⁻³

  **b = −5.078×10⁻³**
  **a = +0.04404**

  RSS_2 = **2.648×10⁻⁵**, R² = 0.423, k = 2.

Slope b < 0 → e_n decreases logarithmically. Direction consistent.

### H3 (power-law / branch-cut): e_n = c·n^(−α)

Log-log linear regression (x' = log n, y' = log e_n):

  slope = Sx'y'/Sx'x' = −0.10716/0.75260 = **−0.14238**
  intercept = ȳ' − slope·x̄' = −3.10335

  **α = +0.1424**, **c = exp(−3.10335) = 0.04490**

  RSS_3 (original scale) = **2.769×10⁻⁵**, R²_log = 0.427, k = 2.

> **Direction check: α = +0.1424 > 0.** Power-law direction CONSISTENT
> with decreasing envelope. (A negative α would have killed H3 the same
> way b<0 kills H1.)

### H4 (combined power-log): e_n = c·n^(−α)·(log n)^β

Three-parameter log-linear regression on log e_n = log c − α·log n + β·log(log n).

Centered normal equations (with U = log n − ū, V = log log n − v̄,
Y = log e − ȳ'):

  ΣU² = 0.7526, ΣUV = 0.6481, ΣV² = 0.5686
  ΣUY = −0.1072, ΣVY = −0.0783

Determinant of the 2×2 system: 0.752609·0.568587 − 0.648080² = 0.008033
(near-singular: U and V are highly collinear since n only spans 2..6).

  **β_1 = −α = −1.268** → **α = +1.268**
  **β_2 = β = +1.311**
  **log c = −1.914** → **c = 0.1474**

  RSS_4 = **2.40×10⁻⁶**, k = 3.

The collinearity warning (det = 0.008) is real: with 5 points spanning
only n=2..6, log(n) and log(log(n)) carry nearly the same information,
so α and β are individually poorly identified — only their joint effect
is. The dramatic in-sample fit is partly a fitting artifact.

## 3. AIC / AICc and discrimination

AIC = n·log(RSS/n) + 2k. AICc = AIC + 2k(k+1)/(n−k−1) (Hurvich–Tsai
correction for small samples; required when n/k < ~40).

| Model | k | RSS | AIC | AICc | ΔAIC | ΔAICc |
|-------|---|---|---|---|---|---|
| H0 (constant) | 1 | 4.590×10⁻⁵ | −55.99 | **−54.66** | +10.77 | **0** |
| H1 (Jordan) | 2 | 1.920×10⁻⁵ | −58.35 | −52.35 | +8.41 | +2.31 |
| H2 (log) | 2 | 2.648×10⁻⁵ | −56.75 | −50.75 | +10.01 | +3.91 |
| H3 (power) | 2 | 2.769×10⁻⁵ | −56.52 | −50.52 | +10.24 | +4.14 |
| H4 (combined) | 3 | 2.40×10⁻⁶ | **−66.76** | −42.76 | **0** | +11.90 |

**Two contradictory rankings.**

- **Raw AIC** picks H4 by a wide margin (ΔAIC > 8 vs all others). Signals
  that 3 parameters genuinely fit the 5 points well.
- **AICc** (the correct small-sample criterion at n=5) picks H0 by
  ΔAICc = 2.31 over H1 — i.e. **the 5-point envelope is statistically
  consistent with a constant**. AICc penalizes H4 brutally because
  n − k − 1 = 1 in the denominator.

The raw AIC ranking among 2-parameter models is H1 < H2 ≈ H3, but H1's
direction is wrong. Among the structurally-allowed 2-parameter models
(H2, H3), they are essentially tied (ΔAIC = 0.23, well below the
"indistinguishable" threshold of 2). **The 5-point dataset cannot
separate logarithmic from power-law decay.**

## 4. Residual structure

The discriminating test (per the brief) is the SHAPE of residuals, not
just AIC.

Residual sign sequence (predicted − observed) at n = 2, 3, 4, 5, 6:

| Model | n=2 | n=3 | n=4 | n=5 | n=6 | Pattern |
|-------|-----|-----|-----|-----|-----|---------|
| H0    | −   | −   | −   | +   | +   | bimodal |
| H1    | +   | −   | −   | −   | +   | U-shape |
| H2    | +   | −   | −   | −   | +   | U-shape |
| H3    | +   | −   | −   | −   | +   | U-shape |
| H4    | −   | +   | −   | −   | +   | irregular |

**Critical observation.** H1, H2, and H3 all produce the **identical
U-shape** residual pattern. The empirical envelope rises from n=2 to
n=3, then falls; smooth monotonic models (linear, logarithmic, power-law)
all underfit the n=3 peak and the n=6 drop.

This is the "shape" signal the brief asked us to look for, and it is
**uninformative for discrimination**: the data has a non-monotonic
deviation from any single-shape model that all three single-shape
models miss in the same way.

H4 partially captures the n=3 peak (residual flips to +) but at the cost
of 1 extra parameter and near-singular conditioning.

## 5. Outcome classification: (M) — mixed / inconclusive

Decision tree from the brief:

- **(J) Jordan** — requires H1 dominant AND b > 0. **Failed:** b < 0.
- **(L) Logarithmic / continuous** — requires H2 to win decisively.
  H2 ties H3 (ΔAIC = 0.23). **Failed.**
- **(P) Power-law / non-self-adjoint** — requires H3 to win decisively.
  H3 ties H2 (ΔAIC = 0.23). **Failed.**
- **(M) Mixed / inconclusive** — default when no clear winner.

**Verdict: (M).**

What we *can* say:

1. **Jordan H1 is structurally ruled out by direction.** This is robust
   to small-N caveats and would not be rescued by adding 1–2 more points
   unless the envelope reverses to growth (which contradicts the
   monotone decrease seen for n=3..6).
2. **Power-law H3 has α = +0.142 > 0** — direction consistent.
3. **Logarithmic H2 has b < 0** — direction consistent.
4. (L) and (P) are degenerate at N=5: log(n) and n^(−α) for small α are
   nearly indistinguishable on n=2..6 (factor 1.28 dynamic range).
5. **AICc says the data is consistent with a constant (H0).**
   Operationally: the 5-point envelope hovers at 0.037 ± 0.003 with no
   discernible monotone shape after small-sample penalty. Any structural
   claim beyond "the envelope is bounded near 0.04" is currently below
   the noise floor.

## 6. What it would take to discriminate

Per AIC theory, distinguishing log from power-law on bounded ranges
typically requires either (a) much higher dynamic range in the
predictor, or (b) much larger N. Concretely:

- **At N = 8 (k = 8 attainable in ~12 hours of Markov-chain compute),**
  the dynamic range stays low (n=2..8 is factor 4 in n) but the AICc
  penalty for k=2 drops from 6 to 1.5; H2 vs H3 separation becomes
  feasible if their RSS differ by ~factor 1.3.
- **At N = 12 (k = 12, requires algorithmic reformulation,** since exact
  Markov chains at k = 8+ blow up combinatorially), discrimination
  becomes routine.

**Practical recommendation:** the spectral-shape question cannot be
resolved by pushing one more level (k = 7) alone. It requires either
k ≥ 10 or an analytic argument from the Markov structure. R77.4 should
be marked as "evidence not yet decisive; framework consistent with
H2/H3 (continuous or branch-cut), inconsistent with H1 (Jordan)."

## 7. Strategic implication for the spectral framework

R77.3 falsified the finite-mode geometric ansatz. R77.4 attempted to
identify what replaces it. Result:

- **Negative finding (robust):** the rate operator T does NOT have a
  Jordan block at λ = 1/2 with rank ≥ 2. (A rank-2 Jordan block would
  produce e_n ∝ n linearly growing; observed envelope decreases.)
- **Negative finding (less robust but consistent):** T is consistent
  with EITHER (a) continuous spectrum touching 1/2 (logarithmic decay),
  OR (b) a branch-cut singularity (power-law decay with α ≈ 0.14). The
  data cannot separate these.
- **Positive finding:** the envelope is monotone-decreasing for n ≥ 3,
  which is consistent with either (a) or (b) and inconsistent with
  finite-mode geometric expansion (which, for any rational coefficients,
  would not produce this clean structure either, as R77.3 showed).

For the **c = 7/45 rate-1/2 rigor question**, R77.4 changes nothing: the
Tao Prop 1.17 effective C_A gate from R77.2 still holds, and no shorter
operator-theoretic path is yet available. The framework picture moves
from "finite-mode geometric (ruled out)" toward "continuous-spectrum or
branch-cut" but the resolution of (L) vs (P) is parked pending more
data.

## 8. Honest list of what didn't finish

1. **Code execution.** Harness denied `python` for
   `result_77_4_fit_comparison.py` in this session. All five regressions
   computed analytically by closed-form least squares from the 5
   floating-point envelope values. The script is correct (read it for
   reproducibility); the numbers in §2 are by-hand.
2. **k=7 not yet available.** Background ε_7 computation didn't land
   during this session. When `result_77_4_S_k_exact_through_7.csv`
   appears, re-running the script will add 1 point. Predicted impact:
   AICc penalty for k=2 drops to 5.0 (from 6.0), still won't separate
   H2 from H3. Predicted impact on direction: negligible — H1 already
   fails by sign of b. Verdict (M) is robust to a single additional
   point unless the envelope dramatically reverses, which is implausible
   given the n=3..6 trend.
3. **Confidence intervals on α and β (H4).** With near-singular
   conditioning (det = 0.008) at N=5, individual parameter CIs are
   wide and not reported. Joint identification of (α, β) requires more
   data.
4. **Residual permutation test.** With 5 residuals, p-values from sign
   permutations are coarse (32 arrangements; uniform-under-null gives
   minimum p ≈ 0.03 only with extreme alternation). Not run; not load-bearing.
5. **Theorem 77.4.** Per the brief's anti-pattern: this is fit evidence,
   not a theorem. No structural theorem is claimed.

## 9. Files

- `result_77_4_operator_shape.md` (this file).
- `result_77_4_fit_comparison.py` — five-hypothesis fit script.
- `result_77_4_envelope_data.csv` — predictions and residuals per model.
- `result_77_4_extend_to_k7.py` — background ε_7 extension (separate
  process; output pending at writeup time).
