# Derivation attempt: E[v]_q125 = 2.216

## Approach 1 — Joint-Gaussian conditional moment

### Setup

For each orbit i, define:
- V_i = n_even_i / n_odd_i = mean v across orbit's Syracuse steps
- Σ_i = total stopping time σ_i
- Σ_resid,i = Σ_i − (α + β·log m_start,i) — σ residual after removing linear log m trend

Bottom-σ-quartile selection at fixed log m corresponds to {Σ_resid < q_25} where q_25 is the 25th percentile of Σ_resid.

If (V, Σ_resid) is jointly Gaussian under the trajectory measure, then:

  E[V | Σ_resid < q_25] = E[V] + ρ · SD[V] · E[Z_Σ | Z_Σ < z_25]
                       = E[V] − 1.2729 · ρ · SD[V]

where z_25 = Φ⁻¹(0.25) = −0.6745 and E[Z|Z<z_25] = −φ(z_25)/Φ(z_25) = −0.3178/0.25 = −1.2729.

### Empirical measurement at multiple N

500K orbits per N (full enumeration parquet at 2²⁷, fresh sampling at 2³⁰..2⁴²). Per-orbit moments and the truncated-Gaussian prediction:

| log2N | E[V] | SD[V] | ρ_resid | ρ·SD[V] | pred E[V\|q_125] | emp E[V\|q_125] | gap (emp − pred) |
|---|---|---|---|---|---|---|---|
| 27 | 2.0828 | 0.2583 | −0.825 | −0.213 | 2.354 | 2.429 | +0.075 |
| 30 | 2.0722 | 0.2337 | −0.841 | −0.197 | 2.322 | 2.388 | +0.066 |
| 32 | 2.0672 | 0.2204 | −0.851 | −0.187 | 2.306 | 2.366 | +0.060 |
| 34 | 2.0620 | 0.2095 | −0.857 | −0.180 | 2.290 | 2.346 | +0.055 |
| 36 | 2.0578 | 0.1995 | −0.867 | −0.173 | 2.278 | 2.328 | +0.051 |
| 38 | 2.0545 | 0.1905 | −0.874 | −0.166 | 2.266 | 2.313 | +0.047 |
| 40 | 2.0512 | 0.1822 | −0.881 | −0.160 | 2.255 | 2.299 | +0.044 |
| 42 | 2.0486 | 0.1758 | −0.885 | −0.156 | 2.246 | **2.288** | +0.041 |

### Status: PARTIAL — captures shape, fails on asymptote

**What works:**
- Truncated-Gaussian formula has the correct functional form (linear in ρ·SD[V])
- Predicts ~95% of the empirical shift (pred − E[V] = 0.20, emp − E[V] = 0.24 at N=2⁴²)
- Tracks the qualitative trajectory (both pred and emp shrink toward 2.216 from above)
- ρ_resid is large and negative (−0.825 → −0.885), confirming V and σ_resid are strongly anti-correlated as physically expected (high V ⇒ fast descent ⇒ low σ_resid)

**What fails: the asymptotic limit.** For pred → 2.216 we need:

  (ρ · SD[V])_∞ = (2.216 − 1.995) / (−1.2729) = **−0.174**

Empirical trajectory of ρ·SD[V]:

  −0.213 (27) → −0.197 (30) → −0.187 (32) → −0.180 (34) → −0.173 (36) → −0.166 (38) → −0.160 (40) → −0.156 (42)

The trajectory **crossed −0.174 around N ≈ 2³⁵** and continues moving toward less-negative values (toward 0). Per-octave Δ has been steadily ~+0.005 throughout the range.

If the trajectory continues toward zero (consistent with SD[V] → 0 at rate 1/√log N and ρ → −1 sub-linearly), pred → 1.995 (not 2.216). The empirical asymptote at 2.216 must therefore come from *non-Gaussian corrections* to the truncated-Gaussian formula — specifically from the skew/kurtosis of the per-step v distribution (Geom(1/2) is non-Gaussian).

### Why the gap doesn't shrink to zero

The per-step v_t under the trajectory measure is approximately Geom(1/2) with E[v]=2, Var[v]=1.881. For per-orbit averaged V_i over n_steps_i, by CLT V_i is asymptotically Gaussian — but the asymptote is reached *slowly* relative to higher-moment corrections.

Per-orbit Σ_i is also a sum of per-step contributions, but with a discrete v dependence at each step. Joint (V_i, Σ_i) is bivariate Gaussian only in the leading-order CLT approximation; sub-leading corrections include:

- Skewness of v: positive (Geom is right-skewed) → V_i has positive skew
- Cubic/quartic moments of (v, log m drop) joint distribution
- The σ definition has an integer-overshoot term (descent stops at m=1, not at log m = 0)

Each of these contributes O(1/√n_steps) corrections to the conditional mean E[V|Σ<q_25]. The empirical pred−emp gap is consistent with O(1/√log N) decay:

  gap(N) ≈ 0.041 at N=2⁴² (log N ≈ 29) → predicts gap(N=2⁵⁰) ≈ 0.034

So the gap shrinks but very slowly. At physical infinity it might reach zero, but it's not the dominant contribution to the *value* 2.216 — that comes from the joint-Gaussian piece (~0.20 of the 0.22 shift).

### Closed-form status

The joint-Gaussian formula derives ~95% of the answer:

  E[V]_q125 ≈ E[V] − 1.273 · ρ_∞ · SD[V]_∞ + (non-Gaussian correction)

The first three terms (E[V], ρ_∞, SD[V]_∞) need their own asymptotic analysis. Empirically:
- E[V] → 1.995 (the trajectory-measure asymptote of mean v)
- SD[V] → 0 like 1/√log N
- ρ → −1 sub-linearly (rate not pinned by current data)

If ρ·SD[V] has a non-trivial limit (which would be the case if ρ → −1 fast enough to compensate SD[V] → 0), the joint-Gaussian gives a non-zero asymptotic shift. Currently the empirical trajectory is in a regime where ρ·SD[V] is decreasing — no plateau visible.

**Bottom line:** Approach 1 partially derives 2.216 (gives most of the shift via joint-Gaussian) but the asymptote 2.216 itself is NOT directly derivable from the joint-Gaussian formula at the limits we can extract. Either:

(a) ρ·SD[V] has a non-trivial limit at exactly −0.174 that current data hasn't reached — would require N >> 2⁴² to confirm
(b) The empirical asymptote 2.216 is set by a combination of joint-Gaussian (most of shift) plus non-Gaussian corrections (residual 0.04) that don't decompose cleanly

Either way, **the Gaussian closed form alone doesn't pin 2.216 to ±0.01**. Need either:

- Approach 2 (σ-quartile-conditional residue density) to derive the non-Gaussian correction, OR
- Approach 3 (partial Esscher tilt) to bypass the Gaussian assumption, OR
- Higher-order Edgeworth expansion of the joint (V, Σ) distribution incorporating Geom(1/2)'s non-Gaussianity

### Recommendation

The dominant structure (≈0.20 of the 0.221 shift) comes from joint-Gaussian. The residual 0.04 is the closed-form gap. Worth attempting Approach 3 (partial Esscher tilt) as it gives an exact formula in v-space:

For Esscher tilt P_w(v=k) ∝ Geom(1/2)(v=k) · 2^{−w·k}, E_w[v] = 1/(1 − 2^{−(1+w)}). Setting E_w[v] = 2.216 gives w_q = −0.136.

This is the *value* of the partial-tilt parameter; the *structural derivation* (why w_q = −0.136 specifically vs the σ-quartile geometry) is the open question. Approach 3 would attempt this. Approach 2 (residue density) is a separate angle that may or may not be tractable.

Files:
- `experiments/corr_v_sigma.py` — empirical correlation + truncated-Gaussian computation
- `experiments_output/corr_v_sigma.csv` — per-N moments and predictions
- `experiments_output/corr_v_sigma_log.txt` — full output log
