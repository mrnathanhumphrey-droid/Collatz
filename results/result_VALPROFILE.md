# RESULT — VALPROFILE (Wilson's norm route): denominator of N is closed-form; numerator rogue; geo-mean ≠ aggregate (2026-07-26)

**Probe:** `probe_valprofile.py`. Wilson's 3-step norm program on `N = Π_{3∤a} π̂_k(a) = N_{ℚ(ζ_{3^k})/ℚ}(π̂(1))`:
(1) v_p(N) at every p, (2) product formula → |N|_∞, (3) geometric mean vs typical √k·3^{−k/2}. Exact k=2..5 (k=6 v₃
running).

## The structural WIN — the denominator of N is fully closed-form (k=2..5)
| k | φ | v₂ | v₃ (pred 1−φ) | v_p for p\|2^M−1, p≠3 |
|---|---|----|----|----|
| 2 | 6 | 0 | −5 ✓ | v₇=−6 |
| 3 | 18 | 0 | −17 ✓ | v₇=v₁₉=v₇₃=−18 |
| 4 | 54 | 0 | −53 ✓ | v₇=v₁₉=v₇₃=v₈₇₂₁₁=v₂₆₂₆₅₇=−54 |
| 5 | 162 | 0 | −161 ✓ | all 11 of {7,19,73,163,2593,71119,87211,135433,262657,97685839,272010961}=−162 |

**`v₂(N)=0`, `v₃(N)=−(φ−1)`, and `v_p(N)=−φ` for EVERY prime `p | 2^M−1` with `p≠3`** (M=ord₂(3^k)=2·3^{k−1}). Since
`2^M−1 = 3^k·(non-3 part)` (LTE: v₃(2^M−1)=k), the **denominator of N is exactly**
`den(N) = 3^{φ−1} · [(2^M−1)/3^k]^φ`.
Every "expected" valuation is pinned in closed form — the p-adic pole structure of N is completely determined.
(λ-adic: `v_λ(π̂(1)) = 1−φ`, i.e. `3·π̂(1)` is an algebraic integer divisible by λ=1−ζ exactly once.)

## The honest NEGATIVES
1. **The numerator is ROGUE and proliferating.** After stripping {2,3,div(2^M−1)} the residual is a pure numerator of
   size **4 → 65 → 740 → 7369 digits** (k=2..5), large unstructured primes (k=2: 19·487; k=3 adds a 29-digit and two
   16-digit primes; …). So N is **not** {2,3,div(2^M−1)}-supported, and **|N|_∞ has no closed arithmetic form** — the
   numerator carries it. The Dwork/Borel shelf pins the p-adic radii (denominator) but the archimedean size stays wild.
2. **The geometric mean does NOT track the aggregate typical** (item-3 not confirmed):
   | k | ln(geo mean)=ln\|N\|/φ | ln(√k·3^{−k/2}) | gap |
   |---|-----|-----|-----|
   | 2 | −1.339 | −0.752 | −0.587 |
   | 3 | −1.932 | −1.099 | −0.834 |
   | 4 | −2.565 | −1.504 | −1.061 |
   | 5 | −3.163 | −1.942 | −1.221 |
   The geo mean sits **below** and the gap **grows ~linearly (−0.21/k)**. This is Jensen (geometric mean ≤ quadratic
   mean = typical, always), so the growing gap = **the spectral log-variance widening ~linearly in k**. Real
   distributional content — but the geometric mean is dominated by the spectral **low tail**, not the bulk/aggregate,
   so it does **not** hand over S_∞.

## Verdict — the norm route gives a valuation theorem, not S_∞
The clean, bankable result is the **closed-form denominator** `den(N) = 3^{φ−1}·[(2^M−1)/3^k]^φ` (v₂=0, v₃=1−φ, v_p=−φ
on all Mersenne primes) — a genuine arithmetic theorem-in-waiting about the norm's poles, and a natural Dwork-type
object (p-adic radii controlled). But the **archimedean side is not closed**: the numerator is rogue and the geometric
mean is a Jensen-below-typical third moment governed by the widening spectral log-variance, not the aggregate. So the
norm route **does not reach S_∞** as item-3 hoped; its deliverable is the pole/valuation structure plus one honest
distributional fact (log-variance grows linearly in k). Not at stake: NORMCHECK (v₃ law), RECENTER, CHANNEL_ID, v₃
HIERARCHY, R1–R30. Cheap (12.7s, k≤5; k=6 v₃ heavy).
