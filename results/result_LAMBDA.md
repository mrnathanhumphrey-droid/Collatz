# RESULT — PROBE LAMBDA: the aligned table; the saving δ is real but SHRINKING toward the ceiling, not a constant (2026-07-26)

**Probe:** `probe_lambda.py`. Wilson's ceiling reframe: at i=k, 3^k≡0 ⟹ π̂=1 ⟹ G(k,j)=2^j; phases≤1 preserve it down
the tree ⟹ **|G(0,m)| ≤ 2^m unconditionally** ⟹ sup|π̂|=2^{−m}|G|≤1 (trivial). So the decay IS the saving of G below
its ceiling — one number, not two competing growths. Deliverable: the aligned one-row-per-k table with explicit
conventions, λ fitted per-unit-m against the ceiling (=2), and the convergence read.

**Conventions (explicit — they crossed earlier this session):** `sup = |π̂(ξ*)|` LINEAR; `R66 = sup²`;
`|G(0,m)| = 2^m·sup` (ceiling 2^m). `|G|/2^m ≡ sup` (identity, checked).

## The aligned table (one row per k, k=3..15)
| k | ξ* | m | dm | sup=\|π̂\| | R66=\|π̂\|² | \|G(0,m)\| | Srate | λ_m | δ=1−λ_m/2 |
|---|-----|---|----|----------|-----------|-----------|-------|-----|-----|
| 3 | 8 | 3 |  | 0.252237 | 0.063623 | 2.018 |  |  |  |
| 4 | 65 | 4 | +1 | 0.176999 | 0.031329 | 2.832 | 0.7017 | 1.4034 | 0.298 |
| 5 | 32 | 5 | +1 | 0.129274 | 0.016712 | 4.137 | 0.7304 | 1.4607 | 0.270 |
| 6 | 64 | 6 | +1 | 0.096106 | 0.009236 | 6.151 | 0.7434 | 1.4869 | 0.257 |
| 7 | 256 | 8 | +2 | 0.075870 | 0.005756 | 19.42 | 0.7894 | 1.7770 | 0.111 |
| 8 | 6049 | 9 | +1 | 0.060891 | 0.003708 | 31.18 | 0.8026 | 1.6051 | 0.197 |
| 9 | 1024 | 10 | +1 | 0.048026 | 0.002307 | 49.18 | 0.7887 | 1.5775 | 0.211 |
| 10 | 54953 | 12 | +2 | 0.038278 | 0.001465 | 156.8 | 0.7970 | 1.7855 | 0.107 |
| 11 | 8192 | 13 | +1 | 0.031944 | 0.001020 | 261.7 | 0.8345 | 1.6690 | 0.166 |
| 12 | 515057 | 14 | +1 | 0.026458 | 0.000700 | 433.5 | 0.8283 | 1.6565 | 0.172 |
| 13 | 65536 | 16 | +2 | 0.022052 | 0.000486 | 1445.2 | 0.8335 | 1.8259 | 0.087 |
| 14 | 131072 | 17 | +1 | 0.019128 | 0.000366 | 2507.1 | 0.8674 | 1.7348 | 0.133 |
| 15 | 262144 | 18 | +1 | 0.016284 | 0.000265 | 4268.9 | 0.8513 | 1.7027 | 0.149 |

**ξ* is a PURE POWER OF 2 at every k=3..15** (no exceptions). **m(k)−k = ⌊(k−4)/3⌋** exactly (m−k = 0,0,0,1,1,1,2,2,2,3,3,3
for k=4..15); the +2 jumps land at k=7,10,13 ⟹ **dm/dk → 4/3**, ceiling `2^{4/3} ≈ 2.52` per k confirmed. The earlier
"λ≈2.0" was per-k dm-averaged and conflated; the clean per-unit-m growth is **λ_m ≈ 1.7**.

## The saving δ is real but SHRINKING — not a constant (the finding)
Per-unit-m λ_m is contaminated at the dm=+2 steps (geometric-mean across a skipped m and a modulus change → reads
high: 1.78, 1.79, 1.83 at k=7,10,13). The clean comparison is the **3-block** (after k=6 every block has Δm=4, Δk=3,
apples-to-apples):
| block | k | \|G\| ratio | λ_m=(ratio)^{1/4} | δ=1−λ_m/2 | sup ratio | vs ceiling 2^{−4}=0.0625 |
|-------|---|-----------|-------------------|-----------|-----------|--------------------------|
| 6→9 | 3→9 | 7.994 | 1.681 | **0.160** | 0.4996 | beats ceiling 8.0× |
| 9→12 | " | 8.815 | 1.723 | **0.139** | 0.5510 | 8.8× |
| 12→15 | " | 9.847 | 1.766 | **0.117** | 0.6156 | 9.85× |

**δ drops ~0.022 per block, monotone.** The pure-dm=1 steps corroborate (δ per block ≈ 0.20 → 0.17 → 0.14). So:
- **λ_m is climbing toward the ceiling 2** (1.45 → 1.58 → 1.66 → 1.72 across blocks); the per-unit-m sup-saving
  `λ_m/2` climbs 0.72 → 0.86 toward 1.
- **Srate climbs 0.70 → 0.85** (R66=sup² rate climbs 0.49 → 0.72). ⚠️ This **overshoots MAXMODE2's "saturation ~0.655"**
  — the max-mode rate did NOT saturate at 0.655; it keeps climbing past it. Correction to that read.
- **δ is k-dependent and shrinking, NOT a fixed number.** 3 clean blocks cannot distinguish δ→δ_∞>0 (geometric decay
  of sup survives) from δ→0 (sub-geometric decay). Linear extrapolation of δ hits 0 near k≈27; a log-linear approach
  would floor at a small positive constant. **Unresolved with k≤15** (k=16 is 3^16≈43M states, heavy).

## Verdict — the target sharpens, and it is UNIFORMITY, not positivity
The closed form is exact (GRECURSION), the ceiling is `2^m` (proved), and the decay is exactly the ⟨2⟩-orbit sum
`Σ_b e(2^b/3^{k−i})` beating its trivial bound, compounded down the tree. That is the Heilbronn/BGK power-saving object
— confirmed. **But the measured saving δ is not a constant; it shrinks with k** (0.30 → 0.12 over k=4..15). So the
precise open question is NOT "is there a power saving δ>0" (there is, at every finite k) but **"is the ⟨2⟩-orbit-sum
power saving UNIFORM in the modulus 3^n (δ bounded below), or does it degrade as n→∞?"** The tree depth grows with k
and the top-level modulus 3^k is the largest, so the downdrift may be a depth-accumulation / large-modulus effect —
Wilson's pen. For the channel budget this matters: a uniform δ gives clean geometric decay of sup|π̂| (hence of R66),
a degrading δ gives sub-geometric decay and changes the S_∞ bookkeeping.

**Hank's sharpened target:** *effective, UNIFORM-in-n power-saving bounds for `Σ_b e(2^b/3^n)` over the ⟨2⟩-orbit,
prime-power modulus 3^n — is the saving exponent bounded below independent of n, or does it →0?* (The distinction
between a fixed and a degrading power saving is the whole question; magnitude-only decay results are silent on it.)

**Not at stake:** GRECURSION (closed form), SINGLEREC, BRIDGE2, CONTRACTION, MEAN1, HIERARCHY, CHANNEL_ID, R1–R30.
⚠️ MAXMODE2's "rate saturates ~0.655" is CORRECTED (climbs to ~0.72 by k=15). Cost: 127s (k≤15; fwd_hat grows ~3^k).
