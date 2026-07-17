# Stopping criteria battery — multiple σ measurements at q=3, N=2³²

**Date:** 2026-05-03. Numerical: `stopping_criteria_battery.py`. CSV: `stopping_criteria_battery.csv`. 200K orbits, 2.5s walk-time, all 200K reached 1.

## 1. Per-criterion firing rates and statistics

| criterion | fired | rate | mean σ | SD | p50 | p95 |
|---|---|---|---|---|---|---|
| σ_to_1 | 200,000 | 100.000% | 224.53 | 74.31 | 218 | 357 |
| σ_prime | 200,000 | 100.000% | 17.39 | 16.60 | 13 | 51 |
| σ_pow2_L10 (≥ 2¹⁰) | 7,932 | **3.966%** | 172.30 | 62.23 | 162 | 288 |
| σ_logdrop_1 | 200,000 | 100.000% | 20.19 | 21.21 | 12 | 61 |
| σ_logdrop_2 | 200,000 | 100.000% | 30.27 | 26.82 | 21 | 84 |
| σ_logdrop_4 | 200,000 | 100.000% | 51.12 | 35.65 | 42 | 120 |
| σ_logdrop_8 | 200,000 | 100.000% | 92.63 | 48.27 | 82 | 185 |
| **σ_res_0mod3** | 66,780 | **33.390%** | **0.00** | 0.00 | 0 | 0 |
| σ_res_1mod8 | 200,000 | 100.000% | 8.99 | 10.69 | 6 | 30 |
| σ_res_5mod8 | 200,000 | 100.000% | 6.98 | 8.10 | 4 | 23 |
| σ_res_0mod5 | 199,998 | 99.999% | 11.86 | 15.50 | 6 | 43 |

## 2. STRUCTURAL FINDING: σ_res_0mod3 reveals a Collatz invariant

**σ_res_0mod3 fires at exactly σ=0 for 33.4% of starts and NEVER fires for the other 66.6%.**

Mechanism: Collatz preserves the set {n : n ≢ 0 mod 3}. Proof:
- If x odd and x ≢ 0 mod 3: T(x) = 3x+1 ≡ 1 mod 3 (always, since 3x ≡ 0 mod 3)
- If x even: T(x) = x/2. Mod 3 arithmetic: x/2 ≡ 2x mod 3 (since inv2 mod 3 = 2). x ≡ 1 → 2; x ≡ 2 → 1. So x/2 ≢ 0 mod 3 whenever x ≢ 0 mod 3 and x even.
- Therefore once an orbit leaves (or never enters) "≡ 0 mod 3", it never visits it again.

**Implication:** the "not divisible by 3" subset is INVARIANT under Collatz iteration (modulo finite transients from initial multiples-of-3). This is a clean structural fact that the standard σ_to_1 framework doesn't surface.

## 3. Cross-criterion correlation matrix (on 2,681 orbits where ALL fired)

|  | σ_to_1 | σ_prime | σ_pow2 | σ_ld1 | σ_ld2 | σ_ld4 | σ_ld8 | σ_1m8 | σ_5m8 | σ_0m5 |
|---|---|---|---|---|---|---|---|---|---|---|
| σ_to_1 | **1.000** | −0.022 | **+0.9999** | +0.344 | +0.412 | +0.551 | +0.761 | −0.080 | +0.241 | −0.029 |
| σ_prime | −0.022 | 1.000 | −0.022 | −0.051 | −0.033 | −0.025 | −0.022 | −0.006 | −0.031 | −0.075 |
| σ_pow2_L10 | +0.9999 | −0.022 | 1.000 | +0.344 | +0.411 | +0.551 | +0.760 | −0.080 | +0.240 | −0.029 |
| σ_ld1 | +0.344 | −0.051 | +0.344 | 1.000 | +0.821 | +0.608 | +0.458 | −0.101 | **+0.633** | +0.032 |
| σ_ld8 | +0.761 | −0.022 | +0.760 | +0.458 | +0.544 | +0.734 | 1.000 | −0.097 | +0.317 | −0.021 |

(σ_res_0mod3 omitted — fires at σ=0 with no variance; correlations undefined.)

## 4. Per-σ_to_1 band ratios (σ_criterion / σ_to_1)

| criterion | q1 (low σ) | q2 | q3 | q4 (hi σ) | q5 (95+) |
|---|---|---|---|---|---|
| σ_prime | **0.133** | 0.090 | 0.071 | 0.054 | **0.043** |
| σ_pow2_L10 | 0.921 | 0.947 | 0.958 | 0.968 | **0.975** |
| σ_logdrop_1 | 0.102 | 0.092 | 0.085 | 0.087 | 0.094 |
| σ_logdrop_2 | 0.149 | 0.137 | 0.127 | 0.132 | 0.145 |
| σ_logdrop_4 | 0.245 | 0.230 | 0.215 | 0.227 | 0.250 |
| σ_logdrop_8 | 0.440 | 0.415 | 0.390 | 0.413 | 0.456 |
| σ_res_1mod8 | 0.076 | 0.047 | 0.036 | 0.026 | 0.020 |
| σ_res_5mod8 | 0.038 | 0.034 | 0.030 | 0.028 | 0.026 |
| σ_res_0mod5 | 0.090 | 0.060 | 0.049 | 0.038 | 0.030 |

## 5. Per-criterion verdicts

| criterion | ρ vs σ_to_1 | per-band ratio pattern | outcome | new info? |
|---|---|---|---|---|
| **σ_pow2_L10** | **+0.9999** | monotone ↑ to ~0.975 | **(a) PROPORTIONAL** | NO — gap = constant ~L steps; σ_pow2 ≈ σ_to_1 − L |
| **σ_prime** | **−0.022** | monotone ↓ (0.13 → 0.04) | **(c) INVISIBLE TO σ_to_1** | YES — captures EARLY dynamics decoupled from total length |
| σ_logdrop_k (k=1..8) | +0.34 → +0.76 | **U-shape** (lowest mid bands) | **(b) STRUCTURED** | partial — increasingly σ-correlated as k grows |
| σ_res_0mod3 | n/a | 0/100% binary | **STRUCTURAL INVARIANT** | YES — Collatz preserves {n ≢ 0 mod 3} |
| σ_res_1mod8 | −0.080 | monotone ↓ (0.08 → 0.02) | (c) early dynamics | YES — early-orbit residue arrival decoupled from total length |
| σ_res_5mod8 | +0.241 | monotone ↓ (0.04 → 0.03) | mostly (c) | partial |
| σ_res_0mod5 | −0.029 | monotone ↓ (0.09 → 0.03) | (c) early dynamics | YES |

## 6. Substantive findings

### A. The "not divisible by 3" Collatz invariant

σ_res_0mod3 reveals that {n ≢ 0 mod 3} is a Collatz-invariant set. This is provable in 2 lines but isn't surfaced by standard σ_to_1 measurement. Two-thirds of all orbits NEVER visit a multiple of 3 (after possibly leaving it on the first halving step from a start ≡ 0 mod 3).

### B. σ_prime is independent of σ_to_1

The first prime visit is essentially uncorrelated with total orbit length (ρ = −0.022). Mean σ_prime = 17 steps regardless of band. **σ_prime captures EARLY-orbit dynamics that are dynamically decoupled from late-orbit dynamics.** This is a clean outcome (c) — σ_prime is a substantively new lens.

Mechanism: primes are dense enough among small integers that any orbit hits one quickly via the descending halving phases. Once the orbit drops to small values, prime visits become near-certain. The "early prime visit" event has its own statistics independent of how the orbit eventually descends to 1.

### C. σ_pow2_L10 ≈ σ_to_1 − constant

When σ_pow2_L10 fires (only 4% of orbits), it fires very close to the end (ratio 0.92 → 0.97 across bands). The gap (σ_to_1 − σ_pow2_L10) is approximately constant ~L step units (since once you hit 2^L ≥ 1024, the orbit just halves down in L steps). Outcome (a) — same dynamics, additive constant offset.

### D. σ_logdrop_k shows U-shape in per-band ratios

For each k, the ratio σ_logdrop_k / σ_to_1 is highest in the lowest and highest σ_to_1 bands, lowest in the middle. Shape echoes the per-σ-quantile band U-shape in slopes (Result 36 follow-up 2). The same structural mechanism underlies both: middle-σ bands have orbits that descend "evenly" with average drift ≈ E[X], while extreme-σ bands have orbits with anomalous v-statistics.

### E. Residue-arrival criteria capture early dynamics

σ_res_1mod8, σ_res_5mod8, σ_res_0mod5 all show monotonically decreasing per-band ratios. Mechanism similar to σ_prime: early-orbit residue arrivals are independent of total length. These residue criteria act as "spread" measurements on early orbit dynamics.

## 7. Implications

**Most criteria (σ_prime, σ_residue_*, σ_pow2_L10) are NOT independent dimensions.** They either:
- Are proportional to σ_to_1 with constant offset (σ_pow2) — outcome (a)
- Capture EARLY-orbit dynamics that decouple from total length (σ_prime, σ_residue_*) — outcome (c) but not new structural constraints, just different time-windowed views

**The substantively new finding is the σ_res_0mod3 invariant** — a clean Collatz dynamics fact that doesn't appear in standard σ analysis. Worth promoting to its own structural observation: the {n ≢ 0 mod 3} subset is closed under T, two-thirds of all orbits live entirely in this subset.

**For trajectory-measure characterization (v3.6+):** the σ_to_1-based analysis (Result 32 family + Result 36 followups) captures the late-orbit / convergence-determining dynamics. Early-orbit dynamics (captured by σ_prime, σ_residue_*) are a separate observational lens but don't add structural constraints — they're independent measurements of a quickly-mixing process.

**For Collatz dynamics structure:** the "not divisible by 3" invariant is a candidate cleaning axiom — most theoretical work could be restricted to {n ≢ 0 mod 3} without loss.

## 8. Files

- `stopping_criteria_battery.py` — walker (~250 lines, includes Miller-Rabin)
- `stopping_criteria_battery.csv` — per-criterion summary
- `stopping_criteria_battery_log.txt` — full diagnostic log
- `stopping_criteria_battery.md` — this document
