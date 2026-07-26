# RESULT — PROBE CONTRACTION: the V<S(1−P) contraction is REFUTED; the sup sits AT the P→1 tight frequency (2026-07-26)

**Probe:** `probe_contraction.py`. Wilson's contraction table: for k=3..10, locate ξ*=argmax|π̂| and read `V / (S(1−P))`
at the frequency carrying the sup. Pre-registered prediction (Wilson): ξ* **odd**, v₂(ξ*)=0, P(ξ*)~median, **ratio<1**
(contraction closes, carries to sup|ρ̂| via BRIDGE2).

## The pre-registration is DEFEATED — cleanly, and the inversion is the finding
| k | ξ* | ξ*/3^k | v₂ | odd? | P(ξ*) | medP | S | V | **V/(S(1−P))** |
|---|-----|--------|----|----|-------|------|---|---|-------|
| 3 | 8=2³ | 0.296 | 3 | EVEN | 0.789 | 0.595 | 0.25224 | 0.215 | **4.04 DEFEATS** |
| 4 | 65=81−2⁴ | 0.802 | 0 | odd | 0.887 | 0.567 | 0.17700 | 0.095 | **4.77 DEFEATS** |
| 5 | 32=2⁵ | 0.132 | 5 | EVEN | 0.944 | 0.568 | 0.12927 | 0.053 | **7.20 DEFEATS** |
| 6 | 64=2⁶ | 0.088 | 6 | EVEN | 0.971 | 0.570 | 0.09611 | 0.034 | **12.24 DEFEATS** |
| 7 | 256=2⁸ | 0.117 | 8 | EVEN | 0.986 | 0.568 | 0.07587 | 0.028 | **25.98 DEFEATS** |
| 8 | 6049=3⁸−2⁹ | 0.922 | 0 | odd | 0.994 | 0.570 | 0.06089 | 0.024 | **60.6 DEFEATS** |
| 9 | 1024=2¹⁰ | 0.052 | 10 | EVEN | 0.997 | 0.570 | 0.04803 | 0.019 | **123 DEFEATS** |
| 10 | 54953=3¹⁰−2¹² | 0.931 | 0 | odd | 0.997 | 0.570 | 0.03828 | 0.016 | **126 DEFEATS** |

**The sup sits at a PURE POWER OF 2** (ξ*=2^m, exponent m≈k drifting slightly above k) — its conjugate mirror N−2^m
is the odd twin (|π̂(ξ)|=|π̂(−ξ)| exactly, so the argmax is the degenerate pair {2^m, N−2^m}). This is EXACTLY the
maximal-v₂, **P→1 tightest** frequency — the OPPOSITE end from Wilson's "odd, v₂=0, P~median" prediction. Wilson's
mechanism (even/high-v₂ ⟹ genuine halvings ⟹ P→1) is CORRECT; his guess about *where the sup lives* is backwards —
the sup sits precisely on the aligned frequencies where the deterministic factor gives **zero room**. The ratio blows
up (4→126) not because V is large but because `1−P ~ 2^{−k} → 0` there.

**Correction to the record:** Wilson's banked a_max (0.259=7/27 at k=3, etc.) is NOT the argmax — 7=2⁻² (odd, v₂=0) is
only **rank-5** at k=3 (|π̂|²=0.0306 vs the max 0.0636). The correctly-signed forward π̂ (= R66 by value:
max|π̂|²=0.06362 vs banked R66 0.06360, matching to 4 digits at every k) peaks at ξ=2^m. The old a_max locations were
mislabelled (pre-sign-fix / different measure); the R66 *values* are right, the *locations* were not.

## Why the bound is the wrong instrument — the feed-values are LARGER than S
The recursion `π̂(ξ*) = Σ_a 2^{−a} e(ξ*2^{−a}/3^k) c_a`, `c_a = π̂(3ξ*2^{−a} mod 3^k)` feeds ξ* (v₃=0) from v₃≥1
frequencies. Those feed-values are **larger** than the sup, not smaller:
| k | max\|c_a\|/S | wmean\|c_a\|/S | \|c\|/S | 1−P |
|---|-----------|--------------|-------|-----|
| 3 | 1.50 | 1.22 | 0.79 | 0.211 |
| 6 | 1.35 | 1.04 | 1.01 | 0.029 |
| 10 | 1.25 | 1.11 | 1.03 | 0.003 |

`max|c_a|/S > 1` at every k (1.50→1.25): the v₃=0 sup is fed by higher-v₃ (near-DC) frequencies of **larger** modulus.
With `|c|/S ≈ 1` and `P ≈ 1`, the fact that `|π̂(ξ*)| = S` (smaller than its inputs) is **phase cancellation among
genuinely large terms** — and a modulus bound `|c|P + V` throws away exactly the phase alignment that produces the
cancellation. So the bound is not merely loose; it is structurally blind to the mechanism. sup|π̂| **does** decay
geometrically (S: 0.252→0.038 over k=3..10, ratio ~0.80/step = R66's rate), but that decay lives in the **phase (arg)
structure of π̂ at the near-DC feed points**, not in the deterministic carrier P (which is ≈1 exactly at the sup).

## Verdict — the contraction route as posed is dead; the seam moves to the phase of π̂
- The `V < S(1−P)` contraction **cannot close** — the sup deliberately occupies the P→1 frequencies where the
  deterministic room `S(1−P) ~ S·2^{−k}` vanishes, while V decays only ~S. Refuted at every k, ratio growing.
- The real decay of sup|π̂| is carried by **phase cancellation among near-DC feed-values of modulus > S** — an
  `arg π̂` phenomenon, not a modulus/`P` phenomenon. (Contrast probe_alpha, which found `arg` "nondescript" globally;
  the live object is `arg π̂` specifically at the feed points `3·2^m·2^{−a}`, near DC.)
- Chain status unchanged upstream: **recursion (SINGLEREC) and bridge (BRIDGE2) still machine-verified.** What's
  refuted is the *closing argument* (deterministic-factor contraction), not the links. The sup→sup|ρ̂| crossing is
  intact; it just needs a phase-aware bound on sup|π̂|, which is Wilson's pen.

**Not at stake:** SINGLEREC, BRIDGE2, MAXMODE2/channels, MEAN1, HIERARCHY, CHANNEL_ID, R1–R30. Cheap (~0.4s).
