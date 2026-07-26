# RESULT — PROBE SINGLEREC: the single-level recursion reproduces (machine precision) + P_k tabulated (2026-07-26)

**Probe:** `probe_singlerec.py`. Wilson's single-level recursion, the independence fork, and the measure-free carrier P_k.

## Independence fork — resolved from the code
`stationary_trunc` transitions `r → (3r+1)·2^{−v}` with weight `(0.5^v)/Z` **independent of r** (not the true
valuation `v₂(3r+1)`). So it is the Syracuse-random-variable model — exactly Wilson's assumption. The object matches.

## The recursion REPRODUCES to machine precision — and the session-long bug was the DFT sign
`π̂_k(ξ) = Σ_{a≥1} 2^{−a} e(ξ2^{−a}/3^k) π̂_k(3ξ2^{−a} mod 3^k)` (weight `0.5^a` real; phase/frequency `2^{−a}` =
`inv(2^a) mod 3^k`):

| k | rel |
|---|-----|
| 3 | 5.5e-16 |
| 4 | 5.7e-16 |
| 5 | 7.3e-16 |
| 6 | 9.0e-16 |

**Two bugs, compounded, explain every reproduction failure this session:**
1. **Wilson's form** — the n-from-(n−1) recursion (GATE_RECURSION, BRIDGE step 1) was wrong; the correct statement is
   **single-level** with the frequency map `ξ ↦ 3ξ2^{−a}` (raising v₃ by one each step ⟹ a finite expansion, not a
   fixed-point iteration — it unrolls to Tao's Syracuse offset map).
2. **My DFT sign** — numpy's `fft` uses the `e^{−2πi}` character, but the derivation's `e(ξX/3^k)` is `e^{+2πi}`. So
   `fft(dense)` was `π̂(−ξ)`; mixing it with the explicit `e^{+}` phases gave the flat rel~1.3 seen in GATE_RECURSION,
   RHOREC, and BRIDGE-1. Fixed by `π̂ = conj(fft(dense))`. (The bridge held regardless because it's self-consistent in
   either convention — which is why *it* reproduced and the recursions didn't, the clue I should have read sooner.)

With both fixed, the recursion is exact. **Correction to the record:** the earlier "recursion doesn't reproduce"
conclusions were these two bugs, NOT object mismatches. The genuine object findings from those probes stand (perpetuity
vs forward measure, additive fft(ν) vs multiplicative fft(ρ), R10's principal-unit coordinate) — but the recursion
itself is correct on the forward measure in the single-level form.

## P_k(ξ) — the measure-free contraction carrier, tabulated
`P_k(ξ) = |Σ_{a≥1} 2^{−a} e(ξ2^{−a}/3^k)| < 1` for ξ≠0 (deterministic, no ν — a geometric-weighted exponential sum
over the ⟨2⟩-orbit mod 3^k):

| k | max P | min 1−P | median P | P(ξ=1) | P(ξ=2) |
|---|-------|---------|----------|--------|--------|
| 3 | 0.789 | 0.211 | 0.595 | 0.595 | 0.232 |
| 5 | 0.944 | 0.056 | 0.568 | 0.634 | 0.245 |
| 6 | 0.971 | 0.029 | 0.570 | 0.586 | 0.303 |
| 8 | 0.994 | 0.006 | 0.570 | 0.499 | 0.376 |

Two facts for the contraction `|π̂(ξ)| ≤ |c|P_k(ξ) + V`, contract iff `V < S(1−P)`:
- **median P ≈ 0.57, stable across k** — at typical frequencies the deterministic room `1−P ≈ 0.43` is healthy and
  k-independent.
- **max P → 1**, with `1−P` shrinking `~2^{−k}` at the worst frequencies (0.21 → 0.006). So the deterministic factor
  alone loses its grip at the tightest frequencies as k grows; there the contraction budget `S(1−P) ~ S·2^{−k}` is
  small, and whether it holds depends on `V` being smaller still. **The whole question is now whether argmax|π̂| lives
  at a typical (P≈0.57) frequency or a worst (P→1) one** — a V-vs-S(1−P) comparison, fully computable and measure-aware.

## Net — both middle links of the chain are now closed
- **Recursion VALIDATED** (single-level, sign-fixed, rel~5e-16). The additive-side dynamics is a machine-precision
  identity.
- **Bridge VALIDATED** (BRIDGE2, rel~1e-13). The additive↔multiplicative crossing is exact.
- So `contraction → sup|π̂| (additive, validated) → Gauss bridge (validated) → sup|ρ̂| (channels)` is a closed chain.
- **P_k is measure-free and tabulated** — the escape/contraction is a `V < S(1−P)` comparison, no `(δ,A)` search. The
  live question: does the sup sit where `1−P` is healthy (~0.43) or where it vanishes (~2^{−k})? That's the next probe.

**Not at stake:** BRIDGE2, MAXMODE2/channels, MEAN1, HIERARCHY, CHANNEL_ID, R1–R30. Cheap (stationary_trunc, ~0.1s).
