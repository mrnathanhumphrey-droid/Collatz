# RESULT — PROBE GATE_RECURSION: the gate caught a convention mismatch (2026-07-26)

**Probe:** `probe_gate_recursion.py`. Wilson's instruction: validate his μ̂ recursion against R66's μ̂_k (k≤6) BEFORE
building the (δ,A) escape scan; if it doesn't reproduce, the recursion is his and must be restated in R66's indexing.
I added a second question: is R66's μ̂ even the same Fourier object as the channel-binding object `fft(ρ)`?

## Finding 1 — R66's μ̂ is `fft(ν)` [raw/additive], NOT the channel object `fft(ρ)` [dlog]

| k | max\|fft(ν)\|² (raw) | max\|fft(ρ)\|² (dlog=channel) | R66 max |
|---|---------------------|-------------------------------|---------|
| 2 | 0.33333 | 0.14283 | 0.14280 |
| 3 | 0.14283 | 0.05171 | 0.06360 |
| 4 | 0.06362 | 0.02004 | 0.03130 |
| 5 | 0.03133 | 0.00778 | 0.01670 |
| 6 | 0.01671 | 0.00450 | 0.00924 |

**R66's μ̂ = `fft(ν)` shifted one level: `R66_k = fft(ν)_{k+1}` exactly** (0.1428=fft(ν)₃, 0.0636=fft(ν)₄, …). So
R66 is the **raw/additive** Fourier transform of ν (build_nu's ν matches R66 up to this one-level offset — the measure
is right).

**But the channel-binding object is `fft(ρ)` [dlog], a THIRD sequence** — MAXMODE2 verified `mean_m(γ_k(m)−1)² =
Σ_a|fft(ρ)(a)|⁴` (the U² identity), so the persistence channels are the Fourier-dual of `|fft(ρ)|²`, not `|fft(ν)|²`.
And `fft(ρ)` decays FASTER (max power ratio ~0.42, amplitude 0.647 — the saturating rate) than `fft(ν)`/R66 (~0.5,
"2^{−k}"). **They are genuinely different transforms, related by the discrete-log (dlog) reindexing** — at k=2 the
value-sets coincide up to the dlog permutation (both max 0.14283), but at k≥3 they diverge.

**⟹ Hank's identification "R66's max mode = the binding mode" conflated two different Fourier objects.** The binding
object is `fft(ρ)` (dlog/multiplicative), max saturating at 0.655; R66's `fft(ν)` (raw/additive) is a slower, distinct
sequence. The lemma must bound **`fft(ρ)`**, not R66's μ̂.

## Finding 2 — Wilson's recursion, faithfully implemented, reproduces NEITHER object

`μ̂_n(ξ) = Σ_{a≥1} 2^{−a} e(ξ2^{−a}/3^n) μ̂_{n−1}(ξ2^{−a} mod 3^{n−1})` (2^{−a} = modular inverse mod 3^n; argument
reduced mod 3^{n−1}; a=1..64):

| n | rel error on fft(ν) | rel error on fft(ρ) |
|---|---------------------|---------------------|
| 3 | 1.57 | 1.55 |
| 4 | 1.43 | 1.20 |
| 5 | 1.38 | 1.41 |
| 6 | 1.36 | 1.32 |

rel ~1.3–1.6 (≫1) on both — the recursion as written does not reproduce either sequence. Since build_nu's ν *is*
R66's measure (up to the one-level offset), the failure is not the measure; it's the recursion's **level indexing or
modulus bookkeeping** (the `R66_k = fft(ν)_{k+1}` offset strongly suggests the recursion's "n" is shifted relative to
build_nu's level, or the `mod 3^{n−1}` vs `mod 3^n` reduction / the `2^{−a}` inverse modulus needs adjustment).

## Verdict — stop at the gate, hand back to the pen

Per Wilson's own criterion ("if it doesn't reproduce, it's mine and must be restated"), **the escape scan is NOT built**
— building it on an unvalidated recursion and/or the wrong Fourier object is exactly the PRODFORM/R29-A reconstruction
failure mode. Two things for the pen before the scan:

1. **Pin the object.** The persistence channels are bound by `fft(ρ)` (dlog/multiplicative Fourier), max saturating at
   0.655. R66's μ̂ and the raw-Syracuse recursion are `fft(ν)` (raw/additive), a slower, different sequence. The escape
   criterion — whose whole content is the "×2-orbit of ξ escapes the high-set" — is a statement about the **raw**
   ξ↦ξ·2^{−a} orbit (fft(ν)). Under dlog, ×2 becomes a *shift*, so the escape structure for `fft(ρ)` is a **different**
   combinatorial object. Decide which object the lemma is written for, and if `fft(ρ)`, restate the ×2 orbit accordingly.
2. **Fix the recursion indexing.** As implemented it's off by rel~1.4; the `fft(ν)_{k+1}=R66_k` offset is the clue.
   Once it reproduces the intended object (k≤6), the (δ,A) scan is safe to run — and it's cheap the moment the object
   and recursion are pinned.

**Not at stake:** MAXMODE2's saturation (fft(ρ), the channel object — solid), MEAN1, HIERARCHY, CHANNEL_ID, R1–R30.
Cheap (build_nu(6), 0.1s). The gate cost nothing and saved building on the wrong object.
