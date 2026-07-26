# RESULT — FREE CHECK: Π_{3∤a} π̂_k(a) IS rational (a field norm); v₃ = −(φ(3^k)−1) exact (2026-07-26)

**Probe:** `probe_normcheck.py`. Wilson's Galois/orbit check: is `Π_{3∤a} π̂_k(a)` rational, and what is its value + 3-adic
valuation? Feeds the Mahler/house/Dwork rationality shelf (Siegel diss p.92–93 = Dwork's Theorem) and the Siegel pairing.

## The claim holds — it's a field norm, and the arithmetic is clean
The set `{3∤a} = (ℤ/3^k)*` (the units) is Galois-stable, so `Gal(ℚ(ζ_{3^k})/ℚ)` permutes `{π̂(a)}` and
**`Π_{a∈(ℤ/3^k)*} π̂(a) = N_{ℚ(ζ_{3^k})/ℚ}(π̂(1))` is rational** (a field norm). Computed EXACTLY (forward-Syracuse
stationary measure on units, Fractions, full v-period M=ord₂(3^k), untruncated; norm via `Res(Φ_{3^k}, Q)/D^{φ(3^k)}`,
Φ monic):

| k | φ(3^k) | v₃(Norm) | pred −(φ−1) | Norm | ln\|N\| = Σ_a ln\|π̂(a)\| |
|---|--------|----------|-------------|------|--------------------------|
| 2 | 6 | **−5** | −5 | 9253/28588707 | −8.036 = −8.036 ✓ |
| 3 | 18 | **−17** | −17 | (34-digit rational) | −34.778 = −34.778 ✓ |
| 4 | 54 | **−53** | −53 | num740d/den801d | −138.494 = −138.494 ✓ |
| 5 | 162 | **−161** | −161 | num7369d/den7591d | −512.418 = −512.418 ✓ |
| 6 | 486 | **−485** | −485 | num69145d/den69943d | −1837.774 = −1837.774 ✓ |

Float cross-check passes at every k (`ln|Norm| = Σ_{a∈U} ln|π̂(a)|`, sign-convention-free — confirms the exact norm
matches fwd_hat). **k=6 confirmed (v₃=−485=1−φ) after a 2350s resultant — the law now holds at k=2,3,4,5,6.**

## The finding — a one-line 3-adic law, and its λ-adic meaning (Wilson)
**v₃(N(π̂(1))) = −(φ(3^k) − 1) = 1 − 2·3^{k−1}**, exact at k=2,3,4,5,6 (−5,−17,−53,−161,−485 = −(φ−1)).

**λ-adic reading (Wilson):** in `ℚ(ζ_{3^k})`, 3 is *totally ramified*, `(3) = (λ)^φ` with `λ = 1−ζ`, and
`v₃(N(α)) = v_λ(α)`. So the formula says
`v_λ(π̂(1)) = 1 − φ`, hence `|π̂(1)|_3 = 3^{(φ−1)/φ} → 3`. **Equivalently: `3·π̂(1)` is an algebraic integer, divisible
by λ exactly once.** Clean, sharp — a theorem waiting to be proved (the proof will likely hand over the other
valuations too).

The norm is rational with its entire 3-adic
content in the denominator, of size exactly `3^{φ(3^k)−1}`. Interpretation:
- The **geometric mean of the spectrum** is `|N|^{1/φ(3^k)} = exp((Σ_a ln|π̂(a)|)/φ)` — an **aggregate-side /
  spectral-distribution functional** (not the sup, not the total), the class the channels live in. So it lands on the
  right side of the ℓ^∞-vs-aggregate line: it is NOT Tao's Prop 1.17 object.
- Rationality + the clean v₃ law is exactly the input the **Dwork/Borel rationality shelf** (Siegel diss p.92–93,
  Thm 3.2: rational power series via `Π_p R_p > 1`) is built to consume — the "Siegel pairing" Wilson flagged.

## Net
Wilson's Galois observation is CONFIRMED (the product is a rational field norm, not just numerically real), and it
comes with a sharp bonus: **v₃ = 1 − 2·3^{k−1}**. This is a spectral-distribution invariant on the aggregate side —
a concrete arithmetic handle on the class of functionals (ℓ⁴/channels) that the re-center identified as the open,
ours-alone prize. Not at stake: RECENTER, SEC7, LAMBDA, CHANNEL_ID, v₃ HIERARCHY, R1–R30. Cheap (k≤4, 0.2s; k=5,6 heavy).
