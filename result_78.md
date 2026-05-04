# Result 78: Cochrane × Kalafatelis attack on eq 190 — outcome (γ), structural obstruction documented

**Date:** 2026-05-04. Closes the Bilinear_Estimates reference-library attack initiated by brief v2.

## Verdict: outcome (γ)

**Cochrane Theorem 2 does NOT directly close Kalafatelis's eq 190** for our specific setup. The obstruction is structural and identified precisely.

The polynomial identification of Kalafatelis's S_{r,ℓ,ε}(m) via binomial expansion of (1+3)^u modulo 3^{r+1} produces a polynomial g(u) whose derivative g'(u) has the property:

> **For every j ≥ 0, coefficient of u^j in g'(u) has v_3 ≥ j+1 − v_3((j+1)!) ≥ 1.**

Factoring g'(u) = 3^τ · H(u) with the largest possible τ = 1 (set by the constant term's v_3), the resulting H(u) mod 3 has **only its constant term non-zero**. Therefore **D := degp H+ = 0 in Cochrane's notation**.

Cochrane Theorem 2 with D = 0 gives "sum vanishes by counting H(a) ≡ 0 condition" — but this is a complete-sum result. Our incomplete sum (range 3^{r-1}, modulus 3^{r+1}) does NOT inherit this trivial vanishing. The mismatch is the structural obstruction.

**Honest characterization:** Cochrane's machinery is designed to extract cancellation from polynomial structure of the character / phase modulo p. The Syracuse map's specific 3-adic structure of (1+3)^u places ALL cancellation at the constant level after scaling by 3^τ — the higher-derivative polynomial behavior is "trivial mod 3." Cochrane's bound on D = 0 sums is the trivial bound, not the non-trivial saving we need.

## What was tested (rigorous calculation)

Levels r = 2, 3, 4. Polynomial identification:
- f(u) = 1 (constant; trivial character on phase-only sum)
- g(u) = c · Σ_{k=0}^r C(u, k) · 3^k − 9m · u
- where C(u, k) = u(u−1)…(u−k+1) / k! and c = c_{ℓ,ε} ∈ (Z/3^{r+1})*

Verified at r = 2, 3, 4 by direct computation of g'(u) coefficients with their 3-adic valuations:

| r | u^0 v_3 | u^1 v_3 | u^2 v_3 | u^3 v_3 | τ | D |
|---|---|---|---|---|---|---|
| 2 | 1 | 2 | — | — | 1 | 0 |
| 3 | 1 | 2 | 3 | — | 1 | 0 |
| 4 | 1 | 2 | 3 | 3 | 1 | 0 |

Each level's H(u) mod 3 is a constant (only u^0 has v_3 = 0 after factoring τ = 1).

## Why this happens (structural explanation)

The binomial expansion (1+3)^u = Σ_{k=0}^r C(u, k) · 3^k mod 3^{r+1} has coefficients C(u, k) at the term 3^k. Differentiating in u: coefficient of u^j in g'(u) gets contributions from k ≥ j+1, each contributing 3^k · (rational with bounded v_3 of denominator).

The dominant contribution to coefficient of u^j is from k = j+1: contribution = 3^{j+1} · (1/j!).

So v_3 of coefficient of u^j ≥ j+1 − v_3(j!) ≥ 1 (using Legendre v_3(j!) < j/2).

After factoring τ = 1, the resulting H(u) has coefficient of u^j with v_3 ≥ j+1 − v_3(j!) − 1 ≥ 0. The constant term has v_3 = 0 (non-zero mod 3); higher-degree terms have v_3 ≥ 1 (zero mod 3).

This is a generic feature of p-adic exponentiation. The "p-adic logarithm" expansion converges p-adically with each successive term gaining a factor of p. The resulting polynomial has 3-adic structure concentrated at the lowest order.

## Step 4 attempt: Cauchy-Schwarz / smoothing — also obstructed

We attempted to bridge the incomplete-to-complete gap via Fourier smoothing:
> Σ_{u=0}^{N−1} f(u) = (1/q) Σ_{ξ mod q} 1̂(ξ) · F̂(−ξ)

where F̂(ξ) = Σ_{u mod q} f(u) e(uξ/q).

Each F̂(ξ) is a complete sum that adds ξ·u to the linear term of g(u). This shifts the constant of g'(u) by ξ but leaves all higher-degree terms unchanged. Therefore **D remains 0 for every ξ** — Cochrane's bound on F̂(ξ) is still trivially "vanishing".

The Pólya-Vinogradov bound |F̂(ξ)| ≤ √q · log(q) ≈ 3^{(r+1)/2} · (r+1) gives:
> |S_{partial}| ≤ (1/q) · q · log(q) · sup |F̂| ≤ log(q) · √q = (r+1) · 3^{(r+1)/2}

Compared to the trivial bound |S_{partial}| ≤ 3^{r-1} (length of sum), Pólya-Vinogradov is **weaker** for our regime: 3^{(r+1)/2} > 3^{r-1} when (r+1)/2 > r−1, i.e., r < 3. So Pólya-Vinogradov gives no saving for r ≥ 3.

For r ≥ 3, neither Cochrane Theorem 2 (D = 0 obstruction) nor Pólya-Vinogradov (trivial bound) yields the η^{1/2+δ} saving Kalafatelis needs.

## Comparison to Kalafatelis's status

Kalafatelis explicitly leaves eq 190 as **open** (Remark 27 of his paper). Our negative result with Cochrane × Kalafatelis is consistent with the difficulty he identifies.

What we've demonstrated: **the most natural application of Cochrane's machinery (polynomial identification via binomial expansion of (1+3)^u) does not close eq 190**, due to the specific structural obstruction D = 0 in Cochrane's notation.

This is a real research-level finding. The identification was the natural first attempt; the obstruction is non-trivial; the path forward requires different machinery.

## Path forward: what would close eq 190

Three plausible alternatives, listed in increasing order of likely difficulty:

1. **Sharper p-adic identification.** Instead of polynomial expansion of (1+3)^u, use a different parametrization that captures the actual 3-adic structure. Candidates: (a) Kalafatelis's own Postnikov-style λ = log_3(4) approach with a specific cancellation analysis; (b) a 3-adic Steinitz / Iwaniec "two-power" decomposition isolating non-trivial cancellation higher up.

2. **Smooth completion of 3^{r+1}**. Cochrane Theorem 1 (the headline result) requires smooth modulus (multiple distinct prime factors). One could try working modulo 3^{r+1} · q for some auxiliary prime q, then averaging q out. This loses some cancellation but might give the η^δ saving.

3. **van der Corput / Weyl differencing**. Apply repeated differencing to reduce the bilinear sum's effective rank, then apply Cochrane Theorem 2 to the reduced form. Standard technique in analytic number theory but typically loses more than it gains for shallow polynomial structures.

None of these are immediate; each is a substantial research project.

## What this means for c = 7/45

Result 78's outcome (γ) does NOT undermine c = 7/45's validity. The empirical evidence (k = 1..6 exact rationals + |ε_n|·2^n stable envelope + sharpened conjecture S_n = 7/15 − (1/30)·(1/2)^n + O((1/4)^n)) remains as strong as before.

What it does mean: **the rate-½ rigorous proof remains open**. The Cochrane reference-library attack identified the exact Theorem 2 + Prop 4 + Cor 6 machinery as the most promising direct route, attempted it, and found the structural obstruction. This narrows the open question and clarifies what new machinery is needed.

## Strategic position

Pre-Result 78: c = 7/45 had 5 rigorous structural identities (R74, R75 Theorems 75.1/75.2, R76 Theorems 76.1/76.3) + 1 derived operator (R77 T_diag). Rate ½ was empirically certified through k=6 with envelope |ε_n|·2^n ≤ 0.04.

Post-Result 78: The Cochrane reference-library route is **closed as a direct attack**. The status of c = 7/45 is unchanged in terms of empirical certification, but the analytical closure path now has an explicit obstruction map: any future approach must either bypass the D = 0 structure of polynomial expansion, or use machinery genuinely different from Cochrane Theorem 2.

The most likely productive direction is probably (1) above — Kalafatelis's own λ = log_3(4) approach with sharper cancellation analysis. His Remark 27 effectively says the same thing: he knows it's a real obstruction and identifies it as outstanding work in his paper.

## Files

- `result_78_polynomial_identification.py` — Step 1 implementation, runs at r = 2, 3, 4
- `experiments_output/r78_polynomial_identification.txt` — saved analysis (UTF-8 issue noted, content verified inline)
- `cochrane2026.txt` — extracted Cochrane paper text (2150 lines)

## Updated open problems for c = 7/45

| Problem | Status |
|---|---|
| Plancherel formula for S_k | ✓ Proved (R75) |
| Tao recursion → diagonal/off-diagonal | ✓ Proved (R75) |
| Conservation law Σ_j M_{n+1}(η_0+j·3^n) = 0 | ✓ Proved (R76) |
| Leading-mode identity S_{n+1} = −2·M_{n+1}(1+3^n) | ✓ Proved (R76) |
| Class collapse P^{+−} = 0 for n ≥ 2 | ✓ Proved (R76) |
| T_diag spectrum {0, 1} | ✓ Proved (R77) |
| (1, 4) deviation eigenvector | ✓ Verified (R76, R77) |
| Off-diagonal rate λ_2 = 1/2 | ◐ Empirical (k=2..6); Cochrane attack obstructed (R78) |
| Coefficient 1/30 = S_∞/14 | ◐ Numerical fit; analytical origin open |
| Kalafatelis eq 190 | ✗ Open (per his Remark 27 + R78) |

**Result 78 conclusion:** Cochrane × Kalafatelis is closed as a direct attack route. The path forward is either Kalafatelis's own approach (his Remark 27 work) or a fundamentally different technique. c = 7/45 remains empirically certified to ≤ 4×10⁻⁴ at k=6, with multiple rigorous structural anchors. The single remaining gap is the off-diagonal rate-½ proof, which is a published open problem in Kalafatelis's 2026 paper.
