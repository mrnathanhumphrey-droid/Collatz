# Result: top 10 eigenvalues of K_k^(q) across q ∈ {3,5,7,11,13}

**Date:** 2026-05-05.  Float64 dense eigensolve (scipy.linalg.eig).

## Configuration

| q | k | states | M = ord_{q^k}(2) | build (s) | eig (s) |
|---|---|---|---|---|---|
| 3 | 4 | 54 | 54 | 0.00 | 0.00 |
| 5 | 4 | 500 | 500 | 0.07 | 0.10 |
| 7 | 4 | 2058 | 1029 | 0.62 | 0.37 |
| 11 | 3 | 1210 | 1210 | 0.39 | 0.32 |
| 13 | 3 | 2028 | 2028 | 1.35 | 0.92 |

## Top 10 |λ_i| across q

| q | k | \|λ_1\| | \|λ_2\| | \|λ_3\| | \|λ_4\| | \|λ_5\| | \|λ_6\| | \|λ_7\| | \|λ_8\| | \|λ_9\| | \|λ_10\| |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 4 | 1.000000 | 0.000029 | 0.000029 | 0.000029 | 0.000029 | 0.000004 | 0.000004 | 0.000004 | 0.000001 | 0.000001 |
| 5 | 4 | 1.000000 | 0.000057 | 0.000057 | 0.000057 | 0.000057 | 0.000036 | 0.000036 | 0.000036 | 0.000036 | 0.000027 |
| 7 | 4 | 1.000000 | 0.000045 | 0.000045 | 0.000045 | 0.000045 | 0.000038 | 0.000038 | 0.000038 | 0.000038 | 0.000002 |
| 11 | 3 | 1.000000 | 0.000001 | 0.000001 | 0.000001 | 0.000001 | 0.000001 | 0.000001 | 0.000001 | 0.000001 | 0.000001 |
| 13 | 3 | 1.000000 | 0.000002 | 0.000002 | 0.000002 | 0.000001 | 0.000001 | 0.000001 | 0.000001 | 0.000001 | 0.000001 |

## λ_2 (rate-controlling) per q

| q | λ_2 (real) | λ_2 (imag) | |λ_2| | Δ from 1/2 | Δ from 1/q |
|---|---|---|---|---|---|
| 3 | -0.000029 | +0.000000 | 0.000029 | 0.5000 | 0.3333 |
| 5 | -0.000040 | +0.000040 | 0.000057 | 0.4999 | 0.1999 |
| 7 | -0.000045 | +0.000000 | 0.000045 | 0.5000 | 0.1428 |
| 11 | +0.000001 | +0.000000 | 0.000001 | 0.5000 | 0.0909 |
| 13 | -0.000002 | +0.000000 | 0.000002 | 0.5000 | 0.0769 |

## Spectral gaps

| q | \|λ_1\| − \|λ_2\| | \|λ_2\| − \|λ_3\| | \|λ_3\| − \|λ_4\| | \|λ_4\| − \|λ_5\| | \|λ_5\| − \|λ_6\| |
|---|---|---|---|---|---|
| 3 | 0.99997 | 0.00000 | 0.00000 | 0.00000 | 0.00003 |
| 5 | 0.99994 | 0.00000 | 0.00000 | 0.00000 | 0.00002 |
| 7 | 0.99995 | 0.00000 | 0.00000 | 0.00000 | 0.00001 |
| 11 | 1.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| 13 | 1.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |

## Verdict

**Headline:** at every q, **|λ_2| is essentially zero** (range 1.2×10⁻⁶ to 5.7×10⁻⁵), not 1/2. K_k^(q) has a Perron eigenvalue at 1 followed by a *huge* gap of essentially-1 down to a near-zero spectrum. This is the **q-universal** generalization of the q=3 finding from [R77.4 erratum](result_77_4_K_spectrum_erratum.md): K_k itself is **not the rate-controlling operator** for the convergence S_k → S_∞^{(q)} at any prime q.

| q | \|λ_2\| | Δ from 1/2 | Δ from 1/q | Δ from 0 |
|---|---|---|---|---|
| 3 | 2.9×10⁻⁵ | 0.4999 | 0.3333 | **2.9e-5** |
| 5 | 5.7×10⁻⁵ | 0.4999 | 0.1999 | **5.7e-5** |
| 7 | 4.5×10⁻⁵ | 0.5000 | 0.1428 | **4.5e-5** |
| 11 | 1.2×10⁻⁶ | 0.5000 | 0.0909 | **1.2e-6** |
| 13 | 2.1×10⁻⁶ | 0.5000 | 0.0769 | **2.1e-6** |

Across all five q tested:
- Δ from 1/2 ≈ 0.500 ± 0.0001 — **|λ_2| is NOT 1/2 at any q**
- Δ from 1/q is comparable in scale to 1/q itself — also not the rate
- Δ from 0 is in the 10⁻⁶ to 10⁻⁵ range, with q-specific value but always tiny

**This destroys the narrative that "1/2 is the rate operator's eigenvalue."** The 1/2 envelope of |ε_n|·2^n at q=3 (which itself was just falsified at k=7) was never reflected in K_3's spectrum. R77.4 erratum noted this at q=3; this probe extends it to q ∈ {5, 7, 11, 13} — the rate is q-universally **not in the per-level Markov spectrum**.

## Implications for the (q-3)/q closed form

The (q-3)/q closed form for `c̃_q = lim S_k^{(q)}/(q/3)^k` (q-sweep test 2, confirmed at q=11, 13, 17 within 1%) does **not** come from K_k^(q)'s spectrum at the algebraic level. Per the brief's decision rule:

> **|λ_2^(q)| varies with q: spectrum is q-dependent. (q-3)/q form is at the Plancherel mass level only, not the spectral level.**

Modified verdict: |λ_2|^(q) is q-universally near-zero, not at 1/2. The (q-3)/q form lives at a different level than K_k^(q)'s spectrum — likely the inter-level renormalization operator R77.4 erratum identified.

## Implications for q=7 anomaly

q=7's |λ_2| = 4.5×10⁻⁵ is *not* anomalous compared to q=5 (5.7×10⁻⁵). Both sit in the same scale band. The q=7 anomaly in c̃_7's deviation from (q-3)/q is **NOT visible at the spectral level**. q=7's irregularity is a Plancherel-mass-level phenomenon, not spectral.

## Cluster structure (group-theoretic)

The top eigenvalues come in clusters of equal magnitude (or near-equal):
- q=3 k=4: |λ_2| = |λ_3| = |λ_4| = |λ_5| ≈ 2.9×10⁻⁵ (4-cluster)
- q=5 k=4: |λ_2..5| ≈ 5.7×10⁻⁵ (4-cluster), |λ_6..9| ≈ 3.6×10⁻⁵ (4-cluster)
- q=7 k=4: |λ_2..5| ≈ 4.5×10⁻⁵ (4-cluster), |λ_6..9| ≈ 3.8×10⁻⁵ (4-cluster)
- q=11 k=3: |λ_2..10| ≈ 1.2×10⁻⁶ (broader cluster)
- q=13 k=3: |λ_2..4| ≈ 2.1×10⁻⁶, |λ_5..7| ≈ 1.1×10⁻⁶

Clusters reflect the group-theoretic structure of (Z/q^k)* — characters of the multiplicative group at non-trivial conductor levels. Cluster sizes are q-dependent (4 at q=3, 4 at q=5/7, larger at q=11/13).

## Files

- `result_qspectrum.py` — script
- `result_qspectrum.csv` — top 10 eigenvalues for each q
- `result_qspectrum.md` — this writeup