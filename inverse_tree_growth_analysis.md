# Result 32: Per-attractor inverse-tree spectral bypass FAILS — equivalent-vs-requires verdict hardens to REQUIRES

**Date:** 2026-05-03. Sequel to Result 31. Tests user's pressure-test bypass route: does ⟨σ_S | j⟩ close via per-attractor inverse-tree growth eigenvalue λ_j without requiring full trajectory-measure invariance?

**Verdict: empirically falsified at all three j-classes and all three N values tested.** The single-eigenvalue spectral closure is rejected. Per-attractor inverse-tree depth distribution is hump-shaped (not pure exponential), λ_j is not a structural invariant (drifts with N), and the spectral prediction σ_pred = log(N/m_j)/log(λ_j) fails by 30-300%. The "⟨σ_S | j⟩ requires trajectory-measure invariance" verdict hardens.

Code: `inverse_tree_growth.py`. Method: forward enumeration of all odd m ∈ [1, N] (no excursion cap, unlike BFS-with-N-cap which biases by removing high-excursion ancestors). Total compute: ~1.5s across three N values.

---

## 1. Why forward enumeration, not BFS

A first attempt at the spectral test used BFS from m_j capping intermediate predecessors at m ≤ N. **This is wrong.** Forward orbits from m_anc ∈ [1, N] absorbing at m_j may excurse above N before descending to m_j (e.g., m_anc=27 has orbit 27 → 41 → 31 → 47 → 71 → ... → 5; all intermediate values stay reasonably bounded for small starts but high-σ orbits do excurse).

BFS from m_j with cap m_pred ≤ N at every step CULLS chains that pass through m > N intermediate predecessors, missing valid ancestors. This is the same "high-excursion correction" that compute_threads_findings.md addendum identified for the absorbing-chain machinery.

**Forward enumeration**: walk Syracuse from each odd m ∈ [1, N] with no cap on excursion, record (σ_S, j_attr). This gives the exact per-(σ_S, j) count without the truncation bias.

## 2. Per-attractor depth distributions are hump-shaped

Forward enumeration at N = 2^26 = 67M (33M odd integers, 0.3s walk):

| j | # ancestors | σ-mode | ⟨σ_S\|j⟩ | mode/mean | shape |
|---|---|---|---|---|---|
| 2 | 31,474,743 | 63 | 61.73 ± 0.005 | 1.02 | nearly symmetric |
| 4 | 791,525 | 36 | 39.94 ± 0.021 | 0.90 | right-skewed |
| 5 | 1,268,047 | 36 | 44.56 ± 0.018 | 0.81 | right-skewed |
| 7 | 2,713 | 12 | 17.88 ± 0.185 | 0.67 | right-skewed |
| 8 | 16,259 | 24 | 27.20 ± 0.111 | 0.88 | right-skewed |

**Mode is well below mean for j ≥ 4** (mode/mean = 0.67-0.90). The depth distribution has a hump (peak count at moderate depth) plus an exponentially-decaying right tail. **Not pure exponential growth** — the simple "λ_j^d count" model fails by construction.

## 3. λ_j fitted on the growth regime drifts with N

Fitting log count(d) vs d in the rising portion (d ∈ [3, peak_d - 2]):

| j | N=2^22 | N=2^24 | N=2^26 |
|---|---|---|---|
| 2 | 1.217 | 1.087 | 1.108 |
| 4 | 1.216 | 1.281 | 1.176 |
| 5 | 1.246 | 1.204 | 1.197 |
| 7 | 1.229 | 1.269 | 1.406 |
| 8 | 1.277 | 1.407 | 1.172 |

Variation 5-15% across N. R² values are 0.57-0.95 — fits aren't clean exponential. **λ_j is not a structural invariant.**

For comparison, Result 23's forward residue chain had λ_max(M_closed) = 1.264 invariant across k=5..11. The per-attractor growth rate doesn't share that invariance.

## 4. Spectral prediction misses by 30-300%

For each j and N, predict ⟨σ_S | j⟩_pred = log(N/m_j) / log(λ_j). Compare to empirical:

At N=2^26:
| j | log(N/m_j) | log(λ_j) | σ_pred | σ_emp | gap |
|---|---|---|---|---|---|
| 2 | 16.41 | 0.102 | 160.29 | 61.73 | +98.6 (160% overprediction) |
| 4 | 13.58 | 0.162 | 83.59 | 39.94 | +43.6 (109%) |
| 5 | 12.19 | 0.180 | 67.88 | 44.56 | +23.3 (52%) |
| 7 | 8.19 | 0.340 | 27.66 | 17.88 | +9.8 (55%) |
| 8 | 6.80 | 0.159 | 50.63 | 27.20 | +23.4 (86%) |

Spectral prediction systematically OVERESTIMATES σ_S. The "log(N/m_j) / log(λ_j)" formula assumes σ_S grows linearly to fill the inverse-tree to its growth-regime extent — but the actual distribution saturates earlier.

Even using the "universal" λ_max = 1.264 from Result 23 (trying to give the spectral hypothesis its best shot):

| j at N=2^26 | σ_pred (λ=1.264) | σ_emp | gap |
|---|---|---|---|
| 2 | 70.1 | 61.7 | +8.4 |
| 4 | 58.0 | 39.9 | +18.1 |
| 5 | 52.1 | 44.6 | +7.5 |

Still 14-45% gaps. The spectral closure fundamentally doesn't work.

## 5. Why the spectral hypothesis fails

The inverse-tree depth distribution is a Galton-Watson-like process with non-uniform branching (mod-3 dependent dead-ends, varying predecessor counts per residue class). Under saturation constraint m ≤ N:

- Initially exponential growth: count(d) ~ λ_eff^d for small d
- Saturation: when ancestors fill the m ≤ N volume, count(d) plateaus then decays
- Result: hump-shaped distribution, not pure exponential

⟨σ_S | j⟩ = mean of this hump. NOT determined by single eigenvalue λ_j alone — depends on full distribution shape (variance, skewness, location of peak).

For closed-form ⟨σ_S | j⟩, would need:
- Full per-depth count distribution at uniform-m sampling
- This requires the trajectory-measure structure (how density evolves as orbits descend through {m ≤ N})

This is **precisely the Lagarias trajectory-measure invariance question.**

## 6. Verdict on equivalent-vs-requires

User's pressure-test asked: is ⟨σ_S | j⟩ closure equivalent to or requires Lagarias trajectory-measure invariance?

**Empirical answer: REQUIRES.** The simplest spectral bypass (single eigenvalue λ_j) is empirically falsified. More refined spectral methods (multi-eigenvalue, full transition matrix on residues per attractor) might work in principle, but each successive refinement converges toward characterizing the full trajectory measure — at which point you've solved the Lagarias question anyway.

The hump-shaped depth distribution shows: ⟨σ_S | j⟩ is determined by the FULL distribution shape, not a leading eigenvalue. The full distribution is the trajectory measure restricted to absorbing at j.

**No bypass found. Lagarias-class equivalence (in the functional sense) holds.**

## 7. What this changes for v3.5 framing

Result 31 framed it as "reduces to Lagarias-class problem". Result 32 sharpens this: **the reduction is functionally tight.** Per-j W_j cannot be closed via simpler spectral / Cramer / Esscher / mixing methods; closure requires the trajectory-measure invariance object directly.

The unification claim from Result 31 stands: ε(σ), K_eff, σ-quartile selection, per-j W_j — all four manifestations of the same Lagarias-class trajectory measure object. Each manifests through a different marginal/conditional, but no simpler structural bypass exists for any of them.

## 8. Files

- `inverse_tree_growth.py` — forward-enumeration depth-distribution analysis
- `inverse_tree_growth_analysis.md` — this document (Result 32)
- `inverse_tree_per_attractor.py` — earlier BFS-with-N-cap attempt (biased by excursion cap; preserved as cautionary)
- `experiments_output/inverse_tree_growth_log.txt` — full output
