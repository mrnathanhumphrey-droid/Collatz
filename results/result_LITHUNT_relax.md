# RESULT — Hank re-task (relaxation-bound lit hunt) + MAXMODE gate (2026-07-26)

## Hank's lit hunt: bound Σ_j|q_j(k)−1/3| for k=3,4 (the persistence inequality's binding channels)

**Bottom line: nothing off-the-shelf is both the right observable and sharp. The decay literature quotes the WRONG
Fourier mode.** The target is a per-level geometric relaxation rate in a narrow window (amplitude ≤ ~0.79, from the
budget: k=3 has Σ|q_j−1/3|≈0.0945, leading term ≈0.0255, sufficient condition Σ≤0.12 ⟹ rate ≤ 1−0.0255/0.12 ≈ 0.79).

Ranked:
1. **BESPOKE — R66/R74 ν-weighted lifting-operator spectral gap (on disk).** The only object of the right type AND
   observable AND rate. R74 has the exact Parseval recursion `S_{k+1}=3^{k+1}‖d_{k+1}‖²`; the binding rate is ν's
   **max (slowest) Fourier mode**, not the average. Missing lemma (for the pen): top eigenvalue of the
   deviation-projected ν-weighted lifting operator on the k-channel ≤ ~0.75, equivalently `max_a|ρ̂(a/3^j)|²` decays
   ≤ ~0.55/level. Finite-dim linear algebra on an explicit chain; the obstacle is π-conservation leakage
   contaminating the finite-level top singular value (R74 open item, "clean basis likely explicit via K₂ rank-2").
2. **Diaconis–Fulman / Holte carries chain** (Holte 1997 Amer.Math.Monthly; DF 2009 arXiv:0806.3583 §5.2 = mult-by-4
   base 3). Eigenvalues = negative powers of base ⟹ rate exactly **1/3**. RIGHT TYPE but WRONG REGIME: it's the
   *uniform-input* (zero-baseline) relaxation; 1/3 is too fast ⟹ **false as an upper bound**. Value = the exact
   transfer skeleton the ν-weighting perturbs.
3. **Tao Prop 1.14/1.17** (arXiv:1909.03562, Forum Math Pi 2022). Superpolynomial `n^{−A}` by design — the read Hank
   was asked for: the proof does NOT secretly give geometric decay (entropy-decrement loses a factor per scale, no
   spectral gap), and the constant `C_A` is non-explicit. Too weak (polynomial-vs-geometric, unbounded margin).
4. **BGK / Heilbronn / Bourgain mod-p^r** — wrong regime (large-subgroup cancellation, non-explicit ε, no per-level
   rate; magnitude-only, cancels in the ratio).

**The key structural claim (Hank's synthesis, flagged as inference):** the binding channel relaxation is governed by
ν's MAX Fourier mode (~0.707 amplitude), not the average (1/3, what every external theorem quotes). Hank asked this
be confirmed numerically first.

## MAXMODE gate — Hank's claim CONFIRMED qualitatively, quantified honestly

`probe_maxmode.py`: ν power spectrum `|ρ̂_j(m)|²` over primitive m (3∤m), max vs average decay, vs channel rates.

- **Average power rate = 0.334 ≈ 1/3, dead flat across all levels j=3..16.** This is exactly the external-theorem /
  DF rate — and it is the WRONG (too-fast) mode. Confirmed clean.
- **Max mode decays SLOWER: power rate median 0.42, amplitude 0.647** (per-level amplitude ranges 0.60–0.76, worst
  spike 0.76 at j=6, next 0.72 at j=13). So Hank's core diagnosis holds: **the max mode is genuinely slower than the
  average, and it is the binding one.**
- **The bespoke lemma's target, quantified:** the max Fourier mode amplitude stays **≤ 0.76 < 0.79 (budget ceiling)
  at every measured level** — empirically inside tolerance, median 0.647, but with a **thin margin at the worst level
  (0.76 vs 0.79)**. A rigorous bound must control the level-to-level spikes, not just the median.

**Honest caveat (the quantitative identification is directional, not clean):** the k=3,4 channel *increment* rates
`|A_j(k)/A_{j−1}(k)|` are NOT cleanly measurable — A_j(3), A_j(4) oscillate/change sign, so single-level ratios
scatter wildly (0.04–5.0). Their smoothed two-step rates (0.66–0.73, from HIERARCHY) sit near the max-mode amplitude
0.647, but there is a power-vs-amplitude ambiguity (A_j is linear in the power spectrum) that this probe can't
cleanly resolve. So: **the "max mode is the binding mode" claim is directionally confirmed (max slower than average,
right ballpark, worst-level amplitude inside the budget), not a proven equality.** Hank was right to flag it as the
first thing to confirm, and it half-confirms: the *diagnosis* (wrong mode) is solid; the *exact rate identity* is
consistent but noisy.

## Net for the pen
- **The relaxation bound is bespoke, not off-the-shelf** — the decay literature (DF 1/3, Tao superpoly, BGK) all
  quote the average mode, which is too fast to be an upper bound. Confirmed: average power rate = 1/3 exactly.
- **The route is R66/R74's ν-weighted lifting operator, bounding the MAX/deviation mode** (amplitude ≤ ~0.75). The
  max mode is empirically ≤0.76 at every level — inside the 0.79 budget with a thin margin.
- **The one lemma:** top eigenvalue of the deviation-projected ν-weighted lifting operator ≤ ~0.75, blocked only by
  π-conservation leakage at finite level (clean-basis extraction via K₂'s rank-2 structure = R74's open item).
- DF/Holte supplies the exact transfer skeleton (the scaffold); the ν-perturbation's max mode is the engine.

**Not at stake:** CHANNEL_ID/CARRYLEMMA, HIERARCHY, R1–R30, R80–R82. Hank web-verified Tao 1.17 + DF/Holte
eigenvalues; MAXMODE cheap (cached ρ + build_nu(11), 5.3s).
