# RESULT — PROBE MEAN1: the mean-1 constraint, sharpened to exact class-means (2026-07-26)

**Probe:** `probe_mean1.py`. Wilson's constraint: since 4 generates G_r, `Σ_k γ_r(k)=3^r` ⟹ `mean_k γ_r(k)=1` exactly.
White is the pinned mean, unconditionally.

## Verified exact — and the class-means are the level-1 seed, exactly

At r=14,15,16: `mean_k γ_r(k) = 1.0000000000` (machine exact). The M±/split is **cleaner than the estimate**:

| | measured (r=14,15,16) | value |
|---|---|---|
| M₊ = mean over 3∣k | 1.66667 | **= 5/3 exactly** |
| M₋ = mean over 3∤k | 0.66667 | **= 2/3 exactly** |
| (1/3)M₊+(2/3)M₋ | 1.00000000 | = 1 ✓ |

Not M₊≈1.5 / M₋≈0.75 — it's **M₊=5/3, M₋=2/3, the level-1 seed values, pinned at every level.** (M₊ excluding the
divergent k=0=X_r term is also 5/3 — the k=0 channel is 1 of ~14M, negligible.)

## Why (two lines) — the mod-3 marginal is frozen

`Σ_{3∣k} C_r(k) = Σ_s ρ(s)·(mass of ρ in class s mod 3) = Σ_{j} f_j²`, where `f_j` = mass of the dlog profile in
residue class j mod 3. **The tower property `ν_r mod 3 = ν_1` freezes `f=(0,1/3,2/3)` at every level** ⟹
`Σ_{3∣k}C_r(k)=5/9` exactly ⟹ `M₊=3·(5/9)=5/3`, `M₋=(3/2)·(4/9)=2/3`. This is the lag-domain dual of R66's
"primitive Fourier sum invariant." So Wilson's mean-1 upgrades to: **the enriched and depleted class-means are exactly
the level-1 collision factors, at every level** — a two-line theorem, and the sentence that anchors the dichotomy
(the channels distribute around a pinned mean whose two class-halves are 5/3 and 2/3).

Wilson's small-k enriched sample (1.237,1.372,1.528,2.112, avg 1.56) undershoots the true M₊=5/3 because low-v₃
channels have relaxed down; the full population averages exactly 5/3.

## Surfaced: R66 and R74 (the operators Wilson asked for) — and the tension in them

**R66 (`notes/_worksheets/decay_law_derivation.md`):** the forward-Syracuse Markov chain on `(ℤ/3^k)*` (coprime-to-3
residues). State m mod 3^k; kernel `K[r→s]=P(T(m)≡s)`, `T(m)=(3r+1)·2^{−v} mod 3^k`, v~Geom(½), 2^{−v} period
`ord_{3^k}(2)=2·3^{k−1}`. Explicit stationary π_r (k=1:(1/3,2/3); k=2:(8,16,11,4,2,22)/63). Fourier: primitive sum
`Σ_a|μ̂(a/3^k)|²≈7/15` invariant; **average `|μ̂|²` decays 3^{−(k−1)} (rate 1/3); MAX decays ~2^{−k} but "SLOWING."**

**R74 (`notes/_worksheets/lifting_operator_spectral.md`):** the lifting operator `L_k: π_k→π_{k+1}` (from
level-(k+1) stationary's conditional sub-cell masses α,β,γ on each level-k cell). PROVEN algebraic Parseval recursion
`S_{k+1}=3^{k+1}‖d_{k+1}‖²`, d = sub-cell deviation from white. Deviation-restricted `M_k=D_{k+1}L_k D_k` (D removes
the π component); finite-level top SVs **0.6901→0.7560→0.8188→0.8850 GROW with k** = π-conservation leakage (retains
spurious π-correlated components). Clean rate would be 1/√3≈0.577; extraction needs "careful basis selection,
likely explicit via K₂'s rank-2 structure."

**THE TENSION (flag before the lemma):**
- **R74's operator governs the AGGREGATE deviation `‖d‖²` — its clean rate is 1/√3 ≈ 0.577 (the AVERAGE mode).**
- **Hank's binding lemma is about the MAX single mode (~0.707, R66's max), which is SLOWER than the aggregate.**
- ⟹ **R74's operator, even with the π-leakage removed, gives 0.577 — TOO FAST to bound the binding max mode (0.707).**
  A single channel k=3,4 relaxes at the slowest coupled mode (max), not the aggregate average. So bounding R74's M_k
  eigenvalue does NOT bound the persistence-binding channel.
- **The lemma actually needs a bound on R66's MAX primitive coefficient `max_a|μ̂(a/3^k)|²`** — and R66 reports that mode
  as ~2^{−k} **"but slowing,"** which is the worry: MAXMODE (probe, k≤16) shows it plateauing at amplitude ~0.647
  median with spikes to 0.76 (vs the 0.79 budget ceiling). The "slowing" appears to saturate below the ceiling, but
  the object is R66's max mode, not R74's aggregate, and the margin is 4%.

**Bottom line for the pen:** R66/R74 are the right files, but they split the work — R74 proves the aggregate recursion
(rate 0.577, the average), R66 has the explicit chain whose MAX coefficient is the binding object (~0.707, slowing).
The lemma is a bound on **R66's max primitive Fourier coefficient** (using R66's explicit stationary π_r), not R74's
aggregate lifting operator. R74's "K₂ rank-2 / careful basis" open item is about extracting the aggregate 7/45
coefficient — a *different* target from the max-mode bound. Don't write the lemma against M_k (wrong mode).

**Not at stake:** HIERARCHY, CHANNEL_ID, R1–R30. MEAN1 cheap (FFT autocorrelation, 6.5s).
