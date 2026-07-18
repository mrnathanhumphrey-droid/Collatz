# Chang ↔ trajectory operator factorization — outcome (γ), structural account of why no factorization

**Status.** Decisive negative for explicit factorization. K (R60 v2 size-stratified
Markov kernel, B=109/log_base=1.5, λ_PF=0.951) does **not** project through to
Chang's P (Definition C.5, depth 13 mod-64 Syracuse cylinder kernel, λ_PF=1.0).

| metric | value |
|---|---:|
| Stationary Pearson(π_chang_32, ρ_K) | 0.5404 |
| Kernel Frobenius rel diff K_residue vs P_chang_32 | 0.577 |
| Tensor separability K vs K_r ⊗ K_b rel Frobenius | 0.710 |
| Per-row Pearson mean | 0.787 (range 0.30–0.9999) |

The two operators answer **different dynamical questions** on related but
distinct systems:
- **Chang's P:** unconditioned 1-step kernel under uniform-lift residue
  distribution; stationary = invariant measure of forward Syracuse without
  absorption
- **K:** survivor-conditioned multi-step kernel on (residue, log-size) joint
  state; stationary = quasi-stationary distribution (QSD) under m=1 absorption

These are not algebraically related operators; they describe complementary
aspects of the same Collatz dynamics. The "different Fourier classes"
finding from Result 62 (Chang σ≈0.92 vs trajectory σ=0) reflects this
structural distinction.

## Step 1 — Chang's P verified, π_chang_32 essentially uniform at mod 32

Built P (32×32, mod 64, 128 lifts/residue per Def C.5). Row-stochastic to
machine precision; leading eigvec π_chang_64 with eigval 1.

**Aggregating π_chang_64 to π_chang_32** (sum pairs r mod 64 → r mod 32):

| r mod 32 | π_chang_32 | π / (1/16) |
|---:|---:|---:|
| 1 | 0.063480 | 1.016 |
| 3 | 0.062408 | 0.999 |
| 5 | 0.062990 | 1.008 |
| 7 | 0.062408 | 0.999 |
| ... | ... | ... |
| 27 | 0.061795 | 0.989 |
| 31 | 0.061795 | 0.989 |

**π_chang_32 varies only 0.0618 to 0.0635** (range 2.7%). At mod 32, Chang's
stationary is essentially uniform. The non-uniform structure underlying
Chang's σ≈0.92 lives at deeper modular resolutions (mod 64, 128, ...).

## Step 2-3 — K_residue construction and ρ_K dramatically non-uniform

Built K from 1.5M orbits at N=2^32, B=109/log_base=1.5. λ_PF=0.951.
Marginalized v_PF to residue: ρ_K(r) = ∑_b v_PF(r, b).

| r | ρ_K | π_chang_32 | ρ/π |
|---:|---:|---:|---:|
| 1 | 0.0935 | 0.0635 | 1.47 |
| 3 | 0.0627 | 0.0624 | 1.00 |
| **5** | **0.1090** | 0.0630 | **1.73** |
| 9 | 0.0547 | 0.0623 | 0.88 |
| **13** | **0.0341** | 0.0625 | **0.55** |
| 19 | 0.0405 | 0.0624 | 0.65 |
| 21 | 0.0402 | 0.0625 | 0.64 |
| 23 | 0.0780 | 0.0624 | 1.25 |
| 25 | 0.0390 | 0.0618 | 0.63 |

**ρ_K varies from 0.034 (r=13) to 0.109 (r=5) — factor of 3.2.** This is the
trajectory measure D_avg (after dividing by π_32 it equals D_avg).

Compared to Chang's near-uniform π_chang_32, ρ_K has **massive structural
variation that Chang's stationary completely lacks at the mod 32 level**.

K_residue self-consistency verified: ρ K_residue = λ_PF ρ to machine precision,
Pearson(ρ, ρ K_residue normalized) = 1.000.

## Step 4 — Element-wise comparison: large structural differences

K_residue (16×16, derived from K via v_PF-weighted projection over bins) vs
P_chang_32 (16×16, projected from mod 64).

| comparison | value |
|---|---:|
| Frobenius ‖K_residue − P_chang_stat‖_F | 1.338 |
| Relative Frobenius (vs ‖P_chang‖_F) | **0.577** |
| Max element gap | 0.588 |
| Mean element gap | 0.032 |
| Per-row Pearson mean | 0.787 |
| Per-row Pearson range | [0.300, 0.9999] |
| Stationary Pearson(π_chang_32, ρ_K) | **0.540** |
| Stationary L1 distance | 0.247 |

**The kernels are NOT equal at element level.** Frobenius difference is
58% of P_chang's norm. Some rows match well (Pearson > 0.99) but others
diverge sharply (Pearson 0.30).

## Step 5 — Hypothesis B (separability) rejected

Tested K ≈ K_r ⊗ K_b tensor approximation. Built K_r as the (r → r')
marginalization, K_b as the (b → b') marginalization, and compared their
outer product to the actual K.

| metric | value |
|---|---:|
| RMS difference K vs K_r ⊗ K_b | 0.0132 |
| Relative Frobenius ‖K − K_tensor‖_F / ‖K‖_F | **0.710** |

**71% Frobenius difference — K is highly non-separable.**

SVD spectrum of K_sub (top 10 singular values):

| rank | σ |
|---:|---:|
| 1 | 1.562 |
| 2 | 1.450 |
| 3 | 1.416 |
| 4-7 | ≈ 1.414 (degenerate cluster) |
| 8 | 1.373 |
| 9 | 1.368 |
| 10 | 1.331 |

**K has effectively flat singular spectrum** for the first ~10 modes
(values 1.33–1.56). Rank-1 approximation captures only 8.1% of K's
Frobenius norm. The kernel is not low-rank; it is irreducibly multi-mode.

Hypothesis B (K factorizes as K_r × K_b separable structure) is **clearly
rejected**.

## Step 6 — Survivor-conditioning explains the gap structurally

By construction:
- Chang's P: empirical P[r → r'] = #{first-step images of uniform lifts of r
  hitting r' mod 64} / 128. **Unconditioned, single-step, uniform-lift sample.**
- K: empirical K[(r,b) → (r',b')] from actual orbit transitions during
  trajectories that complete to m=1. **Survivor-conditioned, multi-step,
  natural-density sample.**

The gap between K_residue and P_chang has structural sources:

1. **Survivor conditioning**: K weights each (r, b) cell by visit frequency
   along surviving orbits. P weights each lift uniformly. Survivor-conditioned
   visit weights differ from uniform lifts because some residues serve as
   transient absorbers (small r near m_j) and others as long-lived corridors.

2. **Multi-step composition**: K accumulates statistics over many Syracuse
   steps within each orbit. P is one step. Even if K were unconditioned,
   K_residue would be a mixture of P, P², P³, ... over the orbit length
   distribution.

3. **Size dimension**: K resolves transitions at log-size resolution. The
   v_PF-weighted projection K_residue effectively integrates over size with
   weights specific to the joint stationary, not the uniform-residue-marginal
   weights P uses.

These three effects compound to explain the 0.57 relative Frobenius gap.

**Per-residue stationary divergence pattern** (see Step 4 table): the
structurally important residues are exactly those flagged in prior work:
- r=5 over-enhanced (1.73× Chang) — descent endpoint to m=1
- r=21 under-represented (0.64× Chang) — m_3=21 attractor in inverse tree
- r=13 most depleted (0.55× Chang) — off main descent paths
- r=1 enhanced (1.47×) — m_1=1 vicinity

These match the trajectory-measure mechanism documented in Results 60-65.
Chang's P loses this structure because uniform lifts + first-step kernel +
no size dimension averages it out.

## Step 7 — Character decomposition: different Fourier classes confirmed

Computed DFT of stationary distributions on Z/16Z (16 odd residues mod 32,
indexed by (r-1)/2).

| k | \|F_chang(k)\| | \|F_K(k)\| | ratio K/chang |
|---:|---:|---:|---:|
| 0 | 0.0000 | 0.0000 | — |
| 1 | 0.0016 | **0.1113** | 70× |
| 2 | 0.0024 | 0.0612 | 25× |
| 3 | 0.0005 | 0.0723 | 145× |
| 4 | 0.0005 | 0.0455 | 91× |
| 5 | 0.0013 | 0.0725 | 56× |
| 6 | 0.0031 | **0.0950** | 31× |
| 7 | 0.0019 | **0.1116** | 60× |
| 8 | **0.0034** | 0.0317 | 9× |
| 9..15 | (mirror) | (mirror) | (mirror) |

**Top 5 dominant Fourier modes:**
- Chang: k ∈ {8, 6, 10, 2, 14} — **even modes** (mostly k=8 and conjugates)
- K: k ∈ {7, 9, 1, 15, 6} — **odd modes** (k = ±1, ±7)

The DFT magnitudes for K are 25–145× larger than for Chang at every
non-trivial k. **Chang's stationary at mod 32 is essentially in the principal
character space** (deviations from uniform are ~10^-3); K's stationary lives
in a much richer character space with substantial spectral content at modes
k=±1 (period 16), k=±7 (period 16/7), k=±6 (period 8/3).

This is the explicit Fourier-class distinction:
- **Chang's character class** = even/period-4 dominated; concentrated at k=8
  (which corresponds to residue period 4 = mod 8 structure, consistent with
  Chang's safe classes {1, 3, 7} mod 8)
- **K's character class** = odd-mode dominated; period-16 (mod 32) and
  period-16/7 (related to 3-adic structure of m_j sequence — k=±7 mod 16
  resonates with 3-divisibility patterns)

## Step 8 — Synthesis verdict

### Outcome (γ): No clean operator factorization

The factorization Hypotheses A and B both fail:
- **A (K_residue ≈ P_chang_projected):** stationary Pearson 0.54, kernel
  Frobenius rel diff 0.58 — too large for "approximate equality"
- **B (K ≈ K_r ⊗ K_b):** rank-1 approx captures 8%, separable approx errs
  by 71% in Frobenius norm — clearly not separable

### What the operators DO have in common

- Both are non-negative, irreducibly Markov-like operators on residue spaces
- Both have well-defined leading eigenvectors (Chang's exists at λ=1, ours
  at λ=0.951 due to absorption)
- Both encode genuine Collatz dynamics structure
- Both have non-trivial Fourier content (though at different character classes)

### What the operators do NOT have

- **Same stationary** (Pearson 0.54)
- **Same kernel structure** (Frobenius diff 58% of norm)
- **Algebraic relation** (no projection, marginalization, or factorization
  recovers one from the other)
- **Same Fourier class** (Chang k=8 even, K k=±1, ±7 odd)

### Why the discrepancy is structural, not numerical

Chang's P answers: *"What is the long-run residue distribution under the
unconditioned Syracuse map starting from uniform lifts?"*

K answers: *"What is the joint (residue, log-size) distribution that
survivor-conditioned trajectories spend time at, accounting for the m=1
absorbing boundary?"*

These are mathematically different objects:
- Chang's = invariant measure of forward Syracuse (without absorption)
- Ours = QSD (quasi-stationary distribution) of forward Syracuse with
  absorption + size-tracking

QSD ≠ invariant measure for absorbing Markov chains. They differ by an
exponential factor in the eigenvalue (Chang's λ=1 vs ours λ=0.951). This is
not a defect; it is a fundamental structural distinction in the
mathematical objects.

### For v3.7 / Chang correspondence

**Lead claim (revised):**

> Chang's transfer operator P (Definition C.5) and the trajectory measure
> kernel K (Result 60 v2) describe complementary observables of the same
> Collatz dynamics. Chang's P encodes the unconditioned forward Syracuse
> kernel on residue cylinders; its stationary π is the invariant measure
> with eigenvalue 1. K encodes the survivor-conditioned forward kernel on
> joint (residue, log-size) states; its leading left eigenvector v_PF gives
> the QSD with eigenvalue 0.951.
>
> The two operators are not algebraically related: K_residue (the v_PF-
> weighted projection of K to residues) differs from P at the kernel level
> (relative Frobenius 0.58, per-row Pearson mean 0.79, range [0.30, 0.9999])
> and at the stationary level (Pearson 0.54, L1 distance 0.25 between
> π_chang_32 and ρ_K).
>
> The character decomposition makes the distinction explicit: π_chang_32 is
> essentially uniform with weak even-mode (k=8) structure; ρ_K has dominant
> odd-mode (k=±1, ±7) structure with Fourier coefficients 25–145× larger
> than Chang's at every non-trivial frequency. The two measures live in
> different character classes of L²(Z/32Z).

### Hypothesis C confirmed

The relation between the two frameworks is **conceptual, not algebraic.**
Both probe Collatz dynamics; they ask different questions and yield
different answers. The σ≈0.92 vs σ=0 distinction from Result 62 reflects
this: smooth Fourier decay applies to deeper-mod refinements of Chang's
nearly-uniform stationary, while atomic Fourier concentration applies to
the survivor-conditioned trajectory measure's m_j attractor structure.

This is a clean negative result that **strengthens** the framework synthesis
chapter rather than weakening it: instead of claiming explicit operator
factorization (which would be wrong), the synthesis can frame the two
frameworks as complementary character decompositions probing distinct
dynamical questions.

## What this opens

1. **QSD vs invariant measure literature**: the K vs P distinction is the
   classical absorbing-Markov-chain QSD distinction. Yaglom-Lloyd /
   Méléard-Villemonais theory provides the right framework.
2. **Comparing at deeper mod resolutions**: at mod 32, Chang's structure
   averages out. Comparing K_residue vs P_chang at mod 64 or mod 128 might
   reveal more nuanced relations (Chang's σ≈0.92 lives there).
3. **Multi-step Chang kernel**: P^t for t > 1 might approach K_residue as
   t grows, since K is multi-step composite. Worth computing P^t for t ∈
   {2, 5, 10, 20} and comparing trajectories.

## Files

- `experiments/85_chang_K_factorization.py` — full analysis script
- `experiments_output/chang_K_factorization.csv` — summary metrics
- `experiments_output/chang_K_kernel_diff.csv` — per-element K vs P table
  (256 entries)
- `experiments_output/K_svd_spectrum.csv` — top 10 singular values
- `experiments_output/chang_K_factorization_log.txt` — full diagnostic log

Compute: ~10s (kernel build + projection + SVD + DFT).
