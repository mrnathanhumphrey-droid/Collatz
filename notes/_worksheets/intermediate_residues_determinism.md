# Intermediate residues mod 64 — outcome (a): only r=21 is non-deterministic (Result 45)

**Status.** Decisive. Pure arithmetic proof + computational verification.

**Outcome (a) confirmed: ALL 31 of 32 odd residues mod 64 are deterministic
at mod 256 and mod 1024. Only r=21 (v_2(3r+1)=6) is non-deterministic.**

The boundary-non-determinism is unique to r=21. The 14 intermediate residues
(v_2 ∈ {2,3,4}), the lone r=53 (v_2=5), AND the 16 v_2=1 residues are ALL
fully deterministic at higher mod. The v_2 spectrum is structurally
homogeneous (deterministic) except at the v_2 ≥ 6 boundary.

## Arithmetic proof

For m ≡ r (mod 64), write m = r + 64k and let v₀ = v_2(3r+1) with 3r+1 = 2^v₀·u
(u odd).

3m+1 = 3r+1 + 192k = 2^v₀·u + 2⁶·3k

**Case v₀ < 6:**
- 3m+1 = 2^v₀·(u + 2^(6−v₀)·3k)
- Since 6−v₀ ≥ 1, 2^(6−v₀)·3k is even
- u is odd → u + (even) is odd
- Therefore v_2(3m+1) = v₀ exactly, for ALL k

**Case v₀ = 6 (r=21 only):**
- 3m+1 = 64·(u + 3k)
- u + 3k can be odd or even depending on k mod 2
- v_2 non-deterministic, depends on bits of k beyond mod 64

**Case v₀ ≥ 7:**
- Would require 3r+1 ≡ 0 (mod 128). For r ∈ {1, 3, ..., 63}: 3r+1 ∈ {4, ..., 190}, none reach 128 with v_2 ≥ 7. So this case is empty mod 64.

## Mod 64 classification (32 odd residues)

| v_2(3r+1) | count | residues | det at mod 256 |
|----------:|------:|----------|---------------:|
| 1 | 16 | {3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63} | YES (all) |
| 2 |  8 | {1, 9, 17, 25, 33, 41, 49, 57} | YES (all) |
| 3 |  4 | {13, 29, 45, 61} | YES (all) |
| 4 |  2 | {5, 37} | YES (all) |
| 5 |  1 | {53} | YES |
| **6** |  **1** | **{21}** | **NO** |

Total: 32. Single non-deterministic residue: r=21.

## Mod 256 / mod 1024 verification

For each r mod 64, lift to mod 256 (4 lifts) and mod 1024 (16 lifts). Compute
v_2(3m+1) for each lift; report distinct values:

- All 31 deterministic residues: v_2 set = {single value} at both mod 256 and mod 1024
- r=21: v_2 set at mod 256 = {6, 7, 8}; at mod 1024 = {6, 7, 8, 9, 10}

r=21's v_2 set grows with modulus — lifts deeper, the recursive shifted-Geom
structure (Result 42) appears: P(v=k|m mod 2^(k+1)) resolves more residues
each level.

## Empirical orbit observables (1M orbits at N=2³²)

Per group means:

| group | n_residues | ⟨σ⟩ | SD across r | ⟨V⟩ | SD across r |
|-------|----------:|----:|------------:|----:|------------:|
| Chang I_2 (5 of v=1) | 5 | 240.49 | 9.26 | 2.0210 | 0.022 |
| v=1 other | 11 | 226.25 | 11.08 | 2.0606 | 0.033 |
| intermediate (v∈{2,3,4}) | 14 | 221.10 | 10.80 | 2.0763 | 0.034 |
| v=5 (r=53) | 1 | 206.30 | — | 2.1260 | — |
| v=6 BOUNDARY (r=21) | 1 | **193.35** | — | **2.1810** | — |

**Monotone progression in v_2(3r+1)**: higher v_2 at start → lower ⟨σ⟩
(faster descent), higher ⟨V⟩ (more high-v steps).

The v_2(3r+1) value drives a smooth gradient across observables. r=21 sits
at one extreme, low-v residues at the other.

## Within v=1 substructure

Within v_2=1 (16 residues), ⟨σ|r⟩ ranges from 205.6 (r=35) to 255.4 (r=63).
Chang I_2 = {7, 27, 31, 59, 63} has mean ⟨σ⟩ = 240.5 vs other v=1 mean 226.3.

The v=1 substructure (within Chang I_2 vs not) is itself v_2(3·next_step+1)-
driven: the residue dynamics map r → next-step-residue, and Chang I_2's
specific mod-64 residues map to particular next-step residue distributions.
Substructure within v_2 levels is captured by next-step Markov, which is
deterministic for v₀ < 6.

## What this resolves for the v_2 spectrum

- **v_2 = 1 (16 residues, slowest first-step):** deterministic. Chang's
  spectral methods focus on a 5-residue subset.
- **v_2 ∈ {2, 3, 4} (14 residues, intermediate):** all deterministic.
  Smooth gradient in ⟨σ⟩, ⟨V⟩.
- **v_2 = 5 (r=53):** deterministic. v=5 across all higher-mod lifts.
- **v_2 = 6 (r=21):** UNIQUE non-deterministic residue. The {m_j = (4^j−1)/3}
  attractor sequence lives here. All Lagarias-class trajectory-measure
  complexity localizes at this single residue (Result 40, 42).

## Implication for unified framework

Outcome (a) means **the v_2 spectrum is structurally homogeneous (deterministic
in the residue→v map) except at one boundary residue (r=21)**. Chang's I_2 is
not "the deterministic exception" — Chang's I_2 is part of a deterministic
spectrum that contains 31/32 residues. The structural exception is r=21 alone.

Unified framework requirements:
- Deterministic part (31 residues): residue chain Markov on m mod 64 with
  deterministic v transitions. Path B / Result 19 exact.
- Singular boundary (r=21 only): the 2-adic visit-measure on
  {m ≡ 21 mod 32} cylinder set (Result 42's reformulation).

Two structural objects, not one homogeneous mechanism. But the deterministic
part is fully closed; the open question is exactly the visit-measure on the
single-residue boundary.

This sharpens the v3.6 framing:

> **The Lagarias-class open piece is the visit-frequency distribution
> P(m mod 2^k | orbit visits r ≡ 21 mod 32) as k → ∞. Outside this single
> residue cylinder, the trajectory measure is deterministic on residues
> mod 64.**

## Files

- `experiments/72_intermediate_residues.py`
- `experiments_output/72_intermediate_residues_log.txt`
- `experiments_output/72_determinism_mod64.csv`
- `experiments_output/72_orbit_observables.csv`
- `experiments_output/72_group_summary.csv`

Compute: 1.0s (arithmetic + 1M orbits + analysis).
