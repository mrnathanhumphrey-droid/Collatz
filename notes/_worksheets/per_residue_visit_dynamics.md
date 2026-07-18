# Per-residue-visit dynamics — all complexity localizes at r=21 mod 32 (Result 39)

**Status.** Decisive structural finding combining outcomes (b) and (c):

- **(c)+** v is **deterministic** mod 32 for 15 of 16 odd residues. Not just
  "deviates from Geom(1/2)" — completely fixed by residue alone.
- **(b)** Visit-number and position dependence exist, but ONLY at r=21
  (the m_j sub-stratum). All hidden-state behavior is concentrated there.

This means the residue-chain Markov framework (Path B / Result 19) is
**structurally exact for 15/16 residues**. The Lagarias-class trajectory
measure question reduces precisely to the dynamics at r=21 mod 32, which
is the m_j = (4^j−1)/3 sub-stratum requiring deeper 2-adic refinement.

## Per-residue P(v) at N=2³⁶ (200K orbits, 16.9M visits)

| r mod 32 | n_visits | ⟨v⟩ | P(v=⟨v⟩) | mathematical reason |
|---:|---------:|------:|----------:|---------------------|
| 1 | 0.997M | 2.000 | 1.000 | 3·1+1 = 4 = 2² |
| 3 | 1.041M | 1.000 | 1.000 | 3·3+1 = 10 = 2·5 |
| **5** | 1.260M | **4.000** | 1.000 | 3·5+1 = 16 = 2⁴ |
| 7 | 1.102M | 1.000 | 1.000 | 22 = 2·11 |
| 9 | 1.042M | 2.000 | 1.000 | 28 = 2²·7 |
| 11 | 1.054M | 1.000 | 1.000 | 34 = 2·17 |
| 13 | 0.980M | 3.000 | 1.000 | 40 = 2³·5 |
| 15 | 1.098M | 1.000 | 1.000 | 46 = 2·23 |
| 17 | 1.130M | 2.000 | 1.000 | 52 = 2²·13 |
| 19 | 0.961M | 1.000 | 1.000 | 58 = 2·29 |
| **21** | 0.960M | **5.924** | **varies** | 3·21+1 = 64 = 2⁶ → v ≥ 5, mod 32 insufficient |
| 23 | 1.104M | 1.000 | 1.000 | 70 = 2·35 |
| 25 | 0.967M | 2.000 | 1.000 | 76 = 2²·19 |
| 27 | 1.093M | 1.000 | 1.000 | 82 = 2·41 |
| 29 | 1.127M | 3.000 | 1.000 | 88 = 2³·11 |
| 31 | 1.011M | 1.000 | 1.000 | 94 = 2·47 |

**For 15/16 residues, v is exactly determined by r mod 32.** No probabilistic
content. The marginal Geom(1/2) baseline P(v=k) = 2⁻ᵏ emerges only after
averaging over residues with the natural visit weights.

## Why r=21 is special

m_j = (4^j − 1)/3: m_1=1, m_2=5, m_3=21, m_4=85, m_5=341, m_6=1365, ...
- m_3 = 21 ≡ 21 (mod 32)
- m_4 = 85 ≡ 21 (mod 32)
- m_5 = 341 ≡ 21 (mod 32)
- m_j ≡ 21 (mod 32) for ALL j ≥ 3

Reason: m_j = (4^j − 1)/3 = (4·4^(j-1) − 1)/3. For j ≥ 3, 4^(j-1) ≡ 0 (mod 16),
so 4·4^(j-1) ≡ 0 (mod 64), hence m_j ≡ −1/3 (mod 64) ≡ 21 (mod 32).

For m ≡ 21 (mod 32), 3m+1 ≡ 64 (mod 96), so v_2(3m+1) ≥ 5 (depending on
higher bits of m). Mod-32 cannot determine v; the m_j sub-stratum dynamics
require mod-2^k for arbitrary k.

## Test 2: visit-number drift (Markov check)

Across 16 residues, only r=21 shows non-zero drift v_1 → v_5:

| r | ⟨v⟩@v=1 | ⟨v⟩@v=5 | drift |
|---:|------:|------:|------:|
| 1–19, 23–31 | (deterministic) | (deterministic) | 0.000 |
| **21** | 5.999 | 5.827 | **−0.173** |

Aggregate: mean drift across all residues = −0.011, SD = 0.042.
Max |drift| = 0.173 (r=21).

**For 15/16 residues, dynamics are Markov on r mod 32 trivially** (since
v is deterministic). For r=21, visits are NOT visit-independent — later
visits have systematically smaller v (mean drops 0.17 from visit 1 to 5).

## Test 3: position-within-orbit drift

Same pattern. Only r=21 shows non-zero position drift:
- Early in orbit ([0, 0.25]): ⟨v|r=21⟩ = 5.937
- Late in orbit ([0.75, 1.0]): ⟨v|r=21⟩ = 5.866
- drift = −0.070

Other residues: drift = 0.000 by construction (deterministic).

## Test 4: autocorrelation across visits to same residue

For 15/16 residues: Cov(v_i, v_{i+1} | same r, same orbit) = 0 trivially
(both v_i and v_{i+1} are constants).

For r=21: Cov = −0.023, corr = **−0.013**. Small negative autocorrelation —
consecutive visits to r=21 within same orbit have slightly anti-correlated
v values.

## Implications for v3.6

### Path B Markov framework is structurally exact for 15/16 residues

For r ∈ {1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31}:
v_t = f(m_t mod 32) is a deterministic function. The residue chain
m_t → m_{t+1} mod 32 is a deterministic map — fully captured by Path B
framework with no hidden state.

### All Lagarias-class complexity = r=21 sub-stratum dynamics

The single residue r=21 mod 32 carries:
- Variable v (Geom-like distribution, ⟨v⟩=5.92)
- Visit-number dependence (drift −0.17 across visits 1→5)
- Position dependence (drift −0.07 early→late)
- Non-zero autocorrelation (−0.013)
- All m_j ≡ 21 (mod 32) for j ≥ 3 — the entire infinite sub-stratum

**The Lagarias-class question reduces precisely to: characterize the
distribution of v_2(3m+1) for m ≡ 21 (mod 32), conditional on the orbit
context (visit number, position within orbit, prior visits).**

This is the sharpest possible localization of the open question. Item A
in Result 34's Lagarias-class catalog (per-j W_j → ⟨σ_S|j⟩) IS precisely
this: the m_j class IS the r=21 mod 32 class (for j ≥ 3, mod 32; for
finer j-discrimination, mod 2^k for larger k).

### r=5 is the m_2 boundary case

r=5 mod 32 has v=4 deterministically (since 3·5+1 = 16 = 2⁴, with no
higher-bit dependence). m_2 = 5 IS this class. Higher m ≡ 5 (mod 32) like
m=37, 69, 101, ... all have 3m+1 ≡ 16 (mod 32 trivially? No — 3·37+1=112=2⁴·7
so v=4; 3·69+1=208=2⁴·13 so v=4). For m ≡ 5 (mod 32), 3m+1 has v=4 always
because 3·5+1=16 and 3·32=96 ≡ 0 (mod 16), so adding 96·k to 16 keeps
v_2 = 4.

So m_2=5 doesn't lift to a sub-stratum the way m_j for j≥3 does. The
Lagarias-class structure starts at j=3 in mod-32 resolution.

### Connection to Result 32 and 33

Result 32 Route B closed constant 4 bulk via Esscher per-step + algebraic
Cov[T,V]. The "Esscher per-step" framework averages over all residues with
their visit-frequency weights, producing the Geom(1/2)-like marginal
P(v=k) = 2⁻ᵏ (since the deterministic-v residues collectively give exactly
this marginal under their occupancy frequencies). The trajectory-measure
"deviation" beyond Geom(1/2) is entirely the r=21 contribution.

For per-band stratification: orbits in different σ-bands visit r=21 with
different frequencies, AND the conditional ⟨v|r=21, band⟩ may differ across
bands. Both effects feed into Result 25's per-band E_band ≠ 2 deviations.

## Files

- `experiments/70_per_residue_visit.py`
- `experiments_output/70_per_residue_visit_log.txt`
- `experiments_output/70_test1_p_v_given_r.csv`,
  `70_test2_visit_number.csv`, `70_test3_position.csv`,
  `70_test4_autocorr.csv`
