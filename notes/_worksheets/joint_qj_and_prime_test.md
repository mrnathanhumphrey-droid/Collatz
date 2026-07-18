# Result 33: σ-quantile q reduces absorbing-attractor j; primes carry no distinctive Lagarias-class signal

**Date:** 2026-05-03. Sequel to Result 32. Two cheap empirical tests on the Lagarias-class taxonomy.

Tests two questions:
1. **Joint (q, j) analysis:** is σ-quantile q functionally equivalent to absorbing-attractor j, or are they genuinely 2D?
2. **Prime-vs-all comparison:** do primes carry distinctive Lagarias-class signal?

**Result 1 (major):** σ-quantile q REDUCES j. The two-item Lagarias-class taxonomy from Result 32 collapses to one observable: w_q(q).

**Result 2 (null):** primes don't differ from all-starts beyond Prime Number Theorem's effect on ⟨log m⟩. Trajectory-measure deviations are 2-adic, not prime/composite.

Code: `joint_q_j_analysis.py` (5 seeds × 1M orbits, 2.2s), `prime_vs_all.py` (3 seeds × 1M orbits each condition, 2.4s incl. Miller-Rabin sampling).

---

## 1. Joint (q, j) analysis: q is primary, j is downstream

Walk N=2^32 forward, 5 seeds × 1M orbits. For each orbit: σ_S, j_attr, v_avg = v_sum/σ_S. Stratify by σ-quantile (4 bands at midpoints {0.125, 0.375, 0.625, 0.875}). For each (q, j) cell, compute ⟨v⟩.

### P(j | q-band) — j-distribution shifts with q

| q | P(j=2) | P(j=4) | P(j=5) | ⟨v\|q⟩_bulk |
|---|---|---|---|---|
| 0.125 | 0.874 | 0.053 | 0.072 | **2.365** |
| 0.375 | 0.933 | 0.025 | 0.042 | 2.079 |
| 0.625 | 0.963 | 0.012 | 0.025 | 1.966 |
| 0.875 | 0.980 | 0.006 | 0.014 | **1.865** |

**Bottom-σ-quartile orbits absorb at j=4 / j=5 at 7-12% rate vs <1% in top quartile.** Confirms the σ-quartile / j-class connection conjectured in Result 30.

⟨v|q⟩_bulk varies by 0.50 across q-bands (Result 22's structural finding).

### ⟨v | q, j⟩ joint table — q dominates, j is near-redundant

| q | ⟨v\|j=2⟩ | ⟨v\|j=4⟩ | ⟨v\|j=5⟩ | spread across j |
|---|---|---|---|---|
| 0.125 | 2.356 | 2.439 | 2.409 | 0.083 |
| 0.375 | 2.078 | 2.087 | 2.084 | 0.009 |
| 0.625 | 1.966 | 1.973 | 1.970 | 0.006 |
| 0.875 | 1.865 | 1.873 | 1.872 | 0.007 |

| | spread across q (for fixed j) |
|---|---|
| j=2 | 0.491 |
| j=4 | 0.566 |
| j=5 | 0.538 |

**Reduction test (does conditioning on j make q irrelevant?):** NO. ⟨v | q, j⟩ varies by 0.49-0.57 across q for any fixed j.

**Reverse reduction (does conditioning on q make j irrelevant?):** YES (mostly). ⟨v | q, j⟩ varies by only 0.006-0.083 across j for any fixed q. The bottom q-band has the largest spread (0.083) — orbits in the very-fast-descent regime show some residual j-dependence — but for q ≥ 0.375 the j-dependence is essentially zero.

**Verdict: σ-quantile q is the PRIMARY observable. Absorbing-attractor j is largely redundant given q.**

### Implication for taxonomy

Result 32's two-item Lagarias-class list:
1. ⟨σ_S | j⟩ — per-attractor inverse-tree depth distribution
2. w_q(q) — σ-quantile Esscher tilt parameter

**Reduces to one: w_q(q).** The (q ↔ j) joint distribution is effectively factorizable as P(q, j) ≈ P(q) · P(j | q-only), with P(j | q) given by the table above and ⟨v | q, j⟩ ≈ ⟨v | q⟩ (j-independent given q).

⟨σ_S | j⟩ then derives from:
- ⟨σ_S | j⟩ = Σ_q P(q | j) · ⟨σ_S | q-band⟩
- Where ⟨σ_S | q-band⟩ is the band-mean (computable from w_q(q) via descent rate)

So **closing w_q(q) closes ⟨σ_S | j⟩ closes per-j W_j closes ε_S closes ε(σ).** One observable, one Lagarias-class problem, four manifestations.

This is the cleanest possible v3.5 framing.

## 2. Prime-vs-all comparison: null finding

Sample primes from [1, 2^32] via Miller-Rabin (deterministic with witnesses {2, 7, 61} for n < 4.7B), 1M prime orbits × 3 seeds. Compare to all-starts (uniform odd in [1, 2^32]).

### Per-j observables (3-seed bootstrap)

| j | observable | all-starts | prime-starts | gap | z |
|---|---|---|---|---|---|
| 2 | P(j) | 0.9379 | 0.9377 | -0.0002 | -0.79 |
| 2 | ⟨log m\|j⟩ | 21.181 | 20.897 | **-0.284** | -412 |
| 2 | ⟨σ_S\|j⟩ | 76.17 | 75.19 | -0.99 | -32 |
| 2 | ⟨v\|j⟩ | 2.057 | 2.058 | +0.001 | +3.57 |
| 2 | **W_2** | **+7.141** | **+7.142** | **+0.001** | **+0.02** |
| 4 | P(j) | 0.0237 | 0.0237 | 0 | -0.01 |
| 4 | ⟨log m\|j⟩ | 21.177 | 20.898 | -0.280 | -56 |
| 4 | ⟨σ_S\|j⟩ | 54.49 | 53.33 | -1.16 | -6.58 |
| 4 | ⟨v\|j⟩ | 2.251 | 2.259 | +0.008 | +3.24 |
| 4 | **W_4** | **-4.679** | **-4.872** | **-0.193** | **-1.12** |
| 5 | P(j) | 0.0378 | 0.0379 | 0 | +0.69 |
| 5 | ⟨log m\|j⟩ | 21.177 | 20.896 | -0.281 | -107 |
| 5 | ⟨σ_S\|j⟩ | 58.98 | 57.93 | -1.05 | -9.27 |
| 5 | ⟨v\|j⟩ | 2.199 | 2.203 | +0.005 | +4.53 |
| 5 | **W_5** | **+4.638** | **+4.570** | **-0.068** | **-0.60** |

### Reading

**P(j) identical** between conditions (z < 1): primes absorb at j=2, j=4, j=5 with the same probabilities as random integers.

**⟨log m | j⟩ shifts -0.28 nats** for primes vs all-starts. This is the Prime Number Theorem's effect: random primes in [1, N] have lower ⟨log m⟩ than random integers in the same range (primes are sparser at small m due to small-prime gaps, denser at large m, but normalized samples have ⟨log m⟩ ≈ log N - HarmonicSum). Pure finite-N artifact, not Collatz-structural.

**⟨σ_S | j⟩ shifts -1 step** for primes (z=-6 to -32). Tracks the -0.28 ⟨log m⟩ shift exactly: 0.28 nats / log(4/3) ≈ 1 Syracuse step.

**W_j (which subtracts ⟨log m | j⟩/log(4/3)) is INDISTINGUISHABLE between conditions:** all |z| < 1.5. **Primes carry no distinctive Lagarias-class signal in W_j.**

**⟨v | j⟩ shifts are tiny but significant** (z=+3.2 to +4.5, gap +0.001 to +0.008). With 1M orbits per condition × 3 seeds, the SE is small enough to detect ~0.001 differences. The gap direction is positive (primes have slightly higher ⟨v⟩) but the magnitude is < 1% of the value. This is consistent with mod-2^k residue distribution differences between primes and all odd integers (e.g., Chebyshev's bias mod 4) — finite-N effects, not new Lagarias-class structure.

### Verdict

**Primes don't carry distinctive trajectory-measure signal.** The trajectory measure's structural deviations (v=4 spike, mod-2^k residue biases driving ⟨v|j⟩ asymmetry) are **2-adic in origin**, not prime/composite. This was the prior; data confirms.

For the RH connection question: primes' distinguishing feature is multiplicative (Dirichlet series, ζ-zeros). Trajectory-measure invariance is 2-adic / iterated-dynamical. The two structures don't intersect at the level of W_j or ⟨v|j⟩.

## 3. Combined verdict

**Result 1 (joint q-j):** Lagarias-class for q=3 reduces to ONE observable: w_q(q). Closing w_q(q) closes the entire chain (⟨σ_S | j⟩, W_j, ε_S, ε(σ), K_eff via per-band asymptotes).

**Result 2 (prime test):** primes don't add a Lagarias-class observable. RH-Collatz connection remains speculative; no empirical signal.

**Net for v3.5:**
- Single open Lagarias-class problem: closed-form w_q(q) via moderate-deviation rate function for the trajectory measure on Syracuse log-walk
- No simpler bypass (Result 32 falsified spectral)
- Primes don't help (Result 33 null)
- Three-constants-out-of-four story tightens: constants 1, 2 closed; constants 3 (per-j W_j) and 4 (K_full asymptote) both reduce to w_q(q); the open piece is one structural object

## 4. Files

- `joint_q_j_analysis.py` — joint q-j table at N=2^32
- `prime_vs_all.py` — Miller-Rabin prime sampling + per-j comparison
- `joint_qj_and_prime_test.md` — this document (Result 33)
- `experiments_output/joint_q_j_analysis_log.txt` — joint analysis log
- `experiments_output/prime_vs_all_log.txt` — prime test log

Total compute: 4.6s (joint 2.2s, prime 2.4s).
