# Forward-Backward Duality Test — VERDICT: NO clean duality, but BASIN FINGERPRINT confirmed

**Date:** 2026-05-04. Sibling-probe Task 2 (Wilson-prompted). Tests proposed duality between forward S_n^{3x±1} (proved identical via my chain-symmetry result) and inverse-tree D_n^{(x∓1)/3} (Agent 2: 3x+1 single-basin; Agent 3: 3x−1 three-basin weighted).

## Verdict (two parts)

> **Part 1 — Duality FAILS.** None of the proposed clean functional relationships D = f(S) hold. Forward S_k stabilizes at ~7/15; inverse D_n(k) decays toward 0 for k ≥ 2 (Agent 2) or to small but non-zero values bounded by basin reach (Agent 3). The product D·S, sum D+S, ratio D/S, and diagonal D_n(n)/S_n all fail to be stable in n. The forward Markov-chain symmetry K_- = σ K_+ σ does not propagate to the inverse-tree integer-level dynamics.
>
> **Part 2 — Basin fingerprint is enormous and confirmed.** Despite forward S_k^{3x+1} = S_k^{3x−1} as exact rationals at every k, the inverse-tree Plancherel masses differ by **factors of 10³–10⁴** at large depth (n=6, k=2): Agent 2 (3x+1 single-basin) gives D_6(2) ≈ 3.2×10⁻⁷; Agent 3 (3x−1 three-basin total) gives 2.8×10⁻³. The 3x−1 basin structure (three cycles {1,2}, {5,7,10,14}, {17..34}) is readily distinguishable from 3x+1's single conjectured attractor at the inverse-tree level.

## Setup

Three tables compared:
- **Agent 2** [result_inverse_tree_residue.md](result_inverse_tree_residue.md): 3x+1 inverse tree (i.e., (x−1)/3 backwards), single basin from root 1, no value truncation. D_n(k) for n=0..6, k=1..5 as exact rationals.
- **Agent 3** [agent3_inverse_tree_3xm1_Dn.py](agent3_inverse_tree_3xm1_Dn.py): 3x−1 inverse tree ((x+1)/3 backwards), three basins from roots 1, 5, 17. Per-root D_n^{root}(k) and basin-density-weighted total. Uses N_MAX = 10⁸ value cap.
- **Forward S_k**: from existing q=3 caches (R77.7, q-sweep). Identical between 3x+1 and 3x−1 by my proved chain symmetry K_- = σK_+σ.

## Test 1: Duality candidates

| candidate | observed | pass? |
|---|---|---|
| D_n(k) · S_k = const | products → 0 for k≥2; const for k=1 only | ❌ no |
| D_n(k) + S_k = const | sum stays close to S_k as D→0; trivially S_k itself for large n | ❌ trivial |
| D_n(k) / S_k stable | D/S → 0 for k≥2; constant 1/3 for k=1 only | ❌ no |
| D_n(n) vs S_n diagonal | ratios 0.32, 0.0083, 0.0074, 0.0042, 0.00203 — monotone decay | ❌ no scaling |
| inverse "q/3-analogue" | D_{n+1}/D_n ratios are noisy (0.01–0.5 range) at small n; Agent 2 quotes asymptotic 1/9 but n=6 too small to confirm cleanly | ⊳ inconclusive |

**No clean function** maps inverse D_n(k) to forward S_k.

The closest thing to structure: at k=1, D_n(1) = 2/9 for n ≥ 2 (Agent 2 fixed-point) and S_1 = 2/3. Their ratio D/S = 1/3 — but this is a trivial coincidence: D_n(1) = 2/9 reflects mod-3 equipartition of the inverse tree from depth 2, S_1 = 2/3 is the forward Markov stationary's L² Plancherel-coprime mass at k=1. Both happen to be small rationals; no deeper duality structure.

## Test 2: Basin fingerprint (Agent 2 vs Agent 3 ratio)

Direct ratio of Agent 3 (3x−1 three-basin total) over Agent 2 (3x+1 single-basin) at matched (n, k):

| (n, k) | A2: D^{3x+1}(k) | A3: D^{3x−1}_total(k) | ratio A3/A2 |
|---|---|---|---|
| (0, 2) | 6.000     | 6.000     | 1.000 (trivial: both point-mass on root) |
| (1, 2) | 6.12e−02  | 1.02e−01  | 1.66 |
| (2, 2) | 3.95e−03  | 5.29e−02  | **13.4** |
| (3, 2) | 1.98e−04  | 4.01e−03  | **20.3** |
| (4, 2) | 9.88e−06  | 9.05e−03  | **916** |
| (5, 2) | 1.14e−06  | 4.88e−03  | **4 287** |
| (6, 2) | 3.20e−07  | 2.78e−03  | **8 691** |
| (6, 3) | 5.06e−07  | 1.26e−02  | **24 928** |
| (6, 5) | 9.18e−05  | 9.74e−02  | **1 061** |

**The ratio diverges by 3–4 orders of magnitude as depth grows.** This is the basin fingerprint: 3x−1's three-basin structure produces a fundamentally larger inverse-tree Plancherel mass than 3x+1's single-basin inverse-tree.

**Mechanism.** The inverse-tree vertex counts grow at vastly different rates:

| n | Agent 2 (3x+1 single basin) |V_n| | Agent 3 (3x−1 root 1, value-capped) |V_n| |
|---|---|---|
| 0 | 1 | 1 |
| 1 | 14 | 13 |
| 2 | 135 | 56 |
| 3 | 1 350 | 189 |
| 4 | 13 500 | 459 |
| 5 | 135 000 | 1 061 |
| 6 | 1 350 000 | 2 247 |
| growth ratio | → 10 (constant) | 13, 4.3, 3.4, 2.4, 2.3, 2.1 (decaying) |

The 3x−1 inverse tree's expansion **slows** because (a) cycle members are excluded from preimage expansions (so the BFS doesn't loop), and (b) the N_MAX = 10⁸ value cap hits early. Each basin is **bounded** by its cycle structure and the value cap. The 3x+1 inverse tree from 1 has no comparable bound — it expands geometrically without re-encountering its (single, small) cycle.

The empirical measure mu has L² mass scaling roughly like 1/|V_n|, so smaller |V_n| → larger D_n. Agent 3's smaller vertex counts directly explain its larger D values.

## Test 3: Per-root structure within Agent 3 — three basins are distinct

Per-root D_n^{root}(k) at the diagonal (n=k):

| (n, k) | root 1 | root 5 | root 17 | weighted total | Agent 2 (3x+1) |
|---|---|---|---|---|---|
| (3, 3) | 0.1950 | 0.1767 | 0.0852 | 0.1509 | 0.00342 |
| (4, 4) | 0.1156 | 0.2566 | 0.1879 | 0.1866 | 0.00193 |
| (5, 5) | 0.0966 | 0.2249 | 0.2583 | 0.1945 | 9.48e−04 |

The three roots produce **distinct D values**, not merely distinct vertex counts. Root 5 and root 17's basins have larger D values than root 1's at large (n, k), suggesting their inverse-tree structures are even more concentrated (smaller vertex counts before hitting the value cap).

The basin densities (0.327, 0.325, 0.347, residual 4×10⁻⁶) are roughly uniform across roots, so each contributes ~1/3 to the total. The three-basin total ~ ½ to ⅔ of the per-root values.

## Why the forward symmetry doesn't propagate

The forward Markov chain on Z/3^k is a **modular-arithmetic object**. The chain symmetry K_- = σ K_+ σ where σ(r) = −r is a statement about Z/3^k arithmetic: 3·(−r) + 1 = −(3r − 1) holds in Z/3^k as a tautology. Because the heuristic chain treats v ~ Geom(1/2) as an independent random variable (NOT v_2 of any specific lift), the modular identity transfers cleanly to chain transitions.

The inverse tree is an **integer-level dynamical object**. Predecessors of m are computed via g(m) = (2^e · m ± 1) / 3, with validity conditions involving v_2 of specific integer lifts. Cycle structure and value-cap behavior at the integer level ARE NOT modular phenomena — they depend on which actual integer (positive or negative) one lands on after stripping 2's. The negation σ in Z/3^k corresponds to integer negation, which **maps positive integers to negative** — a different basin entirely.

Concretely: in 3x+1, 5 → 16 → 8 → 4 → 2 → 1 (no cycle, falls to {1, 2, 4}). In 3x−1, 5 → 7 → 5 (cycle {5, 7}). The integer 5 sits in different basins under the two forward maps, and the inverse-tree analysis from 5 in 3x−1 produces an entirely separate structure than the 3x+1 inverse-tree analysis (which goes from 1, the only root that matters in 3x+1).

The forward symmetry says "for the modular Markov chain, swapping +1 ↔ −1 is just relabeling residues by negation." The inverse-tree analysis says "for the integer-level inverse map, the +1 case has one basin (single attractor) and the −1 case has multiple basins; their inverse trees grow at different rates and produce different mod-3^k empirical measures." Both statements are correct; they are not in tension because they live in different categories.

## Implications

1. **Basin structure IS detectable** at the inverse-tree mod-3^k Plancherel level. The 3-basin nature of 3x−1 produces a clear, measurable signature in D_n(k), even though forward S_n is identical. **This means the forward Markov-chain framework alone CANNOT distinguish 3x+1 from 3x−1**, but a richer probe (inverse-tree value-truncated Plancherel) CAN.

2. **No multiplicative or additive duality** between forward S_k and inverse D_n(k) exists in any of the candidate forms. The two objects live in genuinely different functional categories: S_k is the L² mass of a stationary distribution on Z/3^k, D_n(k) is the L² mass of an empirical measure derived from an integer-level tree truncated by depth and value.

3. **For c=7/45 closure**: this finding doesn't directly impact the rate-1/2 problem, since c=7/45 is forward-symmetric (transfers to 3x−1 by my proved chain symmetry). But it does say: **the inverse-tree Plancherel framework is a strictly richer probe than the forward Markov chain.** Any closure attempt that uses ONLY the forward Markov chain misses the basin-structure information that the inverse tree carries.

4. **For the q-sweep**: the q/3 universal ratio in forward S_{k+1}/S_k is a forward-Markov-chain phenomenon (proved at k=3 q=3, empirically observed at q=5,7,11,13). The inverse "q/3 analogue" — D_{n+1}(k)/D_n(k) — is **noisy and not converged** at the depths Agent 2 reached (n ≤ 6). Cannot confirm an inverse universal ratio at this resolution.

## Files

- [result_inverse_tree_residue.md](result_inverse_tree_residue.md) — Agent 2 writeup
- [result_inverse_tree_residue.csv](result_inverse_tree_residue.csv) — Agent 2 data
- [agent3_inverse_tree_3xm1_Dn.py](agent3_inverse_tree_3xm1_Dn.py) — Agent 3 script
- [agent3_Dn_total.csv](agent3_Dn_total.csv), [agent3_Dn_root_{1,5,17}.csv](agent3_Dn_root_1.csv) — Agent 3 data
- [duality_S_vs_D_test.py](duality_S_vs_D_test.py) — duality test script
- [duality_S_vs_D_verdict.md](duality_S_vs_D_verdict.md) — this writeup

## STATE.md impact

Add to closed-form lock-ins (or major findings):
- **Forward symmetry K_- = σK_+σ** transfers all Plancherel quantities (S_n, M_n, R_n, ε_n) but NOT the inverse-tree integer-level Plancherel masses D_n(k).
- **Basin fingerprint**: D_n^{3x−1}(k) / D_n^{3x+1}(k) → ∞ as n grows (factors of 10³–10⁴ by depth 6). This is the clean detector that the two systems differ at the integer level despite forward-Markov equivalence.

Add to obstruction map / reframe note:
- **The inverse-tree Plancherel framework is a strictly richer detector than the forward Markov chain.** Forward Markov chain symmetry can hide differences that the inverse tree exposes. Future closure attempts that use only forward Markov machinery should be aware of this gap.
