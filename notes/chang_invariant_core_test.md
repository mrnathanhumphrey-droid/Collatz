# Result 41: Chang/Quadrium I_2 and our r=21 mod 32 are INDEPENDENT observables — opposite v_2 extremes

**Date:** 2026-05-03. Empirical test of structural connection between Chang/Quadrium 2603.11066v6's I_2 = {7, 27, 31, 59, 63} mod 64 and our r=21 mod 32 sub-stratum localization (Result 17/19/40).

**Verdict: outcome (b) PARTIAL CONNECTION. Both are residue-selected via v_2(3r+1) arithmetic, but at OPPOSITE ENDS of the spectrum.** Chang's I_2 picks v_2=1 (slowest descent residues); our r=21 picks v_2=6 (fastest descent, boundary). Empirical orbit observables (⟨σ⟩, ⟨V⟩) deviate from generic in OPPOSITE directions.

Code: `chang_invariant_core_test.py`. Compute: ~1.0s.

---

## 1. Step 1: Mod-64 arithmetic

v_2(3r+1) distribution across 32 odd residues mod 64:

| v_2 | residues | count |
|---|---|---|
| 1 | r ≡ 3 mod 4: {3, **7**, 11, 15, 19, 23, **27**, **31**, 35, 39, 43, 47, 51, 55, **59**, **63**} | 16 |
| 2 | {1, 9, 17, 25, 33, 41, 49, 57} | 8 |
| 3 | {13, 29, 45, 61} | 4 |
| 4 | {5, 37} | 2 |
| 5 | {**53**} | 1 |
| 6 | {**21**} | 1 (BOUNDARY) |

**Chang's I_2 ⊂ v_2=1 set (slowest descent).** All 5 elements are r ≡ 3 mod 4 with 3r+1 = 2·(odd).
**Our singular {21, 53} = top-v_2 residues (fastest descent).** v_2(3·21+1) = 6 (boundary, non-deterministic at higher mod), v_2(3·53+1) = 5.

These sit at **opposite extremes of the v_2(3r+1) spectrum**.

## 2. Step 2: Determinism at mod 256

For each r mod 64, test whether v_2(3m+1) is constant across the four lifts m ≡ r, r+64, r+128, r+192 mod 256:

| Residue set | Non-deterministic at mod 256 |
|---|---|
| Chang's I_2 = {7, 27, 31, 59, 63} | **None** — all v_2=1 stable across lifts |
| Our singular = {21, 53} | r=21 yes (v_2 ∈ {6, 7, 8, 9}); r=53 no (v_2=5 stable) |
| Generic (remaining 25) | None |

**Only r=21 mod 64 is the genuine boundary residue** (non-deterministic at higher mod). Chang's I_2 residues are fully deterministic — they don't share the boundary-non-determinism property.

## 3. Step 3: Per-residue orbit observables at N=2^32 (5M orbits)

Walk 5M orbits from uniform odd integers in [1, 2^32], stratify by starting residue mod 64.

| Group | n_residues | ⟨σ⟩ over residues | ⟨V⟩ over residues | n_orbits |
|---|---|---|---|---|
| Chang I_2 ({7, 27, 31, 59, 63}) | 5 | **81.19 ± 3.96** | **2.0208 ± 0.025** | 779,327 |
| Our singular ({21, 53}) | 2 | **65.33 ± 3.37** | **2.1535 ± 0.037** | 312,213 |
| Generic (other 25) | 25 | 74.52 ± 4.38 | 2.0692 ± 0.035 | 3,908,460 |

**Deviations from generic baseline:**

| Comparison | Δ⟨σ⟩ | Δ⟨V⟩ | direction |
|---|---|---|---|
| Chang I_2 vs Generic | **+6.68** | **-0.048** | σ HIGH, V LOW (slow descent) |
| Our singular vs Generic | **-9.19** | **+0.084** | σ LOW, V HIGH (fast descent) |

**OPPOSITE-DIRECTION effects.** The two characterizations are tracking residue-selected orbits with structurally OPPOSITE properties.

P(j=2,4,5) is INDISTINGUISHABLE across all three groups (within sampling noise: range 0.9369-0.9389 for P(j=2) across all 32 residues). The j-class outcome is determined later in orbit by descent path, not by single first-step.

## 4. Mechanistic interpretation

The mechanism:
- **First Syracuse step**: m → (3m+1)/2^v_2(3m+1). The log-step magnitude is log(3) - v·log(2) ≈ 1.099 - 0.693·v.
- For v=1 (Chang's I_2 residues): log-step = +0.405 nats (ASCENT — orbit grows after first step!)
- For v=6 (our r=21 boundary): log-step = -3.06 nats (large DESCENT)

Chang's I_2 picks the **slowest-descent / ascending residues**. Our r=21 picks the **fastest-descent residue** (boundary).

These are structurally **opposite** even though both are "singular" in the sense of being non-generic.

## 5. Verdict per brief outcomes

- **(a) Same structural mechanism:** NO. Opposite-direction effects rule this out.
- **(b) Partial connection (analogous reasons):** PARTIAL. Both indexed via v_2(3r+1) arithmetic, but at opposite extremes.
- **(c) Independent observables:** YES (functionally). Chang's I_2 = slowest-descent, our r=21 = fastest-descent.

## 6. What this means for Chang contact

If the goal is to claim "we've found the same structure": **NO, we haven't.** Chang's I_2 ⊂ v_2=1 residues; our r=21 ⊂ v_2=6 boundary. They're tracking different ends of the same arithmetic spectrum.

If the goal is to share a meta-level observation: **YES, useful.** The arithmetic spectrum v_2(3r+1) per residue has two structural extremes:
- v_2=1 (16 residues, slowest descent) — Chang's spectral methods focus here
- v_2=6 (1 residue, fastest descent, boundary) — our trajectory-measure methods focus here

The remaining residues (v_2 ∈ {2, 3, 4, 5}) are intermediate. A unified framework might characterize the v_2 spectrum as a whole; Chang and we are looking at opposite poles.

## 7. Hausdorff dimension cross-check (Step 5): NOT performed

Chang's reported Hausdorff dim ≈ 0.68 in Z_2 for the divergent-starting-points set requires understanding their specific framework (transfer-operator measure on 2-adic integers). Without seeing the paper carefully, I can't compute the analog quantity for our T-invariant {n ≢ 0 mod 3} subset to test for clean relationship.

Empirically: at N=2^36 (existing data), the {n ≢ 0 mod 3} subset has natural density 2/3 in Z+, so its "dimension" in Z_2 is 1 (full). This doesn't match 0.68. But the relevant comparison would be the dimension of orbits absorbing at {m_j} vs Chang's "divergent set" — these are different objects.

## 8. Files

- `chang_invariant_core_test.py` — full computation
- `chang_invariant_core_test.md` — this document (Result 41)
- `experiments_output/chang_invariant_core_test.csv` — per-residue orbit observables
- `experiments_output/chang_invariant_core_test_log.txt` — full log

Compute: 1.0s walk + analysis.
