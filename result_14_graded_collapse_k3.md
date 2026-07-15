# Result 14 (qx+1 paper) — R11's cell collapse GRADES to k=3, and the grading is TIGHT. The route to Result 1 is live.

**Date:** 2026-07-15. **Verdicts: H_M(gate) ✓ / H_G1 ✓ EXACT (the k=2 fact survives) / ★ H_G2 ✓ EXACT (the collapse GRADES) / H_SHARP ✓ TIGHT both ways / H_CELLS ✓ = d³q³ — and the pre-stated caveat holds.**

**Headline: `v_1` mod `d` · `v_2` mod `dq` · `v_3` mod `dq²` — a triangular grading, confirmed by exhaustive exact-integer equality (zero mismatches in 3.1M checks at q=7). This is the first forward step on Result 1's route since R6 killed the old one.**

Probe: `probe_14_graded_collapse_k3.py`. Log: `result_14_graded_collapse_log.txt`. Runtime: **7.4 s**.

## The derivation (done before running)

Mod `q^k` the `r_0` term dies (R8), leaving `value(v_1..v_k) = Σ_{m=1}^{k} q^{m−1}·2^{−S_m}`, `S_m = v_{k−m+1}+…+v_k`. At k=3, checked against the chain by hand:

> `value(v_1,v_2,v_3) = 2^{−v_3} + q·2^{−(v_2+v_3)} + q²·2^{−(v_1+v_2+v_3)} mod q³`

The m-th term carries `q^{m−1}`, so it needs `2^{−S_m}` only **mod `q^{k−m+1}`**; and `2^{−S} mod q^j` depends on `S` only mod `ord_{q^j}(2) = d·q^{j−1}`. Hence `S_m` matters only mod `d·q^{k−m}`, giving a **triangular grading**:

| coordinate | appears in | matters mod |
|---|---|---|
| `v_1` | `S_k` only (needed mod q) | **`d`** |
| `v_2` | `S_{k−1}, S_k` | **`d·q`** |
| `v_k` | all | `d·q^{k−1} = M` (full range — no tower) |

At k=2 this collapses to R11's exact fact. **k=3 is the first real test.**

## Results — exact integer equality, no tolerances

**H_M (gate):** `ord_{q²}(2) = dq` and `ord_{q³}(2) = dq²` at q=5, 7, 11, 13 — all non-Wieferich ✓

| q | d | M | mode | H_G1 `value(v_1+d,·,·)` | H_G2 `value(·,v_2+dq,·)` |
|---|---|---|---|---|---|
| 5 | 4 | 100 | **exhaustive** (10⁶ triples) | 960,000 checks, **0 mismatches** | 800,000 checks, **0 mismatches** |
| 7 | 3 | 147 | **exhaustive** (3,176,523 triples) | **3,111,696 checks, 0 mismatches** | **2,722,734 checks, 0 mismatches** |
| 11 | 10 | 1210 | sweeps (169 slices × full range) | 202,800 checks, 0 mismatches | 185,900 checks, 0 mismatches |
| 13 | 12 | 2028 | sweeps (144 slices × full range) | 290,304 checks, 0 mismatches | 269,568 checks, 0 mismatches |

- **H_G1 CONFIRMED** — R11's k=2 fact (`v_1` mod `d`) survives at k=3.
- **★ H_G2 CONFIRMED — THE COLLAPSE GRADES.** `v_2` matters only mod `dq`. This is the new claim and it is exact.

**H_SHARP — TIGHT in both directions, 216/216 at every prime:** `v_2+d` **changes** the value (so `d` genuinely does not suffice for `v_2`), and `v_3+dq` **changes** the value (so `dq` does not suffice for `v_3`). The grading is neither a loose upper bound nor secretly stronger than derived.

**H_CELLS:** `cells = d·(dq)·M = d³q³` exactly at all four primes.

## ⚠️ The pre-stated caveat holds — this is STRUCTURE, not the rate

| q | d | cells | states ≈ q²(q−1) | cells/states |
|---|---|---|---|---|
| 5 | 4 | 8,000 | 100 | 80.0 |
| 7 | 3 | 9,261 | 294 | 31.5 |
| 11 | 10 | 1,331,000 | 1,210 | 1,100.0 |
| 13 | 12 | 3,796,416 | 2,028 | 1,872.0 |

`cells/states = d³·q/(q−1) ≈ d³`. **The collapse alone gives NO injectivity** — flagged in the pre-registration before running, and it holds. The weights must do the work (effective `3^k`). **Result 1 remains unproved.** In general: `cells = d^k·q^{k(k−1)/2}` vs `~q^k` states, so `cells/states = d^k·q^{k(k−3)/2}`.

## ★ What the structure hands us for free — the next prediction

The tower ratios are `x_j = 2^{−d·q^{j−1}}`: **`2^{−d}`, `2^{−dq}`, `2^{−dq²}`, …** — **doubly** exponentially separated. At q=41: `2^{−20}`, `2^{−820}`, `2^{−33620}`. Only `v_1..v_{k−1}` carry towers; `v_k` ranges over exactly its modulus and has none.

Summing the towers (same computation as R11, per coordinate) predicts:

> **`ratio_within(k) = [∏_{j=1}^{k−1} (1+x_j)/(3(1−x_j))] / P2^{k−1} − 1` → `∏_{j=1}^{k−1}(1+x_j)/(1−x_j) − 1` ≈ `2·2^{−d}` for EVERY k ≥ 2**

i.e. **the within-cell overlap is essentially k-independent**, because every tower past the first contributes `~2^{−dq}` — not small, *nonexistent* at any precision that exists. At k=2 this reduces to R11's exact identity.

**Consequence, if true:** the within-cell part is *bounded in k* ⇒ it cannot be what breaks domination at q=3. But domination **does** break at q=3 (R8: `ratio_k ≈ 0.4655·k`, linear). ⇒ **the k-growth at q=3 must live entirely in the CROSS-cell (family-b) part** — the same character-sum object R13 isolated. That is a sharp, pre-committable prediction and it is the next probe.

## Not at stake
R10's law, R11, R13, R5's rate, R6, R7, R12, THEOREM_C_745, Th 78.1–78.3, R81b, ε_k. A refutation here would kill only my proposed route to Result 1, not any banked result.

_Reporting discipline: derived on paper first, then tested as exact integer equality — no tolerance to mis-specify (four of my decision rules have been mis-specified this arc). Exhaustive over all M³ triples where affordable rather than sampled. H_SHARP was included specifically so a collapse that is stronger than derived would be caught and reported as good news rather than buried. The "no injectivity" limitation was stated in the pre-registration, before the numbers, not offered afterward as a caveat. Author's structural priors this arc: 13-for-20._
