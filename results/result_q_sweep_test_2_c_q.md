# Q-sweep test 2 — c_q across odd primes (TENTATIVE, EARLY-TERMINATED)

**Date:** 2026-05-04. Tests whether 7/45 is the q=3 instance of a closed-form family c_q across qx+1 systems for q ∈ {3, 5, 7, 11, 13}. Run terminated early at user direction after q=11 k=3 finished and q=13 k=2 confirmed the structural pattern; q=13 k=3 was killed (estimated ~21 hours).

## Verdict (tentative): hybrid (NO-PATTERN literal) / (UNIVERSAL-SHAPE renormalized)

> **The literal hypothesis is FALSIFIED**: S_∞^{(q)} := lim_{k→∞} S_k^{(q)} does NOT exist as a finite limit for any q ≥ 5. The sequence S_k^{(q)} grows asymptotically like (q/3)^k. Therefore c_q = S_∞^{(q)}/q is well-defined ONLY at q=3 (where q/3 = 1 makes the limit finite). 7/45 is q=3-specific.
>
> **A renormalized analog DOES exist**: S_{k+1}^{(q)} / S_k^{(q)} → q/3 universally across all tested q, and c̃_q := lim_{k→∞} S_k^{(q)} / (q/3)^k converges to a finite constant for every q. The c̃_q values across q don't show an obvious closed form at this 5-q resolution, but the existence and stability of the limit is itself a structural finding.

The q-sweep is best read as **two findings**: (1) the original hypothesis is wrong, and (2) the right hypothesis is the renormalized version. Worth follow-up.

## Stage 0 preflight (PASSED before main run; archived in result_q_sweep_test_2_preflight.md)

- All 20 ord_{q^k}(2) values match expected (q=7 critically gives 3, 21, 147, 1029 — NOT 6, 42, 294, 2058).
- q=3 k=3 sanity: generalized chain reproduces S_3 = 31370/67963 over Q.
- Plancherel at q=5 k=2: π_2^{(5)} fiber-projection equals π_1^{(5)} = (1/15, 2/15, 8/15, 4/15) exactly.

## Stage 1+2 results — exact rational X_k^{(q)} = q^k · ‖π_k^{(q)}‖²

Computed: q=3 (k=1..4), q=5 (k=1..4), q=7 (k=1..3), q=11 (k=1..3), q=13 (k=1..2).

| q | k | states | M(q,k) | X_k float | S_k = X_k − X_{k-1} float |
|---|---|---|---|---|---|
| 3  | 1 | 2    | 2    | 1.6667     | 0.6667 |
| 3  | 2 | 6    | 6    | 2.1429     | 0.4762 |
| 3  | 3 | 18   | 18   | 2.6044     | 0.4616 |
| 3  | 4 | 54   | 54   | 3.0686     | 0.4642 |
| 5  | 1 | 4    | 4    | 1.8889     | 0.8889 |
| 5  | 2 | 20   | 20   | 3.2561     | 1.3673 |
| 5  | 3 | 100  | 100  | 5.5230     | 2.2668 |
| 5  | 4 | 500  | 500  | 9.2862     | 3.7633 |
| 7  | 1 | 6    | 3    | 3.0000     | 2.0000 |
| 7  | 2 | 42   | 21   | 7.2560     | 4.2560 |
| 7  | 3 | 294  | 147  | 17.1977    | 9.9417 |
| 11 | 1 | 10   | 10   | 3.6738     | 2.6738 |
| 11 | 2 | 110  | 110  | 13.4723    | 9.7984 |
| 11 | 3 | 1210 | 1210 | 49.3993    | 35.9270 |
| 13 | 1 | 12   | 12   | 4.3354     | 3.3354 |
| 13 | 2 | 156  | 156  | 18.7905    | 14.4550 |

Exact rationals stored in [experiments_output/result_q_sweep_test_2_cache.json](experiments_output/result_q_sweep_test_2_cache.json) and [experiments_output/result_q_sweep_test_2_table.csv](experiments_output/result_q_sweep_test_2_table.csv).

## Stage 3 — convergence analysis: S_k^{(q)} grows like (q/3)^k

Examining successive S-ratios:

| q | S_2/S_1 | S_3/S_2 | S_4/S_3 | q/3 |
|---|---|---|---|---|
| 3  | 0.7143 | 0.9693 | 1.0057 | **1.0000** |
| 5  | 1.5382 | 1.6579 | 1.6601 | **1.6667** |
| 7  | 2.1280 | 2.3359 | —      | **2.3333** |
| 11 | 3.6646 | 3.6666 | —      | **3.6667** |
| 13 | 4.3337 | —      | —      | **4.3333** |

**Reading**: the asymptotic ratio S_{k+1}^{(q)} / S_k^{(q)} → **q/3** universally. q=11, q=13 hit q/3 already at k=2 (4 sig figs); q=5, q=7 still converging at k=3. q=3 (q/3 = 1) is the borderline case where ratio → 1 from below, giving the finite Tao limit S_∞^{(3)} = 7/15.

For q ≥ 5, ratio > 1 means S_k diverges geometrically. **S_∞^{(q)} = ∞ for q ≥ 5.**

## Stage 4 — renormalized constant c̃_q := lim S_k^{(q)} / (q/3)^k

| q | k=1 | k=2 | k=3 | k=4 | apparent limit c̃_q |
|---|---|---|---|---|---|
| 3  | 0.6667 | 0.4762 | 0.4616 | 0.4642 | **7/15 ≈ 0.4667** |
| 5  | 0.5333 | 0.4922 | 0.4896 | 0.4877 | ≈ **0.487** (still settling) |
| 7  | 0.8571 | 0.7817 | 0.7826 | —      | ≈ **0.78** |
| 11 | 0.7292 | 0.7288 | 0.7288 | —      | ≈ **0.7288** (stable) |
| 13 | 0.7697 | 0.7698 | —      | —      | ≈ **0.7698** (only 2 points) |

The **existence and finiteness of c̃_q is universal**. The values themselves don't show an obvious closed form at this resolution:

- q=3 ↔ 7/15 = 0.4667
- q=5 ≈ 0.487
- q=7 ≈ 0.78
- q=11 ≈ 0.7288
- q=13 ≈ 0.7698

The q=3 and q=5 values are close (~0.47, 0.49); q=7, 11, 13 cluster higher (~0.73-0.78). q=7 has ord_7(2) = 3 (not primitive), unlike the others. But this doesn't separate q=3,5 from q=11,13 cleanly. **Closed-form testing on five points is suggestive but unresolved at best.**

Tested simple rational families (`c̃_q = a/q`, `(q-1)/q²`, `(q-1)/(q(q+1))`, etc.) by fitting A in `c̃_q = A · g(q)`; relative-std of A across q exceeded 30% for every g(q) tried. **No single-parameter rational form fits.**

## Why does this happen? Heuristic

Plancherel: X_k^{(q)} = q^k · ‖π_k^{(q)}‖² counts total Fourier mass. For a uniform π on the (q-1)·q^{k-1} q-coprime states, ‖π‖² = 1/((q-1)q^{k-1}), giving X_k^{(uniform)} = q/(q-1) — bounded.

Tao 2019 establishes that for q=3 the actual ‖π_k^{(3)}‖² lies in a regime where X_k^{(3)} grows linearly in k (i.e., increment X_k − X_{k-1} → const = 7/15). Empirically here for q ≥ 5 we observe **‖π_k^{(q)}‖² ~ const · (q/3)^k / q^k = const / 3^k**, so X_k^{(q)} ~ const · (q/3)^k. The factor 1/3^k in ‖π‖² is suggestive of the q=3 "shadow" — possibly because the dynamics retain some 3-adic structure across all q (the multiplier in the iteration is q, but the stripping prime is 2, and 3 enters somehow as an artifact of the qx+1 framework).

**Speculative interpretation:** the (q/3)^k growth rate might reflect the dynamical contraction/expansion factor of qx+1. The mean log-growth per Syracuse step at q is `log(q/2)/log 2 = log_2(q/2)`. For q=3: log_2(3/2) ≈ 0.585. For q=5: log_2(5/2) ≈ 1.322. The ratio (q/3) doesn't directly equal these, so the connection (if any) is via a different mechanism.

This is a **genuine open question**: why is the universal asymptotic ratio exactly q/3?

## Limitations / honest caveats

1. **q=13 only k=2** — c̃_13 estimate is from a single data point. Could be unstable.
2. **q=5, q=7 still settling at k=3 or k=4** — the c̃ values listed are upper bounds; true limit may be 1-3% lower.
3. **No closed-form found** but only 5 data points; a richer family (e.g. q ∈ {3, 5, 7, 11, 13, 17, 19, 23}) might reveal one.
4. **q=7 is the only non-primitive case** in the set — possibly an outlier; need q ∈ {15-th cyclotomic, 23, etc.} to test.
5. **Plancherel identity for general q** verified at q=5 k=2 only; should hold by the same argument as Theorem 75.1, but k≥3 verifications would strengthen.
6. **The (q/3) ratio observation** may have a clean derivation from Tao's framework that we haven't pursued. Worth ~half a day of analytic work.

## What changes for the c=7/45 closure question

**Nothing direct.** The q=3 framework is unchanged. The finding that 7/45 is q=3-specific doesn't help or hurt the rate-1/2 closure — it just clarifies that publishing 7/45 doesn't generalize to a c_q family.

But the **renormalized c̃_q hypothesis** (and the q/3 universal ratio) is publishable on its own as a Plancherel-mass scaling theorem for qx+1 Markov chains, independent of any Collatz closure attack.

## Outcome classification (per brief's decision tree)

- (FAMILY-CLOSED): NO at literal level (S_∞ doesn't exist for q ≥ 5).
- (UNIVERSAL-SHAPE): **YES at renormalized level** — S_{k+1}/S_k → q/3 across all q, and c̃_q exists for all q tested.
- (NO-PATTERN): NO — there IS structure (q/3 universal ratio).

Best classification: **(NO-PATTERN at literal hypothesis) + (UNIVERSAL-SHAPE at renormalized level).** A genuine partial finding: the brief's hypothesis is wrong about which constant has the structure, but a closely-related constant DOES have universal structure.

## Files

- [result_q_sweep_test_2_preflight.py](result_q_sweep_test_2_preflight.py) / output  — Stage 0 verifications
- [result_q_sweep_test_2.py](result_q_sweep_test_2.py) — main run (terminated at q=13 k=2)
- [result_q_sweep_test_2_log.txt](result_q_sweep_test_2_log.txt) — full stdout
- [experiments_output/result_q_sweep_test_2_cache.json](experiments_output/result_q_sweep_test_2_cache.json) — exact-rational X_k^{(q)} cache
- [experiments_output/result_q_sweep_test_2_table.csv](experiments_output/result_q_sweep_test_2_table.csv) — float table
- [result_q_sweep_test_2_c_q.md](result_q_sweep_test_2_c_q.md) — this writeup

## Recommended follow-up (not pursued)

1. **Push q=5 to k=5** (cache to k=5 would take ~30 min) to firm c̃_5 to 4 sig figs.
2. **Run q ∈ {17, 19, 23}** at k=2 (cheap; matter of minutes) to expand the (q, c̃_q) table.
3. **Derive the (q/3) ratio analytically** from Tao's Plancherel framework. Likely a clean argument exists.
4. **Check if c̃_q has a closed form involving ord_q(2)** — q=7's deviation might be the clue.
5. **The "c̃_q = lim S_k^{(q)} / (q/3)^k exists" conjecture** is a clean publishable theorem candidate, independent of Collatz closure status.
