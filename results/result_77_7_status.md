# R77.7 — k=7 ε extension and Padé refresh — STATUS: NOT COMPLETED, SUPERSEDED

**Date:** 2026-05-04. Companion to [result_77_6_generating_function.md](result_77_6_generating_function.md). Records the status of R77.7's Markov-chain extension to k=7.

## Status

> **NOT COMPLETED.** Compute killed at ~8.5 hr wall-time runtime by user direction. R77.7's k=7 stationary distribution at q=3 (N=1458 coprime states, M=1458, O(N³) Gauss elimination over Q with growing-denominator Fractions) was actively grinding (96% CPU, 1.9 GB RAM, no signs of stuck or OOM) but past the 3-6 hr original estimate. ε_7 was not produced; the m+n=5 Padé refresh of E(z) := Σ ε_n z^n was not run; the (G-power) vs (G-log) discrimination on E(z)'s singularity at z=2 remains unresolved at the m+n=4 level from R77.6.

## Why superseded

R77.7's analytical purpose was to distinguish power-law from logarithmic asymptotic behavior of ε_n via deeper Padé approximants of the generating function — an **indirect probe** of analytical obstructions in the trajectory measure's Fourier structure.

The **joint 2-adic / 3-adic Bohr empirical positive** ([result_bohr_probe.md](result_bohr_probe.md), 2026-05-04) **directly observes** the analytical obstruction R77.7 was reaching for: chi²-departures from CRT independence on (Z/2^a)* × (Z/3^b)* growing to z = 16.5 at k=20 (a=5, b=4). This is the same multiplicative joint structure the trajectory measure carries, observed at the empirical level on N=10⁷ Syracuse trajectories. R77.7's analytic-structure-of-E(z) probe is now **substantially superseded** by direct empirical observation of the joint Bohr coupling.

## Re-fire conditions

R77.7 (k=7 + Padé refresh) can be re-fired in the future if:
1. The Bohr finding (`result_bohr_probe.md` + verification chain) requires corroboration from the ε_n side.
2. The m+n ≤ 4 Padé from R77.6 is challenged on its (G-branch-cut, type indeterminate) verdict and a deeper diagonal probe becomes load-bearing.
3. The compute budget admits 10–20 hr of single-threaded Fraction work (e.g., during an idle period).

Until then: R77.7 is parked. R77.6's m+n=4 Padé verdict (branch-cut at z=2, power vs log indeterminate at N=5) stands as the current operator-shape characterization for ε_n.

## Files

- [result_77_7_extend_to_k7.py](result_77_7_extend_to_k7.py) — script (preserved; cache up through k=6 in [experiments_output/result_77_7_eps_exact_through_k7.json](experiments_output/result_77_7_eps_exact_through_k7.json), which still exists and is reusable)
- [result_77_7_log.txt](result_77_7_log.txt) — partial stdout (cache restore only; no k=7 print)
- [result_77_7_status.md](result_77_7_status.md) — this writeup

## Cache retained

The exact-rational ε_1..ε_6 cache from R77.6's Stage 1 is preserved in [experiments_output/result_77_7_eps_exact_through_k7.json](experiments_output/result_77_7_eps_exact_through_k7.json). Future re-fires of R77.7 can resume from this cache (skipping the 7+ minutes of k=1..6 rebuilding) and only run the k=7 Markov chain.

## Pending sibling task — q=7 k=4 (also killed)

The q=7 k=4 probe ([c_tilde_q7_k4_probe.py](c_tilde_q7_k4_probe.py)) was also launched (~20-25 hr estimate, ~30 min into the run) and killed in the same compute snapshot. q=7's anomaly in the c̃_q structure test (δ ≈ 0.21 from (q−3)/q) remains unresolved — could be finite-k transient or q-specific arithmetic. The q=17 probe ruled out the "non-prim-root pattern" explanation, so q=7 is genuinely q-specific or finite-k. Re-fire conditions same as R77.7: only if a follow-up requires the q=7 k=4 datapoint specifically.
