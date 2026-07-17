# c̃_q structure test — VERDICT: SUGGESTIVE clean formula, ambiguous at 5 datapoints

**Date:** 2026-05-04. Sibling-probe Task 3 (Wilson-prompted). Tests structural candidates for c̃_q := lim S_k^{(q)} / (q/3)^k across q ∈ {3, 5, 7, 11, 13}.

## Verdict (one paragraph)

> **Outcome C: SUGGESTIVE but ambiguous.** A clean candidate emerges — **c̃_q ≈ (q − 3)/q** — fitting q=11 and q=13 to within **0.2%** as exact rationals would predict, with the implied multiplicative constant C(q) := c̃_q · q / (q − 3) → 1 (1.0021 at q=11, 1.0007 at q=13). q=3 is in a separate regime entirely (q/3 = 1 makes the (q−3)/q form identically 0; c̃_3 = 7/15 is the unrenormalized Tao limit). q=5 (still settling at k=4, Aitken extrapolation suggests it converges to ~0.482) and q=7 (the only non-(2-prim-root) case in our set, c̃_7 ≈ 0.78 with deviation 0.21 from (q−3)/q) DEVIATE substantially. Cannot determine from 5 points whether q=5's deviation is finite-k transient or structural, or whether q=7's deviation reflects the broader non-prim-root pattern.

> **Recommended next probe: compute c̃_17 at k=2 (cheap, minutes).** q=17 has ord(2 mod 17) = 8 (NOT primitive root, like q=7). If c̃_17 ≈ (17−3)/17 = 14/17 ≈ 0.824 cleanly, q=7's deviation is q-specific finite-k. If c̃_17 has δ ≈ 0.2 like q=7, there's a non-prim-root structural correction.

## Candidate-by-candidate fit quality

### Empirical c̃_q values (last available k)

| q | k_max | c̃_q | ord(2 mod q) | 2 prim root? |
|---|---|---|---|---|
| 3 | 4 | 0.464214 | 2 | YES (2 = q−1) |
| 5 | 4 | 0.487717 | 4 | YES (4 = q−1) |
| 7 | 3 | 0.782583 | 3 | NO (3 ≠ 6 = q−1) |
| 11 | 3 | 0.728798 | 10 | YES (10 = q−1) |
| 13 | 2 | 0.769793 | 12 | YES (12 = q−1) |

### Candidate 1: ord(2 mod q) / (q−1) correlation

Pearson r(ratio, c̃_q) = **−0.48** (weak, wrong sign for the obvious story). The q=7 row has ord/(q−1) = 0.5; all others 1.0 — a binary distinction. But q=7's c̃ (0.78) sits between q=11 (0.73) and q=13 (0.77), not clearly outlying.

**Verdict: ord(2 mod q) alone does NOT predict c̃_q ordering.** The non-prim-root status doesn't simply flag q=7 as the outlier.

### Candidate 2: q=7 outlier, simple-form fits on {3, 5, 11, 13}

| form | fit | RSS |
|---|---|---|
| linear in q: 0.0331·q + 0.348 | R² = 0.98 | 1.30e−3 |
| 1/q: 0.827 − 1.222/q | — | 1.25e−2 |
| log q: 0.182 + 0.224·log q | — | 4.57e−3 |

Linear fit RSS is small, but with 4 points and 2 params, R² = 0.98 is essentially **trivial fit quality** (only 2 degrees of freedom; can't distinguish from random). Residuals 0.017, −0.026, +0.017, −0.008 don't show a clear systematic pattern.

**Verdict: simple linear fit works but isn't constraining at N=4.**

### Candidate 3: Step function vs smooth

Cluster analysis:
- Low cluster {q=3, 5}: c̃ ∈ {0.464, 0.488}, mean 0.476, std 0.012
- High cluster {q=7, 11, 13}: c̃ ∈ {0.78, 0.73, 0.77}, mean 0.760, std 0.023
- Step magnitude: 0.284
- Max within-cluster std / step = **0.081** (well below 0.2 threshold)

**A step function fits.** But the empirical step location (between q=5 and q=7) doesn't have an obvious arithmetic motivation — neither (q−3)/q nor 2-primitivity suggests a step there. Could be coincidental clustering at N=5.

### Candidate 4: Closed-form (q − 3)/q — THE LEADING SIGNAL

c̃_q · q values:

| q | c̃_q · q | (q − 3) | deviation |
|---|---|---|---|
| 3  | 1.393 | 0  | +1.393 (special: q=3 is the boundary) |
| 5  | 2.439 | 2  | +0.439 |
| 7  | 5.478 | 4  | +1.478 (anomalous) |
| 11 | **8.017** | **8** | **+0.017 (0.2%)** |
| 13 | **10.007** | **10** | **+0.007 (0.07%)** |

Implied C(q) := c̃_q · q / (q − 3):

| q | C(q) | dev from 1 |
|---|---|---|
| 5 | 1.219 | +0.22 |
| 7 | 1.370 | +0.37 |
| **11** | **1.0021** | **+0.002** |
| **13** | **1.0007** | **+0.001** |

**The (q−3)/q candidate fits q=11 and q=13 to within 0.2%.** Both cases where 2 is primitive root mod q AND q is "large enough" (q ≥ 11). At q=11 and q=13 we see the universal-shape structure cleanly.

(q−3)/q means c̃_q · q is an INTEGER (= q − 3). Empirically q=11 gives 8.017 and q=13 gives 10.007 — both within 0.2% of integers 8 and 10. **This is the cleanest signal in the data.**

For q=5, 7: deviations 0.22 and 0.37 in C(q). Cannot tell from 5 points whether these are finite-k transients (will converge to 1 with deeper k) or structural deviations.

For q=3: (q−3)/q = 0/3 = 0 identically, but c̃_3 = 7/15. **q=3 is a separate regime**: q/3 = 1 means there's no renormalization, and c̃_3 IS the unrenormalized Tao limit S_∞^{(3)} = 7/15.

### Candidate 5: Unifying formula across q=3 and q ≥ 5

**Two-regime structure:**

- **q = 3 (boundary, no renormalization needed):** c̃_3 = lim S_k^{(3)} = 7/15 ≈ 0.4667. Empirical 0.4642 (k=4, still settling toward 7/15 = 0.4667).
- **q ≥ 5 (renormalized regime):** c̃_q = (q − 3)/q + δ(q) where δ → 0 for large q with 2 prim root.

Empirical δ:

| q | δ = c̃_q − (q−3)/q | reading |
|---|---|---|
| 5  | +0.088 | small — may converge to 0 (Aitken extrapolation gives ~0.482, δ ~ 0.08) |
| 7  | +0.211 | LARGE — only non-(2-prim-root) case |
| 11 | +0.0015 | essentially 0 |
| 13 | +0.0006 | essentially 0 |

The δ(7) ≈ 0.21 is the only deviation that **resists explanation as finite-k transient.** q=7's c̃ at k=2 was 0.7817 and at k=3 was 0.7826 — stable, not still decreasing. So q=7 likely has a real non-prim-root deviation.

**Convergence check on q=5:**
- c̃_5 sequence at k=1..4: 0.5333, 0.4922, 0.4896, 0.4877 (slowly decreasing).
- Aitken Δ² extrapolation: 0.4895 (k=0..2), 0.4822 (k=1..3).
- Aitken predicts c̃_5 → ~0.48, NOT 0.40 (= (q−3)/q).
- So δ(5) ≈ 0.08 is **likely real**, not just finite-k. q=5 also doesn't fit (q−3)/q cleanly.

## Synthesis

What the data IS telling us cleanly:
1. **q=11, q=13: c̃_q = (q−3)/q to within 0.2%.** Strong evidence for the leading-order formula at "large enough" q with 2 primitive root.
2. **q=3: separate regime** (boundary case where renormalization is unnecessary). c̃_3 = 7/15 from forward Tao limit.
3. **q=5, 7: deviate from (q−3)/q.** q=5 modestly (~0.08), q=7 strongly (~0.21).

What the data CANNOT distinguish at N=5:
- Is q=5's deviation finite-k (will → 0 with more k) or structural?
- Is q=7's deviation a "non-2-prim-root" pattern (would predict similar deviation at q=17, 23, etc.) or q-specific?
- Is the convergence to (q−3)/q monotonic in q, or are there oscillations?

## Recommended next probe — q=17

**Cheapest decisive test:** compute c̃_17 at k=2.
- q=17 has ord(2 mod 17) = 8, NOT primitive root (q − 1 = 16).
- Q=17, k=2: N = 17² = 289, M = ord_{17²}(2) = 8·17 = 136. State count = 16·17 = 272 coprime residues. Compute time: ~1 minute (much cheaper than q=13 k=3 which had 2028 states).

**Decision rule:**
- If c̃_17 ≈ (17−3)/17 = 14/17 ≈ 0.8235: q=7's deviation is finite-k transient (q=17 at k=2 already at limit because q is large enough). Promotes the (q−3)/q hypothesis.
- If c̃_17 has δ ≈ 0.2 (i.e., c̃_17 ≈ 1.02): there's a structural non-prim-root correction. Need to characterize.
- If c̃_17 is somewhere else (e.g., 0.7): suggests c̃_q has more structure than (q−3)/q + non-prim-root correction; needs a richer hypothesis.

**Alternative (more expensive but firmer):** extend q=5 to k=5 (~30 min).
- Would firm c̃_5 to 5+ sig figs and resolve whether δ(5) is structural or transient.
- But less leverage than q=17 (one more datapoint at known q vs new datapoint at new q with known structural feature).

## Files

- [c_tilde_structure_test.py](c_tilde_structure_test.py) — analysis script
- [c_tilde_structure_verdict.md](c_tilde_structure_verdict.md) — this writeup

## Honest caveats

1. **5 datapoints with 4 free parameters** in any rational form is overfitting territory. The (q−3)/q candidate has zero free parameters (it's a specific functional form), so the q=11/q=13 fit at 0.2% is genuinely surprising — but two points isn't sufficient to call it confirmed.
2. **q=5 still settling.** k=4 is the deepest we have; the c̃_5 sequence hasn't visibly stabilized. Aitken extrapolation suggests ~0.48 at limit, but this isn't decisive.
3. **q=3 special status.** Cannot fit into the (q−3)/q family by definition (since q − 3 = 0 there). Need a unifying limit framework (e.g. analytic continuation as a function of q) to bridge.
4. **The "step function" reading** (Candidate 3 fitting at the within-cluster-std level) is real but suggests the wrong structural mechanism — there's no arithmetic reason for a step at q ≈ 6 specifically.
5. **The q=7 anomaly** could mean (a) finite-k transient (q=7 only has 3 datapoints; maybe c̃_7 hasn't converged), (b) structural non-prim-root effect, or (c) q-specific arithmetic. q=17 would help distinguish (a)+(c) from (b).

## Verdict, restated

> **C: SUGGESTIVE.** The (q − 3)/q form fits q=11 and q=13 to within 0.2%, suggesting it's the leading-order formula at large q with 2 primitive root. The deviation pattern at q=5, 7 is consistent with (small-q transient + non-prim-root correction) but cannot be confirmed at N=5. Recommend q=17 at k=2 (~1 min compute) as the cheapest discriminating probe.
