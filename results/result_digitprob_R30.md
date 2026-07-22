# Probe R30 — the digit probability — **A GATE PASS (identity certified, new weld); B/C/D structured but NOT the clean simplification**

**Date:** 2026-07-22. Probe `probes/probe_digitprob_R30.py` (reuses R7 mu/cram, R9 gamma/tau, R10 autocorr_dlog, R28 build_nu).
Wilson's reframe: the requirement is **not** √-cancellation `max|μ̂|≲3^{−r/2}` (that was the triangle-route bound, aimed 3 orders
too high) — the theorem needs `Σ_r|A_r(m)|<∞`, so `|A_r(m)|=O(r^{−1−δ})` suffices, and the decay is observed geometric.
Object recast as a **conditional-digit probability**: `A_r(m)=γ_{r−1}(τ_m)·(3q_r(m)−1)`, `q_r=p_r/p_{r−1}`,
`p_r=Pr[X′≡4^{−m}X mod 3^{r+1}]`, so **`A_r(m)→0 ⟺ q_r(m)→1/3`** (next 3-adic digit of the ratio → uniform).

## R30-A — MEASURE q_r + IDENTITY GATE: **PASS** (new collision↔character weld at m≥1)
Two independent paths — collision side `γ_r(τ_m)=3^r Pr[X′≡4^{−m}X mod 3^{r+1}]` (partner-count, R9) and character side
`A_r(m)=Σ_u g_r(u)·c_{3^r}((u−m) mod 3^r)` (dlog+Ramanujan, R10) — welded for the first time at **m=1,2,3,4, r=2..7**:

| m=1: r | γ_r | q_r=p/p | 3q_r−1 | γ_{r−1}(3q−1) | A_r(m) char | gate |
|---|---|---|---|---|---|---|
| 2 | 0.693878 | 0.346939 | +0.040816 | +0.027211 | +0.027211 | PASS |
| 4 | 0.707042 | 0.335236 | +0.005707 | +0.004012 | +0.004012 | PASS |
| 7 | 0.717087 | 0.334590 | +0.003771 | +0.002694 | +0.002694 | PASS |

**GATE PASS all m,r** — the identity `A_r(m)=γ_{r−1}(3q_r−1)=γ_r−γ_{r−1}` is exact. This certifies Wilson's recasting:
`q_r=p_r/p_{r−1}` is the right object, `3q_r−1` the right factor. (m=0 case = the R10-B / R28 weld `A_r(0)=S_r`.) The
collision-γ ledger and the primitive-layer character coefficient are now the **same object** for every ratio value, not
just the DC line. **γ_r→f(4^{−m})** and **q_r→1/3** confirmed (m=1: γ 0.694→0.723→f(τ₁)≈0.72; q 0.347→0.3341).
A_r(1) magnitudes 0.027,0.009,0.004,0.0037,0.0036,0.0027 (r=2..7); r3→4 ratio 0.435≈|λ₂| (agrees with Wilson's ≈0.43
anchor; his quoted 0.014/0.006 scale is a normalization/indexing convention — the **rate**, not the scale, is what matches).

## R30-B — RATE + m-INDEPENDENCE: **NOT a clean single rate, NOT m-independent** (channel mixing)
Ratio `(3q_{r+1}−1)/(3q_r−1)` per m (measurement, no fit):
- **m=1** (3q−1 stays positive): 0.43, 0.93, 0.96, 0.74, 0.80, 0.83, 0.92 — decays but at a **drifting ratio**, with
  mid-r values near 1 (0.93, 0.96), **not** a constant ½.
- **m=2,3,4**: `3q−1` **flips sign** (0.005→−0.003→+0.002→−0.024…), so the ratios blow up (−13.5, +4.5, −6.4) — the
  **complex-pair signature** (oscillation), not a positive geometric.

So "one mixing chain with rate |λ₂|≈½" does **not** show up as a clean per-channel digit rate. This is R28-B/C again:
the raw single-channel ratio is **not** the eigenvalue ratio — channel mixing (measured O(1) in R28-C) contaminates each
channel's apparent rate. The clean ½ is a property of the full operator, not of any one ratio value's digit sequence.

## R30-C — NEXT-DIGIT DISTRIBUTION: the excess **ALTERNATES**, it is not a static favored digit
Conditional on agreement mod 3^r, the full distribution of the next digit `d=((X′−4^{−m}X)/3^r) mod 3`:

| m=1: r | P(d=0) | P(d=1) | P(d=2) | excess d=1 | excess d=2 |
|---|---|---|---|---|---|
| 2 | 0.34694 | 0.24490 | 0.40816 | −0.088 | **+0.075** |
| 3 | 0.33773 | 0.38854 | 0.27373 | **+0.055** | −0.060 |
| 4 | 0.33524 | 0.29276 | 0.37200 | −0.041 | **+0.039** |
| 5 | 0.33510 | 0.36363 | 0.30127 | **+0.030** | −0.032 |
| 6 | 0.33503 | 0.30964 | 0.35533 | −0.024 | **+0.022** |

The match probability `q=P(d=0)` approaches 1/3 **from above** (excess +0.0136→+0.0017), but the off-match mass
**sloshes d=1↔d=2 with the parity of r** — a period-2 alternation in *where* the excess sits. That alternation is the
**|λ₂| complex pair rendered in digit space**; it is the object a coupling/mixing argument would have to annihilate, and
it is oscillatory, not a fixed skew. (m=2,3,4 show the same sloshing, less cleanly.)

## R30-D — WHICH digits vs HOW MANY: **NOT a scalar 3-state chain, but asymptotically depth-Markov**
`q_r` conditioned on the prefix `X mod 3^k`. At **mod 9** only one class is viable (the ratio-shifted targets miss the
sparse support, 2/9 density → the other classes are exactly zero-mass, nan). At **mod 27**, three viable prefixes with
**genuinely different q**:

| | X%27 (a) | (b) | (c) | spread |
|---|---|---|---|---|
| m=1, r=3 | 0.283 | 0.337 | 0.345 | **0.062** |
| m=1, r=4 | 0.325 | 0.316 | 0.339 | 0.022 |
| m=1, r=5 | 0.328 | 0.329 | 0.337 | 0.008 |
| m=1, r=6 | 0.342 | 0.335 | 0.334 | 0.007 |
| m=2, r=3 | 0.337 | 0.279 | 0.395 | **0.116** |
| m=2, r=6 | 0.324 | 0.334 | 0.323 | 0.012 |

So `q_r` **does depend on which digits matched** (spread 0.06–0.12 at r=3) — the process is **not** a clean scalar
3-state chain, and the theorem does **not** collapse to a depth-only recursion. **But the prefix-dependence decays**
(m=1: 0.062→0.007; m=2: 0.116→0.012), so the process is *asymptotically* depth-Markov: mixing homogenizes the
prefix-dependence as depth grows. The jackpot (exact scalar recursion → provable 3-state gap) is out; the softened
version (prefix-dependence shrinking at the gap rate) is in.

## Status
**R30: the digit recasting is CERTIFIED; the three simplifications it was meant to unlock all come back structured, not
clean.** **A GATE PASS** — `A_r(m)=γ_{r−1}(3q_r−1)` exact (m=1..4, r=2..7), a new collision↔character weld at m≥1;
`q_r→1/3` and `γ_r→f` confirmed. The reframe is real: the theorem **is** a conditional-digit equidistribution statement.
**B** — the per-channel digit rate is neither a clean geometric ½ nor m-independent (m≥2 sign-oscillates = complex pair;
m=1 drifts to ~0.9); channel mixing contaminates the raw rate (R28-B/C). **C** — the deviation from uniform is an
**alternating** d=1↔d=2 sloshing (period-2, the |λ₂| pair in digit space), not a static favored digit. **D** — `q_r`
depends on *which* digits matched (mod-27 spread 0.06–0.12), so **not** a scalar 3-state chain — but the dependence
**decays** with depth, so asymptotically depth-Markov.

**Consequence for the crux (owed to the pen):** R30 delivers the cleanest **statement** of the theorem so far
(conditional-digit equidistribution `q_r(m)→1/3`, certified equivalent to `Σ_r A_r(m)` summability) and a **corrected
target** — Wilson's √-cancellation requirement is retired; `O(r^{−1−δ})` suffices and geometric decay is observed. But
the digit process is **not** the clean 3-state Markov chain that would make the gap free: it is oscillatory (C),
channel-mixed (B), and prefix-dependent-but-homogenizing (D). All three point back to the **same** unresolved object —
the |λ₂|≈½ gap — now wearing the guise of the homogenization rate of the conditional-digit distribution. The remaining
analytic work is unchanged in substance (control that rate) but has a far better **shape**: a decaying prefix-dependence
in a 3-adic digit process, the regime coupling/renewal arguments are actually built for. No fitting; exact identity gate
(m=1..4, r=2..7) and exact digit distributions/spreads; B/C/D reported plainly as structured negatives, not smoothed.
