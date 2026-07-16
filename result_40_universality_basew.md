# Result 40 (qx+1 paper) — UNIVERSALITY REFUTED: the boundary is NOT ord_q(w)=2. q=3 does not replicate at other order-2 halving bases; (q=5,w=4) has a (weak) gap, not the q=3 divergence.

**Date:** 2026-07-16. **Verdict: ✗ H_UNIV REFUTED (clean falsifier: (5,4), predicted critical, GAPS). ★ q=3 confirmed UNIQUELY critical. Bonus observation: w=q−1 (q≥5) cases show a clean weak gap r≈1−1/w² (flagged, not claimed).**

**Headline: I predicted a universality class keyed to `ord_q(w)=2` (halving phase collapses to {1,−1}), with a NEW critical point at (q=5, w=4). FALSIFIED. (5,4) shows a geometric gap `r≈0.9375` (X_k CONVERGES), not the q=3-type linear divergence. `ord_q(w)=2` is necessary-looking but NOT sufficient for criticality — q=3 is special because it is the unique q where the halving base is BOTH order-2 AND primitive (`2≡−1 mod q` with `q−1=2`), which is R6's `(q−1)/d=1` criterion, not `d=2`. Honest negative — the seductive universality framing joins the arc's other cross-domain mirages.**

Probe: `probe_40_universality_basew.py` + higher-k confirmation. Runtime: ~1 min.

## The hypothesis (thread 2/3) and why it was worth testing

The qx+1 family (halve by 2, vary q) has one boundary: `d=ord_q(2)=2 ⟺ q=3`. "5x+1" is not a different system — it's q=5, in-family. The genuine universality test varies the **halving base** `w`: map `(qx+1)·w^{−v}`, weights `p_v=(w−1)w^{−v}` (w=2 = Collatz). This is a well-defined q-adic self-similar measure for any `w` coprime to `q`. The mechanism (worksheet §4, FP1) says the gap closes when the halving phase `w^{−v} mod q` collapses to two values — `ord_q(w)=2 ⟺ w≡−1 mod q`. If universal, it predicts a critical point at **(q=5, w=4)** (`4≡−1 mod 5`), away from q=3.

## Method — direct gap diagnostic, no operator

`value(v₁..v_k) = Σ_m q^{m−1} w^{−S_m} mod q^k` (the "+1" factors out as a unit). Unnormalized cell mass `m_cell = Σ_{addr→cell} ∏ w^{−vᵢ}`; with `U2 = Σ_{v≤V} w^{−2v}`, the normalized correlation is `X_k = (Σ_cell m_cell²)/U2^k`, `Δ_k = X_k − X_{k−1}`. **Gap ⟺ `Δ_k` decays geometrically (`Δ_{k+1}/Δ_k → r < 1`, X_k converges); no gap ⟺ `Δ_k → const` (X_k grows linearly, ratio → 1).**

## Results

| q | w | ord_q(w) | w≡−1? | `Δ₄/Δ₃` | higher-k ratio | class | predicted |
|---|---|---|---|---|---|---|---|
| 3 | 2 | 2 | yes | 1.006 | ~1.0 | **NO-GAP** | no-gap ✓ |
| 5 | 2 | 4 | no | 0.508 | ~0.5–0.6 | GAP | gap ✓ |
| 5 | 3 | 4 | no | 0.272 | ~0.3 | GAP | gap ✓ |
| **5** | **4** | **2** | **yes** | **0.9375** | **0.9375 (k=3..6, flat)** | **GAP** | **no-gap ✗ MISMATCH** |
| 7 | 2 | 3 | no | 0.357 | — | GAP | gap ✓ |
| 7 | 3 | 6 | no | 0.614 | — | GAP | gap ✓ |
| 7 | 6 | 2 | yes | 0.972 | 0.9722 (k=3..6, flat) | GAP (weak) | no-gap ✗ |
| 11 | 2 | 10 | no | 0.527 | — | GAP | gap ✓ |
| 11 | 10 | 2 | yes | 0.990 | — | GAP (weak) | no-gap ✗ |

**The decisive falsifier — (5,4):** pushed to k=6, the ratio is dead-constant at **0.9375** (0.9373, 0.9375, 0.9375, 0.9375 for k=3,4,5,6). A constant ratio `< 1` means `Δ_k ~ 0.9375^k` decays geometrically ⟹ **X_k converges ⟹ a GAP.** This is qualitatively unlike q=3 (ratio ~1, X_k grows linearly). So (5,4), predicted critical by the `ord_q(w)=2` mechanism, is **not** critical — it has a weak gap. **H_UNIV is refuted.**

## Why q=3 is actually special (the real boundary criterion)

`ord_q(w)=2` gives a 2-valued halving phase at BOTH (3,2) and (5,4) — yet (3,2) diverges and (5,4) gaps. So the 2-valued-phase mechanism is insufficient. The distinguisher is **R6's criterion `(q−1)/d = 1`** (conservation *determines* the leading off-diagonal mode ⟹ forces `λ₂=λ₁`), NOT `d=2`:
- (3,2): `d=2`, `(q−1)/d = 1` ⟹ mode determined ⟹ **no gap**.
- (5,4): `d=2`, `(q−1)/d = 2` ⟹ mode underdetermined ⟹ **free to gap** (r=0.9375). ✓
- (5,2): `d=4=q−1`, `(q−1)/d = 1` ⟹ but still GAPS (r≈0.62)! So `(q−1)/d=1` alone is ALSO not sufficient.

So neither `d=2` nor `(q−1)/d=1` alone is the boundary in the general (q,w) family. **q=3 needs both simultaneously — `d=2` AND `(q−1)/d=1` — which forces `q−1=2`, i.e. q=3 uniquely.** For the actual Collatz map (`w=2` fixed), the halving base is order-2 only at q=3, so there is exactly one critical point and no class to generalize. **Universality (as an `ord_q(w)=2` class) does not exist.**

## Bonus observation (flagged, NOT claimed)

The `w=q−1` weak-gap ratios are suspiciously clean: (5,4)→15/16, (7,6)→35/36, (11,10)→99/100 = **`1 − 1/w² = 1 − w^{−d}`** (d=2), matching the leading tower ratio `x₁ = w^{−d}`. This may be the **within-cell tower** geometric mode (which I did not subtract from `Δ_k`), not the true cross-cell gap — so it is an observation, not a result. Notable contrast if real: `r_q` for w=2 (Collatz) has NO closed form (R28), but the `w=q−1` gaps land on a clean `1−1/w²`. Not chased (irrelevant to the universality verdict, which stands regardless — X_k converges either way).

## What this means for the arc

- **Honest negative, per the ethos:** universality was the most seductive of the three threads and the one the arc's history (Collatz walls-not-threads, Solar/Cosmology killed) flagged as mirage-prone. Building the falsifier and running it is the win; the finding is that q=3's criticality is genuinely non-generic.
- **q=3's uniqueness is now sharper:** not "the q with 2-valued halving phase" (that's also (5,4)), but "the unique q where the halving base is order-2 AND primitive." This tightens the phase-boundary statement in the L3 briefs.
- **No change to Result 1 / L3:** those are about the w=2 (Collatz) column, where q=3 is the sole boundary. Untouched.

## Not at stake
R10–R39 (all w=2). A refutation here kills only the universality hypothesis, which is what happened — cleanly.

_Reporting discipline: falsifier pre-registered (`(5,4)` predicted no-gap; a gap there = refutation), and it fired. The automated classifier's absolute AMBIG band `[0.85,0.95]` caught (5,4) at 0.9375 on a technicality — an 8th d/V-blind threshold mis-spec this arc — but the higher-k push (constant 0.9375 across k=3–6) resolves it unambiguously as a gap, and the bimodal cluster separation (gap cases 0.27–0.61 vs no-gap-cluster... which turned out to be (3,2) ALONE) is clean. The `1−1/w²` pattern is disclosed as possibly within-cell, not claimed as the cross gap._
