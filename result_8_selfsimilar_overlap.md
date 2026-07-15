# Result 8 (qx+1 paper) — π_k IS a q-adic self-similar measure. The domination step is an OVERLAP estimate, and the "3" is the address measure's participation ratio.

**Date:** 2026-07-15. **Verdicts: H_ADDR (gate) CONFIRMED at machine zero / H_DIAG CONFIRMED / H_LB CONFIRMED / H_DOM SPLIT — holds q≥5, FAILS at q=3 / H_ORD CONFIRMED-monotone but functional form WRONG (→ R9/R10).**

**Headline: the reframing is valid. `‖π_k‖² = diagonal + overlaps`, the diagonal is `(Σ_v p_v²)^k` with `Σ_v 4^{−v} = 1/3` — so R5's named "3" is the L² dimension of a self-similar measure, and "sub-leading characters don't perturb the rate" restates as "collision mass = O(diagonal)".**

Probe: `probe_8_selfsimilar_overlap.py`. Log: `result_8_overlap_log.txt`.

## The reframing (derived, then gated)

Iterating the chain from any `r_0`:
> `r_k = q^k·r_0·2^{−A_k} + Σ_{m=1}^{k} q^{m−1}·2^{−S_m}`, `S_m = v_{k−m+1}+…+v_k`

**Mod `q^k` the `r_0` term vanishes identically** ⇒ `π_k` is exactly the law of `Σ_{m=1}^k q^{m−1}2^{−S_m}`, independent of the start. (Independently consistent with STATE's K_k lemma: the chain mixes in exactly k steps.)

So `π_k` is a **q-adic self-similar measure**: IFS `T_v(x)=(qx+1)/2^v`, weights `p_v = 2^{−v}/Z`. **Every map contracts by exactly `1/q` q-adically** (`|q(x−y)/2^v|_q = q^{−1}|x−y|_q`; the `2^v` is a unit and does nothing). By address `a=(v_1..v_k)`, `p_a = ∏_i p_{v_i}`:

> `‖π_k‖² = Σ_r (Σ_{a→r} p_a)² = **Σ_a p_a²** [DIAGONAL] + **Σ_{a≠a', val=val'} p_a p_{a'}** [OVERLAPS]`
> `DIAGONAL = (Σ_v p_v²)^k → (Σ_v 4^{−v})^k = (1/3)^k`

**The "3" is `1/Σ_v p_v²` — the address measure's own participation ratio; equivalently `D₂ = log3/log q` is the correlation dimension.** This makes R5's constant a *known species* of object rather than a bespoke identity.

**Truncation is EXACT, not an approximation:** `2^{−v} mod q^k` has period `M=ord_{q^k}(2)`, and `P(v ≡ j mod M) = Σ_{i≥0}2^{−(j+iM)} = 2^{−j}/Z` with `Z=1−2^{−M}` — exactly the weights probe_5/probe_6 use. **But therefore `Σ_v p_v² = (1/3)·(1+2^{−M})/(1−2^{−M}) ≠ 1/3` at small M** (13% off at q=5,k=1). All work below uses the exact value.

## H_ADDR (GATE) — CONFIRMED at machine zero

| q,k | addresses | `max|π_addr − π_iter|` | `|diag_enum − (Σp²)^k|` |
|---|---|---|---|
| 3,2 | 36 | **0.000e+00** | 4.2e−17 |
| 3,3 | 5832 | 5.6e−17 | 7.0e−16 |
| 5,2 | 400 | 1.1e−16 | 1.8e−16 |
| 7,2 | 441 | 2.8e−17 | 1.4e−17 |

Two combinatorially unrelated routes (address sum vs power iteration) agree exactly. **The reframing is valid; H_DIAG confirmed.**

## H_LB — CONFIRMED. `C_q ≥ 1` is FORCED, and it explains previously-unexplained empirics.

`offdiag ≥ 0` at every tested (q,k) — Cauchy–Schwarz on each fiber. So **`δ_q > 0` is a theorem, not an observation**, which is exactly why all 8 measured δ in `result_4` are positive. *(This later becomes a kill shot: R9's best linear fit `δ = 0.82/ord − 0.089` predicts δ<0 for ord>9.2, contradicting this theorem.)*

## ★ `ratio_1 = 0` EXACTLY — structural, not luck

At k=1 the value is `2^{−v_1} mod q` with `v_1 ∈ {1..ord_q(2)}` = exactly one full period ⇒ **the coding is a bijection onto `⟨2⟩`** ⇒ zero collisions. Measured `offdiag = 0.00000e+00` at k=1 for q=5,7,11,13 (and later 8 more primes in R9, all ≤3.3e−16). **This is the seed of the whole pillar-3 mechanism** (R9/R10): the first collision cannot occur until a full period shift, which costs `2^{−ord}`.

## H_DOM — SPLIT. Holds q≥5, FAILS at q=3. ⚠️ AND MY DECISION RULE WAS BUGGY.

**⚠️ Rule failure, recorded rather than retro-fixed.** I pre-registered "bounded iff max step < 1.5". **That rule cannot distinguish bounded from linear growth** — a sequence growing like `k` has steps `(k+1)/k → 1` and sails through. The automated line printed "BOUNDED" for q=3. **It is wrong.**

Reading the actual q=3 numbers:
```
ratio:       1.013, 1.604, 2.069, 2.534, 3.000, 3.466, 3.932   (k=2..8)
differences:      0.591, 0.465, 0.465, 0.466, 0.466, 0.466     <- CONSTANT
```
**`ratio_k ≈ 0.4655·k` — LINEAR. And 0.4655 ≈ 7/15 = 0.4667.**

Corrected verdicts:
- **q≥5 CONFIRMED** — converging: q=5→0.210, q=7→0.361, q=11→0.00209, q=13→0.00068 (steps decaying to 1.01–1.03).
- **q=3 REFUTED** — domination fails at the critical case, by exactly a factor of k.

**Not a problem — it is the known critical behaviour.** R7 independently found `‖π_k‖² ~ (7/15)·k·3^{−k}` at q=3. **Fourth independent sighting of the same phase boundary** (R6's `(q−1)/2=1 ⟺ q=3`; R7's geometric-series divergence; R7's X_k linear growth; this).

**Independent cross-check:** at q=3, `X_k = 1 + ratio_k`. Here `ratio_4 = 2.06865` ⇒ `X_4 = 3.06865`. R7's power-iteration `X_4 = 3.06864623`. Different route, same number.

## H_ORD — monotone CONFIRMED, but the functional form is wrong (→ R9/R10)

`ratio/2^{−ord}` spans 1.4×; `ratio·ord` spans 130×. Pillar 3's `δ_q ≈ 0.82/ord_q(2)` is monotone-but-misfit. **Handed to R9 (bake-off) and R10 (conclusive).** Outcome: law A refuted by 2.55e13×; **`δ_q = 2^{1−ord_q(2)}·(q−3)/q`** replaces it.

## Not at stake
THEOREM_C_745, Th 78.1–78.3, R81b, ε_k. R5's rate result stands (q≥5 domination confirmed; q=3's factor-of-k is the critical case, already known from R7).

_Reporting discipline: the gate ran first and passed at machine zero. The H_DOM decision-rule failure is recorded as a failure, with the wrong automated verdict shown, rather than the rule being quietly rewritten. The q=3 refutation is reported as a refutation even though it is benign. H_ORD's monotone "confirmation" is explicitly flagged as insufficient — monotonicity does not identify a functional form — which is what routed R9/R10._
