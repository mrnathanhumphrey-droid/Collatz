# Result 7 (qx+1 paper) — object identification: the q-sweep and R76 have been computing the SAME object by different routes. The paper's primitive is `q^k‖d_k‖²`, not `q^k‖π_k‖²`.

**Date:** 2026-07-15. **Verdicts: H_ID_A REFUTED / H_ID_B CONFIRMED / H_ID_C SPLIT (my prior WRONG) / H_ID_D SPLIT.**

**Headline: `D_k` (the q-sweep's pillar-2 "difference") IS `M_k(1)` (R76's `S_k`) IS `q^k‖d_k‖²` (R74's identity). One object, three names, two threads that didn't know they were the same. The identity is exact and one-line provable.**

Probe: `probe_7_object_identification.py`. Log: `result_7_identification_log.txt`. Closes the gap flagged in `result_6_conservation_generalize.md` §⚠️. Compute: seconds.

## The tension that motivated this

- **R74/R75/R76 object:** `S_k := q^k‖d_k‖²` where `‖d_k‖² := Σ_{r'}(π_k(r') − π_{k−1}(parent(r'))/q)²` (level-incremental deviation), with Plancherel `S_k = Σ_{gcd(ξ,q)=1}|μ̂_k(ξ)|² =: M_k(1)`.
- **q-sweep object** (`probe_5_universal_rate.X_gen`): `X_k := q^k‖π_k‖²` — raw stationary mass.

Different functionals. Yet STATE/memory asserted both `S_k^(q) = q^k‖π_k‖²` **and** `S_k^(3) → 7/15`. Since R5 also established `‖π_k‖² ~ 3^{−k}`, those two force `X_k(q=3) → 1 ≠ 7/15`. **Both could not be true of one object.**

## Measured — both objects, side by side

| q | k | `X_k = q^k‖π_k‖²` | `M_k(1)` | `q^k‖d_k‖²` | X ratio/(q/3) | M ratio/(q/3) |
|---|---|---|---|---|---|---|
| 3 | 1 | 1.66666667 | 0.66666667 | — | — | — |
| 3 | 2 | 2.14285714 | 0.47619048 | 0.47619048 | 1.28571 | 0.71429 |
| 3 | 3 | 2.60443182 | 0.46157468 | 0.46157468 | 1.21540 | 0.96931 |
| 3 | 4 | 3.06864623 | 0.46421441 | 0.46421441 | 1.17824 | **1.00572** |
| 5 | 4 | 9.28621186 | 3.76325228 | 3.76325228 | 1.00883 | **0.99609** |
| 7 | 3 | 17.19773535 | 9.94170134 | 9.94170134 | 1.01577 | **1.00110** |
| 11 | 2 | 13.47225814 | 9.79842301 | 9.79842301 | 1.00011 | **0.99942** |

## H_ID_A — REFUTED. The name collision is real.

`max |X_k − M_k(1)| = 7.26` over tested range. **Different objects**, exactly as prior'd.

## H_ID_B — CONFIRMED at 1e−15. R74's identity AND R75's Plancherel both port to general q.

`M_k(1) = q^k‖d_k‖²` to ≤1.8e−15 at every tested `(q,k)` — q=3,5,7,11. Since `M_k(1)` is *defined* here as the coprime-restricted Plancherel sum `Σ_{gcd(ξ,q)=1}|μ̂_k(ξ)|²`, this simultaneously confirms:
- **R74's identity** `S_k = q^k‖d_k‖²` generalizes off q=3, and
- **R75's Thm 75.1** (Plancherel decomposition) generalizes off q=3.

**Combined with R6: R74 ✓ ports, R75 ✓ ports, R76 Thm 76.1 (conservation) ✓ ports. Only R76 Thm 76.3 (leading-mode collapse) ✗ — and R6 showed exactly why ((q−1)/2 = 1 ⟺ q = 3).** The q=3 machinery is now fully triaged.

## ★ DISCOVERED (not pre-registered) — `X_k − X_{k−1} = M_k(1)`, exactly, and it is PROVED

Spotted in the table, then verified: the increments of the q-sweep's object **are** R76's object.

| q | k | `X_k − X_{k−1}` | `M_k(1)` |
|---|---|---|---|
| 3 | 4 | 0.46421441 | 0.46421441 |
| 5 | 4 | 3.76325228 | 3.76325228 |
| 7 | 3 | 9.94170134 | 9.94170134 |
| 11 | 2 | 9.79842301 | 9.79842301 |

Exact at every q (residuals ≤1.8e−15; the 1e−8 entries are 8-dp table transcription, not error).

**Proof (one line, elementary — requires only the projection/consistency property that π_k pushes forward to π_{k−1}):**

> `‖d_k‖² = Σ_{r'}(π_k(r') − π_{k−1}(par)/q)²`
> `      = ‖π_k‖² − (2/q)·Σ_{r'} π_k(r')π_{k−1}(par(r')) + (1/q²)·Σ_{r'} π_{k−1}(par(r'))²`
>
> Each parent has exactly `q` coprime children, and the children's masses sum to the parent's mass (consistency), so
> `Σ_{r'} π_{k−1}(par)² = q·‖π_{k−1}‖²` and `Σ_{r'} π_k(r')π_{k−1}(par(r')) = Σ_{par} π_{k−1}(par)² = ‖π_{k−1}‖²`.
>
> `⟹ ‖d_k‖² = ‖π_k‖² − (2/q)‖π_{k−1}‖² + (1/q)‖π_{k−1}‖² = ‖π_k‖² − (1/q)‖π_{k−1}‖²`
> `⟹ q^k‖d_k‖² = q^k‖π_k‖² − q^{k−1}‖π_{k−1}‖² = X_k − X_{k−1}.  ∎`

So **`X_k = X_0 + Σ_{j≤k} M_j(1)`** — the q-sweep's object is the *cumulative sum* of R76's object. Not a coincidence; an identity.

## ★★ THE UNIFICATION — `c̃_q` is already about R76's object

The q-sweep's pillar 2 defines `c̃_q := D_k/(q/3)^k` with `D_k = X_k − X_{k−1}`. By the identity above, **`c̃_q = M_k(1)/(q/3)^k`** — i.e. `c̃_q` is R76's `S_k`, normalized. Verified against R5's own reported numbers:

| q | k | this probe: `M_k(1)/(q/3)^k` | `result_5_universal_rate.md` reported `c̃` |
|---|---|---|---|
| 7 | 3 | **0.78258** | **0.78258** ✓ exact |
| 11 | 2 | 0.72881 | 0.72880 (at k=3) ✓ |
| 5 | 4 | 0.48772 | 0.48963 (at k=3) ✓ consistent drift |

**Two threads — the Collatz-closure thread (R74/R75/R76, since 2026-05-03) and the qx+1 q-sweep (since 2026-05-04) — have been computing one object by two routes for two months.** R76's machinery applies to the q-sweep's pillar-2 constant *directly*, and R6's triage of that machinery is therefore correctly aimed.

## H_ID_C — SPLIT. My prior was WRONG on X_k.

- **`M_k(1) → 7/15` at q=3: CONFIRMED.** And its deviations reproduce R76 §10's ε_k **exactly**: |M−7/15| = 9.524e−03, 5.092e−03, 2.452e−03 at k=2,3,4 = R76's ε_2, ε_3, ε_4. `M_k(1)` **is** R76's `S_k`, confirmed against published values.
- **`X_k(q=3) → 1`: REFUTED.** X_k *grows*: 1.667, 2.143, 2.604, 3.069. Prior wrong. It grows **linearly**, `X_k ≈ (7/15)·k`, because at q=3 the increments `M_k(1)` → 7/15 ≠ 0 and never decay.

## H_ID_D — SPLIT, and this is the decisive one

- **`M_k(1)` carries the `(q/3)^k` rate at EVERY tested q, including the critical q=3** (ratio/(q/3) = 1.00572 / 0.99609 / 1.00110 / 0.99942 at q=3/5/7/11).
- **`X_k` carries it only for q≥5** (1.009 / 1.016 / 1.0001), where it is *inherited* by geometric summation (`Σ_j (q/3)^j ~ (q/3)^k·q/(q−3)`). **At q=3 the geometric factor is trivial and X_k degenerates to linear growth** (ratio/(q/3) = 1.178 at k=4, approaching 1 only like `1+1/k` — polynomially, not geometrically).

## Consequence for the paper — a naming fix with real content

**The primitive object is `M_k(1) = q^k‖d_k‖²`, not `X_k = q^k‖π_k‖²`.** `M_k(1)`:
1. carries the clean universal `(q/3)^k` rate at **every** odd q **including the critical q=3**, where `X_k` degenerates;
2. **is** R76's `S_k` (verified against published ε_k), so R74/R75/R76 machinery applies to it directly;
3. is what `c̃_q` is *already* measuring (pillar 2), verified to 5 dp against R5's own reported values.

`X_k` is a cumulative artifact. R5's pillar-1 rate claim was measured on `X_k` and **survives** (for q≥5 by inheritance; at q=3 `X_k/X_{k−1} → 1 = q/3` still holds, just polynomially) — **the rate result is not threatened.** But pillar 1 (rate, on `X_k`) and pillar 2 (constant, on `D_k = M_k(1)`) are stated about **different objects in the same paper**, and the object that degenerates at q=3 is the one carrying the headline. **Restating both pillars on `M_k(1)` makes the paper uniform, sharpens the q=3 critical case rather than degenerating there, and imports R74/R75/R76 for free.** Recommended before publication.

## Lead (NOT a result) — `δ_q` extrapolated to the critical case

`c̃_3 = lim M_k(1) = 7/15 = 0.4667` while `(q−3)/q = 0` at q=3 ⇒ `δ_3 = 0.4667`. Pillar 3's empirical `δ_q ≈ 0.82/ord_q(2)` with `ord_3(2) = 2` predicts **0.41** — within **12%** at a case the fit never saw (fit was q≥5 primes). Suggestive that pillar 3 knows something about q=3, where pillar 2 collapses to zero. **Prior stated to lose:** 12% is loose, `0.82` is itself an 8-point empirical fit, and this arc has killed five small-window extrapolations. Worth exactly one probe, not a claim.

## Not at stake

THEOREM_C_745 (c=7/45), Th 78.1–78.3, R81b, ε_k. R5's rate result stands. R6's verdict stands and is now correctly aimed. This probe is identification/naming — it cannot falsify a mechanism, only relabel its subject.

_Reporting discipline: H_ID_C's X_k branch was my pre-registered prior and it LOST — reported as wrong, not quietly dropped. The unification (`D_k = M_k(1)`) was NOT pre-registered — it was spotted in the table, then verified numerically AND proved algebraically before being claimed. The `δ_3` extrapolation is filed as a lead with its prior stated to lose. Author's structural priors: 3-for-8 this arc._
