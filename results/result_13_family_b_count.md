# Result 13 (qx+1 paper) — the family-(b) collision set is EXACTLY characterized. `ε_q` reduces to a CHARACTER SUM — which is why it is erratic.

**Date:** 2026-07-15. **Verdicts: H_SPLIT CONFIRMED (exact algebra) / ★ H_COND CONFIRMED EXACT (zero false pos AND zero false neg, ~450k cross pairs, 6 primes) / H_J1 CONFIRMED at q=41,47 — FAILED my threshold at q=31 (4th mis-specified rule) / H_EPS EXPLORATORY, no fit (pre-committed).**

**Headline: `2^{−A'} ≡ 2^{−A} + j·s·2^{−v_2} (mod q)` with `s = (2^d−1)/q mod q` is an EXACT iff for family-(b) collisions. Expanding its `⟨2⟩`-membership indicator over characters shows `ε_q` IS A CHARACTER SUM — the erraticness is structural, not noise, and the tool is standard.**

Probe: `probe_13_family_b_count.py`. Log: `result_13_family_b_log.txt`. Runtime: **0.54 s**.

## Setup — the cell collapse (from R11's exact H_V1MOD)

`value(v_1,v_2) = 2^{−v_2} + q·2^{−(v_1+v_2)} mod q²` depends on `v_1` **only mod `d`** (`d = ord_q(2)`). So the `M²` addresses collapse to `M·d` **cells** `(v_2, c)`, `c = (v_1−1) mod d`, with `mass[c,v_2] = G_c·p_{v_2}`. For fixed `v_2` the `d` cells have **distinct** values (`2^{−(c+v_2)} mod q` takes `d` distinct values), so **each `v_2` appears at most once per value-bucket**. Hence:

- **family (a)** = within-cell pairs (same `v_2`, same `c`) — counted exactly in R11
- **family (b)** = cross-cell pairs sharing a value — **necessarily different `v_2`**

**H_SPLIT (gate) CONFIRMED — exact algebra**, no approximation:
`offdiag_total = Σ_cells mass² + crossterms − P2²`, `family(a) = Σ_cells mass² − P2²` ⇒ **`family(b) = crossterms` exactly.** Verified numerically (rel err 1e−13 → 1e−7, float accumulation only).

## ★ H_COND — the derivation, and it is EXACT

**Derivation.** For `v'_2 = v_2 + j·d`, write `2^d = 1 + q·s`. Then mod `q²`: `2^{−jd} ≡ 1 − j·q·s`, so `2^{−v'_2} − 2^{−v_2} ≡ −j·q·s·2^{−v_2}`. The second term carries a factor `q` and needs `A = v_1+v_2` only mod `d`. Collision mod `q²`, divided by `q`:

> ## `2^{−A'} ≡ 2^{−A} + j·s·2^{−v_2} (mod q)`, &nbsp; `s = (2^d − 1)/q mod q`

**Tested as an exact iff on enumerated pairs — zero false positives, zero false negatives, at every prime:**

| q | d | `s` | cross pairs | condition violations | non-multiple-of-d shifts | |
|---|---|---|---|---|---|---|
| 11 | 10 | 5 | 4,950 | **0** | **0** | EXACT |
| 13 | 12 | 3 | 10,296 | **0** | **0** | EXACT |
| 17 | 8 | 15 | 3,808 | **0** | **0** | EXACT |
| 31 | 5 | 1 | 1,550 | **0** | **0** | EXACT |
| 41 | 20 | 32 | 155,800 | **0** | **0** | EXACT |
| 47 | 23 | 22 | 273,493 | **0** | **0** | EXACT |

**The family-(b) collision set is now exactly characterized.** Solvable for `A'` **iff the RHS lands in `⟨2⟩`** (the order-`d` subgroup) — a genuinely arithmetic condition.

## H_J1 — confirmed at large q; ⚠️ my threshold failed at q=31

| q | d | \|j\|=1 share | \|j\|=2 | \|j\|≥3 |
|---|---|---|---|---|
| 11 | 10 | 0.996204 | 0.003793 | 0.000003 |
| 13 | 12 | 0.999985 | 0.000015 | 0.000000 |
| 17 | 8 | 0.984584 | 0.015384 | 0.000032 |
| **31** | **5** | **0.886964** | **0.110870** | 0.002166 |
| 41 | 20 | 0.999996 | 0.000004 | 0.000000 |
| 47 | 23 | **1.000000** | 0.000000 | 0.000000 |

**⚠️ Fourth mis-specified decision rule this arc.** I pre-registered "|j|=1 share ≥ 0.90 for q ≥ 31"; q=31 gives **0.886964** and the rule **fired as failed**. Reported as fired. The rule was `d`-blind: q=31 has the **smallest `d` (=5)**, so `x = 2^{−d} = 2^{−5}` is largest and `|j|=2` (weight `~x²/x = x`) carries 11%. A `d`-dependent threshold was the right specification. Mechanistically the picture is intact — `|j|=1` dominance is governed by `x`, and it is total (1.000000) at q=47. *(Prior rule failures: "bounded iff step<1.5" couldn't see linear growth (R8); ΔR² couldn't separate free fits (R9); a relative tolerance against machine-epsilon noise (R11).)*

## H_EPS — EXPLORATORY. No fit. Pre-committed and honored.

```
q:      11       13       17       31       41       47
eps:  0.0582   0.3822   0.0967   0.0068   0.0066   0.0071
```

**I pre-committed to NOT fitting a law to six points** — that is exactly how pillar 3 went wrong — and I am honoring it. For the record, no obvious candidate tracks it:

| q | d | `ε_q` | `d/q` | `1/(q−1)` | `d/(q−1)` | `(q−1)/d` |
|---|---|---|---|---|---|---|
| 11 | 10 | 0.058185 | 0.909 | 0.100 | 1.000 | 1 |
| 13 | 12 | 0.382195 | 0.923 | 0.083 | 1.000 | 1 |
| 17 | 8 | 0.096711 | 0.471 | 0.063 | 0.500 | 2 |
| 31 | 5 | 0.006819 | 0.161 | 0.033 | 0.167 | 6 |
| 41 | 20 | 0.006595 | 0.488 | 0.025 | 0.500 | 2 |
| 47 | 23 | 0.007107 | 0.489 | 0.022 | 0.500 | 2 |

`d/(q−1)` is flat at 0.5 across q=17, 41, 47 while `ε_q` swings **14×** (0.0967 vs 0.0066). **q=17 and q=41 have identical `(q−1)/d = 2` and `ε_q` differing by 14×.** No smooth function of `d` or `q` can do this.

## ★★ The route — and why erratic is the *expected* answer

With H_COND exact, family (b) is a weighted sum over `(v_1,v_2,j)` of the indicator `1[2^{−A} + j·s·2^{−v_2} ∈ ⟨2⟩]`. Expand the indicator over characters:

> `1[y ∈ ⟨2⟩] = (d/(q−1))·Σ_{χ : χ|_{⟨2⟩} = 1} χ(y)` — a sum over the `(q−1)/d` characters trivial on `⟨2⟩`

> ## ⇒ **`ε_q` is a CHARACTER SUM.**

**This explains the erraticness structurally rather than excusing it.** It is the same species as a Jacobi/Kloosterman sum: the value depends on `q`'s specific arithmetic, not on any smooth function of `d` or `q`. It also means the correct tool is **standard character-sum machinery with square-root cancellation** — not a fit, and not a new theory. The principal-character term gives the "expected" `d/(q−1)` density; the non-principal characters give the fluctuation that makes q=17 and q=41 differ 14× at identical `d/(q−1)`.

**This is the honest next step**, and it is a derivation, not a measurement.

## Where pillar 3 stands after R11 + R13

| piece | status |
|---|---|
| `v_1` determined only mod `d` ⇒ geometric tower | **THEOREM** (R11, exact) |
| leading term `(1+x)/(1−x)·(1−2^{−M})/(1+2^{−M}) − 1 → 2^{1−ord_q(2)}` | **THEOREM** (R11, exact-rational verified) |
| prefactor 2 = tower cross-term | **DERIVED** (R11) |
| `δ_q = ratio_2·(q−3)/q` | **exact** (R10) |
| family-(b) collision set | **EXACTLY CHARACTERIZED** (R13, this) |
| `ε_q` value | **OPEN — reduces to a character sum** (R13) |

## Not at stake
R10's law, R11's family-(a) identity, R5's rate, R6, R7, R12, THEOREM_C_745, Th 78.1–78.3, R81b, ε_k.

_Reporting discipline: H_COND was derived on paper first, then tested as an exact set-equality (both directions) rather than a tolerance — it cannot be gamed by a threshold, and it passed at ~450k pairs. H_J1's rule fired FAILED at q=31 and is reported as fired, with the mis-specification named (d-blind) rather than the threshold widened — fourth such failure this arc. H_EPS carried a pre-committed refusal to fit, and no law is proposed; the character-sum reduction is offered as a ROUTE, not a result. Author's structural priors this arc: 9-for-16._
