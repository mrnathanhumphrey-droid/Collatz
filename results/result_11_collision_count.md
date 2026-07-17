# Result 11 (qx+1 paper) — the collision count: pillar 3's leading term is now an EXACT IDENTITY. The prefactor 2 is derived — and my stated reason for it was wrong.

**Date:** 2026-07-15. **Verdicts: H_V1MOD CONFIRMED (exact) / H_FAMA CONFIRMED as an EXACT IDENTITY (after two errors of mine, below) / H_DOM_A CONFIRMED (99.3%) / H_RESID REFUTED — my prior lost.**

**Headline: `ratio_(a) = (1+x)/(1−x)·(1−2^{−M})/(1+2^{−M}) − 1` with `x = 2^{−ord_q(2)}`, verified by EXACT RATIONAL ARITHMETIC (Fraction equality, not tolerance). Leading order `2^{1−ord_q(2)}`. Pillar 3's leading term is a theorem; the ~0.7% residual is open and is NOT `O(1/q)`.**

Probe: `probe_11_collision_count.py`. Log: `result_11_collision_count_log.txt`. Runtime: **0.18 s**.

## H_V1MOD — the structural fact, EXACT

`value(v_1,v_2) = 2^{−v_2} + q·2^{−(v_1+v_2)} mod q²`. The second term carries a factor `q`, so it needs `A = v_1+v_2` only **mod `d = ord_q(2)`** (since `2^d ≡ 1 mod q`). Therefore:

> **`v_1 → v_1 + d` leaves the value EXACTLY unchanged. `v_1` is only ever determined mod `d`.**

Verified by exact integer equality at q=11, 13, 17, 31, 41, 47 (every row, all M−d of them).

**This overturns the R9/R10 story of "the cheapest collision costs a full period shift."** Collisions are not *cheap* — they are **structural and always present**: every value-bucket contains a whole **geometric tower** in `v_1`.

## H_FAMA — the count. The prefactor 2 falls out of the tower's cross-term.

**Family (a)** := pairs with the same `v_2` and `v_1 ≡ v'_1 (mod d)`. With `x = 2^{−d}`, `tM = 2^{−M}`, `M = ord_{q²}(2)` (`= d·q` at every tested prime, verified):

> `G_c := Σ_{v≡c (d)} p_v = 2^{−c}/(1−x)` &nbsp;(the `Z` cancels exactly)
> `H_c := Σ_{v≡c (d)} p_v² = 4^{−c}(1+tM)/((1−x²)(1−tM))`
> bucket offdiag `= p_{v_2}²·(G_c² − H_c)`

Summing `Σ_{c=1}^{d} 4^{−c} = (1−x²)/3` and `P2 := Σ_v p_v² = (1/3)(1+tM)/(1−tM)`:

> ## `ratio_(a) = (1+x)/(1−x) · (1−tM)/(1+tM) − 1` &nbsp;&nbsp;**[EXACT]**
> → `2x/(1−x)` → **`2^{1−ord_q(2)}`** as `M → ∞`

**★ The prefactor 2 is the cross-term of the geometric tower** — `(1+x)/(1−x) − 1 = 2x/(1−x)` — **not** ordered-pair double-counting. **R9/R10's stated explanation for the 2 was WRONG**; it matched the number for the wrong reason. Corrected here.

### ⚠️ Two errors of mine, recorded

**(1) My decision rule failed — third time this arc.** I pre-registered "H_FAMA CONFIRMED iff |rel err| < 1e−9". It fired **REFUTED** at 1.222e−09 (q=47). The rule was mis-specified: I put a *relative* tolerance on a quantity that shrinks to 2.4e−7 while float noise stays pinned at machine epsilon. The **absolute** errors were 2.9e−16, 2.7e−17, 1.4e−17 — all noise floor. *(Prior rule failures: "bounded iff step<1.5" couldn't see linear growth (R8); ΔR² couldn't separate free fits (R9).)*

**(2) Exact arithmetic caught a real algebra error that float hid.** Fraction equality against my *first* closed form `2x(1−tM)/((1−x)(1+tM))` returned **False** at all six primes — while floats agreed to **18 significant digits**. The discrepancy is exactly `−2·2^{−M}/(1+2^{−M})` (≈1e−33 at q=11, ≈1e−247 at q=41), confirmed by exact match at q=11, 17, 31, 41. **The corrected form gives exact Fraction equality at every prime.** The leading `2x/(1−x)` and the derived prefactor 2 are unaffected — but I would not have found the error in float.

**Method note: run the exact-arithmetic check even when float agrees to 18 digits.**

## H_DOM_A — CONFIRMED

Family (a) share of total offdiag: **0.99345 (q=41), 0.99294 (q=47), 0.99323 (q=31)**.

## H_RESID — REFUTED. My prior lost.

Family (b) (different `v_2`, i.e. `v'_2 = v_2 + jd` with a compensating `v_1`) was predicted to be the `O(1/q)` correction, shrinking with q. **It does not shrink:**

| q | d | resid/family(a) |
|---|---|---|
| 11 | 10 | 0.058185 |
| 13 | 12 | **0.382195** |
| 17 | 8 | 0.096711 |
| 31 | 5 | 0.006819 |
| 41 | 20 | 0.006595 |
| 47 | 23 | 0.007107 |

Positive at every q (consistent with H_LB), **erratic at small q, then flat at ~0.7% across q=31→47.** A plateau, not a decay. Mechanically plausible: family (b) is *also* `O(x)` (its pairs cost `2^{−jd} = x^j`), so `resid/family(a) ~ const` rather than `~1/q`.

**⇒ R9/R10's `(1 + O(1/q))` framing is UNSUPPORTED and should be restated as `(1 + ε_q)` with `ε_q ≈ 0.007` for q ≥ 31, erratic below.** Honest limit: three large-q points cannot distinguish a plateau from slow decay. Either way my stated prior does not survive.

## Where pillar 3 now stands

- **Leading term: THEOREM.** `ratio_(a) = (1+x)/(1−x)·(1−tM)/(1+tM) − 1` is an exact identity, derived and verified in exact rational arithmetic. `→ 2^{1−ord_q(2)}`.
- **Structural fact: THEOREM.** `v_1` determined only mod `d` ⇒ geometric tower ⇒ the 2.
- **Residual (family b): OPEN.** ~0.7% at large q, erratic at small q, **not** `O(1/q)`. Requires the family-(b) count — the compensating-`v_1` solvability condition `2^{−A'} ≡ 2^{−A} + j·t·2^{−v_2} (mod q)`, which is solvable iff the RHS lands in `⟨2⟩`. That is a genuinely arithmetic condition (why small q is erratic) and is the next count.
- **`δ_q = ratio_2·(q−3)/q`** unchanged (R10, exact at q=41, 47).

## Not at stake
THEOREM_C_745, Th 78.1–78.3, R81b, ε_k, R5's rate, R6, R7. R10's **law** stands (measured conclusively); this probe corrects R10's **explanation** of the prefactor and refutes its `O(1/q)` framing.

_Reporting discipline: the decision rule fired REFUTED and is reported as having fired REFUTED, with the mis-specification named rather than the threshold quietly widened. The exact-arithmetic check found an error in my own algebra that float agreed with to 18 digits — reported, located exactly, and corrected. H_RESID was my prior and it LOST; the O(1/q) story it inherited from R9/R10 is withdrawn. The prefactor 2 is now derived, and my previous published reason for it (ordered-pair double-count) is marked WRONG. Author's structural priors this arc: 7-for-13._
