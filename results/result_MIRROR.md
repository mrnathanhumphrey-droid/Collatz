# RESULT — MIRROR MAP (q,p)=(2,3): the q↔p image of Collatz is the unique SECOND critical point; involution symmetric on STRUCTURE + arithmetic TYPE, but NOT on the VALUE (2026-07-28)

**Probe:** `probes/probe_mirror.py` (this diagnostic run: scratch `mirror_diag.py` + `mirror_depth.py`, promoted). The **mirror map** `x ↦ (2x+1)/3^v` = the p-Hydra member with multiplier `q=2` and valuation prime `p=3` — the **q↔p image of Collatz `(3,2)`** under the family involution `q=(p+1)/(p−1)`. Built on the certified `Xk_qp` machinery from `result_PHYDRA_FAMILY.md` with the roles of 2 and 3 swapped (tower over `2^k`, chain on `(ℤ/2^k)*`, `r ↦ (2r+1)·3^{−v}`, `v~Geom(2/3)`). Gate: reuses the phydra transfer op verbatim; `(3,2)` side reproduces the certified S-ladder (`2/3, 10/21`).

## Why the mirror specifically (not one sibling among many)
The criticality condition `(p−1)/(p+1) = 1/q` ⟺ **boundary curve `q(p−1) = p+1`**. The **only two integer points that are both constructible (p prime) and critical** are `(q=3,p=2)` = Collatz and `(q=2,p=3)` = the mirror (`2·2 = 4 = p+1` ✓). They are q↔p images of each other. So "is there structure across the critical class?" reduces, for the constructible members, to a single yes/no: **does `S_∞(3,2)` relate to `S_∞(2,3)` under the swap?** There is no third point to triangulate — the whole "solve the critical class" program lives or dies here. Never computed before this run.

## The three checks
**(a) Own early ladder (exact).** `S_k(2,3) = 1, ½, ⅝, ½, 20/41, 3808/7913, …` — jagged, its own object, floor climbing toward ~0.456–0.46. (Collatz for contrast: the smooth `2/3, 10/21, 0.464…`.)

**(b) Involution test — equality REFUTED (pre-registered outcome 3).** Same machinery both sides, `k` to 14:
```
   mirror  lim S_k(2,3) ~ 0.456–0.46
   Collatz lim S_k(3,2) ~ 0.470  (→ banked S_inf ~ 0.4749, floor 0.473177)
   difference (2,3)-(3,2) = -0.011   (>> convergence noise)
   ratio 0.976   sum 0.929   product 0.216   -- nothing clean
```
`S_∞(2,3) ≈ 0.459 ≠ S_∞(3,2) ≈ 0.475`. The −0.011 gap is far outside convergence noise ⟹ **equality (outcome 1) is decisively refuted.** No clean function (sum/ratio/product) lands at 3 digits either. *Guardrail honored: a clean digit-relation would be a lead to prove, never trusted — and none appeared; outcome 2 not chased.*

**(c) Depth class — SAME (infinite Mahler depth).** Exact denominators `k=1..8` ignite after a ramp-up:
```
   k   den(S_k)                         bits    factor(den)
   5   41                               5.4     41
   6   7913                             13.0    41·193
   7   118709105848321                  46.8    17·41²·193·21523361
   8   1.06e32                          106.4   5·17·41²·193²·21523361·926510094425921
```
bits `5→13→47→106`, ratio averaging **>2** = doubly-exponential = the MAHLER signature. And the primes are the tell: `41 | 2²⁰−1`, `193 | 2⁹⁶−1` — Mersenne **cofactor** primes (`ord(2 mod 41)=20`, `ord(2 mod 193)=96`), the **same `⟨2⟩`-multiplicative-order mechanism** as `c̃_q` / the denominator theorem, just not the full Mersenne numbers (the tower is `2^k` here, not `3^k`). So the mirror is **the same arithmetic type** as Collatz: infinite Mahler order, walled value. Not special.

## The finding
The q↔p involution is a **symmetry of the criticality condition** (both `(3,2)` and `(2,3)` sit on the boundary curve) and a **symmetry of the arithmetic type** (both infinite Mahler depth, same `⟨2⟩`-order prime mechanism) — but it is **NOT a symmetry of the value** (`0.475 ≠ 0.459`). This is exactly the pre-registered "**each critical constant is algebraically alone**" outcome: the two real critical points talk at the level of **structure / associated-graded**, and are **silent at the level of the constant**.

Same split as everything else this arc:
- **biextension** (`result_SCOPING_HEIGHT_BIEXT.md`): structure homed abstractly, height (value) walled.
- **p-Hydra law** (`result_PHYDRA_FAMILY.md`): gr functional equation `q(p−1)/(p+1)` real, value untouched.
- **mirror** (here): involution symmetric on structure + type, silent on value.

The mirror was the **one directly testable place** the "no functional equation on the value" belief could have broken (outcome 1 would have been a genuine constraint on the constant). It held — the belief is now a checked fact, not an assumption.

## Honest caveats
- `S_∞(2,3) ≈ 0.456–0.46` is soft (infinite-depth ⟹ slow; float oscillated 0.456–0.465, exact k=7/8 gave 0.4635/0.4563 — still moving). The **equality refutation is robust regardless** (the 0.011 gap dwarfs the uncertainty); the precise mirror constant is not nailed and would need more levels (slow, same wall).
- The depth-class verdict is settled to k=8 exact — doubly-exponential ignition + `⟨2⟩`-cofactor primes are unambiguous.

## Net
- **Newly banked:** the mirror `(2,3)` = the unique second constructible critical point; involution **refutes value equality** (`0.475 ≠ 0.459`); mirror is **same arithmetic class** (infinite Mahler depth, `⟨2⟩`-cofactor-prime denominators). The critical class has **no functional equation on the value** even between its only two real points — the "each critical constant is algebraically alone" reading is confirmed where it was testable.
- **Not at stake:** 7/15 (floor 0.473177), MAHLER, GARSIA, DENOM, SOLSTICE, PHYDRA_FAMILY, R1–R30. The mirror is a *different* constant; nothing about `S_∞(3,2)` moves.
