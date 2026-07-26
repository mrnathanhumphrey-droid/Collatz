# RESULT — PROBE PERSIST: level-1 dichotomy proof verified + persistence-bound targets (2026-07-26)

**Probe:** `probe_persist.py`. Gates Wilson's level-1 origin of the enriched/depleted dichotomy and measures the
concrete persistence-bound numbers for the pen inequality.

## Level-1 proof — VERIFIED EXACT (Fractions)
- `ord_9(4) = 3` ✓ (4, 16≡7, 64≡1 mod 9) ⟹ `3∣k ⟺ 4^k≡1 mod 9` (ratio trivial at the 2nd tower level).
- `ν_1 (dlog) = (0, 1/3, 2/3)`; `Σν_1² = 5/9` ✓, cross `Σρ(y)ρ(y±1) = 2/9` ✓.
- `γ_1(k) = 3·factor = 5/3 (3∣k) / 2/3 (else)` ✓ — matches the banked seed `γ_1=[5/3,2/3,2/3,5/3]`.
- **Enrichment = Cauchy–Schwarz**: `Σν_1² ≥ Σν_1(y)ν_1(y+c)`, strict since ν_1 is not shift-invariant. The dichotomy
  is born at level 1, exactly.

## Persistence as a two-sided bound (not a sign) — both channels comfortably inside

`γ_r(k) = γ_1(k)·Π_{j=2}^r 3q_j(k)`, `q_j(k)=p_j(k)/p_{j-1}(k)`. Survives to ∞ iff `Π 3q_j(1) < 3/2` (depleted stays
<1) and `Π 3q_j(3) > 3/5` (enriched stays >1):

| k | type | γ_1 → γ_16 | Π 3q_j (to r=16) | threshold | signed Σ(q_j−⅓) | total Σ\|q_j−⅓\| |
|---|------|-----------|------------------|-----------|-----------------|------------------|
| 1 | dep  | 0.667→0.730 | 1.09502 | <1.5 ✓ | +0.0306 | 0.0306 (all +, monotone) |
| 2 | dep  | 0.667→0.473 | 0.71011 | <1.5 ✓ | −0.0996 | 0.1086 (mixed, oscillating) |
| 3 | ENR  | 1.667→1.237 | 0.74224 | >0.6 ✓ | −0.0945 | 0.0945 (falls to white) |

**Room to threshold (the pen's target):**
- k=1 depleted: `log Π 3q_j = +0.0908` vs ceiling `log(3/2)=+0.4055` → **room +0.315** (tail beyond r=16 must not add
  more; the k=1 relaxation is +0.03 total and decaying ~0.9/level ⟹ tail ~0.003, huge margin).
- k=3 enriched: `log Π 3q_j = −0.2981` vs floor `log(3/5)=−0.5108` → **room +0.213** (tail must not subtract more;
  again small and decaying).

The relaxation `Σ|q_j(k)−1/3|` is a small, convergent total (0.03–0.11) — a **bounded** quantity, not a sign.

## The type-match (Wilson's key point, confirmed)
The persistence target is an **upper bound on `Σ_j|q_j(k)−1/3|`** — relaxation-to-white, not sign. That is exactly the
type the "dead" decay shelves supply: Tao Prop 1.14/1.17, BGK, Heilbronn, Bourgain give uniform `|ν̂|` UPPER bounds,
which cancel in a normalized ratio (why they couldn't touch the sign) but directly bound total relaxation. **First
genuine type-match in the arc.** Hank's re-task: bound `Σ_j|q_j(k)−1/3|` for the Syracuse measure.

## Net
- Dichotomy's **base case is proved exact** (level-1, Cauchy–Schwarz); **persistence is a two-sided bounded-relaxation
  inequality** with large margins (rooms +0.315, +0.213), **separable from the rate/turnover question entirely**.
- The decay literature is now the **right type** for the target. The dichotomy is the arc's best-shaped pen target.
- k=2 caution confirmed: it moves *away* from white (0.667→0.473, signed Σ negative) — the bound is genuinely
  two-sided (depleted channels need an upper bound on positive excess; k=2 needs the ceiling; enriched need the floor).

**Not at stake:** CHANNEL_ID/CARRYLEMMA, R1–R30, R80–R82. Cheap (cached ρ + build_nu(11), 2.1s).
