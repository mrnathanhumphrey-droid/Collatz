# RESULT — §7 EFFECTIVE-CONSTANTS EXTRACTION (Hank): C_A = exp(O(A²)), effective for all n≥1, no n₀ barrier (2026-07-26)

**Source:** Hank's full read of Tao 1909.03562 §7 (+ Lemma 2.2) from `scratchpad/tao.txt`. Wilson's T2/T3/T4 (make
Tao's C_A explicit). Bookkeeping extraction, not a proof.

## HEADLINE — the gate answer
**|E e(−2πiξ·Syrac/3ⁿ)| ≤ C_A·n^{−A}, C_A = exp(O(A²)), EFFECTIVE for all n≥1.** Non-triviality window **n ≳ exp(cA)**.
- **No astronomical n₀(A).** Lemma 7.4 (the disjoint-triangle geometry) is **threshold-free** — holds every n≥1,
  triangle separation is the *absolute* constant `(1/10)log(1/ε)`, independent of A and n. The "n sufficiently large"
  worry lives nowhere in the geometry.
- The only "sufficiently large" is on the **induction scale m** (not n): `C_{A,ε}=exp(O(A/ε))=exp(O(A))` from **Case 2
  of Prop 7.8**. Below it the trivial base case `Q_m≤m^A` covers the range.
- **No Baker** (Remark 7.5: periodic structure "we will not exploit … beyond Lemma 7.4" — no linear-forms-in-logs, no
  ineffective Diophantine input). No entropy decrement, no compactness, no unlocated-scale pigeonhole. **Fully
  effective.**
- Combined with T1 (useful regime A≈2, grind A=2,3,4): `C_2~exp(O(4))`, `C_3~exp(O(9))`, `C_4~exp(O(16))` — benign,
  the tower/factorial worry is dead.

## ⚠️ OCR caveat (load-bearing for the algebra)
`pdftotext` garbles Tao's small absolute constant **ε** ("black iff |θ|≤ε") as a bare `3` or drops it. Every
`exp(−3·1_W)` / `exp(−3·#white)` in the dump is **`exp(−ε·…)`** (per-white-point gain), NOT `exp(−3)`. ε is declared a
*sufficiently small absolute constant*, chosen once, **independent of A and n** (implied constants do not depend on ε).
Genuine 3's also occur (Pascal value b=3; log 9; 10^A). Disambiguate per line.

## T4 — Prop 7.3 / 7.8 / the induction (the non-negative large-deviation core)
| Object | Statement (ε restored) | Depends on | Threshold |
|---|---|---|---|
| **Prop 7.3** (L3372) | `E exp(−ε·#{j∈[n/2]: b_j=3,(j,b_{[1,j]})∈W}) ≪_A n^{−A}` | A, ε | none (→7.37 via 7.8) |
| **Q recursion** (7.34) | `Q(j,l)=exp(−ε·1_W(j,l))·E Q((j,l)+Hold)` | ε, Hold law (abs) | — |
| **Base (7.39)** | `Q_m ≤ m^A` (Q≤1) — covers small m | none | all m |
| **Prop 7.8 monotone** (L4738) | `Q_m ≤ Q_{m−1}` for `C_{A,ε} ≤ m ≤ n/2` | A, ε | **m ≥ C_{A,ε}=exp(O(A/ε))** ← THE threshold |
| **Fwd induction** (L4746) | `Q_m ≤ C_{A,ε}^A ≪_{A,ε} 1` all m ⟹ (7.37) ⟹ Prop 7.3 | A, ε | — |

**Prop 7.8 three cases** (induction step `Q(j,l) ≤ m^{−A}Q_{m−1}`, `j=n/2−m`):
| Case | Condition | Threshold on m |
|---|---|---|
| **1** white start | `exp(−ε/2)m^{−A}Q_{m−1}` via (7.42) | `C₁~(A/ε)log(A/ε)` |
| **2** small triangle `s≤m/log²m` | **BINDING**: `(1+O(A/log m))(1−(3ε/4)≫1)≤1` needs `log m ≳ A/ε` | **`C₂~exp(cA/ε)`** |
| **3** large triangle `s>m/log²m` | hit O(A) whites, `P=O_{A,ε}(1)` steps, `R=O(A²)` triangles, union bound `≪A²·4^{−A}` | `C₃(A,ε)`, explicit |

## T3 — the geometry/renewal lemmas (the n₀ hunt) → NO bad n₀
| Lemma | Statement | Depends | Threshold |
|---|---|---|---|
| **7.2** cancellation (L3318) | white ⟹ `|f(3^{2j−2}2^{−l},3)| ≤ exp(−ε)` (=cos θ, |θ|>ε; even `exp(−cε²)`) | ε abs | **none** — per-white gain is an absolute <1 |
| **7.4** black structure (L3393) | `B` = disjoint triangles, each in `[n/2−(1/10)log(1/ε)]×Z`; any two separated ≥`(1/10)log(1/ε)` | **ε only** | **NONE, all n≥1** ← key finding |
| **7.5** (L~4335) | B periodic under `(0,2·3^{n−1})`; Baker could give more but **not used** | — | no Baker |
| **7.6** holding time (L4423) | `Hold` exp tail, mean (4,16), not in any proper-subgroup coset; `j~Geom(4)` | abs | none |
| **7.7** first-passage (L4540) | `P(v_{[1,κ]}=(j,l)) ≪ e^{−c(l−s)}(1+s)^{−1/2}G_{1+s}(c(j−s/4))` | c abs (Lem 2.2, d=2) | none |
| **7.9** many triangles→many whites (L5510) | `E[1_{R≤r}exp(−εΣ1_W+εR)] ≤ exp(ε)` | ε abs | none |
| **7.10** large triangles rare (L5730) | `P(E_{p,s'}) ≪ A²(1+p)/s' + exp(−cA²(1+p))` | **A**, c abs | `s'≥C·A²(1+p)` explicit |

**T3 verdict:** the "triangles not too large / n large" phrasing is discharged entirely inside Prop 7.8's scale
threshold `m ≥ C_{A,ε} ~ exp(cA/ε)`; **Lemma 7.4 is unconditional in n**. n₀ is not astronomical — it is `exp(O(A))` on
the induction scale. To certify `≤ n^{−A}` you need `n ≳ exp(cA)`. **The downstream (channel) budget closes.**

## T2 — assembly + the SINGLE lossy step (Wilson's T6 target)
Symbolic carry (replace every ≪ with these):
```
per-white gain        g        = exp(−ε),  ε absolute < 1/100
Case-2 scale thresh   C_{A,ε}  ≍ exp(c·A/ε)          [BINDING, L4818–5070]
fwd-induction cap     Q_m      ≤ (C_{A,ε})^A ≍ exp(O(A²))
Geom(4) A-th moment   E[j^A]   ≍ A! = exp(O(A log A))
FINAL                 C_A      ≍ exp(O(A²)),   |S(n)| ≤ C_A·n^{−A}
non-triviality        n        ≳ exp(c·A)
```
**The one genuinely LOSSY step** (not merely loose): **Prop 7.1 → 7.3**, `|g|≤1` + reduction to the *non-negative
white-point count* ("does not require capturing cancellation"). This is precisely where the heuristic `exp(−cm)` is
downgraded to `n^{−A}` — **a quantitative loss, not a barrier.** That is the localized lemma for the T6 sharpness audit:
if the forward bound comes in orders worse than the T6-fitted truth, this `|g|≤1` step is where the slack lives.

## Net
Tao's C_A is **effective, `exp(O(A²))`, no n₀ barrier, no Baker, no ineffective step** — the "tedious extraction" is
real and benign. Effective superpolynomial decay of sup|π̂| is now a *citable-and-computable* fact, not a folklore
≪_A. The residual open improvement (exp(−cm), Remark 1.15) is blocked only by the one lossy `|g|≤1` step, which T1
already shows we don't need (A≈2 suffices). **This closes the ℓ^∞ leg completely.** Not at stake: RECENTER, LAMBDA,
CHANNEL_ID, v₃ HIERARCHY, R1–R30. Next: T6 reverse-grind (fit §7 free constants to exact sup through k≤16) + the
aggregate/channel re-center.
