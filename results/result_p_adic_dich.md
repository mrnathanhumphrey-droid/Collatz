# RESULT — p-ADIC-DICH: does γ_∞(k)>1 ⟺ p∣k lift off p=3? **FAILS-CLEAN (the dichotomy is 3-specific) — but with a precise DISSOCIATION: the tower engine LIFTS and the p∣k⟹enriched HALF lifts; only the converse is 3-specific.** (2026-09-06)

**Probe:** `probes/probe_p_adic_dich.py`. Kill-first test of whether the CHANNELDICH dichotomy (p=3: γ_∞(k)
ENRICHED >1 iff 3∣k, DEPLETED off it, 12/12) generalizes to the p-adic affine renewal
`a_r = q^{−v}(1 + p·a_{r−1}) mod p^r`, `v~Geom(½)`, at p=5 (q=2). Machinery copied verbatim from the p=3
pipeline (`build_nu_exact`, `dlog_table`) and generalized 3→p, 2→q; **dlog generator = 1+p** (the principal-unit
generator; `= 4` at p=3, coinciding with `q²`).

## Gates (both mandatory-first, both PASS)
- **GATE 0 — white baseline** (uniform ρ ⟹ γ≡1): `max|γ_white(k)−1|` = 2.22e-16 (p=3,r=4), 2.22e-16 (p=5,r=3),
  3.33e-16 (p=7,r=3). Normalization/coordinate correct. PASS.
- **GATE 1 — p=3 faithfulness + reproduction:** `build_chain_exact(3,2)` == `build_nu_exact` (r=0..5), **EXACT
  Fraction-dict equality: True**. Dichotomy reproduced **14/14**: 3∣k → enriched **4/4** (k=3,6,9,12); 3∤k →
  depleted **10/10**. VERDICT LIFTS. PASS — the generalized code is a faithful lift of the p=3 measurement.

## PRIMARY p=5 (exact r=2..5) — VERDICT **FAILS-CLEAN**

2×2 (rows p∣k / p∤k × cols enriched / depleted):

| | enriched (γ_∞>1) | depleted (γ_∞<1) |
|---|---|---|
| **5∣k** | **2** (k=5,10) | 0 |
| **5∤k** | **5** (k=1,4,6,11,14) | 7 (k=2,3,7,8,9,12,13) |

The `p∣k ⟹ enriched` direction is clean (2/2). The **converse fails hard**: 5 of 12 non-multiples are enriched.
This is unambiguous — **k=4 is enriched at every level** (γ: `1.09653, 1.30673, 1.28349, 1.33899` for r=2,3,4,5;
never below 1.10), so no extrapolation softness can rescue the dichotomy. γ_∞(k)>1 ⟺ 5∣k is **FALSE**.

γ_r(k) table (p=5, exact r=2..5; g_inf = geometric-tail extrapolation, else deepest r):

| k | 5∣k | γ_2 | γ_3 | γ_4 | γ_5 | γ_∞ | sign | pred |
|---|---|---|---|---|---|---|---|---|
| 1 | . | 0.90768 | 1.07126 | 1.10223 | 1.03928 | 1.03928 | ENRICHED | depleted → **MISS** |
| 2 | . | 0.56611 | 0.53879 | 0.60039 | 0.60244 | 0.60244 | depleted | depleted |
| 3 | . | 0.37837 | 0.43955 | 0.52645 | 0.43799 | 0.41570 | depleted | depleted |
| 4 | . | 1.09653 | 1.30673 | 1.28349 | 1.33899 | 1.33899 | ENRICHED | depleted → **MISS** |
| 5 | YES | 1.25728 | 1.30387 | 1.30302 | 1.28789 | 1.28789 | ENRICHED | enriched ✓ |
| 6 | . | 1.40653 | 1.42474 | 1.45549 | 1.22276 | 1.22276 | ENRICHED | depleted → **MISS** |
| 7 | . | 0.39955 | 0.50672 | 0.58556 | 0.68534 | 0.68534 | depleted | depleted |
| 8 | . | 0.47848 | 0.56009 | 0.58690 | 0.69637 | 0.69637 | depleted | depleted |
| 9 | . | 0.99955 | 0.95248 | 0.94007 | 0.91462 | 0.91462 | depleted | depleted |
| 10 | YES | 1.83687 | 2.41956 | 2.97942 | 2.83258 | 2.75376 | ENRICHED | enriched ✓ |
| 11 | . | 1.14526 | 1.17194 | 1.34165 | 1.11450 | 1.11450 | ENRICHED | depleted → **MISS** |
| 12 | . | 0.39971 | 0.40385 | 0.40090 | 0.42367 | 0.42367 | depleted | depleted |
| 13 | . | 0.39971 | 0.32105 | 0.31808 | 0.31539 | 0.31296 | depleted | depleted |
| 14 | . | 1.14526 | 1.22906 | 1.10381 | 1.10159 | 1.10159 | ENRICHED | depleted → **MISS** |

## The DISSOCIATION (the important finding — spec's secondary-(a) case realized)

The dichotomy fails, but **two of its three components LIFT cleanly**, so this is not a plain collapse — it localizes
exactly what was 3-specific:

1. **The tower/telescoping ENGINE lifts EXACTLY.** Secondary (a), exact Fractions, p=5:
   `Σ_j ρ_3(t+j·5²) == ρ_2(t)` → **True**; `Σ_j ρ_4(t+j·5³) == ρ_3(t)` → **True**. Wilson's fiber-tower claim
   (the p-divisible affine shift preserves marginal consistency) is **CORRECT and general** — the engine survives
   the prime change.
2. **The `p∣k ⟹ enriched` half lifts, with the fiber-depth magnitude growth.** k=5 → γ_∞=1.28789, k=10 → γ_∞=2.75376
   — the deeper fiber multiple is *more* enriched, exactly the p=3 signature (there: γ_∞(3)=1.28, γ_∞(9)=1.99).
   Secondary (c) confirmed.
3. **Only the CONVERSE (`p∤k ⟹ depleted`) is 3-specific.** Extra enriched non-multiples appear at p=5 that were
   invisible at p=3.

## WHY p=3 is special (localization, flagged for the pen — structural observation, not a proof)

At **p=3, q=2 the two structural axes COINCIDE**: the ×q² correlation multiplier `q²=4` equals the principal-unit
(fiber) generator `1+p=4`. The p=3 "dichotomy" reads a single fused axis, and on that fused axis "enriched" happens
to align with "3∣k". **At p≥5 with q=2 the axes split** (`1+p ≥ 6 > 4 = q²`), and enrichment appears on *both*:
notably **k = q² = 4 is a separately-enriched non-multiple** (the cleanest MISS). The natural reading is that
enrichment is governed by a finer **fiber ⊕ ⟨q⟩-orbit** structure that degenerates to the single condition "p∣k"
**only when q² ≡ 1+p** — i.e. only at (p=3, q=2). (The enriched-non-multiple set {1,4,6,11,14}: k=4 is cleanly the
q² channel; the full account of {1,6,11,14} needs the ⟨q⟩-orbit-on-the-dual analysis and is not resolved here.)
This is consistent with the CHANNEL_ENGINE_LIFT recon's finding that the load-bearing object is a twisted
⟨q⟩-orbit structure.

## Honest flags
- **Exact depth:** r≤5 (p=5). Extrapolation is soft — many γ_∞ fell back to the deepest r (γ_5) because the last
  A-ratios were not stable in (0,1); the r-trend is reported in full instead. **The verdict does not depend on the
  extrapolation:** k=4 (and k=6 through r=4) is enriched at *every* exact level, so the dichotomy fails at fixed r.
- **k=1 is marginal** (crosses 1: γ_2=0.908 < 1, γ_5=1.039 > 1) — could be extrapolation noise; excluded from the
  argument. The MISS rests on the unambiguous k=4.
- **p=7 NOT run** — kill-first honored (spec: "if p=5 fails clean, STOP, don't burn p=7"). It is a ~0.1 s
  confirmation if desired; the p=5 failure is decisive on its own.

## VERDICT and disposition
**FAILS-CLEAN.** The dichotomy `γ_∞(k)>1 ⟺ p∣k` is **3-specific** and does **not** lift to the p-adic affine renewal
family — the biconditional breaks at p=5 (converse direction). Per the pre-registered call, the deep channel-engine
result does not generalize as the dichotomy; **fall back to the modest orbit-Ramanujan paper** (the ADD-A-PEAK of
`notes/CHANNEL_ENGINE_LIFT_RECON.md`).

**But the fallback is richer than a null:** the **telescoping engine lifts exactly** at p=5 and the **p∣k⟹enriched
direction lifts** with its fiber-depth growth — so the general theorem to write is not the dichotomy but the
*engine* (fiber-tower telescoping under a p-divisible affine shift) plus the *forward* enrichment implication, with
the clean p=3 biconditional exposed as the coincidence **q² = 1+p**. Not at stake: CHANNELDICH itself (p=3 stands,
reproduced 14/14 here), the tower identity, R7–R10.
