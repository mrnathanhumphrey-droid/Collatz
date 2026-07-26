# RESULT — PROBE HIERARCHY: the dichotomy is a v₃(k) hierarchy (2026-07-26)

**Probe:** `probe_hierarchy.py`. Wilson's upgrade: the enriched/depleted dichotomy is not "3∣k" but a **hierarchy in
`v₃(k)`** (the 3-adic valuation). Mechanism, one line: `G_j=⟨4⟩` has order `3^j`, so `4^k=id in G_j ⟺ 3^j∣k ⟺
v₃(k)≥j`. For every level `j≤v₃(k)` the channel-k collision **is literally the m=0 event**, giving an exact base case
at every valuation.

## The exact base case — VERIFIED (Fractions), including Wilson's decisive γ_j(9)=X_j check

`γ_j(k) = X_j` (the m=0 partial sum) for all `j ≤ v₃(k)`, then departs and relaxes. `X_1=5/3`, `X_2=15/7`,
`X_3=2.60437…` (exact big rational), `X_4=3.0686`.

| k | v₃ | γ_1 | γ_2 | γ_3 | γ_4 |
|---|----|-----|-----|-----|-----|
| 3 | 1 | **5/3 = X_1** ✓ | 10/7 ≠ X_2 (departs) | — | — |
| 9 | 2 | **5/3 = X_1** ✓ | **15/7 = X_2** ✓ | 1.912 ≠ X_3 (departs) | — |
| 27 | 3 | **5/3 = X_1** ✓ | **15/7 = X_2** ✓ | **2.6044 = X_3** ✓ | 2.372 ≠ X_4 (departs) |

**Wilson's decisive check `γ_j(9)=X_j` for j=1,2 passes exactly** (`γ_2(9)=15/7` to the Fraction). The hierarchy is
right — and it explains a "coincidence" from CHANNELDICH: `γ_∞(9)≈2.112` sits just below `X_2=15/7=2.143` because k=9
tracks the m=0 channel for 2 levels then barely relaxes.

## The k=27 prediction — HOLDS
v₃(27)=3, so `γ_1..3(27) = X_1,X_2,X_3` exactly (2.6044 to the Fraction), then relaxes: **`γ_16(27)=2.36139`, inside
the predicted [2.3, 2.9]**, enriched as v₃=3 requires.

## The hierarchy grouped by v₃ — γ_∞ increases with v₃, tracks X_{v₃}

| v₃ | base X_{v₃} | side | channels (γ_16) |
|----|-------------|------|------------------|
| 0 | γ_1=0.667 | depleted <1 | 1:0.730, 2:0.473, 4:0.861, 5:0.765, 7:0.428, 8:0.750, 10:0.521, 11:0.592 |
| 1 | X_1=1.667 | ENRICHED >1 | 3:1.237, 6:1.372, 12:1.528, 15:1.451, 21:1.424, 24:1.375 |
| 2 | X_2=2.143 | ENRICHED >1 | 9:2.112, 18:1.962 |
| 3 | X_3=2.604 | ENRICHED >1 | 27:2.361 |

## Why this is the better route — the proof factors into a proved fact + two specific bounds

Budget usage `|log Π_relax| / |threshold|` per channel (relaxation from base X_{v₃}, or γ_1=2/3 for v₃=0):

| binding | k=4 | k=3 | k=6 | k=5 | k=1 | others |
|---------|-----|-----|-----|-----|-----|--------|
| budget% | **63.0** | **58.4** | 38.1 | 33.9 | 22.4 | ≤29, several 0 (safe direction) |

Exactly reproduces Wilson's numbers (63/58/23%). **The dichotomy gets *easier* as v₃ grows** (X_j grows, the floor
`−log X_{v₃}` deepens, more room), so the binding comparison is only **v₃=1 vs v₃=0** — finite and specific, not
uniform-in-k. The proof factors into:

1. **`X_v > 1` for v≥1** — Cauchy–Schwarz / the m=0 linear divergence (the six-sightings result), **already proved**.
2. **Bounded relaxation on the two tightest channels, k=3 and k=4** (58%, 63%) — a specific finite target.

The binding channels' relaxation rates are fast (k=3: 0.73, k=4: 0.66), so the r>16 tail is tiny and the thresholds
hold with margin at the observed values; the remaining work is a *rigorous* bound on that relaxation, now
channel-specific rather than uniform.

## Net
- **The dichotomy is a v₃(k) hierarchy with an exact base case at every valuation** (`γ_j(k)=X_j` for j≤v₃(k)),
  verified exactly for k=3,9,27; Wilson's γ_j(9)=X_j check passes to the Fraction (15/7).
- **k=27 prediction holds** (2.361 ∈ [2.3,2.9]). Twelve+ confirming channels, grouped cleanly by v₃.
- **The proof reduces to one already-proved fact (X_v>1, m=0 divergence) + bounded relaxation on k=3,4** — not a
  uniform-in-k requirement, and not a sign. This is the arc's best-shaped result: exact base case, one line of group
  theory for the mechanism, and a testable prediction that landed.
- Hank's target is now sharp and channel-specific: bound `Σ_j|q_j(k)−1/3|` for k=3 and k=4 (relaxation/upper-bound
  type, the decay shelves are the right type).

**Not at stake:** CHANNEL_ID/CARRYLEMMA, R1–R30, R80–R82. Cheap (cached ρ + build_nu(11), 3.6s).
