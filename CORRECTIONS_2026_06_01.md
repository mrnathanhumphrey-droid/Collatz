# Corrections / Walk-backs (2026-06-01)

**Critical bug found in Fraction-based p_m distribution machinery.** Numerical offset_distribution method confirmed correct via c(m) matching at 30 digits across m = 0..4.

## The bug

`probe_pm_distribution_2026_05_31.py` (depth-2 and any depth-k+ extension) had a structural error in Case A:

When κ_1 = 0 (level-1 dominant, δ_1 = 0 globally) at depth 2:
- d_1 = δ_2^(0) must equal 0 → level-2 must ALSO be dominant
- Hence σ_2 = δ_3^(0), which is the depth-0 distribution of (A_3, B_3)
- Case A contribution should be (1/3) · p_0 of iid copy

The script used `(1/3) · num_p1` (the depth-1 distribution) instead of `(1/3) · num_p0` (the depth-0 distribution). This contaminates the within-coset structure of p_m for m ≥ 2 while leaving QR/NQR mass split nearly intact (why c(2) was off by only 1.4e-4 despite p_2(σ) being off by 0.13 at individual σ).

Same structural error would propagate to p_3 and beyond in `probe_p3_exact_2026_05_31.py`.

## What's WALKED BACK (do not trust)

1. **`pm_distributions_2026_05_31.json`** — DISCARDED, renamed to `DISCARDED_BUGGY_pm_distributions_2026_05_31.json`. The p_1 and p_2 rational distributions inside it are WRONG. Do not load.

2. **"Mass concentration on σ = ±4 ≈ 38%"** — based on buggy p_2. The TRUE numerical p_m looks different (peaks shift across depths in the actual distribution). σ = ±4 may or may not be the peak; need to recompute.

3. **"Character order|4 restriction"** — INCORRECT. The σ↔−σ symmetry restricts μ_∞ to 8 EVEN characters (orders 1, 2, 4, 4, 8, 8, 8, 8), not 4 characters of orders | 4. The order-8 characters DO contribute.

4. **"QR mass = 0.577, NQR mass = 0.423, 15% nonperturbative asymmetry"** — c_∞ = 0.153 is correct (matches to 50 digits), and that equals QR mass − NQR mass, so the SPLIT is right. But the within-coset interpretation tied to ±4 concentration was bug-driven.

5. **"|⟨χ_4, μ_∞⟩|² = c_∞ conjecture"** — FALSE. Actually ⟨ω, μ_∞⟩ = 0. Numerical |⟨ω, p_m⟩|² decays geometrically at ratio (3/5)² = 0.36, exactly matching λ_dom(ω) = −3/5 squared. The conjecture was based on the buggy Fraction p_1, p_2.

6. **"Decay ratio 30× per step" for the conjecture** — artifact of bug.

## What's still CORRECT

1. **Mod-4 dichotomy theorem**: c_∞(q) ≡ 0 ⟺ q ≡ 3 mod 4. Pure algebraic, no buggy data dependency.

2. **c(0) = 19/127 exact**. c(1) = 33-digit num/den exact. c(2), c(3), c(4), c(5) at high numerical precision via offset_distribution method (not the Fraction depth-2 machinery).

3. **c_∞ = 0.15298912060588517527891674877413229926086222622334** (50 digits, from offset_distribution).

4. **σ ↔ −σ symmetry of μ_∞**: forced by D = X−Y → −D = Y−X swap. Pure structural, not bug-dependent.

5. **Per-character dominant eigenvalues over Q(i)**: λ_dom(χ) = 3/(4χ(2) − 1).
   - λ_dom(χ_L) = 1
   - λ_dom(ω, order 4) = −3/5
   - λ_dom(order-8 chars) = (−3 ± 12i)/17, magnitude 3/√17
   These derivations were pure algebra. Still good.

6. **50-digit PSLQ negatives against 720 Q(i) Hecke L-values + 4 frontier frameworks**: These are facts about c_∞ as a number. The negatives stand even though our directional motivation (1+2i hypothesis) was based on shaky inputs.

## Trustworthy method going forward

Use `probe_p3_numerical_2026_06_01.py` style: `offset_distribution_mp(q, n, A_MAX)` at appropriate dps, extract p_m by summing P_X auto-correlation at `d = j · q^m` for j = 1..q−1, normalize.

- p_0 fits in N = q² = 289 entries, instant
- p_1 fits in N = q³ = 4,913 entries, < 1 sec
- p_2 fits in N = q⁴ = 83,521 entries, ~15 sec
- p_3 fits in N = q⁵ = 1.4M entries, ~12 min (computed today)
- p_4 fits in N = q⁶ = 24M entries, ~11 min (computed today as part of sanity)
- p_5 fits in N = q⁷ = 410M entries, RAM-bound (see c(6) note)

## Numerical results that are TRUSTED

| m | c(m) | \|⟨ω, p_m⟩\|² | diff c_∞ |
|---|---|---|---|
| 0 | 0.149606 | 0.159744 | +6.8e−3 |
| 1 | 0.153178 | 0.055145 | −9.8e−2 |
| 2 | 0.153248 | 0.019825 | −1.33e−1 |
| 3 | 0.153005 | 0.007126 | −1.46e−1 |
| 4 | 0.152989 | 0.002564 | −1.50e−1 |

Decay ratios for |⟨ω, p_m⟩|²: 0.345, 0.360, 0.360, 0.366 → asymptotic 9/25 = (3/5)² ✓

## The bigger lesson

The session's downstream conclusions about |⟨ω⟩|² = c_∞ as a "structural identity" came from numerical agreements that were code-bug artifacts. The structural derivations themselves (mod-4 theorem, dominant eigenvalues, σ↔−σ symmetry) survive.

**Previous nulls don't matter much.** The Q(i) Hecke L scan, Fantini-Rella, LSSW, Kriz-Nordentoft falsifications are FACTS about c_∞ as a number at 50 digits, but the FRAMING (directional hypotheses based on the (1+2i) eigenvalue interpretation) was on shaky inputs from the start. Future search should not lean on those directional priors.

We're back at: c_∞ is 0.15298912... to 50 digits, surrounded by structural facts but with no closed-form candidate currently surviving. Fresh slate for hypotheses.
