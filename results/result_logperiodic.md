# Probe LOG-PERIODIC — the bore is LOG-periodic in log₃(r), not periodic in r — **which explains every disagreeing period estimate; and 7/15 is now conditional on an unobserved turnover near r≈27 (walk-back #45)**

**Date:** 2026-07-23. Exact ε chain from `ε_1=1/5, ε_{r+1}=ε_r+2Λ_r` with exact Λ (r≤7), float (r=8..11), and
**ν-Λ (r=12..16, validated 0.00% vs ε-telescoping)**. Prompted by Wilson's walk-back #45 (my "decay ⟹ turnover"
argument in `result_crossing.md` was wrong — decay is not crossing).

## The ε sign chain and its crossings (exact/validated)
`ε_r`, r=1..17: `+ +  − − − − − − −  + + + + + + + +` — **positive r=1,2; negative r=3–9; positive r=10–17+.**
Interpolated zero-crossings:

| crossing | location r | log₃(r) | nearest 3^k |
|---|---|---|---|
| #1 (+→−) | **2.65** | 0.888 | ~3¹ (transient→bore) |
| #2 (−→+) | **9.01** | **2.001** | **3² = 9 (dead-on)** |
| #3 (+→−) predicted | **≈ 27–28** | ~3.0 | **3³ = 27** |

`log₃(9.01) = 2.001` — the second crossing sits **essentially exactly on 3².** Half-periods: **6.36** (crossing
#1→#2), then **≥7 and counting** (crossing #2→now, un-crossed at r=16). The ×3 rule (half-period 6.36 → 19) predicts
crossing #3 at **r ≈ 28**, matching the independent `|Λ_r|` extrapolation (r≈22–28).

## The reframe: log-periodic, not periodic — the period was never fixed
**Half-periods are lengthening** (≈2, 6.36, ≥7→~19) — geometric, ratio ≈3. **A fixed period cannot produce
lengthening half-cycles.** This single fact **explains the entire campaign's period confusion:** every estimate
disagreed (9, 9.06, 9.5, 17.6, 22, 26, 37) and the MOON joint fit railed identically across all channels **because
they were all fitting a fixed frequency to a signal that doesn't have one.** The failures of Probe D, SAT, and MOON
were **misspecification, not resolution.** The correct model is **log-periodicity in log₃(r)** — oscillation
periodic in the *3-adic level's logarithm*, with zero-crossings near powers of 3.

**This is exactly the Lapidus p-adic prediction (F1-E).** A p-adic self-similar string is lattice (Lapidus–Hùng),
and its complex dimensions are periodic in `log(scale)`; the scale here is `3⁻ʳ`, so the geometric oscillation is
periodic in `log₃(r)`, with features near `r = 3^k`. The `9.01 ≈ 3²` crossing is the signature. F1-E located the
oscillation's origin in the p-adic level structure; this locates its *form*.

## 7/15 is now CONDITIONAL on an unobserved turnover (walk-back #45)
**Corrected:** decay ≠ crossing. No-turnover extrapolation (Λ>0 decaying 0.91 from r=17):
`Σ_{r≥17}Λ ≈ +2.37×10⁻³ ⟹ ε_∞ ≈ +9.9×10⁻³ ⟹ S_∞ ≈ 0.4766`, **not 7/15 = 0.4667.** For 7/15 one needs
`Σ_{r≥17}Λ = −2.58×10⁻³` — a substantial negative excursion **not observed** in the data (Λ is positive through
r=16, ε still rising). **So 7/15 hinges on the crossing #3 turnover near r≈27.** The turnover is *well-motivated* —
ε has crossed twice already, the log₃ structure predicts crossing #3 at 3³, and Lapidus theory requires the
oscillation — but it is **unobserved and past the compute wall** (build_nu walls at r≈16–17 on `3^{r+1}` dense
arrays). This is the first time a straightforward extrapolation of the deepest data does not land on 7/15.

## ⚠️ Danger-class caveat (#30): two crossings are two points
The log₃/powers-of-3 structure rests on **two** crossings (#1 partly transient-contaminated; only #2 at 9.01 is a
clean post-transient bore crossing). Two points fit a geometric law trivially — **the literal #30 abundance-trap
shape.** The prediction (crossing #3 at r≈27, half-period ~19) is falsifiable but **requires a third crossing to
confirm**, and r≈27 is ~10 levels past the build_nu wall. So log-periodicity is a **strong, theoretically-motivated
hypothesis, not an established fact.** No fixed-period model is viable (lengthening half-cycles rule them out); a
log-periodic model is consistent with everything and predicted by Lapidus, but under-determined by two crossings.

## Status
**LOG-PERIODIC: the bore is log-periodic in log₃(r) (crossings geometric near powers of 3), which explains why every
fixed-period estimate disagreed — misspecification, not resolution; and 7/15 is conditional on an unobserved turnover
near r≈27.** Exact/validated: crossings at r≈2.65 and **r≈9.01 = 3² (log₃=2.001)**; half-periods lengthening
6.36→≥7 (×3 rule ⟹ crossing #3 at r≈28); ν-Λ validated 0.00% vs ε. **Reframe:** no fixed period exists (lengthening
half-cycles), so Probe D/SAT/MOON failed by fitting a fixed frequency to a log-periodic signal — this is the unified,
correct explanation, superseding "window < ½ period" with "wrong model class." **Correction (walk-back #45):** the
earlier "amplitude decays ⟹ turns over" was wrong (decay ≠ crossing); no-turnover extrapolation gives S_∞≈0.4766,
so **7/15 is conditional** on the r≈27 turnover — motivated by two prior crossings + log₃ structure + Lapidus, but
**unobserved.** **Caveat:** two crossings = #30 two-points; the log₃ law needs a third crossing (r≈27, past the
compute wall) to confirm. Theoretical anchor: Lapidus p-adic complex dimensions periodic in log(scale)=log₃(r)
(F1-E). No fitting beyond linear crossing-interpolation; the log₃ law flagged as hypothesis; 7/15-conditional stated
plainly as the first non-7/15 extrapolation in the record. Not at stake: R1–R30, R80–R82 (the exact structural
identities and M-reality); at stake and now CONDITIONAL: the asymptotic value 7/15 (pending the turnover).
