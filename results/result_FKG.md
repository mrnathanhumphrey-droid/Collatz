# Probe FKG — is the positivity a correlation inequality? — **NO. FIFTH DEATH: the correlation-inequality shelf is CLOSED. ν is NOT log-supermodular (MTP2) in either coordinate — 26–31% of elementary squares violate, worst-ratio 48→233 GROWING with r (group) and 16→178 (dlog); not remotely approximate FKG. So FKG/Holley/Griffiths — the one shelf whose theorems carry a SIGN in the conclusion — does not apply. Two refinements: (C) the ω-triple shows the target is NOT "the ratio distribution prefers its own multiplier 4" — C(4ω) actually EXCEEDS C(4); 4 is the MIDDLE value and beats only the AVERAGE of its two twists, a local CONVEXITY of the ratio distribution at the fiber lag (margin +0.00196→+0.00066, shrinking), not a preference; and C is NOT elevated on the transport's multipliers 𝒮={4^a(−2)^b} (indistinguishable from generic lags), so the "collisions favour the map's own multipliers" mechanism dies too. FKG-B: ×4 scrambles the group order (~52%, chance) and is only ~90% monotone in dlog (s→s+1, carries break 10%) — but with FKG-A failed the bridge is moot.**

**Date:** 2026-07-25. Probe `probes/probe_fkg.py`, log `logs/fkg_run.log`. Exact rationals, r=4..8, banked ν. Reformulation: lag k in dlog = mult by 4^k, so `C(k)=Pr[R=4^k]` (R=X'/X, X,X' iid ~ν); fiber shift N/3 = mult by ω=4^{3^{r−1}} (order-3 gen). Target = `Pr[R=4] > ½(Pr[R=4ω]+Pr[R=4ω²])`. **FKG-A passing is not a result — without FKG-B's ×4-monotone bridge it says nothing about d₁ (stated first, per guardrail).**

## FKG-A — log-supermodularity FAILS in both coordinates, all r
`ν(x∨y)ν(x∧y) ≥ ν(x)ν(y)` on Z/3^r ≅ {0,1,2}^r, elementary squares (MTP2 ⟺ local condition on a product of chains):
| r | group: viol% / worst | dlog: viol% / worst |
|---|---|---|
| 4 | 26.4% / 48.6 | 21.3% / 15.8 |
| 6 | 29.8% / 116 | 28.3% / 96.0 |
| 8 | 30.9% / **233** | 30.5% / **178** |

**Both coordinates fail decisively and the worst-ratio GROWS with r** — this is not "approximate FKG with small violations," it is robustly anti-supermodular at fine scales. Pre-registered threshold was <1% violations / worst <1.05; measured is 30% / hundreds. **Shelf closed in both coordinates.** (Digit order: componentwise ∨/∧ is coordinate-symmetric ⟹ MSB-first and LSB-first give identical elementary-square sets and identical stats — verified.)

## FKG-B — the bridge (moot, but recorded)
Fraction of covering pairs `x ⋖ y` with `4x ≤ 4y` componentwise:
- **group** (a→1+4a): 51.6–53.2% — chance level, **×4 scrambles the order, no bridge.**
- **dlog** (s→s+1): 81.5% → **90.6%** (rising with r) — `s→s+1` preserves digit order except at base-3 carries (~10%), so ×4 is *approximately* monotone in dlog. But **FKG-A fails in dlog anyway**, so the bridge is moot.

## FKG-C — the ω-triple, and the mechanism check (both refine/kill)
`C(k)/C(0)` at the target lags, every r:
| r | C(4) | C(4ω) | C(4ω²) | Pr[R=4] − ½(twists) |
|---|---|---|---|---|
| 4 | 0.2304 | **0.2557** | 0.2012 | +0.00196 |
| 6 | 0.1786 | **0.1894** | 0.1651 | +0.00135 |
| 8 | 0.1458 | **0.1507** | 0.1396 | +0.00066 |

**The target inequality holds at every r (margin >0, shrinking) — but the reformulation's intuition is WRONG: `C(4ω) > C(4)`.** The exact multiplier 4 is the **middle** value, not the peak; it beats only the **average** of its two ω-twists (which straddle it). So the mechanism is a **local CONVEXITY of the ratio distribution at the fiber lag** (4 above the midpoint of its straddling neighbors), NOT "the ratio distribution prefers its own multiplier." This is exactly the AC-LAGS curvature `Δ²C(N/3)`, now read on the ratio distribution.

**𝒮-enhancement DIES:** C on the transport's own multipliers `𝒮={4^a(−2)^b}` has median 0.20–0.29, **indistinguishable from the generic-lag median** (0.16–0.26) at every r. So "collisions favour the map's own multipliers" is **not** the mechanism — C is not elevated on 𝒮 (pre-registered: this kills that explanation).

## FKG-D — not run (FKG-A failed; guardrail: no Gibbs fit without supermodularity)

## Verdict — the fifth death
**KILL condition fired: both coordinates fail FKG-A at every r.** The correlation-inequality shelf — FKG, Holley, Griffiths, the only shelf whose theorems have a SIGN in the conclusion — is **closed**: ν is not MTP2, and getting less so with r. The five deaths now: **R29** (no finite transfer operator), **PRODFORM** (no lacunary factorization), **decay-shelf** (Fourier/Rajchman machinery sign-blind), **CONE** (no forward-invariant cone; positivity is convergence to ν, not self-generated), **FKG** (no correlation inequality; ν not ferromagnetic). Every general route to *proving* the sign is exhausted. What remains is the bare, measured fact: `Pr[R=4] > ½(Pr[R=4ω]+Pr[R=4ω²])` — a local convexity of the Syracuse ratio distribution at the fiber lag, all quantities nonnegative probabilities, margin 0.19%→0.066% shrinking, **with no known theorem-shelf to deliver it and no mechanism (not preference, not 𝒮-enhancement, not ferromagnetic order) surviving.** The pen faces a bespoke convexity statement about the specific ratio distribution of ν.

## Status
**FKG (shelf test, 5th death):** ⚠️**FKG-A FAILS both coords all r** (group 26→31% viol, worst 48→233; dlog 21→30%, worst 16→178; GROWING) ⟹ ν NOT log-supermodular/MTP2 ⟹ **correlation-inequality shelf CLOSED** (FKG/Holley/Griffiths carry a sign but don't apply). **FKG-B:** ×4 scrambles group order (~52%=chance), ~90% monotone in dlog (s→s+1 carries), moot given A fails. ⭐**FKG-C: target `Pr[R=4]>½(Pr[R=4ω]+Pr[R=4ω²])` HOLDS all r (+0.00196→+0.00066, shrinking) BUT C(4ω)>C(4) — 4 is the MIDDLE, beats only the AVG of straddling twists = local CONVEXITY of the ratio distribution at the fiber lag, NOT "prefers its own multiplier."** ⭐**𝒮-enhancement DIES: C on {4^a(−2)^b} indistinguishable from generic ⟹ "collisions favour map's multipliers" mechanism dead.** FKG-D skipped (A failed). ⭐**5th DEATH (R29+PRODFORM+decay+CONE+FKG): every sign-bearing shelf closed. Remaining = bare convexity of the Syracuse ratio distribution at the fiber lag, nonneg probabilities, margin →0.066%, no shelf, no surviving mechanism.** Not at stake: R1–R30, R80–R82, all Thread-3 probes. commit pending.
