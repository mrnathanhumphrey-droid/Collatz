# Probe SAT — the station at infinity — **the "station" is the telescoping boundary datum (not independent tail info); but the tail amplitude disfavors period-9 and favors a LONG period (~17), refining the ≈22 pre-reg down**

**Date:** 2026-07-23. Probe `probes/probe_SAT_station.py`. Uses the global constraint `Σ_{r≥1} Λ_r = −1/10` as a
long-baseline datum and **tests** candidate periods against it (period is an input, not an output — immune to the
under-determination that killed Probe D). Exact `Λ_r` r≤7; float r=8..11 (ε to k=12).

## SAT-A — the station (exact, and an honest correction)
`2Λ_r = ε_{r+1}−ε_r` **telescopes**, so `T_M := Σ_{r≥M} Λ_r = (ε_∞ − ε_M)/2 = −ε_M/2` (given `ε_∞=0`, i.e. S→7/15):

| M | T_M = −ε_M/2 | \|T_M\|/\|Λ_{M−1}\| | sign | ε_M exact? |
|---|---|---|---|---|
| 8 | +3.728e−4 | 1.74 | + | yes |
| 9 | +3.76e−6 | 0.01 | + | float |
| 11 | −7.510e−4 | 1.92 | − | float |
| **12** | **−1.1374e−3** | **2.94** | − | float |

**The station `T = Σ_{r≥12} Λ_r = −ε_12/2 = −1.1374×10⁻³`** (Wilson's rough −1.4e−3, refined). **Honest correction:**
`T_M = −ε_M/2` is the *telescoping identity* — the "station" is the **boundary datum ε_M**, not independent
information beyond k=12; the `−1/10` global constraint is tautologically `−ε_1/2` (`ε_1 = 1/5`) plus convergence, and
carries no info about tail *structure*. **What IS informative:** `|T|/|term| ≈ 2.94` — a tail this large relative to a
single term, from a slowly-damped (`ρ≈0.984`) oscillation, requires a **long period** (small θ ⟹ strong
`1/|1−ρe^{iθ}|` amplification).

## SAT-B — the period scan (fit a,b to r=8..11, ρ=0.984 fixed; predict T; compare to exact)
**No period reproduces the tail with ρ=0.984** — every P *under*-predicts:

| P | 9 | 12 | 18 | 22 | 26 |
|---|---|---|---|---|---|
| T_pred | −6.6e−4 | −6.8e−4 | −6.0e−4 | −5.2e−4 | −4.2e−4 |
| rel err vs T | **42%** | 41% | 47% | 54% | 63% |

Best P=11 at 40% — i.e. a single ρ=0.984 mode fit to r=8..11 **cannot make a tail as large as −1.14e−3** (it gives
~−6e−4, about half). So the tail is *bigger* than the local window predicts — the station sees slow/long structure
the short baseline misses (the VLBI point), and **period-9 fits worst-to-middling, never well.**

## SAT-C — robustness (windows; ρ free)
| window | consistent-P (rel < 20%) |
|---|---|
| r=7..11 | [8.0, 13.5] |
| r=9..11 | [6.5, 20.5] |

The consistent-P set is **broad and window-dependent** (~[8, 20]), not pinned — the 4-point window can't resolve it.

## SAT-D — nearly model-free half-cycle from the tail
Solving `T ≈ −(2/π)·h·|Λ₁₃|/(1+ρ^h)` with `|Λ₁₃| ≈ ρ|Λ₁₁| = 3.80e−4`, `T = −1.137e−3`: **h = 8.8 ⟹ period P = 2h ≈
17.6.** This **refines Wilson's pre-registration**: his rough tail (−1.4e−3) gave P≈22; the exact tail (−1.137e−3, 19%
smaller) gives **P≈17.6.** Either way: **not 9** — the tail is too large for a period-9 oscillation to produce from
the observed term sizes.

## SAT-E — the caveat (load-bearing)
This test **assumes `Σ Λ_r = −1/10` (S→7/15)** — `ε_12` is measured *relative to* 7/15, so the tail `−ε_12/2` is
conditional on 7/15. It tests 7/15 and the period **jointly**: the output is "*if* 7/15, *then* P ∈ [range]." The
range here (~15–20, model-free ~17.6) is plausible and non-empty, so it is **consistent with** 7/15 + a single long
mode — it does not prove either.

## Status
**Probe SAT: the station is the telescoping boundary datum (not new tail info), but the tail amplitude disfavors
period-9 and favors a LONG period ~17.** **A** — `T = −ε_12/2 = −1.137e−3` exact via telescoping; the `−1/10`
constraint is `−ε_1/2` + convergence, carrying no structural tail info (honest correction of the "unreachable-region
datum" framing); the informative content is the *amplitude* `|T|/|term| ≈ 2.94`. **B** — no ρ=0.984 mode fit to
r=8..11 reproduces so large a tail (all under-predict; period-9 40%+ short). **C** — consistent-P is broad/window-
dependent (~[8,20]), unpinned. **D** — model-free `P ≈ 17.6`, refining Wilson's ≈22 down (exact tail is 19% smaller
than his estimate). **E** — theorem-conditional; the plausible non-empty range is consistent with 7/15 + one long mode.

**Consequence for the crux (owed to the pen).** SAT converges with the campaign's other period results into one
picture: **rate ≈0.984 (R81/F1), period NON-integer (F1-B excludes all integers 2–9), and LONG (~15–22; SAT
model-free ~17.6).** Period-9 is a crossover artifact (as R26's ">16" and the retired float-era long-period fits
already hinted). But like Probe D, SAT **cannot pin** the period — the tail is one integral constraint sensitive to
long/slow structure, and it says "long," not a value; the 4-point local window stays under-resolved. The clean
measurement still needs the pure region to span ≥1 period (ε to k≈17–20). No fitting of the period (each was tested,
not fit); the telescoping-identity correction and the window-dependence reported plainly; the ≈22→≈17.6 refinement
follows from using the exact tail. Not at stake: R1–R30, R80–R82, Probe D.
