# Result 37 (qx+1 paper) — the level-3 gate is derived EXACTLY; a next-digit constant `σ₁` DOES enter (chain is a truncation past k=2), BUT its vanishing is not a new boundary (q=5,7 gap at σ₁=0).

**Date:** 2026-07-16. **Verdicts: ★ H_GATE3 ✓ CONFIRMED EXACT (0 mismatch / 1.2M constructed pairs / 6 primes) / ★ H_SIGMA_MATTERS ✓ CONFIRMED (no-σ form fails by exactly `y₁·j₁·σ₁` on every σ₁≠0 pair) / H_SIGMA_INDEP ✓ (σ₁=0 independent of s_R13, 7 primes <130) / H_L3RATE ≈ 1/q (measurement).**

**Headline: the 2nd-order expansion R20 left open is done — `W₂ ≡ Q₁ − y₁(P²−Q) mod q`, `P=j₁·s_R13`, `Q=j₁·σ₁+C(j₁,2)s_R13²`, with `σ₁` = 2nd q-adic digit of `(2^d−1)/q`. σ₁ genuinely enters the level-3 collision gate (coefficient `y₁·j₁`), so the clean k=2 chain is a TRUNCATION. But σ₁'s vanishing does NOT close the gap: q=5 and q=7 have σ₁=0 and healthy gaps (r₅≈0.62, r₇≈0.38, full-cascade). d=2 remains the sole boundary; σ₁ is a fine-structure correction to ω, not a new fixed point.**

Probe: `probe_37_level3_gate_sigma1.py`. Log: `result_37_level3_gate_log.txt`. Runtime: ~30 s.

## Why this probe fired now (the sequencing)

The re-entry recursion worksheet (`PHASE3_REENTRY_RECURSION_WORKSHEET.md`, §8 item 1) flagged that the recursion is exact only for the `k=1→2` step and *structurally extrapolated* beyond. R20 explicitly left level 3 as `W₂ + T₃ ≡ 0 mod q` with `W₂ := (U₁+T₂)/q` "a DEFINITION not a formula — needs a 2nd-order expansion of `2^{−jd} mod q³`." The arc's lesson (the `|ε_n|·2^n` plateau, `⌊r/2⌋+2`, both died exactly one level past their window): **measure the next level before trusting the extrapolation.** A bound proved on the k=2 form could be a bound on the wrong object.

## The derivation (done before running)

Write `2^d = 1 + q·s + q²·σ (mod q³)`, where `s := s_R13 = ((2^d−1)//q) % q` (the R13/R20 level-2 constant) and **`σ := σ₁ = (((2^d−1)//q)//q) % q`** (the NEXT q-adic digit — the object under test). Then `2^{j₁d} = 1 + q(j₁s) + q²(j₁σ + C(j₁,2)s²) mod q³`, so with `P := j₁s`, `Q := j₁σ + C(j₁,2)s²`:
```
    2^{−j₁d} = 1 − qP + q²(P² − Q)  mod q³
    U₁ = 2^{−S₁}[P − q(P² − Q)]  mod q²
    W₂ = (U₁ + T₂)/q  ⟹  W₂ ≡ Q₁ − y₁·(P² − Q)  (mod q)
```
`y₁ = 2^{−S₁} mod q`; `Q₁ = ([ (j₁s)·2^{−S₁} + 2^{−S₂} − 2^{−S′₂} ] mod q²) // q`. **The `σ₁` dependence sits inside `Q`: the coefficient of σ₁ in `W₂` is exactly `+y₁·j₁`.** Dropping σ₁ shifts `W₂` by `y₁·j₁·σ₁ mod q`.

## Ground truth & method

At k=3, `W₂ = (U₁+T₂)/q` is an **exact integer** once levels 1,2 pass — computed by big-int division from residues mod `q³`, no model, no tolerance. `W₂` depends only on `(S₁,S₂,S′₁,S′₂)` and `j₁` (NOT on `v₁` / the full address), so valid level-1&2-passing pairs are **constructed directly** (works at σ₁≠0 primes where full cell enumeration is infeasible — q=11 has 1.33M cells). 200k pairs per prime, tested as exact set-equality.

## H_GATE3 + H_SIGMA_MATTERS — the exact-iff test

| q | d | s_R13 | σ₁ | L2 pairs | with-σ bad | no-σ bad | pred no-σ | shift-ok | match |
|---|---|---|---|---|---|---|---|---|---|
| 5 | 4 | 3 | **0** | 200,000 | **0** | 0 | 0 | ✓ | ✓ |
| 7 | 3 | 1 | **0** | 200,000 | **0** | 0 | 0 | ✓ | ✓ |
| 11 | 10 | 5 | **8** | 200,000 | **0** | **200,000** | 200,000 | ✓ | ✓ |
| 13 | 12 | 3 | **11** | 200,000 | **0** | **200,000** | 200,000 | ✓ | ✓ |
| 23 | 11 | 20 | **3** | 200,000 | **0** | **200,000** | 200,000 | ✓ | ✓ |
| 41 | 20 | 32 | **8** | 200,000 | **0** | **200,000** | 200,000 | ✓ | ✓ |

- **H_GATE3 CONFIRMED:** the with-σ form matches exact `W₂` on all 1.2M pairs, zero mismatch. **The level-3 gate is now derived, closing R20's open item.**
- **H_SIGMA_MATTERS CONFIRMED:** at σ₁≠0 primes the no-σ form fails on **every** pair, and the failure equals the predicted shift `y₁·j₁·σ₁` (shift-ok ✓ = the with-minus-no difference is exactly `y₁j₁σ₁` on all pairs). **σ₁ is a genuine level-3 correction — the clean k=2 chain is a truncation.**
- Note q=5,7 (σ₁=0) show with-σ = no-σ trivially — they verify the *first-order* structure (`P²`, `s²` terms) but not σ₁ itself; the σ₁ content is exercised only at q≥11. (This is why the first run, testing only q=3,5,7, was silently vacuous on the real question — all three have σ₁=0. Caught and fixed.)

## H_SIGMA_INDEP — σ₁ vanishes independently of s_R13

`σ₁` is a separate q-adic digit; fixing `s_R13 ≠ 0` (the s=1 regime, all primes <1093) leaves `σ₁` free over `F_q`. Primes with **σ₁=0 & s_R13≠0**: q = 5, 7, 17, 31, 73, 89, 127 (and q=3). So σ₁ has its own vanishing locus, one digit deeper than s_R13's (`s_R13=0 ⟺ s≥2`, first at q=1093).

## ★ The load-bearing synthesis — σ₁ enters, but is NOT a new boundary

The obvious worry: does σ₁'s vanishing open a *second* gap-closing condition (a new place `r_q → 1`), the way `d=2` does? **No — and there is direct proof:**

> **q=5 and q=7 have σ₁=0 AND fully healthy gaps** (`r₅≈0.62`, `r₇≈0.38`). Those `r_q` were measured with the **full cascade** (R27 direct high-k, R32 operator), so they already include the σ₁=0 correction. **A prime can have σ₁=0 and still gap.** Therefore σ₁=0 is NOT gap-closing.

So the picture is:
- **σ₁ IS a genuine correction to the recursion's fine structure** — the phase/weight factor `ω` (worksheet §3/§7) gains a 2nd-order term. The L3 bound must be proved on the σ₁-inclusive object, not the k=2 truncation. *(This validates firing the probe: the extrapolation was NOT safe.)*
- **σ₁ does NOT add a fixed point.** `d=2` remains the sole gap-closing condition. σ₁ modifies `ω`'s value but not its boundary structure. *(This bounds the damage: no new degeneration axis to characterize.)*
- **σ₁ lives in the phase/weight, not the rate.** The level-3 pass rate stayed `≈1/q` (0.314, 0.198, 0.141 vs 1/q = 0.333, 0.200, 0.143) — σ₁ shifts *which residue* passes gate 3 (the target value `t`), not *how many* pass. Consistent with R36 (rate is boundary-blind; both boundaries live on the weight side).

## What this hands the L3 program

1. **Level 3 gate: DELIVERED** (`W₂` closed form, exact). R20's Phase-1c open item is closed; the cascade is now explicit two levels deep with the second-order carry.
2. **The bound targets the corrected `ω`** — not the k=2 form. The correction is the `+y₁ C(j₁,2)s² + y₁ j₁ σ₁ − y₁ P²` structure in `W₂`.
3. **The boundary characterization SURVIVES** — d=2 is still the only fixed point; σ₁=0 is benign (q=5,7 witnesses). So the worksheet's two-fixed-point structure (d=2 phase collapse, s_R13=0 index shift) is NOT joined by a third at σ₁=0.
4. **Open, one level deeper:** whether the level-4 gate introduces yet another digit (`σ₂`) with the same "enters-but-benign" character, and whether the correction terms telescope into a clean generating function for `ω`. Not fired (level 4 = `d⁴q³` structure; and the k=3 result may already be enough to write the bound).

## Not at stake
R10–R36. A refutation of H_GATE3 would kill only the 2nd-order derivation; the r_q gap measurements (R27/R32) and the boundary (d=2) are independent.

_Reporting discipline: derived on paper first, tested as exact set-equality vs big-int ground truth (no tolerance to mis-specify). The first run was silently vacuous (q=3,5,7 all have σ₁=0) — caught by reading the σ₁ column, fixed by adding σ₁≠0 primes via direct pair construction, and disclosed here rather than buried. H_SIGMA_MATTERS's "the chain survives (inert σ₁)" branch was pre-registered as the falsifier and committed to be reported as good news if it fired; it did not — σ₁ genuinely enters. The "not a new boundary" claim rests on the independent, full-cascade r₅/r₇ measurements, not on this probe. Author's structural priors this arc continue to land; the quantitative σ₁-value is reported, not fit._
