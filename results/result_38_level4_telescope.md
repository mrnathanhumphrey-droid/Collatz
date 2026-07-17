# Result 38 (qx+1 paper) — the digit corrections TELESCOPE: σ₂ enters level 4 with the SAME coefficient σ₁ had, so the whole tower is the base-q expansion of ONE constant `τ = (2^d−1)/q`. ω's digit-dependence has a closed generating function.

**Date:** 2026-07-16. **Verdicts: ★ H_GATE4_SANITY ✓ (order-3 model = exact W₃, 0 mismatch) / ★ H_SIGMA2_ENTERS ✓ (σ₂-dropped model fails on all 120k pairs at σ₂≠0 primes) / ★★ H_TELESCOPE ✓ CONFIRMED (failure = `y₁·j₁·σ₂`, the same coefficient σ₁ carried at level 3) / H_SIGMA2_INDEP ✓ (σ₂ free over F_q).**

**Headline: the level-4 gate needs a 3rd digit σ₂, and it enters `W₃` with coefficient `+y₁·j₁` — IDENTICAL to how σ₁ entered `W₂` (R37). The mechanism: `2^{jd} = 1 + q·j·τ + O(q²τ²)`, `τ := (2^d−1)/q`, so the q-adic digits of the linear term `qjτ` are `j·(digits of τ)` — exactly the measured per-digit coefficient. ⇒ the infinite tower of digit corrections R37 opened is the base-q expansion of a SINGLE constant `τ`, entering ω LINEARLY. No infinite regress; ω's dependence on the tower is closed-form in `τ`. This is the object the L3 bound wants.**

Probe: `probe_38_level4_telescope.py`. Log: `result_38_level4_telescope_log.txt`. Runtime: ~40 s.

## Why this fired (R37's follow-up, worksheet §8 item 3)

R37 derived the level-3 gate exactly and found σ₁ (2nd q-adic digit of `(2^d−1)/q`) genuinely enters — "the clean k=2 chain is a truncation." The open worry: does each deeper level bring a NEW independent digit-constant (σ₂, σ₃, …), forcing the L3 bound to track an infinite tower? Or do the corrections telescope into something closed? This probe decides it at the next level.

## Method — no hand-derivation of `W₃` (the R37 substitution trick, one level up)

The digits `(s, σ₁, σ₂, …)` enter the cascade ONLY through `2^{−j₁d}` (`j₁·d` is the one guaranteed multiple of `d`). So substitute a **truncated model** for that single quantity and read its effect on the exact big-int `W₃`:
- `pow2_true = 2^{−j₁d} mod q⁴` (exact; = inverse of `(1+qs+q²σ₁+q³σ₂)^{j₁}`).
- `pow2_model` = inverse of `(1+qs+q²σ₁)^{j₁} mod q⁴` — **σ₂ dropped** (order-2).

`pow2_model = pow2_true mod q³` (differ only at the `q³` digit, by exactly `j₁σ₂`), so the model still passes levels 1,2,3 cleanly (`W₃` well-defined) and `W₃_model` differs from `W₃_true` only by the σ₂ image. Ground truth `W₃ = (W₂+T₃)/q` is exact integer division at `k=4`. Pairs constructed to pass levels 1,2,3 directly. 120k pairs/prime.

## Results — exact-iff vs big-int ground truth

| q | d | s | σ₁ | σ₂ | pairs | order-3 bad | order-2 bad | pred order-2 | tele bad | telescope |
|---|---|---|---|---|---|---|---|---|---|---|
| 11 | 10 | 5 | 8 | **0** | 120,000 | 0 | 0 | 0 | 0 | ✓ (control) |
| 13 | 12 | 3 | 11 | **1** | 120,000 | **0** | **120,000** | 120,000 | **0** | ✓ |
| 23 | 11 | 20 | 3 | **0** | 120,000 | 0 | 0 | 0 | 0 | ✓ (control) |
| 41 | 20 | 32 | 8 | **15** | 120,000 | **0** | **120,000** | 120,000 | **0** | ✓ |

- **H_GATE4_SANITY:** the full (order-3) model reproduces exact `W₃`, zero mismatch — the level-4 plumbing is correct.
- **H_SIGMA2_ENTERS:** at σ₂≠0 primes (q=13, 41) the σ₂-dropped model fails on **every** pair. σ₂ genuinely enters the level-4 gate. Controls (σ₂=0: q=11, 23) match perfectly, confirming the failure is σ₂ and nothing else.
- **★★ H_TELESCOPE CONFIRMED:** `W₃_true − W₃_model ≡ y₁·j₁·σ₂ (mod q)` on all pairs (tele bad = 0). **σ₂ enters `W₃` with coefficient `+y₁·j₁` — the exact coefficient σ₁ carried into `W₂` at level 3 (R37).**

## ★★ The generating function — why the tower is one constant

`2^{jd} = (1 + (2^d−1))^{j} = 1 + j(2^d−1) + C(j,2)(2^d−1)² + …`, and `(2^d−1) = q·τ` with `τ := (2^d−1)/q`. So:
```
    2^{j₁d} = 1 + q·j₁·τ + q²·C(j₁,2)·τ² + …   (q-adic expansion in τ)
```
The **linear term `q·j₁·τ`** has q-adic digits `j₁·(digits of τ) = j₁·(s, σ₁, σ₂, …)`. That is EXACTLY the per-digit coefficient the probes measure: σ₁ entered `W₂` as `y₁·j₁·σ₁` (R37), σ₂ enters `W₃` as `y₁·j₁·σ₂` (R38) — the `y₁` from `2^{−S₁}`, the `j₁` from the linear q-adic term. So the telescoping is not an inductive guess from two points; it is the q-adic image of a single binomial identity, and the general digit follows.

> **ω's dependence on the LTE tower is a closed generating function: the linear part is `y₁·j₁·τ`, `τ = (2^d−1)/q`, entering through the SINGLE constant τ (its digits are not independent knobs — they are the base-q representation of τ). The nonlinear parts are powers `τ², τ³, …` with coefficients `C(j₁,2), C(j₁,3), …` (the `P²`, `s²` structure of R37). There is NO infinite tower of independent corrections.**

## What this hands the L3 program

- **R37's "truncation" worry is resolved.** The k=2 chain is a truncation of a closed generating function in `τ`, not the head of an unbounded regress. The L3 bound can carry the exact tower constant `τ = (2^d−1)/q` symbolically and get all orders at once.
- **Write ω with `τ`.** Inside the character sum (worksheet §7), the digit-dependence collapses to `y₁·j₁·τ` (linear) + `C(j₁,2)τ²` (…) — a clean object, one constant per prime.
- **The boundary is still `d=2`, and `τ` enters benignly.** `τ` is a q-adic unit iff `s = τ mod q ≠ 0`, i.e. the `s=1` regime (all primes `<1093`). `τ ≡ 0 mod q` is `s≥2` = R35's benign onset-shift. And `τ = s` alone (all higher digits 0) at q=5,7 — which GAP (r₅≈0.62, r₇≈0.38, full-cascade). So the tower constant modifies ω's value, never its boundary. `d=2` remains the sole fixed point; the two-fixed-point structure stands.

## Honest scope
Confirmed at the first two correction digits (σ₁ via R37, σ₂ via R38), with the q-adic binomial giving the general term. Levels ≥5 (σ₃…) not separately run — the generating-function argument makes them the same computation, and the binomial is exact. The nonlinear coefficients `C(j₁,n)τⁿ` were verified only insofar as the order-3 model (which contains them) matches exact `W₃`; their individual telescoping was not isolated (not needed — they are already closed powers of τ).

## Not at stake
R10–R37. A refutation of H_GATE4_SANITY would kill only the level-4 plumbing; the r_q gaps (R27/R32), the d=2 boundary, and R37's level-3 gate are independent.

_Reporting discipline: same substitution method as R37 (no error-prone hand derivation of W₃; big-int ground truth). σ₂=0 controls (q=11,23) included specifically so a spurious "everything fails" bug would be caught — they correctly show 0 mismatch. H_TELESCOPE's falsifier (a DIFFERENT coefficient ⇒ no clean telescoping) was pre-registered and committed to be reported; it did not fire — the coefficient is `y₁·j₁`, matching R37. The generating-function claim rests on the binomial identity (exact), with the probe confirming the first two digits behave as it predicts._
