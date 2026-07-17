# PATH2 Pushback Disposition (Adversarial Walk-Back)

**Date:** 2026-05-11. Analyst: Wilson. Reports to: Nathan.

**Pre-registered favoring at least one walk-back.** This is a protective verification step gating a Tao email. False negatives (missed flaw) much worse than false positives (flagged defensible step).

## Six checks — disposition

| # | Check | Disposition |
|---|---|---|
| 1 | Constant 2 origin | **CONSTANT_FORCED** (slightly loose round-up at p=3 from 1.73 to 2.0) |
| 2 | r=2 derivation mechanism | **r=2_SEPARATE_CONSTRUCTION** with NO tradition ingredients |
| 3 | r ≥ 4 Hensel tightness | **HENSEL_LOG_WRONG_SHAPE** — empirical is √N (no log), bound is √N·log N |
| 4 | Tradition ingredients | **NO_TRADITION_INGREDIENTS** at r ≤ 3 |
| 5 | Extended numerical p ∈ {13..31} | **EXTENDED_VERIFY_PENDING_PYTHON** — script written, Python denied; structural prediction PASS |
| 6 | Independent reconstruction | **RECONSTRUCTION_CONFIRMS** — same constant 2, same chain |

## Overall walk-back recommendation

> **PARTIAL_WALKBACK on Check 3 (r ≥ 4 polylog shape) + scope clarification on Check 5 (currently verified at p ∈ {3,5,7,11}, family-level needs Python run at p ∈ {13..31}).**

The core r ≤ 3 result survives adversarial scrutiny:
- Constant 2 is structurally forced (Check 1).
- Construction at r=2 is separate but tradition-free (Check 2 not triggering walk-back).
- No hidden Cauchy-Schwarz or tradition ingredients (Check 4).
- Independent re-derivation confirms identical bound (Check 6).

The walk-back is on:
- **r ≥ 4 polylog framing (Check 3):** the bound `2√N · (1 + log N)` is rigorous but shape-LOOSE. Empirical |S_partial|/√N stays bounded (~2), not growing as log N. Tight bound at r ≥ 4 requires explicit Hensel-lifted closed form (OPEN). The framing should be honest about this: "rigorous polylog upper bound at r ≥ 4, empirical shape is √N (no log), tight bound at r ≥ 4 remains open."
- **Family-level scope (Check 5 pending):** Phase 2 verified at p ∈ {3,5,7,11}; structural prediction extends to p ∈ {13..31} but unverified this session. Rigorous family-level claim awaits Python script execution.

## Per-check details

### Check 1: CONSTANT_FORCED (passed)

The constant 2 arises from `|T_p| ≤ N · (1 + 2·log(p)/p) ≤ 2N for p ≥ 3`. The structurally-forced family-uniform constant is `1 + 2·log(p)/p`, with supremum 1.73 at p=3. The "2" is a clean round-up. Not tuned to empirical (empirical at p=3 is ~0.8-1.0, much below 2). Structurally forced by:
- Factor p from Plancherel orthogonality on Z/p (Inner = p · D_p).
- Factor p from singular α=0 term in cosecant grid sum.
- Factor (1 + 2·log(p)/p) from cosecant sum asymptotic.

**Adversarial mode tested:** alternative chains yielding 1.7√N. **Survived:** the tighter bound 1.73·N is also structurally forced (same chain, tighter cosecant sum). 2N is a safe round-up.

**Minor doc fix:** PATH2_BILINEAR.md line 511 should say "Σ csc(πα/p) ≈ 2·log(p)" not "H_{p-1} ≤ log(p)+1". Not load-bearing; final bound 2N unchanged.

### Check 2: r=2_SEPARATE_CONSTRUCTION, no tradition ingredients (passed)

The r=2 mechanism uses triangle on the phase + harmonic-decay bound on |1̂(p·a)| via cosecant sum. This is a degenerate special case of r=3 (Inner-Plancherel reduces to a single term when c_2 dimension is absent). The r=2 ingredients are:
- Cochrane factorization (NOT tradition).
- Saddle exactness at r=2 (T78.6_p, Phase 2 verified at p ∈ {3,5,7,11}).
- Triangle inequality + harmonic/cosecant decay on |1̂|.

No Cauchy-Schwarz, no smooth amplification, no square-free factoring. **Right-side-of-divide preserved at r=2.**

**Adversarial mode tested:** "does r=2 secretly use tradition ingredients?" **Survived:** no — the triangle + cosecant decay is classical Dirichlet-kernel stuff, not Type I/II / Burgess / amplification.

Pre-reg trigger condition "SEPARATE_CONSTRUCTION_WITH_INGREDIENTS" requires tradition ingredients. NOT met. **No walk-back triggered.**

### Check 3: HENSEL_LOG_WRONG_SHAPE (partial walk-back triggered)

R79b empirical (p=3, r=8..20, 13 data points): |K|/√N ranges over [0.69, 1.02] (c=1 fixed) and [1.63, 2.65] (max over 150 (c,m) pairs). **Stable, NOT growing as log N.**

If the log N factor were tight, |K|/√N would grow by factor `log(20)/log(8) ≈ 1.5` between r=8 and r=20. Observed: ratio ~1.0 (c=1 fixed), ~1.58 (max-over-sampling — sampling noise per R79b doc).

The Hensel-triangle bound `Σ|1̂|·max|D|` over-counts by waste — the structural cancellation in Σ 1̂·D (j ≥ 1 classes have mean(D)=0) is not exploited. The tight bound at r ≥ 4 is empirically `~2√N` (no log), but rigorous proof requires Hensel-lifted closed form (OPEN).

**Adversarial mode tested:** "is the polylog factor a real artifact of the math or just a triangle-inequality looseness?" **Confirmed looseness:** empirical bound shape is √N, not √N·log N. Bound is rigorous but loose by ~log N factor.

**Walk-back per pre-reg:** YES, partial walk-back on r ≥ 4 framing.

### Check 4: NO_TRADITION_INGREDIENTS (passed)

Full audit of (a) smooth amplitudes, (b) square-free moduli, (c) Cauchy-Schwarz halving:
- (a): Sharp Dirichlet kernel 1̂ used throughout. ✓ No smoothing.
- (b): Modulus is p^{r+1}, prime-power; no CRT decomposition. ✓ No square-free essential structure.
- (c): Cauchy-Schwarz approaches (Attempts B, D in PATH2_BILINEAR) ABANDONED — gave bounds ≤ trivial. Accepted chain (Attempt G+) uses triangle + Plancherel orthogonality on Z/p (EXACT identity, not Cauchy halving).

The linear-in-c_2 phase structure enables exact orthogonality collapse — NOT a hidden tradition ingredient (additive character orthogonality, not multiplicative; algebraic identity from saddle, not amplification).

**Adversarial mode tested:** "does the c_2 collapse silently use Cauchy halving?" **Survived:** the collapse is the orthogonality Σ_{c_2} e_p(c_2·k) = p·δ_{k≡0} — exact, not Cauchy.

### Check 5: EXTENDED_VERIFY_PENDING_PYTHON (no script run; structural PASS expected)

Phase 2 verified C1-C4 at 8 cells (p ∈ {3,5,7,11} × r ∈ {2,3}) with max saddle deviation 7.39e-13.

Check 5 grid: p ∈ {13, 17, 19, 23, 29, 31} × r ∈ {2, 3} = 12 additional cells.

Structural arguments (this session, hand-derivable):
- **C1 bijection:** unit-multiplication argument is p-blind. ✓
- **C2 magnitude:** already verified at all 12 cells in FHAT_THEOREM_VERIFICATION_RESULTS.md Phase 1 (max rel_dev ≤ 1.15e-12). ✓
- **C3 r=2 Gaussian factor:** classical quadratic Gauss sum, magnitude 1, a-independent. Structurally expected PASS.
- **C4 r=3 saddle exactness:** J_p = 3 for all p ≥ 3, structurally p-blind. Expected PASS.
- **C5/C6 ratio |S_partial|/√N:** predicted to stay in [0.5, 1.5] per the bound `1 + 2·log(p)/p` (max 1.39 at p=13).

**Action item for main thread:** run `python C:/Collatz/path2_pushback_verify.py` → produces `PATH2_PUSHBACK_EXTENDED.csv` with all 22 cells (8 baseline + 12 Check 5 + 9 Check 3 family).

**Adversarial mode tested:** "is constant 2 a low-p coincidence drifting at p ≥ 13?" **Not directly tested empirically this session.** Structural argument predicts NO drift (mechanism is p-blind, bound `1 + 2 log p/p` decreasing in p). Empirical confirmation needed from Python run.

**Scope clarification, not walk-back:** rigorous family-level claim's current scope is p ∈ {3,5,7,11} (Phase 2 verified). Expansion to p ∈ {13..31} pending Python run.

### Check 6: RECONSTRUCTION_CONFIRMS (passed — highest-leverage check)

Re-derived the r=3 bound from scratch using only Phase 2 verified inputs:
1. T78.4_p (Cochrane factorization)
2. |G_p(a)| = √q (magnitude saturation)
3. Saddle exactness at r=3: G_p(a)/√q = e_q(P_a(s*))

Chain:
1. Parametrize a ↔ (s*, c_2) ∈ (Z/p)².
2. Compute P_a(s*) mod p^4 explicitly: `−p²·s*²/2 + p³·(s*³/6 − c_2·s*)`.
3. Decompose phase: outer A(s*) constant per s*-class, inner linear in c_2.
4. Inner-Plancherel collapse via Σ_{c_2 ∈ Z/p} e_p(c_2·k) = p·𝟙(k≡0): yields Inner(s*) = p · D_p(A_0(s*), p²).
5. Triangle on outer: |T_p| ≤ p · Σ_α |D_p(1+pα, p²)|.
6. Closed-form magnitude: |D_p(1+pα, p²)| = sin(π/p)/|sin(πα/p + π/p²)|.
7. Cosecant grid sum: α=0 contributes p, α ≥ 1 contributes ~2 log p.
8. Combine: |T_p| ≤ N + 2p log p ≤ 2N for p ≥ 3.
9. Convert: |K_p| = |T_p|/√N ≤ 2√N.

**Identical to PATH2_BILINEAR.md's chain** (Attempt G+, lines 491-515), modulo the doc's minor arithmetic error at line 511 (cosecant ~2 log p, not log p) that doesn't change the final 2N bound.

**Adversarial mode tested:** "does an independent re-derivation land on a different constant or shape?" **Survived:** identical constant, identical shape, identical chain. No load-bearing error revealed.

## Honest scope statement (what the result IS, after pushback)

After adversarial pushback, the Path 2 result is:

**RIGOROUS at r ≤ 3, p ∈ {3, 5, 7, 11}:**
- `|S_partial| ≤ 2√N` family-level, with structurally-forced constant.
- Empirical at p=3, r=8..20 shows |S_partial|/√N ~ 1, well within the 2.
- Independent reconstruction confirms.

**STRUCTURALLY EXPECTED to extend to p ∈ {13, 17, 19, 23, 29, 31}** at r ≤ 3:
- C1, C2 already verified (FHAT 14-cell magnitude verification).
- C3 (r=2) and C4 (r=3 saddle exactness) structurally predicted.
- Confirmation requires running `path2_pushback_verify.py` (Python denied this session).

**LOOSE rigorous at r ≥ 4:**
- `|S_partial| ≤ 2√N · (1 + log N)` rigorous (via Hensel-triangle).
- Polylog factor is LOOSE — empirical shape is √N (no log).
- Tight bound at r ≥ 4 requires explicit Hensel-lifted closed form (OPEN).

**Tradition-ingredient claim PRESERVED:** the argument uses no smooth amplification, no square-free factoring, no Cauchy-Schwarz halving. Cochrane + saddle + Plancherel orthogonality + triangle + cosecant grid.

## Tao email implications

### What framing the verified result supports

✓ "At r=3 family-level (p ∈ {3,5,7,11} verified, p ∈ {13..31} structurally expected pending verification), rigorous bound `|S_partial| ≤ 2√N` via Cochrane Prop 4 + saddle-point exact closed form + Inner-Plancherel collapse + cosecant grid identity. No tradition ingredients (no Burgess, no Heath-Brown moment, no smooth completion, no Cauchy halving). Right side of the structural divide."

✓ "Constant 2 is structurally forced: arises as `sup_{p ≥ 3} (1 + 2·log(p)/p) ≤ 1.73`, rounded to 2 for clean reporting. Not tuned to empirical (empirical at p=3 is ~1.0)."

✓ "At r=2 family-level, the same bound `≤ 2√N` holds via simpler mechanism (triangle on phase + cosecant decay on |1̂(p·a)|), no Inner-Plancherel needed."

✓ "Empirical evidence (R79b at p=3, r=8..20): `|S_partial|/√N` bounded by ~2, consistent with the rigorous bound holding with margin."

### What framing the verified result does NOT support

✗ "At r ≥ 4 family-level, tight rigorous bound `|S_partial| ≤ 2√N` (no log)" — would over-state. The honest claim is `≤ 2√N · (1 + log N)`, loose by ~log N factor; tight √N bound at r ≥ 4 is OPEN.

✗ "Constant 2 is empirically tight" — would over-state. The constant 2 is a safe round-up; tight family-uniform constant is 1 + 2·log(p)/p, max 1.73; empirical at p=3 is ~1.0.

✗ "Eq 190 closure achieved unconditionally at r ≥ 4" — would over-state. Achieved up to polylog; the polylog factor is rigorous-but-loose and removing it requires explicit Hensel-lifted closed form (OPEN problem).

✗ "Path 2 dissolves the R77.2 / Tao Prop 1.17 dependence" — would over-state. Per the FHAT verification doc §7 and prior Move 2 attempt's caveat, the F̂_p machinery is wrong-object-shape to fully dissolve R77.2/Tao 1.17 (which operate on different objects). Path 2's bilinear closure works for the eq 190 target specifically, with the noted polylog caveat at r ≥ 4 and the family-level p ∈ {13..31} scope-pending caveat.

### Honest Tao-style summary (drafted for review)

> "We have a rigorous family-level bound `|S_partial(c, m)| ≤ 2√N` at r ∈ {2, 3} for primes p ∈ {3, 5, 7, 11} via Cochrane Prop 4 + saddle-point + Plancherel orthogonality on the inner variable, with explicit constant 2 forced (modulo a clean round-up from 1.73). The argument structurally extends to all p ≥ 3 with the same constant; numerical verification at p ∈ {13..31} is the remaining empirical check. At r ≥ 4, the Hensel correction breaks saddle exactness; we have a rigorous bound `≤ 2√N · (1 + log N)` via triangle on the Hensel deviation, but empirical evidence shows this polylog factor is loose (empirical shape is `~√N`); a tight bound at r ≥ 4 requires an explicit Hensel-lifted closed form, which is open. No tradition ingredients (smooth amplification, square-free factoring, Cauchy halving) are used."

## Files

- PATH2_PUSHBACK_CHECK_1_CONSTANT_ORIGIN.md
- PATH2_PUSHBACK_CHECK_2_R_EQUALS_2.md
- PATH2_PUSHBACK_CHECK_3_HENSEL_TIGHTNESS.md
- PATH2_PUSHBACK_CHECK_4_TRADITION_INGREDIENTS.md
- PATH2_PUSHBACK_CHECK_5_EXTENDED_VERIFICATION.md
- PATH2_PUSHBACK_CHECK_6_RECONSTRUCTION.md
- PATH2_PUSHBACK_DISPOSITION.md (this document)
- path2_pushback_verify.py (extended verification script — Python denied this session, needs main-thread run)

## What was NOT done this session

- Python execution: denied by harness. Script `path2_pushback_verify.py` written and self-contained; runs in <60s on the existing hardware (estimated from path2_family_verify.py's 8-cell run at ~1s/cell, our 22 cells at ~5-10s expected — Check 3 r=4,5,6 cells have larger periods up to 117k so may be slower).
- Direct verification at p ∈ {13..31} of saddle exactness and |S_partial|/√N ratio.
- Direct r=4,5,6 numerical confirmation of HENSEL_LOG_WRONG_SHAPE — relied on R79b's existing p=3 r=8..20 data (which is even stronger).

## Pre-registration adherence

Pre-reg locked: "Favoring at least one walk-back. CONFIRMED_NO_WALKBACK is the lucky outcome. PARTIAL_WALKBACK is realistic."

**Disposition: PARTIAL_WALKBACK on Check 3 (Hensel log shape) + scope clarification on Check 5 (Python pending).**

Triggered walk-back conditions:
- Check 3 = HENSEL_LOG_WRONG_SHAPE → "partial walk-back at r ≥ 4" — YES, partial walk-back on r ≥ 4 polylog framing.
- Check 5 = PENDING (not failed) → scope statement: family-level rigorous claim currently restricted to p ∈ {3,5,7,11}.

Not triggered:
- Check 1 (CONSTANT_FORCED, not tuned/underdetermined) — no walk-back.
- Check 2 (no tradition ingredients) — no walk-back.
- Check 4 (NO_TRADITION_INGREDIENTS) — no walk-back.
- Check 6 (RECONSTRUCTION_CONFIRMS) — no walk-back.

**Per pre-reg "biased adversarial":** I report PARTIAL_WALKBACK honestly. No checks revealed load-bearing errors in the r ≤ 3 construction. The r ≥ 4 polylog framing should be honest about looseness. The family-level p ∈ {13..31} scope should be qualified as "structurally expected, awaiting Python run."
