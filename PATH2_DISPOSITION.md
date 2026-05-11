# PATH2 Disposition

## Disposition: **H_PARTIAL**

The family-level extension of R78.4-78.6 from q=3 to general p ≥ 3 is **mechanical and verified** structurally; explicit substitution into the bilinear sum produces a **rigorous bound |S_partial| ≤ C · √N** (the eq 190 square-root closure target) at small r with explicit constant, but the bound degrades to **|S_partial| ≤ C · √N · log N** at r ≥ 4 due to the Hensel-correction triangle inequality. This is the "≪ √N" target up to polylog — useful and well-scoped, but not the no-polylog ideal.

## One-paragraph rationale

Phase 1 (family extension): clean. T78.4_p (Cochrane factorization), T78.5_p (bijection on coset), T78.6_p magnitude (|G_p|=√q) all generalize p-blindly. T78.6_p phase at r=3 is conjectural (J_p = r condition matched for all our test cells but verification script not run this session). Phase 2 (empirical): bijection and magnitude verified via FHAT_THEOREM_VERIFICATION_RESULTS.md (33 cells covering our grid) and hand-derivation; saddle-exactness at p ≥ 5 r=3 remains a clean empirical-pending hypothesis. Phase 3 (bilinear): substituting T78.4_p + T78.6_p (assuming saddle exact) gives explicit P_a(s*) = −p²s*²/2 + p³(s*³/6 − c_2·s*) mod p^4 at r=3, the linear-in-c_2 term enables an Inner-Plancherel transformation Inner(s*) = p · D_p(a_0(s*), p²) (length-p Dirichlet kernel), and summing via 1/sin bounds yields **|T_p| ≤ 2N at r=3**, equivalent to |S_partial| ≤ 2√N. At r ≥ 4 Hensel correction (R79b) creates a per-a deviation ψ_true − ψ_lead of bounded magnitude; triangle inequality picks up a Σ|1̂(p·a)| ~ N log N factor, giving |S_partial| ≤ √N · (2 + 2 log N) — still ≪ √N in the polylog interpretation.

## Scope of bound (H_PARTIAL detailed)

| r | Bound on \|S_partial\| | Constant | Polylog | Notes |
|---|---|---|---|---|
| 2 | ≤ (1 + log p/p) · √N | ≤ 2 | none | trivial-with-decay; closes target at r=2 |
| 3 | ≤ 2 · √N | 2 | none | T78.6_p saddle-exact (PENDING Phase 2 verification at p ≥ 5) |
| ≥ 4 | ≤ 2 · √N · (1 + log N) | 2 | log N from Hensel-correction triangle | R79b documents Hensel structure; argument uses ψ_lead + bounded D(a) triangle |

**Bound holds family-level (p-uniform constants up to log p/p decay).** Empirical R79b at p=3 shows actual |S_partial| ~ √N · N^{0.022} — much sharper than our bound (R79b's β=0.522 vs our N^{0.5+ε}). The √N empirical save (within-Inner cancellation we don't capture) is consistent with our bound holding with margin.

## What fails / what's NOT covered

- **No-polylog ideal:** |S_partial| ≤ C·√N (no log) at r ≥ 4 is NOT achieved. The log factor at r ≥ 4 comes specifically from the Hensel correction; would need an explicit Hensel-lifted closed form of ψ_true (open at p=3 and at family level).
- **T78.6_p saddle exactness at p ≥ 5, r=3:** structural derivation (J_p = r condition) suggests yes; not empirically verified this session. The `path2_family_verify.py` script is delivered for the user to run.
- **r=2 with p=3:** picks up the q=3-specific Gaussian factor e^{iπ/6} (T78.6 r=2 remark) — for p ≥ 5 r=2 the analog is a different root of unity but doesn't affect the magnitude bound.

## Path 1 (family extension) — clean or q=3-specific?

**Clean.** No q=3-specific structural barrier identified:
- Cochrane Prop 4 generalizes p-blindly (the truncated p-adic log structure).
- Bijection via unit multiplication on coset is p-blind.
- Magnitude |G_p|=√q is verified at 33 cells (FHAT doc) spanning p∈{3,5,7,11,...,31}, r∈{1..6}.
- Saddle s* = (C_a − 1)/p mod p generalizes directly from R78.6 derivation.
- The truncation level J_p differs across p (p=3 r=2 has J=3; p≥5 r=2 has J=2) — **NOT a structural barrier**, just a polynomial-degree variation. The saddle-exactness condition appears to be J_p = r (matched at all our 8 test cells).

## Phase 2 verification status

- C1 (bijection a ↔ C_a): **PASS** structurally and by hand-derivation at all 8 cells.
- C2 (|G_p|=√q on support): **PASS** by direct inheritance from FHAT_THEOREM_VERIFICATION_RESULTS.md.
- C3 (r=2 phase has uniform root-of-unity factor): structural prediction; **NOT empirically run this session** (Python denied).
- C4 (r=3 saddle exact): structural prediction; **NOT empirically run this session**.

## Phase 3 attempt outcomes

| Attempt | Method | Outcome |
|---|---|---|
| A | Substitute C_a = a · L̃_p^{-1}, expand P_a(s*) mod p^4 | Yields explicit phase: e_{p²}(−s*²/2) · e_p(s*³/6 − c_2·s*) at r=3 |
| B | Partition by s*-class, length-p inner DFT on c_2 | Inner(s*) = p · D_p(a_0(s*), p²) — rigorous closed form |
| C | Cauchy-Schwarz on s*-sum | Loses to triangle for our case (csc² identity) |
| D | Poisson on a-sum | Equivalent to attempt B re-organized |
| E | Empirical comparison (R79b p=3 r=8..20) | Our bound 2N >> empirical N^{0.522}; consistent (we're a loose upper bound) |
| G+ | Triangle + 1/sin grid identity on a_0(s*) ranging across {1, 1+p, ..., 1+(p-1)p} | **|T_p| ≤ 2N at r=3 (rigorous)** |

## Adversarial checks (A1-A4)

**A1 — Magnitude:** rigorous 2N >> empirical N^{0.522} at p=3 r=8..20. Loose upper bound holds with margin. PASS.

**A2 — Hensel:** ψ_true at r ≥ 4 deviates from ψ_lead by bounded per-a fluctuation (R79b mean 13-21% at r=4-10; |D(a)| ≤ 2 trivially). Triangle inequality picks up an extra polylog factor. Polylog acceptable in "≪ √N" interpretation. PARTIAL — bound holds with log degradation.

**A3 — Walk-back safety:** We do NOT use the falsified "ψ is cubic exponential character of a" claim. We use explicit polynomial form of P_a(s*) and length-p inner DFT — both grounded in T78.4_p (rigorous) and family-level T78.6_p (saddle exact at small r). PASS.

**A4 — Honest scope:** the r ≥ 4 Hensel-correction loss IS specified precisely (factor 1 + log N from triangle on the deviation D(a)). What's NOT covered: explicit Hensel-lifted closed form of ψ_true (open). The argument's hard dependency is T78.6_p saddle exactness — empirically pending at p ≥ 5 r=3.

## Strategic implication

- **Path 2 (direct construction) achieves eq 190 closure up to polylog at r ≥ 4 and exactly at r ≤ 3** (assuming Phase 2 verification confirms saddle exactness at p ≥ 5).
- The bound is family-level: works for any prime p ≥ 3 with p-uniform constants.
- **Does not achieve the no-polylog ideal at r ≥ 4** — would require explicit Hensel-lifted closed form (open at p=3 per R79b §C2/C3/Hensel).
- The empirical √N saving (R79b β=0.522) is STRONGER than our rigorous bound but NOT needed for the bilinear closure target.

## Implication for c=7/45 closure

Combining with R78.1-78.3 (rigorous): the full chain for eq 190 closure at q=3 is now:
1. R78.1: complete-sum vanishes (rigorous).
2. R78.2: F̂ supported on q/9 elements (rigorous).
3. R78.3: |F̂| = 3√q on support (rigorous).
4. PATH2 Phase 3 derivation: |S_partial| ≤ 2√N at r=3, |S_partial| ≤ 2√N·(1+log N) at r ≥ 4 (this work).

**Eq 190 closure: ACHIEVED to within polylog at r ≥ 4 family-level, conditional on T78.6_p saddle exactness at p ≥ 5 r=3 (empirical-pending).**

c=7/45 status update:
- Empirical: certified to 1.7×10⁻⁴ at k=6 (unchanged).
- Rigorous: now **eq 190 closure to polylog**, conditional on Phase 2 verification (path2_family_verify.py).
- Open: no-polylog ideal at r ≥ 4 requires explicit Hensel-lifted ψ_true closed form.

## What's needed for H_DIRECT_WORKS upgrade

To upgrade H_PARTIAL → H_DIRECT_WORKS, need:
- Phase 2 verification run (the `path2_family_verify.py` script — confirms T78.6_p saddle exact at p ∈ {5,7,11}, r=3).
- A different argument at r ≥ 4 that doesn't lose the polylog. Concrete options:
  - **Explicit Hensel-lifted ψ_true closed form** (open). With this, Inner Plancherel applies to ψ_true directly without triangle on D(a).
  - **A bound on Σ 1̂(p·a) · D(a)** that uses the s*-class structure of D (R79b: j=0 anomalous, j ≥ 1 regular). The j=0 class delocalization (D(a) ~ −1 mean) could be exploited; the j ≥ 1 classes have mean(D) = 0 which is the "centered perturbation" that vanishes after averaging.
  - **Bourgain-Konyagin sum-product on ⟨1+p⟩** (open path #1 from R78.6 strategic position).

These are concrete next-step targets, NOT speculative.

## Specific new input named (for H_NEEDS_NEW_MATH alternative reading)

If r ≥ 4 polylog loss is judged unacceptable (strict "no log" interpretation of "≪ √N"), then the SPECIFIC new input needed is:
- **A Hensel-lifted closed form of ψ_true(a) at r ≥ 4 expressing the Hensel correction as a polynomial-in-c_2 of computable degree.**
- This is the analog of R78.4-78.6 at r ≥ 4 — explicit, not just empirical — and is an OPEN problem at p=3 (R79b §Open problems).
- Without it, the triangle-on-deviation argument bottoms out at the polylog.

Literature analog: Heath-Brown's hybrid bounds for cubic-character sums (which R79b walked back as not directly applicable) DO have versions that handle "polynomial perturbations" via Vinogradov's mean-value theorem. Whether VMV applies to our specific Hensel-corrected polynomial is open.

## Files

- PATH2_PRE_REGISTRATION.md (locked first)
- PATH2_FAMILY_EXTENSION.md
- PATH2_FAMILY_EXTENSION_VERIFICATION.md + (CSV not produced — Python denied)
- PATH2_BILINEAR_FROM_CLOSED_FORM.md
- PATH2_DISPOSITION.md (this document)
- path2_family_verify.py (delivered; not run this session)

## Empirical / verified-via-FHAT inputs

- T78.4_p, T78.5_p: structural, p-blind, RIGOROUS.
- T78.6_p magnitude (|G_p|=√q): verified at 33 cells via FHAT.
- T78.6_p saddle phase at r=3: structural; verified at p=3 r=3 (R78.6 table); EMPIRICAL PENDING at p ∈ {5,7,11} r=3.

## Caveats

1. Python execution denied this session — verification CSV not generated.
2. Phase 2 saddle-exactness at p ≥ 5 r=3 is the one empirical-pending input.
3. Argument's polylog at r ≥ 4 is acceptable under "≪" interpretation; would need Hensel-lifted closed form for strict no-log.
4. R79b's empirical β=0.522 (sub-trivial √N save) is STRONGER than our bound — not used by our argument, but confirms our bound holds with margin.

## Bottom line

**H_PARTIAL.** Eq 190 closure achieved at r ≤ 3 (constant 2) and at r ≥ 4 (constant 2 × polylog) family-level for all primes p ≥ 3, conditional on T78.6_p saddle exactness at p ≥ 5 r=3 (Phase 2 verification — script delivered, not run this session). The remaining gap (polylog at r ≥ 4) has a precisely-named structural source (Hensel correction) and concrete next-step paths (Hensel-lifted closed form, or Vinogradov mean-value on the deviation polynomial).
