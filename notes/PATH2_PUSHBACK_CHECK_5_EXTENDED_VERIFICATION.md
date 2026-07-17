# PATH2 Pushback — Check 5: Extended numerical verification at p ∈ {13, 17, 19, 23, 29, 31}

**Adversarial frame:** Phase 2 verified `path2_family_verify.py` C1-C4 at 8 cells (p ∈ {3,5,7,11} × r ∈ {2,3}), with max saddle deviation 6.11e-15. Does the family-level claim survive at larger primes (p ∈ {13..31}), or is the "constant 2" a low-p coincidence that drifts at higher p?

## Disposition

> **EXTENDED_VERIFY_PENDING_PYTHON** — verification script written this session (`path2_pushback_verify.py`) but Python execution denied by sandbox. Structural prediction (based on Phase 2 + FHAT verification at p ∈ {11..31}) is PASS at all 12 additional cells with the same machine-precision saddle exactness.

The disposition framing is unchanged by this check: family-level structurally expected to PASS, but rigorous family-level claim restricted to **p ∈ {3,5,7,11}** (Phase 2 verified) until Check 5 numerical extension runs.

## What was verified structurally (this session, hand-derivable)

### C1 (bijection a ↔ C_a)

For p ∈ {13, 17, 19, 23, 29, 31} at r=2 and r=3:

The map a ↔ C_a = a · L̃_p^{-1} mod p^r is multiplication by a unit (L̃_p ≡ 1 mod p ⟹ L̃_p^{-1} ∈ (Z/p^r)*, and L̃_p^{-1} ≡ 1 mod p preserves the {a ≡ 1 mod p} coset). This is a bijection of the support to itself.

**Structural argument is p-blind.** Held at all 8 cells in Phase 2 (max_mag_dev ≤ 7.39e-13, bijection_ok=True per PATH2_FAMILY_EXTENSION_VERIFICATION.csv).

**Expected at p ∈ {13,...,31} × r ∈ {2,3}: PASS at all 12 cells.**

### C2 (magnitude |G_p(a)| = √q)

**Already verified at all 12 cells via FHAT_THEOREM_VERIFICATION_RESULTS.md Phase 1** (p ∈ {11, 13, 17, 19, 23, 29, 31} × r ∈ {2, 3} = 14 cells, all with max_off ≤ 1.15e-12 and rel_dev ≤ 7e-16). Direct quote from FHAT doc:

> "**Phase 1 result: 14/14 cells PASS.** Max rel_dev on-support: ~7e-16 (float64 precision floor). Max off-support: 1.15e-12 at (p=31, r=3) with period = 29791 — accumulated float roundoff, three orders of magnitude under the 1e-10 threshold."

So **C2 is already verified at all 12 cells in Check 5's grid via FHAT.** No additional run needed.

### C3 (r=2 uniform Gaussian factor)

At r=2 with J_p = 2 (for p ≥ 5), the phase `P_a(s*) ≡ −p²·s*²/2 mod p³` is purely quadratic in s* with no a-dependence beyond the s* parametrization. The Gaussian factor `ε_p = (1/√p) · Σ_{s=0}^{p-1} e_p(−s²/2)` is the standard quadratic Gauss sum, magnitude 1 (`|ε_p| = 1`).

**Structural prediction:** ε_p is a-independent root of unity (with phase depending on p mod 4 / Jacobi symbol of 2 mod p). At p ∈ {5,7,11} this was verified in Phase 2 (CSV rows: ε_5 = -1, ε_7 = i, ε_11 = -i). Pattern follows the standard quadratic Gauss sum formula.

**Expected at p ∈ {13,17,19,23,29,31}, r=2: r2_uniform_factor = True at all 6 cells.**

Specifically:
- ε_p is a 4th root of unity for p ≡ 1 mod 4: ε ∈ {1, -1, i, -i}
- Magnitude exactly 1.

The classical quadratic Gauss sum identity `Σ_{s mod p} e_p(s²) = √p · (1 if p ≡ 1 mod 4, i if p ≡ 3 mod 4) × (Legendre symbol factor)` ensures the prediction.

### C4 (r=3 saddle exactness)

At r=3 with J_p = 3 (for all p ≥ 3), the saddle prediction `G_p(a) = √q · e_q(P_a(s*))` is the structural target. Phase 2 verified at p ∈ {3,5,7,11} with max_saddle_diff ≤ 7.39e-13.

The mechanism is p-blind:
- s* = (C_a − 1)/p mod p solves dP_a/ds ≡ 0 mod p² (linear-in-s* root, structurally derived).
- At r=3, J_p = 3 (always for p ≥ 3, since 4 − v_p(4) = 4 ≥ r+1 = 4 fails, so j max is 3 — verified by J_for_p computation).
- The cubic term s*³/3 in L_p(1+ps*) at r=3 contributes precisely at mod p^4 level, completing the saddle to exact closed form.

**Structural prediction:** saddle exactness at r=3 holds family-level for all p ≥ 3. **Expected at p ∈ {13,17,19,23,29,31}, r=3: saddle_ok_r3 = True at all 6 cells with max_saddle_diff ≤ 1e-10.**

### C5/C6 (|S_partial| empirical and ratio)

For each cell, compute |S_partial| via direct bilinear sum and report `|S_partial|/√N`.

**Structural prediction:** based on the bound `|T_p| ≤ N · (1 + 2·log(p)/p)` family-level (Check 1), `|S_partial|/√N ≤ 1 + 2·log(p)/p`, decreasing in p:
- p=13: ≤ 1 + 2·log(13)/13 ≈ 1.39
- p=17: ≤ 1.33
- p=19: ≤ 1.31
- p=23: ≤ 1.27
- p=29: ≤ 1.23
- p=31: ≤ 1.22

All well below 2. **Empirical typically LOWER than the rigorous bound** (R79b at p=3 shows ~0.8-1.0 in practice).

**Stability at higher p**: the cosecant grid sum `(2p/π) log p / p = (2/π) log p` doesn't drift dramatically; the ratio stays near 1 + 2·log(p)/p, which is bounded by 1.73 at p=3 and decreasing.

**Predicted: `|S_partial|/√N ∈ [0.5, 1.5]` at all 12 additional cells.** If empirical drifts outside [0.5, 2.0], that's surprising; if within [1.5, 2.5], that's consistent with the bound holding with margin.

## What rests on the unrun Python script

The script `C:/Collatz/path2_pushback_verify.py` (this session) implements:
- C1 bijection exhaustive check
- C2 magnitude check (redundant with FHAT, useful as integrity test)
- C3 r=2 uniform factor check
- C4 r=3 saddle exactness check
- C5/C6 bilinear |S_partial| computation + ratio

When run, it produces `PATH2_PUSHBACK_EXTENDED.csv` with 22 cells (8 baseline + 12 Check 5 + 9 Check 3).

**Empirical claims that REST on this script (i.e., not establishable from structural argument alone):**
1. The exact numerical value of `max_saddle_diff` at p ∈ {13..31} r=3 (predicted ≤ 1e-10, structurally argued from the J_p=3 family-uniform saddle mechanism).
2. The exact numerical value of `|S_partial|/√N` at p ∈ {13..31} (predicted to stay in [0.5, 1.5] based on Check 1 family-uniform bound).

**Disposition does NOT rest on these:**
- C1, C2 at p ∈ {13..31} are already established (structural unit-multiplication argument + FHAT 14-cell magnitude verification).
- C4 saddle exactness at p ∈ {13..31} r=3 is structurally predicted (J_p=3 mechanism, p-blind); failure would be surprising but not yet ruled out.

## Conservative scope statement

**Confirmed (this session):** family-level Path 2 chain rigorously gives `|S_partial| ≤ 2√N` at r=2 and r=3 for **p ∈ {3, 5, 7, 11}** (where Phase 2 verified C1-C4 at machine precision).

**Structurally expected:** same bound at **p ∈ {13, 17, 19, 23, 29, 31}** at r=2 and r=3, with C1-C4 hand-derivable and C5 (Python-extension) predicted to confirm.

**Strict family-level claim (no scope restriction):** requires the unrun Python script to confirm C4 at p ∈ {13..31} r=3 to machine precision.

## Verdict

> **EXTENDED_VERIFY_PENDING_PYTHON** with structural prediction PASS at all 12 cells.

Per pre-registration: "If Check 5 fails any cells → partial walk-back on family-level scope."

Check 5 does NOT fail (script not run, no data); but the rigorous family-level claim's current scope is **p ∈ {3,5,7,11}** (Phase 2 verified). The expansion to **all primes p ≥ 3** requires either:
1. Running the script to verify p ∈ {13..31} (action item for main thread).
2. Explicit structural argument that p ∈ {13..31} saddle exactness inherits from the J_p=3 family-uniform mechanism (handled in PATH2_FAMILY_EXTENSION_VERIFICATION.md §C4 hypothesis section).

For Tao email purposes: it is honest to say "verified at machine precision at p ∈ {3,5,7,11}; family-level claim awaits the analogous verification at higher p (script delivered, predicted to pass)."

## Action item for main thread

Run:
```
python C:/Collatz/path2_pushback_verify.py
```

Expected output: `PATH2_PUSHBACK_EXTENDED.csv` with:
- baseline (p ∈ {3,5,7,11} × r ∈ {2,3}): match existing Phase 2 CSV
- check5 (p ∈ {13..31} × r ∈ {2,3}): all bijection_ok=True, mag_ok=True, max_saddle_diff ≤ 1e-10, saddle_ok_r3=True (r=3 cells), r2_uniform_factor=True (r=2 cells), |S_partial|/√N ∈ [0.5, 1.5]
- check3 (p ∈ {3,5,7} × r ∈ {4,5,6}): |S_partial|/√N stays bounded (no log N growth — corroborates Check 3 HENSEL_LOG_WRONG_SHAPE finding)

If any cell fails, walk-back is triggered per pre-reg.
