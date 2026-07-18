# TAUBERIAN_SCOPING_MATCH — Phase 3: matching singularity to theorem

**Date:** 2026-05-12. Wilson. Phase 3 of the Tauberian scoping probe.

## Purpose

Cross-reference Phase 1's empirical singularity detection (R77.6 + cached ε_n) against Phase 2's theorem hypotheses. For each candidate Tauberian theorem applicable in principle, determine whether its singularity-type hypothesis matches what R77.6 detected, what asymptotic it predicts, and whether the prediction is empirically supported through k=6.

---

## Matching table

Candidate theorems from Phase 2 (Chevalier 2507.15394 and Flajolet-Sedgewick Ch. VI passed initial filter as power-series, no non-negativity required). Below: each theorem's singularity-type hypothesis vs R77.6's detected structure.

| Theorem | Hypothesis on E(z) | Matches R77.6 (branch cut at z=2)? | Predicted ε_n asymptotic | Empirical support (n=2..6) |
|---|---|---|---|---|
| **Chevalier Thm 1.14** (square-root branch) | g(z) = h(√(1−z/2)), h holomorphic in neighborhood of D̄(1,1)^{1/2}, h'(0) ≠ 0 | **TYPE-CONSISTENT**, branch order specifically α=1/2 | ε_n ~ C·(1/2)^n · n^{-3/2} | **NOT SUPPORTED**; see below |
| **Chevalier Thm 1.16** (meromorphic h with pole of order M at 0) | g(z) = h_p(√(1−z/2)), h_p meromorphic, pole at 0 of order M ≥ 1 | TYPE-CONSISTENT; with M=1 gives the boundary case | ε_n ~ D·(1/2)^n · n^{M − 3/2} | M=1: ε_n ~ (1/2)^n · n^{-1/2}; **also not supported** |
| **Chevalier Rem 1.15** (general α-branch) | g(z) = h((1−z/2)^α), α ∈ (0,1) | TYPE-CONSISTENT for any α ∈ (0,1) | ε_n ~ C·(1/2)^n · n^{-α-1} | Could fit ANY α; not informative at small N |
| **FS Ch. VI** (singularity analysis, full table) | f(z) has singular expansion at z=2 from the FS singularity table | TYPE-CONSISTENT; framework, not single theorem | Determined by the specific singular expansion fit | Depends on the specific expansion |
| **FS Ch. VI: simple pole** | f(z) ~ c/(1−z/2) | NO — R77.6 explicitly rules out simple pole at z=2 | ε_n ~ c·(1/2)^n (constant prefactor) | **PARTIALLY** matches leading-order |
| **FS Ch. VI: log singularity** | f(z) ~ c·log(1/(1−z/2)) | TYPE-CONSISTENT (R77.6 cannot rule out log) | ε_n ~ c·(1/2)^n · 1/n | Mixed support |
| **FS Ch. VI: pole + log** | f(z) ~ c·log(1/(1−z/2))/(1−z/2) | TYPE-CONSISTENT (R77.6 cannot distinguish at N=5) | ε_n ~ c·(1/2)^n · log(n)/n? | Not directly testable |

---

## Empirical verification at n=2..6 (the key disambiguation)

The cached ε_n values give a clean empirical test. From Phase 1 §(e):

| n | |ε_n|·2^n |
|---|---|
| 2 | 0.038095 |
| 3 | 0.040736 |
| 4 | 0.039236 |
| 5 | 0.036856 |
| 6 | 0.031866 |

The empirical |ε_n|·2^n is **nearly constant**, fluctuating in [0.032, 0.041] around the value 1/30 = 0.03333 (R76 §10's conjectured leading coefficient).

### Test against Chevalier Thm 1.14 (α=1/2 square-root): prediction ε_n ~ C·(1/2)^n · n^{-3/2}

If this prediction were correct, |ε_n|·2^n · n^{3/2} should be approximately constant. Computed values:

| n | |ε_n|·2^n · n^{3/2} |
|---|---|
| 2 | 0.1077 |
| 3 | 0.2117 |
| 4 | 0.3139 |
| 5 | 0.4121 |
| 6 | 0.4683 |

These are **NOT** constant — they grow roughly linearly (range = 4.3× across the data). Chevalier α=1/2 prediction **FALSIFIED** as a description of the leading-order ε_n behavior.

### Test against Chevalier Thm 1.16 (meromorphic h, M=1): prediction ε_n ~ D·(1/2)^n · n^{-1/2}

|ε_n|·2^n · n^{1/2} should be approximately constant.

| n | |ε_n|·2^n · n^{1/2} |
|---|---|
| 2 | 0.0539 |
| 3 | 0.0706 |
| 4 | 0.0785 |
| 5 | 0.0824 |
| 6 | 0.0781 |

Range = 1.5× across data. **Better than n^{-3/2}**, but still not constant. M=1 not strongly supported either.

### Test against simple-pole reading (α=0 effectively): prediction ε_n ~ c·(1/2)^n (no n-correction)

|ε_n|·2^n should be approximately constant.

From the table above: range [0.032, 0.041], a factor of 1.28×. This is the **tightest** fit among the tested hypotheses.

### Test against logarithmic singularity (FS Ch. VI log entry): prediction ε_n ~ c·(1/2)^n / n

|ε_n|·2^n · n should be approximately constant.

| n | |ε_n|·2^n · n |
|---|---|
| 2 | 0.0762 |
| 3 | 0.1222 |
| 4 | 0.1569 |
| 5 | 0.1843 |
| 6 | 0.1912 |

Range = 2.5× — **not constant**.

### Implied exponent β fit: |ε_n|·2^n ~ C/n^β

From log-log slopes on consecutive pairs (from Phase 1 §(e)):

| n_1, n_2 | implied β |
|---|---|
| 2 → 3 | −0.165 (negative: would grow with n) |
| 3 → 4 | +0.130 |
| 4 → 5 | +0.281 |
| 5 → 6 | +0.798 |

The implied β grows from negative to nearly 1 across the data — no stable exponent emerges. This is **inconsistent with any pure-power-law asymptotic**, including α=1/2 (which would require β = 3/2 stably).

---

## Refined interpretation

The empirical pattern (constant ~1/30, slow upward drift then downward decay) is consistent with the **decomposition**:

ε_n = (leading-rate component, simple-pole-like) + (subleading branch-cut correction at the same rate).

Specifically, if E(z) has both a **simple-pole-like behavior** at z=2 with residue (1/30)·2 = 1/15 (matching R76 §10's 1/30 leading coefficient), PLUS a **branch-cut endpoint** at z=2 contributing a slowly-decaying correction, then:

- Leading: ε_n ≈ −(1/30)·(1/2)^n. This dominates and explains the near-constancy of |ε_n|·2^n at small n.
- Subleading: deviation (|ε_n|·2^n − 1/30) ≈ correction term decaying with n.

This is consistent with FS Ch. VI's **mixed singularity** picture: E(z) = c/(1 − z/2) + branch-cut term at z=2 of unknown precise form.

R77.6's "no simple pole" detection is reading the *subleading* branch structure (the Padé approximant's pole drift is sensitive to the branch part), not denying the dominant simple-pole behavior.

### Test: subleading correction analysis

Define δ_n := |ε_n|·2^n − 1/30. Then δ_n is the subleading correction normalized to (1/2)^n.

| n | δ_n |
|---|---|
| 2 | +0.00476 |
| 3 | +0.00740 |
| 4 | +0.00590 |
| 5 | +0.00352 |
| 6 | −0.00147 |

δ_n is **non-monotone** and changes sign between n=5 and n=6. This non-monotonicity is the strongest empirical evidence we have, and it points to:

1. The subleading correction is **not a single decaying power**. A simple n^{-α-1} would be monotonically decaying and same-signed throughout.
2. The non-monotonicity suggests either:
   - **Multiple subleading terms** of competing signs (e.g., one term in Chevalier's expansion at n=2..6 dominated by c_1/n correction, the next dominated by c_2/n^2, etc.).
   - **Oscillatory contributions** from secondary singularities (e.g., a pole or branch on the second sheet of E(z), or a complex-conjugate pair at z = 2·e^{±iθ}).
   - **Subleading log structure** (FS log·power entry) where leading c_k·log(n)^k changes character at small n.

3. The empirical sign at n=6 going negative is consistent with an asymptotically-negative correction term that emerges after a transient. This is the standard pattern when an asymptotic series has multiple competing terms at the boundary of where the leading-order approximation breaks down.

### N=5 limitation acknowledged

R77.6 itself states: *"The 5-coefficient budget cannot separate (G-power) from (G-log)."* Phase 3 confirms this is **also true** for the empirical asymptotic verification:

- Six data points cannot robustly distinguish n^{-3/2} from n^{-1} from n^{-2} from constant when the leading-order term dominates.
- The non-monotonicity of δ_n means we cannot do a clean log-log fit to the subleading; the fit is *unstable* across consecutive pairs.

---

## Disposition for matching

### Strict reading: Which Chevalier or FS hypothesis hold for E(z)?

**Chevalier Thm 1.14 (pure α=1/2 branch):** Empirically **falsified** as a description of the leading-order ε_n behavior. The prediction n^{-3/2} decay is not seen.

**Chevalier Thm 1.16 (meromorphic h with M=1):** Predicts n^{-1/2}; better fit than M=0 (Thm 1.14) but still not matching the near-constancy. Could be reconcilable if the leading constant from the residue dominates the n^{-1/2} subleading term.

**Chevalier Rem 1.15 (general α):** The "best fit" α would have to be near 0 (giving near-constant behavior); but α=0 falls outside the (0,1) range Remark 1.15 covers — at α=0, the result degenerates to a simple pole.

**FS Ch. VI: simple pole + branch-cut correction:** Best qualitative match. Leading-order simple pole (residue 1/30 in the rate-1/2 normalization) explains the near-constancy; subleading branch-cut at z=2 of indeterminate type explains the small deviations and Padé pole drift.

### Honest reading: We need more data

The fundamental issue is that **at N=5 (six ε_n values), the leading-order behavior is constant, and the subleading correction is too small and non-monotone to identify its precise form.**

This matches R77.6's own honest verdict: *"The branch type (power-law vs logarithmic) is not resolvable from N = 5 coefficients."*

### Choice of best match

If forced to pick one Tauberian theorem as "the right tool for E(z) at z=2":

**Flajolet-Sedgewick Ch. VI's singularity-analysis framework**, with the specific expansion being a *mixed* singularity at z=2 — leading-order simple-pole-like, with a subleading branch-cut correction whose precise form (power, log, or hybrid) is indeterminate from N=5 data.

This is **NOT** a single Tauberian theorem in the form "hypothesis → conclusion"; it is the **framework** that handles the mixed structure. Within that framework:

- **Chevalier 2507.15394 Thm 1.16 (meromorphic h)** is the cleanest single-theorem candidate if we read R76 §10's "+O((1/4)^n)" as "+leading pole-like residue at z=2 contributing main rate, plus a subleading branch correction at z=2 captured by the meromorphic h with M ≥ 1 pole at 0."

- **Chevalier Thm 1.14 (pure α=1/2)** is the **wrong single theorem** because the leading order is not √-branch.

---

## What additional data would disambiguate

To clearly distinguish among candidates, the following would help:

1. **Extending to ε_7 (k=7 Markov chain).** Adds one diagonal Padé point [3/3]. Predicted by R77.6: pole at z ≈ 2.030–2.040.
2. **Extending to ε_8 (k=8 Markov chain).** [4/4] diagonal point.
3. **Differential approximants (D-Padé).** Could distinguish power-law exponent at lower N, but requires additional infrastructure.
4. **Direct analytic continuation of E(z) past z=2.** If achievable in closed form (unlikely given Collatz dynamics, but worth flagging), the branch structure would be directly visible.
5. **Numerical evaluation of E(z) on the second sheet near z=2.** Would reveal whether the cut is "real" branch (power-law) or logarithmic.

---

## Output for Phase 4

Phase 4 will:
- Document the **falsification** of the H_SQUARE_ROOT_MATCHES_PLUS_EMPIRICAL hypothesis: Chevalier Thm 1.14's pure α=1/2 prediction is empirically not seen at n=2..6.
- Document the **structural reading**: FS Ch. VI's mixed-singularity framework is the right level of abstraction; the exact theorem that closes the analysis requires either ε_7+ data or a different structural input.
- Surface H_AMBIGUOUS as the honest disposition given the N=5 limitation.
