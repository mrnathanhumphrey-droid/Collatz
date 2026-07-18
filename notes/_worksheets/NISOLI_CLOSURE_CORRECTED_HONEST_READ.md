# NISOLI_CLOSURE_CORRECTED_HONEST_READ — Phase 4 on A's extractability

**Date:** 2026-05-12. Phase 4 of Route A.

---

## 1. The question

Phase 3 found the closure inequality fires at (r=3, A=3, K=6) under the optimistic `C_A = 1` parameterisation. But the inequality's load-bearing parameter is `A` (and its companion constant `C_A`) from Tao Prop 1.17. **Can A — with a usefully-bounded C_A — be extracted from the project's machinery?**

This is the honest gate. Even if (r=3, A=3, K=6) is mathematically a firing cell, the closure is rigorous only if a specific A and C_A are deliverable.

---

## 2. Sources of A — Phase 4 inventory

### Source (a): Tao Prop 1.17 effective bookkeeping

**Status:** **INFEASIBLE** per BOOKKEEPING_PHASE1_DISPOSITION (2026-05-11).

The iterated-cubic recursion in Case 3 of Tao Prop 7.8's proof (the `p_{i+1} ≤ 40A(1+p_i)³ + O(A)` iterated R = A²/ε times) forces:

- Optimistic floor: `C_A ≥ A^{O(A)}` (super-polynomial in A).
- Faithful bookkeeping: `C_A ≥ exp(exp(A²))`.

At K = 10, no A value produces a C_A small enough to satisfy Nisoli `η < 1` with M_3'' ≈ 24, **even under the optimistic A^{20A} floor** (requirement fails by 40+ orders of magnitude per BOOKKEEPING_PHASE1 §"override check 3").

**Source (a) is structurally closed**, not "open with effort." Tao's renewal-process method does NOT produce a polynomial-in-A bound, regardless of bookkeeping quality.

### Source (b): Cross-freq machinery directly

The cross-freq derivation (CROSS_FREQ_DISPOSITION + T_LEAD_CORRECTED_PHASE1) gives `W_+(g) = 2^{−g+1}/15` for g ∈ {2, 4, 6, ...}, with geometric tail. Σ_g W_+(g) = 2/45 (the Off_lin contribution to T_lead's eigenvalue).

**Is this an "A"?**

No. The cross-freq A_eff in the sense of geometric decay in `g` is the decay of the off-diagonal cross-frequency weight as g grows — i.e., how the bilinear sum's cross-frequency contributions decay as the geometric-step difference grows. This is **structurally a different quantity** from Tao's A:

- **Cross-freq decay:** in `g = v' − v`, where v, v' are geometric step indices. Lives in T_lead's matrix structure.
- **Tao's A:** in `n`, where n is the Syracuse modulus level. Lives in `|μ̂_n(ξ)|`'s decay.

The Nisoli closure inequality `|K_bil(r)| · K^{−A} · M_3'' < 1` uses Tao's A, not cross-freq's geometric-in-g decay. **The two are not interchangeable.**

**Source (b) does NOT deliver A**, but it does deliver M_3'' (which it did, exactly: 24.43).

### Source (c): T_lead's eigenvalue 0 on (9, −7)

T_lead is rank-1 with spectrum {43/45, 0}. The zero eigenvalue means: any component along (9, −7) is killed in one step. This is sometimes informally described as "infinite A for this mode."

**Is this a useful Tao A?**

No, for two reasons:

1. The zero eigenvalue is on (9, −7), the **null direction** orthogonal-like to the (1, 4) eigenvector. The (1, 4) direction (which carries the rate-bearing data) has eigenvalue 43/45 — i.e., decays at rate 43/45 ≈ 0.956 per step, NOT geometrically fast.
2. T_lead's spectral structure is **within-level** (one-step iteration of the 2×2 operator). Tao's A is **across-level** (asymptotic decay of `|μ̂_n(ξ)|` as n → ∞). The operator-iteration semantics of T_lead at one n doesn't translate to Tao's across-n decay (cf. T_V_DISPOSITION's H_M_RECURSION_UNDERSPECIFIED finding).

**Source (c) does NOT deliver a usable A.**

### Source (d): New Tao Prop 1.17 analogue from cross-freq framework

Could cross-freq's framework provide a NEW bound on `|μ̂_n(ξ)|` independent of Tao §7?

The cross-freq machinery establishes:
- Plancherel decomposition of `|μ̂_n(ξ)|²` summed by frequency class. (This is fact-of-Plancherel, not a decay bound.)
- Closure of the bilinear recursion on V_M = span{M_n^{ab}(g, c) : g ∈ {0, 2, 4, ...}}.
- T_lead's spectrum {43/45, 0} on the (1, 4)-projection of V_M.

None of this implies a per-frequency upper bound `|μ̂_n(ξ)| ≤ C · n^{−A}`. The within-level Plancherel sum `S_n = Σ_ξ |μ̂_n(ξ)|² → 7/15` gives a per-class total but **not a per-frequency decay rate** like Tao's A.

A polynomial-decay `|μ̂_n(ξ)| ≤ n^{−A}` is genuinely different machinery from cross-freq's V_M closure; in principle one could imagine a result deriving it from cross-freq's bilinear structure, but **this is not in the project today** and is exactly the "novel technique outside Tao's method" gap flagged in BOOKKEEPING_PHASE1_DISPOSITION §"Alternative routes."

**Source (d) is the only conceivable unblocker, but it's a research project, not an extraction.**

---

## 3. The honest disposition

> **A cannot be extracted from the project's machinery in any form that delivers the Nisoli closure inequality at the corrected rate λ = 43/45.**

This is consistent with — and a generalisation of — the rate-1/2 R77.2 finding (Tao C_A INFEASIBLE) and BOOKKEEPING_PHASE1_DISPOSITION (the Tao-bookkeeping obstacle is structural, not a quality-of-extraction problem).

**Specifically:**
- Source (a) is closed (BOOKKEEPING_PHASE1's INFEASIBLE).
- Source (b) doesn't deliver A (cross-freq decays in g, not in n).
- Source (c) doesn't deliver A (T_lead's null eigenvalue is within-level, not asymptotic).
- Source (d) is a research project, not an extraction (and would require novel Fourier-analytic technique outside Tao's renewal-process method).

**Therefore the corrected-rate Nisoli closure inequality is unevaluable under any rigorous A.**

---

## 4. What this means for the corrected-rate framework

The T_LEAD_CORRECTED probe landed H_T_LEAD_CARRIES_DIFFERENT_RATE — a positive structural finding that c=7/45 has a within-level cross-freq-derived eigenvalue 43/45 over Q. **This finding stands independently of the Nisoli closure status.**

What changes:

- The seven-probe trajectory's "Nisoli structurally inapplicable at rate-1/2" framing **does not flip** at the corrected rate. The structural obstruction is now at a different layer:
  - At rate-1/2: no discrete eigenvalue exists (R77.3 falsified the 3-mode model; R77.4 K_k has no λ near 1/2).
  - At rate-43/45: discrete eigenvalue exists (T_lead's 43/45 is clean Q-rational) — BUT Tao Prop 1.17's A is INFEASIBLE for the closure to fire rigorously.

- **The bequest of the eighth probe**: "T_lead has a clean Q-rational eigenvalue 43/45 with explicit closed-form origin, and the resolvent norm M_3'' = 24.43 is small enough that the closure inequality would fire at (r=3, A=3, K=6) IF Tao Prop 1.17 delivered polynomial-in-A constants. Tao's method does not deliver this; the structural Nisoli closure remains unavailable."

This is **honestly halfway** — the within-level structural finding stands, but the rigorous closure does not follow.

---

## 5. Comparison: where the rate-1/2 vs corrected-rate probes land

| Layer | rate-1/2 (R77.3 + M3_DISPOSITION) | corrected rate 43/45 (this probe) |
|---|---|---|
| Discrete eigenvalue exists | NO (3-mode model falsified) | YES (T_lead's 43/45 over Q) |
| M_3 finite and computable | YES (≈ 100-1000 vs T_3) but operator-theoretically hollow | YES (= 24.43 exact, anchored on T_lead) |
| Closure inequality fires at small (A, K) under `C_A = 1`? | YES at (r=3, A=3, K=6) if M_3 = 100 | YES at (r=3, A=3, K=6); product 0.679 |
| Tao Prop 1.17 C_A extractable? | NO (BOOKKEEPING_PHASE1 INFEASIBLE) | NO (same obstruction) |
| Rigorous Nisoli closure achievable? | NO | NO |

The rate-1/2 disposition has an extra failure mode (operator-theoretic anchor falsified by R77.3). The corrected-rate disposition has the operator-theoretic anchor (T_lead is a clean structural object) but **the same Tao-A blocking obstruction**.

So the corrected rate IMPROVES the structural anchor but does NOT improve the closure feasibility — the load-bearing obstruction is shared.

---

## 6. Honest verdict

> **H_A_EXTRACTION_HARD.**
>
> The corrected-rate Nisoli closure inequality at λ = 43/45 fires at (r=3, A=3, K=6) under the optimistic `C_A = 1` parameterisation — structurally the SMALLEST firing cell, and a clean within-level marker of where the inequality wants to go. But the **A required (= 3) is not extractable from any project source**:
>
> - Tao Prop 1.17 effective bookkeeping: INFEASIBLE (BOOKKEEPING_PHASE1 §C_A super-exponential).
> - Cross-freq machinery: delivers decay in `g`, not in `n` (structurally different objects).
> - T_lead's null eigenvalue: within-level, not the asymptotic A Nisoli needs.
> - Novel polynomial-in-A bound outside Tao's method: research project; not on file.
>
> The closure inequality is **mathematically parameterised** (the (r=3, A=3, K=6) cell satisfies the form) but **structurally unevaluable** (no rigorous A delivery).

This is the honest disposition. The corrected-rate finding is a **positive structural anchor (43/45 over Q exists)** but **not a rigorous closure (Tao-A still infeasible)**.

---

## 7. Files

- BOOKKEEPING_PHASE1_DISPOSITION.md (Tao C_A INFEASIBLE)
- NISOLI_CLOSURE_CORRECTED_FIRING.md (Phase 3 firing cells)
- CROSS_FREQ_DISPOSITION.md (cross-freq machinery; doesn't deliver A)
- NISOLI_CLOSURE_CORRECTED_HONEST_READ.md (this file)
- NISOLI_CLOSURE_CORRECTED_STRUCTURAL.md — Phase 6 implication for seven-probe trajectory
