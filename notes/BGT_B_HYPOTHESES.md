# BGT_B — Karamata's Tauberian theorem

## Phase 0 — verbatim statement

**Source:** Kevei 2019 Theorem 18 (= BGT 1.7.1; Feller XIII §5). Lines 639-649 of `Kevei_2019_Regularly_Varying_Functions_Notes.txt`.

> **Theorem 18.** Let U be as above [U nondecreasing right-continuous on R, U(x)=0 for x<0], c ≥ 0, ρ ≥ 0, ℓ ∈ SV. The following are equivalent:
> (i) U(x) ~ c · x^ρ · ℓ(x) / Γ(1+ρ) as x → ∞;
> (ii) Û(s) ~ c · s^{−ρ} · ℓ(1/s) as s ↓ 0.

Where Û(s) = ∫_[0,∞) e^{−sx} dU(x) is the Laplace-Stieltjes transform.

## Hypothesis types

- h_1: U nondecreasing, right-continuous, vanishes for x < 0 (sign constraint).
- h_2: U(x) is regularly varying with index ρ ≥ 0 (or equivalently Û(s) is RV at 0).
- h_3 (load-bearing): KNOWN Laplace-side asymptotic Û(s) ~ c · s^{−ρ} · ℓ(1/s) as s ↓ 0.

## Phase 1 — hypothesis × ε_k matrix

For ε_k, the natural U is the partial sum S_n = Σ_{k≤n} ε_k. S_n converges (numerically near 0.198), so U(n) is a *bounded* sequence, not RV with positive index.

| hyp | check | verdict |
|---|---|---|
| h_1 | S_n must be nondecreasing — ε_k has negative signs for k ≥ 3, so S_n is non-monotone (max at n=2 of 0.2095, then declining) | **FAILED**. |
| h_1 (variant) | take Σ |ε_k| — this is nondecreasing, but the absolute partial sum doesn't deliver the chain-side quantity needed for closure | SATISFIED for variant; but conclusion-shape doesn't deliver closure. |
| h_2 | S_n is bounded (converges to ~0.198), so RV with ρ ≥ 0 requires ρ = 0 (slow variation of bounded function); but slow variation of a bounded function requires lim S(n) = const, satisfied | SATISFIED at ρ = 0 trivially. |
| h_3 | Laplace-side asymptotic Û(s) ~ c · ℓ(1/s) as s ↓ 0 with ρ = 0 — this is essentially constant = lim S_n = 7/45·3 = 7/15 (if S_n target identified) | **MODE_H_CIRCULAR** — the Laplace-side asymptotic IS the target c=7/45 result we're trying to prove. |

**Phase 1 verdict: FAILED (h_1 non-monotone) + MODE_H_CIRCULAR (h_3 = target).** Even if we work with |ε_k|, the conclusion at ρ=0 is just the trivial "constant ~ constant" — no Fourier-decay content.

## Phase 2 — conclusion shape

The Tauberian conclusion at ρ = 0 with ℓ slowly varying gives U(x) ~ c·ℓ(x). For S_n converging to a constant, ℓ → c trivially. No rate information — the Tauberian theorem doesn't deliver the (1/2)^n decay rate to S_n from any side. To get a *rate* (i.e., S_n − 7/15 = O((1/2)^n)), one would need a Tauberian theorem with REMAINDER, which Theorem 18 does not supply.

**Phase 2 verdict: SHAPE_MISMATCH.** Conclusion is bare-asymptotic, no rate.

## Phase 3 — multi-regime check

The k=7 jump is irrelevant if we're at ρ=0 trivial Tauberian — but it IS relevant if one tries to bootstrap to a remainder theorem (where ε_k itself is the object). For an "abel-summability + remainder" theorem (Hardy-Littlewood + 2RV refinement), the k=7 jump again kills the second-order slow-variation hypothesis (same reason as candidate A Phase 3).

**Phase 3 verdict: STRUCTURALLY_BLOCKED** (for the rate / remainder version) and **N/A** (for the bare ρ=0 version since the bare version doesn't deliver closure).

## Disposition: MODE_H_CIRCULAR + NO_FIT

Bare Karamata Tauberian doesn't deliver rate. Refined-remainder version (2RV) has the same k=7 obstacle as candidate D. Either way, the Laplace-side asymptotic IS the closure target after Plancherel — Mode H circular.
