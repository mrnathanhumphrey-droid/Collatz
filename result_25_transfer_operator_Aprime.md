# Result 25 — ROUTE 1A′ (weighted cascade operator): the operator is GATE-VALIDATED; `λ₁ = 1/3` confirms R5 as the top eigenvalue; q=3 divergence = DEGENERATE top eigenvalue. But `r_q` for q≥5 is NOT cleanly `λ₂/λ₁` (tower-entangled + unconverged in L).

**Date:** 2026-07-16. **Type:** probe (`probe_25_transfer_operator_Aprime.py`). **Verdict: the transfer operator is real and validated; it confirms R5 and the q=3 mechanism; it does NOT yet pin `r_q` for q≥5.**

## The operator

State `s = (a, b, γ)`, `a = 2^{-S} mod q^L`, `b = 2^{-S'} mod q^L` (both in `H = ⟨2⟩ mod q^L`), `γ ∈ Z/q^L`. One transfer step prepends a coordinate to each address (multipliers `2^{-Δ}, 2^{-Δ'} ∈ H`, folded-geometric weights), keeps the branch iff the digit condition `(γ+T) ≡ 0 mod q` holds (`T = a'−b'`), and updates `γ' = (γ+T)/q`. Built empirically (one step enumerated exactly — no hand-derived matrix entries). `M` is nonnegative sub-stochastic; `v₀ = 𝟙(a=1,b=1,γ=0)`; `sum(Mᵏ v₀) = P(collide to depth k) = ‖π_k‖²`, exact for `k ≤ L`.

## G_GATE — PASS (the operator is correct)

`sum(M_L^k v₀)` equals float `‖π_k‖²` (from `stationary`) to <1e-9 rel for **all k ≤ L, every q** (machine precision). E.g. q=3 L=3: k=1,2,3 match to ≤2.9e-16. The construction is faithful — this is the load-bearing validation.

## Win 1 — `λ₁(M) = 1/3` confirms R5 as the top eigenvalue

| q | L | λ₁ |
|---|---|---|
| 3 | 3 | 0.333336 |
| 5 | 2 | 0.333334 |
| 7 | 2 | 0.333334 |

R5's `‖π_k‖² ~ 3^{-k}` (the universal "3", `D₂ = log3/logq`) **is literally the Perron eigenvalue** of this validated operator. Independent confirmation of R5 from a completely different construction.

## Win 2 — q=3 divergence = DEGENERATE top eigenvalue (mechanism confirmed)

`λ₂/λ₁` at q=3 converges cleanly in L: **0.600 (L=1) → 0.9916 (L=2) → 0.9997 (L=3) → 1**. It reaches 1 because **`λ₂ → λ₁ = 1/3`** — the top eigenvalue is (asymptotically) degenerate. A degenerate Perron eigenvalue / Jordan block gives `k·(1/3)^k` growth ⇒ `cross(k)` linear ⇒ divergence. **`r₃ = 1 ⟺ λ₁ degenerate`, exactly as pre-committed.** The q=3 mechanism is now an eigenvalue statement.

## The block — `r_q` for q≥5 is NOT cleanly `λ₂/λ₁`

Two independent problems:
1. **Unconverged in L.** q=7: `λ₂/λ₁ = 0.655 (L=1) → 0.981 (L=2)`; q=5: `0.728 → 0.979`. Two truncation levels is not convergence, and L=3 for q≥5 exceeds the state-count cap (q=5 L=3 ≈ 1.25M).
2. **Tower entanglement.** The subdominant spectrum is a *cluster* near 1/3 (q=5 L=2 top |eigs|: 0.333, 0.326, 0.326, 0.308, 0.308, 0.285). These are the within/tower modes (R24: the tower dominates the near-top spectrum). The cross-specific rate is buried among them, not isolable as plain `λ₂`. Extracting it needs an **amplitude-resolved** decomposition (project `v₀` and the summation functional onto eigenvectors) or a within-subtracted operator.

## The revision this forces

The operator hints `r_q` for q≥5 may be **much closer to 1** than the "~0.6" implied by low-k `ρ_k`. This is **consistent with R22's `ρ₅` climbing** (0.534 → 0.508 → 0.624, still rising). If `r₅ ≈ 0.9+`, then q≥5 converges but *slowly*, and the `3/q = 0.6` prior was wrong on magnitude, not just value. **Honest status: `r_q` (q≥5) ∈ (0.62, ~1), rising, NOT pinned.**

## Plan status after R25

| phase | status |
|---|---|
| 3 — the bound | **operator BUILT + GATE-VALIDATED**; `λ₁=1/3` = R5; q=3 divergence = degenerate top eigenvalue (`r₃=1` confirmed as an eigenvalue statement). **`r_q` for q≥5 OPEN** — `λ₂/λ₁` unconverged in L + tower-entangled; needs amplitude-resolved extraction, higher L, more exact `ρ_k` (q=5 k=6), or route B (Fourier, mode-separated). |

## Not at stake
R10–R24 (R5 confirmed, R22 `ρ` rising vindicated), R5's rate, R6, R7, R12, THEOREM_C_745.

_Reporting discipline: the gate pass is stated as the load-bearing win (the operator is real); `λ₁=1/3` and the q=3 degenerate-eigenvalue mechanism are reported as confirmed; the failure to extract `r_q` for q≥5 is stated plainly with both causes (unconverged L + tower entanglement), NOT spun as a near-miss; and the "r_q may be ≈1" revision is flagged as a hint consistent with R22's rising ρ, not asserted as pinned._
