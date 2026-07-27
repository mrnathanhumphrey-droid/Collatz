# RESULT — P6J: ENRICHED hemisphere — the 3|k channels CONVERGE; all residual is depleted (2026-07-26)

**Probe:** `probe_p6j.py`. Split the channel sum `T_i = Σ_k4⁻ᵏγ_i(k) = ½S_{i+1}` by parity of the base-4 lag:
`T⁺` (same-parity, 3|k, enriched, M₊=5/3) vs `T⁻` (cross-parity, 3∤k, depleted, M₋=2/3 — where the cascade/d₁ live).
`γ_i(k) = 3ⁱ⟨ρ_i, shiftₖρ_i⟩` from the certified numerator profile (built cheaply, no dense matrix).

## GATE + E-A — all exact
`build_rho4 == build_level['rho']` to **0.0** (i=2..6). Pinned mean `(1/3)M₊ + (2/3)M₋ = 1.00000000` every level
(M₊=5/3, M₋=2/3 exact — the average of γ over all k is exactly 1, forced). Enriched channels reproduce the certified
constants: `γ_∞(3) → 1.2370` (cert 1.2372), `γ_∞(6) → 1.3716` (1.3717), `γ_∞(9) → 2.1115` (2.112).

## E-B / E-D — the hemispheres separate; enriched converges, depleted is the whole residual
`T⁺ + T⁻ = 0.235676 = T₁₅` (full) — the redundancy check passes (E-D). But the two behave oppositely:

| | i=8 | i=11 | i=13 | i=15 | Λ (i=15) |
|---|---|---|---|---|---|
| **T⁺ (enriched)** | 0.019709 | 0.019678 | 0.019672 | **0.019672** | **−1e-7 (≈0, converged)** |
| **T⁻ (depleted, ours)** | 0.213620 | 0.214793 | 0.215456 | 0.216004 | **+2.6e-4 (rising)** |

- **T⁺ CONVERGES fast** — `Λ⁺ → 0` by i≈13 (Λ⁺: −5e-5 at i=8 → +2e-7 at i=13), limit **0.019672**, monotone from above.
  (= `Σ_{3|k}4⁻ᵏγ_∞(k) = 1.2372/64 + 1.3717/4096 + 2.112/4⁹ + … ≈ 0.01967`.)
- **T⁻ carries the entire slow residual** — rising, `Λ⁻ ≈ +2.6e-4`, deparitied two-step rate settling at **0.91** (i=13,14,15
  → 0.905, 0.914, 0.906), identical to P6I's whole-system rate. So the P6I residual **is** the depleted hemisphere.

Decomposition of the limit: `S_∞/2 = T⁺_∞ + T⁻_∞ ≈ 0.01967 (enriched, done) + 0.2186 (depleted, rate 0.91) = 0.2383`
→ `S_∞ ≈ 0.4766`. Consistent with P6H/P6I.

## E-C — looking south shows no turn: the enriched side is closed
`T⁺` has **already converged** to its limit (0.019672) with `Λ⁺ ≈ 0` and a monotone, turn-free approach — it is not "nearer
a turn," it is **done**. So the enriched hemisphere gives **no early-warning signal** of a system-wide turnover; whatever
happens to 7/15 must happen in the **depleted 3∤k channels alone**, which is exactly where the slow rising residual lives.
This *isolates* the open sign question: the enriched half is settled and positive-contributing; the depleted half is the
entire unresolved object. Not at stake: P6D–P6I identities, S_{i+1}=2T_i, the value bracket [0.4714, 0.478], R1–R30.
