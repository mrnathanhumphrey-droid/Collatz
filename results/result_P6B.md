# RESULT — P6B: parity mechanism CORRECTED (3∤m→cross, not odd-lag); collapse inconclusive (2026-07-26)

**Probe:** `probe_p6b.py`. (1) Test the load-bearing parity mechanism via the odd/even-lag corollary — from data, not
assertion. (2) Test Wilson's cross-term collapse `ν_o = ½(×2⁻¹)_*ν_e`. Reuses P1/P6 machinery, per level.

## (1) MECHANISM — CORRECTED: the rule is 3∤m, not odd/even lag
Parity decomposition `A_j(m) = A[even] + A[odd] + cross` for m=1..6 (j=4 and j=6, identical pattern):
| m | verdict | | m | verdict |
|---|---------|--|---|---------|
| 1 | **CROSS** | | 4 | **CROSS** |
| 2 | **CROSS** | | 5 | **CROSS** |
| 3 | SAME | | 6 | SAME |

**`3∤m → cross-parity` (same-parity shells = 0); `3\|m → same-parity` (cross = 0).** Both prior mechanism claims are
wrong:
- ⚠️ **RETRACTED (Claude, result_P6):** "m=1 is odd-lag; odd-lag = cross, even-lag = same." FALSE — m=2 (even) is CROSS,
  m=3 (odd) is SAME. The odd/even-lag framing is dead.
- **Wilson's "lag-1 wants same-parity" (a′−a≡2)** also doesn't match — m=1 is cross.
- **The correct empirical law: A_j(m) is cross-parity iff 3∤m.** d₁ (m=1) is 3∤1 → cross-parity, so the P6 *finding*
  (d₁ = cross-parity cross-term) STANDS; only the *reason* was wrong.

**This aligns with the enriched/depleted dichotomy:** the same-parity channels are exactly `3|k` (enriched, class mean
5/3); the cross-parity channels are `3∤k` (depleted, class mean 2/3), where d₁ lives. So the branch-parity split IS the
dichotomy — the enriched channels are carried within-parity, the depleted (d₁ included) across-parity. Structurally:
`|ρ̂_e(a)|²+|ρ̂_o(a)|²` has m-channel support only on `3|m` (a constraint on the same-parity power spectrum). **The
mechanism is Wilson's pen — I will not assert a reason again; the empirical law is 3∤m → cross.**

## (2)(3) COLLAPSE — INCONCLUSIVE (wrong coordinate + a≥1 convention), not refuted
Testing `ρ_o[s] = c·ρ_e[(s−δ)%N]` (Wilson: c=½, δ=dlog(2⁻¹)) by roll-scan: best δ ≠ 2⁻¹ mod N, c ∈ {−1.12, 0.67,
−0.31}, **residual ~0.8** (huge). And the collapsed form `½·3^j⟨ρ_e,shift_{2⁻¹}ρ_e⟩` vs A_j(1): ratios −8.98, 0.53,
−4.44 (≠1). So the roll-based test **fails**. But two reasons it is likely the wrong test, not a refutation:
- **Coordinate:** ×2⁻¹ is a shift in the **base-2** dlog, but `2⁻¹ ≡ 2 mod 3` so it **flips the ⟨4⟩ coset** — it is
  NOT a pure roll in the bridge's base-4/R10 coordinate. A roll-scan cannot represent the coset-flip, so it can't see
  the identity even if it holds.
- **Convention:** `stationary_trunc` uses **a≥1** (v≥1), while the clean collapse needs **a≥0** (odd a=1 ↦ a′=0, which
  is excluded, leaving the ½(m₁)_*ν boundary term Wilson flagged). The a=0 term is genuinely absent here.

So the collapse is **not confirmed and not refuted** — it needs either the base-2 dlog coordinate (where ×2⁻¹ is a
plain shift) or the a≥0 iid measure. I hand this back rather than build a new coordinate transform (reconstruction
guardrail). What IS clean and re-usable: the parity sub-profiles ρ_e, ρ_o are computable (partial_pihat→bridge→ifft),
and `ρ_e+ρ_o = ρ_full` to the SINGLEREC truncation (~2^{−Amax}).

## Net
The load-bearing sentence is corrected from data: **A_j(m) is cross-parity iff 3∤m** (d₁ = m=1 = cross-parity stands;
the "odd-lag" reason is retracted). The parity split = the enriched/depleted dichotomy. Wilson's `ν_o=½(×2⁻¹)ν_e`
collapse — and hence the reduction of d₁ to a single-sub-measure autocorrelation at the ratio-(−2) lag — is untested by
this probe (needs base-2 coordinate / a≥0); it's the right next pen step, but the mechanism it rests on is `3∤m→cross`,
not odd/even. Not at stake: P6 finding (d₁=cross-term), P1LVL, BRIDGE2, CHANNEL_ID, MEAN1, dichotomy, R1–R30. Cheap.
