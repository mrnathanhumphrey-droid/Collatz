# Phase 1c — Projected C_A Looseness vs Nisoli Requirement

**Companion to:** Phase 1a (`TAO_PROOF_CONSTANT_MAP.md`), Phase 1b (`TAO_BOOKKEEPING_TRACTABILITY.md`).
**Goal:** assuming Phase 2 fires and successfully extracts each constant at a conservative upper bound, project what C_A would look like and check against Nisoli's η < 1 requirement at K = 6 and K = 10.

## 1. Setup — what number we need

From R77.2 §3.3 (`result_77_2_nisoli_certification.md`):

- M_3 := sup_γ ‖R(z, T_3)‖ ≈ 800–1000 (conservative; sharper rationals give ≈ 50–100 for ‖V⁻¹‖ alone but we keep M_3 ≈ 10³).
- γ = circle of radius 1/8 around λ_2 = 1/2.
- ℓ(γ) = π/4 ≈ 0.7854.
- Nisoli Lemma 2.9 requires η = ε_K · M_3 < 1 ⇒ **ε_K < 1/M_3 ≈ 10^{-3}**.

The translation from Prop 1.17 to the operator-norm bound ‖T − T_K‖ used in Nisoli:

```
‖T − T_K‖_op ≤ (Plancherel constant) · sup_ξ |μ̂_K(ξ)| ≤ Plancherel · C_A · K^{-A}
```

(The "Plancherel constant" depends on how the Syrac-distribution Fourier transform feeds into the rate operator; conservatively bounded by O(1) — say ≤ 4 — by the explicit T_3 construction in R77.2. **NOT** the squared form C_A² · K^{−2A}; the squaring in the R77.2 cert §3.4 was an artifact of an over-cautious reading. The cleanest bound is linear in C_A · K^{−A}.)

So the Nisoli closure requirement becomes:
```
4 · C_A · K^{-A} · M_3 < 1
⇔ C_A < K^A / (4 · M_3)
⇔ C_A < K^A / 4000  (with M_3 ≤ 10³)
```

At **K = 6**:
- A = 1: C_A < 6/4000 = 1.5 × 10^{-3}. **Vanishingly tight; practically infeasible.**
- A = 2: C_A < 36/4000 = 9 × 10^{-3}. **Still very tight.**
- A = 3: C_A < 216/4000 = 0.054. Tight but conceivable.
- A = 4: C_A < 1296/4000 = 0.324. Plausibly achievable.
- A = 6: C_A < 6^6/4000 = 46656/4000 ≈ 11.7. **Very achievable.**

At **K = 10**:
- A = 1: C_A < 10/4000 = 2.5 × 10^{-3}.
- A = 2: C_A < 100/4000 = 0.025.
- A = 3: C_A < 1000/4000 = 0.25.
- A = 4: C_A < 10000/4000 = 2.5. **Achievable for reasonable C_A.**
- A = 6: C_A < 10^6/4000 ≈ 250. **Very achievable.**

So the question reduces to: **for fixed A (large enough that K^A/4000 > what we'd get for C_A), what does the bookkeeping yield for C_A?**

## 2. Projected C_A as a function of A — conservative upper bound

We assemble C_{A,ε} ≡ C_A (since C-3 ε is absorbed into the absolute constants once we fix ε = 1/100 at the end) by walking through §7.4's three cases.

### 2.1 Case 1 contribution

The Case 1 threshold m₀^{(1)}(A, ε) ≈ exp(C₁ · A/ε) where C₁ is a small absolute constant — from C-20/C-21 of Phase 1b. Conservative C₁ ≤ 10 (very loose). With ε = 1/100, this gives:

```
m₀^{(1)} ≤ exp(10 · 100 · A) = exp(1000 · A).
```

This is the m below which Case 1 cannot fire and we must use the base case Q_m ≤ m^A. So the **effective C_A from Case 1 alone** is at most `(m₀^{(1)})^A = exp(1000 · A²)`.

**This is enormous.** At A = 1, this is exp(1000) ≈ 10^{434}. The looseness here is dominated by the ε^{-1} factor inside the threshold, which traces directly back to Tao's choice ε ∈ (0, 1/100) and to the white-point gain factor exp(−ε) (C-19). The (−ε) gain is *small* and so you need *many* white-point crossings to compete with the m^{−A} factor degradation — equivalently, m has to be very large.

### 2.2 Case 2 contribution

Case 2 threshold m₀^{(2)}(A, ε) is similar shape: requires `O(A/log m)` correction to be < ε/2, i.e., `m ≥ exp(C₂ · A/ε)` for absolute C₂. Conservative C₂ ≤ 20 (looser than Case 1 because we go through Lemma 7.7 / C-14).

```
m₀^{(2)} ≤ exp(2000 · A).
(m₀^{(2)})^A ≤ exp(2000 · A²).
```

### 2.3 Case 3 contribution

Case 3 is where the local-CLT (C-14) and large-deviations (C-28) compound. The threshold P = P(A, ε) must satisfy:
1. P > 10A/3 + 1 (from C-30).
2. The recursion `p′ ≤ p + 40·A·(1+p)³ + 10A/3 + 1` must iterate R = A²/ε times within P steps. Iterating `p_new = 40·A·(p+1)³ + O(A)` from p = 0 for R iterations, the bound grows hyper-exponentially: after one step `p ≈ 40A`, after two `p ≈ 40A · (40A)³ = 40^4 A^4`, after three `≈ 40^{1+3+9} A^{1+3+9} = 40^{13} A^{13}`, and after R = A²/ε ≈ 100 A² iterations the exponent has hit `3^{100 A²}` in the worst case before saturation.

**This is the dominant looseness source.** Even at A = 1, ε = 1/100, R = 100, the iterated cubic recursion gives P at least double-exponentially large in A.

Tao's exposition does not need to track this because Prop 1.17 only claims n^{−A} for fixed A. The implicit constant in `P = O_{A,ε}(1)` can be absorbed without effort. **In the bookkeeping pass, this iterated-cubic loosening is the headline cost.**

Conservative estimate:
```
P(A, 1/100) ≤ exp(exp(C₃ · A²))   with C₃ ≤ 5 absolute,
C_A ≤ (P(A, 1/100))^A ≤ exp(A · exp(5 A²)).
```

At A = 1: C_A ≤ exp(exp(5)) = exp(148) ≈ 10^{64}.
At A = 2: C_A ≤ exp(2 · exp(20)) ≈ exp(10^9) — astronomical.

### 2.4 Plus accumulated Vinogradov drift

The 18 named-unspecified constants in Phase 1a, under Tao's "vary line to line" convention, contribute a multiplicative drift factor of order 2^{18} ≈ 2.6 × 10^5 at worst-case interpretation, or O(1) if every ≪ is the same fixed constant. We use 10^5 as a conservative absorber.

### 2.5 Total projected C_A (conservative upper bound)

```
C_A^{conservative} ≤ 10^5 · exp(A · exp(5 A²))
```

(Case-3 dominates the Case-1, Case-2 thresholds.)

### 2.6 Tabulated against Nisoli requirement

| A | C_A^{cons} | Nisoli req @ K=6 | Nisoli req @ K=10 | Pass @ K=6? | Pass @ K=10? |
|---|---|---|---|---|---|
| 1 | 10^{64} | 1.5×10^{-3} | 2.5×10^{-3} | NO | NO |
| 2 | 10^{10^9} | 9×10^{-3} | 0.025 | NO | NO |
| 3 | 10^{10^{20}} | 0.054 | 0.25 | NO | NO |
| 4 | 10^{10^{34}} | 0.324 | 2.5 | NO | NO |
| 6 | 10^{10^{78}} | 11.7 | 250 | NO | NO |
| 10 | 10^{10^{217}} | ~10^7 | ~10^{10} | NO | NO |

**No value of A produces a C_A^{cons} small enough to satisfy Nisoli η < 1.** The growth of C_A in A massively outpaces the gain from K^A.

## 3. Sharpened (less conservative) bookkeeping — does it close the gap?

The C-3 (Case 3 iterated cubic) is the dominant source. **Even if Phase 2 produced a SHARP bookkeeping** — replacing the worst-case `exp(C₃·A²)` recursion bound with the tightest plausible `O(A^{20})` polynomial (this is a fantasy floor, not justified by anything other than wishing), we would still have:

```
C_A^{optimistic} ≈ A^{20·A}
```

| A | C_A^{opt} | Nisoli req @ K=10 | Pass @ K=10? |
|---|---|---|---|
| 1 | 1 | 2.5×10^{-3} | NO |
| 2 | 2^{40} ≈ 10^{12} | 0.025 | NO |
| 3 | 3^{60} ≈ 10^{29} | 0.25 | NO |
| 4 | 4^{80} ≈ 10^{48} | 2.5 | NO |
| 6 | 6^{120} ≈ 10^{93} | 250 | NO |

**Even under the optimistic floor, Nisoli closure fails at K = 10 for every A.**

## 4. Why? Structural diagnosis

The reason is fundamental to the proof shape:

- Tao's bound is `n^{−A}` for any fixed A. The implicit constant `C_A` is **never claimed to grow slowly in A**; in fact the proof structure (iterated cubic in Case 3) is consistent with C_A growing **at least exponentially in A²**.
- Nisoli's η < 1 requirement asks for `C_A · K^{−A} < 1/M_3`. Equivalently: `A log K > log C_A + log M_3`, i.e., **A must outpace log C_A divided by log K.** If log C_A grows like A² or worse, **no finite A solves the inequality** at any finite K.
- The only way out is `C_A` polynomial in A — which Tao's proof does NOT deliver.

**The proof is structurally inadequate for any effective C_A application at the project's verified K range.** This is not a bookkeeping defect; it is a feature of the technique. Tao's renewal-process argument optimizes for "any fixed A" qualitatively, not for any specific A with a small constant. To get a polynomial-in-A bound, you would need to **rewrite the proof using a different method** (e.g., direct Fourier-analytic L² bound, sum-product / Bourgain-Konyagin technique, or Heilbronn-type ζ-function arguments). That is **out of scope** for what was framed as a bookkeeping pass.

## 5. Sensitivity check — could M_3 shrink to make it work?

The R77.2 certification listed M_3 ≈ 800–1000 conservatively, with sharper rationals giving ≈ 50–100 for ‖V‖·‖V⁻¹‖. Even using M_3 = 50 (4× lower than the conservative 800 estimate):

- K = 10, A = 4: req C_A < 10^4 / (4·50) = 50. C_A^{opt} = 10^{48}. Still fails by 46 orders of magnitude.

**No tightening of M_3 within plausible bounds saves the calculation.** The C_A growth is the dominant obstruction.

## 6. Sensitivity check — could a different operator-norm route avoid the Plancherel-style ‖T − T_K‖ ≤ C_A · K^{−A} link entirely?

This is a different question from "is bookkeeping tractable." If the project found a way to bound ‖T − T_K‖_op directly (e.g., via the spectral / Riesz representation of T at finite N truncations without going through Prop 1.17), Nisoli could fire without using C_A at all. **That is the R77.3 route** (`result_77_3_nisoli_bypass.md`) and is the project's already-identified alternative path. Phase 1c does not address that route; it remains the most promising overall direction.

## 7. Conclusion of Phase 1c

> **Even if Phase 2 bookkeeping were performed perfectly and produced C_A in optimal closed form, the resulting bound would not satisfy Nisoli η < 1 at any K ≤ 10 for any A.** The looseness is structural to Tao's renewal-process proof technique, not to the bookkeeping. The Case-3 iterated-cubic recursion makes C_A grow at least like exp(exp(A²)) under faithful bookkeeping, and even an optimistic A^{O(A)} growth (not justified by anything in the proof) fails the Nisoli closure at K = 10.

**Implication for Phase 1 disposition.** Phase 1b's three positive tractability checks pass (no expertise blocker), but Phase 1c is unambiguously negative. This locks the disposition at **TRACTABLE_BUT_LOOSE** or **INFEASIBLE** depending on interpretation. See disposition file.

---

End of Phase 1c — Looseness Projection.
