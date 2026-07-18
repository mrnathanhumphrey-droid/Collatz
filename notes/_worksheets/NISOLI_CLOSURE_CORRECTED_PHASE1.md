# NISOLI_CLOSURE_CORRECTED_PHASE1 — closure inequality articulated at λ = 43/45

**Date:** 2026-05-12. Wilson, Route A of the T_lead Nisoli-bypass re-evaluation. Phase 1 of the Nisoli closure-inequality test at the corrected rate.

---

## 1. The Nisoli closure inequality (verbatim from R77.2)

From `result_77_2_nisoli_certification.md` §3.1 (Nisoli 2026 Lemma 2.9 / Assumption 2.8):

> **η := ε_K · sup_{z∈γ} ‖R(z, T_K)‖ < 1**

where
- `ε_K := ‖L − L_K‖_op` is the truncation error (project-translation: `|K_truncated| / √q`, with `q = 3^K` the modulus).
- `γ` is a closed contour in ρ(T_K) isolating the target eigenvalue λ_2.
- `M := sup_{z∈γ} ‖R(z, T_K)‖` is the resolvent norm on γ.

For the project's parameterisation (used in the brief and in `M3_CLOSURE_TABLE.md`):

> **|K| · K^{−A} · M_3 < 1**

where
- `|K|` is the bilinear bound on `|S_partial|` at level r (the Plancherel-truncation level; `r = K − 1` in some conventions).
- `K^{−A}` is the Tao Prop 1.17 algebraic decay on `|μ̂_n(ξ)|` (the Syracuse Fourier coefficient).
- `M_3 = ‖(I − T)^{−1}‖` is the resolvent norm at z = 1 (the target eigenvalue's location).

At the corrected rate, the target eigenvalue is **λ = 43/45** (T_lead's leading eigenvalue), and the resolvent of interest is `(I − T_lead)^{-1}` evaluated at z = 1.

---

## 2. M_3'' verified

T_lead = (1/45) · [[7, 9], [28, 36]] (from `T_LEAD_CORRECTED_SPECTRUM.md`).

`I − T_lead = (1/45) · [[38, −9], [−28, 9]]`

`det(I − T_lead) = (38·9 − (−9)(−28))/45² = (342 − 252)/2025 = 90/2025 = 2/45.`

By Cramer:

> **(I − T_lead)⁻¹ = (1/2) · [[9, 9], [28, 38]] = [[4.5, 4.5], [14, 19]]**

Operator ℓ² norm via SVD: with A = (I − T_lead)⁻¹,

  A^T A = (1/4) · [[9² + 28², 9·9 + 28·38], [..., 9² + 38²]]
        = (1/4) · [[865, 1145], [1145, 1525]]

  trace(A^T A) = 597.5
  det(A^T A) = (865·1525 − 1145²)/16 = (1,319,125 − 1,311,025)/16 = 8100/16 = **506.25** (exact)

Eigenvalues: σ² = (597.5 ± √(597.5² − 4·506.25))/2 = (597.5 ± √354,981.25)/2 = (597.5 ± 595.803)/2.

  σ²_max = 596.6515
  σ_max = √596.6515 ≈ **24.426**

> **M_3'' = ‖(I − T_lead)⁻¹‖₂ ≈ 24.43** (exact via the closed-form computation; the slight gap above the spectral radius 45/2 = 22.5 is the κ(V) ≈ 1.086 condition factor of T_lead's non-orthogonal eigenbasis (1, 4) ⊥̸ (9, −7)).

This **matches** T_LEAD_CORRECTED_CLOSURE.md §2's ≈ 24.4 reading.

---

## 3. Bilinear bound `|K_bil|` at level r (PATH2 + HENSEL)

From `PATH2_DISPOSITION.md` (eq 190 closure family-level at q = 3) and `HENSEL_DISPOSITION.md` (Hensel-lifted closed form at r ≥ 4):

| r | N = 3^{r−1} | √N | |K_bil| bound | regime |
|---|---|---|---|---|
| 2 | 3 | 1.732 | 2·√3 ≈ 3.464 | strict 2·√N (PATH2 r ≤ 3) |
| 3 | 9 | 3 | **6** | strict 2·√N (PATH2 r ≤ 3) |
| 4 | 27 | 5.196 | 2·√3·√27 = 2·9 = **18** | polylog-free 2·√p·√N (HENSEL r ≥ 4) |
| 5 | 81 | 9 | 2·√3·9 ≈ **31.18** | polylog-free |
| 6 | 243 | 15.59 | 2·√3·15.59 ≈ **54.0** | polylog-free |
| 8 | 2187 | 46.77 | ≈ 162.0 | polylog-free |
| 10 | 19683 | 140.30 | ≈ 486.0 | polylog-free |

(p = 3 throughout; q = 3 in the c = 7/45 context.)

---

## 4. The A constant (Tao Prop 1.17 algebraic decay)

Tao Prop 1.17 (Tao 2022, p. 12, eq 1.25):

> For ξ ∈ Z/3ⁿZ not divisible by 3, `|E e^{−2πi ξ Syrac(Z/3ⁿZ)/3ⁿ}| ≪_A n^{−A}` for any fixed `A > 0`, with implied constant uniform in n, ξ.

Two facts about A:

**(a)** A is a **free parameter** in Tao's statement. It is not a derived quantity for our specific operator. Tao's proof works for any A > 0 — picking larger A in principle gives faster decay but with a worse `C_A` constant.

**(b)** The **effective `C_A`** is the load-bearing question. `BOOKKEEPING_PHASE1_DISPOSITION` rules: **C_A extraction from Tao §7.2-7.4 is INFEASIBLE.** The iterated-cubic recursion in Case 3 of Prop 7.8's proof forces super-exponential A-dependence regardless of bookkeeping quality. At K = 10, no value of A produces a `C_A` small enough to satisfy Nisoli η < 1 with the corrected M_3'' ≈ 24.

The closure inequality as parameterised — `|K_bil| · K^{−A} · M_3 < 1` — assumes `C_A = 1`. This is the **same parameterisation used in `M3_CLOSURE_TABLE.md`**, with the explicit caveat: "conditional on C_A = O(1)." Phase 1c of BOOKKEEPING rules out the C_A = O(1) regime.

---

## 5. Closure inequality at the corrected rate — the threshold

Substituting:

> `|K_bil(r)| · K^{−A} · 24.43 < 1`
>
> ⟺ `K^A > |K_bil(r)| · 24.43`

The thresholds at p = 3:

| r | |K_bil| | |K_bil| · M_3'' (threshold) |
|---|---|---|
| 2 | 3.464 | 84.6 |
| 3 | 6 | **146.6** |
| 4 | 18 | 439.7 |
| 5 | 31.18 | 761.7 |
| 6 | 54.0 | 1,319 |
| 8 | 162.0 | 3,957 |
| 10 | 486.0 | 11,873 |

For closure to fire at level r and Nisoli truncation parameter A:

> **K > (|K_bil(r)| · M_3'')^{1/A}**

This is the load-bearing inequality for Phase 2-3 tabulation.

---

## 6. Comparison to rate-1/2 case (R77.3 falsified)

| Quantity | rate-1/2 (R77.3 falsified) | corrected rate 43/45 |
|---|---|---|
| Target eigenvalue λ | 1/2 | 43/45 ≈ 0.956 |
| 1/(1−λ) (spectral radius of (I−T)⁻¹) | 2 | 45/2 = 22.5 |
| κ(V) condition factor | depends on T_3's eigenbasis (V is companion matrix, ‖V‖·‖V⁻¹‖ ≈ 1.84·1000 → 944 loose; ~50-200 tight) | T_lead's (1,4)⊥̸(9,−7), κ ≈ 1.086 |
| M_3 (operator norm) | R77.2 quoted 800-1000 loose; M3_DISPOSITION estimated 50-200 numerical | **24.43 exact** |
| |K_bil| at r=3 | 6 (strict) | 6 (strict) |
| Threshold |K_bil|·M_3 at r=3 | 300-1200 (loose); 300-600 (tight) | **146.6** |

**Counterintuitive finding:** even though λ = 43/45 sits much closer to 1 than λ = 1/2 does, the **exact M_3'' = 24.43 is SMALLER than the rate-1/2 case's M_3 in absolute terms** because:
- At rate-1/2, T_3 is a 3×3 non-normal companion matrix with κ(V) ≈ 50-1000.
- At rate 43/45, T_lead is a 2×2 rank-1 matrix with κ(V) ≈ 1.086 (eigenbasis is nearly orthogonal in Frobenius sense for this specific (1,4) vs (9,−7)).

So at the **same r = 3**, the corrected-rate threshold (146.6) is *easier* than the rate-1/2 threshold (300-1200) — but the entire question is structural-meaning, not numerical.

---

## 7. Closure inequality structure — summary

For the Nisoli closure test at the corrected rate λ = 43/45:

- **|K_bil(r)|**: explicit, rigorous from PATH2 (r ≤ 3) and HENSEL (r ≥ 4). At p = 3, |K_bil(r=3)| = 6, |K_bil(r=4)| = 18.
- **M_3'' = 24.43**: exact (closed-form 2×2 SVD), verified above and in T_LEAD_CORRECTED_CLOSURE.md.
- **K^{−A}**: this is Tao Prop 1.17's decay. The exponent A is a free parameter; the constant C_A is the load-bearing question.

The inequality fires iff **K^A > |K_bil(r)| · M_3''** at some plausible (A, K, r). Phase 2-3 tabulates this.

The **honest gate** (Phase 4): even if Phase 2-3 finds firing cells at e.g. A=2, K=20, r=3, **A is not extractable** in any quantitatively useful form from the project's machinery (BOOKKEEPING_PHASE1_DISPOSITION's INFEASIBLE finding). So the closure inequality is **structurally parameterised** but not **evaluated** against any specific A delivered by Tao.

---

## 8. Files

- T_LEAD_CORRECTED_CLOSURE.md (M_3'' ≈ 24.4 derivation, source)
- T_LEAD_CORRECTED_SPECTRUM.md (T_lead spectrum {43/45, 0})
- result_77_2_nisoli_certification.md (R77.2 Nisoli framework; original closure setup)
- BOOKKEEPING_PHASE1_DISPOSITION.md (Tao C_A INFEASIBLE)
- PATH2_DISPOSITION.md + HENSEL_DISPOSITION.md (bilinear bounds)
- M3_CLOSURE_TABLE.md (analog rate-1/2 tabulation for comparison)
- NISOLI_CLOSURE_CORRECTED_PHASE1.md (this file)
- nisoli_closure_corrected.py (main-thread verification script)
