# Result — PROBE F2: LEMMA E-FORM is GATE-CONFIRMED (closed-form compressed-chain entries, exact 32,886/32,886 at L=3); the LTE ladder is EXACT ({6,12} at L=3); the QSD braid PASSES; the coalescence 2×2 is LOWER-TRIANGULAR (c₀ protected, tower→kinematic the only coupling).

**Date:** 2026-07-16. Judge/instrument work: gate Wilson's closed-form entry formula for the frozen compressed chain, settle the LTE cascade, confirm the QSD braid, and dump the 2×2 effective-pair data for the coalescence derivation. No proof authored here, no rate fit. Direct/dense at q=3 (INSTRUMENT LAW). Probe `probes/probe_phase2b_F2.py`, log `logs/probe_phase2b_F2_log.txt`, dumps `outputs/eform_gate_q3_L2.tsv`, `outputs/lte_ladder_q3.tsv`, `outputs/effective_2x2_q3.tsv`.

**Headline: the closed-form entry formula `entry[(e,γ)→(e′,γ′)] = R(e′−e)·N(e′,γ,γ′)/D` reproduces EVERY nonzero entry of the frozen compressed chain to 2.8e-16 — 414/414 at L=2, 32,886/32,886 at L=3, support identical (every zero is exactly N=0, since R(s)>0 always). LEMMA E-FORM is a gate-confirmed identity. The LTE ladder `v₃(2^{e′}−1)` is exact — the valuation-2 targets at L=3 are EXACTLY {6,12} as pre-registered, and from γ=0 the carry valuation into target e′ is deterministically t(e′)−1. The QSD braid PASSES (partner above c₀ at L=2, below at L=3 — matching the true partner). The 2×2 effective-pair matrix is LOWER-TRIANGULAR: c₀'s kinematic mode never feeds the tower (B[kin,tow]=0 to machine precision = Real-T1 protection made 2×2-visible), the only coupling is tower→kinematic, and the detuning changes sign (the braid) while the eigenvector-coalescence ratio |coupling|/|detuning| grows 17→189 (the EP approach).**

## F2-1 — LEMMA E-FORM: GATE-CONFIRMED (the judge)
Claimed identity (exact for E1's frozen source-side uniform average):
```
entry[(e,γ)→(e′,γ′)] = R(e′−e) · N(e′,γ,γ′) / D
   R(s) = Σ_δ w_δ w_{δ−s}          (circular autocorrelation = Real-T1's R_k(0) object, k=0)
   N(e′,γ,γ′) = #{units u : (γ + u·(1−2^{e′})) ≡ 0 mod 3  AND  ((γ+u·(1−2^{e′}))//3) mod 3^L == γ′}
   D = 2·3^{L−1}
```
| L | dim | nonzeros | max\|entry_actual − entry_formula\| | zero pattern | min R(s) |
|---|---|---|---|---|---|
| 2 | 54 | 414 | **2.78e-16** ✅ | identical ✅ | 0.0847 (>0) |
| 3 | 486 | 32,886 | **2.78e-16** ✅ | identical ✅ | 1.30e-3 (>0) |

- **Every entry matches to machine precision, both L. Pre-registered exact (≤1e-12) — PASSED with 4 orders of margin.**
- **Zero pattern proven structural:** min R(s) > 0 at both L ⟹ R(s) > 0 for all s (all weights positive) ⟹ `entry = 0 ⟺ N = 0` exactly. The support is carried entirely by the carry-count N, never by the move algebra R.
- **The factorization mechanism (why it holds), confirmed:** T = a·2^{−δa}·(1−2^{e′}) = u·(1−2^{e′}) with u = a·2^{−δa} a unit; as `a` ranges over the source class, u ranges over ALL D units, so the per-move carry-sum `Σ_a[gate·carry=γ′] = N(e′,γ,γ′)` is δa-INDEPENDENT. Summing the free (δa,δb) with δa−δb≡s gives the circular autocorrelation R(s). **Move algebra ⟹ R; dynamics ⟹ N; they factor.** Hand-checkable on row (0,0): self-loop `= R(0)·6/6 = Σw² = 65/189` (T=0, N=6); `(0,0)→(2,1) = R(2)·3/6` (T=6u mod 9 ∈ {3,6}, all pass, carry splits 3/3). Full per-entry audit with (s, R(s), N) in `outputs/eform_gate_q3_L2.tsv` (414 rows).
- **Cross-link:** the R in E-FORM is the SAME autocorrelation object as Real-T1's R_k at k=0 (`R(s)=Σ_δ w_δ w_{δ−s}`). E-FORM's move-algebra factor is literally the Real-T1 eigenvalue machinery restricted to the compressed chain — the kinematic and dynamical halves share the one autocorrelation.

## F2-2 — the LTE cascade (retro-explains C2's class-dependence)
Since T = u·(1−2^{e′}) with u a unit, **v₃(T) = v₃(2^{e′}−1) = 1 + v₃(e′/2)** for even e′ (LTE). Odd e′ give a unit T ⟹ die from γ=0.
- **(a) Exact integer ladder (L=3):** valuation-2 targets **EXACTLY {6,12}**; valuation-1 even = {2,4,8,10,14,16}; valuation-0 odd = {1,3,5,7,9,11,13,15,17}. **Pre-registered {6,12} — CONFIRMED.** (9 | 2^{e′}−1 ⟺ ord₉(2)=6 divides e′ ⟺ e′∈{6,12} in range; 27 | · ⟺ 18|e′ ⟺ only e′=0 ⟹ max finite valuation is 2 at L=3.) Full ladder `outputs/lte_ladder_q3.tsv`.
- **(b) Deterministic γ=0 carry law (theorem-shaped):** from the Δ-adjacent γ=0 sector, a move into even target e′ produces a carry of valuation **exactly t(e′)−1**, uniformly over all units u (v₃(T)=t(e′) for every unit ⟹ v₃(carry)=t−1). Verified deterministic at L=2 and L=3 for every even target (e.g. e′=6 → carry v₃=1). The "deep jumps" (carry v₃=1) are produced by exactly the t=2 targets {6,12}.
- **(c) C2 reproduction + attribution:** C2's per-(θ=e mod 3, γ) cascade deviation reproduced — max **0.317 at L=2 / 0.106 at L=3** (matches C2's reported 0.32 / 0.11). The target-side generator of the carry valuations is the LTE ladder; it pins WHERE the deep carries come from ({6,12}) and is exact/deterministic in the γ=0 sector. **Honest scope:** the full per-source-class deviation is a convolution of the source-γ valuation with this target ladder; F2-2 pins the target factor exactly (ladder + γ=0 determinism), not every per-class number. "Class-dependent cascade" is now target-side = the LTE ladder of the phase shifts.

## F2-3 — QSD admissibility, condition two (the braid) — PASS
The committed novel pre-registration: the QSD-compressed partner is ABOVE c₀ at L=2 and BELOW at L=3 (the crossing the uniform average provably missed — uniform is below-then-above, the wrong braid).
| L | c₀ | QSD partner (side) | true partner (side) |
|---|---|---|---|
| 2 | 0.343915 | 0.346827 (**above**) | 0.346827 (above) |
| 3 | 0.333336 | 0.333231 (**below**) | 0.333236 (below) |

- **BRAID CONFIRMED: above@L2, below@L3 — and it matches the TRUE partner's braid at both L.** ⟹ **QSD carries proof weight** (it reproduces the crossing uniform missed, condition two of admissibility met). Uniform (E's baseline) gives the wrong braid (below@L2, above@L3); QSD gives the right one.

## F2-4 — the 2×2 effective-pair (the coalescence-derivation data)
Spectral projection onto the 2-dim invariant subspace span{c₀-mode, partner} of the FULL operator (direct dense eig + EP-robust orientation), basis (kinematic = γ=0-aligned, tower). Dump `outputs/effective_2x2_q3.tsv`.

| L | c₀ = B[kin,kin] | partner = B[tow,tow] | B[kin,tow] | B[tow,kin] (coupling) | detuning c₀−partner | discriminant (l₀−lₚ)² |
|---|---|---|---|---|---|---|
| 2 | 0.34391534 (=65/189) | 0.34682666 | **0** (≤5e-9) | 0.05053138 | −2.911e-3 (partner **above**) | 8.476e-6 |
| 3 | 0.33333588 | 0.33323630 | **0** (≤5e-9) | 0.01882012 | +9.958e-5 (partner **below**) | 9.915e-9 |

- **The 2×2 is LOWER-TRIANGULAR.** `B[kin,tow] = 0` to machine precision at both L: **the kinematic (c₀, γ=0) mode never feeds the tower.** This is Real-T1's protection lemma made 2×2-visible — c₀ has an exact LEFT eigenvector (ℓ₀ on γ=0), so the restriction is triangular; the eigenvalues are exactly the diagonal (c₀ and partner, recovered to 8 digits). Consistent with P's left/right split (c₀ left localizes on γ=0; partner right lives in the tower, γ=0-weight 0.000).
- **The only coupling is tower→kinematic** (B[tow,kin] = 0.0505 → 0.0188). The partner is a dressed tower mode that leaks UP into c₀'s sector; c₀ does not leak DOWN. This is the exact asymmetry P measured.
- **The detuning changes sign** (−2.9e-3 → +9.96e-5) = the braid, now read directly off the diagonal. **Discriminant → 0** monotone (8.48e-6 → 9.92e-9).
- **EP-approach diagnostic (observation, two points, NOT a fit):** the eigenvector-coalescence ratio |coupling|/|detuning| = 0.0505/0.00291 = **17.4** (L=2) → 0.0188/0.0000996 = **189** (L=3). The detuning shrinks faster than the coupling ⟹ the two eigenvectors become more parallel ⟹ the finite-L operator approaches a defective Jordan block (R39) at the EP. **This is the coalescence: both diagonal entries → 1/3 while the ratio → ∞.** The 2×2 values are the JUDGE for Wilson's L-trend derivation (from E-FORM + the LTE ladder); the trend is NOT fitted here.
- **L=4: SIZED + DEFERRED.** The QSD-weighted L=4 2×2 needs μ = |dominant right eigenvector| of the 236,196-state operator (nnz ~2.3e8, ~5.5GB CSR). μ is the LARGEST mode ⟹ power iteration (pure sparse matvec, GPU-friendly, NOT the banned interior-near-EP LU), but the near-degenerate top (gap < 1e-4) makes subspace separation slow. Heavy compute ⟹ Lambda/greenlight. **Confirmatory only** (rate-fit banned); L=2,3 are the load-bearing points.

## Adjudication
| item | verdict |
|---|---|
| **F2-1 LEMMA E-FORM** | **GATE-CONFIRMED** — every entry exact to 2.8e-16, both L; zero pattern structural (entry=0 ⟺ N=0). The compressed chain's entries are the closed form R·N/D. |
| F2-2 LTE ladder | CONFIRMED exact — {6,12} at L=3 (pre-reg); γ=0 carry valuation = t(e′)−1 deterministic; C2 deviations reproduced (0.317/0.106), target-side = LTE. |
| F2-3 QSD braid | PASS — above@L2/below@L3, matches true partner ⟹ QSD admissible (proof weight). |
| F2-4 2×2 | LOWER-TRIANGULAR (c₀ protected, B[kin,tow]=0); tower→kinematic the only coupling; detuning sign-flips (braid); coalescence ratio 17→189. Values dumped as Wilson's judge; no fit. |

**⟹ The compressed-chain program now has a gate-confirmed closed form for its entries (E-FORM), an exact arithmetic ladder for the carry structure (LTE), a proof-weight braid (QSD), and a triangular 2×2 whose off-diagonal is the single tower→kinematic coupling — the exact object the coalescence derivation acts on. The kinematic side (c₀ = exact eigenvalue, protected) and the dynamical partner (dressed tower mode) are cleanly separated in the effective matrix; only the L-trend of the detuning/coupling remains for Wilson's pen.**

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, H, J, L-A, L-B, Real-T1, R1-STOP, Probe P/C/E/W. No `r_q` value changes; **no rate-law fit** (the true 2.9e-3, 1.0e-4 sequence untouched; the 2×2 detuning agrees with it but is NOT fitted to it — the discriminant L-trend is Wilson's to derive, the dumped values judge).

_Reporting discipline: E-FORM's pass is reported with the exact residual (2.8e-16) AND the structural zero-pattern argument (R>0 ⟹ zero ⟺ N=0), not just "matches." The LTE {6,12} prediction is reported as CONFIRMED with the exact ord₉/ord₂₇ reasoning. C2's deviations are reproduced and the attribution is scoped honestly (target-ladder pinned; source-γ factor acknowledged). The 2×2 triangularity is reported as machine-precision-zero (≤5e-9), tied to Real-T1's exact ℓ₀, not asserted as algebraic-zero without the numeric caveat. The EP-approach ratio is flagged as a two-point observation, not a fit. L=4 is sized and deferred (heavy compute, confirmatory). c₀-masquerade criterion applied at every extraction._
