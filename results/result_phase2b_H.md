# Result — Real q=3 Δ-channel operator: H_EXACT + H_CIRC CONFIRMED (family complete to L=3, 18/18), H1 refined, instrument catch (shift-invert→LU). The toy is the parity quotient.

**Date:** 2026-07-16. Nathan's Δ-channel derivation of the real q=3 operator + Session-Two verdict (S refuted, co-invariance mechanism, braid, Real-T1 posed). This file banks that verdict and my gate results (H_EXACT, H_CIRC via H2, H1 amplitudes). Probe `probes/probe_phase2b_H.py`, log `logs/probe_phase2b_H_log.txt`.

**Headline: VALIDATION PASSES (gap L=2=2.911e-3, L=3=9.958e-5; c₀ closed form exact; braid anchors reproduced) — the build is confirmed identical to Nathan's. H_EXACT ✅ (c₀=Σw² exact eigenvalue at L=1,2,3). H_CIRC ✅ AT EVERY LEVEL — the FULL circulant family is in the spectrum: 2/2 (L=1), 6/6 (L=2, 2e-15), 18/18 (L=3, LU pivots). ⚠️ The 18/18 required correcting a shift-invert artifact: iterative eigensolvers FAIL on the defective q=3 operator (R39 Jordan) — same lesson as ESPRIT/G0c; LU-pivot singularity is the reliable tool. H1: the amplitude bulk sits on c₀ + the dynamical partner (confirmed), but the zero-amplitude cluster is the ODD-k subfamily (not all k≠0), and it is mid-spectrum — NOT R32's near-λ₁ tower (that identification not confirmed).**

## The object (Nathan's Δ-channel derivation)
Real operator `M = build_M_gen(3, L, 2, [λ^δ])` — phase group `⟨2⟩ mod 3^L`, order `D = 2·3^{L−1}`; folded weights `w_δ = λ^δ/Σ` (`= 2^{−δ}/(1−2^{−D})` at λ=½). Δ = {(a,a,0)} (D states). From (a,a,0), move pair (δ_a,δ_b) gives T = a(2^{−δ_a}−2^{−δ_b}); gate from γ=0 is `δ_a ≡ δ_b mod d` (d=ord₃(2)=2). Three channels:
- **Equal (δ_a=δ_b), weight w_δ²:** T=0, carry stays 0, stays in Δ ⇒ a **circulant convolution on ℤ/D** with eigenvalues `c_k = Σ_δ w_δ² χ_k(δ)`, `c₀ = Σw²`.
- **Same-parity unequal:** gate passes, T≠0, carry jumps ⇒ leaks out of Δ (per-step leak ≈ 5/9−1/3 = 2/9).
- **Opposite parity:** gate fails from γ=0 ⇒ dies.
- **Identity making the toy inevitable:** total Δ-survival = (Σ_odd w)² + (Σ_even w)² = the toy's λ₁. **The toy is the parity shadow (quotient) of the real Δ-dynamics — not an analogy, a quotient.**

## Validation gate (must pass) — ✅ PASSED
| L | D | dim | gap | banked | c₀ = Σw² | closed `(1/3)(1−4^{−D})/(1−2^{−D})²` |
|---|---|---|---|---|---|---|
| 1 | 2 | 12 | 0.2222 | — | 0.555556 | 0.555556 (diff 0) |
| 2 | 6 | 324 | 2.911e-3 | 2.9e-3 ✓ | 0.343915 | 0.343915 (diff 6e-17) |
| 3 | 18 | 8748 | 9.958e-5 | 1.0e-4 ✓ | 0.333336 | 0.333336 (diff 0) |

L=1 family = {5/9, −1/3} (the toy's spectrum exactly). L=2 all six c_k are eigenvalues (max dist 1.94e-15). Braid anchors reproduced: L=2 dynamical **0.346827 > c₀=0.343915** (dynamical on top); L=3 **c₀=0.333336 > dynamical 0.333236** (kinematic on top). **Build ≡ Nathan's build.**

## H_EXACT — ✅ CONFIRMED
`Σw² = c₀ = (1/3)(1−4^{−D})/(1−2^{−D})²` is an **exact eigenvalue** of the full M at L=1,2,3 (machine-zero match to the closed form). The 6-digit "suspicious precision" was exact, as pre-registered. (Left-eigenvector concentration on the a=b sector: Nathan-verified in the Session-Two read — character on Δ to 3e-16, exactly zero on diagonal-with-carry, dressed off-diagonal.)

## H_CIRC — ✅ CONFIRMED AT EVERY LEVEL (family complete)
| L | members confirmed | tool | residual / pivot |
|---|---|---|---|
| 1 | 2/2 | dense eig | 2.7e-17 |
| 2 | 6/6 | dense eig | ≤1.9e-15 |
| 3 | **18/18** | **LU pivots** | minpiv 2e-16…1e-13 (controls c_k+0.05: 0.02–0.34) |

**⚠️ INSTRUMENT CATCH (durable):** the first H2 pass used shift-invert `eigs(M, sigma=c_k)` and reported only 2/18 "hits" (c₀ and c_{D/2}), the rest "missed" (residual 0.3). **That was an artifact** — the q=3 operator is *defective* (R39 exceptional point), and ARPACK shift-invert fails on it exactly as ESPRIT failed on the Jordan (Q6-B) and eig+inv failed at G0c. The robust **LU-pivot singularity test** `min|diag(U)| of LU(M−c_kI)` (with a c_k+0.05 non-singular control) confirms **all 18 members are exact eigenvalues.** Nathan's "the rest window-limited-**not refuted**" call was correct; the window limitation is now resolved. *Rule reinforced: on the defective q=3 operator use direct/LU methods, never iterative eigensolvers (ESPRIT, ARPACK shift-invert).*

## H1 — amplitude in the true-v₀ mass sequence (v₀ = indicator(1,1,0)); dense exact at L=2
Exact modal amplitudes `A_i` in `‖π_k‖² = 1ᵀMᵏv₀ = Σ_i A_i μ_i^k` (left/right eigenvectors, NOT ESPRIT — mass-blindness lesson respected by using the exact decomposition, not a fit).
- **Bulk amplitude on c₀ and the dynamical partner — ✅ as pre-registered.** |A(c₀)|=42.7, |A(dynamical 0.346827)|=44.0 (magnitudes inflated by the near-defect `l·r→0`, so unreliable in size but unambiguously the carriers).
- **Zero-amplitude cluster = the ODD-k circulant subfamily** (robust ~5e-16): c₁, c₃, c₅ carry EXACTLY zero. **Deviation from "all k≠0 zero":** the EVEN-k members c₂, c₄ carry small but nonzero amplitude (|A|=0.039 each). So the split is even-k (carry) vs odd-k (zero), not 0 vs k≠0.
- **R32 identification NOT confirmed (deviation).** The pre-reg hoped the zero-amplitude family = R32's near-λ₁ "mystery tower cluster." But the zero-amplitude odd-k members are **mid-spectrum** (moduli 0.286, 0.206), not near λ₁=1/3. So they are a genuine zero-amplitude subfamily but **not** R32's near-λ₁ cluster; that identification stands open.

## Session-Two verdict (Nathan; banked as reported)
- **S (rotation symmetry) REFUTED:** commutator 0.44/0.26, nowhere near zero — the carry's integer division breaks the rotation symmetry, so the exact family is **not symmetry-protected.** One line, one dead pretty hypothesis.
- **Mechanism = CO-INVARIANCE (not right-invariance, not symmetry):** each left eigenvector ℓ_k is an exact eigen-observable, `E[ℓ_k(next)|·] = c_k·ℓ_k(·)` — character on Δ, exactly zero on diagonal-with-carry, dressed off-diagonal (right eigenvectors spread, 92% off Δ). The zeros on diagonal-with-carry are a genuine **consistency condition** (a theorem to earn), not bookkeeping.
- **Braid, as a statement:** the coalescing pair = one **kinematic** member (c₀ = Σw², closed form, tower-clocked 1/3 + (2/3)2^{−D}) + one **dynamical** (non-family); the **Perron role swaps** between them (dynamical on top at L=2, kinematic at L=3). gap(L) = their separation, dominated by the dynamical member's approach to 1/3.
- **Toy placement:** L=1 family = {5/9, −1/3} = the toy's spectrum. **D1 is the ground floor of the circulant family;** the tower refines D = 2 → 6 → 18. One proof template should climb the whole thing.
- **Real-T1 (posed, Nathan's pen next session):** construct ℓ_k as the c_k-discounted Δ-return functional (χ_k on Δ, 0 on diagonal-with-carry, expected discounted character value at return off-diagonal); verify Mᵀℓ_k=c_kℓ_k by the gate algebra. Honest content: the diagonal-with-carry zeros impose a real consistency condition; the discounted-return sum's convergence needs the carry-excursion weights to beat c_k. Build L=1 (12 states) fully by hand, induct up.

## Standing requests (updated)
- **H1/H2 fired (this file).** H_EXACT ✅, H_CIRC ✅ 18/18. H1 refined (odd-k zero; R32 id open). *H2's shift-invert retired in favor of LU pivots for this operator.*
- **G — still standing:** L=4 partner via localization reduction (braid point 3 + rate-law point 4); the full 236196-state solve is out of reach — needs the reduced operator on the coalescing subspace.

## Bookkeeping flag (Nathan) — addressed
The banked **gap(L=1) = 0.889** (Q2 recon table, `result_phase2a_recon.md`) does NOT map onto this build's L=1 spectrum ({5/9, 1/3, 1/9}; top-pair gap = 5/9−1/3 = 2/9 = 0.2222). Different normalization in the earlier Q2 extraction at L=1 only. **L=2/L=3 validate exactly (gaps 2.9e-3 / 1.0e-4), so nothing downstream is at risk.** The Q2 L=1 entry is annotated in STATE as a superseded L=1 normalization (this Δ-channel build is canonical for L=1).

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, THEOREM D1 (the toy, now placed as the family's ground floor). No `r_q` value changes. The dynamical partner (the non-family coalescing member) is untouched = the remaining unknown for D3.

_Reporting discipline: validation passed before H1/H2 were read. H_CIRC's 18/18 required catching a shift-invert artifact (defective operator) and re-confirming with LU pivots + controls — reported as a caught artifact, not buried. H1's deviations (even-k carry small amplitude; zero-amplitude set is mid-spectrum not R32's near-λ₁ cluster) are reported AS deviations, not "consistent with." The S-refutation, co-invariance, braid, and Real-T1 are Nathan's Session-Two results, banked as reported (left-eigenvector structure Nathan-verified; I confirmed the spectral/amplitude facts). The gap(L=1)=0.889 discrepancy is annotated, not silently overwritten._
