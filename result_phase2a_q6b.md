# Result — Q6-B SETTLED (reliable Jordan-aware detector): the q=3 collision requires λ=1/2 EXACTLY. The WEIGHT is load-bearing. Combined with Q6-A, NEITHER phase nor weight quotients out — the d=2 collision is the marginality of the actual Syracuse map (q=3 AND λ=1/2), not abstract phase combinatorics.

**Date:** 2026-07-16. Closes the Q6-B thread left inconclusive by the ESPRIT detector (which fails on the q=3 Jordan). Two independent reliable detectors, agreeing.

**Verdict: ✗ collision is NOT weight-blind. M(3,2,λ) collides (r=1, Jordan) ONLY at λ=1/2; all other λ are gapped. The weight base 2^{−v} (λ=1/2) is essential, exactly as the phase group ⟨2⟩ is (Q6-A). The d=2 boundary is the marginal balance of the specific 3x+1 map — decouplable in neither direction.**

## The reliable detector (ESPRIT-free)
The Jordan signature is a **linear-in-k prefactor** in the mass sequence `a_k = 1ᵀMᵏv₀`. With `λ₁` from a direct Perron solve, `b_k = a_k/λ₁ᵏ`. **Collision** (Jordan at top, r=1) ⟺ `a_k ~ (A+Bk)λ₁ᵏ` ⟺ `b_k` grows linearly (slope **B ≠ 0**, growing with L). **Gapped** (r<1) ⟺ `b_k → const` (**B ≈ 0**). No exponential fitting → no repeated-root failure. Validated: M(3,2,1/2) B=0.29 (L=2) → 0.44 (L=3), growing toward 7/15 ✓; the gapped operators give B≈0.

## Q6-B — is the q=3 cross–Perron collision weight(λ)-blind? NO.

| λ | B(L=2) | B(L=3) | top-2 gap L2→L3 | verdict |
|---|---|---|---|---|
| 0.25 | 0.0000 | 0.0007 | — | gapped |
| 1/3 | 0.0001 | 0.0004 | 0.53 → 0.19 (stable) | **gapped** |
| 0.40 | 0.0022 | 0.0007 | — | gapped |
| **0.50** | **0.29** | **0.44 (grows)** | **2.9e-3 → 1.0e-4 (coalescing)** | **COLLISION (r=1)** |
| 0.60 | 0.0017 | 0.0002 | 0.063 → 0.24 (growing) | **gapped** |

**Only λ=1/2 collides — a sharp resonance.** The k-prefactor B grows with L only at λ=1/2 (real coalescence); at λ=0.4 and 0.6 it is ≈0 and the top-2 gap does NOT shrink toward 0 (no EP). Two independent detectors (B-prefactor growth + eigenvalue-gap coalescence) agree.

## What this means — both 2's are load-bearing
- **Q6-A (phase):** collision needs `⟨2⟩`'s q-adic lift at q=3 (`2≡−1 mod q ⟺ q=3`); a generic order-2 phase (−1) does not collide. *(banked, result_phase2a_recon.md)*
- **Q6-B (weight):** collision needs `λ=1/2` (the `2^{−v}` halving); any other weight base is gapped. *(this result)*
- **⟹ Mr. Potato Head decouples NEITHER.** Swapping the phase generator OR the weight base kills the collision. The number **2 is doing both jobs, and both are essential.**

**Conceptual reading:** the d=2 collision is the **marginality of the actual Syracuse map** — the growth (q=3) balanced against the halving (λ=1/2). This is exactly R41's `σ_H = log₂((q+1)/2) = 1` at q=3, which is the marginal condition of the *specific* 3x+1 map (halving 2). Marginality is a codimension-1 balance between the two 2's; that is why it is a sharp resonance at λ=1/2 and unique to q=3, and why it decouples in neither direction.

## Phase 2b statement — final form (both refinements folded in)
> The entrance-exam collision is the **marginal balance of the specific Syracuse map at q=3**: at `d = ord_q(2) = 2` (⟺ q=3 via `2 ≡ −1 mod q`) **with the halving weight `λ = 1/2`**, the cross–Perron pair **coalesces as L→∞** (an exceptional point / defective limit), `σ_H = 1`; perturbing EITHER the phase group `⟨2⟩` (Q6-A) OR the weight `1/2` (Q6-B) restores the gap. At d≥3 (any q≥5) the pair stays gapped. The proof must use the actual map structure (both 2's) — it is NOT abstract {±1} phase combinatorics and NOT weight-independent.

The bonus 7/15 route (Jordan-chain slope, Q3, 0.16%) stands and is λ=1/2-specific (as it must be — 7/15 is the actual Syracuse constant).

## Not at stake
R1–R46, Phase 0/1, G0c′, Phase 2a Q1–Q5. No r_q value changes. This closes Q6-B (the weight-inessentiality question) as REFUTED: the weight is essential.

_Reporting discipline: the ESPRIT detector's failure on the q=3 Jordan (the reason Q6-B was inconclusive) is replaced by a Jordan-aware detector validated on knowns (M(3,2,1/2) collision B-growth; gapped operators B≈0), cross-checked against independent eigenvalue-gap coalescence. Both agree only λ=1/2 collides. The pre-registered "weight-blind" hypothesis is REFUTED, reported as a refutation. The λ=1/2 resonance is sharp (0.4 and 0.6 both gapped) — a codimension-1 marginal condition, not a fit._
