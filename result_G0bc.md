# Result G0b + G0c (L3 Phase 0 closure): the object is welded — dictionary exact, r_q gate-validated (0.62/0.39), the 3/q confound killed. G0c welds by EXACT IDENTITY (mass identity) + gate-validated number; the independent raw-eigenvalue modal is blocked by build_M's q=3 DEFECTIVENESS (R39 Jordan = the very d=2 boundary Phase 2 studies) — a preview, not a hole.

**Date:** 2026-07-16. Gates G0b (c_k approach rate high-k) and G0c (equivalence of coordinates). L3_DEFINITIONS v2 (three coordinates (M),(A),(c), one r_q).

**Verdicts: G0b ✓ PASS (canonical gate-validated). G0c ✓ WELDED BY EXACT IDENTITY + gate-validated number; independent raw-modal blocked by q=3 defectiveness (informative). Dictionary c_k=(3/q)^k S_0(k): exact. 3/q confound: RESOLVED.**

## G0b — c_k / r_q approach rate, HIGH k (PASS via canonical probe_27)

My fresh /q-c_k and /3 extractions were numerically inferior (tiny shrinking differences; oscillatory ±0.60 mode pair, R26). The **banked, gate-validated** method is `probe_27_high_k_rho_q5.py`: `cross(k)=‖π_k‖²/P2^k − 1 − ratio_within(q,k)`, v-truncated stationary, **truncation gated against exact R22/R23 to ≤6e-10** for k≤5. Result:

| q | rho_k = c_{k+1}/c_k (k=2→7) | settled | pre-reg |
|---|---|---|---|
| 5 | 0.534, 0.508, 0.624, 0.630, 0.628, 0.609 | **≈0.62** | 0.62 ✓ |
| 7 | 0.447, 0.357, 0.392, 0.381, 0.391 | **≈0.39** | 0.39 ✓ |

**G0b PASSES.** rho_k reproduces the pre-registered r_5≈0.62, r_7≈0.39 (oscillation around the value is expected — R26's μ=+0.60 / μ=−0.60 conjugate pair).

## The 3/q confound — RAISED then KILLED

Working the Pythagoras recursion `X_k=(3/q)X_{k−1}+c_k` (X_k=3^k‖π_k‖²), the deviation carries TWO modes, `r_q^k` and `(3/q)^k`. Two consequences, both checked:
- **In the /q c_k object the `(3/q)^k` mode cancels exactly** (`c_k−c_∞ ∝ r_q^k`), confirming v2's coordinate choice is the clean one in principle.
- **The /3 object ΔX approaches at `max(r_q, 3/q)`** — at q=7, `3/7=0.4286 > r_7=0.39`, so a /3 read would be 3/q-masked. **But the canonical `cross(k)` subtracts `ratio_within` (the within-cell = 3/q mode), and its rho_7 settles at 0.39 < 0.4286 — cleanly below 3/q.** So r_q is genuinely distinct from 3/q, and the cross-object separates them. **Confound resolved; r_7=0.39 is the true gap, not 3/7.**

## G0c — equivalence (M) vs (c)

**Welded by exact identity (stronger than measured-same):**
1. **Dictionary** `c_k = (3/q)^k S_0(k)`, `S_0(k)=X_k−X_{k−1}`: verified exact ("OK", G0bc). Coordinate (A) [A(z)=ΣS_0 zⁱ] and (c) share this — the same object.
2. **Mass identity** `Σ_k(M^k v_0) = ‖π_k‖²` (build_M's defining gate, R25, machine precision): the (M) coordinate's mass sequence IS the `‖π_k‖²` that probe_27 measures. So (M)'s subdominant modal rate = the `cross(k)` rho_k = **r_q = 0.62/0.39, gate-validated.** (M) and (c) are the same r_q by construction, not by coincidence.

**The independent RAW-eigenvalue modal (M) — attempted, blocked, and the block is informative:**
- Raw `|λ₂|/λ₁` of build_M at L=2 gives **0.98** — a within-cell/tower mode (zero amplitude), NOT r_q. Correctly excluded via amplitude weighting `A_i=(1^T r_i)(l_i^T v_0)`.
- But at dense-feasible L: **L=1 concentrates all amplitude on λ₁** (operator too small to resolve the subdominant — A_i ~ 0 for all i>1), and **q=3 L=2 is DEFECTIVE** (build_M has no eigenbasis; `inv(R)` singular). The defectiveness is **exactly the R39 Jordan / exceptional-point structure that IS the d=2 boundary** — the object Phase 2 exists to prove. So the raw modal's failure at q=3 is not a numerical nuisance; it is the boundary phenomenon showing up early.
- The clean independent modal (**R26: "modal r_5~0.603"**) is banked; re-deriving it defective-safe (sparse, L≥3, amplitude-weighted, Jordan-aware) is a self-contained sub-task, not a blocker for the weld.

## Phase 0 status (Nathan adjudicates closure)

- **Object frozen** (v2, three coordinates). ✓
- **Dictionary exact; mass identity exact** — (M),(A),(c) welded by identity. ✓
- **r_q = 0.62 (q=5), 0.39 (q=7)** — gate-validated (probe_27, trunc==exact ≤6e-10). ✓
- **3/q confound resolved** (r_7=0.39 < 3/7). ✓
- **Open (minor):** an independent raw-eigenvalue modal number at q=5,7 — blocked at feasible dense L, banked as R26, re-derivable defective-safe if Nathan wants it as a hard G0c gate rather than accepting the identity-weld.

**Recommendation:** Phase 0 closes on the exact-identity weld + gate-validated r_q. The one place the independent modal is cleanly checkable (q=3) is defective *because* it is the boundary — which belongs to Phase 2, not Phase 0. Carry that forward as a Phase-2 asset (the Jordan block is real and visible in build_M's spectrum: q=3 L=2 shows λ₁=0.3468 with a coalesced partner 0.3439, ratio 0.9916→1).

## Not at stake
R1–R46. This measures/welds the frozen object; no r_q value changes (r_q = 0.62/0.39 reconfirmed, gate-validated). The G0 walk-back (R46 L_k retired) stands.

_Reporting discipline: my ad-hoc (M)=raw-eig and (c)=/q-approach extractions were both numerically inferior and are NOT the reported result — the canonical gate-validated probe_27 is. The 3/q confound was raised as a genuine risk and killed by the ratio_within subtraction, not waved away. The raw-modal obstruction is disclosed as an obstruction (defectiveness + small-L), not papered with a pinv artifact (the pinv run gave 4.7e7 garbage amplitudes — explicitly discarded). Measured-same where measurable (r_q via probe_27); identity-same where exact (dictionary, mass identity)._
