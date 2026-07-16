# Result G0 (L3 Phase 0 — FREEZE THE OBJECT): the gate CAUGHT a definitional drift in R46. The uniform-lift refinement operator L_k does NOT carry the spectral gap — its gain is 1/√3 for EVERY q (gapless, first-moment). r_q is second-moment (build_M / R42), exactly where it always lived. R46's "L_k realizes r_q, σ_max<1/√3" is REFUTED. Phase 0 did its job.

**Date:** 2026-07-16. **Gate G0** of the L3 campaign (Phase 0). Fresh code of `L3_DEFINITIONS.md` §1–§7 (`probe_G0_phase0.py`), importing NOTHING from prior probes, computing BOTH the operator norm σ_max(L_k) (full SVD) and the stationary gain γ_k, at q=3,5,7 per k.

**Verdict: ✗ G0 FAILS its pass bar — and the failure is the point. R46 fused two different objects; G0 separated them.**

## The three facts (fresh implementation)

| q | REFINE `‖L_k d_k − d_{k+1}‖` | γ_k (= L_k's gain on d_k) | σ_max(L_k) (operator norm) | R46 claimed √(r_q/3) |
|---|---|---|---|---|
| 3 | 1.3e−16 … 1.1e−15 | 0.5684, 0.5790, 0.5782, 0.5778, 0.5769 → **≈1/√3** | 0.634, 0.665, 0.670, 0.6706, 0.6706 | 0.5774 |
| 5 | 2e−16 … 6e−16 | 0.5758, 0.5762, 0.5775 → **≈1/√3** | 0.6070, 0.6072, 0.6072 | 0.4546 |
| 7 | 3e−16 | 0.5777, 0.5768 → **≈1/√3** | 0.6456, 0.6459 | 0.3606 |

1. **REFINE holds** (1e-16): `L_k d_k = d_{k+1}` — the identity is real. ✓
2. **γ_k ≈ 1/√3 for EVERY q** (0.577 at q=3, 5, and 7). The uniform-lift refinement gain has **NO gap**. It does **not** equal √(r_q/3) — R46's q=5,7 values (0.455, 0.361) are wrong for this operator.
3. **σ_max(L_k) > 1/√3 for all q** (0.61–0.67), non-monotone — not r_q, not below the boundary.

## The drift G0 caught

R46 reported "σ_max(L_k) = √(rate_k/3)" where `rate_k = cB_{k+1}/cB_k`, `cB_k = 3^k(‖π_k‖² − (1/3)‖π_{k−1}‖²)`. But:

- **L_k propagates the `/q` deviation** `d_k = π_k − π_{k−1}/q` (§4 lift is the uniform fiber spread, weight 1/q). Its gain γ_k = ‖d_{k+1}‖/‖d_k‖ → **1/√3, for all q** (the first-moment deviation decays at exactly the Perron rate 1/3; **no gap**).
- **R46's `cB_k` is the `/3` Perron object** `‖π_k‖² − (1/3)‖π_{k−1}‖²`. The 1/3 = Σ_v p_v² is the **pair second-moment** eigenvalue, NOT a fiber lift. Its rate → r_q (0.62, 0.39) — the real gap.
- **For q=3, 1/q = 1/3, so the two objects coincide** (R45 and R46-at-q=3 are unaffected). **For q=5,7 they differ**, and R46 measured the gap from the `/3` object while labeling it the `/q` operator L_k's gain. Category error, now separated.

## Walk-back (visible, per campaign protocol: gate fails → phase stops → walk-back → STATE entry)

**REFUTED (R46 over-claim):**
- "L_k = P_W∘lift∘K realizes r_q as its gain / is a cleaner injection-free L3 target than build_M."
- `L3_DEFINITIONS.md` §9 as written ("σ_max(L_k) < 1/√3 for d≥3, = 1/√3 at d=2"): FALSE for this L_k (σ_max = 0.67 > 1/√3 even at q=3=d=2; γ_k = 1/√3 for all q, so no d≥3 vs d=2 distinction).

**SURVIVES:**
- The three identities FORGET / ONE-STEP / REFINE (G0 re-confirms REFINE at 1e-16). They are real structural facts about the **first-moment** refinement — just not where the gap is.
- The q=3 2/3-death (R45 / R46-at-q=3): at q=3, /q = /3, so those results stand.
- r_q itself (0.62, 0.39): unchanged, banked, living in build_M / R42 exactly as before.

**SHARPENED (the genuine gain from the catch):**
- **The spectral gap is intrinsically second-moment.** No single-address / first-moment refinement operator can carry it — G0 proves the single-vector L_k has gain 1/√3 for every q. This is the platform map (R44) made quantitative: first moment = Tao = density-1 = gapless; **second moment = Nathan = r_q = the gap.**
- **The L3 object is the pair / second-moment operator** (build_M's (a,b,γ), or R42's renewal A(z)) — where r_q has always lived. There is no simpler single-vector home. R46's hope for one was the drift; G0 closed it.

## Re-freeze required (Phase 0 not yet complete — Nathan to adjudicate)
`L3_DEFINITIONS.md` §7/§9 must be re-authored on the **second-moment object**. Candidates: (a) build_M pair operator, r_q = |λ₂|/λ₁ (banked, gate-validated); (b) R42 renewal, r_q = subdominant singularity of A(z) = Σ S_0(i) z^i. The clean single-vector L_k is retired as the target (kept only as the first-moment structural scaffold).

## Not at stake
R1–R45 (all unaffected). R46's identities stand; R46's "L_k carries r_q" claim is walked back. No r_q value changes — this corrects an operator LABEL, and relocates the L3 target back to the (correct) second-moment object.

_Reporting discipline: G0 was coded fresh from the written page precisely to catch page-vs-build drift, and it did — the σ_max/γ separation plus the /q-vs-/3 object distinction exposed an R46 conflation that the q=3-only R45 could not have surfaced (1/q=1/3 there). The pass bar failed honestly; the walk-back is banked before any Phase-1 work builds on the wrong object. This is the campaign's Phase-0 control working as designed._
