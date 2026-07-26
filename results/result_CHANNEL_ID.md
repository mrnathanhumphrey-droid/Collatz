# PEN RESULT + Probe CHANNEL_ID — the dominant mode IS the m=1 collision channel — **EXACT IDENTITY (derived by pen, verified exact-rationally): `d1_r = A_r(1)/A_r(0)` for r≥2 — the MODES frontier and the charledger are THE SAME OBJECT. Corollary (same proof): `Re δ̂_r(n) = A_r(n)/A_r(0)` for every 3∤n — the whole mode ladder is the channel ledger normalized by the m=0 increment (the 7/15-slope object). Equivalences: d1>0 ⟺ A_r(1)>0 (γ_r(1) strictly increasing = MOON's measured m=1 monotone+) ⟺ q_r(1) = Pr[X'≡4X mod 3^{r+1} | X'≡4X mod 3^r] > 1/3 (the conditional lifting of the ×4-collision beats the neutral cost) ⟺ adjacent dlog positions have positively correlated new-digit biases. CONSEQUENCE: the S2 splitting run's survival fractions q_s ≈ 0.3335 (m=1 channel, plateau r=17..38) ARE this object — S2 was unknowingly measuring sign(d1_r) statistically out to r≈38, and it reads POSITIVE throughout. Exact side: q_16 = 0.333751 > 1/3, seamlessly continuing into S2's ~0.3335. Three previously separate banked measurements (MODES d1 ladder, MOON m=1 monotonicity, S2 survival plateau) are one measurement.**

**Date:** 2026-07-25. Pen derivation (Claude, at Wilson's request: "try with your pen"). Verification `probes/probe_channel_id.py`, log `logs/channel_id_run.log` + `channel_id_v5.log`.

## The derivation (3 steps, all exact)
Setup: `C^{(r)}(k) = Pr[X'≡4^k X mod 3^{r+1}] = p_r(k)` (X,X' iid ~ν_r; the FKG ratio-distribution reformulation, which survives the shelf's death). Target lags {1, m−1, m+1}, m=3^{r−1}.
1. **The three lags are ONE lag in the three fiber channels:** `m−1 ≡ 1+2m (mod 3m)` and C even ⟹ target lags = `1+δm, δ=0,1,2`.
2. **The fiber sum telescopes to the level below:** for fixed s, `{s+1+δm}` is exactly the full fiber over (s+1 mod m) — no carries (the 3 points differ by m, 3m≡0) — and the tower gives `Σ_j ρ_r(t+jm) = ρ_{r−1}(t)` **exactly** (ν_r mod 3^r = ν_{r−1}, induction on the renewal; dlog commutes with reduction). Hence `Σ_δ C^{(r)}(k+δm) = C^{(r−1)}(k)` for any k.
3. **Collapse:** `Num = 2C(1)−C(m−1)−C(m+1) = 3p_r(1)−p_{r−1}(1)`; `Den = 2[C(0)−C(m)] = 3p_r(0)−p_{r−1}(0)`. Multiply by 3^{r−1} and recognize the banked charledger objects `γ_r(m')=3^r p_r(m')`, `A_r(m')=γ_r−γ_{r−1}`:
$$\boxed{\ d_1^{(r)} = \frac{A_r(1)}{A_r(0)}\ }\quad(r\ge2);\qquad \operatorname{Re}\hat\delta_r(n)=\frac{A_r(n)}{A_r(0)}\ \ (3\nmid n,\ \text{same proof}).$$
(r=1 degenerate: m=1, lag m−1=0. A_r(0)=X_r−X_{r−1} → S_∞ — the established linear-divergence increment.)

## Verification (all falsifiers passed)
- **V1 EXACT (Fractions, r=2..5):** tower fold `ρ_r → ρ_{r−1}` EXACT equality (True all r); identity == 5-lag formula **as exact Fractions** (True all r); `d1_2 = 2/35` exactly.
- **V2 float r=2..16:** identity vs banked d1 ladder, worst rel **2.0e−13**.
- **V3 corollary n=2:** vs banked d2 (r=12..16), rel ≤2.3e−7 (= banked print precision).
- **V4 the Λ weld:** `Σ_{m=1..13} 4^{−m}A_r(m)` reproduces the banked ν-route `Λ_r` (r=12..16) to ≤1.7e−5 (= 5-digit banked precision) — **the identity is welded to the ε ladder.**
- **V5 the q-ladder (the payoff table):**

| r | γ_r(1) | A_r(1) | q_r(1) | 3q−1 | A_r(0) |
|---|---|---|---|---|---|
| 1 | 0.66667 | −3.33e−1 | 0.2222 | −3.3e−1 | 0.667 |
| 2 | 0.69388 | +2.72e−2 | 0.34694 | +4.1e−2 | 0.476 |
| 4 | 0.70704 | +4.01e−3 | 0.33524 | +5.7e−3 | 0.464 |
| 8 | 0.71926 | +2.17e−3 | 0.33434 | +3.0e−3 | 0.466 |
| 12 | 0.72566 | +1.39e−3 | 0.33397 | +1.9e−3 | 0.469 |
| 16 | 0.73001 | +9.14e−4 | **0.333751** | +1.25e−3 | 0.471 |

**The seed (hand-exact):** ν₁ = (0, 1/3, 2/3) in dlog (X≡7 mod 9 iff v odd, prob 2/3) ⟹ γ₀(1)=1, γ₁(1)=2/3 — **the channel FALLS first** (A₁(1)=−1/3), then grows monotonically forever after (A_r(1)>0 all r=2..16). d1>0 for r≥2 is a recovery-and-growth statement.

## The unification (the big consequence)
`q_r(1)` is **exactly the m=1 per-level conditional survival the S2 splitting run measured** (`q_s ≈ 0.3335`, plateau r=17..38, N=4×10⁵×20 replicas). So:
- **S2 was a direct statistical measurement of sign(d1_r) for r=17..38 — and it reads positive (q>1/3) throughout.** The exact ladder (q₁₆=0.333751) continues seamlessly into the statistical plateau (~0.3335). The dominant mode's positivity is now *observed* to r≈38, not 16. (Statistical caveat: per-level 3q−1 ≈ 1e−3 is at the per-level noise scale; the *cumulative* rise was the z≈5 statement. Consistent-with-positive, not per-level 3σ.)
- **Three banked measurements are one object:** MODES' d1 ladder ≡ MOON's "m=1 channel MONOTONE(+)" ≡ S2's survival plateau. The 0.19% margin, the m=1 monotonicity, and the 0.3335 plateau are the same number in three normalizations.
- **The digit-bias form:** `p₀ − 1/3 =` weighted mean of `⟨α_t, α_{t+1}⟩` where α_t = (conditional top-digit distribution at coarse position t) − uniform (since `Σ_j a_j b_j = 1/3 + ⟨a−u,b−u⟩` for probability 3-vectors). d1>0 ⟺ **adjacent dlog positions have positively correlated new-digit biases** — a local smoothness of the digit field.

## What this does and does not prove
- **Proved (pen + machine-verified):** the identity, the corollary, the equivalences. The 7/15-vs-0.477 decider is now: **does γ_r(1) (the ×4-collision mass) keep growing — is q_r(1) > 1/3 forever?**
- **Not proved:** that q_r(1) > 1/3 persists. But the question now sits on the channel ledger with structure: (i) the **m=0 sibling is ESTABLISHED** — A_r(0)>0 with increments →~7/15 is the six-sightings linear divergence; the target is the m=1 version of a proven-at-m=0 phenomenon; (ii) **not universal in m** — MOON: m=3 decays, m=2 oscillates — so it is genuinely channel-specific; (iii) **new proof surface:** upstairs on (v,v') the measure is PRODUCT-geometric (FKG-perfect — the FKG failure was a property of the pushforward ν, not the source); all difficulty lives in the non-monotone digit functionals of X'−4X = Σ_k 3^k(2^{−S'_k}−4·2^{−S_k}); (iv) the survival/QSD framing (per-level survival → Perron eigenvalue of the killed pair-chain; target λ>1/3) is available **with the R29 caveat attached** (the pair-tail state is infinite-dimensional; truncations may not converge — same disease as R29; do not construct an operator).
- **Criticality reading:** BOTH channels sit at the 3q=1 critical edge — m=0 approaches from above with excess ~1/r (linear growth), m=1 with excess ~1.3e−3 shrinking. "Marginal by construction" (POINCARE) is this criticality, seen a third way.

## Status
**CHANNEL_ID (pen + verification):** ⭐⭐⭐**EXACT IDENTITY d1_r = A_r(1)/A_r(0)** (r≥2), corollary Re δ̂_r(n)=A_r(n)/A_r(0) (3∤n) — MODES ≡ charledger, verified exact-rationally (tower fold EXACT, identity==5-lag as Fractions, d1_2=2/35) + float to 2e-13 (r≤16) + d2 corollary + **Λ-weld to the ε ladder** (1.7e-5). ⭐⭐**UNIFICATION: S2's survival plateau q_s≈0.3335 IS q_r(1) ⟹ S2 statistically measured sign(d1) to r≈38 = POSITIVE; exact q₁₆=0.333751 continues seamlessly. MODES ladder ≡ MOON m=1 monotone(+) ≡ S2 plateau = one object.** ⭐**Seed: γ₀(1)=1 → γ₁(1)=2/3 (falls), then grows forever (A_r(1)>0 r=2..16). Decider = does the ×4-collision channel keep growing (q_r(1)>1/3)?** NOT proved: persistence. Proof surface: m=0 sibling established (7/15 slope); channel-specific (m=3 decays); product-geometric upstairs (FKG failure was the pushforward's); QSD framing carries the R29 caveat. Criticality: both channels at the 3q=1 edge (m=0 excess ~1/r, m=1 excess ~1.3e-3). Not at stake: R1–R30, R80–R82, all Thread-3 probes — this identity WELDS them. commit pending.
