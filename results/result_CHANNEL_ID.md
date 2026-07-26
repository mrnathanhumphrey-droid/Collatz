# PEN RESULT + Probe CHANNEL_ID — the dominant mode IS the m=1 collision channel — **EXACT IDENTITY (derived by pen, verified exact-rationally): `d1_r = A_r(1)/A_r(0)` for r≥2 — the MODES frontier and the charledger are THE SAME OBJECT.

> **⚠️ AUDIT CORRECTIONS (Wilson, 2026-07-25; quantified in `probe_channel_audit.py`, see the CORRECTIONS section at
> the bottom).** (1) The headline **"S2 measured sign(d1) positive throughout to r≈38" is RETRACTED** — q−1/3 decays
> ~tenfold across the window ("plateau" was flat at reporting precision only); S2's per-level SE on q is
> 1.7–3.5×10⁻⁴ vs a needed ~3×10⁻⁵ at r=38 (per-level z = 0.09–0.18σ, unresolvable; window-aggregate 1.9–3.9σ).
> Honest claim: **S2 and the exact ladder agree where both are informative** — consistency, not extension. The 2.3σ
> r≈36 rollover sits exactly where excess is smallest and resolution worst. (2) **The m=0-inheritance line is CUT**
> — γ_r(0) diverges (increments forced positive); γ_r(1) converges (~0.738), so its increment-positivity = monotone
> approach from below, never overshooting. Different asymptotics; nothing to inherit. (3) **A_r(0) = S_r exactly**
> (= banked R7: X_r is the cumulative sum of S), so the identity's final form is **`d1_r = A_r(1)/S_r`**. (4) The Λ
> weld is **identical** (2e-13 vs exact-ε Λ, K=40); V4's 1.7e-5 was the banked LAM_NU 5-digit rounding (2.1e-5
> half-ulp), truncation 6.7e-10. (5) γ_∞(1) ∈ [0.7361, 0.7382]; low-complexity candidates 14/19, 31/42, 45/61,
> 59/80; (1+S_∞)/2 = 0.7383 at the band edge given S_∞'s ±0.001 — underdetermined, no claim. (6) Wilson's blind
> anchors (d1_2 = 2/35, p_2(1) = 102/1323 = 0.077098) landed independently — seed doubly confirmed. Surviving lead:
> **upstairs FKG sandwich** (monotone events bracketing {X'≡4X} on the product (v,v') lattice, gap < 4×10⁻⁴) —
> Wilson's pen. Caveat carried: "FKG-perfect measure, non-monotone observable" is the *generic shape of an FKG
> failure* — a lead, not a route.

**(Original headline, as corrected by the audit):** Corollary (same proof): `Re δ̂_r(n) = A_r(n)/A_r(0)` for every 3∤n — the whole mode ladder is the channel ledger normalized by the m=0 increment (= S_r, the R7 cumulative-sum weld). Equivalences: d1>0 ⟺ A_r(1)>0 (γ_r(1) strictly increasing = MOON's measured m=1 monotone+) ⟺ q_r(1) = Pr[X'≡4X mod 3^{r+1} | X'≡4X mod 3^r] > 1/3 (the conditional lifting of the ×4-collision beats the neutral cost) ⟺ adjacent dlog positions have positively correlated new-digit biases. The S2 splitting run's survival fractions q_s ≈ 0.3335 (m=1 channel, r=17..38) ARE this object — exact q₁₆ = 0.333751 continues seamlessly into them; per the audit, this is **consistency where both are informative** (per-level sign at the deep end unresolvable), not an extension of positivity to r≈38. Three previously separate banked measurements (MODES d1 ladder, MOON m=1 monotonicity, S2 survival) are one object.

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

---

## AUDIT CORRECTIONS - quantified (probe_channel_audit.py, logs/channel_audit_run.log, 2026-07-25)

**A) A_r(0) = S_r (the R7 weld) => d1_r = A_r(1)/S_r.** Exact Fractions at r=2: A_2(0) = 10/21 = S_2 EXACTLY
(= Wilson's blind anchor A_2(1) = (2/35)(10/21) = 4/147). Float: rel 3.6e-16 (r=6) to 9.5e-15 (r=12); r=14/16 at
1.3e-8/6.4e-9 = the precision of the EPS chain itself (built from 5-digit LAM_NU). The r=3..5 exact-Fraction
comparison FAILS as rationals - this is the vmax=120 truncation in build_nu_exact (2^-120 mass defect: visible to
exact rationals, invisible to float), NOT a failure of the identity; the banked exact-eps json used exact geometric
resummation. Structurally the weld is banked R7: A_r(0) = X_r - X_{r-1} = S_r (X = cumulative sum of S).

**B) The Lambda identity is IDENTICAL.** Lam_r = sum_{m>=1} 4^-m A_r(m) vs the exact-eps (eps_{r+1}-eps_r)/2, K=40:
rel 1.0e-14 to 2.0e-13 for r=2..7 (exact-eps rows); ~1e-12 for r=8..10 (float-eps). V4's 1.7e-5 decomposed: banked
LAM_NU half-ulp rel = 2.1e-5 (dominates); K=13 channel truncation = 6.7e-10. "Reproduces to 1.7e-5" becomes
"identically (limited only by the banked constants' print precision)."

**C) S2's error bar ON q - the plateau downgrade, with numbers.** SE(q) per level: naive sqrt(q(1-q)/N)/sqrt(20) =
1.67e-4; ESS-discounted (uniq-frac 0.223) = 3.53e-4. Propagated excess 4.23e-4 * 0.887^(r-16): at r=38 = 3.0e-5 =>
per-level z = 0.18/0.09 - UNRESOLVABLE. Window-aggregate (17..38): mean excess 1.40e-4 => z = 3.9 (naive) / 1.9
(ESS). => "positive throughout to r~38" RETRACTED; honest = consistency where informative; the strong statistical
statement remains the telescoped eps-hat rise (z~3.4); the 2.3-sigma r~36 rollover lives exactly where excess is
smallest + resolution worst - the two must be read together.

**D) gamma_inf(1) scan.** gamma_16(1) = 0.730013, A_16(1) = 9.14e-4; extrapolation band over rho in [0.87,0.90]:
[0.73613, 0.73824]. Low-complexity residents: 14/19 (0.736842), 31/42 (0.738095, in-family denominator 2*3*7),
45/61 (0.737705), 59/80 (0.737500) (+8 more q<=141). Named: (1+S_inf)/2 = 0.73830 at S_inf=0.4766 sits just
OUTSIDE the band but S_inf's own +-0.001 reaches the edge - not excludable. Underdetermined; banked as a scan, no
claim. A second closed form on the ledger remains open.

**The surviving lead (Wilson takes the pen):** the upstairs FKG sandwich - the (v,v') measure is product-geometric
(FKG for free); the whole content is whether {X'=4X mod 3^(r+1)} admits a monotone reformulation or can be
sandwiched between monotone events with gap < the 4e-4 excess. Caveat carried: "FKG-perfect measure, non-monotone
observable" is the generic shape of an FKG failure - a lead, not a route.
