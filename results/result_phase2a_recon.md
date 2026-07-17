# Result — PHASE 2a BOUNDARY RECON (gates Phase 2b's proof target; no proofs attempted)

**Date:** 2026-07-16. Recon per the Phase-2a worksheet (Q1–Q6 + Q2b). Instruments: build_M sparse (L≤3), ESPRIT-on-mass-sequence, exact-ish stationary (v=64), generalized operator M(q,gen,λ). Committed pre-registrations; **deviations reported AS deviations.**

**Headline: Q2's fork is (b) — the q=3 pair is distinct-but-coalescing, so Phase 2b's theorem is an L→∞ COALESCENCE statement, not finite-L algebra. And Q6 REFUTED the pure-phase hypothesis: `|⟨g⟩ mod q|=2` alone does NOT collide (M(q,−1) stably gapped at 0.831); the q=3 collision is specific to `⟨2⟩`'s q-adic lift (`2≡−1 mod q ⟺ q=3`, `⟨2⟩ mod 9` has order 6, not {±1}). Phase 2b's statement changes: the boundary is d=ord_q(2)=2 via the lifted `⟨2⟩` structure — NOT decouplable to abstract {±1} sign combinatorics.**

Probes: `probe_phase2a_recon.py` (Q1–Q5), `probe_phase2a_q2b_q6.py` (Q2b, Q6). Logs: `result_phase2a_recon_log.txt`, `result_phase2a_q2b_q6_log.txt`.

## Q1 — Mode census at the collision (q=3, L=2,3)
build_M q=3 L=3 top-6 |λ|/(1/3): **1.000, 1.000**, 0.974, 0.974, 0.910, 0.910.
- **TWO eigenvalues at ratio 1.000** (within ≡ Perron, a Jordan pair — since `1/q = 1/3` merges the within mode onto Perron). ✓ pre-reg "TWO".
- **DEVIATION:** the third mode is **0.325 (ratio 0.974)** — this is the **cross mode climbing toward 1/3** (r₃→1 as L→∞), NOT R25's unrelated "0.273". The "something else" is the cross mode mid-collision, not a foreign mode.
- **Adjudication:** at finite L, within+Perron are already merged (Jordan pair at 1/3); the cross is a distinct third mode climbing to 1/3.

## Q2 — Exact Jordan vs exceptional-point limit (the fork; NAMES the theorem)
| L | dim | gap \|λ₁−λ₂\| | eigvec overlap \|⟨v₁,v₂⟩\| |
|---|---|---|---|
| 1 | 12 | 0.889 | 0.000 |
| 2 | 324 | 2.91e-3 | 0.998 |
| 3 | 8748 | **9.96e-5** | **0.99999** |

- **Fork (b) CONFIRMED:** eigenvalues stay *distinct* at finite L (gap ≠ 0) but the gap → 0 AND eigenvectors → parallel (overlap → 1), both **monotonically in L**. The exceptional point is reached only at **L = ∞**. ✓ pre-reg (b).
- **★ THEOREM SHAPE (Phase 2b):** an **L→∞ coalescence-rate statement** — the cross–Perron gap → 0 as L→∞ at d=2 (forced), stays bounded away at d≥3 (forbidden). NOT finite-L algebra.

## Q2b — EP Puiseux signature (physics-frame, weight-deform knob f(v)=v)
- q=3 (L=3): gap(ε) is **V-shaped** — decreases to a minimum ~3.3e-5 at ε≈3e-4, then rises (finite-L detuning floor). Log-log regression slope **0.575** (√ε-ish). The V-shape = a trajectory passing *near* an EP.
- q=7 (L=2): gap(ε) **flat ≈0.13**, slope ≈0 — the top pair is robustly separated (no EP).
- **DEVIATION / partial:** the clean pre-registered slopes (½ vs 1) did **not** cleanly materialize — q=3 is non-monotone (finite-L detuning), and the q=7 gap-metric measured the complex top-pair distance (robustly flat), not a linear splitting. **The CONTRAST (q=3 EP-proximity V-shape vs q=7 flat) corroborates (b), but Q2b is suggestive, not a clean √ε confirmation.** Q2's primary data (gap→0 + overlap→1, monotone) already settles the fork.

## Q3 — Derive-7/15 hook (q=3, feasibility)
`c_k = X_k − X_{k−1}` (the k·(1/3)^k Jordan slope): 0.476, 0.462, 0.464, 0.4655, 0.4662, 0.4655, **0.4659** (k=8). vs 7/15 = 0.46667 → **0.16% — PASS (<1%).** The top-2 Jordan-block restriction has off-diagonal coupling **0.019** at ~equal eigenvalues (a genuine chain). **Phase 2b's bonus route (Jordan chain → slope 7/15) is feasible.**

## Q4 — Separation control at d≥3 (q=7 primary; q=5 in (c) only)
- q=7 ESPRIT modes: Perron **0.333 (ratio 1.000)**, cross **0.127 (ratio 0.3798)**, sub 0.052.
- **DEVIATION:** ESPRIT does NOT isolate the within mode 1/7 (ratio 3/7=0.4286) — within (0.143) and cross (0.127) are only **0.016 apart (9%)**, so ESPRIT merges them (returns 0.38). **within = 3/q is EXACT by Lemma 5** (the Pythagoras homogeneous rate), a *theorem*, not an ESPRIT-measurable mode.
- Three theoretical modes {1/3, 1/7, 0.127}: **min pairwise gap = 0.016** (within–cross); Perron well-separated (gap 0.19). So d≥3 is "separated" but the within–cross gap is modest.
- q=5 (c)-coords only [concealing prime, no modal]: cross-rho `0.534, 0.508, 0.624, 0.630, 0.628, 0.609 → ~0.62`, the ±0.60 oscillatory pair as expected.

## Q5 — Wieferich spot (q=1093)
`d = ord_1093(2) = 364` (≥3 ⇒ gapped by theorem-shape); `2 ≢ −1 mod 1093` (q=1093 is NOT a boundary — only q=3 has 2≡−1). k=1: ‖π₁‖²=1/3. **k=2 hit the heavy-compute guard (n·vmax ≈ 7.6e7 > 30M) — SKIPPED**, honest. No collision at 1/3 in reach; the gap onset is index-shifted to higher k (R35), beyond the cheap exact route.

## Q6 — Decoupling / Mr. Potato Head (recon on Phase 2b's HYPOTHESES)
Generalized `M(q, gen, λ)`: phases ⟨gen⟩ mod q^L, weights ∝ λ^v. `M(q,2,1/2) = build_M`.

**GATE G6:** `M(3,2,1/2)` collision (banked r₃=1), `M(7,2,1/2)` gapped. ⚠️ **The ESPRIT-on-mass collision detector is UNRELIABLE for the q=3 Jordan** (gave 0.831/1.000/0.559 at L=1/2/3 — the repeated root + k-prefactor confuse ESPRIT). It IS reliable for *gapped* operators (see below). So q=3's collision is taken from banked values, not this detector.

**TEST A (swap phase gen, keep weight) — REFUTED.**
| operator | \|⟨gen⟩ mod q\| | r (L=1,2,3) | verdict |
|---|---|---|---|
| M(7,−1,½) | 2 | 0.831, 0.831, 0.831 | **GAPPED (stable)** |
| M(5,−1,½) | 2 | 0.831, 0.831, 0.831 | **GAPPED (stable)** |
| M(13,−1,½) | 2 | 0.831, 0.831, 0.831 | **GAPPED (stable)** |
| M(13,3,½) | 3 | 0.573 | gapped |

- **`|⟨g⟩ mod q| = 2` is NOT sufficient for collision** — every gen=−1 case is stably gapped at 0.831 (flat across L, so a real gap). The pre-registered "collision follows ord(gen)=2, weight-blind + 2-blind, pure {±1} combinatorics" is **REFUTED.**
- **Why q=3 is different (the mechanism):** `⟨−1⟩ mod q^L = {1,−1}` at every L (never lifts), but `⟨2⟩ mod 9 = {1,2,4,5,7,8}` has **order 6** — it lifts q-adically. And −1 is the *only* order-2 element in F_q*, and it never lifts ⇒ **no genuine order-2-lifting analog of q=3 exists at any other prime.** The q=3 collision is specific to `2 ≡ −1 mod q ⟺ q=3`, carried by the *lifted* `⟨2⟩ mod 3^L` structure.

**TEST B (swap weight λ, keep gen=2 at q=3) — INCONCLUSIVE.** r came out 1.0/1.0/0.865/1.0 for λ=1/3, 2/5, 3/5, 1/2 — but this is exactly the q=3 Jordan case where the ESPRIT detector is unreliable, so the values (esp. the λ=3/5 "unglue") are not trustworthy. Weight-blindness is *suggested* (3/4 give r=1) but **not established.** The dilution Σp_v² moves as designed (0.501 → 0.432 → 0.274 → 0.344).

## ADJUDICATION (per Q; Q2 + Q6 set Phase 2b's target)
| Q | one-line adjudication |
|---|---|
| Q1 | TWO at 1/3 (within≡Perron Jordan pair); 3rd = cross climbing (0.325), not R25's 0.273 |
| Q2 | **(b): distinct-but-coalescing → theorem is L→∞ coalescence, not finite-L algebra** |
| Q2b | q=3 EP-proximity (V-shape) vs q=7 flat corroborates (b); clean √ε slope NOT clean (partial) |
| Q3 | 7/15 hook feasible (slope 0.4659, 0.16%); Jordan coupling present |
| Q4 | 3 modes; within=3/q EXACT by Lemma 5 (not ESPRIT-isolable); within–cross gap 0.016 |
| Q5 | d=364≥3, not a boundary; gap beyond cheap reach (index-shifted), no collision seen |
| Q6 | **A REFUTED: order-2 phase alone ≠ collision (M(q,−1) stably gapped 0.831); B inconclusive** |

**★ PHASE 2b STATEMENT (revised by this recon):** the entrance-exam theorem is **not** "d=2 alone, pure {±1} phase combinatorics, weights quotient out." It is: **at `d = ord_q(2) = 2` (⟺ q=3, via `2 ≡ −1 mod q`), the cross–Perron pair COALESCES as L→∞ (an exceptional point / defective limit), forced by the q-adic lift of `⟨2⟩ mod 3^L`; at d≥3 the pair stays gapped.** The proof must use `⟨2⟩`'s lifted structure (a genuine order-2 phase that lifts to order `2·q^{L−1}`), which is realizable *only* at q=3 — not abstract sign combinatorics. Whether the weight base is inessential is left OPEN (Test B inconclusive; a reliable q=3 collision detector is needed to settle it).

## Not at stake
R1–R46, Phase 0/1, G0c′. No r_q value changes. This is recon: it names Phase 2b's theorem shape (L→∞ coalescence) and corrects its mechanism hypothesis (⟨2⟩ q-adic lift, not generic {±1}).

_Reporting discipline: the pre-registered Q6-A ("order-2 ⟹ collision") and the clean Q2b slopes were REFUTED/not-clean and reported as such, not as "consistent with." The ESPRIT detector's unreliability on the q=3 Jordan is disclosed (it invalidates Test B and the G6 q=3 read — banked values used instead), rather than trusting the flaky 0.559/0.865. The Test A refutation rests only on the RELIABLE gapped-case measurement (flat 0.831 across L). The q-adic-lift mechanism (⟨2⟩ mod 9 order 6 vs ⟨−1⟩ order 2) is exact arithmetic, not a fit._
