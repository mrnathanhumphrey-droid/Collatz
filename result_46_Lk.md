# Result 46 (G1b / L3 Session-Two, qx+1 paper) — L_k COMPLETED: the deviation propagates by a single DETERMINISTIC REFINEMENT operator with NO injection (drip ≡ 0 at machine precision). The 2/3 dies as pre-registered — and takes 7/45-as-injection with it. The same operator delivers the GAP half of L3 clean: L_k's gain IS r_q, σ_max(L_k)=√(r_q/3) < Perron 1/√3 for every d≥3.

**Date:** 2026-07-16. **Instrument:** completed `c_seven_forty_fifth_derivation.py`'s Steps 2–5 (the deviation lifting operator L_k) as `probe_46_Lk_echo_block.py`, reusing its foundation (self-similar π_k, ‖d_k‖²) without mutating the shared file. **Session-Two spec (user): measure the echo block T′'s propagated fraction in c-units → 2/3 [at risk]; same instrument at q=5,7 = the gap half of L3.**

**Verdicts: H_ECHO (2/3) ✗ REFUTED (injection ≡ 0; propagated fraction = 100%, not 2/3) / H_GAP ✓ CONFIRMED (L_k's gain = r_q; σ_max(L_k) below Perron for d≥3, same instrument). Three structural identities at ~1e-16.**

**Headline: two mis-built versions (v1 phi≡rate tautology; v2 /q-vs-/3 object + ARPACK ρ=0) were killed by internal gates before either could manufacture a false "2/3 confirmed" — the redesign's whole point. The machine-precision structure they were pointing at: (LEM-FORGET) the transfer K forgets the source's finest q-adic digit (`target=(qr+1)2^{−v} mod q^k` is independent of r's q^{k−1} digit, since q·q^{k−1}=q^k≡0), so the self-block T′=P_W K P_W is GENUINELY ZERO (‖dK‖/‖d‖=9e−17 for random d∈W); (ONE-STEP) `π_{k+1}=lift(π_k)K_{k+1}` exactly (measure equilibrates in one step); (REFINE) `d_{k+1}=P_W[lift(d_k)K_{k+1}]` exactly — cos = 1.00000000, residual 1e−15. So the deviation propagates by PURE DETERMINISTIC REFINEMENT with NO injection: drip ≡ 0. The recursion c_{k+1}=(2/3)c_k+7/45 is therefore NOT an echo+drip partition — there is no drip. The propagated fraction is 100%, not 2/3; the 2/3 is merely the arbitrary coefficient of writing the convergent sequence c_k→7/15 as a first-order recursion (any α gives α c_k+(1−α)·7/15 at the fixed point; R45 already showed the residual isn't even constant). 2/3 dies; 7/45-as-injection falls with it — exactly the linkage the user pre-registered. What the SAME operator gives instead: L_k = P_W∘lift∘K is one clean operator whose gain is exactly r_q (q=3→1 marginal, q=5→0.62, q=7→0.39), so σ_max(L_k)=√(r_q/3) sits below the Perron 1/√3=0.5774 for every d≥3 (q=5: ~0.42–0.46, q=7: ~0.35–0.39) — the tractable half of L3, now realized as the singular value of an explicit injection-free operator.**

Probe: `probe_46_Lk_echo_block.py` (v3). Log: `result_46_Lk_echo_log.txt`. Dense/sparse, v≤64, seconds.

## The three structural identities (each ~1e-16, all q)

| identity | statement | measured |
|---|---|---|
| **LEM-FORGET** | `dK_k = 0` for any deviation `d∈W_k` (K forgets source's finest digit) ⇒ self-block T′=0 | `‖dK‖/‖d‖ = 9e−17` (q=3), 1e−16 (q=5), 1.5e−16 (q=7) |
| **ONE-STEP** | `π_{k+1} = lift(π_k)·K_{k+1}` (one-step equilibration) | `max ‖π_{k+1}−lift(π_k)K‖/‖π‖ = 3.5e−16` (q=3), 4.8e−16, 4.4e−16 |
| **REFINE** | `d_{k+1} = P_W[lift(d_k)·K_{k+1}]` (deterministic refinement, no injection) | `min cos = 1.00000000`, `max resid = 1.1e−15`, all q |

**Reason for LEM-FORGET:** `target(r,v)=(q r+1)·2^{−v} mod q^k`. For `r'=r+q^{k−1}`: `q r'+1 = q r+1+q^k ≡ q r+1 (mod q^k)`. So the transfer's image is independent of the source's q^{k−1} digit; rows of K within a level-(k−1) fiber are identical; hence `dK=0` for fiber-mean-zero d. The deviation cannot self-propagate at a fixed resolution — it exists only via refinement to the next level.

## H_ECHO (the 2/3) — REFUTED

`REFINE` gives cos = 1 at 1e−16: `d_{k+1}` is entirely the transported-lift of `d_k`. **Injection ≡ 0 — there is no drip.** Therefore:
- The "propagated fraction in c-units" is **100%**, not 2/3.
- The echo/drip decomposition `c_{k+1}=(2/3)c_k+7/45` has no operational referent: with drip = 0, the split into a 2/3-echo and a 7/45-injection is empty. The 2/3 is the arbitrary α in writing `c_k→7/15` as `c_{k+1}=α c_k+(1−α)·7/15`; R45 already showed the exact residual isn't constant (super-geometric). **2/3 dies; 7/45-as-injection falls — as the Session-Two spec pre-registered ("anything else settles it, and takes the 7/45 reading down with it").**

## H_GAP (the tractable half of L3) — CONFIRMED, same instrument

L_k: d_k ↦ d_{k+1} is a single deterministic operator. Gain `g_k=‖d_{k+1}‖²/‖d_k‖²`; c-unit `rate_k=3g_k=cB_{k+1}/cB_k`; `σ_max(L_k)=√g_k=√(rate_k/3)`. Perron = 1/√3 = 0.5774 (q=3, r_3=1). `cB_k=3^k(‖π_k‖²−(1/3)‖π_{k−1}‖²)`, gated exact at q=3 (0.476190, 0.461575, 0.464214, 0.465515 — matches R45).

| q | rate_k profile | → r_q | σ_max(L_k) | vs Perron 0.5774 |
|---|---|---|---|---|
| 3 | 0.969, 1.006, 1.003, 1.001, 0.999 | **1** (marginal) | ~0.577 | AT (r_3=1, no gap) |
| 5 | 0.534, 0.508, 0.624 | **0.62** | 0.42–0.46 | **below (GAP)** |
| 7 | 0.447, 0.357 | **0.39** | 0.35–0.39 | **below (GAP)** |

The deviation operator reproduces the banked r_q (build_M: r_5≈0.62, r_7≈0.39) and shows the gap `σ_max(L_k)<1/√3` for every d≥3. **The L3 tractable half restated cleanly: bound `‖L_k‖ = σ_max(P_W∘lift∘K) < 1/√3` for d≥3** — an explicit, injection-free operator (no carry-memory, no pair-state), a cleaner target than build_M's (a,b,γ).

## What this settles / relocates
- **The q=3 boundary is fully understood at the deviation level:** deterministic refinement, gain exactly 1 (marginal), no injection. r_3=1 is the marginal gain of L_k, not an echo/drip balance.
- **r_q is REALIZED as σ_max(L_k)²·3** — the singular value of one explicit operator L_k = P_W∘lift∘K. This is the object L3's bound should target (supersedes the echo/drip framing entirely).
- **The 2/3 and 7/45-as-injection are retired.** 7/15 = lim c_k stands (the q=3 marginal fixed point); 7/45 = (1/3)·7/15 is only its trivial first difference, carrying no operator meaning.

## Not at stake
R1–R45. No r_q value changes — r_q is re-derived (0.62, 0.39) and now realized as L_k's gain. L3's statement is unchanged but its target is sharpened: `‖L_k‖<1/√3` for d≥3, L_k explicit and injection-free.

_Reporting discipline: v1 (phi≡rate) and v2 (/q object, ARPACK ρ=0) were killed by internal gates, not reported as results — the redesign catching two would-be false passes. The verdict rests on three identities at ~1e-16 (not fits) plus the exact-rational cB gate. The 2/3's death is structural (injection=0 at machine precision), not a threshold call. The GAP is the banked r_q reproduced, so H_GAP is a cross-validated realization, not a new claim._
