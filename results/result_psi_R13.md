# Probe R13 — lobe constant / ψ-existence / transport — **D gate PASS; A/B measured; C decider OPEN**

**Date:** 2026-07-21  Reuses R7/R9/R10; exact where marked. Probe `probes/probe_psi_R13.py`. Gates Wilson's
closed-form weight, the ψ-uniformity discriminator, the ψ-existence decider, and the renewal-is-the-chirp identity.

**Weight closed form verified:** Re w(x) = 15/(2D) − 1/2, D = 17−8cos2πx; Re w(0)=+1/3, Re w(½)=−1/5 ✓; ⟨Re w⟩
numeric = −1.3e−17 (**mean exactly 0**); ∫₀¹|Re w|dx numeric = **0.160861** = K0 = 1−(2/π)arccos(1/4) (Wilson's
0.160867 to 5 digits; arccos rounding). ‖w‖²=1/15 = 1/(4²−1), the N=2 Mersenne member — 15 is the weight's own L²
norm, the denominator of 7/15.

## R13-D — RENEWAL-IN-ORBIT GATE (forced): **GATE PASS**
Both routes reproduce the frozen Syrac(ℤ/3^{n+1}) exactly from Syrac(ℤ/3^n) + Geom(2), n = 1…5:
- frequency recursion **X_{n+1} = 1 + 3·2^{−v}·X_n** → μ_{n+1} exact (True ×5);
- orbit-coordinate **t′ = β(2^{−v}·4^t)** (via R12's β-table) → same transition exact (True ×5).

**The chirp β is not a coordinate change — it is the second half of the dynamics** (exponentiate, multiply by the
unit 2^{−v}, apply β). The corpus's R78/R81/R81b/R12-A chirp structure was the renewal's transfer function the whole
time. Two routes, one transition, exact. Walk-back #32 **not** incurred.

## R13-A — LOBE CONSTANT vs UNIFORM (discriminator; measurement, NO verdict)
| r | L+M measured | S_r·K0 [uniform-ψ] | abs diff | ratio meas/unif |
|---|---|---|---|---|
| 3 | 0.0723313 | 0.0742495 | −0.00192 | 0.97417 |
| 4 | 0.0744961 | 0.0746741 | −0.00018 | 0.99762 |
| 5 | 0.0726717 | 0.0748833 | −0.00221 | 0.97047 |
| 6 | 0.0715506 | 0.0749885 | −0.00344 | 0.95415 |
| 7 | 0.0738668 | 0.0748795 | −0.00101 | 0.98648 |

**The ratio is NOT a flat 0.96 — it oscillates 0.954…0.998** (mean ≈0.976). ψ is **not uniform** (measured mass is
below the uniform-ψ prediction at every r), but the deficit **oscillates in r** (0.2%–4.6%) rather than sitting at
Wilson's pre-registered fixed ~4.1%. The deficit and its oscillation are the same period-9 structure seen in the Λ
ledger, now in the total |Re w|-mass. No verdict; the shape question is R13-B.

## R13-C — ψ-EXISTENCE DECIDER: **OPEN at r≤7** (the differences do not uniformly shrink)
γ_r(τ_m) across r = 1…7 (exact), with successive differences:

| m | v₃ | γ_r (r=1…7) | successive diffs | trend |
|---|---|---|---|---|
| 1 | 0 | 0.667→0.717 | +.027,+.009,+.004,+.0037,+.0036,+.0027 | **monotone shrinking** (converging) |
| 2 | 0 | 0.667→0.476 | −.177,+.005,−.002,+.001,**−.012**,−.007 | non-monotone (jump at r=5→6) |
| 3 | 1 | 1.667→1.242 | −.238,−.094,−.010,−.044,−.035,−.003 | non-monotone |
| 4 | 0 | 0.667→0.868 | +.150,+.062,−.017,−.001,+.009,−.002 | oscillating, small |
| 9 | 2 | 1.667→2.089 | +.476,−.231,+.031,+.027,**+.073,+.046** | **tail GROWING** (drifting up) |
| 27 | 3 | 1.667→2.377 | +.476,+.462,−.232,−.004,+.057,−.048 | non-monotone |

**Verdict: ψ-existence is neither confirmed nor cleanly killed at r≤7.** Only m=1 shows clean monotone convergence.
The rest are dominated by a **common oscillation with a turnover near r=6** (m=2, m=9, m=27 all jump there) — the
period-9 structure, inherited by the γ's. m=9 in particular has *growing* tail differences (+0.073 at r=5→6). So the
differences do **not** uniformly shrink: convergence, if it holds, is **oscillatory** and its damping cannot be
established from r≤7. This is R9-C's boundedness with the convergence question still open — the decider needs r≥8
(the exact wall, μ₈) or an analytic damping argument, not a deeper numerical read at this depth. Reported plainly per
pre-registration: **the ψ programme is not vindicated here, and not dead — it is genuinely undecided at reachable
depth.** The pen adjudicates whether the r=6 turnover is a damping oscillation (ψ exists) or not.

## R13-B — WHERE THE DEPLETION SITS (measurement, NO fit)
Layer profile |θ̂(k)|² binned into 12 angle-bins, measured mass ÷ uniform (S_r·bin-fraction). At r=7 (best
resolved; symmetric by conjugation):

| x (bin center) | Re w(x) | meas/unif (r=7) | (r=6) | (r=5) |
|---|---|---|---|---|
| 0.042 (≈0) | +0.309 | **0.921** | 0.838 | 0.861 |
| 0.125 | +0.161 | 1.100 | 1.073 | 1.119 |
| 0.208 | +0.002 | 1.069 | 1.305 | 1.235 |
| 0.292 | −0.107 | 0.944 | 0.835 | 0.765 |
| 0.375 | −0.169 | 0.978 | 0.833 | 0.881 |
| 0.458 | −0.197 | 0.989 | 1.104 | 1.176 |

**The depletion is concentrated near x≈0 (the trivial-character end): the x≈0 bin is below uniform at all three r
(0.86, 0.84, 0.92), the most robust deviation.** Per Wilson, a depletion at θ≈0 is a **decay statement for ν̂ near
the trivial character — a mixing statement in the most standard form.** The intermediate bins oscillate with r (the
x≈0.208 bin swings 1.24→1.31→1.07), the same period-9 modulation; by r=7 the whole profile is flatter (ratios
0.92–1.10) than r=5,6. The near-x≈0 depletion is the candidate signal; no verdict, table verbatim.

## R13-E — R85 rung-1 debt (feasibility only, NOT run)
R85 rung-1 (operator-DFT chirp identity, PASSED but owed the r=5/n=8 extension) is the **same object as this
transport**. With the β/U tables now built (R11/R12) and the support law (R12-A) pruning U to ~1/6 density
(k≡ξ mod 3, block-diagonal by v₃), the n=8 (N=3⁸=6561) extension is a dedicated FFT/Bluestein-route probe
(O(N log N) per block) — **cheaper than July but not free.** Deferred, as instructed.

## Status
**[UPDATE 2026-07-21: R13-C RESOLVED in R14.** Wilson's follow-up argument — γ_n(τ_m) is a partial sum of the
frozen A-series (A_r=C_{r+1}(m)/3), bounded ⟹ no nonzero A_r-limit ⟹ **no non-uniform ψ exists** (limiting-shape
reading dead; the object is the deviation field δ_r). The "open decider" below is closed in the *ruled-out*
direction. See result_deviation_field_R14.md.]**

**R13: D gate PASS** (the chirp β *is* the renewal — X_{n+1}=1+3·2^{−v}X_n and t′=β(2^{−v}4^t) both exact, n=1…5;
#32 not incurred), **A/B measured** (ψ is not uniform — an oscillating ~0.2–4.6% deficit — with the depletion
concentrated near the trivial character x≈0, a candidate mixing signal), **C decider OPEN** (γ_r(τ_m) convergence is
undecided at r≤7: only m=1 monotone, the rest dominated by a period-9 oscillation with a common r≈6 turnover and
m=9 tail-growing — ψ-existence needs r≥8 or an analytic damping argument). The R12-C "ψ-existence" language is
scoped down (correction applied to result_lobes_R12.md): what stabilizes is the |Re w|-moment, and whether ψ *exists*
is exactly this open decider. **Still owed (pen):** ψ-existence (the r≈6 turnover's damping), then Σ_{r≥1}Λ_r=−1/10
as its transient — via the now-explicit transport operator θ̂′(k)=𝔼_v𝔼_t[e(k·β(2^{−v}4^t)/3^r)] (ψ its fixed
point). No fitting; exact γ and renewal gate, labeled numeric lobes/bins; the decider's openness reported as
openness, not smoothed.
