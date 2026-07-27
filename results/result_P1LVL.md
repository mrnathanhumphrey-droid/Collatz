# RESULT — P1 (per-level): the cascade / d₁ is PHASE-carried — inverts the top-level P1 verdict (2026-07-26)

**Probe:** `probe_p1lvl.py`. Wilson's correction to P1: apply the bridge PER LEVEL. The tower telescopes
`γ_r(k) = 1 + Σ_{j≤r} A_j(k)`, `A_j(k) = Σ_{3∤a}|ρ̂_j(a)|²e(ak/3^j)` = the primitive shell at level j (since
`ρ̂_r(3^{r−j}b)=ρ̂_j(b)`). The bridge reaches each shell AT ITS OWN LEVEL, where 3∤a and |τ_a|=√q (gate-verified). So
scramble arg π̂_j at each level, bridge, accumulate. My earlier "disconnected/amplitude" verdict was an artifact of
applying the bridge only at the top level (where the increment is 10⁻³ of the total by construction). **SUPERSEDES the
verdict of result_P1** (its GATE and Gauss-sum facts stand; its "amplitude/disconnected" conclusion is inverted).

## GATES — both pass
- **GATE1** per-level bridge `A_j` == direct `A_j`: max abs **3.3e-15** (bridge valid at each level on 3∤a).
- **GATE2** telescoping `1 + Σ_{j=1}^7 A_j(k)` == certified `3^7⟨ρ,shift_k ρ⟩`: max abs **4.4e-16**. The tower shell
  decomposition is exact. (Confirms `ρ̂_r(3^{r−j}b)=ρ̂_j(b)` and that per-level `stationary_trunc` is tower-consistent.)

## The shells, and the d₁ re-derivation
Real per-shell `A_j(1)` (= `d₁^{(j)}·S_j`):
`A_1=−0.33333` (→ 1+A_1 = 2/3 = class mean exactly), then **all positive**: `A_2=+0.0272, A_3=+0.0092, A_4=+0.0040,
A_5=+0.0037, A_6=+0.0036, A_7=+0.0027`. cascade `c_1 = Σ_{j≥2}A_j(1) = +0.0504`. Each shell >0 is `d₁^{(j)}>0` — the
arc's positivity, shell by shell. **"Is the cascade positive" = "is d₁ positive."**

## THE FINDING — the cascade is PHASE-carried (class mean is amplitude)
Keep `|Ŵ_j|` exactly, scramble arg Ŵ_j at each level (positivity flag: inverse spectrum SIGNED both modes — formal):
| k | 3\|k | c_k real | c_k ZERO-phase | c_k RANDOM-phase |
|---|------|----------|----------------|-------------------|
| 1 | no | +0.0504 | +0.1414 | +0.1498 |
| 2 | no | −0.1906 | −0.0108 | +0.0558 |
| 3 | yes | −0.4243 | −0.1304 | −0.2979 |
| 4 | no | +0.2009 | +0.0015 | −0.1860 |
| 5 | no | +0.1018 | −0.0600 | −0.0629 |
| 6 | yes | −0.2821 | −0.3914 | −0.4415 |

- **The class mean A₁ = −1/3 survives every scramble** (phase-free frozen mod-3 marginal) — as pre-registered.
- **The cascade does NOT survive.** Keeping |π̂| and changing only the phase changes the shells wildly (A₂(1):
  +0.027 → −0.007 / +0.089; signs flip), so `c_k` is **not a function of |π̂| alone**. Ratios scatter 0.01–2.97; the
  dichotomy sign-structure (c: +,−,−,+,+,− real) is **not preserved** under scramble. **The cascade — the
  channel-distinguishing content, = Σ_{j≥2}A_j = the d₁ sequence = where 7/15-vs-0.4737 lives — is PHASE-carried.**

## Verdict — arg π̂ carries the cascade; |π̂| carries only the class mean
The Syracuse measure's **amplitude** |π̂| fixes the class mean (the symmetric/marginal part, 10/21). Its **phase**
arg π̂ (the measure's ASYMMETRY — natural, since d₁ is a directional lag-1 autocorrelation) carries the **cascade**
(the d₁ content, the deviation from 10/21, the whole S_∞ question). The zero-phase (symmetric) measure gives the wrong
cascade. So:
- **The day's |π̂| AMPLITUDE work (LAMBDA sup, HOMOG shape, VALPROFILE valuations) does NOT determine the fine S_∞
  structure** — it's the class-mean/amplitude side, already known (10/21).
- **P6 is now well-aimed:** the cascade lives in arg π̂ = `arg ∫_{ℤ₂} e(−{t·χ₃(z)}₃)dz`, the 2-adic domain integral.
  Analyzing the phase is exactly what's needed. The P1→P6 order is vindicated: P1 says phase carries the cascade, P6
  reads the phase off the domain.
- **"Prove d₁>0" is a PHASE statement:** each real shell A_j(1)>0, but generic phases give negatives — the positivity
  is a property of the actual arg π̂, not of |π̂|.

Wilson's exact split stands with the better address: `S_∞ = 10/21 + 2Σ_k 4^{−k}Σ_{j≥2}A_j(k)`, 10/21 = level-1
truncation (= S₂ — possibly not a coincidence, one look owed), correction = accumulated primitive increments (phase).
Not at stake: BRIDGE2, P4, GATE facts of result_P1, CHANNEL_ID, MEAN1, R1–R30. Cheap (0.3s).
