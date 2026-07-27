# RESULT — P1: the bridge is BLIND to the channel; arg π̂ cannot carry γ (amplitude, but NOT the |π̂| work) (2026-07-26)

**Probe:** `probe_p1.py`. Wilson's amplitude-vs-phase gate. Refactor BRIDGE2's exact χ/DL/τ/offset transfer, gate
bit-for-bit, then test whether scrambling arg π̂ (= arg Ŵ, the forward additive FT) can move γ. **The gate
discriminated — by revealing the bridge cannot reach γ's support.**

## GATE — refactor is faithful to BRIDGE2's actual domain (3∤a)
Factored bridge `ρ̂(a)=(1/τ_a)Σ_t χ̄_a⁺(t)μ_Y(t)` (μ_Y from the offset formula on Ŵ) vs BRIDGE2 direct `Σ_Y ν_Y χ_a(Y)`:
- **3∤a (primitive, BRIDGE2's validated domain — it only ever tested a∈{1,2,4,5,7}): max rel = 7.5e-13, BIT-FOR-BIT.**
  `|τ_a| = 81 = √q` for all 1458 primitive characters.
- **3|a (a>0): bridge FAILS** (max rel 1.02). `|τ_a| ≈ 1.9e-11` — the Gauss sum **vanishes for imprimitive characters**.
  This is NOT a refactor bug: BRIDGE2 never covered 3|a.
- **Dominant carrier a=N/3=3^{n−1}: bridge gives |ρ̂|²=0.0001, truth is 1/3; τ=9e-12.** The bridge is blind to it.

## THE STRUCTURAL FINDING — bridge domain ∩ γ support = ∅
`γ_n(k)−1 = Σ_{a≠0}|ρ̂(a)|²e(ak/3^n)` is carried **entirely by 3|a** (the coarse, imprimitive frequencies — P4):
| k | γ−1 | from 3∤a (bridge domain) | from 3\|a (channel) | 3\|a share |
|---|-----|--------------------------|---------------------|-----------|
| 1 | −0.283 | +0.0027 | −0.286 | **100.9%** |
| 2 | −0.524 | −0.0065 | −0.518 | **98.8%** |
| 3 | +0.242 | −0.0033 | +0.246 | **101.4%** |
| 6 | +0.385 | −0.0010 | +0.386 | **100.3%** |

The bridge's domain (3∤a, fine, primitive, τ=√q) is **disjoint** from γ's support (3|a, coarse, imprimitive, τ=0). The
3∤a part — the only thing arg π̂ maps to through the bridge — contributes **~0%** to every channel.

## VERDICT — arg π̂ cannot carry the channel, and neither can |π̂|
`γ = Σ|ρ̂|²e(ak/N)` is amplitude-only in ρ̂ (Wiener–Khinchin — arg ρ̂ never enters). Its support is the coarse 3|a
frequencies = the **frozen mod-3^j marginals** (tower-fixed: |ρ̂(N/3)|²=1/3 is the mod-3 marginal exactly, P4; the
cascade is the mod-9, mod-27, … marginals, r-stable). The bridge maps arg π̂ ONLY to 3∤a, where it's blind to γ.
So:
- **All three of Wilson's readings resolve to AMPLITUDE:** the +cascade, the dichotomy, and the v₃-hierarchy all live
  in the 3|a frozen marginals, which arg π̂ (via the bridge) cannot touch — so they survive any phase scramble trivially.
- **⚠️ BUT the amplitude is NOT the day's |π̂| work.** The fine |π̂| structure (LAMBDA sup at 2^m, HOMOG diffuse tail) is
  all 3∤a — the same γ-irrelevant part the bridge covers. The channel is carried by a **THIRD object**: the coarse ρ̂
  marginals (tower-frozen), which neither |π̂| amplitude, arg π̂, nor the bridge reaches. Wilson's "if amplitude, the
  |π̂| work comes back into play" does **not** follow — π̂ (amplitude AND phase) is *disconnected* from γ, because the
  bridge only connects it to the fine part. The seam holds.
- **The scramble is moot** by construction: the phase-carrying operator (the bridge) is blind to γ's support.

## Wilson's exact S_∞ split (pen, banked here)
Since the class mean is identical for all 3∤k channels (class-mean theorem), the dominant carrier carries zero
channel-distinguishing info; everything separating channels lives in the cascade `c_k := γ_∞(k) − M_class(k)`. With
Σ_{3|k}4^{−k}=1/63, Σ_{3∤k}4^{−k}=20/63: **S_∞^class = 2[⅔·20/63 + 5/3·1/63] = 10/21 = 0.476190 (= banked S₂)**, and
**S_∞ = 10/21 + 2Σ_k 4^{−k}c_k.** Measured S_∞≈0.4737 ⟹ Σ4^{−k}c_k ≈ **−0.00125** (while c₁ alone gives +0.0167 —
strongly cancelling cascades). **7/15 requires Σ4^{−k}c_k ≈ −0.00476 — ~4× more negative** — the sharpest statement of
what 7/15 demands, now a single number about a single object.

## Net
The gate did its job: the bridge is structurally blind to the channel (τ=0 on 3|a = γ's entire support). γ is
amplitude-carried by the frozen coarse marginals — a third object, disconnected from both the fine |π̂| spectrum and
its phase. This **closes the P1 gate as AMPLITUDE** and, importantly, does NOT resurrect the |π̂| work: P2–P6 (all on
π̂/arg π̂) are aimed at a part the bridge shows is γ-disconnected. Bears directly on P6 (which analyzes arg π̂). Not at
stake: BRIDGE2 (reproduced 3∤a bit-for-bit), P4, CHANNEL_ID, MEAN1, v₃ HIERARCHY, R1–R30. Cheap (0.3s).
