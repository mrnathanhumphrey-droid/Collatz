# Probe R20 — the thin window — **the whole route collapses to R7; the j≥2 tension is an oscillation; two R19 framings corrected**

**Date:** 2026-07-21  Reuses R7/R9/R10; exact rationals where marked. Probe `probes/probe_window_R20.py`. Per-r
shape statistics (no fit). Certifies the weight identity `Λ_r = Σ_{m≥1} 4^{−m}A_r(m)` (the R7 engine), settles the
m=9 tension, and vacates two of R19's headline framings.

## R20-E — WEIGHT IDENTITY (forced): **PASS** — the deviation-field route *is* the R7 channel engine
`w(k) = 1/(4χ_k(4)−1) = Σ_{m≥1} 4^{−m} e(−mk/3^r)`, so w's orbit-Fourier coefficients are exactly `4^{−m}`, and
**Λ_r = Σ_{m=1}^{3^r} 4^{−m} A_r(m) == Λ_r(R10) == OffDiag_{r+1}/2 exact, r=2…5** (the m=3^r≡0 term carries the DC
`A_r(0)=S_r` at weight `4^{−3^r}`; drop it and the rational equality fails by ~4^{−3^{r−1}}). With
`A_r(m)=γ_r(τ_m)−γ_{r−1}(τ_m)=C_{r+1}(m)/3`: **the deviation-field route, the C-table engine (R7), and the
collision-γ identity (R9) are one computation, not three.** (Third time this campaign walks back to R7 — Wilson's
note; "Re w smooth" in #39 was a rediscovery of the geometric weight.)

## R20-A — THE m≲r WINDOW (measurement, NO fit): the running sum **saturates by m≈3–4**; window coefficients are O(1)
Running `Σ_{m=1}^{M} 4^{−m}A_r(m)` vs Λ_r, all m in the window:

| r | Λ_r | running sum at m=r | gap | terms beyond m=4 |
|---|---|---|---|---|
| 3 | +1.3199e−3 | +1.1526e−3 | −1.7e−4 | (window m≤3) |
| 4 | +6.5026e−4 | +6.8513e−4 | +3.5e−5 | ≤6e−5 |
| 5 | +3.2692e−4 | +3.2132e−4 | −5.6e−6 | ≤3e−5 |
| 6 | −3.3867e−4 | −3.3838e−4 | +2.8e−7 | ≤2e−6 |
| 7 | +2.1485e−4 | +2.1462e−4 | −2.2e−7 | ≤6e−6 |

**The sum is dominated by m=1,2,3; every term with m≥4 is ≤6e−5 and the running sum reaches Λ_r to within ~1e−7 by
r≥6.** The window ratios |A_r(m)|²/typical are O(1) across the whole window (max ≈3.3 at m=3, r=6; no growth trend).
**This is the thin window the proof needs, confirmed:** the crude tail bound `|Σ_{m>M} 4^{−m}A_r(m)| ≤ (4/3)4^{−M}S_r`
is realized — the coefficients that matter are m ≲ r (really m ≲ 4), a logarithmically thin, O(1) set.

## R20-B — m=9 and m=27 SETTLED (measurement): **j≥2 is an oscillation, not exceptional** (R13-C read the oscillation)
Successive differences A_r(m)=γ_r(τ_m)−γ_{r−1}(τ_m):

| r | A_r(9) | sign | | A_r(27) | sign |
|---|---|---|---|---|---|
| 3 | −0.2308 (=−S₃/2 DEF) | − | | +0.4616 | + |
| 4 | +0.0314 | + | | −0.2321 (=−S₄/2 DEF) | − |
| 5 | +0.0269 | + | | −0.0040 | − |
| 6 | +0.0731 | + | | +0.0571 | + |
| 7 | +0.0459 | + | | −0.0482 | − |

**m=9: same-sign but the magnitude turns over at r=6 (0.073 → 0.046) — non-monotone.** R13-C reported
+0.027, +0.073, +0.046 (r=5,6,7): that is the up-swing *and the turnover* — **an oscillation sampled three times,
not exceptional growth.** **m=27: signs alternate (−,+,−)** — plainly oscillating. So the j≥2 coefficients are
bounded/oscillating, not exceptional. And decisively: **m=9, 27 sit outside the window** (weight 4⁻⁹≈4e−6, 4⁻²⁷≈6e−17),
so even the |A|²/typ enhancement (4.8, 5.9 at m=9) contributes ~2e−7 to Λ_r — negligible either way. **R13-C's
decider is answered: no exceptional harmonic; the tension is dissolved.**

## R20-C — ADDITIVE ARGMAX TRAJECTORY (measurement, NO fit): the slow mode **migrates to the trivial character**
| r | argmax ξ | x=ξ/N | dist-to-0 | (2/3)^r |
|---|---|---|---|---|
| 2 | 4 | 0.444 | 0.44444 | 0.44444 |
| 3 | 19 | 0.704 | 0.29630 | 0.29630 |
| 4 | 16 | 0.198 | 0.19753 | 0.19753 |
| 5 | 32 | 0.132 | 0.13169 | 0.13169 |
| 6 | 64 | 0.088 | 0.08779 | 0.08779 |
| 7 | 256 | 0.117 | 0.11706 | 0.05853 |

**The distance from the argmax to the trivial character matches (2/3)^r to all printed digits, r=2…6** (r=3's ξ=19
is not a non-monotone excursion — its conjugate N−19=8=2³ is the near-trivial member, dist 8/27=(2/3)³). **The
additive max is a near-trivial-frequency slow mode migrating to the trivial character geometrically.** r=7 is a
near-degenerate excursion (the argmax jumps to 2⁸=256 while 2⁷=128, dist 0.0585≈(2/3)⁷, is a near-tie in the top-5).
**⚠️ R19-A correction:** "pinned to ⟨2⟩ / 2-adic resonance" is **vacuous and retracted** — 2 is a primitive root mod
3^r, so ⟨2⟩ is the entire unit group and "the argmax is a power of 2" is content-free. The real, measured content
is this: **x → 0 at rate (2/3)^r — the slow mode lives near the trivial character on the *additive* side** (the
coordinate where #38's intuition belonged).

## R20-D — THE 0.233 CHECK (forced): **VACATE** — A_r(3^{r−1}) = −S_r/2 exactly
A_r(3^{r−1}) == −S_r/2 exactly, r=2…7 (the coarsest non-DC stratum is two conjugate members fixed by the DC-mass
split). **True and vacuous** — as R9-D vacated P_n=S_n. The "≈0.233 constant" is just −S_r/2 → **−S_∞/2 = −7/30**,
a restatement of S_∞=7/15, **not a second constant.** **⚠️ R19-D correction:** the "A-side dominant spike at
m*=3^{r−1}, magnitude ≈0.233, *not* 7/30" was reading this artifact — it *is* −S_r/2 → −7/30 exactly. The claimed
"orthogonal ⟨3⟩-orbit spike" is vacated. (R19-D's *conclusion* — the additive max-spike does not touch the b_r
route — still stands, but now on R19-B + R20-A/C, not on an A-side argmax that was definitional.)

## Status
**R20: the route collapses to R7 and the tension dissolves.** **E PASS** — Λ_r = Σ 4^{−m}A_r(m) = OffDiag/2 exact;
deviation-field = C-table engine = collision-γ, one computation. **A** — the running sum saturates by m≈3–4; the
window m≲r coefficients are O(1); the thin-window structure is confirmed. **B** — m=9 (turnover at r=6) and m=27
(sign-alternating) are **oscillations, not exceptional**; R13-C sampled the m=9 oscillation; and both are
weight-suppressed (4⁻⁹, 4⁻²⁷) out of Λ_r regardless. **C** — the additive slow mode migrates to the trivial
character at (2/3)^r (measured, r≤6); "⟨2⟩ resonance" retracted as vacuous. **D** — A_r(3^{r−1})=−S_r/2 vacated;
R19-D's A-side spike was this artifact.

**Consequence for the crux (owed to the pen):** the whole deviation-field program is the **R7 channel engine**, and
the owed estimate is now sharp and thin: a **summable-in-r bound on the O(1) window coefficients A_r(m), m ≲ r** (the
tail m>r is ≤(4/3)4⁻ʳS_r, and the j≥2 modes oscillate + weight-suppress). R19-B and R20-A have measured the
dominant part (m≤4) as O(1) through r=7. **Conditional on that O(1) persisting, Σ_r 4^{−m}A_r(m) converges and Λ_r
is summable — the theorem — for the harmonics that carry weight.** Measured, not proved; five to seven points. No
fitting; exact weight/vacate gates, labeled numeric window/trajectory; two R19 framings corrected as corrections to
my statements, not to R19-B (which stands).
