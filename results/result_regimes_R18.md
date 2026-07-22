# Probe R18 — regimes and roughness — **A settles the conflation, B corrects "high-frequency", C PASS, D confirms the Prop-1.17 gap**

**Date:** 2026-07-21  Reuses R7/R10; exact where marked. Probe `probes/probe_regimes_R18.py`. Settles the
rate-conflation (ρ ↔ √(1/3)) and measures the roughness / branch structure / max-coefficient profile of the
deviation field δ_r.

## R18-A — EXACT RATIO TABLE (measurement, NO fit): the exact-regime ratio is **~1/2, not √(1/3), and not geometric**
√(1/3) = 0.57735 is the generic √-cancellation **magnitude** (= √⟨|D|²⟩, R17). The exact term-ratios are a
*different* thing and they are **not** near it:

| object | exact-regime \|ratio\| (pre-oscillation run) | then |
|---|---|---|
| d_k = ε_k | 0.535, 0.482, 0.470, 0.432 (k=3→6) — **declining, ≈0.43–0.53** | k=7 → **2.36** (magnitude turns up) |
| Λ_r | 0.181, 0.493, 0.503 (r=3→5) → **≈0.50** | r=6 → **1.036** (sign flip) |
| bulk S_r·b_r | 0.374, 0.492, 0.503 (r=3→5) → **≈0.50** | r=6 → **1.036**, r=7 → 0.634 |

The exact-regime ratio settles toward **1/2** (Λ_r, bulk b_r both hit 0.493/0.503 at r=4,5), **below** √(1/3)=0.577,
and it is **not a clean geometric** — it declines through r=5 and then the **sign oscillation** begins at r=6
(ratio jumps to 1.036, sign flips −). Six exact terms do not complete one oscillation period, so **no term-ratio
rate is establishable** (consistent with the standing "stop measuring the period").

**The conflation, settled — three distinct numbers, three distinct objects:**
- **√(1/3) = 0.577** — the mean **L²-contraction magnitude** ⟨|D|²⟩=1/3 (R17). A per-step *amplitude* of the field.
- **≈0.50** — the exact **term-ratio** of the signed sequence in its short pre-oscillation run (not geometric).
- **ρ ≈ 0.988** — a **2-mode envelope-decay fit** of the *signed* d_k over k=7..13 (`phase_routeB_prime_eps_fit.json`,
  period 37.16, sign-change near k=9). This is the decay of an **oscillation's amplitude envelope**, not a term ratio.

**⚠️ Correction to my own record (banked):** I had used ρ≈0.984 as if it were the same object as √(1/3) (in the
depth steer and the dead-#38 argument). It is not: an L²-contraction magnitude (0.577), a signed term-ratio (~0.5),
and a signed-envelope decay (0.988) are three different quantities. The "6 terms vs ~9-period" framing was measuring
a *signed-envelope* regime with a *magnitude* number — a category error. Corrected here.

## R18-B — ROUGHNESS OF δ (EXACT): δ_r is **broad / equipartitioned across orders — NOT high-frequency**
|ν̂(χ_k)|² binned by v₃(k) (R8 W_j strata), exact via **Σ_{v₃(k)=j}|ν̂|² = Σ_u g_r(u)·c_{3^{r−j}}(u mod 3^{r−j})**
(Ramanujan; cross-checked against N·g(0)−1 = **exact match, all r**). j=0 = finest (order-3^r, roughest);
j=r−1 = coarsest (order-3).

| r | ‖δ‖²(dlog) | frac v₃=0 (finest) | per-stratum fractions (j=0 … r−1) |
|---|---|---|---|
| 4 | 2.0686 | 0.2244 | .2244 .2231 .2302 **.3223** |
| 5 | 2.5342 | 0.1837 | .1837 .1832 .1821 .1879 **.2631** |
| 6 | 3.0003 | 0.1554 | .1554 .1552 .1547 .1538 .1587 **.2222** |
| 7 | 3.4658 | 0.1343 | .1343 .1345 .1343 .1339 .1332 .1374 **.1924** |

**The mass is essentially EQUIPARTITIONED across order-strata** — each v₃=j stratum carries ≈1/r of ‖δ‖², *flat in
j*, with a **mild ~1.3× enhancement at the coarsest (order-3) stratum** (bold), which is the single largest bin at
every r. The finest (v₃=0) fraction declines only as ~1/r — i.e. every *fixed* stratum's fraction → 0 simply because
the *number* of strata grows, not because mass flees to high frequency. **δ_r is a broad field, not a
high-frequency-dominated one; if anything it is mildly coarse-loaded.**

**⚠️ Honest negative on the "better crux" adjective:** the proposed picture "δ_r is a broad, **non-decaying,
high-frequency** field" — the *broad* and *non-decaying* parts hold (‖δ‖²(dlog) grows ≈ +0.466/step, ~linear in r;
consistent with R17-C's broadly-spread, R14's no-non-uniform-ψ). But **"high-frequency" is not supported**: the fine
strata do not grow, the coarse order-3 stratum stays the most-loaded. So the correct roughness statement is
**"b_r small = smooth (low-harmonic) Re w paired against a broad, equipartitioned-across-orders field"**, not "smooth
× high-frequency." The word is corrected; the R17-C broad-field picture stands, now stratified.

## R18-C — BRANCH FACTORIZATION (forced): **PASS** — T = U₊D₊ + U₋D₋ certified
The v-parity split of the R16-A-certified renewal (ord = 2·3^{k−1} is **even**, so v-parity is constant within each
residue class mod ord; the two β-branches are exactly **2^{−v} ≡ ±1 mod 3**, i.e. v even/odd):

- **μ_even + μ_odd == μ_k EXACT**, k=2…5 (byte-identical to `R7.build_mu`).
- **DC split** (ξ=0 branch weights): mass(+, v even) = **1/3**, mass(−, v odd) = **2/3** — exact rationals, all k.
- **Circle-average branch weights** ⟨|D_±|²⟩ = Σ_{v even/odd} 4^{−v} = **1/15** (even) + **4/15** (odd) = **1/3** exact
  (closed form).

The branch factorization **replaces the dead U∘D composition** and is certified as computation. **#39 not incurred.**

## R18-D — MAX-COEFFICIENT PROFILE (measurement): the max **spikes above typical — the Prop-1.17 gap, made concrete**
max_{3∤ξ}|μ̂_r(ξ)| (additive) vs typical √(S_r/(2·3^{r−1})):

| r | max\|μ̂\| | typical | **max/typ** |
|---|---|---|---|
| 2 | 0.3779 | 0.2817 | 1.34 |
| 3 | 0.2522 | 0.1601 | 1.58 |
| 4 | 0.1770 | 0.0927 | 1.91 |
| 5 | 0.1293 | 0.0536 | 2.41 |
| 6 | 0.0961 | 0.0310 | 3.10 |
| 7 | 0.0759 | 0.0179 | **4.25** |

**max/typical GROWS (accelerating): 1.34 → 4.25.** The max coefficient itself decays at rate ~0.70–0.79/step —
**slower than the generic √-cancellation rate √(1/3)=0.577/step** — so the peak (the slow mode, R17) outruns the
typical coefficient by an ever-widening margin. **A triangle-inequality / max-coefficient bound cannot close the
estimate:** it would need max_ξ|μ̂_r| ≲ typical (constant factor), and the data shows the opposite — a growing spike.
This is the concrete Prop-1.17 quantity, measured for the first time: Tao's proposition bounds every coefficient
small **in the level** (r^{−A}); this campaign needs no coefficient larger than typical **in the shell**, and the
max **is** larger than typical, growingly. **Different demands, not a stronger version of the same one** — confirmed
numerically. (The slow-mode rate ~0.75 sits above the mean-L² rate 0.577 of R17; the growing separation is the
localization of the difficulty into the slow mode.)

## R18-E — R85 RUNG-1 FEASIBILITY (one line): **exact DEAD, float CHEAP**
Exact-rational rung-1 at r=8 is **DEAD** two ways: Λ_8 = (ε_9−ε_8)/2 needs ε_9, which is **float-only** (k≥9 wall),
and μ_8's autocorr is ~19M-Fraction (support 4374). Float: **CHEAP** — Bluestein/FFT autocorr O(N log N), N=6561.
**Verdict: only a numeric rung-1 is available, no exact term.** (Stop asking for the exact one — the ε_9 wall
forecloses it.)

## Status
**R18: A settles the conflation** (exact-regime ratio ≈ **1/2**, below √(1/3)=0.577 and not geometric — declines
then oscillates at r=6; ρ≈0.988 is a *signed-envelope* decay of an oscillation, a different object; my prior
0.984-as-√(1/3) usage corrected), **B corrects "high-frequency"** (δ_r is **broad/equipartitioned across orders**,
mildly coarse-loaded — NOT high-frequency-dominated; the *broad, non-decaying* claim holds, the *high-frequency*
adjective does not), **C PASS** (branch factorization T = U₊D₊ + U₋D₋ certified — exact renewal split, DC weights
1/3 & 2/3, circle-avg weights 1/15 & 4/15; #39 not incurred), **D confirms the Prop-1.17 gap** (max/typical grows
1.34→4.25; the slow-mode peak decays slower than √(1/3), so a max-coefficient bound cannot close — the level-bound
r^{−A} is the wrong species for the shell demand), **E exact-dead/float-cheap**.

**The crux, refined (owed to the pen):** the difficulty is **not a norm contraction** — ‖δ_r‖² does not decay (R18-B)
and the max coefficient spikes (R18-D). It is a **low-frequency-decay / regularity statement**: b_r = ⟨δ_r, Re w⟩ is
small because a smooth (low-harmonic) Re w is paired against a broad, equipartitioned field, so the owed estimate is
**summable decay of the fixed-frequency coefficients A_r(m) = C_{r+1}(m)/3** (the C-tables — full circle to R12-B's
"no fifth coordinate"), **with the triangle-inequality/max route confirmed dead**. Same object as R5's open qx+1
step, now sharpened from "tower contraction" to "fixed-m C-table decay with phase cancellation." No fitting; exact
ratio/branch/Ramanujan gates, labeled numeric max/roughness; two honest adjective corrections banked.
