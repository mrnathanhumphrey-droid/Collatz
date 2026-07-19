# Probe D2-a — the ½-flux gate + injection tables (D-2's factor-two raw material)

**Date:** 2026-07-18  CPU, exact rationals. Probe `probes/probe_trackD2a.py`. Judge data (C) HELD until
Wilson's factor-two derivation posts.

## A — THE ½-FLUX GATE: **PASS** (exactly ½ each, exact rationals, both L)
On the mean-field 4-cell chain's stationary Perron flow (left Perron eigenvector = QSD occupation; T[src][dst],
row-sums are survivals):

| cell | (ev,v0) | (ev,v≥1) | (od,v0) | (od,v≥1) |
|---|---|---|---|---|
| occupation π | 1/3 | 1/6 | 1/3 | 1/6 |
| survival | 2/9 | 5/9 | 5/18 | 4/9 |
| flux π·surv | 2/27 | 5/54 | 5/54 | 2/27 |

- Total surviving flux = **1/3** (= the Perron eigenvalue).
- **v₃=0 (unit-carry) flux = 1/6 → share 1/2**;  **v₃≥1 (divisible) flux = 1/6 → share 1/2**. EXACT.
- Mechanism confirmed: occupation splits by population (2/3, 1/3) and survival compensates (v₃=0 survives
  ~half as often as v₃≥1) so the two fluxes are equal. The QSD occupation π = [1/3,1/6,1/3,1/6] = population
  weights by class (2/3 unit-carry, 1/3 divisible).

**⟹ the symbol's first factor (1+e^{iθ})/2 is theorem-grade** (the ½/½ source split is exact and L-invariant).

## B — INJECTION TABLES (raw material only; W := T̂ + c = γ′ − ⌊γ/3⌋, flux-weighted along the stationary flow)
Delivered to `outputs/injection_tables_q3.tsv` (q=3, L=3, exact rationals). Contents, no interpretation:

**W mod 3 by (source v₃-class, e′ mod 6):**
- v₃=0 (unit-carry, channels e′≡1,3,5): W mod 3 **uniform** {0:1/3, 1:1/3, 2:1/3} in every channel.
- v₃≥1 (divisible): e′≡0 → **W≡0 deterministically**; e′≡2 and e′≡4 → **{W≡1: 1/2, W≡2: 1/2}** (never 0).

**(W mod 3) × (destination v′-class) joint (survival-conditioned injection):**
- W≡0 → {v′=0: 17/24, v′≥1: 7/24}
- W≡1 → {v′=0: 31/48, v′≥1: 17/48}
- W≡2 → {v′=0: 31/48, v′≥1: 17/48}

The W mod 9 tables (finite-L correction order) are in the TSV as well.

## D2-b RIDER — W's top-scale histogram (fork discriminator): **UNIFORM → boundary-sourced**
Flux-weighted along the stationary flow, exact rationals (`probes/probe_trackD2b.py`), L=2,3.
- **W mod 3 marginal = exactly {1/3, 1/3, 1/3}** (uniform), both L — Wilson's pre-registration confirmed.
- **Full W histogram = uniform interior with half-weight endpoints**: L=3 → W=0:1/18, W=1..8: 1/9 each, W=9:1/18;
  L=2 → 1/6, 1/3, 1/3, 1/6. The endpoint (boundary) weight contracts 1/6 → 1/18 (ratio 1/3) from L=2→3. This is
  the finite-section-of-uniform (Lebesgue-on-the-tape) signature.
- **Top-scale leading trit ⌊W/3^{L−2}⌋**: L=3 → {t=0: 5/18, t=1: 1/3, t=2: 1/3, t=3(overflow): 1/18}. NOT the
  {½,½,0} Bernoulli; it is **uniform {⅓,⅓,⅓} with vanishing boundary corrections** (the t=0 depletion and the
  t=3 overflow are the two half-weight endpoints, → 0 as L→∞).

**⟹ Fork verdict: factor two is BOUNDARY-SOURCED (uniform top scale), not a bulk top-scale Bernoulli — the
derivation routes through the section-edge analysis.** Consistent with D-1's H2 (essential-curve finite
sections): the injection is equidistributed on the tape, and the winding/factor-two structure lives at the
truncation edge, not in a bulk coin flip.

## D2-c — three parallel checks (judge still HELD). `probes/probe_trackD2c.py`
**C1 — endpoint contraction: DERIVED, L=4 = 1/54.** The W=0 endpoint atom is the flux of the e′=0 channel (the
only W=0 route: v≥1 sources, T=0), which equals `(1/3D)·Σ_e AC(e) = (1/3D)·(Σw)² = 1/(3D)`, normalized → **1/D
= 1/(2·3^{L−1})**. Sequence 1/6 → 1/18 → **1/54** (×1/3 per level), exact. The Lebesgue-restriction law is not
just extrapolated — it falls out of Σ_e AC(e) = (Σw)² = 1.

**C2 — the m=2 seat EXISTS (claim confirmed; exact values not).** Seats σ(θ)=(1/3)((1+e^{iθ})/2)², θ_m=m·2π/3^{L−1}:
- L=3 m=1 seat (0.2943, 0.698): occupied by the m=1 doublet 0.2979–0.3000 / 0.656–0.662 (phase 0.940 of seat, per the march).
- **L=3 m=2 seat (0.1956, 1.396): OCCUPIED** — nearest pair 0.02024+0.18363j (|·|=0.1847, arg=1.461) and 0.00406+0.19035j (|·|=0.1904, arg=1.550). Modulus ~0.185–0.19 (below the seat 0.196, comparable-corrections as flagged), phase ~1.41–1.55 (seat 1.40).
- **Ordering ✓**: m=2 (|·|≈0.19) is further from 1/3 than m=1 (|·|≈0.30). **Phase ratio ✓**: m=1:m=2 ≈ 0.66:1.41 ≈ **2:1**.
- L=4 m=1 seat (0.3288, 0.2327): **nearly dead-on** — the doublet 0.32895–0.32914 / 0.2306–0.2308 (phase 0.993 of seat).
- L=4 m=2 seat (0.3156, 0.4654): **below D1-C block-6 depth** — not resolved in existing data; needs a deeper block.

**C3 — doublet precision (clean target for the ladder-solution prediction).**
- L=3 (dense, 12 digits): 0.237639959367+0.183030417014j and 0.234998609841+0.183154982890j; **splitting 2.644285e-3**.
- L=4 (D1-C block, ~6 digits): 0.320423+0.075242j and 0.320223+0.075252j; **splitting 2.002e-4**.

## D-2 LADDER JUDGE ITEM — the k=±4 pair (committed before looking): **CONFIRMED**
The ladder theorem identifies the m-index as the coprime class; coprime classes mod 9 = {±1, ±2, ±4}, so
**three pair families at the top of the L=3 spectrum, no more**. Committed seat for k=±4: modulus
(1/3)cos²(2.793/2) = **0.01005**, phase **2.793**. Hunt in the dense L=3 spectrum (`probes/probe_trackD2_k4hunt.py`):

| family | committed seat (mod, phase) | nearest measured pair |
|---|---|---|
| k=±1 | (0.2943, 0.698) | 0.29995 / 0.656 & 0.29794 / 0.662 (the doublet = one ladder's top two internal modes) |
| k=±2 | (0.1956, 1.396) | 0.18474 / 1.461 |
| **k=±4** | **(0.01005, 2.793)** | **−0.00955+0.00287j → 0.00998 / 2.849** |

**The k=±4 pair EXISTS where committed**: modulus 0.00998 vs seat 0.01005 (0.7% — nearly dead-on), phase 2.849
vs seat 2.793 (within the dressing march). Three pair families present, no more at the top. The ladder
theorem's identification (m-index = coprime class, phases 2πk/3^{L−1}) is **vindicated**, and the doublet's
1.009:1 phase ratio is explained as the k=±1 ladder's *two internal* top modes (not two integer rungs).

## D2-d — the k=±4 seat + sector ID + census closure (`probes/probe_trackD2d.py`, dense L=3, committed)
**T1 — seat occupied.** Nearest pair to the committed seat (mod 0.01005, phase 2.793): **−0.009554+0.002875j →
modulus 0.00998, phase 2.849, distance 0.00057.** (The modulus band [0.003, 0.05] holds ~21 pairs; this is the
unique one at the seat.)
**T2 — sector identification CONFIRMED.** The nearest pair's gauge-character mass leads in **k=14 → ±4 class**
(next k=5, also ±4), so the mode is ±4-dominant. D1-A's "k=4 tail" of the leading pair was this mode's shadow —
identification confirmed. (Sector-broad, ~8% in the leading k, as all these modes are.)
**T3 — census closes; NO fourth coprime family (warning shot did not fire).** Among the 354 conjugate pairs
with modulus > 0.003, dominant gauge-k → coprime class:
- **Unit coprime families = {±1, ±2, ±4} — exactly THREE, as predicted (MATCH).** No fourth unit class (and
  none is possible mod 9).
- The remaining dominant-k classes are all **divisible by 3** (r=0, 3, 6 → k ∈ {9}, {3}, {6,15}): 118 pairs =
  the **internal ladder rungs** (frequencies 3k, 9k of the same three ladders), exactly what the ladder theorem
  predicts as the deeper rungs — not a fourth coprime family.
Caveat (honest): the dominant-k grouping is coarse (each class spans a wide modulus range), so this is a
coprime-class **completeness** census, not a clean ladder-top count; the sharp confirmations are T1 (seat
exists) and T2 (seat is ±4). **The ladder theorem's coprime indexing holds — three families + internal rungs +
real modes, no unexplained fourth family.**

## Status
½-flux gate PASS (first symbol factor theorem-grade). Injection tables (factor two's raw material) delivered.
D2-b rider: W top-scale UNIFORM ⟹ boundary-sourced (section-edge route). W mod 3 marginal exactly uniform.
D2-c: C1 endpoint 1/54 derived; C2 m=2 seat exists (right ordering, 2:1 phases); C3 doublets dumped.
D2-d: k=±4 seat occupied (0.00998, 2.849), sector-ID ±4 confirmed, census closes at 3 coprime families (no 4th).
σ-vs-spectra judge SEALED until write-up.
**Held per protocol** — Wilson derives factor two blind from these tables; the assembled σ(θ) vs dense L=2,3
spectra is the C-stage judge, same choreography as τ and the v₁ trit.

Probe `probes/probe_trackD2a.py`; tables `outputs/injection_tables_q3.tsv`.
