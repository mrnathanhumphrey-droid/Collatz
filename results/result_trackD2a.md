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

## Status
½-flux gate PASS (first symbol factor theorem-grade). Injection tables (factor two's raw material) delivered.
**Held per protocol** — Wilson derives factor two blind from these tables; the assembled σ(θ) vs dense L=2,3
spectra is the C-stage judge, same choreography as τ and the v₁ trit.

Probe `probes/probe_trackD2a.py`; tables `outputs/injection_tables_q3.tsv`.
