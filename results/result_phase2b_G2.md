# Result — MICRO-PROBE G2: gating Wilson's fine-structure lemmas. D2 (the v≥1 digit-tape law) is CONFIRMED entirely, to the last edge. U0 (v=0 mean-field) SPLITS: the carry MAP is exactly count-uniform on ℤ/3^{L−1} (Wilson's proof correct for the map), but the OPERATOR's geometric weighting breaks per-state MASS uniformity (dev up to 0.47) — so U0-as-stated is refuted and its "v=0 contributes zero sub-cell correction" consequence fails at the mass/operator level.

**Date:** 2026-07-16. Per-state / per-edge gate of the two limit-worksheet lemmas. Direct/exact at q=3 (INSTRUMENT LAW). Probe `probes/probe_phase2b_G2.py`, log `logs/probe_phase2b_G2_log.txt`. No proof authored, no rate fit; combinatorial + machine-precision.

**Headline: D2 is CONFIRMED to the last edge — t≥2 targets preserve g mod 3 deterministically (139968/139968 at L=3), t=1 from g≡0 states land v′=0 with weight 1 (69984/69984), and t=1 from g≡1,2 splits exactly [1/2,1/2]. The v≥1 sector's digit-tape law is real and exact. U0 SPLITS into a confirmed half and a refuted half: the carry MAP equidistributes γ′ uniformly on ℤ/3^{L−1} in COUNT (edge-count spread = 0 exactly, count v′-split = geometric to 0.00e+00 — Wilson's units-equidistribute proof is correct for the map), BUT the transfer operator weights edges by w_δa·w_δb, so the per-state MASS distribution is NOT uniform (weight spread 0.22/0.26) and the per-state mass v′-split deviates from mean-field by up to 0.27 (L=2) / 0.47 (L=3). U0 as literally stated ("the surviving MASS's new-carry distribution is uniform") is REFUTED, and its architectural consequence ("v=0 sources zero sub-cell correction; all deviation is in v≥1") FAILS in the operator's own measure.**

## LEMMA D2 — the v≥1 fine law (γ = 3g). CONFIRMED (all signatures, both L where testable).
For a v₃(γ)≥1 tower source, surviving moves reach only even targets (gate needs T≡0 mod 3); split by LTE class t(e′)=v₃(2^{e′}−1):
| signature | claim | result L=2 | result L=3 | verdict |
|---|---|---|---|---|
| **D2(a)** t≥2 (e′≡0 mod 6) | γ′ ≡ g mod 3 (digit-shift, deterministic) | 432/432 | **139968/139968** | ✅ deterministic |
| **D2(b)** t=1 (e′≡2,4 mod 6), g≡0 mod 3 | γ′ lands v′=0 with weight 1 | untestable* | **69984/69984** | ✅ weight 1 |
| bonus t=1, g≡1,2 mod 3 | split [v′=0, v′≥1] = [1/2, 1/2] | [0.500000, 0.500000] | [0.500000, 0.500000] | ✅ exact 1/2 |

\* g≡0 mod 3 tower states need v₃(γ)≥2, which first appears at L=3 (γ=9,18).

- **D2(a) and D2(b) are COMBINATORIAL (weight-independent in a,b) and hold to the last edge.** The digit-shift channel (t≥2: the second 3-adic digit of γ promotes to the first of γ′, γ′≡g mod 3) is exact; the g≡0 kill (t=1 from 9|γ states lands entirely on v′=0) is exact.
- **The bonus [1/2,1/2] is w-weighted and comes out exactly 1/2** — the two reachable residues (g+1, g+2 mod 3) receive symmetric mass. The v≥1 cascade IS driven by γ's digit tape, in closed form, exactly as D2 states. **This is where C2's 0.32-deviations come from — now pinned as the LTE-class digit law.**

## LEMMA U0 — v=0 sources "exactly mean-field per state". SPLITS: map ✅ / mass ✗.
| measure | per-state v′-split vs geometric {2/3,2/9,1/9} | γ′ support | within-state spread | verdict |
|---|---|---|---|---|
| **COUNT** (unweighted, the carry MAP) | max dev **0.00e+00** (both L) | 3^{L−1} (3 / 9) | edge-count spread **0** | ✅ **exactly count-uniform** |
| **MASS** (w-weighted, the OPERATOR) | max dev **0.27** (L=2) / **0.47** (L=3) | 3^{L−1} (3 / 9) | weight spread 0.22 / 0.26 | ✗ **REFUTED** |

- **CONFIRMED (map): the carry map equidistributes γ′ uniformly on ℤ/3^{L−1} in counting measure — exactly, per state.** The number of gate-passing (δa,δb) pairs landing on each of the 3^{L−1} carries is identical (spread 0), so the COUNT v′-split is exactly the geometric {2/3,2/9,1/9}. **Wilson's proof (passing units equidistribute, γ′=(γ+um)/3 = uniform h) is correct — for the MAP.**
- **REFUTED (mass, as literally stated): the transfer OPERATOR weights each edge by w_δa·w_δb (geometric).** So the per-state MASS distribution over the 3^{L−1} carries is NOT uniform (weight spread 0.22/0.26), and the per-state MASS v′-split deviates from mean-field by up to 0.47. Worst state (a,b,γ)=(1,4,22) at L=3: mass split {v′=0: 0.401, v′=1: 0.019, v′=2: 0.580} vs mean-field {0.667, 0.222, 0.111} — heavily skewed to v′=2 (the heavy-weight small-δ edges cluster on the deep carries).
- **Consequence for the architecture (flagged): U0's stated consequence — "the v=0 sector contributes ZERO sub-cell correction; all deviation from mean-field is sourced in v≥1 states" — FAILS in the operator's measure.** The v=0 sector's geometric reweighting is itself a per-state sub-cell correction (up to 0.47), even though the underlying carry map is measure-preserving. The eigenfunction-expansion architecture (corrections sourced only in v≥1, U0 making v=0 correction-free) needs restating: v=0 is correction-free in COUNTING measure but not under the w-weighted mass the operator actually propagates. Whether a w-adapted reference measure restores the clean U0-consequence is the open question (Wilson's call — I do not resolve it).

## Adjudication
| lemma | verdict |
|---|---|
| **D2(a) digit-shift** | ✅ CONFIRMED deterministic (139968/139968 @L3). |
| **D2(b) g≡0 kill** | ✅ CONFIRMED weight-1 (69984/69984 @L3). |
| **D2 bonus [1/2,1/2]** | ✅ CONFIRMED exact. |
| **U0 map (count-uniform)** | ✅ CONFIRMED exact (spread 0, dev 0.00e+00) — Wilson's proof correct for the map. |
| **U0 mass (per-state mean-field, as stated)** | ✗ REFUTED (dev 0.27/0.47) — operator w-weighting breaks per-state mass uniformity. |
| **U0 consequence (v=0 correction-free)** | ✗ FAILS in the operator measure — v=0 geometric reweighting IS a sub-cell correction source. |

**⟹ D2 is solid — build the v≥1 fine structure on it (digit-shift + g≡0 kill + [1/2,1/2], all exact). U0's mechanism (carry-map equidistribution) is real and exact in counting measure, but the lemma's "mass uniform / v=0 correction-free" form does not survive the operator's geometric weighting. The contraction estimate (the worksheet's stated crux) cannot assume the v=0 sector is correction-free at the mass level; it must carry the geometric reweighting or work in a measure where the map's exact count-uniformity is the operative structure.**

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, H, J, L-A, L-B, Real-T1, R1-STOP, Probe P/C/E/W/F2/G0/G1. No `r_q` value changes; no rate-law fit (the true 2.9e-3, 1.0e-4 untouched). G1's cell-AGGREGATE cascade {2/3,2/9,1/9} stands — G2 shows it holds per-state only in COUNT, not in w-weighted mass (the aggregate averages the geometric mass skew back to the count value). Real-T1 / PARTNER-CHAR unaffected.

_Reporting discipline: U0 is reported as a genuine SPLIT — the count/map half CONFIRMED to machine precision, the mass half REFUTED as stated, with the architectural consequence's failure named explicitly rather than smoothed. The map-vs-operator distinction (count-uniform carry map vs w-weighted operator mass) is the precise diagnosis, not a hedge. D2's confirmations are reported as deterministic edge-counts (weight-independent), separating them from the w-weighted bonus. Worst-case states are dumped, not just maxima. No "consistent with mean-field" language for the refuted mass version._
