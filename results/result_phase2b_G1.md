# Result — MICRO-PROBE G1: the coarse tower cell transfer. The 4×4 is EXACT and L-INVARIANT (byte-identical L=2/3), row-sums match {2/9,5/18,5/9,4/9}, dest-parity is set by source v₃-class, and the cascade is universal geometric {2/3,2/9,1/9}. Perron(4×4)=1/3 EXACTLY (not the partner) — the (parity,v₃) cells resolve the survival/cascade structure but are too coarse for the partner eigenvalue.

**Date:** 2026-07-16. Raw transfer dump (the derivation's judge). Direct/exact at q=3. Probe `probes/probe_phase2b_G1.py`, log `logs/probe_phase2b_G1_log.txt`, dumps `outputs/tower_cell_transfer_q3_L{2,3}.tsv`. No proof authored, no rate fit. Cells over the CLOSED carry tower (γ≠0; no return to γ=0 since γ+T≥q ⟹ γ′≥1 — the L-A/Real-T1 protection).

**Headline: the 4×4 cell transfer (source-uniform within cells, exact rationals) is BYTE-IDENTICAL at L=2 and L=3. Its row-sums are exactly {2/9,5/18,5/9,4/9} (pre-reg, gate-matched). The new content — the cell-to-cell split — factorizes cleanly: (i) destination PARITY is determined entirely by the SOURCE v₃-class (v₃=0→odd, v₃≥1→even; source parity irrelevant); (ii) the destination v₃-CASCADE is the universal geometric law {2/3,2/9,1/9}=2·3⁻⁽ʲ⁺¹⁾, identical for all four source cells. The cell-dependence lives ONLY in the survival magnitude, not the cascade shape. Perron(4×4)=1/3 EXACTLY at both L (spectrum {1/3,1/27,0,0}) — NOT the partner (0.346827/0.333236): the 4-cell lumping is not exact (within-cell spread 0.35), so it collapses to the mean 1/3 and the partner's deviation from 1/3 lives entirely below the 4-cell resolution.**

## (1) The 4×4 cell transfer — exact, L-invariant
`T[c,c′] = (1/|c states|)·Σ_{src∈c}Σ_{dst∈c′} M[dst,src]`; cells `(e_ρ mod 2, v₃(γ)∈{0,≥1})`. **Identical at L=2 and L=3:**

| src \ dst | even,v₃=0 | even,v₃≥1 | odd,v₃=0 | odd,v₃≥1 | row-sum |
|---|---|---|---|---|---|
| **even,v₃=0** | 0 | 0 | 4/27 | 2/27 | **2/9** ✅ |
| **even,v₃≥1** | 10/27 | 5/27 | 0 | 0 | **5/9** ✅ |
| **odd,v₃=0** | 0 | 0 | 5/27 | 5/54 | **5/18** ✅ |
| **odd,v₃≥1** | 8/27 | 4/27 | 0 | 0 | **4/9** ✅ |

- **Row-sums match the pre-registered {2/9, 5/18, 5/9, 4/9} exactly (all four, both L).** (Guaranteed given G0-2's cell survivals; confirmed here as a consistency check.)
- **L-INVARIANCE:** the 4×4 is byte-identical at L=2 and L=3. The coarse survival/parity structure is scale-free; the partner's L-flow (0.346827→0.333236) lives entirely in the finer (e_ρ,γ) structure that this cell-coordinate discards.
- **Destination parity is set by the SOURCE v₃-class:** every v₃(γ)=0 source (rows 1,3 in the parity sense — the even,v₃=0 and odd,v₃=0 rows) sends all mass to ODD destination cells; every v₃(γ)≥1 source sends all mass to EVEN destination cells. **The source's own parity does not affect where the mass goes** — only its v₃-class does. (Mechanism: the carry γ′ shifts the phase by δa−δb, and the gate-passing move-pairs from a v₃=0 vs v₃≥1 source select opposite phase-shift parities.)

## (2) The cascade split — universal geometric, class-independent
Per source cell, the destination v₃(γ′) distribution (flow to v′∈{0,1,≥2}), as a fraction of that cell's survival:

| source cell | v′=0 | v′=1 | v′=2 (L=3) |
|---|---|---|---|
| all four cells | **2/3** | **2/9** | **1/9** |

- **The cascade fraction-of-survival is EXACTLY {2/3, 2/9, 1/9} = 2·3⁻⁽ʲ⁺¹⁾ (tail piled at the top) — identical for ALL four source cells.** At L=2 it is {2/3, 1/3} (no v′≥2 level); at L=3 the ≥1 mass splits geometrically into {v′=1: 2/9, v′=2: 1/9}. The raw exact flows (survival × cascade): e.g. even,v₃=0 → {4/27, 4/81, 2/81}; even,v₃≥1 → {10/27, 10/81, 5/81}.
- **The cell-dependence is entirely in the survival MAGNITUDE (row-sum) and the dest-PARITY, not the cascade SHAPE.** The v₃-cascade is the universal geometric ladder regardless of source cell. (This is the cell-aggregate ratio; C2's finer per-(θ=e mod 3, γ) cascade deviations are a sub-cell effect — at this coarser grouping the cascade shape is exactly uniform.)

## (3) Judge (reported, not in the raw dump)
- **Within-cell spread of the per-state dest-cell flow = 0.35 (both L) ⟹ the 4-cell lumping is NOT exact** — states within a cell have very different dest-cell distributions; the 4×4 is a source-uniform AVERAGE, not a Markov lumping.
- **Perron(4×4) = 1/3 EXACTLY at both L** (spectrum {1/3, 1/27, 0, 0}: trace = 10/27 = 1/3+1/27 ✓, rank 2). **This is NOT the partner.** ρ(M_tower) = 0.346827 (L=2) / 0.333236 (L=3); |Perron − ρ| = 1.35e-2 (L=2) / 9.73e-5 (L=3). The shrinking |Perron−ρ| is just the partner APPROACHING 1/3 (the 4×4 Perron is pinned at 1/3, L-invariant), not the coarse operator tracking the partner.
- **Interpretation:** the (parity, v₃) cell coordinate captures the survival bookkeeping (row-sums), the parity-routing, and the geometric cascade — all exactly and L-invariantly — but it is TOO COARSE to carry the partner eigenvalue. The partner's O(10⁻²→10⁻⁴) offset from 1/3 lives in the sub-cell (e_ρ,γ) structure the 4×4 averages away. This is the exact-rational, coarse-grained analog of E's "uniform compression too crude" and consistent with C (partner needs full (e_ρ,γ)) and W (needs quasi-stationary weight).

## Adjudication
| item | verdict |
|---|---|
| 4×4 transfer | EXACT rationals, **L-invariant** (byte-identical L=2/3); row-sums {2/9,5/18,5/9,4/9} ✅. |
| cell-to-cell split (new content) | dest-parity set by source v₃-class; cascade universal geometric {2/3,2/9,1/9}; cell-dependence only in survival magnitude. |
| Perron / L-trend | Perron(4×4)=1/3 exactly (both L), NOT the partner; 4-cell lumping non-exact (spread 0.35) ⟹ coarse coordinate collapses to the mean. |

**⟹ The raw 4×4 + cascade is delivered for Wilson to derive blind from E-FORM + the corrected {2/9,5/18,4/9,5/9} cell law and judge entry-by-entry. Its structure is clean and L-invariant (parity-routing + universal geometric cascade), but its Perron is exactly 1/3 — the (parity,v₃) cells are the right survival/cascade bookkeeping and the WRONG resolution for the partner. The partner's finite-L flow toward 1/3 is a sub-cell phenomenon (full (e_ρ,γ) + quasi-stationary weight, per C/E/W); G1 pins exactly what the coarse cells DO and DON'T carry.**

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, H, J, L-A, L-B, Real-T1, R1-STOP, Probe P/C/E/W/F2/G0. No `r_q` value changes; no rate-law fit (2.9e-3, 1.0e-4 untouched; the 4×4 Perron is 1/3, explicitly NOT identified with the partner or the rate). E-FORM/G0 unaffected (the transfer entries are aggregates of the same gate-confirmed structure; row-sums reproduce G0-2's cell survivals).

_Reporting discipline: the row-sum match is reported as a consistency check (guaranteed by G0-2), not a new pass. The Perron=1/3 is reported as NOT the partner (the coarse coordinate's limitation), with the shrinking |Perron−ρ| explicitly attributed to the partner approaching 1/3 rather than the 4×4 resolving it — no "consistent with the partner" smoothing. The cascade uniformity is scoped as a cell-AGGREGATE ratio (C2's finer per-class deviation is a sub-cell effect, not contradicted). Within-cell spread (0.35) is reported so the 4×4 is understood as a source-uniform average, not an exact lumping. All entries exact rationals._
