# Result — PROBE 2c1: zeroth bracket + parity-sector tables + scale ledger. (A) Bracket [0.1184, 0.6041] ∋ 1/3 & partner, L-invariant — and its entire width is carried by the ONE (odd,v₃=0) cell (every cell's mean is exactly 1/3). (B) Exact parity-sector R_{D/2}(s) and N_{D/2} tables delivered for Wilson's corrector. (C) HONEST NEGATIVE — the raw defect-mass scale-ledger does NOT cleanly separate q=3 from q=7 (mean/sector ratios overlap: q=3 0.93/0.64 vs q=7 0.72); the marginality reading gets no cheap confirmation from bare ‖D_k‖², it must come from the tax-weighted / selection-graded combination.

**Date:** 2026-07-16. Zeroth Collatz–Wielandt bracket + parity-sector inputs for 2c-2 + scale-ledger marginality discriminator. Direct/entry algebra, no eigen-solves. q=3 (+q=7 control). Probe `probes/probe_phase2c1.py`, log `logs/probe_phase2c1_log.txt`, dumps `outputs/parity_{R,N}_q3_L*.tsv`. Claude gates; ALGEBRAIC/STATISTICAL labeled.

**Headline: (A) the zeroth bracket [0.11836735, 0.60408163] (width 0.486, L-invariant) contains 1/3 and the partner — and the histogram localizes the ENTIRE width to the single (odd,v₃=0) cell (ratios [0.118, 0.604], 6 distinct values); the other three cells are tight ([0.29, 0.43]) and every cell's mean ratio is EXACTLY 1/3. (B) The parity sector k=D/2 is delivered in closed form: exact R_{D/2}(s) tables (L=2: {−13/63, 2/21, −4/63, 0, −4/63, −2/21}; L=3 full 18-entry) and the exact N_{D/2}(e′,γ,γ′) table at L=2 (62 nonzero, values ±1,±3). (C) HONEST NEGATIVE on the marginality discriminator: the raw defect-mass ledger does NOT cleanly separate q=3 from q=7 — count-normalized mean/sector ratios overlap (q=3: 0.934 shallow, 0.642 deep; q=7: 0.718), and q=3's larger deep-scale mass share (26% vs 7%) is mostly the 3-adic sector count (12/4/1) not per-sector intensity. The cheap raw-mass test is INCONCLUSIVE; if marginality is real (and P4/E4/G0 say it is — partner coalesces at q=3, gapped at q=7), it lives in the tax-weighted / selection-graded combination, not the bare ‖D_k‖².**

## A) The zeroth bracket — and where its width lives
h⁰ = coarse Perron right-eigenvector lifted by cell (2/3, 5/3, 5/6, 4/3); operator Mᵀ (flow convention). Bracket = [min, max] of (Mᵀh⁰)/h⁰.
| L | bracket | width | 1/3 in | partner in |
|---|---|---|---|---|
| 2 | [0.11836735, 0.60408163] | 0.4857 | ✅ | ✅ |
| 3 | [0.11836735, 0.60408163] | 0.4857 | ✅ | ✅ |

Ratio histogram by cell (L=3; identical structure at L=2):
| cell | n | [min, max] | mean | # distinct ratios |
|---|---|---|---|---|
| (even, v₃=0) | 2916 | [0.2959, 0.3776] | **1/3** | 3 |
| (even, v₃≥1) | 1296 | [0.2857, 0.4286] | **1/3** | 4 |
| **(odd, v₃=0)** | 2916 | **[0.1184, 0.6041]** | **1/3** | **6** |
| (odd, v₃≥1) | 1296 | [0.2857, 0.4286] | **1/3** | 3 |

- **Every cell's mean ratio is EXACTLY 1/3** — h⁰ makes the mean-field exact at cell level; the bracket width is pure within-cell defect.
- **The entire bracket width is carried by the (odd, v₃=0) cell** (ratios spanning [0.118, 0.604]); the other three cells sit tight in [0.286, 0.429]. This is the SAME cell that is the sole STATISTICAL survival cell (G4 / 2c-0a). **The zeroth defect scale is concentrated in one cell** — a sharp localization for the corrector: 2c-2 need only tame (odd, v₃=0).
- Bracket is L-invariant (identical L=2/L=3) — the width 0.486 is the L-independent baseline every correction order must beat.

## B) Parity-sector closed-form tables (inputs for Wilson's 2c-2 corrector)
Parity sector k = D/2 (χ_{D/2}(2^j) = (−1)^j), the deepest 3-adic scale (v₃(D/2)=L−1).
- **R_{D/2}(s) = Σ_m w_m w_{m−s} (−1)^m, exact (dump `outputs/parity_R_q3_L{2,3}.tsv`):**
  - L=2 (D=6): s=0..5 → **{−13/63, 2/21, −4/63, 0, −4/63, −2/21}** (note R(3)=R(D/2)=0; antisymmetry R(s)·sign vs R(D−s)).
  - L=3 (D=18): full 18-entry table (R(0)=−52429/262143, …, R(9)=0, …), dumped.
- **N_{D/2}(e′,γ,γ′) = Σ_{u:gate,carry→γ′} (−1)^{dlog u}, exact (dump `outputs/parity_N_q3_L2.tsv`):** L=2, **62 nonzero of 64**, integer values ∈ {−3, −1, 1, 3}. These are the twisted unit-counts of the parity block (the k=D/2 case of G2's N_κ).
- The parity block entry is R_{D/2}(s)·N_{D/2}(e′,γ,γ′)/D (the D/2 diagonal sector, G2 with k_in=k_out=D/2). ALGEBRAIC, exact.

## C) Scale ledger — the marginality discriminator. INCONCLUSIVE (honest negative).
Defect mass ‖D_k‖² bucketed by scale a = v_q(k):
| q, L | a=0 | a=1 | a=2 | mean/sector ratios |
|---|---|---|---|---|
| q=3, L=2 | 4 sec, 86.2%, mean 1.606 | 1 sec, 13.8%, mean 1.032 | — | 0.642 |
| q=3, L=3 | 12 sec, 73.5%, mean 17.00 | 4 sec, 22.9%, mean 15.88 | 1 sec, 3.7%, mean 10.20 | 0.934, 0.642 |
| **q=7, L=2** | 18 sec, 92.6%, mean 15.47 | 2 sec, 7.4%, mean 11.11 | — | 0.718 |

- **The raw defect-mass ledger does NOT cleanly discriminate.** Pre-registered "q=3 flat, q=7 decaying" is NOT borne out: q=3's count-normalized mean/sector ratios (0.934 shallow, **0.642** deep) OVERLAP q=7's (0.718) — indeed q=3's deep ratio (0.642) is *steeper* than q=7's (0.718). Both show mild, comparable decay.
- **q=3 does spread more mass into deep scales (26% at a≥1 vs q=7's 7%), but that is mostly the sector COUNT structure** — q=3 (3-adic) has 12/4/1 sectors per scale, q=7 (7-adic) has 18/2 — not a higher per-sector intensity. Count-normalized, the two are similar.
- **Verdict: the cheap raw-mass test is INCONCLUSIVE.** It neither confirms nor cleanly kills the marginality reading. Marginality IS established independently (P4/E4/G0: the partner coalesces with c₀ at q=3, stays gapped at q=7), so it is real — but it does NOT show up in the bare ‖D_k‖²-per-scale. **If Wilson's marginality argument is to rest on a scale-ledger, it must be the tax-weighted / selection-graded combination (‖D_k‖² × cascade-tax, graded by the confirmed 3^{L−1−j}|k selection rule), which is his 2c-3 pen — not this bare-mass diagnostic.** This saves the session the discriminator was meant to protect: don't build the marginality case on raw defect mass.
- Finite-L / small-bucket caveats: the deepest q=3 scale (a=2 at L=3) is a single sector (k=9, parity) and the most truncation-sensitive; its lower value (10.20) may be partly a boundary effect. q=7 has only 2 scales at L=2. No extrapolation.

## Adjudication
| part | verdict |
|---|---|
| A bracket | [0.1184, 0.6041] ∋ 1/3 & partner, L-invariant; width entirely in the (odd,v₃=0) cell (all cells mean exactly 1/3). |
| B parity tables | exact R_{D/2}(s) (L=2,3) + N_{D/2} (L=2) delivered (dumps); the corrector's closed-form targets. ALGEBRAIC. |
| C scale ledger | INCONCLUSIVE — raw ‖D_k‖²-per-scale does not separate q=3 (0.93/0.64) from q=7 (0.72); marginality must come from the tax-weighted combination, not bare mass. |

**⟹ A and B hand Wilson exactly what 2c-2 needs: the L-invariant baseline (width 0.486, localized to one cell) and the parity sector's closed-form R/N. C is an honest steer: the raw defect-mass ledger is NOT the marginality diagnostic — q=3 is not flatter than q=7 in bare ‖D_k‖². The marginality (real, per the spectra) is carried by the selection-graded, cascade-taxed combination the 2c-3 pen must assemble, not by amplitude and not by raw per-scale mass. Both routes that "carry no decay" (amplitude 2c-0b, raw mass 2c-1C) are now ruled out as the mechanism; the exact 3^{L−1−j}|k selection rule (2c0) × cascade tax 3^{−j} is the only remaining candidate.**

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, H, J, L-A, L-B, Real-T1, R1-STOP, Probe P/C/E/W/F2/G0/G1/G2, 2c(0+1), 2c0. No `r_q` value changes; no rate-law fit (the bracket contains the partner but is not fitted; the ledger is reported, not fit). Marginality of q=3 (partner→c₀) stands on P4/E4/G0, not on part C.

_Reporting discipline: part C is reported as an HONEST NEGATIVE — the auto-"flat" criterion that wrongly labeled q=7 flat was removed; the actual finding (q=3 and q=7 mean/sector ratios overlap) is stated plainly, and the deep-scale mass-share difference is attributed to sector count, not intensity, rather than being read as confirmation. The marginality reading is not smoothed into "consistent with flat"; it is explicitly declared UNsupported by raw mass, with the real diagnostic (tax-weighted) named as Wilson's pen. The bracket's width-localization to one cell is a new sharpening (not in 2c-01). Parity tables are exact (Fraction / integer). Caveats on finite-L/small-buckets stated._
