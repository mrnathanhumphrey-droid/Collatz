# Result — PROBE P: the dynamical partner lives in the CARRY TOWER (γ≠0). Recon complete; the carry-tower transfer table is the handoff.

**Date:** 2026-07-16. Recon on the dynamical partner's home in the carry coordinate (no proof, no rate fit). Instrument law honored: dense eig + 2D invariant-subspace restriction for q=3 (near the EP); ARPACK used only for the gapped q=7 control (justified). Probe `probes/probe_phase2b_P.py`, log `logs/probe_phase2b_P_log.txt`, data `outputs/partner_*.tsv`.

**Headline: the partner's HOME is the carry tower (γ≠0) — its RIGHT eigenvector carries zero mass at γ=0 (L=3: 67%/22%/11%/0% over carry levels v=0..3). Its LEFT eigenvector localizes on γ=0 (98.85%→99.98%) — the pre-registered left/right split near the EP, robustly confirmed (subspace extraction), NOT an extraction error, and consistent with Real-T1. The partner mixes twist sectors (FORK b) with coupling concentrated in k≡0 mod D/3. The carry-tower transfer table (P3) is L-stable and is the hand-derivation handoff. q=7 control: same tower home, the difference is RATE not LOCATION.**

## P1 — γ-profile of the partner (the core deliverable)
Partner confirmed distinct from all family c_k (L=2: 2.9e-3 from nearest c_k; L=3: 9.96e-5 from c₀ — the coalescence gap — and 0.15 from the nearest *other* c_k). Extracted via 2D invariant-subspace restriction (EP-robust; a raw dense-eig column collapses onto c₀'s pure-γ=0 support at L=3).

| L | eigenvalue | LEFT profile m(v), v=0..L | RIGHT profile (the HOME) |
|---|---|---|---|
| 2 | 0.346827 | [0.0056, 0.0060, **0.9885**] | [**0.6505, 0.3495**, 0.0000] |
| 3 | 0.333236 | [0.0001, 0.0001, 0.0000, **0.9998**] | [**0.6724, 0.2210, 0.1066**, 0.0000] |

- **HOME = carry tower (γ≠0).** The RIGHT eigenvector (where the mode's mass sits) is tower-graded with **exactly zero mass at γ=0** — consistent with Real-T1's corollary (the partner must carry mass off zero-carry, since the family exhausts the γ=0-supported left eigenvectors).
- **The LEFT eigenvector localizes on γ=0** (98.85% → 99.98% as L grows) — the pre-registered **left/right split**. This is the coalescence signature: the partner's *co-observable* concentrates toward c₀'s γ=0 support as the EP is approached, with an **essential, shrinking non-zero tower tail** (0.56%+0.60% at L=2) that keeps it a distinct non-family mode. **Robustly confirmed via subspace restriction — not an extraction artifact.**
- **Pre-registered "pure-γ=0 ⟹ stop, extraction suspect": addressed.** The pure-γ=0 is LEFT-only and robust; the RIGHT home is γ≠0, so Real-T1 is *upheld*, not contradicted. (A mode is characterized by its right home here; the near-γ=0 left is non-normality/EP localization, and it is not exactly γ=0-supported.)

## P2 — gauge structure: FORK (b), coupled, but sub-structured
Does ℓ_partner factor as ω^{−e_a}·f(e_ρ, γ)? **No exact factorization at any twist k** (best residual 1.9e-2 at L=2, 9.6e-4 at L=3 — both ≫ 1e-9). **FORK (b): the partner mixes twist sectors ⇒ D3 is intrinsically a coupled-tower object.** But the coupling is *not* uniform: the near-factorization sits at **k ≡ 0 mod D/3** (L=3: k∈{0,6,12} all at 9.6e-4; the subgroup 6ℤ/18 ≅ ℤ/3). So the partner lives in a **ℤ/3 twist sub-family**, mixing those three sectors — a sharper target than "all D sectors."

## P3 — carry-tower transfer (THE HANDOFF)
Raw operator flow aggregated by carry level `v = v₃(γ)` (γ=0 ↦ level L). Row-normalized (source level → destination), L-stable:

| src level | → dest levels (L=2) | → dest levels (L=3) |
|---|---|---|
| v=0 | [2/3, 1/3, 0] | [2/3, 2/9, 1/9, 0] |
| v=1 | [2/3, 1/3, 0] | [2/3, 2/9, 1/9, 0] |
| v=2 | — | [2/3, 2/9, 1/9, 0] |
| v=L (γ=0) | [2/3, 0, 1/3] | [2/3, 2/9, 0, 1/9] |

- From every tower level: **2/3 to v=0**, the remainder geometric down the tower (L=3: 2/9, 1/9 — ratios 6:2:1). The γ=0 level feeds the same profile with its self-slot moved to γ=0. **This is the explicit carry-tower recursion whose Perron root the partner should be** — the object next session's hand-derivation opens on. Full tables (raw + partner-weighted): `outputs/partner_transfer_{raw,wt}_q3_L{2,3}.tsv`.

## P4 — q=7 control: same HOME, the q=3 specialness is the RATE
| quantity (q=7, L=2) | value |
|---|---|
| subdominant *eigenvalue* (by modulus) | 0.301675−0.126026j, \|·\|/c₀ = 0.9808 — **a FAMILY member** (dist to c_k = 3e-16), amplitude-invisible |
| top **non-family** tower-partner | 0.158414 (real), \|·\|/c₀ = 0.4752, at **rank 22** by modulus (= R26's L=2 value 0.475; c₀≈1/3 at q=7 too) |
| its γ-profile m(v) v=0..2 | [0.6476, 0.3190, 0.0334] — **tower-graded, same as q=3** |

**Pre-registration UPHELD (after correcting search depth — the top is dense with amplitude-invisible family members up to 0.98·c₀).** The q=7 tower-partner has the **same tower-graded home** as q=3 (67/32/3 ≈ q=3's 67/22/11). **The q=3 specialness is NOT where the mode resides — it is the RATE: at q=3 the partner rises to coalesce with c₀ (gap 1e-4); at q=7 it sits far below (rank 22, 0.475).** Exactly the pre-registered expectation ("specialness in the approach to the family ceiling, not the home").

## Adjudication
| probe | verdict |
|---|---|
| P1 | HOME = carry tower (RIGHT, γ≠0); LEFT localizes γ=0 (split, robust). Real-T1 upheld. |
| P2 | FORK (b): mixes sectors, coupling concentrated in k≡0 mod D/3 (ℤ/3 sub-family). |
| P3 | L-stable carry-tower transfer emitted — the recursion the partner is the Perron root of (handoff). |
| P4 | same tower home at q=7 (rank-22 partner); q-specialness = rate (approach to c₀), not location. |

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, H, J, L-A, L-B, Real-T1, R1's STOP. No `r_q` value changes; **no rate-law fit** (2.9e-3, 1.0e-4 stay a two-point sequence until derived).

_Reporting discipline: the L=3 LEFT "pure-γ=0" pre-registered flag fired and was resolved (robust subspace extraction confirms it is a genuine left/right split, not an extraction error; the RIGHT home is γ≠0, upholding Real-T1). P2's fork is reported as (b) with the sub-structure (k≡0 mod D/3) as a refinement, not a fitted claim. P4's initial "no non-family in top-12" was corrected by searching to rank 30 — the partner exists at q=7 with the same home; reported as pre-reg-upheld, with the honest note that the near-top is family. No rate-law extrapolation._
