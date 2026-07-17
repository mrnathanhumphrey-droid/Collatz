# Result — PROBE 2c0: THE DEFECT OBJECT IS FROZEN AND GATE-CONFIRMED. D-FORM reconstructs M_tower exactly (G1), the sector form closes to machine precision with the twist pinned as R_{k_in}·N_{k_in−k_out} (G2), and the D-DEPTH selection rule HOLDS — with the index pinned to 3^{L−1−j}|k (G3), an exact off-by-one correction to Wilson's stated 3^{L−j}|k from the carry's /3 level-drop. The fluctuation hierarchy is graded by ord₃(k) via EXACT coset orthogonality, not amplitude.

**Date:** 2026-07-16. Gate of Wilson's 2c-0 (D-FORM sector decomposition + D-DEPTH selection rule). Pure entry algebra, no eigen-solves. Direct/exact at q=3 (INSTRUMENT LAW). Probe `probes/probe_phase2c0.py`, log `logs/probe_phase2c0_log.txt`. Claude gates; per-state claims labeled ALGEBRAIC vs STATISTICAL.

**Headline: the defect freeze is confirmed end-to-end. (G1) The single-state D-FORM kernel — from (a,e,γ), target gauge u=a′, shift s, weight w_m w_{m−s}, gate/carry reading u — rebuilds M_tower with max|diff| = 0.000e+00, nnz identical, both L. (G2) The gauge-Fourier sector form B̂[k_out,k_in] = R_{k_in}(s)·N_{k_in−k_out}(e′,γ,γ′) closes to 6.3e-14 over 10.4M sector entries at L=3 — the twist convention is pinned (R on the input index k_in, N on the difference k_in−k_out), k=0 recovering the gate-confirmed E-FORM. (G3, the crown) The D-DEPTH selection rule HOLDS EXACTLY, but with the index the gate pins to 3^{L−1−j}|k (NOT the stated 3^{L−j}|k): resolving γ′ to mod 3^j excites precisely the k with 3^{L−1−j}|k, forbidden k giving |N_k| ≤ 2.3e-14. The off-by-one is the carry's /3 dropping γ′ one 3-adic level. The grading is by ord₃(k) through exact orthogonality on the principal-unit tower — amplitude decay is not invoked and not needed (|R_k|/R_0 → 1 stands).**

## G1 — RECONSTRUCTION: the D-FORM single-state kernel IS M_tower. EXACT.
Built M_tower independently from the D-FORM parametrization (source (a,e,γ); loop target gauge u∈U and shift s; e′=e+s, m=dlog a − dlog u ∈{1..D}, weight w_{m}w_{m−s}, gate `(γ+u(1−2^{e′}))≡0 mod 3`, carry `γ′=(γ+u(1−2^{e′}))/3`), diffed against `build_M_gen` (which loops δa,δb).
| L | max\|D-FORM − build_M_gen\| | nnz (form / true) | verdict |
|---|---|---|---|
| 2 | **0.000e+00** | 3240 / 3240 | ✅ EXACT |
| 3 | **0.000e+00** | 892296 / 892296 | ✅ EXACT |

- **The (a,e,γ)+(u,s) gauge parametrization is verified against the raw (δa,δb) build to the last bit.** u=a′ is the target gauge, m=dlog a − dlog u the halving depth, and the gate/carry read the TARGET gauge u — all confirmed. The single-state D-FORM is exact (ALGEBRAIC).

## G2 — SECTOR FORM: the gauge-Fourier of M_tower. Twist pinned; E-FORM recovered.
For each source block (e,γ)→(e′,γ′), the D×D gauge submatrix B[a′,a] gauge-Fourier-transformed to B̂[k_out,k_in], compared to the D-FORM sector prediction.
| L | sector entries checked | max residual | closing convention |
|---|---|---|---|
| 2 | 13,824 | **2.7e-15** | ✅ B̂[k_out,k_in]=R_{k_in}(s)·N_{k_in−k_out} |
| 3 | 10,415,952 | **6.3e-14** | ✅ same |

- **The sector form closes to machine precision at both L.** Twist convention pinned (as Real-T1's s=+1 was): `R_k(s)=Σ_m w_m w_{m−s} ω^{km}` carries the INPUT sector index k_in; `N_κ(e′,γ,γ′)=Σ_{u:gate,carry→γ′} ω^{κ·dlog u}` carries the sector DIFFERENCE κ=k_in−k_out. **k_in=k_out=0 gives R_0(s)·N_0 = the gate-confirmed E-FORM (F2-1) verbatim** — the defect algebra IS the kinematic algebra, character-dressed, exactly as claimed. The R_k are Real-T1's twisted family; the N_κ are the twisted unit-counts.

## G3 — D-DEPTH SELECTION RULE (the crown): HOLDS, with the index PINNED to 3^{L−1−j}|k.
Aggregating N_k over γ′ at depth-j resolution (γ′ mod 3^j), the twisted count vanishes on all forbidden k:
| L | j (γ′ mod 3^j) | PINNED excited k (3^{L−1−j}\|k) | observed nonzero k | max\|N_k\| forbidden | Wilson-stated (3^{L−j}\|k) |
|---|---|---|---|---|---|
| 3 | 0 | {0, 9} | **{0, 9}** | 2.3e-14 | {0} |
| 3 | 1 | {0,3,6,9,12,15} | **{0,3,6,9,12,15}** | 2.3e-14 | {0, 9} |
| 3 | 2 | all (0..17) | **all** | 0.0 | {0,3,6,9,12,15} |
| 2 | 0 | {0, 3} | **{0, 3}** | 4.7e-15 | {0} |
| 2 | 1 | all (0..5) | **all** | 0.0 | {0, 3} |

- **The selection rule HOLDS EXACTLY** — observed excited-k = pinned set at every (L,j), forbidden k give machine-precision zeros. **The gate PINS the index to 3^{L−1−j}|k**, an exact off-by-one below Wilson's stated 3^{L−j}|k. Mechanism (confirmed): fixing γ′ mod 3^j fixes u mod 3^{j+1} (the carry's `/3` means γ′-depth j needs u-depth j+1), i.e. a coset of the principal-unit subgroup `U_{j+1}=⟨2^{2·3^j}⟩`; the character χ_k sums to zero over that coset unless χ_k is trivial on U_{j+1} ⟺ **3^{L−(j+1)}|k = 3^{L−1−j}|k**. The stated 3^{L−j}|k missed the single level-drop of the carry — exactly "the part the gate must pin."
- **The hierarchy is graded by ord₃(k), through EXACT orthogonality — not amplitude.** Coarse observables (small j, reading few γ′-digits) excite only the sparse high-order sectors (j=0 → only {0,9} at L=3); finer resolution unlocks more sectors. This is the sparsity the contraction needs: `|R_k|/R_0 → 1` (the flagged obstruction stands — 2c-0b), so the defects carry no amplitude decay, but they are **exactly orthogonal to shallow carry structure**. The contraction candidate — sparse selection × cascade tax 3^{−j} — is now on confirmed ground: amplitudes are not asked to decay.

## G4 — scope number (bookkeeping, closes the G0-2 record tension)
Per-STATE survival (column sums of M_tower): **only cell (odd, v₃=0) is STATISTICAL** — within-cell spread 1/3, half-range **0.1667 ≈ Wilson's 0.17**; the other 3 cells are ALGEBRAIC (spread 0). G0-2's "spread 0" was the class-average. (Same as 2c-0a; folded in here for the record.)

## Adjudication
| gate | verdict |
|---|---|
| G1 reconstruction | ✅ EXACT (0.000e+00, nnz match) — single-state D-FORM kernel = M_tower. |
| G2 sector form | ✅ CLOSES (6.3e-14) — twist pinned R_{k_in}·N_{k_in−k_out}; k=0 = E-FORM. |
| G3 D-DEPTH rule | ✅ HOLDS with index PINNED to **3^{L−1−j}\|k** (off-by-one from stated, carry level-drop); forbidden k = machine zeros. |
| G4 scope | (odd,v₃=0) sole statistical cell, half-range 0.167 ≈ 0.17. |

**⟹ The defect object is frozen and gate-confirmed on every axis: exact position-space kernel (G1), exact sector form with pinned twist (G2), and an exact depth selection rule graded by ord₃(k) (G3). The one correction the gate supplies is the index: 3^{L−1−j}|k, not 3^{L−j}|k — the carry's /3 level-drop, precisely the bookkeeping Wilson flagged as the gate's job. The mechanism is confirmed to be exact coset orthogonality on the principal-unit tower, so the contraction (2c-2/2c-3) rests on sparsity-by-selection × cascade-tax 3^{−j}, with the non-decaying amplitudes (2c-0b) no longer load-bearing. D-FORM is ALGEBRAIC and now gate-confirmed; nothing statistical anywhere on the page except the single (odd,v₃=0) survival cell.**

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, H, J, L-A, L-B, Real-T1, R1-STOP, Probe P/C/E/W/F2/G0/G1/G2, 2c stages 0+1 (result_phase2c_01). No `r_q` value changes; no rate-law fit. E-FORM (F2-1) is recovered as the k=0 sector (consistency). The 2c-3 amplitude obstruction (|c_k|/c₀→1) stands and is now explicitly bypassed by the selection rule, not resolved.

_Reporting discipline: G3 is reported as HOLDS-with-pinned-index, not as a pass of the stated rule — the off-by-one is named as a correction (3^{L−1−j}|k) with the exact mechanism (carry level-drop / U_{j+1} coset), and Wilson's stated 3^{L−j}|k is shown alongside as what was missed. The forbidden-k zeros are reported at machine precision (≤2.3e-14). G2's closing convention is stated explicitly (which index carries R, which carries N) rather than asserting "the formula closes." The amplitude obstruction is reiterated as standing (bypassed, not solved). Algebraic/statistical labels applied per the standing rule._
