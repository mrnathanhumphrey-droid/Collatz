# Probe D2-e — the ladder matrices (raw material for the (★) self-consistency; σ-judge SEALED)

**Date:** 2026-07-18  CPU, dense/direct at q=3 (INSTRUMENT LAW). Probe `probes/probe_trackD2e.py`,
log `logs/probe_trackD2e_log.txt`, matrices `outputs/ladder_matrices_q3.tsv`. **σ-vs-spectra judge stays
SEALED** — this delivers the per-ladder effective matrices and localizes each banked mode; it does not judge σ.

## Construction (ALGEBRAIC, exact tensor structure)
The tower factorizes as **block (e_ρ,γ) ⊗ gauge (dlog a)** — every (e_ρ,γ) block is full-D over the gauge circle
(verified: L=3 tower dim 8424 = 468 blocks × D=18; L=2: 288 = 48 × 6). Gauge-Fourier is therefore `I ⊗ F`
(unitary DFT over the gauge circle), an **independent** construction of the sector operator (parallel to 2c0-G2's
`Bhat = R·N`). Reduction to the effective gauge-frequency operator `A[k_out,k_in]` uses the **E-form (k=0)
block-QSD** (left/right Perron of the k=0 sector block operator) to collapse the (e_ρ,γ) index. Ladder rungs =
`k₀·3^j mod D`, j=0..L−1.

## [3] COUPLING STRUCTURE vs the Decimation Lemma — pre-registered SHAPE
Raw Frobenius coupling `C[k_out,k_in] = ‖block sub-block‖_F` (L=3, 18×18):
- **Corner-dominance ABSENT (CONFIRMED).** No single dominant corner; C is **diagonal-dominant** (survival
  transport on the k_out=k_in diagonal, entries 2.30/2.24/2.09…) **plus a broad background** — the opposite of a
  cyclic companion. This corroborates session-two's negative result (companion route dead) from the *coupling
  matrix*, independently of the 2π/3^{L−1} phase-spacing argument.
- **×3 ladder-step present, not sparse.** Feed *into* the ladder-frequency rungs (multiples of 3) is enhanced:
  `C[:,0]` (into DC bottom rung) rows at k=3,6,9,12,15 = 1.52/1.52/**1.95**/1.52/1.52 vs off-rung ~0.91. So the
  ×3 decimation delta rides on a broad tail rather than being an isolated subdiagonal.
- **Tail feeding the bottom rung (CONFIRMED).** The bottom rungs receive the strongest broad in-feed: k=±1
  bottom (gf=9) row max off-diag **1.95**, k=±2 bottom (DC) row max **1.27** — the tail-fed closure of the
  otherwise-nilpotent shift, as the Decimation Lemma predicts.
- **After QSD reduction the rungs DECOUPLE** (reduced A off-diagonals ~1e-4…1e-3 rel): the effective ladder
  matrix is **near-diagonal**, so the ladder eigenvalues are essentially the diagonal rung multipliers with a
  *small* dressing correction — the worksheet's `λ ≈ σ(θ_k)(1+dressing)` leading balance, made literal.

## [1]/[2] PER-LADDER EFFECTIVE MATRICES + banked-mode locator (L=3)
Reduced 3×3 spectra (mod): **k=±1 → {0.3003, 0.1939, 0.1128}**; **k=±2 → {0.3343, 0.2408, 0.1274}**.

| banked mode | modulus | LOCATED in | value | dist |
|---|---|---|---|---|
| doublet #1  0.23764+0.18303j | 0.3000 | **k=±1 reduced TOP rung (gf=1)** | 0.23747+0.18377j | **0.0008** |
| doublet #2  0.23500+0.18315j | 0.2979 | **k=±1 block-resolved** (rung-1 split) | 0.23523+0.18315j | **0.0002** |
| m=2 seat #1  0.02024+0.18363j | 0.1847 | **k=±1 2nd rung (gf=3)** block-resolved | 0.00330+0.19398j | 0.0199 |
| m=2 seat #2  0.00406+0.19035j | 0.1904 | **k=±1 2nd rung (gf=3)** block-resolved | 0.00330+0.19398j | 0.0037 |

**Findings:**
1. **The doublet IS the k=±1 ladder's top rung (VINDICATED).** The reduced top-rung eigenvalue 0.3003 = the
   doublet centroid (d=0.0008); the block structure **splits** it into the two banked members (d=0.0002, 0.0008).
   Leading balance: 0.3003 = σ(θ₁)·(1+0.020), σ(θ₁)=(1/3)cos²(π/9)=0.2943 — a **+2.0% dressing** on the top rung,
   the small correction (★) supplies. Wilson's "doublet = the k=±1 ladder's two internal top modes" is confirmed
   *operationally*: same rung, split by the block (e_ρ,γ) fabric.
2. **The m=2 seat is a DEEPER RUNG of the SAME k=±1 ladder (refinement of D2-c), not the k=±2 top.** The mod-0.19
   occupants sit at the k=±1 gf=3 rung (reduced 0.00094+0.19385j; block-resolved 0.00330+0.19398j), heavily
   dressed up from σ(θ₃)=(1/3)cos²(π/3)=0.0833 to 0.194. This relabels D2-c's "m=2 seat" as an **internal ladder
   rung**, consistent with the ladder theorem's own internal-rung census (D2-d) — the naive θ_m=2πm/9 seat conflates
   the m-index with gauge frequency; the actual gauge-frequency ladder places the doublet (gf=1) and the m=2 modes
   (gf=3) as successive rungs of the **k=±1** ladder.
3. **The 0.267/0.244 accounting is the k=±2 ladder's solution set (CONFIRMED).** k=±2 block-resolved top carries
   0.2718 / 0.2667 (≈0.267) and 0.2417 / 0.2385 (≈0.244), exactly the modes Wilson's worksheet assigns to (★) on
   the k=±2 ladder. The reduced k=±2 top rung (0.2408) matches the 0.244 branch; the k=±2 mean-field seat (m=2,
   0.196) is **not** its reduced value — **this is where the leakage/dressing is largest and where (★) must do the
   most work** (localized, as deliverable-2 asks).

## Two honest caveats
- **The E-form (k=0) Perron is the finite-L partner, not weight-free 1/3.** 0.3410 (L=2) → 0.3343 (L=3),
  converging to 1/3 — it is the DC/partner rung (the k=0 gauge sector = the E-form), i.e. the finite-section
  partner eigenvalue that the braid tracks, **not** the coarse mean-field S of D0.3 (a different, weight-free
  object). No exactness is broken; the label in the raw log ("DEV from 1/3") is this reframing.
- **The block-resolved subspace is large** (Nb·L = 1404 at L=3), so matches to *small*-modulus modes are density
  coincidences, not captures. The apparent k=±4 (mod 0.010) match into the k=±1 block-resolved pool is one such
  artifact — **discarded**; k=±4 was not in the pre-registered ladder set and is not claimed here. Only the
  **dominant** modes (doublet, 0.267/0.244) are meaningful subspace captures.

## Status
Ladder decomposition **captures the dominant modes**: doublet = k=±1 top rung (both members, d≤0.0008),
m=2 occupants = k=±1 gf=3 rung (d≤0.02), 0.267/0.244 = k=±2 ladder solution set. Coupling shape:
corner-dominance ABSENT (companion route dead, corroborated), tail feeds the bottom rung, ×3 delta rides a broad
background, reduced rungs near-decoupled (small dressing = leading balance λ≈σ(θ_k)(1+dressing), +2.0% on the top
rung). **Largest leakage = the k=±2 top-rung seat**, flagged for (★). Matrices in `outputs/ladder_matrices_q3.tsv`.
**σ-vs-spectra judge SEALED until Wilson's D-2 write-up posts.**
