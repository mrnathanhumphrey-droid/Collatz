# Probe J — THE JUDGE, unsealed (the campaign's verdict table)

**Date:** 2026-07-18  CPU. Probe `probes/probe_judge_J1.py` (dense L=2,3), `probes/probe_judge_J2.py`
(SpMV L=4). Everything pre-committed. Model = the QSD-reduced gauge-frequency operator (D-2e machinery):
`A[k_out,k_in] = Lp·M̃[·,k_out,·,k_in]·R`, block-QSD from the E-form (k=0) sector. **σ unsealed.**

---

## J1 — COMPLETENESS AT L=2,3 — **SPLIT VERDICT (honest)**

Every tower mode with |λ|>0.05 (L=3: **626 modes**; L=2: 20) assigned to its dominant gauge frequency →
(ladder k₀, rung gf); model = the per-frequency reduced diagonal `A[gf,gf]`.

### ✅ What the pre-registration got right
- **Totality / no orphan FAMILY (CONFIRMED).** All 626 modes assign to gauge frequencies whose ladders are
  **exactly {±1, ±2, ±4} + DC** — no fourth coprime family, none possible mod 9. The gauge-frequency labeling is
  total; the three ladders + internal (div-3) rungs + DC/Nyquist reals exhaust the |λ|>0.05 spectrum. This
  vindicates D2-d's census at the level of *every* mode, not just the top pairs.
- **The JURY modes ARE captured by the rung-diagonal (1e-3–1e-2):**

  | banked jury mode | measured | gf (ladder,rung) | model A[gf,gf] | resid |
  |---|---|---|---|---|
  | partner (→1/3) | 0.33324 | 0 (DC) | 0.33431 | **0.0011** |
  | doublet #1 | 0.29995+0.6563i(ph) | 1 (±1,r0) | 0.23747+0.18377j | **0.0008** |
  | doublet #2 | 0.29794 | 1 (±1,r0) | " | **0.0025** |
  | m=2 rung | 0.20169 / 0.19247 | 3 (±1,r1) | 0.00094+0.19385j | **0.008 / 0.004** |
  | 0.244 branch | 0.24407 / 0.23539 | 2 (±2,r0) | 0.09305+0.22214j | **0.004 / 0.006** |

### ❌ Where the STRONG pre-registration FAILS (reported straight)
- **"Residuals at the block-structure scale, not the mode scale" does NOT hold for the subdominant modes.** The
  rung-diagonal `A[gf,gf]` is a **single number per frequency**, but each gauge frequency hosts *many* modes
  (gf=1: 21 modes from 0.30 down to ~0.05; gf=9: 55 modes). Residual-to-diagonal reaches **0.484** — comparable
  to the mode modulus itself for the low modes. 20 modes exceed 25%-of-|λ| residual.
- **Gauge frequency is NOT a sharp quantum number.** Purity (mass in the dominant gf) is **0.07–0.32** for almost
  all modes — the modes are broadly spread across the gauge circle. The "dominant-gf" assignment is a *labeling*,
  not a projection onto an invariant subspace.
- **Structural reason (not a bug):** the sector form couples k_in→k_out via `N_{k_in−k_out}` across the *whole*
  gauge circle (2c0-G2), with only a graded/triangular selection rule (2c0-G3: aggregated N vanishes unless
  3^{L−1−j}|k). So the gauge-frequency ladders are a **triangular (nested) grading, not a block-diagonalization**.
  The rung-diagonal predicts the *top* mode of each rung (the jury); the ~600 subdominant modes are O(1)
  block-hybridized and the diagonal cannot reach them.

### J1 verdict
**The ladder decomposition is a COMPLETE LABELING and captures the JURY, but is NOT a complete spectral
reduction.** Completeness holds as "every mode has a ladder/DC home, no fourth family" (strong, confirmed);
it FAILS as "every mode = a lightly-dressed rung eigenvalue" (the diagonal only predicts the dominant mode of
each rung; subdominant modes are heavily block-mixed, gauge freq not sharp). The symbol σ's predictions live on
the **dominant modes**, which are confirmed at 1e-3–1e-2; the full-spectrum reduction would need the block
structure carried explicitly (most of the tower), not the ½-page symbol. Dump: `outputs/judge_completeness_L23.tsv`.

---

## J2 — THE L=4 ARM (SpMV, no dense eig) — **PASS ×3, pre-committed directions**
L=4 tower **233 280 = 4320 blocks × D=54**, nnz 225M, built in 123s (cached `~/j2_L4`). Block-QSD Perron
ρ=0.33345094 (vs g4 partner ρ₄=0.33350, agree 5e-5). Reduced k=±1 ladder (rungs [1,3,9,27]) top-rung
eigenvalue = **d(1,4)**, all via SpMV (no dense eig). σ(θ₁)=(1/3)cos²(π/27)=0.328841, θ₁=2π/27=0.232711.

| quantity | value |
|---|---|
| **d(1,4)** | 0.320294 + 0.075340j (\|d\|=0.329035, arg=0.231023) |
| d(2,4) (k=±2 dominant non-DC) | 0.284822 + 0.138760j (\|d\|=0.316825) |

**Pre-registered directions — all PASS:**
- **(i) |d(1,4)|/σ(θ₁) = 1.00059 ∈ (1, 1.0201) — PASS.** Dressing collapses **+2.01% (L=3) → +0.059% (L=4)**,
  approaching 1 **from above**, exactly the ledger's law. The doublet converges onto the symbol σ(θ₁).
- **(ii) phase ratio arg(d₁)/θ₁ = 0.99275 ∈ (0.9434, 1) — PASS.** Marches **0.940 (L=3) → 0.9928 (L=4)** toward
  1 from below, matching the banked block-level 0.993 to 4 digits.
- **(iii) off-diagonal couplings = 3.6e-5 (k=±1), 4.7e-5 (k=±2) — PASS.** Far below the L=3 scale (~1e-3);
  diagonal dominance not only persists but **strengthens** with L. The rungs decouple faster as L grows.

**Independent cross-check:** the SpMV reduced-diagonal d(1,4)=0.320294+0.075340j reproduces the banked L=4
doublet (D1-C block-6 Rayleigh-Ritz: 0.320423+0.075242j) to **1.6e-4** — the block-splitting scale. Two
disjoint instruments (dense block subspace iteration vs SpMV gauge-Fourier collapse) agree on the L=4 doublet.

### RIDER — 12-digit L=4 doublet split (block-splitting ledger) — **LANDED, 12 digits**
Within-block subspace iteration (block-28 Rayleigh-Ritz, INSTRUMENT-LEGAL), **converged at it 259, Δtop=7.5e-14**
— full 12-digit resolution (the projected operator diagonalizes exactly within the subspace, beating the naive
matvec-floor estimate). Partner **0.333499901322** (= g4 ρ₄=0.333499901324 to **12 digits**).

| L=4 doublet member | value |
|---|---|
| p₀ | 0.320422712770 + 0.075242317692j |
| p₁ | 0.320222549235 + 0.075251807019j |
| **splitting \|p₀−p₁\|** | **2.003883e-4** |

**Block-splitting ledger:** splitting contracts L=3 → L=4 as 2.644285e-3 → 2.003883e-4 = **ratio 0.0758 ≈ ×0.076**
— the pre-registered doublet-contraction ratio, confirmed to 3 digits. The subspace iteration also resolved the
k=±2 family (0.290584+0.144542j, |·|=0.32455) and d(2,4)=0.284976+0.138565j (matching J2's 0.284822+0.138760j).
Dump `outputs/judge_L4_doublet_split.json`. **§6 block-splitting entry now complete.**

