# Probe D1 — the dark state's name — **walk-back #29: it's a readout zero, not a symmetry**

**Date:** 2026-07-20  CPU, dense L=3 banked vectors. Probe `probes/probe_darkstate_D1.py`, log
`logs/probe_darkstate_D1_log.txt`. Discovery-first, as specified. **Result: no involution commutes with M₊
(P dies as predicted), and the dark state is annihilated by the AGREEMENT FUNCTIONAL — a readout/DC-orthogonality
selection, not a Hamiltonian symmetry.** Wilson's P-conjecture as posed is dead; the true selection rule is named.

## D1-A — commutator census: **NONE commute** (P dies, #29)
‖[M₊,S]‖_F/‖M₊‖_F on the L=3 full operator (dim 8748):

| involution S | ‖[M,S]‖/‖M‖ | verdict |
|---|---|---|
| **P (swap a↔b)** | **1.374** | **NOT a symmetry** — exactly as B1 forced: P is the *intertwiner* M₋=PM₊P, so [M₊,P]=(M₊−M₋)P≠0 |
| G (a,b↦a⁻¹,b⁻¹) | 1.414 | no relation (√2) |
| PG, Gm, PGm, Sig(−a,−b,−γ), Z(γ↦−γ) | ~1.414 | no relation |

**No global commuting involution exists among the natural candidates.** P failing is not a surprise — it is
_forced_ by Theorem 7: since M₊ ≠ M₋ (B1 measured them 0.25 apart), P cannot commute with M₊. **Walk-back #29:
the pair-swap is the +↔− intertwiner, not a symmetry of M₊, so the dark state is not "P-odd."**

## The mechanism (diagnosis): the dark state is a READOUT annihilation ⟨1|r⟩ = 0
The dark doublet member is killed by the **agreement functional** (all-ones readout), not by the init and not by
a symmetry:

| doublet member | λ | **⟨1\|r⟩** (readout) | ⟨ℓ\|v₀⟩ (init) | \|A\| |
|---|---|---|---|---|
| **m0 (dark)** | 0.235+0.183j | **2.7e-14 ≈ 0** | 2.4e-2 (nonzero) | 2.5e-15 |
| m1 (bright) | 0.238+0.183j | 0.175 | 4.7e-4 | 3.2e-4 |

The init **excites** the dark mode (⟨ℓ|v₀⟩=0.024); the **readout cannot see it** (⟨1|r⟩=0). The dark state is a
**mean-zero / DC-free eigenvector** — orthogonal to the uniform all-ones functional. The functional and init
*are* both P-symmetric (verified `<1|P1>`✓, `v₀` P-fixed ✓ — the mechanism's "other half" survives), but since P
is not a symmetry of M, the eigenvectors are **not** P-eigenstates (measured P-"parity" of the doublet = +0.237 /
−0.149, nowhere near ±1). So the parity framing does not apply; the selection is functional-orthogonality.

## D1-C (corrected) — the READOUT-resolved census (the corrected raw material)
The selection rule is **visibility = |⟨1|r_j⟩|** (readout/DC coupling), and it splits the band cleanly. Raw dump,
L=3 (product = A_j(λ_j−1/3), R4-D conventions):

| mode | \|λ\| | phase | \|⟨1\|r⟩\| | channel |
|---|---|---|---|---|
| doublet m0 | 0.298 | 0.662 | 2.7e-14 | **DARK (DC-free)** |
| doublet m1 | 0.300 | 0.656 | 0.175 | visible |
| k2 m0, m1 | 0.185, 0.174 | 1.46, 1.41 | 0.45, 0.54 | visible |
| partner | 0.333 | 0 | 29.6 | visible (couples maximally) |
| band0 | 0.325 | 0.460 | 4.2e-14 | **DARK** |
| band4 | 0.277 | 1.290 | 8.0e-15 | **DARK** |
| band7 | 0.244 | 1.163 | 8.1e-15 | **DARK** |
| band1,2,5,6 | 0.30–0.25 | — | 0.03–0.21 | visible |

**~Half the band is dark (⟨1|r⟩ ~ 1e-14, systematic across the band — not a single fragile mode).** The dark
modes still have nonzero init overlap (⟨ℓ|v₀⟩ up to 0.18) — they are populated but unread. **Constructive note for
the pen:** the selection rule *removes ~half the band* from the agreement observable's shell sum — the naive band
count overcounts the effective density by ~2×, and the dark (DC-free) half is exactly the null column the
corrected edge density must drop.

## D1-D (corrected) — the selection in ladder/readout coordinates
There is no commuting-S action to report (D1-A). The corrected characterization: the agreement functional ⟨1| is
the **DC projector** — ⟨1|r⟩ = Σ_{a,b,γ} r = the k=0 (gauge-DC) component summed over (e_ρ,γ) blocks. **A band
mode is dark iff its DC-gauge/block-sum content vanishes** (pure-AC in the readout). The k=±1 doublet splits:
m1 carries DC leakage through the block coupling (visible), m0 is DC-free (dark). So the visibility lives in the
DC (k=0) sector-sum, not in a symmetry parity. This + R4-D's frozen conventions is the pen's corrected input:
the edge-density integral runs over the **visible (⟨1|r⟩≠0) half** of the band, with the DC-free half projected out.

## Status
**Walk-back #29 logged:** P is the +↔− intertwiner, not a symmetry of M₊; no natural involution commutes. **The
dark state is a readout zero** ⟨1|r⟩=0 (DC-free mode), not a symmetry eigenstate — a functional-orthogonality
selection, with ~half the band dark to the agreement observable. Init excites the dark modes; only the readout
kills them. The corrected raw material is the readout-resolved census (visible ⟨1|r⟩≠0 vs DC-free null column),
frozen for the pen's symbolic Ĝ. No near-EP extraction; deviations reported; the null column is the evidence.
