# Result — REQUEST F: the e=−1 sub-block is NILPOTENT (ρ₋ = 0). D1's maximality closes with maximal margin.

**Date:** 2026-07-16. Pre-registered probe (`probes/probe_phase2b_F.py`, log `logs/probe_phase2b_F_log.txt`). Closes the one owed piece of **THEOREM D1** (`BRIEF_D1_TOY_GAP.md`): that no off-diagonal `e = ab = −1` mode beats the diagonal-ray subdominant `|λ₂| = (1−λ)/(1+λ)`.

> **⚠️ CORRECTION (2026-07-16, from Request-F-corollary gate `result_phase2b_Dmax.md`):** the "nilpotency index = 2L" claimed below was a **two-point extrapolation** (L=1→2, L=2→4). At L=3 the index is **5, not 6** — the true sequence is **{2, 4, 5}** for L={1,2,3}. The proof's load-bearing bound **`index ≤ 2L` HOLDS** (and `Msub^{2L}=0`), so **ρ₋=0 / nilpotence / maximality are unaffected** — only the "= 2L exactly" flourish is retracted. Read "≤ 2L" wherever "= 2L" appears below.

**Headline: ρ₋ = 0 EXACTLY at all 18 grid points — the e=−1 sub-block is genuinely NILPOTENT (nonzero entries, but no cycle), nilpotency index ≤ 2L. This is stronger than P1 (`ρ₋ < |λ₂|`): the off-diagonal sector cannot sustain ANY mode. Claude's "two-step kill" is confirmed as literal nilpotence — at L=1 the block satisfies Msub² = 0 (mass exits in ≤2 steps).**

## Method (per the guards)
From the frozen `build_M_gen(q, −1, λ)` (frozen weights `λ^δ → [u,s]`), the `e=−1` sector = states `(a,b,γ)` with `b ≡ −a mod q^L`. ρ₋ = spectral radius of the **principal submatrix** `M[S,S]` — which **drops** the `us`-exit branches (flip-stay / stay-flip → `e=+1`) **without renormalizing**, so the block is strictly sub-stochastic by construction. Dense exact eig (`np.linalg.eig`); no ESPRIT, no mass sequences. Full-`M` sanity: `build_M_gen` reproduces D1's `λ₁ = (1+λ²)/(1+λ)²` and `|λ₂| = (1−λ)/(1+λ)` to 6 digits (right object).

## ρ₋ table — grid q ∈ {5,7,13} × L ∈ {1,2} × λ ∈ {0.30, 0.50, 0.70}
**Every one of the 18 points: ρ₋ = 0.000000.**

| λ | \|λ₂\|=(1−λ)/(1+λ) | us=λ/(1+λ)² | u²=1/(1+λ)² | ρ₋ (all q, all L) | P1 |
|---|---|---|---|---|---|
| 0.30 | 0.5385 | 0.1775 | 0.5917 | **0.000000** | PASS |
| 0.50 | 0.3333 | 0.2222 | 0.4444 | **0.000000** | PASS |
| 0.70 | 0.1765 | 0.2422 | 0.3460 | **0.000000** | PASS |

`q`-spread of ρ₋ over {5,7,13} = **0.00e+00** at every (λ,L); `|L2−L1|` = **0.00e+00** at every (λ,q).

## The block is NONZERO but nilpotent (not an empty-extraction artifact)
| q | L | λ | dim | nnz(Msub) | max column-sum | nilpotency index |
|---|---|---|---|---|---|---|
| 5 | 1 | 0.5 | 10 | 4 | 0.4444 = u² | **2** |
| 7 | 1 | 0.5 | 14 | 4 | 0.4444 = u² | **2** |
| 5 | 1 | 0.3 | 10 | 4 | 0.5917 = u² | **2** |
| 7 | 2 | 0.5 | 98 | 28 | 0.4444 = u² | **4** |
| 13 | 2 | 0.7 | 338 | 52 | 0.3460 = u² | **4** |

**The flip-flip candidate is present and carries exactly its `u²` weight** (max column-sum = u² at every λ — the `u² = 4/9 > 1/3` "dangerous" transition Claude flagged is really in the block). It simply **forms no cycle**: the block is a nilpotent (strictly upper-triangularizable) transient operator. Nilpotency index ≤ **2L** (measured {2,4,5} at L={1,2,3} — see correction banner; the ≤2L bound is what the lemma uses).

## Pre-registration adjudication
- **P1 — `ρ₋ < |λ₂|` strictly, everywhere. ✅ PASS (18/18), with the maximal possible margin: ρ₋ = 0.** No off-diagonal mode approaches the ray's `|λ₂|`; D1's committed theorem statement stands.
- **P2 — `ρ₋` at `us`-class or below, not `u²`-class. ✅ CONFIRMED and SHARPENED.** ρ₋ is not merely `us`-scale — it is **exactly 0**. The `u²`-weight flip-flip transitions exist (max column-sum = u²) but are non-recurrent, so the sector's spectral radius vanishes. **Claude's two-step-kill mechanism is CORRECT** and is literally the nilpotency index (2 at L=1). The bounding constant for the proof is therefore trivial: the `e=−1` reachable subgraph is a **DAG** (acyclic) — mass reaches a `us`-exit within 2 steps per carry level.
- **P3 — q-flat and L-stable. ✅ CONFIRMED (q-flat: spread exactly 0).** Spectral radius is L-stable (0 at every L). **Refinement (observation, not a P-deviation):** the *nilpotency index* grows with L (measured {2,4,5} at L={1,2,3}, ≤2L; the earlier "=2L" was a 2-point fit — see banner) while ρ₋ stays 0 — the once-for-all-q form is justified (q-independent structure).

## What this hands the proof
D1's maximality is now a **confirmation, not an exploration**: the finite reachable-graph argument reduces to *"the `e=−1` reachable carry-subgraph is acyclic"* — provable once for all odd q (P3 q-flatness), with the carry gate `γ ± 2a ≡ 0 mod q` supplying the exit in a bounded number of steps (nilpotency index ≤ 2L). **D1 fully closes** on the numerical side; the pen-and-paper acyclicity proof (Nathan) meets a ρ₋ = 0 target, not a `< |λ₂|` inequality. **→ Proof written and gated: `result_phase2b_Dmax.md`.**

## Not at stake
R1–R46, Phases 0/1/2a, Phase 2b Sessions 1–2, D1 (a)–(d). No `r_q` value changes. This closes D1's owed maximality piece numerically.

_Reporting discipline: ρ₋ = 0 is disclosed as STRONGER than the pre-registered `ρ₋ < |λ₂|` (P1), and verified NOT to be an empty-extraction artifact (nnz > 0, max column-sum = u², nilpotency index 2L). P2's "us-class" pre-registration is reported as under-shot (actual = 0), with the mechanism (two-step carry kill) confirmed as the cause. P3's q-flatness is exact (spread 0); the L-dependence of the nilpotency index is flagged as a refinement, not a deviation. Full-M subdominant cross-checked against D1's closed forms._
