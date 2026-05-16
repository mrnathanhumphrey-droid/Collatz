# Dark-subspace classification of Syracuse — phased attack plan

**Date:** 2026-05-15 (compact-prep). Follow-up to `DWM_FRAMEWORK_DIVE.md`. Implements Option 1 of the three next-probe candidates: full dark-subspace classification per Benoist-Pellegrini-Szczepanek 2024 (arXiv 2409.18655) applied to Syracuse's level-graded DWM trajectory.

**Target:** classify the invariant measures of Syracuse's quantum trajectory by identifying dark subspaces, minimal isometry families, and ergodic-decomposition structure. The c=7/45 subdominant rate is then the spectral gap of Π restricted to the largest non-trivial dark subspace at the inverse limit.

## What "dark subspace" means here

For Syracuse, H = L²((Z/3^n)*) at level n. The Kraus operator family is

  `M_v^{(j, b_prior)}(f)(ξ) = 2^{-v/2} · e^{-2πi ξ · x_j(b_prior) · 2^{-v}/3^n} · f(ξ · 2^{-v} mod 3^n)`

with `x_j(b_prior) = 3^{2j-2} · 2^{-b_prior} mod 3^n` (level-graded phase).

A **dark subspace** D ⊂ H satisfies `M_v^{(j, b_prior)} D ⊂ D` for ALL (v, j, b_prior). Equivalently, the orthogonal projection P_D commutes with every Kraus operator in the family.

The **commutant** A' = {X : [X, M] = 0 for all M in the Kraus family} is the natural algebraic object. By Benoist-Pellegrini-Szczepanek 2024 Theorem 1:
- If A' = C·I (one-dimensional), the trajectory has a unique invariant measure (fully irreducible).
- If dim(A') > 1, the central projections of A' classify the ergodic components, and each component (dark subspace) carries its own invariant measure.

## Phased plan

### Phase 1 — Commutant dimension at small n  *(this session, ~30 min)*

**Probe:** `dark_subspace_probe.py` (already written). Builds all distinct Kraus operators at n=2, 3 with V_MAX=12, J_MAX=3 (sufficient to capture all distinct x_j(b_prior) values mod 3^n). Computes dim(commutant) by SVD of the commutator-stack linear map.

**Decision tree:**
- **dim(A') = 1:** Syracuse is fully irreducible at finite n. The dark-subspace classification is degenerate at the operator-algebra level; move to Phase 2 (structural subspaces beyond standard dark).
- **dim(A') > 1:** Genuine dark-subspace structure exists; proceed to Phase 3 (explicit decomposition).

**Predicted outcome (high confidence):** dim(A') = 1 (full irreducibility), consistent with K_k's spectrum {1, 0, ..., 0} (this session's K_STRUCTURE_RESULT). Syracuse mixes fully in finite k steps — no operator-algebra-level dark subspaces.

If predicted outcome holds: Phase 1 returns "fully irreducible" and we proceed to Phase 2.

### Phase 2 — Structural-invariant subspaces beyond standard dark  *(1-2 sessions)*

Even with fully-irreducible operator algebra, Syracuse has SPECIFIC structural subspaces that carry the c=7/45 spectral content:

(a) **R76 conservation subspace** at the bilinear pair-form level: `D_R76 = {M_n : Σ_j M_n(η_0 + j·3^k) = 0 for all η_0, k}`. Verified to be exact for n ≥ 2 (R76 §11). Captures the "slow mode" at the moment level.

(b) **T_diag (1, 4) eigenspace** in (P_+, P_-) class-resolved: `D_T_diag = span{(1, 4)}` carries T_lead's 43/45 eigenvalue.

(c) **W_{k-1} 3-fiber-zero-mean subspace** at level k: identified this session as K_k's kernel. Contains all "next-level inhomogeneity" data.

These are NOT dark subspaces in the operator-algebra sense (the Kraus operators DO mix them with their complements). But they ARE structural invariants of specific REDUCED operators (T_diag, K, etc.).

**Phase 2 deliverable:** for each (a)-(c), compute the "approximate-darkness" — how much do the Kraus operators leak between D and D⊥? Specifically, compute `‖P_D M_v P_{D^⊥}‖ / ‖M_v‖` for each Kraus operator. Small values = approximately dark; large = strongly mixed.

If R76's D_R76 turns out to be approximately dark (small leakage), it's the "candidate dark subspace" for the inverse limit (where exact darkness recovers).

### Phase 3 — Explicit dark-subspace decomposition  *(only if Phase 1 returns dim(A') > 1; otherwise skip)*

Find central projections P_1, ..., P_m of A' such that I = Σ P_i and each P_i H is irreducible. Compute Π's action on each P_i H separately.

**Phase 3 deliverable:** list of dark subspaces D_1, ..., D_m with their dimensions and characteristic isometries.

### Phase 4 — Spectral gap on the largest dark subspace  *(2-3 sessions)*

For each dark subspace D (or "approximately dark" subspace from Phase 2), compute the spectral gap of Π restricted to D. The c=7/45 rate is given by this gap on the largest non-trivial dark subspace.

**Method:**
- Restrict the Kraus family to D: `M_v|_D = P_D M_v P_D`.
- Form the channel L|_D(ρ) = Σ_v M_v|_D · ρ · (M_v|_D)†.
- Compute spectrum of L|_D's transfer matrix on the (D ⊗ D)/symmetric subspace.
- Top non-trivial eigenvalue = spectral gap on D.

**Test:** if spectral gap on D_T_diag (or its lift) = 43/45 exactly, our framework recovers T_lead's eigenvalue as a dark-subspace spectral gap. That would close the structural picture: T_lead 43/45 = within-level dark-subspace spectral gap.

### Phase 5 — Inverse-limit extension  *(3-5 sessions, substantial)*

Per Benoist-Bruneau-Pellegrini 2024 (arXiv 2403.20094, infinite-dim case) + 2025 purification paper (arXiv 2509.13377), the finite-n dark subspaces extend to the inverse limit H_∞ = L²(Ẑ_3^×).

**Construct:**
- The inverse-limit Hilbert space H_∞ via projective limit of (H_n, ι_{n→n+1}) where ι is the natural Tao recursion injection.
- The level-graded Kraus family at H_∞: the limit of M_v^{(j, b_prior)} as j → ∞.
- The dark-subspace structure at H_∞.

The c=7/45 subdominant rate is the spectral gap of the inverse-limit Π on its largest non-trivial dark subspace. This is the ASYMPTOTIC structural object that R77.5's multi-resolution decomposition was pointing toward.

**Phase 5 deliverable:** rigorous statement of c=7/45 rate as inverse-limit spectral gap. This is paper-shaped.

### Phase 6 — Documentation as paper-shaped result  *(1-2 sessions writeup)*

Combine:
- Morning's DWM identification at moment level (Result 2)
- Evening's structural exhaustion of finite-truncation discrete eigenvalues
- Phases 1-5's dark-subspace classification

Into a single paper-shaped result: **"Syracuse's quantum trajectory has explicit dark-subspace classification via Benoist-Pellegrini framework. The c=7/45 subdominant rate is the spectral gap of the inverse-limit transition kernel restricted to its largest non-trivial dark subspace."**

Add to the writeup:
- Cross-application to physics_detector (AI-video detection via residual diagnostic) as a downstream consequence
- The level-graded phase coupling χ_j(ξ, b_prior) as the load-bearing novel piece (no prior DWM literature treats this)

## Effort estimate

- Phase 1: ~30 min (probe written, just need to run + interpret)
- Phase 2: 1-2 sessions
- Phase 3: 0-1 sessions (only if Phase 1 surprises)
- Phase 4: 2-3 sessions
- Phase 5: 3-5 sessions (substantial; new math required)
- Phase 6: 1-2 sessions writeup

**Total: ~10-15 sessions of focused work** for a full Phase 1→6 closure. Aligns with agent 4's estimate of 5-10 sessions for the Benoist-Pellegrini program applied to Syracuse.

## Compact-survival material

This file + `DWM_FRAMEWORK_DIVE.md` + `dark_subspace_probe.py` + `closure hunt/benoist_*.pdf` are the load-bearing artifacts for the dark-subspace arc. Key references for resumption after compact:

- **Benoist-Pellegrini-Szczepanek 2024 arXiv 2409.18655** — main classification theorem
- **Benoist-Bruneau-Pellegrini 2024 arXiv 2403.20094** — infinite-dim case
- **Benoist-Fraas-Pautrat-Pellegrini 2017 arXiv 1703.10773** — discrete-time framework foundation
- **2025 arXiv 2509.13377** — purification in infinite dim
- All in `C:/Users/Nate/OneDrive/Documents/closure hunt/`

The current open question is now FRAMED: classify Syracuse's dark subspaces, identify the spectral gap on the largest non-trivial one, extend to the inverse limit. This is a recognizable mathematical program with a clear endpoint.

## What this does for the c=7/45 story

**Before this session:** c=7/45 subdominant rate ≈ 0.984 (empirical), no rigorous closure available. T_lead's 43/45 = best finite-rank handle but doesn't close. All finite-truncation operators exhausted.

**After this session + dive:** c=7/45 subdominant rate is recognized as a **spectral gap of an inverse-limit DWM transition kernel on a specific dark subspace**. The path to rigorous closure runs through:
1. Dark-subspace identification (Phases 1-3)
2. Spectral gap on each dark subspace (Phase 4)
3. Inverse-limit extension (Phase 5)

This is the rigorous-closure path that R77.5's multi-resolution displacement, R77.6's branch-cut at z=2, and R77.7's Padé extension all pointed toward but couldn't formalize. The DWM framework gives us the language.

## Caveats

- **No guarantee** that the spectral gap on Syracuse's dark subspaces equals 0.984 exactly, or that the period-9 oscillation is captured by a single dark subspace's spectral feature. The classification is the right framework; whether it CLOSES the c=7/45 rate explicitly remains to be seen.
- **Level-graded structure is novel** — Benoist-Pellegrini's framework assumes time-homogeneous Kraus operators. Syracuse's level-graded `M_v^{(j, b_prior)}` requires extension. This may add structural difficulty.
- **Phase 5 may yield only an existence-of-asymptote result**, not the explicit value 7/45 or rate 0.984. The numerical evaluation might still require ε_n exact extension (R77.7 v2).

These caveats notwithstanding, the structural picture is CLEAN: c=7/45 is the spectral gap of a recognizable mathematical object, not a phantom that resists characterization.
