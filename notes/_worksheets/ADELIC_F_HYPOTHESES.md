# ADELIC_F — Anker-Schapira-Trojan 2013 (heat kernel on affine buildings)

**Source:** C:/tmp/adelic/Anker_Schapira_Trojan_2013_Heat_Kernel_Affine_Buildings.txt. Main result Theorem 3.

## Verbatim Theorem 3 (p. 15, lines ~1139-1162)

> "Theorem 3. Let J ⊊ I_0. Suppose that (ω_n : n ∈ ℕ) is a sequence of co-weights such that V_{ω_n}(o) is contained in the support of p(n; ·). We assume that δ_n = n^{-1} ω_n satisfies
> (18a) lim_{n→∞} n^{-1} dist(δ_n, ∂M)^{-2η} = 0,
> (18b) lim_{n→∞} ⟨δ_n, α⟩ dist(δ_n, ∂M)^{-2η} = 0, for all α ∈ Ψ_+,
> (18c) ⟨δ_n, α⟩ ≥ ξ, for all α ∈ Φ_+ \ Ψ_+,
> for some ξ > 0. Then for any sequence of good vertices (v_n : n ∈ ℕ) such that v_n ∈ V_{ω_n}(o),
> p(n; v_n) = n^{-r/2 - |Ψ_++|} ρ^n e^{-n φ(δ_n)} P_Ψ(ω_n) Q_Ψ(t_n) (1 + E_n(δ_n))
> with |E_n(δ_n)| ≤ C Σ_{α ∈ Ψ_+} (⟨δ_n, α⟩ + n^{-1}) dist(δ_n, ∂M)^{-2η}
> where t_n = (I − T_Ψ) o_n, o_n = ∇ φ(δ_n), and φ(δ) = max{⟨u, δ⟩ − log κ(u) : u ∈ 𝔞}."

## Hypotheses isolated

- **h1 (SPACE):** Affine building Ã_r (or its 1-dimensional case = (q+1)-regular tree). The paper covers root systems Φ of type Ã_r with Weyl group W.
- **h2 (RANDOM WALK):** Nearest-neighbor isotropic random walk on the building, with transition kernel p(n; x, y). Specifically the *W-invariant* random walk associated to a probability measure on co-weights.
- **h3 (CRAMÉR ZONE):** δ_n = n^{-1} ω_n lies in M = int(conv hull V), the Cramér zone (interior of convex hull of vertex set V).
- **h4 (BOUNDARY-DISTANCE CONDITIONS):** (18a)-(18c) regularity of δ_n's approach to ∂M.
- **CONCLUSION:** Sharp asymptotic for p(n; v_n) with explicit exponential decay e^{-n φ(δ)} and polynomial-in-n prefactor n^{-r/2 - |Ψ_++|}. Macdonald spherical function P_ω as combination.

## Hypothesis × input check

| Hyp | (1) μ_n | (4) R78 (1+3)^u |
|---|---|---|
| h1 (affine building Ã_r) | The relevant building for ℚ_3 is Ã_1 = (3+1)-regular tree = T_3 (per Lubotzky 2013, T8 in BT). Same tree object as BT. | — |
| h2 (W-invariant isotropic random walk) | **FAILED.** Tao's recursion is *NOT* a W-invariant isotropic random walk on T_3. Per BT_CANDIDATE_CONSTRUCTIONS Candidate B, Tao's walk preserves ω = ∞, hence is *anisotropic* on T_3 (one direction is "vertical", up toward ω). ASTrojan's framework is for *spherical* (K-bi-invariant) walks; Tao's walk is parabolic-invariant, not K-bi-invariant. | — |
| h3 (Cramér zone) | UNVERIFIABLE for Tao — δ_n = n^{-1} ω_n with ω_n = E[trajectory's tree-coordinate]. Tao stays at *fixed depth n* on the encoded tree (per BT_CANDIDATE_CONSTRUCTIONS Candidate C: "dynamics is lateral motion at fixed depth-n"); doesn't escape to a Cramér point. | — |
| h4 (boundary-distance conditions 18a-c) | N/A given h2 failure | — |

## Disposition for F

**NO_FIT (categorical).**

ASTrojan's framework is for *spherical (K-bi-invariant)* isotropic random walks on affine buildings. Tao's walk is *parabolic-invariant* (preserves ω = ∞ via M_v all sharing this fixed point). These are *different categories* of random walk:

- K-bi-invariant ↔ Macdonald spherical functions ↔ sharp heat-kernel asymptotic (ASTrojan Theorem 3).
- Parabolic-invariant ↔ horocycle-mixing dynamics ↔ no Cramér zone, no spherical-function expansion.

Per BT_PHASE0_THEOREMS T7: "the Tao map is not (visibly) K-bi-invariant" — already noted at BT phase.

**Asymmetry problem (noted at user pre-reg):** the asymmetry between the multiplicative-by-3 step (3-adic shrinking) and divide-by-2 step (2-adic shrinking) is the structural reason Tao isn't W-invariant. The two steps use different primes; W-invariance would require a single-prime symmetric structure.

**Adelic factorization tag:** **NON_ARCH_ONLY** — the affine building is over a single non-archimedean local field. ASTrojan doesn't produce archimedean-place output. Even if Tao's walk had been K-bi-invariant (which it isn't), the conclusion would close only F_3, not F_∞.

## Note on Riemannian-symmetric-space analog

Anker has parallel work on heat kernels on noncompact Riemannian symmetric spaces (e.g., Anker-Ostellari, ref [1] in 2013 paper). For Tao, the natural Riemannian analog would be a heat kernel on (something acting on) ℝ — i.e., a "1-attractor" analysis at the archimedean place ∞. But Tao's measure μ_n has NO archimedean component (it's profinite); the Riemannian heat-kernel framework operates on a smooth manifold, and Tao provides none.

This is the same observation as candidate B's verdict for the archimedean local factor: even though BT says the attractor lives at the archimedean place, μ_n provides no archimedean factor to plug in.
