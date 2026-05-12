# Phase 2: Candidate Syracuse-derived δ_a constructions

For each candidate, the IFS is φ_a(x) = (x + a + δ_a(x))/3 for a ∈ {0,1,2}. Verification proceeds against:
- C² on [0,1]
- Uniform contraction: |φ'_a(x)| = |1 + δ'_a(x)|/3 ∈ (0, 1) with sup < 1
- UNI (Claim 2.2 / cond. (10)): ∃ codings ξ, ζ with c < |d/dx(log f'_ξ − log f'_ζ)(x)| ≤ m'
- **Derivation fidelity (A2):** δ_a's formula explicitly uses Syracuse data (parity / v_2 / Tao mod-3 / T_lead eigenstructure), not just "perturbation labeled Syracuse"

**A note on what "Syracuse data" means.** Syracuse-native data:
- The parity sequence (a_n) ∈ {0,1}^ℕ of iterates n_k → n_{k+1} = T(n_k)
- The 2-adic valuation v_2(3n+1) of the odd-step's denominator
- The residue-mod-3 transition kernel of x ↦ Syr(x) on Z_3
- The T_lead = (1/45)·[[7,9],[28,36]] action on V_M ⊆ ℂ² (class-mass deviation)

Anything else dressed up to look "Syracuse" is ad hoc.

---

## Candidate (C1): 2-adic valuation perturbation

**Construction attempt.** Let h: ℝ → ℝ be a C² interpolation of v_2 restricted to positive integers. Set

> δ_a(x) := ε · h(3·(3x + a) + 1)

where 3·(3x+a)+1 is the Syracuse-step numerator applied to digit a in front of x. For x ∈ [0,1] and a ∈ {0,1,2}, the argument ranges over [3a+1, 9+3a+1] ⊆ [1, 19].

**Derivation-fidelity check (A2).** v_2(3n+1) is genuinely Syracuse-native data — it's the geometric-step exponent in Syr(n) = (3n+1)/2^{v_2(3n+1)}. PASS.

**C² check.** **FAIL.** v_2: ℤ_+ → ℤ_+ is integer-valued; on the reals, v_2 has no canonical smooth extension. Possible smoothings:
- Mahler/Volkenborn interpolation: continuous on ℤ_2 but **not on ℝ**. v_2 lives on ℤ_2; embedding ℤ_2 ↪ ℝ via base-2 expansion produces a measurable, not continuous, function. (Lemma: a function on ℤ_2 continuous in the 2-adic topology is generically discontinuous in the archimedean topology pulled back to [0,1] via base-2 expansion.)
- Empirical smoothing: h(y) = E[v_2 | y ± δ] for some real averaging window. This is C² by construction but the smoothing parameter δ is not Syracuse-derived.
- Heuristic interpolation: h(y) = log(1 + |sin(πy/2)|^{−1})/log 2 or similar. C² off integer points, blows up at integers — NOT globally C² on [0,1].

**Verdict on C1: derivation is Syracuse-native, but there is no canonical C² interpolation. Any C² choice introduces a non-Syracuse smoothing parameter, breaking A2 — this is precisely Probe 3 Candidate (b)'s ad-hoc-bump failure mode in disguise.**

---

## Candidate (C2): Cross-frequency perturbation from T_lead eigenvector

**Construction attempt.** T_lead = (1/45)·[[7,9],[28,36]] has eigenvector (1, 4) on V_M ⊆ ℂ²; the (1, 4) ratio reflects asymptotic class-mass ratio P_−/P_+ = 4 between Syracuse parity classes (R76, R77).

Let τ(x) be the class-indicator function on Z_3: τ(x) = +1 if x lies in the "+ class" (v_2(3x+1) even), τ(x) = −1 if "− class". Set

> δ_a(x) := ε · g(τ(κ^{-1}(x))) · χ(a)

where κ: ℤ_3 → [0,1] is the base-3-expansion map, g: {+1, −1} → ℝ is constant-valued (so essentially δ_a is a step function on K_Φ), and χ(a) ∈ {1, 4} reflects the (1, 4) eigenvector.

**Derivation-fidelity check (A2).** PASS — uses T_lead's (1, 4) eigenvector directly.

**C² check.** **FAIL.** τ ◦ κ^{-1}: [0,1] → {+1, −1} is a measurable indicator function on a Cantor-like partition of [0,1]. It is **not even continuous**, let alone C².

Smoothing τ via Gaussian convolution τ_σ = τ ∗ G_σ produces a C^∞ function — but the smoothing parameter σ has no canonical Syracuse value. Worse: as σ → 0, τ_σ → τ pointwise but τ_σ → constant (= 0 or = average) **as iteration of the IFS folds K_Φ densely into itself**. The smoothed function loses its Syracuse-derived structure on the attractor in the small-σ limit.

**Verdict on C2: Syracuse-derived in source but the natural function is discontinuous; smoothing reinstates the ad-hoc parameter problem.**

---

## Candidate (C3): Tao recursion residue-density perturbation

**Construction attempt.** Tao's recursion has a transition kernel on residues mod 3: in Z_3 the Syracuse map x ↦ (3x+1)/2^{v_2(3x+1)} permutes the residue classes via 3x + 1 (mod 9) determining v_2 mod 2, etc. Let k: Z/3 → ℝ encode the stationary distribution of this kernel; let κ: [0,1] → Z/3 send x to floor(3x) mod 3.

Set

> δ_a(x) := ε · k(floor(3x) mod 3) · k(a)

**Derivation-fidelity check (A2).** PASS in source, but...

**C² check.** **FAIL.** floor(3x) mod 3 is a step function on [0,1] with jumps at x = 1/3 and x = 2/3. Three problems:

1. δ_a is C⁰ off the jumps, **discontinuous** at x = 1/3 and 2/3 (unless k is constant on Z/3 — in which case δ_a is constant in x and the IFS is affine).
2. The Tao recursion's stationary measure on Z_3 in Haar measure is **uniform** by the (3, 2) coprimality and dynamics of the Syracuse map (this is essentially the standard result that the Syracuse iteration is ergodic with Haar invariant measure on Z_3). So k(x mod 3) = 1/3 constant — gives δ_a constant in x — back to affine.
3. Even if we replace floor(3x) mod 3 with a smoothed version, the canonical "uniform stationary" answer kills the x-dependence.

**Verdict on C3: PINCER. Either (a) k is non-constant (then δ_a is discontinuous, fails C²) or (b) k = 1/3 (true stationary on Z_3, then δ_a is constant in x, fails UNI by giving an affine IFS). No middle path that is both Syracuse-canonical AND non-vacuous.**

---

## Candidate (C4): Parity-branch density smoothed via Gaussian

**Construction attempt.** Joint density of (x, branch parity) under Syracuse iteration: for x ∈ Z_3 in Haar, the parity of v_2(3x+1) is a Bernoulli(1/2) variable in the limit (the v_2 distribution is geometric P(v_2 = k) = 2^{-k} so parity is 50/50). For x lifted to [0,1] via κ, joint density factorizes: ρ(x, parity) = 1_{[0,1]}(x) · 1/2.

Set

> δ_a(x) := ε · ∫ G_σ(x − y) · b(y) dy,   b(y) = parity-indicator function

**Derivation-fidelity check (A2).** Source is Syracuse but smoothing σ is again non-canonical.

**C² check.** PASS — Gaussian convolution gives C^∞.

**Uniform contraction check.** |φ'_a(x)| = |1 + δ'_a(x)|/3 = |1 + ε · b'_σ(x)|/3 where b'_σ is the smoothed parity derivative. Bounded for small ε. PASS for small ε.

**UNI check.** This is the load-bearing check. The cocycle for two codings ξ, ζ:

  d/dx log φ'_{ξ|n}(x) = d/dx Σ_k log|1 + δ'_{a_k}(x_k)|/3

where x_k = φ_{a_{k-1}} ◦ ⋯ ◦ φ_{a_0}(x). Since δ'_a(x) = ε · (G_σ' ∗ b)(x) is **the same function** for each a ∈ {0,1,2} (the smoothed parity derivative does not depend on the digit a in any Syracuse-canonical way — parity is a property of v_2, not of the digit a we just emitted), the per-step derivative log|1 + δ'_a(x)| is **independent of a**. 

For two σ-periodic codings ξ = (a, a, a, ...) and ζ = (b, b, b, ...) (a ≠ b), the limit derivative d/dx log f'_{ξ|n}(x) converges to a value that depends only on the orbit limit point of φ_a iteration (a fixed point x_a*), and similarly for ζ. Since δ'_a doesn't depend on a, the two limits **coincide if the fixed points x_a*, x_b* happen to be in the same parity class**. Generically they aren't — so UNI fires marginally.

But the UNI strength m(ε): as ε → 0, δ'_a → 0, log|1 + δ'_a(x)|/3 → log(1/3), and the cocycle becomes constant — m(ε) → 0.

Quantitative estimate: m(ε) = O(ε · sup|δ''_a|) = O(ε/σ²) for Gaussian width σ. UNI strength is linear in ε.

**Verdict on C4: C² and UC pass for small ε. UNI fires non-vacuously but with strength m ∝ ε → 0 as ε → 0. This is exactly H_DELTA_EXISTS_BUT_UNI_DEGENERATE in the limit.**

But more importantly: even when C4 fires UNI for ε > 0, the parity-indicator function b is itself **not derived canonically** — it depends on a choice of which Bernoulli(1/2) realization to use. The "average" choice (where b ≡ 1/2) makes b constant and δ_a constant. A non-constant b means selecting a specific almost-sure realization, which is not a canonical Syracuse object.

---

## Candidate (C5, additional): Class-symmetry-broken IFS via T_lead

A more ambitious construction: use T_lead's structure to define an IFS that's not even of the form (x + a + δ_a(x))/3 but rather two distinct branch maps reflecting the (+/-) class split.

> Φ' = {φ_+(x) = (x + 0)/3 + ε·g_+(x), φ_-(x) = (x + 1)/3 + ε·g_-(x), φ_2(x) = (x + 2)/3}

with g_+, g_- chosen so the (1, 4) class-mass ratio is preserved under iteration of Φ'.

**Derivation-fidelity (A2).** Honest attempt — T_lead's two-class structure motivates a two-branch perturbation.

**C² check.** Need g_+, g_- ∈ C²([0,1]). Constructible.

**UNI check.** If g_+ ≠ g_- as C² functions, then d/dx log φ'_+ ≠ d/dx log φ'_-, UNI fires non-vacuously.

**Canonical-relation check (deferred to Phase 3).** The self-conformal measure ν of Φ' lives on K_{Φ'} ⊆ [0,1]; the candidate "Syracuse" piece is the choice g_+(x) − g_-(x) reflecting the (1, 4) eigenvector. But the question is: does ν have a canonical relation to μ_n on Z_3?

**Verdict on C5: PASSES C² + UC + UNI (potentially). Canonical-relation to μ_n is the open gate — analyzed in Phase 3.**

This is the most promising candidate. But the choice of g_± is still constrained to satisfy:
1. Self-conformal measure ν of Φ' projects to μ_n via base-3 expansion + class-decoration
2. Real-Fourier decay of ν transfers to 3-adic Fourier decay of μ_n

Phase 3 examines whether either is achievable.

---

## Summary table

| Candidate | Syracuse-derived (A2)? | C² | UC | UNI | Notes |
|---|---|---|---|---|---|
| C1: v_2 perturbation | YES | **FAIL** | — | — | No canonical C² interpolation of v_2 on ℝ |
| C2: T_lead (1,4) eigenvector | YES | **FAIL** | — | — | Class indicator τ is discontinuous; smoothing is non-canonical |
| C3: mod-3 residue density | YES (in source) | **FAIL** | — | — | Pincer: non-constant → discontinuous, constant → affine |
| C4: parity-density Gaussian-smoothed | partial | YES | YES (small ε) | marginal, m∝ε | Degenerate UNI limit |
| C5: class-broken two-branch IFS | YES (via T_lead) | YES | YES | non-vacuous | **Canonical-relation to μ_n is the deferred gate (Phase 3)** |

Only C5 advances past the entry gate. C4 advances structurally but with degenerate UNI strength; flag for limit-rate analysis in Phase 4 if C5 fails Phase 3.
