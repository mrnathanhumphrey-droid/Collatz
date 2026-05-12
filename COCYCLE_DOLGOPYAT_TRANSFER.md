# Transfer mechanism: from ARHW conclusion to μ_n Fourier decay

ARHW Thm 1.1 conclusion: **|F_q(ν)| = O(|q|^{−α})** where ν is the self-conformal measure on the IFS attractor K ⊆ [0,1] ⊆ ℝ.

Syracuse stationary measure μ_n lives on ℤ_3 (or ℤ/3^n ℤ in finite-dimensional truncations) — a different topological space.

## (a) Relationship between μ_n and ν for the surviving candidate

The only surviving candidate from `COCYCLE_DOLGOPYAT_CANDIDATES.md` is (b): twisted base-3 IFS φ_a(x) = (x + a + δ_a(x))/3 with bolted-on δ_a making it non-conjugate to linear.

Map κ : ℤ_3 → [0,1] via the base-3 digit expansion: an element of ℤ_3 has expansion Σ a_i 3^i with a_i ∈ {0,1,2}; the corresponding point in [0,1] is Σ a_i 3^{−i−1}. This is a continuous surjection, **measure-isomorphic up to a countable boundary set** (the rationals with terminating base-3 expansions). Push-forward κ_∗ μ_n gives a measure on [0,1].

For untwisted base-3 (candidate a), κ_∗ μ_n IS the self-similar measure ν_p with weights from μ_n's marginal on {0,1,2}.

For twisted base-3 (candidate b), the conjugacy between the two IFS attractors is not the identity; ν of the twisted IFS lives on the twisted attractor K_Φ ⊆ [0,1] which differs from [0,1]. The push-forward κ_∗ μ_n need NOT equal ν.

**Conclusion: the relationship κ_∗ μ_n ↔ ν is only canonical for the linear case Φ_3 (a), and that case is precisely the one ARHW EXCLUDES. For the twisted IFS (b), the relation breaks.**

## (b) Does polynomial decay on ν imply decay on μ_n?

Suppose we accept candidate (b) and grant |F_q(ν)| = O(|q|^{−α}). To get decay on μ_n, we need to relate Fourier transforms across the κ map. The Syracuse Markov chain's natural Fourier object is on ℤ_3 — characters of ℤ_3 are χ_q for q ∈ ℚ_3 (in the form q = j/3^k), not characters of ℝ. So "F_q(μ_n)" must be defined p-adically.

The Tao 1909.03562 §7.4 iterated-cubic obstruction (the ultimate target) is stated for 3-adic Fourier analysis on ℤ_3. The polynomial-in-A bound needed for Prop 1.17 is a 3-adic-Fourier-decay statement, not a real-Fourier-decay statement on a measure on [0,1].

There is no general principle saying real-Fourier decay of κ_∗ μ_n implies 3-adic-Fourier decay of μ_n. They are different transforms: real characters x → e^{2πiqx} don't restrict to 3-adic characters under κ.

**The transfer mechanism doesn't exist generically.** Even granting the framework's hypotheses fire, the resulting decay is on the wrong Fourier transform.

This is exactly the Probe 2 T1-transfer failure mode (A3 adversarial check, called out in pre-registration): "Framework gives decay on ν living on IFS attractor in ℝ. μ_n lives on ℤ_3. Different topological spaces. Probe 2 fell on T1 transfer; this one might too." Confirmed.

## (c) Polynomial-in-A bound

The ARHW exponent α(ν) is implicit through the cocycle data: spectral-gap parameters (ε, α, γ) in Theorem 2.8 depend on UNI constants m, m' (Claim 2.2 / Theorem 2.4 part 5), contraction parameters ρ, ρ_min, distortion constant L. The paper does NOT extract a polynomial-in-A constant in the form needed by Tao Prop 1.17 (which requires C · A^{−k} type bounds with explicit k as a function of the iteration count).

Even WERE the framework to apply, the conversion to a polynomial-in-A bound would require an additional certification step (extracting how α(ν) varies with the number of iterations / the size of the cyclic structure A), which is not done in ARHW.

**No polynomial-in-A bound can be named from this framework as it stands.**
