# Result 77 (sketch): explicit T_lead operator + spectral analysis for rate-1/2

**Status:** sketch / partial — full execution deferred. This document records the operator construction and outlines the rigorous-rate proof path.

## 1. The transfer operator picture

The trajectory measure π on Z_3 admits a **transfer operator** structure inherited from Tao's recursion:
> Syrac(Z/3^{n+1}) = (3·Syrac(Z/3^n) + 1) · 2^{−Geom(2)} mod 3^{n+1}

This induces an operator T on functions on Z_3 (or more precisely on its dual). The level-n stationary distribution π_n is the projection of the inverse-limit fixed point.

For our problem:
- The eigenvalue **λ_1 = 1** corresponds to the stationary direction (Perron mode, gives mean-1 normalization).
- The subdominant **λ_2 = 1/2** (conjectured) governs the rate of S_n → 7/15.

## 2. The leading deviation mode

From R76's leading-mode identity S_{n+1} = −2·M_{n+1}(1+3^n), the rate of S → 7/15 equals the rate of R_n := M_n(1+3^{n−1}) → −7/30.

Define the **deviation sequence**
> δ_n := R_n − (−7/30) = R_n + 7/30

with δ_n → 0 at rate ½ (empirical, R75/R76).

Equivalently, δ_n = −ε_n / 2 where ε_n = S_n − 7/15.

## 3. Mod-3 class decomposition

Split μ̂_n by class of r mod 3:
> μ̂_n^+(ξ) := Σ_{r ≡ 1 mod 3} π_n(r) e^{−2πi r ξ/3^n}
> μ̂_n^−(ξ) := Σ_{r ≡ 2 mod 3} π_n(r) e^{−2πi r ξ/3^n}

So μ̂_n = μ̂_n^+ + μ̂_n^−.

Phase identity (derived from r·ξ mod 3 ∈ {1, 2}):
> μ̂_n(ξ(1 + 3^{n−1})) = ω^ξ · μ̂_n^+(ξ) + ω̄^ξ · μ̂_n^−(ξ)

where ω = e^{2πi/3}, and ω^ξ depends only on ξ mod 3.

## 4. Class-resolved bilinear moments

Define (for each pair (a, b) ∈ {+, −}² and class c ∈ {1, 2}):
> P_n^{a,b}(c) := Σ_{ξ ≡ c mod 3, ξ ∈ (Z/3^n)*} μ̂_n^a(ξ) · μ̂_n^b*(ξ)

This gives 8 quantities at each level (4 type-pairs × 2 ξ-classes), but Hermitian symmetry P^{ba} = (P^{ab})* reduces independent ones to 6 (real(P^{++}), real(P^{−−}), real(P^{+−}), imag(P^{+−}), each at c=1 and c=2; with constraints).

R_n is a specific linear combination:
> R_n = ω̄·[P^{++}(1) + P^{+−}(2) + P^{−+}(1) + P^{−−}(2)] + ω·[P^{++}(2) + P^{+−}(1) + P^{−+}(2) + P^{−−}(1)]

(For real R_n: the ω̄- and ω-coefficients are conjugates, R_n = 2·Re(ω · α_n) with α_n the ω-coefficient.)

Similarly S_n = M_n(1) = Σ_{ξ: 3∤ξ} |μ̂_n(ξ)|² has decomposition:
> S_n = (P^{++}(1) + P^{++}(2)) + (P^{−−}(1) + P^{−−}(2)) + 2 Re(P^{+−}(1) + P^{+−}(2))

## 5. Tao recursion induces T_lead

Tao recursion gives μ̂_{n+1} from μ̂_n. The class-conservation rule (from R66 + the chain dynamics):
- v even → r' ≡ 1 mod 3: μ̂_{n+1}^+ contribution
- v odd → r' ≡ 2 mod 3: μ̂_{n+1}^− contribution

So:
> μ̂_{n+1}^+(ξ) = Σ_{v even, v≥2} 2^{−v} A_v(ξ) μ̂_n(ξ·2^{−v} mod 3^n)
> μ̂_{n+1}^−(ξ) = Σ_{v odd, v≥1} 2^{−v} A_v(ξ) μ̂_n(ξ·2^{−v} mod 3^n)

where A_v(ξ) = e^{−2πi ξ 2^{−v}/3^{n+1}}.

The level-n+1 class-resolved moments P_{n+1}^{ab}(c) are bilinear in μ̂_{n+1}^a, μ̂_{n+1}^b*, hence quadratic in μ̂_n. Substituting Tao's recursion expresses them as quadratic forms in {P_n^{ab}(c)}.

This defines the **transfer operator T** acting on the 6-dim vector P_n = (P_n^{++}(1), P_n^{++}(2), P_n^{−−}(1), P_n^{−−}(2), Re P_n^{+−}(1), Im P_n^{+−}(1), Re P_n^{+−}(2), Im P_n^{+−}(2)) (8 reals, with constraints from total mass).

S_n = linear combination of P_n entries; R_n similarly.

## 6. Spectral gap conjecture

> **Conjecture 77.1:** The transfer operator T has eigenvalue 1 (Perron, captures Plancherel total mass) and subdominant **eigenvalue exactly 1/2**, with eigenvector that projects onto the (R_n − R_∞) deviation.

If T has finite truncations T_N (e.g., level-N restriction), then:
- T_N's eigenvalues converge to T's eigenvalues (Nisoli Theorem 2.15)
- The 1/2 eigenvalue of T_N is computable exactly over Q (since T_N's matrix has rational entries given Tao's recursion + truncated Geom(2))
- Verifying T_N's eigenvalue = 1/2 at small N + Nisoli perturbation bound certifies T's eigenvalue = 1/2

## 7. Implementation outline

1. **Build T_N matrix exactly over Q at level N = 1, 2, 3:**
   - For each pair (a, b) and class c, compute Σ_{ξ ∈ Z/3^N coprime, c-class} μ̂_N^a(ξ) μ̂_N^b*(ξ) symbolically.
   - This gives a 6 × 6 (or 8 × 8 with constraints) matrix at each level.

2. **Find T_N's eigenvalues:**
   - Compute characteristic polynomial over Q.
   - Verify eigenvalue 1/2 appears.

3. **Nisoli perturbation bound:**
   - Estimate ‖T − T_N‖ via Tao Prop 1.17 (super-poly bound on |μ̂_n(ξ)|).
   - Apply Nisoli Lemma 2.9 / 2.12 to get certified eigenvalue bound.

4. **Conclude rate ½ rigorously:**
   - λ_2 = 1/2 ± δ_N → 0 as N → ∞ certifies rate ½.
   - Combined with R75's |ε_n|·2^n empirical bound, gives full certification of c = 7/45.

## 8. Connection to existing results

This program closes the loop on c = 7/45:
- **R74**: Algebraic identity S_{k+1} = 3^{k+1}·‖d_{k+1}‖² (rigorous)
- **R75**: Plancherel decomposition S_k = Σ |μ̂_k|² over 3∤ξ (rigorous)
- **R76**: Conservation law + leading-mode identity (rigorous structural)
- **R77 (this)**: T_lead spectrum gives rate ½ rigorously (pending implementation)

## 9. Why this should work

The key observation is that **rate ½ = P(v=1) for v ~ Geom(2)**. The geometric distribution's mass at v=1 is 1/2, and at each level k → k+1, the "v=1 contribution" perturbs the chain at the rate 1/2 of the previous level's perturbation.

This is the structural origin of rate ½ identified in R71/R73 ("rate ½ matches P(v=1) = 1/2").

The transfer operator T, applied to the deviation subspace, has the dominant eigenvalue 1/2 because the v=1 contribution (which dominates the chain's first-step dynamics) carries weight 1/2 in the Geom(2) distribution.

## 10. Outstanding work

- Implement steps 1–4 above
- Verify spectral structure on T_2, T_3 explicitly
- Apply Nisoli to certify λ_2 = 1/2

Estimated effort: 1-2 hours of focused implementation given the existing infrastructure.
