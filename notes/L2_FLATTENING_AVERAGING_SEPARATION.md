# Phase 4 — BKS averaging and separation hypotheses for Syracuse μ_n

**Date:** 2026-05-12. L²-flattening structural-compatibility probe, Phase 4.

---

## 1. Scope of this phase

Phase 3 concluded that the BKS L²-flattening hypothesis itself does not transfer to Syracuse μ_n in any of the three plausible translations. Strictly speaking, Phase 4 (verifying BKS's averaging step and separation hypothesis for Syracuse) is moot — if the framework doesn't apply at the flattening step, the downstream hypotheses are irrelevant.

However, the procedure asks for the averaging and separation checks **specifically** so the disposition documents what fails and where. For completeness:

- §2 documents what the BKS averaging step (Step 1 in the unified strategy) requires, and whether it is satisfied for Syracuse μ_n.
- §3 documents what BKS separation (Step 3) requires, and whether it is satisfied.
- §4 articulates the **partial** transfer: Tao's renewal-process recursion μ̂_{n+1}(ξ) = Σ 2^{−v} A_v(ξ) μ̂_n(...) **is** an averaging structure, but it's the wrong one for BKS to flatten.

## 2. Averaging step

### 2.1 What BKS requires

(See L2_FLATTENING_BKS_HYPOTHESIS.md §3.) The averaging step expresses
> μ̂(ξ) = ∫ μ̂_x (D_x f · ξ) dν(x)

where:
- {μ_x} is a **disintegration** of μ along the dynamical phase space (e.g., μ_x = the local geometric piece of μ at point x in the dynamical fiber).
- D_x f is the **derivative cocycle** of the underlying smooth dynamics f at x.
- ν is the SRB / Patterson–Sullivan / Gibbs / IFS-invariant measure on the dynamical phase space (typically the **base measure**, e.g., the equilibrium state of the expanding map).

The key analytic feature: as ξ varies in frequency space, the **image frequencies D_x f · ξ** sweep out a measure-on-frequency-space (namely (D_x f)\_* ν), and the L²-flattening hypothesis is applied to this measure.

### 2.2 Syracuse μ_n's natural averaging structure

The Tao recursion for Syracuse μ̂_n is (R76, R77; see also Tao 1909.03562 §7):
> μ̂_{n+1}(ξ) = Σ_{v ∈ Geom(1/2)} 2^{−v} A_v(ξ) μ̂_n(ξ · 2^{−v} mod 3^n)

where A_v(ξ) is a phase factor from the "+1" step and 2^{−v} is the geometric weight from the v-step decomposition.

This **is** an averaging step (literally: a weighted average over v of values of μ̂_n at scaled frequencies). It has the same shape as the BKS averaging step:
> [continuous BKS] μ̂(ξ) = ∫ μ̂_x(D_x f · ξ) dν(x)
> [discrete Syracuse] μ̂_{n+1}(ξ) = Σ_v p_v · A_v(ξ) · μ̂_n(D_v · ξ)
>   with p_v = 2^{−v}, D_v = 2^{−v} mod 3^n.

So the **averaging step has a clean analog** for Syracuse: the geometric distribution v ~ Geom(1/2) plays the role of the SRB measure ν, and multiplication by 2^{−v} mod 3^n plays the role of the derivative cocycle D_x f.

### 2.3 Where the averaging step fails for BKS purposes

Three problems:

1. **The cocycle is discrete.** D_v is multiplication by 2^{−v} **mod 3^n**, an action on the finite group Z/3^n Z, not a smooth derivative on ℝ^d. BKS Step 2's flattening is applied to the **distribution** of D_x f images, and the relevant non-concentration hypothesis is **affine non-concentration in ℝ^d** — which doesn't make sense for D_v on Z/3^n Z.

2. **The cocycle distribution is discrete-multiplicative-by-2 mod 3^n.** The orbit of ξ under {D_v} is {2^{−v} ξ mod 3^n : v ≥ 1}, which is contained in (Z/3^n Z)\* and has order **2·3^{n−1}** (the order of 2 mod 3^n). The distribution of D_v · ξ (with v ~ Geom(1/2)) is supported on this 2·3^{n−1}-element set. There is no "affine hyperplane" structure on this set; the BKS non-concentration hypothesis doesn't have a direct analog.

3. **The flattening would have to be applied to (D_v) ξ distribution, not π_n.** This is the key conceptual point: in BKS, L²-flattening is applied to the **pushforward of the SRB measure ν by the derivative cocycle**, not to μ itself. For Syracuse this would mean applying L²-flattening to the distribution (D_v)\_* Geom(1/2) on Z/3^n Z — which is the **measure on Z/3^n Z** that puts mass 2^{−v}/Z_v on each 2^{−v} mod 3^n. This distribution is **non-uniform on the cyclic subgroup ⟨2⟩ ⊂ (Z/3^n Z)\***, and it is **NOT** affinely non-concentrated in any continuous-space sense.

### 2.4 Partial transfer assessment

The averaging structure of Tao's recursion is a **structural twin** of BKS's averaging step. This is a positive observation: the abstract architecture of "Fourier transform expressed as an average over cocycle-scaled frequencies" is present in Syracuse. **But** the flattening hypothesis that BKS applies downstream **requires continuity** (affine non-concentration in ℝ^d) and **discrete-finite-group structure doesn't deliver it**.

The averaging step is **structurally analogous** but **technically incompatible** with BKS's flattening machinery.

## 3. Separation / non-linearity hypothesis

### 3.1 What BKS requires

(See L2_FLATTENING_BKS_HYPOTHESIS.md §3, Step 3.) The non-linearity hypothesis requires the **derivative cocycle D_x f to be non-linear** in a quantitative sense — typically expressed as a non-concentration estimate on the distribution of log D_x f (the **Lyapunov spectrum** structure), or equivalently, the cocycle not being cohomologous to a constant.

Concretely, in the Patterson–Sullivan / Gibbs cases, this is the **non-arithmeticity** or **non-integrability** of the cocycle (Sahlsten-Stevens-Sahlsten thesis material, also Bourgain–Dyatlov-style discretized sum-product). In the IFS self-similar case, this is a **Diophantine condition** on the IFS contraction ratios (log r_i and the rotation parts R_i).

### 3.2 Syracuse separation candidate

The Syracuse cocycle is D_v · ξ = 2^{−v} ξ mod 3^n. The "cocycle exponent" is **log 2 mod 3^n**, which is a single value (not a spectrum). The "Diophantine" question is then whether 2 is multiplicatively a primitive root mod 3^n — and **yes**, 2 is a primitive root mod 3^n for all n ≥ 1 (this is classical; ord_{3^n}(2) = 2·3^{n−1} = φ(3^n)).

So **the cyclic group ⟨2⟩ ⊂ (Z/3^n Z)\* generates the full unit group**, and the orbit {2^{−v} ξ : v ≥ 0} for ξ coprime to 3 covers all of (Z/3^n Z)\* with multiplicity exactly 1.

### 3.3 Where the BKS separation hypothesis maps to (or fails)

This **looks** like a non-trivial Diophantine fact (2 is a primitive root, non-arithmetic in the relevant sense). However, the BKS separation hypothesis requires the cocycle to have **continuous non-arithmeticity** — i.e., the log of the derivative cocycle to be "spread out" in a quantitative way on ℝ. The fact that 2 is a primitive root mod 3^n is a **discrete** non-arithmeticity statement; it says nothing about the continuous log structure.

For the BKS framework, the relevant input is more like: the **distribution of log |D_v|** (over the SRB / Gibbs measure on the dynamical phase space) is non-concentrated. For Syracuse, log |D_v| = log |2^{−v}| = −v log 2, with v ~ Geom(1/2) → log |D_v| is concentrated on the discrete arithmetic progression {−k log 2 : k ≥ 1} with geometrically decaying weights. This is **maximally arithmetic** — concentrated on a single arithmetic progression generated by log 2.

So in BKS's Step-3 sense, **the Syracuse cocycle is arithmetic (concentrated on a single log-progression), not non-arithmetic**. The BKS separation hypothesis would NOT hold even if Steps 1 and 2 could be made to work in some weakened form.

The "2 is a primitive root mod 3^n" fact is a **different** non-arithmeticity statement that lives in the discrete-multiplicative group, not in the continuous log-derivative spectrum.

### 3.4 Partial transfer assessment

The Syracuse cocycle has **discrete-arithmetic non-trivial structure** (2 is a primitive root mod 3^n) but is **continuous-arithmetic trivial** (log derivative concentrated on a single AP). BKS Step 3 requires the second kind of non-arithmeticity, which Syracuse doesn't have.

## 4. Bottom-line picture: what transfers and what doesn't

| BKS framework component | Syracuse μ_n analog | Status |
|---|---|---|
| Ambient space ℝ^d | Z/3^n Z (finite discrete) | Mismatch |
| Measure μ on ℝ^d | π_n on Z/3^n Z | Discrete (not continuous) |
| Smooth dynamics f, derivative cocycle D_x f | Multiplication by 2^{−v} mod 3^n | Discrete cocycle (not smooth) |
| Averaging step (Step 1) | Tao recursion μ̂_{n+1} = Σ p_v A_v(ξ) μ̂_n(D_v ξ) | **Structural twin (positive)** |
| L²-flattening hypothesis (Step 2) | Affine non-concentration in Z/3^n Z | **Not well-posed** (mismatch §2 above) |
| Separation / non-linearity (Step 3) | Continuous non-arithmeticity of log|D_x f| | **Fails** (discrete-arithmetic) |
| Output | Polynomial Fourier decay |c=7/45| **Not delivered by this route** |

The averaging step has a structural twin. The flattening step is not well-posed in the discrete setting. The separation hypothesis fails in the form BKS uses.

**Net assessment:** BKS's unified strategy has one structural twin component (averaging) but two structurally incompatible components (flattening, separation). The framework as a whole does not transfer.

## 5. Implication for the polynomial-in-A Fourier-decay obstruction

The polynomial-in-A Fourier bound on |μ̂_n(ξ)| that c = 7/45 closure requires is NOT delivered by Baker–Khalil–Sahlsten 2407.16699 because the framework's L²-flattening hypothesis does not transfer to the Syracuse setting, and even if it did, the separation hypothesis fails.

The closest plausible BKS-style transfer would route through the embedded 3-adic Cantor measure (Translation B in Phase 2), which fails the Diophantine condition and delivers only polylog decay anyway. The discrete-group version (Translation A) is equivalent in difficulty to the polynomial-in-A bound itself (Phase 3 §2.2).

## 6. What's needed for closure now

The BKS route being closed means the polynomial-in-A Fourier bound has to come from elsewhere. Promising directions (already on the INDEX shortlist):

1. **Algom–Baker–Sahlsten cocycle Dolgopyat** (arxiv:2306.01275) — different framework, applies to substitution-shift / cocycle settings. Worth a separate compatibility probe (call it Probe 2).
2. **Furstenberg-measure Rajchman / Hochman-Solomyak** (arxiv:2108.06006 / 1610.02641) — closest setup-match (SL_2 random walk with algebraic entries), Rajchman decay. Worth a separate probe.
3. **Glynn-Zeevi + Lyapunov-Foster-Poincaré** (Glynn-Zeevi / arxiv:2005.08145) — drift-condition fallback, parallel framework. Reasonable next probe if Fourier-decay routes all close.
4. **Bourgain–Konyagin sum-product on (Z/3^n Z)\*** — direct Translation-A attack, but very hard for the structural reason in Phase 3 §2.2 (it's equivalent in difficulty to the target).
5. **Tighter Tao §7.4 bookkeeping** — already INFEASIBLE per BOOKKEEPING_PHASE1_DISPOSITION.

The natural next move per the procedure's "If H_L2_FLATTENING_FAILS, route to Probe 2 (SL_2-embedding) or drift-condition fallback" is therefore **Probe 2 on Furstenberg-Rajchman (arxiv:2108.06006)** as the strongest setup-match, or **Probe 3 on the cocycle Dolgopyat framework (arxiv:2306.01275)** as the strongest generality.

---

End Phase 4.
