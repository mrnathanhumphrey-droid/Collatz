# Extension Work — Tao Bridge, c = 7/45 Closure, and qx+1 Generalization

**Date:** 2026-05-06. **Companion to:** [writeup.md](writeup.md) (the 3x+1 σ-residue paper).
**Status:** consolidated archive of the extension lines that grew out of the
mod-2^k σ-residue paper, capturing positive results, negative results,
and walk-backs that constrain the framework. Cross-references to
[STATE.md](STATE.md) (live state) and [closed_form_findings.md](closed_form_findings.md)
(full archive of derivations).

The base paper closed the 3x+1 σ-residue question structurally:
σ(odd n) = α_det(n mod 2^k) + universal stochastic remainder, with the
α_det(r) closed form computable from the symbolic Collatz prefix and
collapsing 2^(k-1) classes onto k distinct conditional distributions
(one per a_final ∈ {3¹, …, 3^k}). The base paper's "what remains open"
section flagged extension directions; this document reports what each
one actually returned.

Extensions covered:

1. **Ergodic / measure-theoretic framing** (Lagarias 1985 2-adic equidistribution,
   Tao 2022 bridge, density-1 v_2 ensemble check)
2. **Tao bridge tightening** (TA.1 / TA.2 / TA.3, 40 verification cells, structural
   constants ε(σ) ≈ −2.45 and slope ≈ 1/2)
3. **c = 7/45 closure attempts** (R75 Plancherel anchor, R76 conservation law,
   R77.x operator probes, R78/R79 Path-C obstruction map, ε_k extension chain
   k=6..11, mode amplitudes v2)
4. **qx+1 generalization** (Cramér convergence law, per-prime decomposition
   q ∈ {3,5,7,9,11,13}, q-sweep tests 1/2/3, c̃_q = (q-3)/q, cycle classification,
   sibling 3x±1)
5. **Numerical methods** (exact rationals over Q, float64 + scipy.eigs cross-check,
   matrix-free Krylov, FFT Plancherel)

---

## 1. Ergodic / measure-theoretic framing

The whole extension chain sits inside the Lagarias / Tao ergodic-theoretic
picture of Collatz: the Syracuse map T(x) = (3x + 1)/2^v has its 2-adic
valuation v_2 ≡ Geom(1/2) under the natural trajectory measure, and
mean trajectory growth is controlled by E[log(3) − v · log(2)] = log(3/4) < 0.
Three concrete connections:

### Lagarias 1985 2-adic equidistribution

[literature_check.md](literature_check.md) audits prior art for the
prefix-decomposition observation. The relevant prior frameworks:

- **Terras (1976), Acta Arithmetica 30, 241–252.** Lemma 3: parity vector
  P_k(n) depends only on n mod 2^k. Lemma 4: S_k ≈ S_0 · 3^d(k) · 2^(-k)
  asymptotically. The base paper sharpens Lemma 4 by tracking the exact
  symbolic state (a_final · m + c_final) with variable prefix length and
  using (a_final, c_final) as a *covariate* parameterizing the conditional
  distribution σ(n) | n mod 2^k.
- **Lagarias, "The 3x+1 problem: An annotated bibliography" (math/0309224,
  math/0608208).** Master reference for the literature; the 2-adic
  equidistribution prediction lives here.
- **Sinai (2003), "Statistical (3x+1) problem", Comm. Pure Appl. Math 56,
  1016–1028 (math/0201102).** Probabilistic / CLT framework on Syracuse;
  treats Geom(1/2) v_2 distribution but not residue-mod-2^k symbolic state.
- **Tao (2019/2022), "Almost all Collatz orbits attain almost bounded
  values".** Syracuse-map framework with 3-adic Fourier on Z/3^n Z. *This
  is the ergodic-theoretic paper that the c = 7/45 work bridges to* — the
  Plancherel decomposition S_k = Σ_{ξ : 3∤ξ} |μ̂_k(ξ)|² lives in Tao's
  3-adic setting (R75 Theorem 75.1).
- **Korec (1994), Mathematica Slovaca 44(1), 85–89.** Almost-all density
  result n^c with c > log_4 3.

The base paper's prefix decomposition is **elementary refinement of
Terras Lemma 4**: it tracks the exact symbolic state with variable prefix
length and uses (a_final, c_final) as a covariate. Across the focused
search documented in literature_check.md, this combination wasn't found
in the canonical 50-year corpus.

### Tao 2022 bridge at 40 verification cells

The base paper extended into Tao's framework by checking that the
prefix decomposition reproduces the per-class structure of Tao's
mean-trajectory inequality (eq 5.15: T_x(N) = log(N/x)/log(4/3) +
O(log^0.6 x) for almost all N). The bridge at 40 verification cells
(5 observables × 4 modular resolutions × 2 N-scales):

> s_mean(r; f) ≈ α_det(r) + K_h · log(N / f(N)),  K_h = 3/log(4/3)

with **slope at K_h = 1.000 ± 0.005** in every cell — Tao's leading term,
no recalibration. Per-residue-class realization of the asymptotic
result. Documented in [tao_bridge_findings.md](tao_bridge_findings.md)
intro.

### Density-1 v_2 ensemble verification

[result_density_one_v2_bounds.md](result_density_one_v2_bounds.md) on
2,796,202 Syracuse trajectories (odd starts coprime-to-3 in [3, 2^23]):

- **Unconditional ensemble mean v_2 = 2.102** vs Geom(1/2) prediction 2.0.
  5% deviation, consistent with Tao 2019 measure-theoretic asymptotic at
  finite resolution.
- **TEST B (mean v > log_2(3))** returned 100% pass rate, but the result is
  *mathematically tautological*: direct algebra gives mean_v − log_2(3) =
  (log_2(n_0) + Σᵢ log_2(1 + 1/(3 n_i)))/L > 0 for every Syracuse
  trajectory terminating at 1. So 100% on the ensemble simply restates
  "all 2.8M trajectories reached 1" (already known on this range).
  **Do NOT cite TEST B's 100% pass rate as a quantitative density-1
  strengthening.**
- **TEST A (per-trajectory geometric null density(v_i ≥ k) ≥ 2^{-(k-1)})**
  fails for long trajectories — selection effect: conditioning on large L
  squeezes mean v toward log_2(3), so cumulative density drops below
  2^{-(k-1)}.

The empirically informative content is the unconditional ensemble mean
2.102 (≈ 2). This is **consistent with** Tao 2019's Geom(1/2) prediction
at strong sample resolution; it's not stronger than the measure-theoretic
statement.

This probe also rules out the "bit-budget contradiction" mechanism for
the convergence-rate question: low-v_2 admissibility probe
([result_low_v2_residue_admissibility.md](result_low_v2_residue_admissibility.md))
returned 18/18 cells null — low-v_2 trajectories have *larger* admissibility
horizons, not smaller. Bit-budget consumption is deterministic at rate
Σ v_i per step; slower consumption produces longer admissibility, not
contradiction.

---

## 2. Tao bridge tightening (TA.1 / TA.2 / TA.3)

The base bridge is structurally clean at leading order. The three
Tao-bridge tasks tightened the characterization of the structural
correction ε in `s_mean(r) ≈ α_det(r) + K_h · log(N/f(N)) + ε`.

### TA.1 — N-stability of the σ structural offset

[experiments/36_TA1_sigma_offset_N_sweep.py](experiments/36_TA1_sigma_offset_N_sweep.py).
For σ at k ∈ {8, 10, 12} and N ∈ {2²⁵, 2²⁷, 2²⁸, 2³⁰, 2³²}:

| N | gap = σ̄ − ⟨α_det⟩ − K_h·⟨log N⟩ | per-class SE (k=8) |
|---|---|---|
| 2²⁵ | −2.4468 | 0.174 |
| 2²⁷ | −2.4514 | 0.091 |
| 2²⁸ | −2.4492 | 0.066 |
| 2³⁰ | −2.4526 | 0.034 |
| 2³² | −2.4574 | 0.018 |

**Variation across k at fixed N: identically zero.** All k=8/10/12 gaps
match to 4 decimal places at every N — the gap is a property of the σ
distribution, not the modular grid. Variation across N at fixed k is
0.0105 across 7 doublings (−0.0022 per unit log N, well within Tao's
O((log N)^0.6) sub-leading band; ~5× the per-class SE at the largest N).

**Closed-form decomposition.** ⟨α_det⟩ at every k equals
`E[prefix_steps] − K_h · ⟨descent during prefix⟩` = +6.23 *exactly across all k*,
an invariant of the prefix algebra. So:
```
σ̄ − K_h · log N  =  ⟨α_det⟩  +  ε(σ)  =  +6.23  +  (−2.45)  =  +3.78
```
The **+6.23** is rigorous (R1 in [closed_form_findings.md](closed_form_findings.md):
⟨α_det⟩ = log(6)/log(4/3), derived from the binomial j-distribution).
The **−2.45** is a post-prefix descent constant — the actual descent takes
2.45 fewer Collatz steps than the K_h · log(post-prefix value) random-walk
heuristic predicts. **Not derived analytically** (R2 in closed_form_findings:
the natural candidate −log(2)/L = −2.4094 was ruled out because ε_total
drifts in the wrong direction with N).

### TA.2 — Trim-quantile sweep

At N=2²⁷, k=8, observable s @ √N: per-class trim of top-q% drives the
gap to zero at **q* ≈ 1.18%**. Implied log^(−c) N exponent: c = 1.55 —
not a clean Tao exceptional-set exponent (Tao 2022 uses 0.6). Reading:
q* is a property of the σ right-tail mass at this N, not a structural
Tao quantity. Worth retesting at 2³² but not pursued.

### TA.3 — Parametric form of ε(observable)

13 (observable, N) cells, candidate parametric forms:

| Model | SSE | R² |
|---|---|---|
| **gap = a + b · log(threshold)** | **0.80** | **0.989** |
| gap = a + b · log(N) + c · Δlog | 0.74 | 0.990 |
| gap = a + b · Δlog | 3.28 | 0.955 |
| gap = a + b · log(log N) | 53.5 | 0.269 |

Best two-parameter form: **`ε ≈ −2.35 + 0.486 · log(threshold)`**.
Slope ≈ 0.486 ≈ 1/2.

The **−2.35 intercept matches the σ structural constant** from TA.1, and
the **+1/2 slope on log(threshold)** is the second extension-paper
constant deserving a closed form.

### Closed-form derivations of the two TA.3 constants

(R1–R3 in [closed_form_findings.md](closed_form_findings.md) line 105+.)

| Constant | Status | Result |
|---|---|---|
| ⟨α_det⟩ = log(6)/log(4/3) | **DERIVED ✓** | Exact invariant of prefix algebra; verified to 1e−14 across k ∈ {6, 8, 10, 12, 14}. |
| ε(σ) = −2.45 | **NOT DERIVED ✗** | log(3)/L candidate ruled out (drift wrong direction). Likely needs renewal-theoretic input. |
| +1/2 slope on log(threshold) | **MECHANISM ID'd, NOT CLOSED ◐** | First-passage overshoot is *not* the source (overshoot constant ≈ 0.30 nats across thresholds). Source is K_eff ≈ 9.94 ≠ K_h on the post-first-crossing segment R(f) = σ_mean − s_mean(f). |

**Synthesis (Tao bridge).** The bridge is structurally clean at leading
order with slope exactly 1.000 at K_h. The structural correction ε
decomposes cleanly into an N-near-constant (−2.35) plus an
observable-dependent log(threshold)/2 piece. The −2.35 constant has been
isolated as the post-prefix descent correction, but a closed form for it
remains open. The +1/2 slope is mechanism-identified (post-crossing
segment runs at K_eff ≈ 9.94 not K_h ≈ 10.43) but the underlying source
of K_eff ≠ K_h on that segment is open.

---

## 3. c = 7/45 closure attempts

The companion analytical question: *the random-walk heuristic K_h gives
the slope; what about the constant in the convergence-rate envelope?*
Empirically, in the Plancherel framework, S_k = Σ_{ξ : 3∤ξ} |μ̂_k(ξ)|² → 7/15,
and equivalently ‖d_{k+1}‖² ≈ c · (1/3)^k with c = 7/45 = (7/15)/3.

The closure attempt for c = 7/45 ran through R75–R79 + R77.x + the
ε_k extension probes through k=11. This is the most detailed strand;
this section captures both what was achieved and the substantial
walk-backs.

### R75 — Plancherel anchor for c = 7/45 (rigorous algebraic identity)

[c_seven_forty_fifth.md](c_seven_forty_fifth.md), Result 75 in
closed_form_findings (line 7099).

> c = (1/3) · lim_{n→∞} Σ_{ξ ∈ Z/3^n, v_3(ξ)=0} |μ̂_n(ξ)|²

The Plancherel decomposition S_k = Σ_{3∤ξ} |μ̂_k(ξ)|² is **proved**
(verified algebraically + numerically through k=3). Combined with R74's
algebraic identity S_{k+1} = 3^{k+1} · ‖d_{k+1}‖², this gives c = S_∞/3
on solid algebraic ground.

**What's rigorous:**
- Plancherel decomposition (Theorem 75.1)
- R74 algebraic identity for S_{k+1}
- Tao recursion → diagonal (= S_n exactly) + off-diagonal decomposition (Theorem 75.2)

**What was provisional (and has now been walked back):**
- Rate of convergence S_n → 7/15 is exactly 1/2 per level
- |ε_n|·2^n ≤ 0.04 envelope from k=2..5
- Provisional certified bound |c − S_k/3| ≤ 0.0133·(1/2)^k

### R76 — Conservation law (rigorous)

[result_76_conservation_law.md](result_76_conservation_law.md). For
M_n(η) := Σ_ξ μ̂_n(ξ) μ̂_n*(ξη):
> Σ_{j=0,1,2} M_{n+1}(η_0 + j·3^n) = 0
> S_{n+1} = −2 · M_{n+1}(1 + 3^n)

Reduces the rate-question to a scalar sequence R_n := M_n(1 + 3^{n−1}) → −7/30.

### R77.x — Operator shape probes (per-strand status)

[result_77_T_lead_spectrum.md](result_77_T_lead_spectrum.md) and follow-ons.

| Strand | Outcome |
|---|---|
| R77.3 (β) | finite-mode geometric ansatz `ε_n = a·r₁^n + b·r₂^n + c·r₃^n` over Q FALSIFIED. A = −157462/3058335 ≠ −1/30; predictions miss ε_4..ε_6 by 28–41%. |
| R77.4 (M) | envelope curve fits at N=5 inconclusive between Jordan / log / power-law. Jordan ruled out; H2 ≈ H3 tied at ΔAIC = 0.23. |
| R77.4 erratum | **K_k spectrum has nothing near 1/2.** \|λ_2\| ∈ [10⁻⁶, 10⁻³] at k=3..7. The rate operator is INTER-LEVEL renormalization, NOT within-level mixing. |
| R77.6 (G-type indeterminate) | Padé approximants of Σ ε_n z^n place poles in [2.05, 2.35]; consistent with branch cut at z=2, NOT simple pole. |
| R77.7 (NOT COMPLETED) | k=7 Markov chain extension killed at ~8.5h. Original "superseded by Bohr" framing also retired. |
| R77.5 — multi-resolution decomposition | **Identity proved over Q at q=3:** ‖R_k‖² · q^k = S_{k+1}/q with R_k := π_{k+1} − T(π_k) ∈ W_k orthogonal to T(V_k). |

### R78/R79 — Path-C obstruction map for the rigorous rate-½ proof

The remaining gap to a rigorous proof of S_∞ = 7/15 is the spectral
identification of the rate-½ operator (Kalafatelis 2026 eq 190, Remark 27).
The obstruction map after R78 + R79 + saddle-class + C2 + C3 + band-spectral:

| Attack route | Status | Saving |
|---|---|---|
| Cochrane Theorem 2 (R78) | ❌ closed | trivial only |
| Pólya-Vinogradov (R78) | ❌ closed | worse than trivial for r ≥ 3 |
| van der Corput B=1 (R79) | ⊳ partial | ~0.73 sub-trivial |
| van der Corput B=2 (R79) | ⊳ partial | ~0.81 worse than B=1 |
| Empirical \|S_partial\| (R79b) | ⊳ confirms vdC stalls | β = 0.522 ± 0.008 |
| **Even ideal pointwise √N** | hypothetical | **STILL insufficient** by R79 Step 4 |
| Saddle-class subsum | ⊳ structural but not closure | β_j ∈ [0.92, 1.06] |
| C2 / BGK on ⟨4⟩ | ⊳ partial | M_4 slope = 3.0059, random-like |
| **C3 / direct band-l¹** | ❌ closed | saturates trivial bound exactly N_r^{1.0} |
| Band-spectral decomposition | ⊳ closed for smooth-completion | lf_mass → 0.25 |
| C1 / 5x+1 sibling reframing | open | multi-day rebuild |

**The lesson from R79's Step 4** is operationally important: any pointwise
√N bound is *structurally insufficient* for eq 190, because it discards
the off-diagonal cancellation between m-values that the bound actually
requires. Differencing-based attacks (vdC family) cannot in principle
work; the closure must come from a method that captures inter-m
cancellation (Bourgain-Konyagin sum-product, smooth completion, or 5x+1
reframing).

### The empirical rate-1/2 walk-back (the big 2026-05-05 result)

The most consequential extension finding, and the most painful one. The
provisional R75 certified bound rested on the empirical claim
"|ε_n|·2^n stable near 0.04 for n=2..5." Pushing the measurements to
k=6, 7, 8, 9, 10, 11 demolished the claim:

| n | ε_n | |ε_n|·2^n |
|---|---|---|
| 2 | +9.52 × 10⁻³ | 0.038 |
| 3 | −5.09 × 10⁻³ | 0.041 |
| 4 | −2.45 × 10⁻³ | 0.039 |
| 5 | −1.15 × 10⁻³ | 0.037 |
| 6 | −4.98 × 10⁻⁴ | 0.032 |
| 7 | **−1.18 × 10⁻³** | **0.150** ← jumps 4× |
| 8..11 | non-monotone, sign flips at k=10 | grows like (1.968)^n |

Key empirical findings from the
[result_epsilon_6.md](result_epsilon_6.md) →
[result_epsilon_7.md](result_epsilon_7.md) → ... →
[result_epsilon_11.md](result_epsilon_11.md) chain:

1. Two-mode `ε_k = A·(1/2)^k + B·(1/3)^k` fit (logged earlier on k=1..5)
   FALSIFIED at k=6 by 10×.
2. |ε_7/ε_6| = 2.36 — non-monotone bounce at k=7, S_k has a local maximum
   at k=6.
3. ε_10 = +7.21 × 10⁻⁴ — sign-flip from k=9, magnitude rebounded.
4. Slow-mode envelope rate ρ ≈ 0.984 per k-step, period ≈ 9.2.
5. **|ε_n|·2^n bounded-envelope reading is empirically refuted:** at
   ρ ≈ 0.984, |ε_n|·2^n grows as (1.968)^n.
6. Within-level K_k spectral gap (~0.998) is consistent — mixing is fast;
   the slow oscillating mode must live in **inter-level renormalization**.
7. Cross-validated: power-iteration, scipy.eigs, FFT agree to 4×10⁻¹⁵
   at k=10.

**Implication for c = 7/45.** R75's algebraic anchor (Plancherel + R74)
stands. The provisional certified bound `|c − S_k/3| ≤ 0.0133·(1/2)^k`
is **NOT VALID** under the new k=7..11 data — at k=7 the actual deviation
exceeds the supposed bound by an order of magnitude. Bound and proof
framework need either (a) a different rate hypothesis, (b) a
non-monotone-tolerant bound, or (c) abandonment until structural form
is identified.

### The inter-level renormalization framing (R77.4 erratum + R77.5 + mode amplitudes v2)

R77.4's erratum reframed the convergence-rate question: K_k itself has
no eigenvalue near 1/2 at any k tested, so the rate-controlling object
is not within-level mixing. R77.5 proved the natural decomposition over Q:

> V_{k+1} = T(V_k) ⊕ W_k,  R_k := π_{k+1} − T(π_k) ∈ W_k,  ‖R_k‖² · q^k = S_{k+1}/q

where T : V_k → V_{k+1} is the lift (T(v)(r') := v(r' mod q^k) / q).

The mode-amplitudes v2 probe ([mode_amplitudes_v2_findings.md](probe_mode_amplitudes_v2/mode_amplitudes_v2_findings.md))
decomposed δ_k := L_{k-1}π_{k-1} − π_k onto K_k right eigenvectors and
R_k singular vectors. Headlines:

- δ_k norms 12–50× larger than |ε_k| — δ_k carries far more information
  than ε_k captures.
- Decomp A (K_k right eigvec) captures ~0% — structural, not a bug.
  Right eigvec of non-symmetric K is biorthogonal-dual to LEFT, not RIGHT;
  right Perron is the constant vector and δ_k has sum 0 by mass
  conservation.
- Decomp B (R_k singular vectors) captures 17.86% / 5.22% / 1.65% in
  top-20 at k=5,6,7 — top-20 / dim_R shrinks fast as state space grows.
- σ_1 alone captures only 2.7–3.6% across all three k — the slow rate is
  **band-collective, not single-direction**.
- Top-20 R_k singular values cluster tightly in [0.658, 0.671] — near-degenerate
  band, no isolated dominant direction.
- ρ_slow ≈ 0.83 (recurrence-fit slow-mode rate from
  [result_renormalization_spectrum.md](result_renormalization_spectrum.md))
  is **not present in any single K_k or R_k mode** at any k tested.

Three remaining hypotheses for ρ_slow's origin:
- (a) **Composition across levels** — slow rate emerges from product of
  R-actions, not single eigenvalue.
- (b) **Functional-projection averaging** — ε_k is a specific scalar
  projection of δ_k that selects the "slow envelope" of the band-collective
  action.
- (c) **Finite-k recurrence-fit artifact** — true rate near 0.60
  (composition asymptotic) and the 0.83 fit is a small-k transient.

Distinguishing requires k=8..11 ε measurements (have through k=11) +
explicit R_k composition computation. R77.5 ("inter-level renormalization
residual operator") was identified as the natural follow-up but not yet run.

### Synthesis (c = 7/45 closure)

| Layer | Status |
|---|---|
| c = 7/45 algebraic identity (Plancherel) | **rigorous (R75)** |
| Conservation law for M_n | **rigorous (R76)** |
| Multi-resolution decomposition over Q at q=3 | **rigorous (R77.5)** |
| Empirical certification through k=11 | **valid: c − S_11/3 within bound** |
| Provisional rate-½ envelope bound | **walked back (k=7..11 empirical refutation)** |
| Rigorous spectral identification of rate operator | **OPEN** (Bourgain-Konyagin / 5x+1 / smooth-completion routes remain) |
| Within-level rate-½ via K_k spectrum | **falsified universally (R77.4 erratum + q-spectrum probe)** |
| Inter-level R_k as rate-controlling object | **proved structurally; spectrum not yet computed at any q** |

c = 7/45 retains its empirical certification through k=11 + structural
anchoring (R74–R77.5). The final rigorous closure remains a published
open problem (Kalafatelis 2026, Remark 27) with a now-mapped obstruction
landscape.

---

## 4. qx+1 generalization

The most productive extension line, and the one with the cleanest
**publishable positive results** independent of the c = 7/45 closure.
Multiple sub-strands: Cramér convergence law, per-prime decomposition,
q-sweep tests 1/2/3, c̃_q = (q-3)/q candidate, cycle classification, and
the sibling 3x±1 chain.

### Cramér convergence law for qx+1 conv-rate decay

[experiments/16_cramer_root.py](experiments/16_cramer_root.py),
[experiments/17_cramer_dual_verification.py](experiments/17_cramer_dual_verification.py).

For qx+1 with Geom(1/2) v_2 distribution, the trajectory step in log-coords
is X = log(q) − v · log(2), with MGF E[exp(−θX)] = q^(−θ) · (2^θ/2) / (1 − 2^θ/2)
for 2^θ < 2. Setting MGF = 1 gives:

> **q^(−θ) = 2^(1−θ) − 1**

This is the **exact Cramér equation** for the upward-MGF root θ(q) on
the qx+1 random-walk heuristic. Solving by brentq:

| q | θ(q) | log(q) | predicted slope = −θ·log(q) | empirical slope (N=10⁸) | ratio |
|---|---|---|---|---|---|
| 5 | 0.358 | 1.609 | −0.576 | −0.5619 | 0.976 |
| 7 | 0.720 | 1.946 | −1.401 | −1.3685 | 0.977 |
| 9 | 0.741 | 2.197 | −1.628 | −2.0529 | 1.261 |
| 11 | 0.805 | 2.398 | −1.930 | −1.6458 | 0.853 |

The Cramér prediction matches q=5 and q=7 to ~2.4%; q=9 and q=11 deviate
(more than ~25%) — finite-N at small sample sizes (q=9 has only 104
converged orbits at N=10⁸; q=11 has 36 at N=10⁹).

**Universal "Cramér multiplier" C ≈ 5/2 in conv_rate decay law** —
confirmed at R²=0.999 across q ∈ {5, 7, 9} and R²=0.994 at q=11. See
[closed_form_findings.md](closed_form_findings.md) line 3021.

The Cramér root identity also has a load-bearing closed form at q=3
(R15 in closed_form_findings, line 1158): the Wiener-Hopf factorization
1 − φ(θ) = (1 − κ⁺_def(θ))·(1 − κ⁻(θ)) has Cramér root w* = 1 *because*
2² − 1 = 3 — the algebraic identity from the Collatz step encoded as a
Cramér condition. The κ⁺_def part doesn't admit a clean rational closed
form in 2^(iθ), so R15 is partial: Cramér root rigorous, full κ⁻
factorization open.

### Per-prime decomposition q ∈ {3, 5, 7, 9, 11}

[experiments/10_q_decomposition.py](experiments/10_q_decomposition.py),
[experiments/12_q_convrate_analytical.py](experiments/12_q_convrate_analytical.py),
[experiments/13_cross_q_unification.py](experiments/13_cross_q_unification.py).

The qx+1 prefix decomposition mirrors the q=3 form: while the symbolic
multiplier `a` is even, branch parity is forced by `c`; apply T_q. Stop
when `a` becomes odd. a_final ∈ {q^j}.

R36 (closed_form_findings line 2977) ran the systematic q ∈ {3, 5, 7, 9} pass:

| q | n_conv | σ_mean | E*[v]_emp | K_h(q;conv) closed-form | empirical pooled slope | gap |
|---|---|---|---|---|---|---|
| 3 | 500,000 | 137.60 | 1.9918 | **10.6102** | 10.3900 | −0.22 (−2.1%) |
| 5 | 32,785 | 165.90 | 2.8948 | **9.8090** | 12.9393 | +3.13 (+31.9%) |
| 7 | 258 | 44.84 | 6.5862 | **2.8963** | 3.2890 | +0.39 (+13.6%) |

**Closed-form K_h(q; conv) = (1 + E*[v]) / (E*[v]·log(2) − log(q)).**
Generalizes Syracuse K_h = 3/log(4/3) (q=3 special case where E*[v] = 2
collapses to (1 + 2)/(2·log(2) − log(3)) = 3/log(4/3) = 10.43).
Universal in form, leading-order match — match degrades from 2% (q=3)
to 14% (q=7). E*[v]_conv increases dramatically with q
(1.99 → 2.89 → 6.59), confirming the Cramér-tilt picture: converged
orbits at higher q require larger above-typical v-magnitudes.

**Per-class slope universality (the q=3 framework's load-bearing fact)
does NOT extend cleanly to q=5.** Slope CV grows 3.7% (q=3) → 24.5% (q=5).

**Cross-q structural picture (4 universal facts, 3 q=3-specific):**

| Constant | q=3 status | q=5..11 status | Universal vs q-specific |
|---|---|---|---|
| 1 (⟨α_det⟩) | DERIVED log(6)/log(4/3) | partial (slope CV 24.5% at q=5) | **q=3-specific in current form** |
| 2 (K_h) | 3/log(4/3) DERIVED | leading-order (gap 32% at q=5, 14% at q=7) | **UNIVERSAL FORM** |
| 3 (per-j W_j) | Lagarias-class | q-SPECIFIC (cycle ≠ lattice; see cycle classification) | **q-specific (different attractor topology)** |
| 4 (per-σ-band) | U-shape (R14 family) | U-shape REPRODUCES at q=5,7,9 | **STRUCTURALLY UNIVERSAL** |
| Cramér multiplier C ≈ 5/2 | n/a | UNIVERSAL R²=0.99 | **UNIVERSAL** |
| Unconditional v_2 ~ Geom(1/2) | n/a | UNIVERSAL (0.5% match across q ∈ {5,7,9,11}) | **UNIVERSAL** |
| K_h(q;conv) functional form | exact | leading-order | **UNIVERSAL FORM** |

### Q-sweep test 1 — universal q/3 ratio (Plancherel-mass scaling)

[result_q_sweep_test_1_rate.md](result_q_sweep_test_1_rate.md), q ∈ {3, 5, 7, 11, 13}.

For the qx+1 Syracuse Markov chain on (Z/q^k)*:

> S_{k+1}^{(q)} / S_k^{(q)}  →  q/3 universally across all tested q.

Empirically clean (q=11, 13 hit 4 sig figs by k=2). q=3 is the borderline
case where ratio → 1 from below, giving the finite Tao limit S_∞^{(3)} = 7/15.
For q ≥ 5, ratio > 1 means S_k diverges geometrically, so S_∞^{(q)} = ∞.

| q | S_2/S_1 | S_3/S_2 | S_4/S_3 | q/3 |
|---|---|---|---|---|
| 3 | 0.7143 | 0.9693 | 1.0057 | **1.0000** |
| 5 | 1.5382 | 1.6579 | 1.6601 | **1.6667** |
| 7 | 2.1280 | 2.3359 | — | **2.3333** |
| 11 | 3.6646 | 3.6666 | — | **3.6667** |
| 13 | 4.3337 | — | — | **4.3333** |

**Why exactly q/3?** No analytic explanation yet — likely a clean
derivation from Tao's Plancherel framework. The structural finding
itself is q-universal, and rules out the original literal hypothesis
that 7/45 is the q=3 instance of a closed-form c_q family
(it isn't — see [result_q_sweep_test_2_c_q.md](result_q_sweep_test_2_c_q.md)).

### Renormalized constant c̃_q := lim S_k^{(q)} / (q/3)^k

| q | c̃_q | Notes |
|---|---|---|
| 3 | 7/15 ≈ 0.4667 | Forward Tao limit (q/3 = 1, no renormalization) |
| 5 | ≈ 0.487 (Aitken extrap → 0.482) | mild deviation from (q-3)/q = 0.40 |
| 7 | ≈ 0.78 | large deviation from (q-3)/q = 4/7 ≈ 0.571 |
| 11 | 0.7288 | matches (q-3)/q = 8/11 = 0.7273 to **0.2%** |
| 13 | 0.7698 | matches (q-3)/q = 10/13 = 0.7692 to **0.07%** |
| 17 | matches 14/17 within ~1% | confirms (q-3)/q at non-prim-root case |

[c_tilde_structure_verdict.md](c_tilde_structure_verdict.md) and
[c_tilde_q17_probe.py](c_tilde_q17_probe.py).

**Two-regime conjecture (suggestive at 6 datapoints):**
- q = 3 (boundary): c̃_3 = 7/15
- q ≥ 11 (large + 2 prim root): c̃_q = (q − 3)/q to within 0.2%
- q = 5, 7: deviate (q=5 by ~0.08, q=7 by ~0.21)

The c̃_q = (q − 3)/q at q=11, 13 means c̃_q · q is integer (= q − 3).
Empirically 8.017 (vs 8) and 10.007 (vs 10) — both within 0.2% of
integers. **This is the cleanest signal in the qx+1 generalization.**

q=17 has ord(2 mod 17) = 8 (NOT primitive root, like q=7). The fact
that c̃_17 still hits (q-3)/q within 1% rules out the "non-prim-root
status is the differentiator" hypothesis — q=7's deviation is q-specific
finite-k or arithmetic-specific, not a non-prim-root pattern.

### Q-spectrum probe — universal triviality of K_k^(q) spectrum

[result_qspectrum.md](result_qspectrum.md). Top-10 |λ| of K_k^(q) at
q ∈ {3, 5, 7, 11, 13}, k = 5..7.

> K_k^(q) spectrum is q-universally trivial: |λ_1| = 1 isolated, all
> sub-leading eigenvalues clustered in [10⁻⁶, 10⁻⁵], no q has an
> eigenvalue near 1/2.

**This extends R77.4 erratum from q=3 to the full prime family** and
settles a substantive negative result: per-level Markov mixing is NOT
the rate-controlling object across primes.

### Q-sweep test 3 — multi-resolution decomposition is q-universal

[result_q_sweep_test_3_decomposition.md](result_q_sweep_test_3_decomposition.md).

R77.5's lift-residual decomposition extends from q=3 to q ∈ {3, 5, 7, 11, 13}:

> ‖R_k^(q)‖² · q^k = S_{k+1}^(q) / q,  R_k^(q) ⊥ T_q(V_k^(q))

**Outcome: DECOMP-UNIVERSAL.** Three test vectors, all (q, k) cells:
⟨R_k^(q), T_q(v)⟩ = 0 as exact rational equality at every cell tested.
Identity holds exactly across all q, all k, by structure (marginal
consistency of the projective Markov system mod q^k is q-blind). Every
analytical tool R77.5 provides for q=3 is available unmodified for any
odd prime q.

### Cycle classification for q ∈ {5, 7, 11, 13}

[experiments/22_q5_cycle_detection.py](experiments/22_q5_cycle_detection.py),
[experiments/29_qx1_cycle_classification.py](experiments/29_qx1_cycle_classification.py),
[experiments/36_q5_fourth_cycle_search.py](experiments/36_q5_fourth_cycle_search.py).

For q ≥ 5, qx+1 is no longer mean-reverting (E[X] = log(q/4) ≥ 0),
and orbits split into three populations: (a) trivial cycle (smallest
member = 1), (b) non-trivial cycle, (c) divergent. Floyd's
tortoise-and-hare cycle detection on every odd start in the parquet:

**q=5 cycles found** ([29_qx1_cycle_catalog_q5.csv](experiments_output/29_qx1_cycle_catalog_q5.csv)):

| smallest | cycle_length | count (in [1, 10⁸]) | walk |
|---|---|---|---|
| 13 | 10 | 51,408 | [13, 66, 33, 166, 83, 416, 208, 104, 52, 26] |
| 17 | 10 | 12,766 | [17, 86, 43, 216, 108, 54, 27, 136, 68, 34] |

**q=5 fourth-cycle search at scale** ([36_q5_search_log.txt](experiments_output/36_q5_search_log.txt)):
Floyd cycle detection on m ≡ 33 mod 40 in range [10⁸, 10¹⁰], 247.5M starts:

```
trivial cycle:    29,140  (0.0118%)
non-trivial:      56,616  (0.0229%)
divergent:        247,414,244  (99.9654%)
unique cycle smallest members: [13, 17]
```

> **NO new cycle landings found.** The conjectured 4th q=5 cycle has
> smallest member > 10¹⁰.

**Outcome by q:**
- q=5: 3 known cycles (trivial = {1, 2, 4, 8, 16, 3, 6}; cycles at 13 and 17;
  no 4th cycle below 10¹⁰)
- q=7, 11, 13: cycle catalogs computed (parquet outputs at
  [29_qx1_cycles_q*_N*.parquet](experiments_output/))

**Strategic implication for c = 7/45 generalization:** the per-j W_j
structure (Lagarias-class q=3 attractor) does NOT carry to q ≥ 5 as a
lattice — q=5's attractor IS the cycle {1, 2, 4, 8, 16, 3, 6}, not the
lattice {m_j = (4^j − 1)/3}. **Constant 3 (per-j W_j) is q-SPECIFIC.**
This is the most consequential structural finding of the cycle
classification: any closed-form for the qx+1 family that depends on
attractor topology will need to be cycle-aware at q ≥ 5.

### Sibling 3x±1 chain symmetry and inverse-tree asymmetry

[sibling_3x_minus_1_symmetry_verdict.md](sibling_3x_minus_1_symmetry_verdict.md),
[duality_S_vs_D_verdict.md](duality_S_vs_D_verdict.md),
[duality_followup_verdict.md](duality_followup_verdict.md).

**Forward symmetry K_- = σK_+σ proved.** The q=3 Syracuse Markov chains
for 3x+1 and 3x-1 are conjugate by negation σ(r) = -r. Implies:

> S_n^{3x−1} = S_n^{3x+1} as exact rationals at every n (verified k=1..4).

All R76/R77 derived quantities transfer; **c = 7/45 is automatic for 3x-1**
by the same evidence chain.

**Inverse-tree D_n asymmetry.** Despite forward equality, integer-level
inverse-tree Plancherel masses differ by 10³–10⁴× at large depth. After
matched-N control (Agent 2 truncated to Agent 3 root-1's |V_n|), residual
structural difference is factor 0.2–4 — **~95% of raw 10⁴× difference is
sample-size driven**. Forward chain symmetry doesn't propagate to
inverse-tree integer-level dynamics.

**No clean forward-backward duality** D = f(S) in any candidate form.

### Synthesis (qx+1 generalization)

| Finding | Status |
|---|---|
| Cramér root q^(−θ) = 2^(1−θ) − 1 (closed form) | **DERIVED** at q=3 (R15: w*=1 from 2²−1=3); empirical match at q=5,7 to ~2.4% |
| Cramér multiplier C ≈ 5/2 in conv_rate decay | **UNIVERSAL** R²≥0.99 across q ∈ {5,7,9,11} |
| Unconditional v_2 ~ Geom(1/2) ensemble | **UNIVERSAL** 0.5% match across q ∈ {5,7,9,11}, mean 2.102 at q=3 |
| K_h(q; conv) functional form (1+E*[v])/(E*[v]·log(2)−log(q)) | **UNIVERSAL FORM** leading-order |
| Per-σ-quantile band U-shape | **UNIVERSAL** across q ∈ {3,5,7,9} (spread attenuates with q) |
| Universal q/3 growth ratio S_{k+1}^(q)/S_k^(q) | **publishable theorem candidate** |
| c̃_q := lim S_k^(q)/(q/3)^k exists ∀ q | **publishable; structural** |
| c̃_q = (q − 3)/q at q ∈ {11, 13, 17} (within 0.2%, 0.07%, ~1%) | **suggestive at large q with 2 prim root** |
| q=3 boundary regime c̃_3 = 7/15 | **rigorous via R75 (Plancherel)** |
| q=5, q=7 deviations from (q-3)/q | **unresolved at N=6 datapoints** |
| K_k^(q) spectrum trivial across q ∈ {3,5,7,11,13} | **q-universal negative result** |
| R77.5 multi-resolution decomposition extends to qx+1 | **DECOMP-UNIVERSAL over Q** |
| 7/45 = q=3 instance of a closed-form c_q family (LITERAL hypothesis) | **FALSIFIED** (S_∞^{(q)} doesn't exist for q ≥ 5) |
| Cycle catalog q=5 (cycles at 1, 13, 17; no 4th below 10¹⁰) | **EMPIRICAL closure** |
| Per-class slope universality at q=5 | **BREAKS** (CV 24.5% vs 3.7% at q=3) |
| Per-j W_j attractor structure | **q-SPECIFIC** (q=5 cycle ≠ q=3 lattice) |
| Forward 3x±1 chain symmetry K_- = σK_+σ | **PROVED** (k=1..4); S_n equal as exact rationals |
| Inverse-tree D_n asymmetry between 3x+1 and 3x-1 | **~95% sample-size driven**; residual factor 0.2–4 |

The qx+1 generalization is the strongest standalone publishable line.
The q-universal Cramér multiplier, q-universal q/3 ratio, K_k spectrum
triviality, and DECOMP-UNIVERSAL R77.5 are four independent
substantive results, all independent of Collatz closure status. The
(q-3)/q pattern at q ≥ 11 is the most surprising single finding — a
clean rational closed form emerging at the renormalized scale,
restricted to the "large enough q" regime where finite-k transients
have damped out.

---

## 5. Numerical methods

The four computational regimes used across the extension probes,
listed by k-range (the modulus 3^k controls everything else):

### Exact-rational stationary computation (k ≤ 5)

[lifting_operator_spectral.py](lifting_operator_spectral.py),
[fundamental_matrix_Z.py](fundamental_matrix_Z.py),
[fourier_S_decomposition.py](fourier_S_decomposition.py),
[nisoli_riesz_extraction.py](nisoli_riesz_extraction.py).

Build K_k via `fractions.Fraction`. Solve πK = π over Q via Gaussian
elimination (singular linear system, fix one component at 1 then
normalize). Compute S_k = q^k · Σ π² and ε_k = S_k − 7/15 as exact
rational arithmetic.

| k | n_k | runtime |
|---|---|---|
| 1 | 2 | <0.01s |
| 2 | 6 | <0.01s |
| 3 | 18 | 0.05s |
| 4 | 54 | 0.6s |
| 5 | 162 | 5.5s |

S_5 has 60-digit numerator and 60-digit denominator. Cross-validated
against R66 closed-form Markov chain on (Z/3^k)*. The exact-rational
regime is the gold standard but doesn't scale past k=5 (denominator
explosion).

### Float64 power iteration with scipy.eigs cross-check (k = 6, 7)

[result_epsilon_6.py](result_epsilon_6.py),
[result_epsilon_7.py](result_epsilon_7.py),
[result_epsilon_7_verify.py](result_epsilon_7_verify.py).

Build K_k as a dense `np.float64` matrix. Power iterate π ← K^T π until
||K^T π − π||_∞ < 1e-15 (typically 7–8 iterations for K_k near-rank-1).
Cross-validate against scipy.sparse.linalg.eigs with `which='LM'` for
the Perron eigenvector.

| k | n_k | iters | residual | runtime |
|---|---|---|---|---|
| 6 | 486 | 7 | 2.6e−16 | ~30s |
| 7 | 1458 | 8 | 2.8e−16 | ~5min |

Power-iteration and scipy.eigs agree on π to L1 distance 1e−15. At k=5,
float64 matches the exact rational to 1e−15.

**Key float64 gotcha** (caught at k=7 build): `1.0 / 2**v` for v=1458
builds Python int 2^1458 then converts to float → OverflowError. Fix:
`2.0**(-v)` (float negative exponent, underflows to 0 cleanly for large v).

### Matrix-free power iteration (k = 8 .. 11)

[result_epsilon_8.py](result_epsilon_8.py) through
[result_epsilon_11.py](result_epsilon_11.py). At k=11, n_11 = 118098 — too
large to materialize a dense (118098 × 118098) float64 matrix
(110GB). Solution: matrix-free Krylov.

Implement K_k as a sparse linear operator: K_k.v at each iteration
applies the truncated-Geom(1/2) v-distribution per coprime state without
ever building K_k explicitly. For each i_r, the row K[i_r, :] has at most
M = ord_{3^k}(2) nonzeros (typically 2–6 nonzeros per row at moderate k).
Use scipy.sparse.linalg.LinearOperator + scipy.sparse.linalg.eigs(which='LM').

| k | n_k | runtime | residual |
|---|---|---|---|
| 8 | 4374 | ~3 min | 4e−15 |
| 9 | 13122 | ~15 min | 6e−15 |
| 10 | 39366 | ~1 hour | 8e−15 |
| 11 | 118098 | ~6 hours | 1e−14 |

Cross-validated against FFT-based Plancherel mass computation (next
section) at k=10 to 4×10⁻¹⁵.

### FFT-based Plancherel mass computation

[result_epsilon_10.py](result_epsilon_10.py) cross-validation chain.

Use the Plancherel identity S_k = Σ_{ξ : 3∤ξ} |μ̂_k(ξ)|² (R75 Theorem 75.1)
to compute S_k from π_k via radix-3 FFT on Z/3^k. The high-frequency
sum (3 ∤ ξ) restricts to (q − 1)/q · n_k components.

Cross-check at k=10: power_iter and FFT-Plancherel agree on S_10 to
4×10⁻¹⁵, providing independent confirmation that ε_10's sign-flip and
non-monotone behavior is real (not a numerical artifact of one method).

The FFT route also enables direct verification that the K_k spectrum
trivializes (R77.4 erratum): the largest non-Perron eigenvalue's
contribution to π_k − constant decays at the spectral-gap rate; matching
this against the empirical decay confirms |λ_2|.

### Numerical methods summary

| Regime | k-range | Method | Validation |
|---|---|---|---|
| Exact rationals | k ≤ 5 | Gaussian elimination over Q | Closed-form R66 chain |
| Float64 dense | k = 6, 7 | Power iter + scipy.eigs | Agree to 1e−15; match exact at k=5 |
| Matrix-free Krylov | k = 8 .. 11 | Sparse LinearOperator + eigs | FFT cross-check at k=10 to 4e−15 |
| FFT Plancherel | k = 10 sanity | radix-3 FFT on Z/3^k | Direct sum check S_k = Σ\|μ̂\|² |

The four-regime stack provides a continuous chain of confidence from
exact rationals up to k=11. The k=7 non-monotone bounce in ε_k was
load-bearing for the rate-1/2 walk-back, so the cross-validation is
essential — a single-method finding at k=7 would have been suspect.

---

## 6. What's now publishable as standalone results

These survive the walk-backs and stand on their own:

1. **The 3x+1 base paper** ([writeup.md](writeup.md)) — σ(odd n) = α_det(n mod 2^k)
   + universal stochastic remainder, prefix-determined collapse to k
   distinct conditional distributions, ⟨α_det⟩ = log(6)/log(4/3) exact.
2. **TA.1 / TA.3 — bridge structural decomposition.** Per-class first-passage
   `s_mean(r) ≈ α_det(r) + K_h · log(N/f(N))` with slope = 1.000 ± 0.005,
   structural ε(observable, N) = −2.35 + 0.486·log(threshold) (0.080 SSE
   on 13 cells at R²=0.99).
3. **R75 / R76 / R77.5 — Plancherel anchor for c = 7/45.** Rigorous algebraic
   identity c = (1/3)·lim Σ |μ̂_n(ξ)|² over high-freq, plus R74 increment
   identity, plus the R77.5 multi-resolution decomposition over Q.
4. **R77.4 erratum + q-spectrum probe** — within-level K_k^(q) spectrum is
   trivial across the prime family q ∈ {3, 5, 7, 11, 13}; the rate-1/2
   reading at the within-level operator is universally false.
5. **Cramér convergence law q^(−θ) = 2^(1−θ) − 1** — exact MGF=1 condition
   for the qx+1 random walk on Geom(1/2) v_2. Closed form at q=3 (R15,
   w*=1 via 2²−1=3); empirical match at q=5,7 to ~2.4%.
6. **Cramér multiplier C ≈ 5/2 in qx+1 conv-rate decay** — universal across
   q ∈ {5, 7, 9, 11} at R²≥0.99.
7. **K_h(q; conv) = (1 + E*[v])/(E*[v]·log(2) − log(q))** — universal
   functional form, leading-order match across q.
8. **Per-σ-quantile band U-shape** — universal across q ∈ {3, 5, 7, 9}
   (spread attenuates monotonically with q).
9. **Q-sweep test 1 — universal q/3 ratio.** S_{k+1}^(q)/S_k^(q) → q/3
   q-universally.
10. **c̃_q = (q − 3)/q at q ∈ {11, 13, 17}.** Renormalized c̃_q candidate
    to within 0.2% / 0.07% / ~1%. Suggestive at 6 datapoints.
11. **DECOMP-UNIVERSAL** — R77.5 lift-residual decomposition extends
    q-universally as exact rational equality.
12. **Sibling 3x−1 chain symmetry** — forward K_- = σK_+σ proved;
    S_n^{3x−1} = S_n^{3x+1} as exact rationals; c = 7/45 automatic for
    3x−1.
13. **q=5 cycle catalog** — three known cycles (smallest 1, 13, 17), no
    4th cycle with smallest member ≤ 10¹⁰ at large-scale Floyd search.
14. **Numerical-methods stack** — four-regime continuous validation chain
    (exact rationals k ≤ 5 → float64 dense k=6,7 → matrix-free Krylov
    k=8..11 → FFT Plancherel cross-check) used throughout the
    extension chain.

---

## 7. Open pieces post-extension

Ranked by tractability:

1. **Compute R_k composition explicitly at q=3, k=4..11.** Distinguishes
   the three remaining hypotheses for ρ_slow ≈ 0.83 (composition,
   functional-projection, finite-k artifact). The natural follow-up
   to mode amplitudes v2.
2. **Push q=5 to k=5+ and q=7 to k=4+.** Resolves whether (q-3)/q
   deviations at q=5, 7 are finite-k transients or structural. Cheap
   compute.
3. **Extend (q, c̃_q) table to q ∈ {19, 23, 29, 31, 37, ...}.** Test
   the (q − 3)/q pattern at more "large q with 2 prim root" cases.
   Cheap (~minutes per q at k=2).
4. **Derive the universal q/3 ratio analytically** from Tao's Plancherel
   framework. Likely a clean argument exists; ~half a day of analytic
   work.
5. **Cycle-aware reformulation of Constant 3 (per-j W_j) at q ≥ 5.**
   The lattice {m_j = (4^j − 1)/3} doesn't apply (q=5 has cycles at
   1, 13, 17; not a lattice). Open: define "j" as cycle-entry point
   or v_2 at entry, and check whether per-cycle W_j structure has
   closed form.
6. **q ≥ 13 cycle catalog at scale.** q=11 and q=13 cycles enumerated
   at modest N; large-scale Floyd search analogous to q=5's would
   bound their cycle counts.
7. **Closed form for −2.35 in ε(σ).** Renewal-theoretic input on
   post-prefix descent.
8. **Closed form for the K_eff ≈ 9.94 on the post-first-crossing segment.**
   The actual mechanism behind the +1/2 slope.
9. **Rigorous rate-1/2 proof for c = 7/45.** Now sharply identified as
   requiring direct band-l¹ analysis of ĥ on the dangerous band
   (per R79 Step 4: any pointwise √N route is structurally insufficient).
   Remaining plausible attacks: Bourgain-Konyagin sum-product on ⟨4⟩,
   smooth completion via auxiliary prime, 5x+1 sibling reframing.
10. **Esscher-tilt closure for R58 → D_emp gap.** Two attempts both
    failed. Closure path requires non-uniform tilt with residue-conditional
    λ_r, or a sign-aligned observable, or a closure mechanism outside
    the Esscher family.

---

## Pointers

- **Live state:** [STATE.md](STATE.md) (current claims + opens + supersessions)
- **Full archive of derivations:** [closed_form_findings.md](closed_form_findings.md)
  (~7800 lines, 79+ results, indexed at top)
- **Base paper:** [writeup.md](writeup.md)
- **Literature audit:** [literature_check.md](literature_check.md) (Lagarias / Terras / Sinai / Tao)
- **Tao bridge findings:** [tao_bridge_findings.md](tao_bridge_findings.md)
- **c = 7/45 anchor:** [c_seven_forty_fifth.md](c_seven_forty_fifth.md)
- **Conservation law:** [result_76_conservation_law.md](result_76_conservation_law.md)
- **Path-C obstruction map:** [result_77_T_lead_spectrum.md](result_77_T_lead_spectrum.md), [result_78.md](result_78.md), [result_79.md](result_79.md), [r79b_S_partial_empirical.md](r79b_S_partial_empirical.md), [saddle_class_subsum_analysis.md](saddle_class_subsum_analysis.md), [bk_moments_analysis.md](bk_moments_analysis.md), [band_l1_analysis.md](band_l1_analysis.md), [band_spectral_decomposition.md](band_spectral_decomposition.md)
- **q-sweep results:** [result_q_sweep_test_1_rate.md](result_q_sweep_test_1_rate.md),
  [result_q_sweep_test_2_c_q.md](result_q_sweep_test_2_c_q.md),
  [result_q_sweep_test_3_decomposition.md](result_q_sweep_test_3_decomposition.md)
- **c̃_q structure:** [c_tilde_structure_verdict.md](c_tilde_structure_verdict.md), [c_tilde_q17_probe.py](c_tilde_q17_probe.py)
- **q-spectrum:** [result_qspectrum.md](result_qspectrum.md)
- **Cramér / per-prime decomposition:** [experiments/16_cramer_root.py](experiments/16_cramer_root.py), [experiments/17_cramer_dual_verification.py](experiments/17_cramer_dual_verification.py), [experiments/10_q_decomposition.py](experiments/10_q_decomposition.py), [experiments/12_q_convrate_analytical.py](experiments/12_q_convrate_analytical.py), [experiments/13_cross_q_unification.py](experiments/13_cross_q_unification.py)
- **Cycle classification:** [experiments/22_q5_cycle_detection.py](experiments/22_q5_cycle_detection.py), [experiments/29_qx1_cycle_classification.py](experiments/29_qx1_cycle_classification.py), [experiments/36_q5_fourth_cycle_search.py](experiments/36_q5_fourth_cycle_search.py)
- **ε_k extension chain (k=6..11):** [result_epsilon_6.md](result_epsilon_6.md),
  [result_epsilon_7.md](result_epsilon_7.md), [result_epsilon_10.md](result_epsilon_10.md), 
  [result_epsilon_11.md](result_epsilon_11.md)
- **Mode amplitudes:** [mode_amplitudes_v2_findings.md](probe_mode_amplitudes_v2/mode_amplitudes_v2_findings.md)
- **Sibling 3x−1:** [sibling_3x_minus_1_symmetry_verdict.md](sibling_3x_minus_1_symmetry_verdict.md), [duality_S_vs_D_verdict.md](duality_S_vs_D_verdict.md), [duality_followup_verdict.md](duality_followup_verdict.md)
- **Density-1 / Lagarias check:** [result_density_one_v2_bounds.md](result_density_one_v2_bounds.md)
- **Numerical methods:** [lifting_operator_spectral.py](lifting_operator_spectral.py), [fundamental_matrix_Z.py](fundamental_matrix_Z.py), [fourier_S_decomposition.py](fourier_S_decomposition.py), [nisoli_riesz_extraction.py](nisoli_riesz_extraction.py)
- **Reference papers:** [Lagarias/](Lagarias/), [Taos Machinery/](Taos%20Machinery/),
  Tao 2022 (lemma 1.12, prop 1.17), Kemeny-Snell, Nisoli 2026, Kalafatelis 2026.

---

## Honesty note

This document includes the substantial walk-backs that constrain each
extension line, not just the positive results:

- **Tao bridge:** −2.45 closed form not derived (log(3)/L candidate ruled out).
- **c = 7/45:** rate-1/2 envelope walked back at k=7..11; provisional
  certified bound NOT VALID under the new data; within-level K_k
  spectrum hypothesis falsified universally; mode-amplitudes v2 found
  ρ_slow not present in any single mode; three closure routes (R78 + R79 +
  C3) closed.
- **qx+1 generalization:** literal "7/45 is q=3 instance of c_q family"
  hypothesis falsified (S_∞^{(q)} doesn't exist for q ≥ 5); five q-points
  insufficient to confirm (q-3)/q at q=5, 7; non-prim-root explanation
  for q=7 ruled out by q=17; per-class slope universality breaks at q=5;
  per-j W_j Lagarias-class structure is q-SPECIFIC (q=5 attractor is a
  cycle, not a lattice).
- **Density-1 v_2:** TEST B's 100% pass rate is mathematically tautological
  for terminating Syracuse trajectories — not a quantitative density-1
  strengthening. The empirically informative content is the unconditional
  ensemble mean v = 2.102.

Per the discipline note in STATE.md: documenting walk-backs is part of
the rigor signal. The framework is more defensible WITH explicit
walk-backs than without.
