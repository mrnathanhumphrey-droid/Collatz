# Wiener-Hopf attempt for W_j on the Syracuse log-walk

**Date:** 2026-05-02. Bibliography read: Vidmar 2015, Kuznetsov-Kyprianou-Pardo 2010, Alsmeyer-Buckmann 2018, Denisov-Wachtel 2026 (intros + key definitions).

This document delivers Deliverable 1 (classification) of the Wiener-Hopf attempt brief. Deliverables 2 and 3 are not yet attempted; per the brief, classification ships first.

---

## Setting: the Syracuse log-walk

Steps:
> X_t = log(3 + 1/m_t) − v_t · log(2),    v_t = ν_2(3 m_t + 1)

In the large-m limit, X_t ≈ log(3) − v_t · log(2). Under iid Geom(1/2) approximation for v:

> E[X] = log(3) − 2 log(2) = log(3/4) ≈ −0.288 < 0
> Var[X] = log(2)² · Var[v] = log(2)² · 2 ≈ 0.961
>
> P(X = log(3) − k · log(2)) = 2^(−k),  k = 1, 2, 3, ...

Step distribution support:
- Single positive value: X = log(3/2) ≈ +0.405 (when v = 1, P = 1/2)
- Countable negative values: X = log(3) − k log(2) for k ≥ 2 (P = 2^(−k))
  - X = log(3/4) ≈ −0.288  (k=2, P=1/4)
  - X = log(3/8) ≈ −0.981  (k=3, P=1/8)
  - X = log(3/16) ≈ −1.674 (k=4, P=1/16)
  - …

**Characteristic function** (closed form in the iid Geom(1/2) approximation, valid in the strip Im(θ) < 1):

> φ(θ) = E[e^(iθX)] = 3^(iθ) / (2 · 2^(iθ) − 1)

This is **rational in 2^(iθ)** (denominator-degree-1 in 2^(iθ)) but transcendental in θ.

**The first-passage problem we want:** start at log(m_start), descend (mean drift negative), first-hit a small target log(m_j) where m_j = (4^j − 1)/3 ∈ {5, 85, 341, 1365, ...}. The Wald-Lorden boundary residue at the entry lattice is W_j.

---

## Classification by candidate framework

### 1a. Oscillating walk (Denisov-Wachtel 2026, main case)

**DOES NOT FIT.**

Denisov-Wachtel define oscillating walks by `lim sup S(n) = ∞ AND lim inf S(n) = −∞ a.s.` (DW p.1). Our walk has E[X] = log(3/4) < 0, so by SLLN S(n)/n → log(3/4), hence lim S(n) = −∞ a.s. and lim sup S(n) is finite a.s. **NOT oscillating; descending.**

This puts us outside the main DW results which are tailored to the oscillating case.

DW also covers Wiener-Hopf factorisation in their Section 3 (for general iid walks, not only oscillating); that machinery is generic and applies. But their *universality method* (the more robust modern technique) is specifically for the oscillating + scaling-to-Brownian-motion case, which we don't have.

**Verdict:** main DW machinery doesn't fit. Their general WH-factorisation review *does* apply but is just a pointer to Spitzer/Kyprianou.

---

### 1b. Upwards skip-free Lévy chain (Vidmar 2015)

**DOES NOT FIT.**

Vidmar Definition 2.1: X is an upwards skip-free Lévy chain if (i) X is a compound Poisson process, AND (ii) for some h > 0, supp(ν) ⊂ ℤh = {hk : k ∈ ℤ}, AND (iii) supp(ν|_{(0,∞)}) = {h}.

Apply to our walk:
- (i) Discrete-time random walk, not compound Poisson — but the discrete-time analog (right-continuous random walk on ℤh) is what Vidmar 2015 references via [25] (Avram-Kyprianou) and works "by analogy". So this is not the obstacle.
- (iii) Single positive jump value at +log(3/2) ✓ — this is satisfied
- (ii) **Lattice obstruction.** Negative jumps lie at sizes log(3) − k·log(2) for k ≥ 2:
  - log(3) − 2 log(2) = log(3/4)
  - log(3) − 3 log(2) = log(3/8)
  - …
  Differences between consecutive negative jumps: log(2). For all jumps to lie on a single lattice ℤh with h = log(3/2), we would need both log(3/2) and log(2) to be integer multiples of h. Then log(2)/log(3/2) ∈ ℚ. But:
  > log(2) / log(3/2) = log(2) / (log(3) − log(2)) = 1 / (log_2(3) − 1)

  Since log_2(3) is irrational (Gelfond-Schneider, equivalently, 2 and 3 are multiplicatively independent), log(2)/log(3/2) is irrational. **No lattice contains both jump sets.**

**The walk has the "non-random overshoot upward" property** (every upward jump is exactly +log(3/2), so first-passage-upward overshoot is bounded and computable directly) but lacks the lattice structure that Vidmar's machinery requires.

Furthermore, **we need DOWNWARD first-passage** (to small targets m_j on the log-scale). Even if Vidmar applied, his framework treats *upward* exit. The natural dual walk X' = −X has unrestricted positive jumps (the original negative jumps reflected) and a single negative jump value at −log(3/2). This is "downward skip-free", but Vidmar's results are for "upward skip-free". The dual-walk framing is partial: downward skip-free IS its own Vidmar-style class (Avram-Kyprianou treat the spectrally negative case [3]) — but again only on lattices.

**Verdict:** the strong "single positive jump value" property is real and useful, but the lattice obstruction blocks Vidmar's scale-function recursion (his Eq. 4.10–4.11) from applying directly. Some adaptation could be possible (lifting his linear recursion to the continuous-state "non-random overshoot" Lévy class via Bertoin VII), but is not a direct citation.

---

### 1c. Meromorphic Lévy class (Kuznetsov-Kyprianou-Pardo 2010)

**DOES NOT FIT** the strict M-class definition. Suggestive structural feature, no direct theorem.

KKP Theorem 1 characterises M-class processes by equivalent conditions including:

> (vi) The Laplace exponent ψ(z) is a *real meromorphic function* satisfying Im(ψ(z)/z) > 0 for Im(z) > 0.

For our walk, the "discrete-time Laplace exponent" is
> κ(z) := −log E[e^(zX)] = −log(3^z / (2 · 2^z − 1)) = log(2 · 2^z − 1) − z log 3

This involves `log(2 · 2^z − 1)`, which is **not meromorphic in z** — the inner expression `2 · 2^z − 1` is entire in z, but log of an entire function with zeros is NOT meromorphic (it has logarithmic branch points at the zeros). So κ(z) has logarithmic branch singularities at z values satisfying 2 · 2^z = 1, i.e., z = −1 + 2πi k / log 2 for integer k.

**However:** the *characteristic function itself*, φ(θ) = 3^(iθ) / (2 · 2^(iθ) − 1), is meromorphic in θ within the strip Im(θ) < 1 (with poles on Im(θ) = 1). So the **probability generating function is meromorphic in the variable 2^(iθ)** — a substitution to an "exponential of θ" variable.

This is a softer property than KKP's M-class (which requires meromorphic Laplace exponent). KKP's machinery — partial-fraction decomposition into pole/root sequences, Eq. (6) — relies on the rational structure of `q − ψ(z)` in z. Our `q − κ(z)` is transcendental in z.

**No direct theorem in KKP applies.** A generalization to the "meromorphic in exp(z)" case might exist in the literature (this is related to "lattice" or "exponential-functional" Lévy processes), but it is not in KKP 2010 as written.

**Verdict:** strict M-class doesn't apply. The rational-in-2^(iθ) structure is suggestive — there *may be* an analogous theory for processes with rational-in-exp(z) characteristic function, related to "exponential Lévy processes" — but I don't have a citation in this bibliography for it. KKP's residue-analysis style would have to be adapted, not directly invoked.

---

### 1d. Markov-modulated random walk (Alsmeyer-Buckmann 2018)

**FITS structurally.** This is the right framework if iid is the wrong baseline.

Alsmeyer-Buckmann set up: a Markov chain M = (M_n) on a countable state space S with transition matrix P and stationary distribution π, and a Markov-modulated walk
> S_n = X_1 + ... + X_n,    P((X_1, ..., X_n) ∈ · | M_0 = i_0, ..., M_n = i_n) = K_{i_0 i_1} ⊗ ... ⊗ K_{i_{n−1} i_n}

Empirical confirmation that we need this:
- v_t = ν_2(3 m_t + 1) under the trajectory measure has Cov[v_0, v_k] ≈ −0.4 at lags 50–100 (Result 8 / exp 43) — *non-iid by ~3× the variance*.
- E[v]_traj = 1.995 ≠ 2.000 (the iid Geom(1/2) value).

For our setting, the natural driving chain is the orbit value m_t mod some modulus (or an equivalent Markov decomposition of the Syracuse map's residue dynamics). The compute_threads_findings.md absorbing-Markov-chain analysis already constructs such a chain on m ∈ [3, M] with deterministic transitions — that's the M_n process; the X_t are the corresponding log-step increments.

**A-B's Theorem framework** addresses fluctuation-type trichotomy (positive divergence / negative divergence / oscillation) for MRWs, with explicit conditions in terms of moments under the embedded ordinary-RW restricted to return times of M to a fixed state.

**Wiener-Hopf for Markov-additive processes:** A-B reference Asmussen Ch XI, Prabhu et al. [44], and others for the matrix-Wiener-Hopf factorization in the MRW setting. The factorization extends, but with matrix-valued objects (the state-space dimension introduces matrix structure on the WH factors).

**For W_j specifically:** the absorbing-Markov chain machinery in compute_threads_findings.md already captures the GEOMETRIC invariant P(j) exactly (P(j=2) = 0.938 to ±0.005). What it FAILS to capture is the METRIC invariant W_j — by 1.0+ unit. The A-B framework gives a route to the missing metric piece: the *matrix-WH factorization* of the MRW gives the conditional joint law of (τ_j, S(τ_j) | M absorbs at j), which determines W_j directly.

**Verdict:** A-B is the **correct theoretical framework**. The challenge is operational — identifying the driving chain at the right granularity and computing the matrix-WH factors explicitly. Section 4 of A-B (ladder variables) and the references [12, 24, 44] for matrix-WH would be the next reads.

---

### 1e. Generic Spitzer-Baxter (Spitzer 1956)

**FITS in the iid approximation.** Universal but doesn't immediately give closed form.

Spitzer's identity: for any iid random walk and |s| < 1,
> exp(Σ_{n≥1} (s^n / n) · E[e^(iθ S_n) · 1_{S_n > 0}])
>   = E[s^(τ⁺) · e^(iθ S(τ⁺))]   where τ⁺ = inf{n ≥ 1 : S_n > 0}

This decomposes 1 − s · E[e^(iθX)] = (1 − s · h⁺(θ)) · (1 − s · h⁻(θ)) into ascending and descending ladder factors. For our walk, the ascending factor h⁺ is conditioned on the (rare) positive-jump events with single value +log(3/2); the descending factor h⁻ is the dominant component.

**Closed-form W_j from S-B:** in principle, the joint distribution of (τ_descent, S(τ_descent)) at first downcrossing of any level x < 0 is determined by h⁻(θ). For our walk the descending ladder has a clean structure because of the negative jump law:
> P(S(τ⁻) ≤ −y | τ⁻ < ∞, no positive intermediate excursion) determined by overshoot at first crossing

The Lorden-style upper bound for first-passage to a half-line (W_j_iid_Lorden ≈ E[X²] / (2|E[X]|) ≈ 1.67 in our walk) is universal but doesn't capture the conditional-on-hitting-specific-target boundary residue.

**Verdict:** Spitzer-Baxter gives the universal apparatus and an iid baseline. It cannot deliver closed-form W_j for first-passage to a *specific small target m_j on a discrete set* without additional structure — the target set {5, 85, 341, ...} is determined by Collatz number theory, not random-walk geometry, and Spitzer-Baxter doesn't see this structure.

---

### 1f. The right combined framework

The two operationally-relevant pieces:

1. **Spitzer-Baxter (iid baseline)** for `W_j_iid_first-passage-to-half-line`. Gives the universal Wald-Lorden constant ≈ 1.67. Empirical W_2 = 7.156 is **far from** this baseline; the gap is dominated by:
   (a) **Conditional-on-hitting-m_j**: lattice-target restriction multiplies the residue substantially.
   (b) **Markov dependence**: empirical autocorrelation `Cov[v, v_lag]` integrates non-trivially.

2. **Alsmeyer-Buckmann (Markov correction)** for converting (a) and (b) into the correct empirical W_j. The matrix-Wiener-Hopf factorization on the driving-chain state space gives the discrete-target conditional distribution explicitly.

The KKP and Vidmar frameworks are **structurally aligned** (both produce explicit WH factors for processes with special-structure characteristic functions) but **don't fit our walk's hypotheses strictly**:
- Vidmar's lattice obstruction is real (jumps in log(3)ℤ + log(2)ℤ, dense in ℝ).
- KKP's meromorphicity-in-θ doesn't hold (log of entire function is not meromorphic).

The Syracuse log-walk's characteristic function being **rational in 2^(iθ)** is a striking structural feature that doesn't fit any of the named frameworks directly. Whether this admits an *adapted* meromorphic-in-exp(z) framework (perhaps in a generalized KKP class) is open from this bibliography alone.

---

## Deliverable 1 verdict

| Framework | Fits? | Key obstruction (if any) | Path to W_j |
|---|---|---|---|
| Oscillating (DW main results) | NO | descending walk, not oscillating | n/a |
| Upwards skip-free (Vidmar) | NO | lattice obstruction (log(2)/log(3/2) irrational) | not directly; requires continuous-state generalization |
| Meromorphic Lévy (KKP) | NO | char fn not meromorphic in θ (only in 2^(iθ)) | requires generalization to "rational in exp(z)" class — not in this bibliography |
| Markov-modulated (Alsmeyer-Buckmann) | **YES** | requires identifying driving chain | matrix-WH factorization for MRW |
| Generic Spitzer-Baxter | YES (iid baseline) | doesn't see lattice target {m_j} | Wald-Lorden + corrections |

**Recommendation for Deliverable 2:** combine Spitzer-Baxter (iid Geom(1/2) baseline, exact computation feasible because φ(θ) is rational in 2^(iθ)) with Alsmeyer-Buckmann's matrix-WH factorization (Markov dependence correction). The first half should give a concrete numerical W_2_iid baseline, and the gap to empirical 7.156 will quantify how much of the residue is iid-WH versus how much requires Markov correction.

The non-trivial structural insight — **rational-in-2^(iθ) characteristic function** — is suggestive and should be explored: it may give *exact closed-form* Spitzer-Baxter ladder factors via explicit summation of geometric series in 2^(iθ), avoiding the need for numerical contour integration. This is the most promising route in the bibliography to a genuinely closed-form W_j.

**Honest scope statement:** none of the cited frameworks gives a "drop-in" closed-form W_j. The bibliography supplies the right machinery (Spitzer-Baxter for iid, Alsmeyer-Buckmann for Markov modulation) but the application requires non-trivial work to (a) compute the iid baseline using the rational-in-2^(iθ) structure, (b) identify the right driving chain for Markov modulation, (c) carry out matrix-WH factorization. Each is a multi-page calculation; not a 1-line citation.

---

## Deliverable 1 status

Deliverable 1 complete. User directed: take Path A.

---

# Deliverable 2: iid Wald-Lorden baseline via Spitzer-Baxter

## 2.1 The Spitzer-Baxter / Sparre-Andersen identities for the descending ladder

For an iid random walk S_n with E[X] < 0, the strict descending ladder time σ⁻ = inf{n ≥ 1 : S_n < 0} is a.s. finite, and the strict descending ladder height L⁻ = −S(σ⁻) > 0 has the Wiener-Hopf representation

> 1 − ∫ e^(iθ x) F_X(dx) = (1 − κ⁺(θ)) · (1 − κ⁻(θ))

where the factors are the characteristic functions of (a) the strict ascending ladder height with killing (since E[X] < 0, ascent has positive probability of never happening), and (b) the strict descending ladder height (a.s. finite). The factors are uniquely determined by:

- κ⁻(θ) is bounded and continuous in the closed lower half-plane Im(θ) ≤ 0
- κ⁺(θ) is bounded and continuous in the closed upper half-plane Im(θ) ≥ 0
- κ⁻(0) = 1, κ⁺(0) < 1

For our walk under iid Geom(1/2):

> 1 − φ(θ) = (2 · 2^(iθ) − 1 − 3^(iθ)) / (2 · 2^(iθ) − 1)

The denominator `2 · 2^(iθ) − 1` is *rational in 2^(iθ)* — meromorphic with a single pole at 2^(iθ) = 1/2 (Im(θ) = 1). The numerator `2 · 2^(iθ) − 1 − 3^(iθ)` is **not** rational in 2^(iθ): the term 3^(iθ) = (2^(iθ))^(log_2 3) has irrational exponent log_2 3, hence no rational expansion.

This is the precise reason why the *full* Wiener-Hopf factorization isn't cleanly closed-form via the methods in the bibliography:
- KKP (meromorphic Lévy class) needs meromorphic Laplace exponent in the spectral variable — fails here
- Vidmar (skip-free) needs lattice support — fails here (irrational lattice ratio)

The denominator's rational structure DOES help on the κ⁻ side (the descending-direction analyticity is governed by the pole at Im(θ) = 1, which is the rational-in-2^(iθ) part). The numerator's irrationality affects the ascending κ⁺ side. This is structurally consistent with the asymmetry of our walk: positive jumps are at a single value (+log(3/2)), so the ascending ladder is "simple"; negative jumps span a countable set with weights 2^(−k), so the descending ladder is "rich" but governed by the rational-2-side.

## 2.2 Numerical computation of E[L⁻] (`wh_numerical_check.py`)

Direct simulation under iid Geom(1/2), 10⁶ orbits, max 20K steps each (zero failures):

> **E[L⁻]_simulated = 1.006 ± 0.002 nats = 3.497 step units** (95% CI)

Two analytic candidates from elementary moments:

| Candidate | Value (nats) | Value (step units) | Interpretation |
|---|---|---|---|
| E[X²] / (2 \|E[X]\|) | 1.814 | 6.305 | Lorden's *asymptotic residual life* upper bound |
| Var[X] / (2 \|E[X]\|) | 1.670 | 5.805 | (alternate moment formula) |
| E[X⁻] | 0.490 | 1.705 | Mean magnitude of single-step negative jump |
| **E[L⁻] (simulated)** | **1.006** | **3.497** | **Strict descending ladder mean (this work)** |

**The Lorden value 6.305 is NOT the strict descending ladder height mean.** It is the limit of E[overshoot at first crossing of −y] as y → ∞ — the *asymptotic residual life* in the ladder renewal process. The strict descending ladder height L⁻ at first crossing of 0 has a smaller mean (3.50 step units), since L⁻'s distribution has more mass near 0 than the equilibrium residual distribution.

This refines compute_threads_findings.md's reported "Wald-iid Lorden = 6.305" — that number is the Lorden upper bound, not the strict ladder mean. Both values are correct iid quantities; they measure different things.

## 2.3 Comparison to empirical W_j

The empirical W_j is structurally:

> W_j = ⟨σ_S | absorbed at j⟩ − (⟨log m_start⟩ − log(m_j))/log(4/3) − 1
>     = E[overshoot at first hit of log(m_j) | absorbed at j] / log(4/3)

This is a **conditional** quantity: conditional on the orbit absorbing exactly at m_j (a specific small target). For iid Wiener-Hopf, the analogous *unconditional* first-passage overshoot at log(m_j), starting from log(m_start) ≈ log N → ∞, asymptotes to Lorden = 1.814 nats = 6.305 step units.

| Quantity | Value (steps) | Value (nats) |
|---|---|---|
| Empirical W_2 (m=5)             | **+7.156** | +2.059 |
| Empirical W_4 (m=85)            | **−4.755** | −1.368 |
| Empirical W_5 (m=341)           | **+4.590** | +1.320 |
| Iid Lorden (asymptotic residual)| 6.305 | 1.814 |
| Iid strict ladder mean (sim)    | 3.497 | 1.006 |

**Key findings:**

1. **The iid Wiener-Hopf gives a single universal value** (whatever you take — Lorden 6.305 or strict ladder 3.497). It cannot vary by j-class. Empirical W_j varies dramatically across j (sign-changing, magnitudes 4–7 step units).

2. **Empirical W_2 EXCEEDS even the Lorden upper bound** by +0.85 step units (+0.244 nats). This is striking: the conditional overshoot at first-hitting m=5 is *larger than the asymptotic residual life of the unconditional iid walk*. Conditioning on landing at the smallest m_j selects exactly those orbits with anomalously deep negative excursions before absorption — pushing the conditional residue above the iid asymptote.

3. **W_4 = −4.76 is negative** — the orbit absorbing at m_4 = 85 has σ_S systematically *below* the linear prediction. Conditioning on landing at m_4 (a less common attractor, P=0.024) selects orbits with shallower-than-average descent at the boundary.

4. **Per-j Markov structure dominates iid Wiener-Hopf.** The full empirical residue ε_S = Σ P(j) · [W_j − log(m_j)/log(4/3) + 1] requires the *signed-and-class-specific* W_j values. iid WH cannot produce these.

Per-class breakdown of how much iid captures vs how much is Markov correction:

| j | log(m_j)/log(4/3) | Lorden W_iid (steps) | W_j_emp − Lorden | Markov correction |
|---|---|---|---|---|
| 2 (m=5)     | 5.59 | 6.305 | **+0.851** | small over Lorden |
| 4 (m=85)    | 15.45 | 6.305 | **−11.06** | large negative correction |
| 5 (m=341)   | 20.34 | 6.305 | **−1.71** | moderate negative correction |

The per-j Markov correction **dominates** the residue structure (range |correction| up to 11 step units, vs iid baseline 6 step units). This is conditional-on-target structure that requires Markov machinery (Path B) — pure iid Wiener-Hopf cannot get there.

## 2.4 Decision-criteria verdict

The brief's success criteria for Deliverable 2:
- ±0.05 step match: framework fits cleanly → success
- ±0.1 match: approximate fit → partial success
- > ±0.5 gap: framework doesn't fit, document the gap

**Result:** the iid Wiener-Hopf framework is **±0.85 step units off** for W_2, and **±11.06 off** for W_4. Both fail the ±0.5 success threshold. **iid framework documented as insufficient** for closed-form W_j.

What the iid framework DID deliver: a clean universal baseline at 6.305 step units (Lorden) or 3.497 (strict ladder), serving as the reference against which Markov corrections can be measured. The cross-class variation in W_j is now isolated as a Markov-modulation phenomenon worth ~5–11 step units per class.

---

## 2.5 Path C connection (the user's observation, 2026-05-02)

The rational-in-2^(iθ) structure of the *denominator* `2 · 2^(iθ) − 1` is exactly what makes the asymmetric pole/branch structure of φ(θ) tractable: the descending ladder factor κ⁻ is governed by this rational piece. The numerator's `3^(iθ)` term carries the irrational ratio log_2 3 — the algebraic obstruction to KKP-style meromorphic factorization.

If the closed-form program for W_j is to advance via a **Path C "exponential-meromorphic" class**, the natural research direction is exactly this asymmetry:

- Define a class of Lévy / random-walk processes whose characteristic function `φ(θ)` is rational in `e^(α · iθ)` for some α > 0, with potentially-irrational additional factors.
- For our walk, α = log(2) makes the denominator rational; the numerator's `3^(iθ)` factor is the "non-rational" piece.
- **The Wiener-Hopf factor κ⁻(θ) is determined entirely by the denominator structure** in the lower half-plane (Im(θ) ≤ 0), which IS rational in 2^(iθ) in our case.
- This suggests κ⁻(θ) may admit a closed form via residue analysis at the (single) lower-half-plane pole of `1/(2 · 2^(iθ) − 1)`, weighted by the numerator's value at that pole.

**Concrete conjecture for Path C:** the Spitzer-Baxter ladder factor κ⁻(θ) for the Syracuse log-walk is

> κ⁻(θ) = a closed-form rational expression in 2^(iθ), modulated by 3^(iθ),

derivable by the residue-style method of KKP applied to the rational-in-2^(iθ) part. If true, this gives W_iid in closed form (without need for simulation), AND extends to a structural framework for the per-j conditional W_j by combining with Markov modulation (Path B).

**Half the work is done:** Deliverable 1 + 2 have established that the rational-in-2^(iθ) structure is the right algebraic invariant; the remaining work is to (a) carry out the residue computation explicitly and (b) couple it to the matrix-WH for Markov modulation. These are concrete next steps, not hand-waving.

---

## Status

Deliverable 1 ✓ (classification: Markov-modulated is the right framework, generic Spitzer-Baxter for iid baseline; rational-in-2^(iθ) structure is the algebraic invariant)

Deliverable 2 ✓ (iid Wald-Lorden baseline 6.305 step units / strict-ladder 3.50 step units; empirical W_j 7.16 / −4.76 / +4.59; per-j Markov correction dominates iid by 0.85–11 step units per class)

**Stopping per brief instruction. Deliverable 3 (ε_S aggregation + comparison to empirical asymptote 1.375) conditional on user direction.**

Next-step options:
- **Continue to Deliverable 3** with the iid Lorden W_iid as-is (universal across j) and aggregate via Σ_j P(j) · [W_iid − log(m_j)/log(4/3) + 1]. This will reproduce a iid prediction for ε_S and compare to empirical 1.375.
- **Pivot to Path B** (Alsmeyer-Buckmann matrix-WH) to derive per-j W_j with Markov modulation.
- **Pursue Path C** (residue-style closed form for iid κ⁻ via rational-in-2^(iθ)) — clean derivation now that the structure is identified.
