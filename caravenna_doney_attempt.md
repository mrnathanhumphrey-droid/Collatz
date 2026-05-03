# Caravenna-Doney + Borovkov local-LD attack on W_j

**Date:** 2026-05-02. Sequel to Results 15 (Path C), 16 (Esscher duality), 17 (Path B framework setup).

This document is the parallel attack on the conditional Wald-overshoot constants W_j via Caravenna-Doney 2019 (`references/CaravennaDoney2019_local_large_deviations.pdf`) plus Borovkov-style adaptation for the Gaussian-domain Syracuse log-walk. Numerical: `cd_numerical.py`.

---

## 1. Framework adaptation: CD α ∈ (0,1) ∪ (1,2) → Gaussian-domain α = 2

CD's main theorems explicitly require domain-of-attraction index α ∈ (0,1) ∪ (1,2):

- **Theorem 1.1 (LLD):** for F satisfying (1.2) with α ∈ (0,1) ∪ (1,2), bound P(S_n ∈ x + J, M_n ≤ γx) ≤ C₀ · (n/(a_n · A(x)))^⌈1/γ⌉. The Cauchy α=1 case is left out (handled by Berger 2017).
- **Theorem 1.12 (SRT for random walks):** α ∈ (0,1) only, with positivity index ρ ∈ (0,1).

The Syracuse log-walk has step X = log(3 + 1/m) − v · log(2) with v ~ Geom(1/2). This walk has **finite all-order moments** (Geom(1/2) MGF E[e^(sv)] = (1/2)/(1 − (1/2)e^s) is finite for s < log 2). Therefore X is in the domain of attraction of a normal distribution (α = 2), **explicitly outside CD's main theorems**.

**Borovkov adaptation.** For α = 2 (Gaussian domain), the analogous local-LD machinery comes from Stone (1965) for the local CLT and Borovkov & Mogulskii (2008, *Asymptotic Analysis of Random Walks: Heavy-Tailed Distributions*) for the renewal-theoretic refinements. The relevant statements:

1. **Stone's local renewal theorem (Gaussian case):** For a non-arithmetic random walk with finite variance and negative drift, the renewal density of the descending ladder process u_↓(y) satisfies u_↓(y) → 1/μ_L as y → ∞ (Blackwell, Stone). The convergence is *exponentially fast* if F has finite exponential moments on a neighborhood of zero.

2. **Local large deviation upper bound (Borovkov):** for finite-variance walks, P(S_n ∈ x + J) is governed by a Gaussian local CLT with sub-asymptotic corrections that are exponentially small in x²/n (large-deviation regime) or O(1/√n) corrections (CLT regime).

3. **Asymptotic conditional first-passage overshoot:** for IID negative-drift walk, E[overshoot at first crossing of −y] → E[L⁻²]/(2·E[L⁻]) as y → ∞ (Lorden 1970, applied to the ladder process). This is the **ladder-Lorden** formula.

The CD machinery's *spirit* — that local renewal density at a discrete target governs conditional first-passage statistics — adapts to the Gaussian case via this Borovkov framework. The α=2 case is in some sense *easier* than α ∈ (0,1) ∪ (1,2) because exponential moments give exponentially-fast convergence to asymptotic limits.

**For our problem:** W_j is the conditional Wald-overshoot residue at first-passage to the discrete absorbing target m_j ∈ {5, 85, 341, ...}. The CD/Borovkov framework predicts:

> W_j_iid → E[L⁻²]/(2·E[L⁻]·log(4/3)) as log(m_j) → ∞, **independent of j**.

This is the load-bearing prediction we test.

---

## 2. The Lorden formula: ladder vs walk-onestep

**Substantive correction to Path C / Path A.** Path C (Result 15) and Path A (`wh_numerical_check.py`) reported "Wald-iid Lorden = 6.305 step units" computed via E[X²]/(2·|E[X]|). For our Geom(1/2) Syracuse log-walk:

- E[X] = log(3) − 2·log(2) = −log(4/3) = −0.2877 nats
- E[X²] = log²(3) − 2·log(3)·log(2)·E[v] + log²(2)·E[v²] = 1.043 nats² (using E[v]=2, E[v²]=6)
- **E[X²]/(2·|E[X]|) = 1.814 nats = 6.305 step units**  (one-step / X-renewal formula)

**This is the wrong formula for first-passage overshoot.**

For a random walk with negative drift, the asymptotic conditional overshoot at first crossing of level −y (as y → ∞) is given by the **renewal process whose interarrivals are the descending ladder heights L⁻**, not the walk's one-step increments X. The correct formula is:

> E[overshoot at first crossing of −y, y → ∞] = E[L⁻²]/(2·E[L⁻])  (ladder-Lorden)

For our walk (1M-orbit Path A re-simulation in `cd_numerical.py`):
- E[L⁻] = 1.005 nats = 3.49 step units
- E[L⁻²] = 1.992 nats²
- **E[L⁻²]/(2·E[L⁻]) = 0.991 nats = 3.45 step units**  (correct ladder formula)

Direct simulation in Step 4 of `cd_numerical.py` (500K-orbit IID Syracuse walk descending from log(2³⁶) = 25.0 nats, measured first-crossing overshoot at log(m_j)) confirms:
- Mean overshoot at log(5)   = 0.996 nats ± 0.005 = 3.46 step units
- Mean overshoot at log(85)  = 0.989 nats ± 0.005 = 3.44 step units
- Mean overshoot at log(341) = 0.993 nats ± 0.005 = 3.45 step units

All three within sampling SE of the ladder-Lorden value 0.991 nats. **The correct i.i.d. baseline is 3.45 step units, not 6.305.** The 6.305 value Path C reported is the asymptotic residual life of the X-renewal process — a distinct quantity that doesn't apply here.

**Implication:** the +0.85 step-unit "conditional-on-target Markov correction" Path C identified for W_2 (= W_2_emp − 6.305 = 7.156 − 6.305 = 0.85) was based on the wrong baseline. The correct Markov correction is W_2_emp − 3.45 = **+3.71 step units**, four times larger.

---

## 3. Borovkov finite-y correction to ladder Lorden

The Lorden value is the y → ∞ asymptote. At finite y, there's a sub-asymptotic correction Δ(y) → 0 governed by the renewal equation:

> m(y) = E[(L − y)·1_{L > y}] + ∫₀^y m(y − L) dF_L(L)

where m(y) := E[overshoot at first crossing of level y]. As y → ∞, m(y) → Lorden.

**Numerical solution** (`cd_numerical.py` Step 3): discretize on uniform grid Δy = 0.05 nats, use empirical L⁻ histogram from 1M-orbit Path A simulation, iterate the renewal equation forward.

Results at the j-class targets:

| j | m_j | log(m_j) (nats) | log(m_j) / μ_L | m(log m_j) (nats) | m(log m_j) (steps) | Δ vs Lorden (steps) |
|---|---|---|---|---|---|---|
| 2 | 5   | 1.609 | 1.60 | 1.011 | 3.515 | +0.070 |
| 4 | 85  | 4.443 | 4.42 | 1.024 | 3.561 | +0.116 |
| 5 | 341 | 5.832 | 5.80 | 1.020 | 3.547 | +0.102 |

The finite-y correction Δ(log m_j) is **at most +0.12 step units** across all j tested. Even at log(m_2) = 1.6 μ_L (just above the asymptotic regime), Δ is only +0.07 step units. log(m_4), log(m_5) > 4 μ_L are deeply in the asymptotic regime.

**Effective i.i.d. prediction:**

> **W_j_iid ≈ 3.45 step units for all j ∈ {2, 4, 5}** (Lorden + Borovkov correction; j-independent within ±0.1 step units)

---

## 4. Comparison to empirical W_j and Markov correction

Empirical W_j (50M-orbit measurement at N=2³⁶, from compute_threads_findings.md):

| j | W_j_empirical | W_j_iid (CD-Borovkov) | Markov correction | |Markov| / W_iid |
|---|---|---|---|---|
| 2 | **+7.156** ± 0.006 | +3.46 | **+3.69** step units | 107% |
| 4 | **−4.755** ± 0.06  | +3.44 | **−8.20** step units | 239% |
| 5 | **+4.590** ± 0.06  | +3.45 | **+1.14** step units | 33% |

**Per the brief's decision criteria:**

- "Caravenna-Doney + Borovkov closed form for W_j matches empirical to ±0.05" → **NO.** All three j miss by 1+ step units.
- "Closed form matches W_2 within ±0.1 but misses W_4, W_5 by 1+" → **NO.** W_2 misses by 3.69, not within 0.1.
- "Closed form gives a structural framework for W_j but residual gap is large across all j" → **YES.**

**Verdict:** the i.i.d. local-LD framework gives a clean closed-form baseline (W_j_iid = 3.45 step units, j-independent), but the empirical W_j has substantial Markov-modulation that the i.i.d. framework cannot capture. Markov is structurally necessary at all three j tested.

**Why j=2 is not "captured by i.i.d." despite being dominant.** The brief anticipated that the dominant class j=2 (P = 0.938) might be largely i.i.d.-explained because the "selection bias" of conditioning on landing at the most-common target should be minimal. Empirically this is not the case: W_2 = 7.16 vs i.i.d. = 3.46, ratio 2.07×. The Markov selection effect on j=2 is small in *probability* (Q matrix mixes to uniform-on-32-residues fast, per `path_b_derivation.md`) but large in *Wald residue*. Conditional on landing at m_2 = 5, the orbit's path-properties differ substantially from a typical IID walk's path.

---

## 5. ε_S decomposition

ε_S = Σ_j P(j) · [W_j − log(m_j)/log(4/3) + 1]   (in step units)

Computed from the i.i.d.-baseline-only W_j and from empirical W_j separately:

| component | i.i.d.-only contribution | empirical contribution |
|---|---|---|
| j=2 (P=0.938) | −1.061 | +2.402 |
| j=4 (P=0.024) | −0.261 | −0.455 |
| j=5 (P=0.038) | −0.600 | −0.556 |
| **ε_S total** | **−1.92**  | **+1.39**  |

Per CTF empirical: ε_S asymptote = 1.375 ± 0.005 (sign positive). Empirical from this calculation: 1.39 (matches CTF within 0.02; rounding).

**The i.i.d. baseline alone gives ε_S = −1.92 (negative)**, far from log(4) = 1.386 ≈ empirical. The Markov contribution to ε_S is **+3.31 step units** (= empirical 1.39 − i.i.d. baseline (−1.92)) and effectively delivers the empirical value's proximity to log(4).

**Quantitative target for Path B:** the Markov-modulated matrix Wiener-Hopf framework on the residue chain Q must produce, when aggregated over j with P(j) weights, a per-j Markov correction summing to **+3.31 step units = +0.953 nats** to bring ε_S from −1.92 (i.i.d.) to the empirical +1.39. The per-j Markov corrections that need to be reproduced:

| j | P(j) | Markov W_j (step units) | P(j) × Markov (step units) |
|---|---|---|---|
| 2 | 0.938 | +3.69 | +3.46 |
| 4 | 0.024 | −8.20 | −0.20 |
| 5 | 0.038 | +1.14 | +0.04 |
| **sum** | | | **+3.31** |

The dominant Markov contribution is from j=2 (94% of the sum). **Path B's primary target is closing the +3.69 step-unit Markov correction for W_2.** The j=4 and j=5 corrections, while large in magnitude per orbit, contribute negligibly to ε_S because P(j) is small.

---

## 6. Where the framework succeeds

1. **Universal i.i.d. baseline.** The Lorden ladder formula gives W_j_iid = 3.45 step units **independent of j**, with finite-y Borovkov correction < 0.12 step units. This is a clean closed form: W_j_iid = E[L⁻²]/(2·E[L⁻]·log(4/3)) where E[L⁻] and E[L⁻²] are themselves not closed-form (Path C falsified all natural-constant candidates for E[L⁻]) but are operational empirical constants from Path A simulation.

2. **Quantification of Markov correction.** The Markov correction = empirical W_j − i.i.d. baseline is a sharply defined target for Path B's matrix WH, broken down by j and weighted by P(j) for the ε_S sum.

3. **Identification of the dominant Markov channel.** P(j=2) × Markov_2 = 3.46 step units accounts for 94% of the Markov contribution to ε_S. Path B's success on j=2 alone would close most of ε_S.

4. **Substantive math correction.** Path C's "Lorden = 6.305 step units" baseline was based on the wrong renewal formula (one-step X-residual life vs ladder-overshoot). The correct baseline (3.45 step units) makes the Markov correction four times larger than previously estimated and reframes the physical picture: Markov modulation, not "small conditional-on-target adjustment," is the dominant effect setting W_2.

---

## 7. Where the framework fails

1. **No cross-class differentiation.** i.i.d. local-LD predicts W_j_iid ≈ 3.45 for all j. Empirically W_j varies wildly (range −4.76 to +7.16). The framework is structurally blind to j-dependence; the j-dependence lives entirely in the Markov chain.

2. **No closed form for E[L⁻] itself.** Path C falsified all natural-constant candidates for E[L⁻] at 10⁷-orbit precision (Result 15). The Lorden value 0.991 nats is therefore semi-closed — i.i.d.-derivable from a non-closed-form E[L⁻] and E[L⁻²]. To get a fully closed-form W_j_iid we would need closed forms for both ladder moments, which the same Gelfond-Schneider obstruction (irrational log₂ 3 truncation in Sparre-Andersen) blocks.

3. **No handle on the Markov correction beyond i.i.d. − empirical decomposition.** Quantifying the +3.31 step-unit Markov contribution to ε_S is a description, not a derivation. The Markov correction itself requires the matrix WH framework of Path B (Result 17 partial).

---

## 8. Verdict

**Brief's third outcome:** "Closed form gives a structural framework for W_j but residual gap is large across all j; Markov is structurally necessary; Caravenna-Doney provides the i.i.d. baseline for Path B to extend." — **CONFIRMED.**

The CD-Borovkov framework delivers:
- a clean i.i.d. baseline W_j_iid = 3.45 step units (j-independent)
- a corrected Lorden formula (3.45, not Path C's 6.305)
- a sharp quantitative Markov-correction target for Path B (+3.31 step units in ε_S, 94% from j=2)

It does NOT deliver:
- closed-form W_j matching empirical at any j
- closed-form E[L⁻] (still open)
- the j-dependence of W_j (purely Markov phenomenon, requires Path B)

**Two-framework status (CD parallel + Path B):**
- CD: i.i.d. baseline closed (this document, Result 20)
- Path B: matrix WH framework set up (Result 17), full W_j extraction queued

If Path B succeeds at recovering the +3.31 Markov contribution from the residue chain, the two frameworks together close ε_S in semi-closed form (i.i.d. baseline + Markov correction). If Path B also succeeds, ε_S closure reduces to closing E[L⁻] (Path C / Esscher route, currently blocked by Gelfond-Schneider).

---

## Files

- `cd_numerical.py` — numerical verification (1M L⁻ simulation, Borovkov renewal-equation finite-y correction, 500K direct first-crossing simulation, Markov-correction quantification)
- `caravenna_doney_attempt.md` — this document
- `closed_form_findings.md` — Result 20 entry

## Citations

- Caravenna F. & Doney R., "Local large deviations and the strong renewal theorem," EJP 24 (2019) — the framework paper, α ∈ (0,1) ∪ (1,2) main theorems
- Lorden G., "On excess over the boundary," Ann. Math. Statist. 41 (1970) — asymptotic overshoot formula via ladder process
- Stone C., "On local and ratio limit theorems," Proc. 5th Berkeley Symp. (1965) — local CLT for Gaussian-domain
- Borovkov A. & Mogulskii A., *Asymptotic Analysis of Random Walks*, Cambridge (2008) — Gaussian-domain renewal-theoretic adaptations
- Berger Q., "Strong renewal theorems and local large deviations for multivariate random walks and renewals," EJP 24 (2019) — Cauchy α = 1 extension
- Asmussen S., *Applied Probability and Queues* Ch IV — descending ladder process formalism

## Honest scope statement

The brief proposed Caravenna-Doney 2019 as the framework. CD's main theorems (Thm 1.1 LLD, Thm 1.12 SRT) are for stable α ∈ (0,1) ∪ (1,2), which excludes our α=2 walk. The adaptation via Borovkov (or equivalently, classical Stone-LCLT + Cramér-Lundberg) is what actually delivers the Lorden ladder formula in our setting. CD's *spirit* (renewal density at discrete target governs conditional first-passage statistics) is what carries through; their specific theorems do not directly apply.

The substantive output is therefore: (a) the Lorden-ladder correction to Path C's mistaken 6.305 baseline, (b) the j-independent universal i.i.d. baseline 3.45 step units, (c) the quantification of the per-j Markov correction targeted by Path B. The closed-form W_j the brief sought is *partial* — only the i.i.d. baseline, not the empirical value.
