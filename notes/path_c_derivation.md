# Path C: closed-form attempt for κ⁻(θ) via rational-in-2^(iθ) residue analysis

**Date:** 2026-05-02. Sequel to `wiener_hopf_attempt_v1.md` (Deliverables 1, 2).

This document attempts the closed-form derivation of the descending ladder factor κ⁻(θ) for the iid Geom(1/2) Syracuse log-walk, using the rational-in-2^(iθ) algebraic invariant identified in Path A. Numerical verification: `path_c_numerical.py` (10⁷-orbit simulation + structural calculations).

---

## 1. Setup and notation

Walk: X = log(3) − v·log(2), v ~ Geom(1/2) on {1, 2, 3, ...}.

**Characteristic function** (in the strip Im(θ) < 1):

> φ(θ) = E[e^(iθ X)] = 3^(iθ) / (2 · 2^(iθ) − 1)

In the variable u = 2^(iθ): φ(u) = 3^(iθ) / (2u − 1) = u^(log_2 3) / (2u − 1).

**Wiener-Hopf factorization** for the walk with negative drift E[X] = -log(4/3) < 0:

> 1 − φ(θ) = (1 − κ⁺_def(θ)) · (1 − κ⁻(θ))

where κ⁻ is analytic in the closed lower half-plane Im(θ) ≤ 0 with κ⁻(0) = 1, and κ⁺_def is analytic in the closed upper half-plane Im(θ) ≥ 0 with κ⁺_def(0) = q < 1 (defective, since ascent has positive probability of never happening for negative-drift walks).

The numerator of 1 − φ(θ) is

> N(θ) := 2 · 2^(iθ) − 1 − 3^(iθ)

which is entire in θ. Its zeros distribute between upper and lower half-planes; the WH factorization splits them.

## 2. The Cramer root: an algebraic identity 2² − 1 = 3

**Real zeros of N(θ).** On the imaginary axis θ = iw with w ∈ ℝ, N(iw) = 2·2^(-w) − 1 − 3^(-w) = 0 gives:

> 2^(1−w) = 1 + 3^(−w)

| w | LHS = 2^(1-w) | RHS = 1 + 3^(-w) | match |
|---|---|---|---|
| 0 | 2 | 1 + 1 = 2 | ✓ |
| -1 | 4 | 1 + 3 = 4 | ✓ |

**w = -1 is the Cramer root**, equivalently θ = -i (in the lower half-plane Im(θ) = -1). It corresponds to the **algebraic identity 2² − 1 = 3** — a Collatz-specific number-theoretic fact:

> The Cramer root w* = 1 (so that E[e^X] = 1 at w = 1) is exactly the identity that makes the Collatz step "+1 then halve" have the right algebraic balance: each unit of Cramer-MGF on log 3 is matched by 2² of MGF on log 2, with the "−1" corresponding to the discrete shift.

The Cramer root being **exactly 1** (not an irrational number) is unusual for general iid random walks. It's a structural consequence of the Collatz map's "3x+1" scaling.

**Consequence:** the walk's supremum M (for our negative-drift walk) has tail

> P(M > x) ~ C · e^(-x) as x → ∞

with **decay rate exactly 1 in nats per unit log-displacement**. This is a closed-form result, courtesy of the algebraic identity.

## 3. Esscher tilt at w* = 1 — closed-form parameters

Under the exponential change of measure dP*/dP = e^X (which is a probability measure because E[e^X] = 1 by the Cramer identity), the walk has positive drift.

**Tilted step distribution.** Under P*, P*(v = k) ∝ e^(-k·log 2) · P(v = k) = (1/2)^k · 2^(-k) = 4^(-k). Normalizing: Σ_{k≥1} 4^(-k) = 1/3, so:

> P*(v = k) = 3 · 4^(-k)   (k = 1, 2, 3, ...)

**Tilted moments** (closed form):

> E_P*[v] = 3 · Σ k · 4^(-k) = 3 · (4/9) = **4/3**
>
> Var_P*[v] = E_P*[v²] − (E_P*[v])² = 3·(20/27) − (4/3)² = 20/9 − 16/9 = **4/9**

**Tilted drift μ\*** = log(3) − E_P*[v]·log(2) = log(3) − (4/3)·log(2):

> μ* = log(3 · 2^(-4/3))   ≈ +0.1744 nats per step (positive, as expected)

These are all closed-form expressions in {log 2, log 3} via the Cramer-root identity.

## 4. Wiener-Hopf factorization attempt with rational-in-u ansatz

Goal: factor N(θ) = 2·2^(iθ) − 1 − 3^(iθ) into upper/lower half-plane parts.

Substitute u = 2^(iθ):
> N(u) = 2u − 1 − u^(log_2 3)

The exponent log_2 3 ≈ 1.585 is irrational, so u^(log_2 3) does **not** admit a rational expansion in u. The function N(u) is *not* a polynomial or rational function in u — it's a "generalized polynomial" with one term having an irrational exponent.

**Lower-half-plane zeros of N(θ).** Two real-line / imaginary-axis zeros: θ = 0 and θ = -i. The latter is in the open lower half-plane.

**Are there other zeros?** N is entire (as a function of θ ∈ ℂ), so by general entire-function theory it has either none or infinitely many complex zeros. Off the imaginary axis, the equation 2·u − u^(log_2 3) = 1 with u = 2^(iθ) has complex solutions wherever the two terms balance. The structure depends on the irrational exponent and is not amenable to simple residue analysis.

**Ansatz attempt (rational in u for κ⁻).** Try:

> 1 − κ⁻(θ) = A · (u − 1)(u − 2) / Q(u)

where Q is a polynomial with no zeros in the lower half-plane (i.e., zeros must be in u ∈ {|u| ≥ 1, arg(u) ∈ [0, π]}, corresponding to upper θ-half-plane). Normalization: as θ → -i∞ (u → ∞ along positive real), 1 − κ⁻ → 1, so Q(u)/((u-1)(u-2)) → A·u^(-(degree balance)) — for finite limit, deg Q = 2 with leading coeff A.

But this ansatz forces 1 − κ⁻(θ) to be rational in u, while the WH identity says (1 − κ⁻) · (1 − κ⁺_def) = N(u)/(2u − 1). The numerator carries the irrational power u^(log_2 3); a rational 1 − κ⁻ forces 1 − κ⁺_def to absorb the irrational structure. Then κ⁺_def is NOT rational in u, breaking the symmetry of the rational-in-u program.

**Conclusion:** the rational-in-u ansatz for κ⁻ alone is consistent but PUSHES the irrationality into κ⁺_def. It doesn't give a fully closed-form factorization in the rational-in-u algebra. The algebraic obstruction is the irrational exponent log_2 3 in u^(log_2 3) = 3^(iθ).

## 5. Numerical verification of E[L⁻] (10⁷ orbits)

High-precision simulation of the iid Geom(1/2) walk's strict descending ladder height:

| Quantity | Simulated value | 95% CI |
|---|---|---|
| E[L⁻] | **1.00456 nats** = **3.4919 step units** | ±0.00061 nats |
| E[σ⁻] | **3.4917 steps** | ±0.0033 |
| Wald check | E[σ⁻]·log(4/3) = 1.0045 vs E[L⁻] = 1.0046 | diff = -6×10⁻⁵ ✓ |
| Var[L⁻] | 0.9818 nats² | — |
| E[L⁻²]/(2·E[L⁻]) | 0.9910 nats | — |

**Wald's first identity confirmed at 10⁻⁵ precision.** E[σ⁻]·log(4/3) = E[L⁻] holds tightly.

## 6. Closed-form conjecture testing

| Candidate | Value (nats) | Diff from sim 1.00456 | Verdict |
|---|---|---|---|
| (7/2)·log(4/3) | 1.00689 | -0.00233 (-7.4σ) | **falsified** |
| log(e) = 1 | 1.00000 | +0.00456 (+14.5σ) | falsified |
| log(8/3) | 0.98083 | +0.02373 (+75σ) | falsified |
| log(2) + log(3/4)/2 | 0.54931 | +0.45525 (+1453σ) | falsified |
| log(3) | 1.09861 | -0.09406 (-300σ) | falsified |

**None of the simple closed-form candidates in {log 2, log 3, log(4/3), log(8/3), e, log(e)} matches.** The strict descending ladder height mean for the iid Syracuse log-walk does not have a simple closed form in these natural Collatz-related constants.

The (7/2)·log(4/3) conjecture (equivalently E[σ⁻] = 7/2) is closest at 0.2% off but **rejected at 7.4σ**. E[σ⁻] ≈ 3.4917, and looking for a clean form: 3.4917/log(4/3)·log(4/3) → 3.4917 isn't obviously rational or simply expressible in {log 2, log 3}.

## 7. Lorden asymptotic vs strict ladder mean — clarification

Path A's Section 2.2 noted that Lorden's formula E[X²]/(2|E[X]|) = **1.8139 nats = 6.305 steps** is the *asymptotic residual life* (limit of overshoot at first crossing of -y as y → ∞), NOT the strict descending ladder height mean E[L⁻] = 1.0046 nats.

Path C confirms this distinction. The two values differ by a factor of ~1.8 in our walk:

> Lorden / strict ladder = 1.8139 / 1.0046 ≈ 1.806

For comparison to empirical W_2 = 7.156 step units, **the Lorden value 6.305 is the relevant iid baseline** (since log m_start ≫ log m_j puts us in the asymptotic regime where Lorden applies).

## 8. Spitzer formula — mechanical and computational obstacles

The Wiener-Hopf factor 1 − κ⁻(θ) has the formal Spitzer representation:

> 1 − κ⁻(θ) = exp(−Σ_{n=1}^∞ (1/n) · E[e^(iθ S_n); S_n ≤ 0])

For our walk, S_n = n·log 3 − T_n·log 2 with T_n ~ NegBin(n, 1/2). The condition S_n ≤ 0 becomes T_n ≥ ⌈n · log_2 3⌉, which involves the **irrational threshold log_2 3** truncating at different integer values for each n.

**Why this doesn't close in closed form:**

1. The full untruncated sum Σ_{k≥n} P(T_n = k) · e^(-k·log 2 · iθ) = (some rational expression in 2^(iθ)) — closed-form thanks to NegBin generating function.

2. The truncated sum Σ_{k ≥ ⌈n·log_2 3⌉} ... lacks closed form because the truncation index is *irrational-determined* (the integer ceiling of an irrational multiple of n).

3. As n varies, the truncation rounds {n·log_2 3} mod 1 over a sequence whose distribution by Weyl's equidistribution theorem is uniform on [0, 1) — no algebraic structure in the truncation pattern.

So the Spitzer formula gives a **formal series** for κ⁻ but no closed-form summation. Numerical evaluation works (truncate the series at large N, sum), but doesn't yield an analytic expression.

This is the precise mechanical obstacle that Path C cannot surmount. The KKP-style residue analysis would be applicable IF the characteristic function were meromorphic in u = 2^(iθ). Our φ(u) = u^(log_2 3)/(2u-1) is meromorphic only modulo the irrational power u^(log_2 3); the residue at the single pole u = 1/2 doesn't capture the full WH structure.

## 9. What Path C achieves and what it doesn't

### Achievements

1. **Cramer root w* = 1 derived in closed form** via the algebraic identity 2² − 1 = 3.
2. **Esscher tilt parameters** computed in closed form: μ* = log(3·2^(-4/3)), P*(v=k) = 3·4^(-k), E_P*[v] = 4/3.
3. **Tail decay of supremum M** is e^(-x) exactly (rate 1 in nats), a closed-form consequence of w* = 1.
4. **Wald's identity** verified: E[σ⁻]·log(4/3) = E[L⁻] at 10⁻⁵ precision.
5. **Closed-form conjectures for E[L⁻] ruled out** at simulation precision (none of the natural candidates match within sampling error).

### Not achieved

1. **κ⁻(θ) as a closed-form expression** in u = 2^(iθ) and 3^(iθ). The irrational log_2 3 exponent in u^(log_2 3) = 3^(iθ) blocks rational factorization.
2. **E[L⁻] in closed form.** Empirical 1.00456 nats, no clean expression in {log 2, log 3, log(4/3), log e}.
3. **Conditional W_j formula.** The renewal measure of the descending ladder process at target log(m_j) requires κ⁻ explicitly, which we don't have.

### Decision criteria from the brief

- "Closed-form κ⁻ matches simulation E[L⁻] to ±0.5%": **not achieved** (no closed form derived).
- "Closed form gives Lorden = 6.305 exactly": vacuously fails (Lorden formula in closed form is well-known: E[X²]/(2|E[X]|), reproduces 6.305 trivially, but isn't *derived from κ⁻*).
- "Conditional W_j formula reproduces iid portion": not achieved (no closed-form κ⁻ to extract conditional).
- **"If any of these fails: the rational-in-2^(iθ) structure is incomplete":** confirmed. The structural feature missing is the irrational exponent log_2 3 carried by 3^(iθ) = u^(log_2 3), which is exactly the obstruction we identified in Path A and is now confirmed by direct attempt.

## 10. Verdict

**Path C identifies the algebraic structure (Cramer-root identity, Esscher-tilt closed forms, supremum-tail decay rate) but does not deliver a closed-form κ⁻(θ).** The rational-in-2^(iθ) algebraic invariant from Path A is *partial* — it governs the denominator structure but the numerator carries irrational exponent log_2 3 that resists rational factorization.

The Cramer-root identity 2² − 1 = 3 is itself the **load-bearing closed-form result** of Path C: it provides the exact tail decay rate of the supremum and the closed-form Esscher-tilt parameters. These are genuinely closed in terms of {log 2, log 3}.

**For the W_j program:** the iid baseline is now characterized at three levels of precision:
- Lorden formula 6.305 step units (asymptotic residual life, applies to large-log-target first-passage)
- Strict ladder mean 3.492 step units (per-descending-step structural)
- Empirical W_2 = 7.156 step units (conditional on target = m_2, requires Markov)

The **conditional-on-target Markov correction** (the gap between Lorden 6.305 and empirical 7.156) is +0.85 step units. This is the structural piece Path C cannot capture — it requires Path B (Alsmeyer-Buckmann matrix-WH).

## 11. Path B is the next move

Why: the irrational-log_2-3 obstruction in the iid framework means closed-form κ⁻ via "rational in exp(z)" alone won't work. The Markov-modulated framework, by contrast, captures the cross-class W_j variation directly through state-space structure (residues mod 2^k or absorbing-Markov chain on m). The matrix-WH factorization for MRWs gives the conditional joint law of (τ_j, S(τ_j)) at first-hitting of target class j explicitly via the chain's transition kernel.

**Concrete next step:** apply Alsmeyer-Buckmann's matrix-WH factorization to the absorbing Markov chain on odd m ∈ [3, M] with absorbing classes at m_j = (4^j − 1)/3. The driving chain is the residue process on the Syracuse map; the Markov-additive part is the log-m random walk modulated by it. The chain machinery from compute_threads_findings.md (which already derives P(j) exactly) extends to W_j via the metric structure.

Note: the bibliography's Alsmeyer-Buckmann reference is the right tool, and compute_threads_findings.md already has 50% of the chain machinery built. Path B is operationally tractable.

---

## Files

- `wh_numerical_check.py` — Path A iid baseline simulation (10⁶ orbits, brief)
- `path_c_numerical.py` — Path C high-precision (10⁷ orbits) + Cramer/Esscher analysis + closed-form conjecture tests
- `wiener_hopf_attempt_v1.md` — Deliverables 1, 2 (classification + iid baseline)
- `path_c_derivation.md` — this document (Path C closed-form attempt)

## Citations

- **Spitzer 1956**, "A combinatorial lemma" — the Spitzer-Baxter identity and its consequences for ladder factor representation. Directly used: Sparre-Andersen formula for ladder height MGF.
- **Kyprianou 2014**, Ch 6 — modern formulation of Wiener-Hopf factorization. Directly used: factorization 1 − sφ = (1 − sκ⁺)(1 − sκ⁻) at s = 1.
- **Kuznetsov-Kyprianou-Pardo 2010** — meromorphic Lévy class via residue analysis at poles. Considered: rational-in-2^(iθ) structure as candidate analog. Falsified: numerator carries irrational log_2 3 power, breaks rational factorization.
- **Asmussen, "Applied Probability and Queues"** Ch VIII — Esscher transform and Cramer-root analysis. Directly used: Esscher tilt at Cramer root w* = 1, tilted moment computation.
