# DRIFT_SCOPING_HYPOTHESES — Phase 1, verbatim framework hypotheses

**Date:** 2026-05-12. Scoping probe of Foster–Lyapunov / Glynn–Zeevi / Lyapunov-Foster-Poincaré frameworks for Syracuse μ_n on Z/3^n Z. First Markov-chain-native scoping after three closed continuous/smooth probes.

## 1. Taghvaei–Mehta (arxiv:2005.08145) — the explicit "drift → spectral gap" theorem

**Setting.** Discrete-time time-homogeneous Markov process {X_n}_n on Polish state space (X, B), Markov operator Pf(x) = E[f(X_1) | X_0 = x], invariant measure π. (p. 1)

**Assumption 1 (reversibility, load-bearing).** "P admits a unique reversible invariant measure π" — i.e., P self-adjoint on L²(π):
> ⟨f, P g⟩_π = ⟨P f, g⟩_π, ∀ f, g ∈ L²(π).

**Assumption 2 (Foster–Lyapunov "(v4)" condition).** A positive function V: ℝ^d → [1, ∞), constants b < ∞, α, λ > 0, set K ⊆ X, probability measure ν such that:
> **(4) drift:**  P V ≤ (1 − λ) V + b · 1_K
> **(5) minorization:**  P 1_A(x) ≥ α ν(A) 1_K(x), ∀ A ∈ B.

**Theorem 1 (Poincaré inequality).** Under Assumptions 1–2,
> "P admits a Poincaré inequality (2) with constant β_+ = λ / (1 + 2b/α)."

**Corollary 1 (spectral gap, positive-semidefinite case).** Under same assumptions, if P is positive semi-definite,
> "P admits a spectral gap β = β_+ = λ / (1 + 2b/α)."

**Counter-example 1 (p. 3, §III.B).** Without additional structure, Theorem 1 yields a bound on β_+ only — NOT on β_−. The explicit counter-example (P = [[ε, 1−ε],[1−ε, ε]]) shows constants in (4)–(5) can be ε-independent while the true gap β_− = 2ε goes to 0. So drift + minorization alone do not give a spectral-gap bound in general.

**Extensions.**
- §IV.A (Prop 1, "stronger condition" Assumption 3): if K = {V ≤ R}, R > 2b/λ, the spectral gap for P² gives ‖P‖ ≤ (1 − β_+)^{1/2}.
- §IV.B (Prop 2, non-reversible): if BOTH P and its L²(π) adjoint P† satisfy (4)–(5) with the same V, K, λ, b, α, then P†P satisfies a drift inequality. Spectral gap conclusion is in terms of ‖P‖_{L²₀(π)}.

**Stated conclusions of spectral gap β > 0:**
- §II.A consequence 1 (p. 2): "Geometric convergence of the moments in L²(π)"
  > ‖P^n f − π(f)‖_{2,π} ≤ (1 − β)^n · ‖f − π(f)‖_{2,π}.
- §II.A consequence 2: "Geometric convergence of the probability distribution in the total-variation distance" ‖μP^n − π‖_{TV} ≤ (1 − β)^n · ‖h − 1‖_{2,π}.

**What the framework bounds (verbatim, scope-limited):** convergence rate of P^n to π in L²(π) and in TV. **It does NOT bound Fourier coefficients of π itself.** π̂(ξ) = ⟨exp(2πi ξ·), 1⟩_π is the L²-inner product of the trivial function with the character — but the framework gives convergence of ⟨P^n f, g⟩ → ⟨π(f), g⟩, which is a different object.

## 2. Glynn–Zeevi (2008) — bounds on stationary expectations

**Setting (p. 195).** Markov jump process X = (X(t))_{t≥0} on discrete S, rate matrix Q, stationary distribution π. (Also: DTMC, SDE, jump diffusion — Theorem 1 phrased in unified extended-generator form.)

**Theorem 1 (p. 199).** Suppose g ∈ D(A) is a **non-negative** function with sup_x (Ag)(x) < ∞. Then:
> (i) −E_x ∫_0^t (Ag)(X(s)) ds ≤ g(x);
> (ii) ∫_S π(dx) (Ag)(x) ≥ 0.

**Corollary 1 (jump process), 2 (SDE), 3 (jump diffusion), 4 (DTMC) (pp. 200):** all phrased for **non-negative f: S → R⁺**.

E.g., **Corollary 4 (DTMC):** "Let X = (X_n)_{n≥0} be a discrete-time S-valued Markov chain with kernel P, and suppose f: S → R is **non-negative**. If there exists a **non-negative** function g: S → R and constant c with ∫ P(x, dy) g(y) ≤ g(x) − f(x) + c, then ∫_S π(dx) f(x) ≤ c."

**Lower bound theorems (§3, Prop 2/3 and Cor 5/6):** require π-integrability of g, give lower bounds on πf for **non-negative f**.

**What the framework bounds (verbatim, scope-limited):** stationary expectations π(f) of **real-valued non-negative cost functions f** via the drift inequality A g ≤ −f + c. Returns πf ≤ c (upper bound) or πf ≥ c̃ (lower bound). Throughout, f is non-negative (the "cost function" / "instantaneous cost rate"; cf. p. 195 "cost accrues when X is in state x").

**Critical scope observation:** The paper's bounds f ≥ 0 / g ≥ 0 are not stylistic — the non-negativity is load-bearing in the proof (Theorem 1.i uses Fatou; Cor 1 uses Q g ≤ −f + ce and integrates against π exploiting (Ag)(x) ≥ 0 a.e.). The framework **does not naturally accommodate complex-valued oscillatory test functions** like χ(x) = exp(2πi q x / 3^n).

## 3. Hairer — Convergence of Markov Processes (lecture notes)

**Setting (p. 1).** Discrete- and continuous-time Markov processes on Polish X with transition kernel P, invariant measure µ.

**§1.2 Foster–Lyapunov criteria (Proposition 1.3, p. 4) — countable state space.** Discrete generator L = P − I.
> - **Transient** iff ∃ V: X → R⁺, A ⊆ X non-empty, L V(x) ≤ 0 for x ∉ A, ∃ x ∉ A with V(x) < inf_A V.
> - **Recurrent** iff ∃ V: X → R⁺ with {V ≤ N} finite ∀ N, L V(x) ≤ 0 for all but finitely many x.
> - **Positive recurrent** iff ∃ V: X → R⁺ with L V(x) ≤ −1 for all but finitely many x.

**§3 Harris-type theorems (V-norm spectral gap).**
> "Level sets [Has80, MT93]. If the Lyapunov function is strong enough, one then has a spectral gap **in a weighted supremum norm** [MT92, MT93]."

The weighted norm is ‖φ‖ = sup_x |φ(x)| / V(x) (p. ~16). Harris-type result is: ‖P^n φ − π(φ)‖ ≤ C γ^n ‖φ‖ for some γ < 1 — **convergence in V-weighted supremum norm**, NOT Fourier-coefficient decay of π.

**§4 Total-variation convergence (p. ~22).** Weighted-TV convergence of P^n µ_0 to π, expressed as ∫|D_µ(x) − D_ν(x)| G(x) dx for weight G.

**What the framework bounds (verbatim, scope-limited):**
- Recurrence/transience dichotomy (Prop 1.3)
- Convergence of P^n to π in V-weighted supremum norm (§3)
- Convergence in weighted total-variation norm (§4)

**It does NOT bound Fourier coefficients of π itself.** Hairer's notes are about **chain mixing time to π**, not about **arithmetic properties (Fourier decay) of π**.

## 4. The shared structural conclusion across all three references

All three frameworks bound the **rate at which the chain reaches stationarity** (in L²(π), TV, or V-weighted sup norm), parameterized by the drift constants (λ, b, α, K). None of the three bounds **Fourier coefficients π̂(ξ) of the stationary distribution itself**.

The corpus INDEX framing — "polynomial-in-A bounds on stationary expectations" — is technically accurate but rhetorically misleading for the c = 7/45 unblocking question:
- π̂_n(ξ) = E_π[exp(2πi ξ · / 3^n)] IS a stationary expectation.
- BUT the test function exp(2πi ξ x / 3^n) is **complex-valued and not non-negative**, falling outside Glynn–Zeevi's natural scope (f ≥ 0 throughout).
- Glynn–Zeevi Cor 4 with f = exp(2πi ξ x / 3^n) would require BOTH f ≥ 0 and existence of a non-negative drift function g with P g ≤ g − f + c. The first is violated by definition.

The framework's natural output is bounds on E_π[V(X)] for the Lyapunov function V (Theorem 1.i with f = some non-negative function dominated by Ag): "π V is finite, with bound c." This bounds **moments of π under V**, not characters of π.

## 5. The transfer step (the load-bearing question)

For the framework to deliver the polynomial-in-A Fourier bound, one of two transfer mechanisms must exist:

**(T-A)** A literature theorem connecting L²-spectral gap β > 0 of P to Fourier-coefficient decay |π̂(ξ)| ≤ C(β, ...) · ψ(ξ) on Z/3^n Z. **Not found in any of the three cited references.** Roberts–Rosenthal [13], Cattiaux–Guillin [4], Kontoyiannis–Meyn [10] in 2005.08145's bibliography all give convergence-rate results — not Fourier decay of π.

**(T-B)** A character-by-character argument: for each character χ_ξ, define a chain-related quantity f_ξ such that π(f_ξ) factors through π̂(ξ), find a Lyapunov function g_ξ controlling drift in f_ξ, derive a Glynn–Zeevi-style bound. **No general mechanism** for this in the three references — characters χ_ξ are complex-valued and oscillatory, breaking the non-negativity hypothesis.

The transfer step is **not addressed by any of the three references**. The frameworks were built for non-arithmetic questions (mixing time, queueing, MCMC convergence, stochastic stability) where Fourier coefficients of π are not the object of interest.
