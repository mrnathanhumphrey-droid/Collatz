# DRIFT_SCOPING_DISPOSITION — top-level

**Date:** 2026-05-12. Probe 4 (Foster–Lyapunov / Glynn–Zeevi / Lyapunov-Foster-Poincaré) of the c=7/45 framework-compatibility chain. First **Markov-chain-native** scoping after three closed continuous/smooth-dynamical probes (L²-flattening / SL_2(ℝ) / cocycle-Dolgopyat — all H_*_FAILS).

Reporting to Nathan.

---

## DISPOSITION: **H_DRIFT_EXISTS_BUT_SPECTRAL_GAP_DOESNT_TRANSFER**

The Foster-Lyapunov drift framework applies cleanly to Chain B (Syracuse Markov chain on (Z/3^n Z)\ {0}) — drift functions exist trivially on the finite state space, minorization holds via Doeblin's condition, an L²-spectral gap β_n > 0 follows. **But spectral gap → Fourier-coefficient decay of π_n is structurally broken**, for three independent reasons:

1. **Object mismatch (T1).** Spectral gap bounds the rate at which P_n^k μ_0 → π_n (transient mixing). Fourier decay of π_n is a statement about the **stationary** distribution itself, not its rate of approach. The quantification (∀ k ≥ 1) gives no information at the limit k → ∞.

2. **Scope mismatch (T2).** Glynn-Zeevi's framework is for **non-negative** test functions (verbatim: f: S → R⁺ in Theorem 1, Cor 1–4). Fourier characters χ_ξ = exp(2πi ξ ·/3^n) are complex-valued oscillatory functions, not non-negative. The framework's hypothesis excludes them.

3. **Constants are n-dependent (T3).** Even granting T1 + T2 surmountable, on Chain B the Foster-Lyapunov constants (λ, b, α, ν(K)) for any natural V scale with n; spectral gap β_n ≲ poly(n)/3^n at best. **The framework does not deliver a uniform-in-n spectral gap** without additional smooth/IFS structure (which is exactly what Probes 1-3 ruled out for Syracuse).

The corpus INDEX framing ("polynomial-in-A bounds on stationary expectations") is **technically correct but rhetorically misleading**: the framework bounds π(f) for non-negative cost functions f via drift, not Fourier coefficients of π. The pre-registered most-likely hypothesis is confirmed.

---

## Pre-registration outcome

| Hypothesis | Status |
|---|---|
| H_DRIFT_FRAMEWORK_GIVES_FOURIER_DECAY | **NO** — no literature theorem found connecting L²-spectral gap of Syracuse Chain B to Fourier decay of π_n; the three references are about chain mixing, not stationary-distribution arithmetic. |
| **H_DRIFT_EXISTS_BUT_SPECTRAL_GAP_DOESNT_TRANSFER** | **CHOSEN** — drift functions trivially exist on finite Chain B; spectral gap follows; transfer to Fourier decay broken by T1+T2+T3. |
| H_DRIFT_FUNCTION_DOESNT_EXIST | NO — bounded V on finite state space trivially satisfies drift. The drift inequality is not the obstruction. |
| H_FRAMEWORK_GIVES_WRONG_RATE | partial overlap — the framework gives n-dependent (exponentially small) constants, but the primary issue is the wrong-flavored object (mixing vs Fourier), not the rate. |
| H_AMBIGUOUS | NO — three reinforcing structural arguments converge on the same disposition. |

Pre-registration favored H_DRIFT_EXISTS_BUT_SPECTRAL_GAP_DOESNT_TRANSFER. **Confirmed.**

---

## Adversarial checks resolved

**(A1) Pre-registration honesty.** H_DRIFT_FRAMEWORK_GIVES_FOURIER_DECAY was the optimistic outcome. Searched bibliographic context (Kontoyiannis-Meyn, Meyn-Tweedie, Roberts-Rosenthal, Rosenthal, Bobkov-Tetali, Diaconis-Saloff-Coste, BFLM, He-de Saxcé) for a literature theorem connecting L²-spectral gap on Syracuse-type chains to Fourier decay of stationary measure. **Not found.** The Diaconis-Saloff-Coste Plancherel-character-on-Cayley-graph mechanism would be the closest if Chain B were a Cayley-graph walk, but Chain B's kernel is state-dependent (multiplier 2^{−v(3x+1)} depends on x) — not a Cayley graph walk.

**(A2) Reversibility.** Chain B is not reversible. Taghvaei-Mehta Theorem 1 requires reversibility (Assumption 1); must use Proposition 2 (non-reversible) which is strictly stronger. This is a stacked technical obstruction on top of the structural one.

**(A3) State space.** The natural Markov chain on Z (Chain A) has no proper stationary measure (Collatz orbits absorb to {1,2,4}). The chain on Z/3^n Z (Chain B) has the right stationary measure but is a different chain at each n, with no uniform-in-n control from the framework. The limit chain on Z_3 has μ_∞, but its "Fourier coefficients" are precisely the level-n characters of Z/3^n Z — same object, same transfer problem.

**(A4) Transfer is load-bearing.** Did not paper over T1 (mixing-rate quantification gives no information at k = ∞), T2 (non-negativity exclusion of complex characters), T3 (n-dependent constants). Each is a separate, sufficient obstruction.

**(A5) §5 inherited-claim discipline.** The INDEX claim "drift conditions can give polynomial-in-A bounds on stationary expectations" is verified for **non-negative cost functions** f (Glynn-Zeevi Cor 4 verbatim) — accurate as stated. The **slippage** is in interpreting "stationary expectations" to include Fourier coefficients π̂(ξ) = π(χ_ξ); the characters χ_ξ are not non-negative, so the bound's hypothesis is violated. The INDEX framing is technically correct, the rhetorical application to Fourier decay is not.

---

## Trajectory placement

| Probe | Framework | Disposition |
|---|---|---|
| 1 | Baker-Khalil-Sahlsten L²-flattening 2407.16699 | **H_L2_FLATTENING_FAILS** (continuous/smooth) |
| 2 | Furstenberg/Hochman-Solomyak SL_2(ℝ) 2108.06006 / 1610.02641 | **H_SL2_EMBEDDING_DOESNT_EXIST** (continuous/smooth) |
| 3 | Algom-Baker-Sahlsten cocycle Dolgopyat 2306.01275 | **H_COCYCLE_DOLGOPYAT_LINEAR_EXCLUSION** (continuous/smooth) |
| **4 (this)** | **Foster-Lyapunov / Glynn-Zeevi / 2005.08145** | **H_DRIFT_EXISTS_BUT_SPECTRAL_GAP_DOESNT_TRANSFER** (Markov-chain-native) |

**Four structural negatives.** Three from continuous/smooth-dynamical frameworks (category-of-object mismatch: discrete arithmetic chain in a continuous-smooth setting). One from the Markov-chain-native framework (wrong-flavored conclusion: chain-mixing bounds, not Fourier-of-π bounds).

The unifying meta-pattern: **modern Fourier-decay frameworks deliver their conclusion for objects that Syracuse μ_n on Z/3^n Z is not** — smooth-IFS measures, Furstenberg measures on P^1, Markov-chain time-to-stationarity. The polynomial-in-A bound on |π̂_n(ξ)| is a **discrete-arithmetic Fourier-decay statement** about the **stationary distribution itself**, not its rate of approach, and existing frameworks are not built around this exact object.

---

## Routing — what becomes the natural next move

The corpus's three Fourier-decay framework probes (1, 2, 3) and the one Markov-chain-native probe (4) have all returned structural negatives. The routing options narrow:

### Primary recommendation: Tauberian arc (Flajolet-Sedgewick / Chevalier 2507.15394)

Per the MEMORY entry `project_collatz_r78_bilinear_cracked.md`: "Tauberian arc opened: Flajolet-Sedgewick Ch. VI / Chevalier 2507.15394 Thm 1.16 candidate, single-theorem selection pending ε_7 exact-rational compute. 73 PDFs literature bundle at burgess/literature/."

The Tauberian framing is **structurally different** from the four probes:
- Probes 1-3 ask: does a Fourier-decay framework apply to π_n? (NO across the board.)
- Probe 4 asks: does a Markov-chain mixing framework give Fourier decay of π_n? (NO — wrong object.)
- Tauberian asks: does the asymptotic structure of Σ_k |π̂_k(ξ)|² (or related generating series) give a tail bound on |π̂_n(ξ)|?

This is a **generating-function / complex-analytic** route, NOT a dynamical or chain-mixing route. The four structural failures of frameworks built around the dynamical operator T_lead suggest the next move should be **away from T_lead direct compatibility** and toward a **generating-series / Plancherel-trace** mechanism. This is exactly the Tauberian arc.

The Chevalier 2507.15394 Thm 1.16 candidate, with ε_7 exact-rational compute pending, is the immediately actionable next probe.

### Secondary recommendation: Bourgain-Konyagin discrete sum-product on Z/3^n Z

Per Cocycle-Dolgopyat disposition routing #3: "Bourgain-Konyagin on sum-product in Z/p^k Z." This is the most natural framework for discrete-arithmetic Fourier decay (it's literally designed for the question "Fourier decay of measures on Z/p^k Z"). Per MEMORY, the literature bundle is at `C:/Collatz/Bourgain-Konyagin`.

This is a HARDER probe (Bourgain-Konyagin in characteristic p is technical and the Syracuse-specific sum-product structure is non-trivial) but it's the categorically correct framework: discrete-arithmetic question, discrete-arithmetic answer. No category mismatch.

### Tertiary: argument-portability scoping

If both Tauberian and Bourgain-Konyagin close, the remaining option is to step back to Tao's Prop 1.17 directly and ask whether the iterated-cubic obstruction can be circumvented **at the level of the recursion** (not by importing an external Fourier-decay framework). This is the most ambitious of the three.

---

## Synopsis (one paragraph)

The Foster–Lyapunov drift / Glynn–Zeevi stationary-expectation / Lyapunov-Foster-Poincaré spectral-gap framework family — the corpus's "Markov-chain-native" candidate, in contrast to the three continuous/smooth-dynamical probes already closed — provides bounds on the rate at which the Syracuse Markov chain on Z/3^n Z approaches its stationary distribution π_n (in L²(π_n), TV, V-weighted sup norm), parameterized by drift constants. Drift functions exist trivially on the finite state space; minorization holds via Doeblin; spectral gap β_n > 0 follows. But the framework does NOT bound the Fourier coefficients π̂_n(ξ) of π_n itself: spectral gap quantifies transient mixing not stationary arithmetic (object mismatch), Glynn–Zeevi's hypothesis is non-negative test functions excluding complex characters (scope mismatch), and the constants on Chain B scale with n giving at best exponentially small β_n (no uniform-in-n control). The corpus INDEX framing "polynomial-in-A bounds on stationary expectations" is technically correct for non-negative cost functions f but does not deliver Fourier-coefficient decay of π_n. **The natural next move is the Tauberian arc (Flajolet-Sedgewick Ch. VI / Chevalier 2507.15394 Thm 1.16), which is structurally different from the four exhausted dynamical/chain-mixing probes — it operates on generating series / Plancherel traces, not on the dynamical operator T_lead directly. Secondary route: Bourgain-Konyagin discrete sum-product on Z/3^n Z, which is categorically correct (discrete-arithmetic question, discrete-arithmetic framework) though technically harder.** Four framework probes closed; the next move is OUT of the framework-compatibility scoping pattern and INTO the generating-series / sum-product arena.

---

## Deliverables

- `DRIFT_SCOPING_HYPOTHESES.md` — Phase 1, verbatim Taghvaei-Mehta / Glynn-Zeevi / Hairer hypotheses, with §4-§5 noting the framework's natural scope (TV/L²-mixing, non-negative test functions) does not include Fourier-of-π
- `DRIFT_SCOPING_CANDIDATES.md` — Phase 2, four candidate V's on Chains A/B/Z_3, each fails uniformly: A has no stationary measure for c=7/45, B has trivial bounded V's but n-dependent constants, Z_3 has stationary measure but the framework's conclusion is wrong-flavored
- `DRIFT_SCOPING_TRANSFER.md` — Phase 3 (load-bearing): three independent obstructions (T1 object mismatch, T2 non-negativity scope, T3 n-dependent constants) to spectral-gap → Fourier-decay transfer; no literature theorem connecting them found
- `DRIFT_SCOPING_DISPOSITION.md` — this file
