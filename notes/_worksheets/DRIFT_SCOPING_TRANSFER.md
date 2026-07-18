# DRIFT_SCOPING_TRANSFER — Phase 3, the load-bearing question

## The question

Given: Foster–Lyapunov drift on a Syracuse Markov chain ⇒ L²-spectral gap β > 0 (Taghvaei–Mehta Theorem 1 / Hairer Harris-type result).

Does β > 0 imply a polynomial-in-A bound on |π̂_n(ξ)| for ξ ∈ (Z/3^n Z)\* uniform in n?

**Pre-registration position:** H_DRIFT_EXISTS_BUT_SPECTRAL_GAP_DOESNT_TRANSFER as most likely. This phase tests it.

## What the L²-spectral gap gives, precisely

Taghvaei–Mehta (p. 2, consequence 1):
> ‖P^n f − π(f)‖_{2,π} ≤ (1 − β)^n · ‖f − π(f)‖_{2,π}

For f ∈ L²(π), this bounds the deviation of the n-step average from the stationary expectation π(f). The bound depends on:
- (1 − β)^n — number of chain iterations from a non-stationary start
- ‖f − π(f)‖_{2,π} — L²-norm of f's deviation from its π-mean

**Critical interpretation:** π(f) is the stationary expectation; the bound is on **convergence to π(f), not on π(f) itself**.

## The Fourier-coefficient question

π̂_n(ξ) = Σ_{x ∈ Z/3^n Z} π_n(x) · exp(2πi ξ x / 3^n) = E_{π_n}[χ_ξ] = π_n(χ_ξ)

where χ_ξ(x) = exp(2πi ξ x / 3^n) is the additive character.

**For ξ ≠ 0**: by orthogonality of characters, Σ_x χ_ξ(x) = 0, hence the χ_ξ are mean-zero with respect to the **uniform** distribution on Z/3^n Z. The DEVIATION of π_n(χ_ξ) from uniform-distribution's χ_ξ-expectation (which is 0) IS exactly π̂_n(ξ).

This is the key observation. For ξ ≠ 0:
> π̂_n(ξ) = π_n(χ_ξ) − 0 = π_n(χ_ξ) − unif_n(χ_ξ)

So |π̂_n(ξ)| measures how far π_n is from uniform along the character χ_ξ.

**Spectral gap → Fourier decay candidate transfer:**
> ‖P_n^k μ_0 − π_n‖_{2,π_n} ≤ (1 − β_n)^k · ‖dμ_0/dπ_n − 1‖_{2,π_n}

If μ_0 = uniform_n (the uniform distribution on (Z/3^n Z)\*), then:
> ‖P_n^k unif_n − π_n‖_{2,π_n} ≤ (1 − β_n)^k · ‖unif_n / π_n − 1‖_{2,π_n}

Now P_n^k unif_n → π_n as k → ∞. The Fourier-coefficient question reverses time: as k → ∞, we want |⟨χ_ξ, π_n⟩ − ⟨χ_ξ, unif_n⟩| = |π̂_n(ξ)|.

In the **limit k → ∞**, the spectral gap controls **how fast P_n^k unif_n approaches π_n**, but at k → ∞ exactly we have the limit π_n itself; the bound (1 − β_n)^k → 0 gives no information about π_n's deviation from unif_n.

**Re-formulated correctly:** the spectral gap bound applied at any finite k gives convergence of P_n^k μ_0 → π_n; choosing μ_0 = unif_n means at k = 0 we have the deviation ‖unif_n − π_n‖, which IS what we want to bound, BUT the inequality bounds the RATE of decay starting from this initial deviation — it does NOT bound the initial deviation itself.

> **The spectral gap β bounds the decay of ‖P_n^k μ_0 − π_n‖ in k, NOT ‖μ_0 − π_n‖ at k = 0.**

## Where the transfer goes wrong (the structural mismatch)

The drift framework asks: **how quickly does the chain mix to π?** Answer: rate β.

The Fourier-decay question asks: **how arithmetically uniform is π itself?** Answer: needs |π̂_n(ξ)| ≤ poly(1/n).

These are **different objects**:
1. Mixing rate β is about the **transient regime** (P^k → π).
2. |π̂_n(ξ)| is about the **stationary regime** (π itself).

There is **no general theorem** in the three references (Glynn-Zeevi, Hairer, Taghvaei–Mehta) connecting these. The literature reviewed in Taghvaei–Mehta's bibliography ([10] Kontoyiannis–Meyn, [12] Meyn–Tweedie, [13] Roberts–Rosenthal, [14] Rosenthal) all give convergence rates — not Fourier decay of π.

## Three negative arguments

### (T1) Order-of-quantifiers mismatch

L²-spectral-gap consequence: ∀ k ≥ 1, ‖P^k f − π(f)‖_2 ≤ (1 − β)^k · ‖f − π(f)‖_2.

Fourier-decay target: ∀ ξ ≠ 0, |π̂_n(ξ)| ≤ poly(1/n, 1/|ξ|, ...).

The spectral gap quantifies over (k, f); the Fourier target quantifies over ξ at fixed (k = ∞, f = π). The spectral-gap quantification gives no information at k = ∞ because (1 − β)^∞ = 0 trivially.

### (T2) Non-negativity scope of Glynn-Zeevi

Glynn-Zeevi's bound πf ≤ c requires f ≥ 0. The Fourier character χ_ξ = exp(2πi ξ ·/3^n) is **complex-valued and not non-negative**. Decomposition χ_ξ = cos + i sin into real/imaginary parts gives **two real-valued oscillatory functions, not non-negative ones**. The framework's natural scope does not accommodate them.

Could one write |π̂_n(ξ)|² = E_π[χ_ξ] · conj(E_π[χ_ξ]) and use bilinear drift? Not within Glynn-Zeevi's framework — bilinear drift inequalities are a strictly different machinery (Lyapunov for bilinear forms, e.g., variance Lyapunov).

### (T3) Constants are n-dependent for Chain B

Even if one accepts T1 + T2 are surmountable, on Chain B at level n with finite state space:
- β_n = λ/(1 + 2b/α) with α scaling at least as 3^{−n} (Doeblin on a 2·3^{n−1}-state finite chain has Doeblin's α at most |minorizing set|/|state space| ≤ |K|/2·3^{n−1}).
- Resulting β_n ≤ poly(n) / 3^n at best.

So **even if** the transfer T1 + T2 were resolved, the spectral gap **itself** is exponentially small in n for the natural choices of V, K, ν. **Uniform-in-n spectral gap is not provided by the framework** for Chain B.

Achieving a **uniform-in-n** spectral gap on Z/3^n Z would require a smooth-IFS or continuous-group structure that the drift framework explicitly doesn't use — which puts us back in the territory of Probes 1-3 (BKS, Furstenberg, ARHW), all of which failed.

## Can the literature beyond the three references rescue the transfer?

Searched mental literature for any known theorem connecting Markov-chain spectral gap to stationary-distribution Fourier decay:
- **Bobkov–Tetali (modified log Sobolev)**: gives entropy/KL decay, not Fourier decay.
- **Diaconis–Saloff-Coste (Cheeger / spectral gap on finite groups)**: gives ℓ²(unif) mixing time bound, expressed in terms of |1 − λ_2(P)| ≥ β. For random walks on finite **groups**, mixing time IS related to character decay via the Plancherel formula: ‖P^k unif − π‖_{2,unif}² = Σ_{χ ≠ trivial} |χ̂(P)|^{2k} · |something|. **This is for random walks on groups (e.g., Cayley graphs).** Syracuse Chain B on Z/3^n Z is NOT a Cayley graph walk — the kernel depends on the state (multiplication by 2^{−v(3x+1)} where v = v_2(3x+1) depends on x). So Diaconis–Saloff-Coste's character-eigenfunction structure does not apply.
- **Bourgain–Furman–Lindenstrauss–Mozes (random walks on torus)**: requires the random walk to be by a fixed measure on SL_d(Z); gives quantitative equidistribution. Syracuse chain doesn't have this structure (the multipliers 2^{−v} are scalars not in SL).
- **He–de Saxcé (linear random walks on torus)**: extends BFLM; same structural mismatch.

**No literature connection found** between Foster-Lyapunov drift on Syracuse Chain B and Fourier decay of π_n.

## Conclusion of Phase 3

The transfer from L²-spectral gap to Fourier decay of the stationary measure π_n is **structurally broken**:

1. (T1) Spectral gap is a statement about **transient mixing**, Fourier decay is a statement about **stationary arithmetic**. Different objects.
2. (T2) Glynn-Zeevi's framework is for non-negative test functions; characters χ_ξ are not non-negative.
3. (T3) Even granting T1 + T2, the spectral gap β_n on Chain B is at best exponentially small in n (no uniform spectral gap from Foster-Lyapunov data alone).

**H_DRIFT_EXISTS_BUT_SPECTRAL_GAP_DOESNT_TRANSFER is the disposition.**

The pre-registration was correct: the framework gives a wrong-flavored bound (TV/L²-mixing of P^k to π, not Fourier of π itself), and the constants are n-dependent.
