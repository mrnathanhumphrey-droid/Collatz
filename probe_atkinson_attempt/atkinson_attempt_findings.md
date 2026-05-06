# Atkinson Rota-Baxter construction attempt — findings

**Date:** 2026-05-06.
**Time:** ~3 hours of the 9-hour budget.
**Outcome:** mixture of C and D (see disposition below).

This is the load-bearing deliverable for the probe. The construction did not succeed in a non-trivial sense at k=3, so Phase 3 (k=5) was not run per the pre-registered gate (brief: "if construction succeeds at k=3, repeat at k=5"; trivial-projector successes don't count as "construction succeeds" in the sense the brief asked about).

## 4.1 Summary of constructions attempted

Phase 2 tested 7 candidate operators T : C^18 → C^18 for the Rota-Baxter axiom of weight θ ∈ {0, 1}, where C^18 = functions on the 18 coprime classes of (Z/27)* under **pointwise multiplication**.

| Candidate | Description | RB axiom (θ=1) max residual | Notes |
|---|---|---|---|
| T_a | Diagonal projection onto first 9 indices | 0 | Trivial pass (projector, subset-indicator) |
| T_c | Diagonal projection onto residues r < 14 | 0 | Trivial pass |
| T_d | Diagonal projection onto r ≡ 1 (mod 3) | 0 | Trivial pass |
| T_K | K_3 Markov kernel (row-stochastic) | 2.44 | FAIL — not a projector, not RB |
| T_resK | (I − K_3)⁺ pseudoinverse | 5.27 | FAIL |
| T_K2 | K_3² | 1.26 | FAIL |
| T_chi | Fourier projector onto even-index characters of (Z/27)* | 1.62 | FAIL — Fourier subspaces are NOT pointwise-multiplicatively closed |

(Full per-pair residuals and spectra in `construction_attempts_k3.csv`.)

## 4.2 What "passes the axiom" actually means here

On a finite commutative semisimple algebra C^N (with pointwise product), the only subalgebras are the *coordinate subalgebras* `R_A = {f : f|_{A^c} = 0}` for subsets A ⊆ {1, ..., N}. Each such subset gives a diagonal projector T_A that automatically satisfies the Rota-Baxter axiom with θ = 1 (because both T_A·R = R_A and (Id − T_A)·R = R_{A^c} are subalgebras, and a projector with both image and kernel as subalgebras automatically satisfies eq. (5.2) per Atkinson).

So **every subset of {1,...,18} gives a valid weight-1 Rota-Baxter projector**. There are 2^18 = 262,144 such projectors. The axiom is informationally void on this algebra: it doesn't pick out a canonical operator.

The empirical "passes" of T_a, T_c, T_d are three arbitrary instances of this trivially-large family. None of them connects to the framework's empirical structure in any non-trivial way.

## 4.3 Spectrum analysis

For every passing T (all are projectors), the spectrum is exactly **{0 (mult 9), 1 (mult 9)}**. No geometric structure, no band, no 1−q^n sequence.

| Reference | Empirical / claimed spectrum |
|---|---|
| Atkinson Section IX 1−q^n (different RB regime) | infinite geometric sequence |
| Framework R_k singular value band | continuous band on [0.488, 0.671], k-stable |
| Framework K_k spectrum | {1} ∪ near-zero cluster (~10^-3 magnitude) |
| Atkinson Section VII projector | {0, 1} |
| **This construction's T_a, T_c, T_d** | **{0, 1}** |

The construction's spectrum matches Atkinson's projector spectrum exactly, which was already evident before computation. It does not match anything in the framework's empirical operator spectra.

**Important qualification on the 1−q^n claim** (from reading Atkinson Section IX more carefully than the brief did):

The Section IX statement "eigenvalues will form sequences of the form 1−q, 1−q^2, ..." cites equations (6.3)/(6.4), which appear in **Section VI**, not Section VII. Section VI's example uses T = (E−V)⁻¹ where V is a homomorphism (Vf(ξ) = f(qξ)). T is **not a projector** and is NOT the discrete Wiener-Hopf operator. The 1−q^n eigenvalue sequence is a feature of the homomorphism-resolvent class of RB operators, not the Wiener-Hopf projector class.

The brief conflated these. The framework's empirical R_k band σ ∈ [0.488, 0.671] cannot match Atkinson Section VII (whose projectors give only {0, 1}); and it doesn't fit 1−q^n geometric sequence either (the band is dense, not a discrete geometric sequence; band-width 0.18 is not consistent with successive 1−q^n where q < 1 forces the n=1 term to be the smallest).

## 4.4 Repunit position — a category error in the brief

The "repunit residue (4^k − 1)/3" target set {1, 5, 21, 85, 341, ...} is a **2-adic** feature of the framework — these are the residues mod 2^k that take the maximum number of deterministic-prefix steps. They are residues in (Z/2^k)*, not in (Z/3^k)*.

At k=3, the brief's expected repunit is 21. But:
- 21 mod 3^3 = 21
- 21 = 3 · 7, so **21 is NOT coprime to 3**
- Therefore 21 is **NOT a state** in the K_k Markov chain on (Z/3^k)*
- The character algebra of (Z/27)* doesn't contain an indicator for residue 21

So the construction's natural algebra (3-adic character algebra) cannot index the repunit's special role. The repunit lives in a different algebra (2-adic) that this Atkinson construction doesn't reach. This is a structural mismatch, not a failure of the construction technique.

## 4.5 Failure modes for non-projector candidates

T_K, T_resK, T_K2, T_chi all fail the RB axiom with residuals of order 1-5. The reasons:

- **K_3 (Markov kernel):** stochastic, not multiplicative on functions. K(f·g) ≠ K(f)·K(g) generically. Markov kernels are linear operators on the function algebra but not algebra homomorphisms; they don't fit the homomorphism-resolvent class either (since K is not invertible — eigenvalue 1 — and T_resK = (I−K)⁺ doesn't preserve the algebra structure).

- **T_chi (Fourier projector):** projects onto a subspace closed under *convolution*, not under pointwise multiplication. C^N as a pointwise algebra and ℓ²((Z/27)*) as a convolution algebra are dual via Fourier, but the subalgebra structures are different. T_chi's image is a convolution-subalgebra, not a pointwise-subalgebra, so the RB-with-θ=1 axiom fails on the pointwise algebra.

These failures aren't accidents — they reflect real algebraic mismatches. To get a non-trivial RB structure, one would need either:
1. A different algebra structure (convolution instead of pointwise), where Fourier projectors become subalgebras
2. A non-projector RB operator with weight 1, of the form T = (E−V)⁻¹ for a multiplicative homomorphism V on the chosen algebra
3. A weight-θ ≠ 1 setup where K_k or its variants fit naturally

## 4.6 Honest scope

This probe is a Phase 2 construction attempt at k=3 only, with Phase 3 (k=5 generalization) skipped because Phase 2's "successes" are informationally trivial. The conclusions:

- The framework's natural operators (K_k, R_k) do **not** appear to be Rota-Baxter operators on the natural character algebra of (Z/3^k)* under pointwise multiplication — at least not in any way detectable by the Phase 2 candidate sweep.
- The Atkinson Section VII discrete Wiener-Hopf structure applies trivially (every subset projector works) but produces only {0, 1} spectrum, which doesn't match the framework's empirical operator spectra.
- The Atkinson Section IX "1−q^n" eigenvalue sequence comes from a different RB regime (Section VI's homomorphism-resolvent), not from the discrete Wiener-Hopf projector setup. The brief's hope of matching the empirical band [0.488, 0.671] to 1−q^n was based on a misreading of which Atkinson construction produces that spectrum.
- The repunit target set {1, 5, 21, 85, ...} is a 2-adic feature, not indexable in the 3-adic character algebra of (Z/3^k)*.

## 4.7 Disposition vs pre-registered outcomes

- **Outcome A (construction succeeds at k=3 and k=5):** NOT obtained.
- **Outcome B (works at k=3, fails at k=5):** N/A (Phase 3 not run).
- **Outcome C (no viable T):** PARTIAL — diagonal subset projectors trivially satisfy axiom, but they're informationally void. Non-trivial candidates (K_3, K_3², (I−K)⁺, Fourier projector) all fail.
- **Outcome D (reading reveals incompatibility):** PARTIAL — the brief's expected eigenvalue match (1−q^n band) confuses two distinct RB regimes in Atkinson. Section VII (Wiener-Hopf projector) gives {0, 1}, not 1−q^n. Section VI (homomorphism resolvent) gives 1−q^n but doesn't apply to row-stochastic Markov kernels. The framework's empirical band [0.488, 0.671] doesn't fit either Atkinson regime.

**Net disposition:** Outcome C with strong Outcome-D flavor. The Atkinson Rota-Baxter algebraic frame, as the brief proposed it, doesn't directly apply to the framework's K_k / R_k operators on the natural character algebra of (Z/3^k)*.

## 4.8 Recommended next directions

If the goal of finding an algebraic frame for the empirical band-supported R_k spectrum and the repunit's special structural role remains live, the leads are:

1. **Try the convolution algebra rather than the pointwise algebra.** The character group (Z/27)* ≅ ℤ/18 is cyclic. Functions on (Z/27)* under convolution form a ℤ/18-graded commutative algebra (Pontryagin dual structure). Fourier projectors onto subsets of ℤ/18 are convolution-subalgebra projectors → automatically RB. This brings the Atkinson Section VII frame back as a non-trivial structure with a natural ℤ/18 "winding number" index. Whether the framework's R_k operators are RB on this algebra is a separate test.

2. **Rota-Baxter weight θ ≠ 1 on the convolution algebra.** Could surface a homomorphism-resolvent class with 1−q^n eigenvalues for some q. Match attempt to band [0.488, 0.671] — under the corrected reading, it's still unlikely to match a discrete geometric sequence, but worth one shot.

3. **2-adic side, not 3-adic.** The repunit target set lives in (Z/2^k)*, not (Z/3^k)*. A separate RB construction on the 2-adic character algebra might capture the repunit's role. The framework's deterministic-prefix algorithm naturally lives on the 2-adic side, which is the correct algebra for "repunit indexing."

4. **Pollicott-Ruelle resonance theory** (suggested as Outcome-D fallback in the brief). Liverani-Faure-Sjöstrand's anisotropic Banach space approach for transfer operators with band spectra is a more natural algebraic frame for K_k's empirical structure. The R_k band on [0.488, 0.671] is suggestive of an essential spectrum band in that theory rather than a discrete RB eigenvalue sequence.

5. **Non-commutative RB (EFP Section 4.7).** The framework's K_k acts on a function algebra; a non-commutative RB formulation might handle Markov-kernel-type operators that the commutative Rota-Baxter setup doesn't accommodate.

## Files

- [`reading_summary.md`](reading_summary.md) — Phase 1 reading notes, identifying the Section VII vs Section IX RB regime distinction
- [`construction_attempts_k3.csv`](construction_attempts_k3.csv) — per-candidate axiom verification at k=3 (7 candidates × 2 θ values)
- [`spectrum_comparison.csv`](spectrum_comparison.csv) — full spectra of axiom-passing candidates (all are {0, 1} projectors)
- [`probe_phase2_k3.py`](probe_phase2_k3.py) — script
- [`atkinson_attempt_findings.md`](atkinson_attempt_findings.md) — this writeup

Time spent: ~3 hours of 9-hour budget. Phase 3 (k=5) skipped per gate.
