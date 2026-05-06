# Reading summary — Atkinson 1963 + EFP 2019

Phase 1 of the Atkinson Rota-Baxter probe. Sources read:
- **Atkinson 1963**, "Some Aspects of Baxter's Functional Equation", J. Math. Anal. Appl. 7, 1-30. Located at `C:\Collatz\Baxter Spitzer\1-s2.0-0022247X63900751-main.pdf`. Sections II, V, VI (eq. 6.3/6.4 example), VII (discrete Wiener-Hopf), and IX (operator identities) read carefully.
- **Ebrahimi-Fard & Patras 2019**, "From Iterated Integrals... to Hopf and Rota-Baxter Algebras", arXiv:1911.08766v1. Located at `C:\Collatz\Baxter Spitzer\1911.08766v1.pdf`. Section 4 (Rota-Baxter algebras), specifically Definition 44 and Lemma 45.

## Rota-Baxter axiom (concrete form)

Atkinson eq. (1.1), equivalent to EFP Definition 44: a linear operator T : R → R on an associative algebra R is **Rota-Baxter of weight θ** if for all u, v ∈ R,

  T(u)·T(v) = T( T(u)·v + u·T(v) − θ·u·v )

The framework's natural target is θ = 1 (matches Atkinson eq. (5.2)/(7.2), the discrete Wiener-Hopf case).

**Key lemma (EFP 45):** if T satisfies the weight-θ axiom, so does T̃ := Id − T (companion operator). They satisfy mixed identities (EFP eq. 16). Atkinson denotes T as (+) and T̃ = U as (−).

## Discrete Wiener-Hopf structure (Atkinson Section VII)

The canonical example — the algebraic frame the brief asks us to test:

- **Algebra R**: doubly-infinite power series Σ_{n=−∞}^{∞} α_n ζ^n with absolutely convergent coefficients on the unit circle |ζ| = 1.
- **Operator T**: projection onto nonnegative powers, T(Σ α_n ζ^n) = Σ_{n≥0} α_n ζ^n. T² = T (projector).
- **Identity e**: the constant function 1 (= ζ^0).
- **Critical structural fact:** T·R and (E−T)·R are both **subalgebras** of R (since the product of two power series in nonneg powers is again in nonneg powers; same for nonpositive). This + projector property automatically gives the weight-1 RB axiom (Atkinson eq. (5.2)).

This is the algebraic frame to map onto the framework's character algebra: a projector T splitting R into two subalgebras under pointwise multiplication.

## Winding-number index and components

For z ∈ R invertible (Atkinson's S = {z ∈ R : z⁻¹ exists}, i.e., z ≠ 0 on unit circle by Wiener's theorem),

  ind z = (2π)⁻¹ { arg z(e^{2πi}) − arg z(e^0) } ∈ ℤ

is the change in arg z(ζ) as ζ traces a positive circuit of the unit circle. The connected components of S are exactly S_m = {z ∈ S : ind z = m} for m ∈ ℤ, and each S_m contains ζ^m as a representative.

**Theorem 5 application:** the homogeneous Wiener-Hopf equation T(φz) = 0 with φ ∈ T·R has −m linearly independent solutions if ind z = m < 0; no nontrivial solutions if m ≥ 0.

For finite-dimensional analogs: the winding number takes only finitely many values, so any "Atkinson-style" components on a finite character group must reduce to a finite index set (e.g., ℤ/N for the dual of a cyclic group of order N).

## Eigenvalue claim "1−q, 1−q², …" — IMPORTANT QUALIFICATION

Section IX (final paragraph): "Assuming no elements to be nilpotent, eigenvalues will form sequences of the form 1−q, 1−q², ...(cf. (6.3) and (6.4))."

Crucially, the references (6.3)-(6.4) are NOT in Section VII. They appear in **Section VI**, where the construction is:

  V = a homomorphism on R defined by Vf(ξ) := f(qξ) for fixed real q ∈ (0,1)
  T = (E − V)⁻¹ when defined

This is Atkinson's example in eq. (5.3)-(5.4) carried into Section VI. T here is **not a projector** — it's a resolvent of a homomorphism. Acting on x^n (where x is the function f(ξ) = ξ), one has Tx^n = (x^n − q^n·x^n)/(1−q^n) ... etc, with eigenvalue structure 1−q^n on (E−T)-related operators.

**The 1−q^n eigenvalue sequence is therefore a feature of the homomorphism-resolvent class of RB operators (Atkinson Section VI), NOT of the projector class (Atkinson Section VII / discrete Wiener-Hopf).** Mapping the framework's empirical spectrum onto 1−q^n requires the RB operator to be of the (E−V)⁻¹ form, not a Wiener-Hopf projector.

For the framework's R_k singular-value band σ ∈ [0.488, 0.671]: if these were 1−q^n for n=1,2,..., the band-bottom σ_min ≈ 0.488 would require q ≈ 0.512 with n=1 giving 0.488 — but then n=2 gives 1−0.262 = 0.738, n=3 gives 0.866, etc. The empirical band doesn't fit the 1−q^n sequence. (The band is densely populated, not a discrete geometric sequence.)

## Implication for Phase 2 construction

Two distinct RB-operator regimes from Atkinson:

1. **Projector regime (Section VII):** T² = T, T·R and (E−T)·R both subalgebras. On a finite commutative semisimple algebra C^N (such as the character algebra of (Z/3^k)*), pointwise-multiplicatively closed subalgebras correspond exactly to subsets A ⊆ {1,...,N}: R_A = {f : f|_{A^c} = 0}. Every such subset gives a valid RB projector T_A. The axiom is automatic; the structure is not interesting because spectrum of T_A is {0, 1} only.

2. **Homomorphism-resolvent regime (Section VI):** T = (E−V)⁻¹ for a homomorphism V. This is where 1−q^n eigenvalues appear. For the framework, a candidate is K_k itself (which is Markov-stochastic and NOT a homomorphism in general) or a specifically-chosen homomorphism on the character algebra.

The framework's K_k is row-stochastic, not a multiplicative homomorphism on the function algebra. The map K_k(f·g) ≠ K_k(f)·K_k(g) in general. So K_k is unlikely to fit the homomorphism-resolvent regime cleanly.

## Outcome-D candidate flagged

Atkinson's discrete Wiener-Hopf (Section VII) presupposes a ℤ-graded algebra (the power-series grading), with the projector T cutting the grading at 0. The framework's character algebra on (Z/3^k)* has a finite cyclic dual group ℤ/N where N = 2·3^{k-1}, which is a *finite* analog of ℤ but with no canonical "positive half." Any choice of half-set A ⊂ ℤ/N gives a valid RB projector by the trivial argument above, but the resulting structure is informationally void (spectrum {0,1}, no winding number beyond {0, ..., N−1}).

This suggests **Outcome D is plausible**: Atkinson's discrete Wiener-Hopf doesn't directly apply to finite cyclic state spaces in a non-trivial way. The structure is too easy to achieve (every subset works) and too poor (only {0,1} spectrum). The brief's hope of getting a 1−q^n eigenvalue match would require a different RB regime — Section VI's homomorphism-resolvent — and the framework's K_k / R_k aren't of that type.

Phase 2 will test this empirically: try the projector candidates (which should all trivially work) and try K_k as a candidate non-projector RB (which probably won't work), and report.
