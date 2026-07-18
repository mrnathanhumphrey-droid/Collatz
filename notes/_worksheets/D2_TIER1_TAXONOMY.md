# D2 — Tier 1 Taxonomy (monotone variants)

**Date:** 2026-05-15
**Mode:** E. Verbatim quotes with page references. Single-source citations from PDFs pulled to `C:/Users/Nate/OneDrive/Documents/closure hunt/`.
**Reads with:** `H1_PRIME_DISPOSITION.md`, `H1_PRIME_STRUCTURAL_ARGUMENT.md`, `D1_DISPOSITION.md`.

## 0. Scope

Three Tier 1 candidates were prescribed:

1. **Anti-monotone independence** (Muraki 2002; flip of monotone).
2. **Bi-monotone independence** (Gu-Hasebe-Skoufranis 2017, arXiv:1708.05334).
3. **Indented independence** (Hasebe 2010 three-state, arXiv:1009.1505; the former name of what is now called α-free product).

For each: extract the verbatim definition and check whether the iid-copies-of-A_λ assumption that H1' identified as load-bearing also appears.

---

## 1. Anti-monotone independence

### 1.1 Verbatim sources

**Source 1: Hasebe monograph 2026 (in closure-hunt corpus), p28.**

The antimonotone product is defined as a sibling of monotone:

> "a) The antimonotone product
>   (φ_1 ◁ φ_2)(a_1 a_2 · · · a_n) := [ Π_{k: a_k ∈ A_1} φ_1(a_k) ] · φ_2( →Π_{k: a_k ∈ A_2} a_k ),
> which is just the flip of the monotone product and is essentially the same."

(Hasebe, *Monotone Probability Theory* monograph, p28, Section 2.5 "Notes" — antimonotone listed as item (a) among four universal products that satisfy associativity besides monotone.)

**Source 2: Hasebe 2010 three-state (arXiv:1009.1505), p2 introduction.**

> "Boolean independence [7, 35], monotone independence [26] and antimonotone independence (essentially equivalent to monotone independence)."

(Hasebe 2010, p2.)

**Source 3: Hasebe 2010 three-state (arXiv:1009.1505), p5 Defn 1.9(ii).**

> "(ii) The antimonotone product of {(A_i, φ_i)}_{i∈I}, denoted ⊳_{i∈I}(A_i, φ_i), is the algebraic probability space ⊲_{i∈I^op}(A_i, φ_i). In other words, the defining rule for (A, φ) = ⊳_{i∈I}(A_i, φ_i) is the same as (1.2) except that 'local maximum' is replaced with 'local minimum'."

(Hasebe 2010, p5 Defn 1.9(ii).)

**Source 4: Gu-Hasebe-Skoufranis 2017 (bi-monotone paper), p1.**

> "the rule for anti-monotonic independence is essentially the same as the one for monotonic independence upon reversing the order structure on random variables."

(Gu-Hasebe-Skoufranis 2017, p1.)

### 1.2 iid-copies status

Anti-monotone is **the same algebraic structure as monotone, restated on the opposite-ordered index set**. The HS 2014 Defn 2.2 peak-rule (with peaks defined relative to the order) carries over verbatim under the order-reversal; HS 2014 Defn 2.3 (the iid-copies dot operation) carries over verbatim as well.

**iid-copies requirement: YES, identical to monotone.** Anti-monotone is monotone under the opposite order — it inherits HS 2014 Defn 2.3's iid-copies construction wholesale. No new degree of freedom relative to monotone.

### 1.3 Definitional features

- **Algebraic family:** linearly ordered family (A_λ)_{λ ∈ Λ}, but with the opposite order Λ^op.
- **Peak rule:** at index i where λ_{i−1} > λ_i < λ_{i+1} ("valley" in original order = "peak" in opposite order), substitute X_i → φ(X_i).
- **iid copies:** required for the dot operation N.X = X^(1) + ... + X^(N), where {X^(j)} are antimonotone-independent identically distributed copies. (Mirror of HS 2014 Defn 2.3.)

---

## 2. Bi-monotone independence

### 2.1 Verbatim source

**Gu-Hasebe-Skoufranis 2017 (arXiv:1708.05334), p1 abstract.**

> "In this article, the notion of bi-monotonic independence is introduced as an extension of monotonic independence to the two-faced framework for a family of pairs of algebras in a non-commutative probability space."

**Defn 2.1 (p5, type I), Defn 2.4 (p6, bi-monotonic independence).**

> "Definition 2.1. Let (A_{1,ℓ}, A_{1,r}) and (A_{2,ℓ}, A_{2,r}) be two pairs of algebras with linear functionals φ_k: A_{k,ℓ} ⊔ A_{k,r} → C. The bi-monotonic product of φ_1 and φ_2, denoted φ = φ_1 ⊲⊲ φ_2, is the linear functional on (A_{1,ℓ} ⊔ A_{1,r}) ⊔ (A_{2,ℓ} ⊔ A_{2,r}) defined as follows. Let Ã_{k,ℓ} and Ã_{k,r} be the unitizations of A_{k,ℓ} and A_{k,r} respectively, let φ̃_k be the unique unital linear extension of φ_k to Ã_{k,ℓ} * Ã_{k,r} ≃ Ã_{k,ℓ} ⊔ A_{k,r}, let δ_1 be the delta functional on Ã_{1,ℓ} ⊔ A_{1,r}, and let (φ̃, ψ̃) = (φ̃_1, δ_1) ** (φ̃_2, φ̃_2) be the c-bi-free product on Ã_{1,ℓ} ⊔ A_{1,r} * Ã_{2,ℓ} ⊔ A_{2,r} ≃ (̃A_{1,ℓ} ⊔ A_{1,r}) ⊔ (A_{2,ℓ} ⊔ A_{2,r}). Then φ is defined to be the restriction of φ̃ to (A_{1,ℓ} ⊔ A_{1,r}) ⊔ (A_{2,ℓ} ⊔ A_{2,r})."

(Gu-Hasebe-Skoufranis 2017, p5 Defn 2.1.)

**Defn 2.4 (p6):**

> "A linearly ordered family {(A_{k,ℓ}, A_{k,r})}_{k∈K} of pairs of algebras in a non-commutative space (A, φ) is said to be bi-monotonically independent (of type I) with respect to φ if the joint distributions of {(A_{k,ℓ}, A_{k,r})}_{k∈K} with respect to φ and ⊲⊲_{k∈K} φ|_{A_{k,ℓ} ⊔ A_{k,r}} coincide."

(Gu-Hasebe-Skoufranis 2017, p6 Defn 2.4.)

### 2.2 iid-copies status

**Gu-Hasebe-Skoufranis 2017, p13 Section 3.1 ("The dot operation"):**

> "A crucial ingredient in defining the monotonic cumulants is the dot operation introduced in [12, 13]. Roughly speaking, if a_1, ..., a_n are random variables in a non-commutative probability space (A, φ), then for N ≥ 1, (N.a_1, ..., N.a_n) denotes the tuple (a_1^(1) + · · · + a_1^(N), ..., a_n^(1) + · · · + a_n^(N)), where {a_1^(i), ..., a_n^(i)}_{i=1}^N are **identically distributed and monotonically independent** with respect to φ. ... In the pairs of algebras setting, we introduce a dot operation as follows.
>
> Let (A, φ) be a non-commutative space and let A_ℓ and A_r be subalgebras of A. Take copies A^(i) = A_ℓ ⊔ A_r and φ^(i) = φ|_{A^(i)} for i ≥ 1, and let Ã = ⊔_{i≥1} A^(i) and φ̃ = ⊲⊲_{i≥1} φ^(i). ... By construction, for every χ: {1,...,n} → {ℓ, r}, the **two-faced families {({a_p^(i)}_{p∈χ^{-1}(ℓ)}, {a_q^(i)}_{q∈χ^{-1}(r)})}_{i=1}^∞ are identically distributed and bi-monotonically independent** with respect to φ̃."

(Gu-Hasebe-Skoufranis 2017, p13.)

**iid-copies requirement: YES, identical structure to HS 2014 Defn 2.3.** Bi-monotone requires iid copies of the pair-of-algebras structure (A_ℓ, A_r) to define cumulants. Inherits the iid-copies architecture from monotone.

### 2.3 Definitional features

- **Algebraic family:** linearly ordered family of **pairs of algebras** (A_{k,ℓ}, A_{k,r})_{k ∈ K}. Two-faced framework — each k indexes a LEFT and a RIGHT algebra.
- **Peak rule:** generalized through the c-bi-free product realization; effectively the same peak-substitution mechanic as monotone, but with the bi-free convolution structure on left/right pairs.
- **Single algebra vs paired algebra:** Syracuse has ONE algebra A_j per step j, not a pair (A_{j,ℓ}, A_{j,r}). To use bi-monotone, one would need to artificially split each A_j into a left/right pair — but Syracuse's X̃_j is a single operator, not a left/right pair.
- **iid copies:** required for cumulant theory (p13 dot operation).

---

## 3. Indented independence (= α-free product = "3-state independence")

### 3.1 Verbatim source — historical naming

**Hasebe 2010 three-state (arXiv:1009.1505), p7 footnote 1:**

> "This product was formerly named indented product."

(Hasebe 2010, p7, footnote 1 to the α-free product definition (1.13).)

So **"indented"** is the historical name (pre-2010) for what is now called **α-free product** or **α-freeness**. The α-free product is defined on triplets of states.

### 3.2 Verbatim definition (α-free / former-indented)

**Hasebe 2010, p7 equation (1.12)–(1.13):**

> "The new product can be described as
>   (φ_1, ψ_1, θ_1) ⋋ (φ_2, ψ_2, θ_2) = (φ_1 θ_1 *^{ψ_2} φ_2, ψ_1 θ_1 *^{ψ_2} ψ_2, θ_1 θ_1 *^{ψ_2} θ_2),    (1.12)
> which will be called the α-free product.¹ Although this product is made out of c-free products, the associativity of α-free product is not a consequence of the associativity of the c-free product. In particular, the product for unital 2-algebraic probability spaces
>   (ψ_1, θ_1) ⋋ (ψ_2, θ_2) = (ψ_1 θ_1 *^{ψ_2} ψ_2, θ_1 θ_1 *^{ψ_2} θ_2),    (1.13)
> ¹ This product was formerly named indented product."

(Hasebe 2010, p7.)

**Hasebe 2010 abstract, p1:**

> "We define a new independence in non-commutative probability, called α-freeness, with respect to a triplet of states. This concept unifies several independences in non-commutative probability, in particular, free, monotone, antimonotone and Boolean ones as well as conditionally free, conditionally monotone and conditionally antimonotone independences."

(Hasebe 2010, p1 abstract.)

### 3.3 iid-copies status

α-freeness is the **maximal unification** of all five natural independences (tensor, free, Boolean, monotone, antimonotone) plus their conditional variants (c-free, c-monotone, c-antimonotone). It is defined on **triplets of states (φ, ψ, θ)** per algebra.

α-freeness inherits cumulants and convolution via the same dot-operation iid-copies construction Hasebe & Saigo set up in [86, 87] for natural independence cumulants (referenced in Hasebe 2010 p1 + Section 4.5 "α-free cumulants").

**Hasebe 2010, p23 (used for c-monotone, a direct restriction of α-free):**

> "Definition 8.2. In a C*-algebraic probability space (A, φ, ψ) equipped with two states, we say that a self-adjoint operator X has a c-monotone infinitely divisible distribution if for any n ≥ 1 there exist a C*-algebraic probability space (A_n, φ_n, ψ_n) and **identically distributed, c-monotone independent random variables X_1, ..., X_n ∈ A_n** such that X has the same distribution as X_1 + · · · + X_n with respect to the two states."

(Hasebe 2011 "Conditionally monotone independence I", arXiv:0907.5473, p23 Defn 8.2.)

**iid-copies requirement: YES**, in the same shape as monotone — the cumulant theory uses identically distributed independent copies. Furthermore, α-free uses a triplet of states, all required to be identically distributed across copies.

### 3.4 Definitional features

- **Algebraic family:** linearly ordered family of triplets (A_λ, φ_λ, ψ_λ, θ_λ).
- **Peak rule:** generalized; reduces to c-monotone (or to free/Boolean/monotone) under specialization of the triplet.
- **iid copies:** required for cumulant theory.
- **Number of states per algebra:** 3 (vs monotone's 1). Syracuse has 1 conditional expectation ϕ = E_{B_∞} — there is no natural triplet of states.

---

## 4. Out-of-scope candidates noted during scan

These came up in the literature scan but are not in Tier 1; documented for reference:

| Framework | Reference | Status |
|---|---|---|
| BMT independence (Boolean/Monotone/Tensor mixture) | Arizmendi-Mendoza-Vázquez-Becerra 2023 (arXiv:2309.04123) | Tier 2/3 candidate; uses pairwise independence graph; iid-copies built into model |
| Bigraph independence (mixes ALL 5 natural) | Gilliers-Jekel 2026 (arXiv:2601.15215) | Tier 2/3 candidate; uses iid-copies via Hilbert space construction |
| Conditionally monotone (c-monotone) | Hasebe 2011 (arXiv:0907.5473) | Restriction of α-free; uses iid copies + triplet/pair of states |
| BM-independence (Wysoczański 2007) | Wysoczański 2007 | Tier 2/3; precursor to BMT |
| Λ-monotone, Λ-Boolean (Lenczewski) | Lenczewski [10] cited in BMT | Tier 2/3; mixes tensor with monotone/Boolean |

All listed candidates inherit the iid-copies architecture from their parent natural-independence frameworks. None resolve the structural single-X̃_j-per-step issue identified in H1'.

---

## 5. Summary table — iid-copies status

| Tier 1 candidate | Algebraic family | Number of states per algebra | iid-copies required? | Source quote |
|---|---|---|---|---|
| **Anti-monotone** (Muraki 2002) | (A_λ)_{λ ∈ Λ^op}, linearly ordered (opposite of monotone) | 1 | **YES** (inherits HS 2014 Defn 2.3 under order flip) | Hasebe 2010 p2; Hasebe monograph p28 |
| **Bi-monotone** (GHS 2017) | (A_{k,ℓ}, A_{k,r})_{k ∈ K}, linearly ordered pairs | 1 (two-faced via L/R structure) | **YES** (GHS 2017 p13 explicit dot operation) | GHS 2017 p13 |
| **Indented / α-free** (Hasebe 2010, formerly indented) | (A_λ, φ_λ, ψ_λ, θ_λ)_{λ ∈ Λ}, linearly ordered triplets | 3 | **YES** (via c-monotone restriction; Hasebe 2011 Defn 8.2) | Hasebe 2011 p23 Defn 8.2 |

**All three Tier 1 candidates require iid copies of A_λ-elements at different positions in a word.**

This is structurally the same load-bearing assumption that H1' identified as failing for Syracuse: Syracuse has a SINGLE X̃_j per step j, used at every occurrence of step j in a word — not independent copies (HS 2014 Defn 2.3 type). See `H1_PRIME_STRUCTURAL_ARGUMENT.md §3.5 STRUCTURAL-1`.

---

## 6. Files

- Anti-monotone source: `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_monotone_probability_theory_monograph.pdf` (p28)
- Anti-monotone source: `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_2010_three_state_independence.pdf` (p2, p5)
- Bi-monotone source: `C:/Users/Nate/OneDrive/Documents/closure hunt/bi-monotonic_gu_hasebe_skoufranis_2017.pdf` (p1, p5, p6, p13)
- Indented / α-free source: `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_2010_three_state_independence.pdf` (p1, p7)
- c-monotone iid Defn 8.2: `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_2011_conditionally_monotone.pdf` (p23)
- BMT (out of Tier 1): `C:/Users/Nate/OneDrive/Documents/closure hunt/bmt_independence_2023.pdf`
- Bigraph (out of Tier 1): `C:/Users/Nate/OneDrive/Documents/closure hunt/bigraph_independence_mixture_2026.pdf`
- HS 2014 Defn 2.3 (the iid-copies reference def): `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_saigo_2014_operator_valued_monotone.pdf` (p3)
