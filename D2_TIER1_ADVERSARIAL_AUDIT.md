# D2 — Tier 1 Adversarial Audit

**Date:** 2026-05-15
**Mode:** E (verbatim PDF quotes, page references)
**Method:** `pypdf` direct extraction from PDFs in `C:/Users/Nate/OneDrive/Documents/closure hunt/`. No D2 agent summary trusted; all quotes verified.
**Subject:** `D2_TIER1_DISPOSITION.md`, `D2_TIER1_TAXONOMY.md`, `D2_TIER1_FIT_CHECK.md`.

---

## Claim 1 audit — anti-monotone iid-copies

### 1.1 Verbatim Hasebe 2010 Defn 1.9(ii) (p5)

> "(ii) The antimonotone product of {(A_i, ϕ_i)}_{i∈I}, denoted ⊳_{i∈I}(A_i, ϕ_i), is the algebraic probability space ⊲_{i∈I^op}(A_i, ϕ_i). In other words, the defining rule for (A, ϕ) = ⊳_{i∈I}(A_i, ϕ_i) is the same as (1.2) except that 'local maximum' is replaced with 'local minimum'."

(Hasebe 2010 arXiv:1009.1505v3, p5.)

**Verdict:** D2 agent's quote is **verbatim accurate**. The definition is "monotone on Λ^op" — same operational mechanic with peak ↔ valley relabel.

### 1.2 Order-flip equivalence to monotone — verbatim source

**Hasebe 2010, p2 introduction:**
> "Boolean independence [7, 35], monotone independence [26] and antimonotone independence (essentially equivalent to monotone independence)."

**GHS 2017 p1 (additional source the D2 agent flagged):**
> "the rule for anti-monotonic independence is essentially the same as the one for monotonic independence upon reversing the order structure on random variables."

**Hasebe monograph p28 (verified verbatim by my own pypdf pull, p28 §2.5 "Notes"):**
> "a) The antimonotone product
>   (φ_1 ◁ φ_2)(a_1 a_2 · · · a_n) := [ Π_{k: a_k ∈ A_1} φ_1(a_k) ] · φ_2( →Π_{k: a_k ∈ A_2} a_k ),
> which is just the flip of the monotone product and is essentially the same."

**Verdict:** Three independent verbatim sources confirm "essentially the same." D2 agent's framing is accurate.

### 1.3 Possible loopholes: dot-operation analog, phase direction

**Loophole A — dot operation under flip.** Hasebe monograph p32 Theorem 3.10 / Remark 3.11 explicitly construct **monotonically iid copies** via the natural embedding ι_i : A → ▷_{i∈N}(A_i, φ_i). The order flip just replaces ▷ with ⊳; iid copies still get built via the SAME mechanism on Λ^op. **No relaxation.** Antimono has its own dot-operation analog with the same iid-copies requirement.

**Loophole B — phase direction.** The user's adversarial question (Q3) asks whether Syracuse's phase x_{j_2} = 3^{2j_2−2} · 2^{−b_{[1,j_2]}}, depending on PRIOR accumulators b_{[1,j_1]} for j_1 < j_2, matches anti-monotone (LATER-accumulator dependence) rather than monotone. **Answer:** the dependence direction in monotone vs anti-monotone is a **labeling convention** on Λ, NOT a structural distinction — if Syracuse natively builds phase from PRIOR accumulators, flip the index labeling and you have anti-monotone. But the row (f) failure is at the iid-copies LEVEL (row j), which is order-independent. So even the correct directionality assignment can't rescue row (f).

**Audit finding:** Claim 1 **SURVIVES** intact. The "essentially the same" framing is solid across three published sources. No dot-operation loophole. No phase-direction loophole. The iid-copies failure is genuinely order-flip-invariant.

---

## Claim 2 audit — bi-monotone iid-copies

### 2.1 Verbatim GHS 2017 p13 (verified)

> "A crucial ingredient in defining the monotonic cumulants is the dot operation introduced in [12, 13]. Roughly speaking, if a_1, ..., a_n are random variables in a non-commutative probability space (A, φ), then for N ≥ 1, (N.a_1, ..., N.a_n) denotes the tuple (a_1^(1) + · · · + a_1^(N), ..., a_n^(1) + · · · + a_n^(N)), where {a_1^(i), ..., a_n^(i)}_{i=1}^N are **identically distributed and monotonically independent** with respect to φ. ... In the pairs of algebras setting, we introduce a dot operation as follows.
>
> Let (A, φ) be a non-commutative space and let A_ℓ and A_r be subalgebras of A. Take copies A^(i) = A_ℓ ⊔ A_r and φ^(i) = φ|_{A^(i)} for i ≥ 1, and let Ã = ⊔_{i≥1} A^(i) and φ̃ = ⊳⊳_{i≥1} φ^(i). ... By construction, for every χ: {1,...,n} → {ℓ, r}, the **two-faced families {({a_p^(i)}_{p∈χ^{-1}(ℓ)}, {a_q^(i)}_{q∈χ^{-1}(r)})}_{i=1}^∞ are identically distributed and bi-monotonically independent** with respect to φ̃."

(Gu-Hasebe-Skoufranis 2017, p13 §3.1.)

**Verdict:** D2 agent's quote is **verbatim accurate**.

### 2.2 Two-face structure

**Defn 2.4 (p6):**
> "A linearly ordered family {(A_{k,ℓ}, A_{k,r})}_{k∈K} of pairs of algebras in a non-commutative space (A, φ) is said to be bi-monotonically independent (of type I) with respect to φ if the joint distributions of {(A_{k,ℓ}, A_{k,r})}_{k∈K} with respect to φ and ⊲⊲_{k∈K} φ|_{A_{k,ℓ} ⊔ A_{k,r}} coincide."

Each k indexes a PAIR (A_{k,ℓ}, A_{k,r}). Syracuse has a single A_j = B⟨X̃_j⟩, not a pair.

### 2.3 Possible loopholes

**Loophole A — embed Syracuse as a two-face system with conditional B-action as the "right" face.** User's Q7. Tested: GHS Defn 2.1 (p5) uses c-bi-free convolution on the unitization, which has δ_1 (delta state) playing a specific load-bearing role for the LEFT algebra. The "right" face plays no analogous role — both faces receive ϕ̃_2 = ϕ̃_2 (the SAME state) in the c-bi-free realization. Setting A_{j,r} = B-action does not generate new structural content — it collapses back to monotone over A_{j,ℓ}.

**Loophole B — moments-based definition not requiring iid copies.** User's Q8. The DEFINITION of bi-monotonic independence (Defn 2.4) is "joint distributions coincide" — purely a moment-relation between distinct subalgebra-pairs. **The iid-copies requirement only enters when defining CUMULANTS via the dot operation.** This is a real distinction: the relation itself can be checked between distinct algebra pairs without iid copies. **BUT Syracuse has only ONE algebra-pair per step, with a single X̃_j as generator — to apply Defn 2.4 to Syracuse you need a family {(A_{j,ℓ}, A_{j,r})}_{j∈K} where each pair has a single fixed operator. The Defn 2.4 relation between distinct pairs at distinct j is the level-graded reading already used in H1, which gave row (f) = 0 via monotone peak rule (the bi-monotone Defn 2.1 collapses to monotone for trivial right algebra).** So the moments-based bi-monotone relation reduces to monotone for Syracuse — row (f) fails identically.

**Audit finding:** Claim 2 **SURVIVES** with a refinement: the iid-copies failure is a CUMULANT-theoretic obstruction (Defn 3.1 (3) extensivity uses N.a). The INDEPENDENCE-RELATION failure for Syracuse is structurally upstream (no natural two-face structure → degenerates to monotone, row (f) fails by the peak rule). Both modes confirmed.

---

## Claim 3 audit — indented / α-free iid-copies

### 3.1 Verbatim Hasebe 2010 p7 fn 1 (verified)

> "(ψ_1, θ_1) ⋋ (ψ_2, θ_2) = (ψ_1 θ_1 *^{ψ_2} ψ_2, θ_1 θ_1 *^{ψ_2} θ_2), (1.13)
> ¹ This product was formerly named indented product."

(Hasebe 2010 arXiv:1009.1505v3, p7, footnote 1 to equation (1.13).)

**Verdict:** D2 agent's quote is **verbatim accurate**. "Indented = α-free product" historical identification is sourced correctly.

### 3.2 c-monotone Defn 8.2 (Hasebe 2011)

> "Definition 8.2. In a C*-algebraic probability space (A, ϕ, ψ) equipped with two states, we say that a self-adjoint operator X has a c-monotone infinitely divisible distribution if for any n ≥ 1 there exist a C*-algebraic probability space (A_n, ϕ_n, ψ_n) and **identically distributed, c-monotone independent random variables X_1, ..., X_n ∈ A_n** such that X has the same distribution as X_1 + · · · + X_n with respect to the two states."

(Hasebe 2011 arXiv:0907.5473, p23 Defn 8.2.)

**Verdict:** D2 agent's quote is **verbatim accurate**.

### 3.3 Two-state loophole assessment (Q11)

User's adversarial question: could Syracuse embed with ϕ = E_{B_∞} and ψ = E_{B_marginal} as two states, working around iid-copies?

**Hasebe 2010 p23 Defn 4.5 — spreadability cumulants** explicitly use the dot operation N.a = a^(1) + ... + a^(N) on the spreadability system S = (U, ϕ̃, (ι^(i))_{i=1}^∞), where ι^(j) are "i.i.d. copies" (Remark 4.4 p21: "{a^(i)}_{i∈N} is meant to be 'i.i.d. copies of a'"). The α-free spreadability system §4.5 p27 explicitly constructs: "We take countable copies (A_n, ϕ_n, ψ_n, θ_n) := (A, ϕ, ψ, θ) for n ∈ N and set (U_F, ϕ̃α, ψ̃β, θ̃γ) := ⋋_{n∈N} (A_n, ϕ_n, ψ_n, θ_n)."

**Even with TWO or THREE states**, the cumulant theory and dot operation use **identically distributed copies of the FULL TUPLE (A, ϕ, ψ, θ)**. Adding a second/third state does NOT relax the iid-copies architecture — it just requires the copies to be iid in EACH state simultaneously.

**Audit finding:** Claim 3 **SURVIVES**. The two-state escape hatch fails: Hasebe 2010 §4.5 explicitly constructs iid copies even with multiple states. Adding states doesn't help with the single-X̃_j-per-step issue.

---

## Claim 4 audit — cross-cutting "outside universal-product framework"

### 4.1 Is the classification complete?

**Hasebe monograph p29 (verified):**
> "Speicher [146], and then Ben Ghorbal and Schürmann [28], showed that the associative products of linear functionals on the free product algebra with some conditions are only tensor, free and Boolean products. Muraki dropped one assumption of Ben Ghorbal and Schürmann, and as a result, the classification list contained two more products: monotone and antimonotone [123]. Muraki [124], Gerhold and Lachs [75] gave further results in this direction."

**Verdict:** Classification IS complete in its category (associative universal products on the free product algebra, with the Muraki-relaxed axioms). The five natural independences are exhaustive within that category.

### 4.2 Is iid-copies intrinsic to the classification?

**This is where the D2 agent's claim deserves a refinement.** The iid-copies architecture is intrinsic to the **dot operation / cumulant theory** built atop universal products (Defn 4.2 spreadability system, Defn 4.3 dot operation, Defn 4.5 cumulants). The **universal-product axioms themselves** (associativity + the centering rule) operate on a family (A_i, φ_i)_{i∈I} of DISTINCT algebraic probability spaces — they don't presuppose iid copies *as such*; they presuppose that elements at different positions in a word come from DIFFERENT algebras in the indexed family (with `φ(a_p) → 0` substitution at peaks rendering peak-elements scalar).

**Syracuse failure mode, refined:** the issue is NOT just "no iid copies" — it's a **shared randomness** issue. The single X̃_j operator appears at MULTIPLE positions in a word (e.g., positions 1 and 3 in row (d), positions 1+3 and 2+4 in row (f)). When two positions share the same X̃_j operator, the peak substitution X̃_j → φ(X̃_j) at one position destroys the correlated factor at the other position — but the actual moment computation reflects shared randomness, not factored randomness.

This shared-randomness issue is a property of how Syracuse's A_j = B⟨X̃_j⟩ is generated, NOT a property of the universal-product classification. Universal products allow A_j to be ANY subalgebra — but the standard cumulant theory (HS 2014 Defn 2.3, Defn 4.5) assumes distinct iid copies at distinct positions, which Syracuse doesn't supply.

### 4.3 B-valued extension status

The user's Q13 asks whether the universal-product classification extends to B-valued probability with relaxed iid copies. The Speicher 1998 memoir (`free-random-variables-1nbsped-082186999x.pdf` in the corpus) and the Hasebe-Saigo 2014 operator-valued monotone paper extend to B-valued, but **iid copies are RETAINED as a B-homomorphism condition (HS 2014 Defn 2.3 conditions (1)-(3))**. B-valued doesn't relax iid copies — it strengthens them.

**Audit finding:** Claim 4 **SURVIVES WITH CAVEAT**. The "outside universal-product framework" framing is technically too strong — Syracuse's A_j subalgebras CAN be embedded into a universal-product family (just don't supply iid copies). The **precise** characterization is:

> "Syracuse's X̃_j family is outside the universal-product CUMULANT framework, not outside the universal-product MOMENT-RELATION framework. The CUMULANT theory of all five natural independences (and their conditional/two-faced/three-state extensions) requires iid copies via the dot operation. Syracuse's single fixed X̃_j per step prevents application of this cumulant theory."

This is a real refinement — the cumulant level vs the moment level is a distinction the D2 disposition partially conflates.

---

## Claim 5 audit — BMT + bigraph spot-check

### 5.1 BMT (Arizmendi-Mendoza-Vázquez-Becerra 2023 arXiv:2309.04123)

**Defn 3.4 (p9, verified):**
> "Let (A, ϕ) be a non-commutative probability space. Suppose (A_i)_{i∈I} is a family of sub-algebras of A and G = (I, E) is a digraph on the set of indices I. The family (A_i)_{i∈I} is said to be BMT independent with respect to the pair (ϕ, G) if for every integer m ≥ 1 and variables a_1 ∈ A_{i_1}, ..., a_m ∈ A_{i_m} we have ϕ(a_1 a_2 · · · a_m) = Π_{B ∈ ker_G[i]} ϕ[(a_k)|_B]."

**Important loophole:** BMT independence is defined DIRECTLY via a moment formula on a family (A_i)_{i∈I} of subalgebras with a digraph G. **No iid-copies assumption appears in the DEFINITION.**

**However,** the CLT (Section 5, p20) imposes:
> "we are given a non-commutative probability space (A, ϕ) together with a sequence of random variables (a_i)_{i=1}^∞ that are identically distributed with zero mean and unit variance. Each finite sequence a_1, a_2, ..., a_N is assumed to be BMT independent with respect to a digraph G_N."

So iid-copies enter when one wants a CLT (or cumulant theory by analogy). The INDEPENDENCE RELATION itself is iid-free.

**Refinement to D2 agent's verdict:** BMT independence as a MOMENT RELATION between distinct subalgebras (A_i)_{i∈I} doesn't require iid copies. BUT Syracuse's structure has one A_j = B⟨X̃_j⟩ per step j, with X̃_j a fixed generator. The moment formula Defn 3.4 evaluates ϕ(a_1 ... a_m) as a product over kernel-blocks of factors ϕ[(a_k)|_B]. For Syracuse row (f) = ϕ(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}), with i = (j_1, j_2, j_1, j_2), kernel = {{1,3}, {2,4}}.

Under the monotone digraph (j_1 < j_2 with edge j_2 → j_1), ker_G[i] for this i requires checking: are positions 1 and 3 in the same kerG-block? Defn 2.8 requires that (i_ℓ, i_k) is an edge of G for all ℓ between them. At ℓ=2, i_ℓ = j_2 > j_1 = i_1, so we need (j_2, j_1) ∈ E. Under monotone digraph (j_2 → j_1), YES. So kerG[i] = {{1,3}, {2,4}}, and the moment factors as ϕ(X̃_{j_1}²) · ϕ(X̃_{j_2}²).

**This is potentially NON-ZERO** if X̃_{j_1}² and X̃_{j_2}² have non-zero ϕ-moments. **This is exactly what Syracuse row (f) shows non-zero (6.09×10⁻¹).**

**POSSIBLE TIER-3 HIT.** BMT under the monotone digraph predicts row (f) = ϕ(X̃_{j_1}²) · ϕ(X̃_{j_2}²) ≠ 0, matching Syracuse's non-zero observation. This contradicts the D2 verdict, which claims BMT "predicted to fail identically."

**Caveat — does it actually match Syracuse's value 6.09×10⁻¹?** Need to verify ϕ(X̃_{j_1}²) · ϕ(X̃_{j_2}²) numerically against the D1 data. If they match, BMT under monotone digraph **fits**. If they don't, the framework predicts a different number but at least the qualitative non-zero matches.

**More importantly:** the moment formula gives a DETERMINISTIC factorization at kernel-blocks, while Syracuse's actual moment depends on the SHARED X̃_j operator. The BMT formula assumes ϕ((a_k)|_B) at block B with all positions sharing index i is well-defined — for Syracuse this is ϕ(X̃_{j_1}²) for B={1,3}, which IS well-defined (single second moment of X̃_{j_1}). So the BMT formula CAN BE EVALUATED on Syracuse without needing iid copies. **The D2 agent's claim that BMT "inherits the same iid-copies architecture" is too strong for the INDEPENDENCE RELATION.**

### 5.2 Bigraph (Gilliers-Jekel 2026 arXiv:2601.15215)

**Defn 1.4 (p4, verified):**
> "In a non-commutative probability space (A, φ), a family of subalgebras (A_v)_{v∈V} is G-independent if for every k ∈ N, labeling c: [k] → V, and elements a_j ∈ A_{c(j)}, φ(a_1 ··· a_k) = Σ_{π ∈ P(c,G)} K_π^free(a_1, ..., a_k), where K_π^free denotes the partitioned free cumulant."

**Same situation as BMT.** Bigraph independence is a MOMENT RELATION between distinct subalgebras (A_v)_{v∈V}, defined via a SUM of partitioned free cumulants K_π^free over a bigraph-constrained partition set P(c, G). **No iid-copies in the definition.**

The D2 agent quoted Gilliers-Jekel p3: "enabling the construction of independent copies of arbitrary non-commutative probability spaces." This is about the **Hilbert space MODEL** (Proposition 1.9 p4) that CONSTRUCTS bigraph-independent copies — not about whether the independence relation REQUIRES iid copies to be defined.

**Application to Syracuse row (f):** with V = {j_1, j_2}, c = (j_1, j_2, j_1, j_2), partition with i_1 ∼ i_3 and j_1 ∼ j_4 → π = {{1,3}, {2,4}} (this is a CROSSING partition). Defn 1.2 condition (3): if i_1 < j_1 < i_2 < j_2 and i_1 ~π i_2 (1~3) and j_1 ~π j_2 (2~4) and i_1 ≁π j_1 (1 ≁ 2), then (c(i_2), c(j_1)) ∈ E_2. So {{1,3},{2,4}} ∈ P(c, G) iff (j_1, j_2) ∈ E_2 (E_2 = tensor edges).

Under a chosen bigraph G with appropriate E_2, the crossing pair-partition is admitted and contributes K_2^free(X̃_{j_1}, X̃_{j_1}) · K_2^free(X̃_{j_2}, X̃_{j_2}). For centered X̃_j with non-zero second moment, this is **non-zero**. **Potentially matches Syracuse row (f) qualitatively** under a tensor-edge bigraph between j_1 and j_2.

But Syracuse's row (b) is ZERO (cross-step second moment), which means in a bigraph fit, the pair partition {{1,2}, {3,4}} at the same word should NOT contribute. Under bigraph with E_2 ⊃ (j_1, j_2), the non-crossing pair partition {{1,2},{3,4}} ALSO contributes — and would give K_2^free(X̃_{j_1}, X̃_{j_2}) = 0 (cross-step free 2-cumulant assumed zero). Wait — K^free here is the partitioned free cumulant on (a_1, a_2) where a_1 ∈ A_{j_1} and a_2 ∈ A_{j_2}; under free independence of A_{j_1} and A_{j_2}, this would vanish, but under TENSOR independence (which is what E_2 ⊃ (j_1, j_2) encodes), this K_2^free might NOT vanish — it could equal Cov_{tensor}(X̃_{j_1}, X̃_{j_2}).

This needs careful evaluation, but **the qualitative loophole is real**: bigraph independence's moment formula CAN produce non-zero 4-alternating moments under tensor-edge bigraphs, without invoking iid copies. **The D2 agent's claim that bigraph "inherits iid-copies architecture explicitly" is too strong** — only the constructive Hilbert space MODEL requires iid copies, not the moment-relation DEFINITION.

### 5.3 Other Tier 2/3 candidates not requiring iid copies?

Both BMT (Defn 3.4) and bigraph (Defn 1.4) define independence directly as moment relations, NOT as iid-copies constructions. Similarly:

- **ε-independence / Λ-independence** (Mlotkowski 2004): mixture of tensor and free, defined as moment relation on a graph (bigraph paper references).
- **Tree independence** (Jekel-Liu 2020, Jekel-Ohta-Wendell 2025; in bigraph paper refs JL20, JOW25): free-Boolean-monotone mixtures defined on trees, moment-formula-based.

**Audit finding:** Claim 5 **FAILS** as stated. The D2 agent's spot-check of BMT + bigraph as "inheriting iid-copies architecture" is incorrect for the INDEPENDENCE RELATIONS themselves; iid copies enter only in cumulants/CLT or in explicit Hilbert space constructions of independent copies. **There is a genuine Tier 2/3 loophole the D2 agent missed:** BMT and bigraph as MOMENT RELATIONS may admit Syracuse with single fixed X̃_j per step. Whether they actually fit Syracuse's row (f) value 6.09×10⁻¹ needs a direct moment-formula check (one row of arithmetic).

---

## Claim 6 audit — next-arc recommendations

### 6.1 Spreadability iid-status

**Hasebe 2010 §4.3 (verified):**
> "Definition 4.2. Let (A, ϕ) be a unital algebraic probability space. A unital spreadability system for (A, ϕ) is a triplet S = (U, ϕ̃, (ι^(j))_{j=1}^∞) satisfying the following properties: (i) (U, ϕ̃) is a unital algebraic probability space, (ii) ι^(j): A → U is a unital homomorphism such that ϕ = ϕ̃ ∘ ι^(j) for each j ≥ 1. For simplicity, ι^(j)(a) is denoted by a^(j), a ∈ A. (iii) [spreadability identity (4.1)] ..."
>
> "Definition 4.3. (i) The dot operation for a ∈ A is the sum of copies of a: N.a := a^(1) + a^(2) + ··· + a^(N) ∈ U, N ∈ N; 0.a := 0."
>
> "Remark 4.4. Roughly, {a^(i)}_{i∈N} is meant to be 'i.i.d. copies of a' ..."

**Verdict:** Spreadability systems explicitly DO use iid copies via Defn 4.2(ii)+(iii). The "spreadability" terminology hides the iid-copies requirement under the embedding ι^(j) — but ι^(j) are unital homomorphisms A → U with the same φ-image, and the spreadability identity (4.1) is **exactly** an exchangeability-like condition on those copies.

**D2 agent recommendation to chase Hasebe-Lehner 2023 spreadability as iid-free is OPTIMISTIC at best.** The spreadability framework as published in Hasebe 2010 §4.3 still requires the iid-copy embeddings.

### 6.2 Popa 2008 / operator-valued without B-homomorphic copies

This is in the closure-hunt corpus only as a future-arc reference (Popa 2008 not pulled). Cannot audit without the PDF.

**Audit finding:** Claim 6 **SURVIVES_WITH_CAVEAT**. Spreadability is NOT iid-free (D2's claim that it "may or may not require iid copies — needs verification" is now verified: it DOES require them). Popa 2008 unaudited.

---

## Verdict: **SURVIVED_WITH_CAVEAT**

The D2 Tier 1 disposition's per-candidate verdicts are **largely correct** for the iid-copies-of-A_λ assumption embedded in Tier 1 CUMULANT theories (Hasebe 2010 Defn 4.5 spreadability cumulants, HS 2014 Defn 2.3 monotone cumulants, GHS 2017 §3.1 bi-mono cumulants, Hasebe 2011 Defn 8.2 c-mono iid).

**Per-claim verdicts:**

| Claim | Verdict | Notes |
|---|---|---|
| 1 (anti-mono iid) | SURVIVED | 3 verbatim sources concur; no loopholes |
| 2 (bi-mono iid) | SURVIVED | GHS p13 verbatim; two-face structure & moment-only loopholes fail |
| 3 (indented/α-free iid) | SURVIVED | Hasebe 2010 p7 fn1 + Defn 8.2 verbatim; two-state escape closed by §4.5 multi-state iid construction |
| 4 (cross-cutting "outside universal-product") | **SURVIVED_WITH_CAVEAT** | Refinement: Syracuse is outside the universal-product CUMULANT framework, not outside the moment-relation framework |
| 5 (BMT + bigraph spot-check) | **FAILED** | BMT Defn 3.4 + bigraph Defn 1.4 are MOMENT relations not requiring iid copies in the relation itself; iid copies only enter in cumulants/CLT/Hilbert-space-realization. **Genuine loophole missed by D2 agent.** |
| 6 (spreadability + Popa) | SURVIVED_WITH_CAVEAT | Spreadability NOT iid-free (Hasebe 2010 §4.3 explicit); Popa 2008 unaudited |

---

## Specific errors / loopholes found

### Error 1 — D2 conflates cumulant-level and moment-level iid-copies requirements

**Quote from D2_TIER1_DISPOSITION.md §2:**
> "Cross-cutting prediction. Any Tier-N framework constructed by universal-product axioms (= any framework derived from the Muraki-Speicher-BenGhorbal classification) will share the iid-copies feature and therefore fail at row (f) for the SAME reason."

This is technically inaccurate. The universal-product axioms (associative product on the free product algebra) operate on a family of DISTINCT subalgebras (A_i, ϕ_i)_{i∈I}, NOT on iid copies of one algebra. The iid-copies architecture is a CUMULANT-THEORETIC tool (Hasebe 2010 Defn 4.3 dot operation, HS 2014 Defn 2.3) that DERIVES cumulants from the universal-product. **The independence RELATIONS themselves (Defn 1.9, Defn 2.4, Defn 3.4 BMT, Defn 1.4 bigraph) are moment formulas that do NOT presuppose iid copies.**

### Error 2 — BMT + bigraph spot-check too dismissive (Tier 2/3 loophole)

**Quote from D2_TIER1_DISPOSITION.md §2:**
> "BMT uses pairwise independence graph atop the five natural primitives ... All listed candidates inherit the iid-copies architecture from their parent natural-independence frameworks. None resolve the structural single-X̃_j-per-step issue."

**Verbatim BMT Defn 3.4 (p9):** moment formula with NO iid copies in the definition.
**Verbatim bigraph Defn 1.4 (p4):** moment formula with NO iid copies in the definition.

For row (f) ϕ(X̃_{j_1} X̃_{j_2} X̃_{j_1} X̃_{j_2}) under BMT monotone digraph (j_2 → j_1): kerG[i] = {{1,3}, {2,4}} → moment factors as ϕ(X̃_{j_1}²) · ϕ(X̃_{j_2}²) which is **non-zero**, qualitatively matching Syracuse's 6.09×10⁻¹. This **may or may not** quantitatively fit Syracuse, but the D2 agent's "predicted to fail identically" is **not justified by reading BMT Defn 3.4 directly**.

### No errors in the per-Tier-1-candidate iid-copies analysis

The Tier 1 candidates (anti-mono, bi-mono, indented/α-free) **do** all use the iid-copies dot operation for their cumulant theories, and **do** fail Syracuse row (f) when one tries to apply the cumulant-level analysis. The D2 agent's per-candidate Tier 1 verdicts are correct.

The errors are scope/conflation errors at the Tier 2/3 boundary and at the cross-cutting recommendation.

---

## Recommended follow-up

1. **Direct moment-formula check of BMT Defn 3.4 against Syracuse row (f).** Compute ϕ(X̃_{j_1}²) · ϕ(X̃_{j_2}²) numerically and compare to 6.09×10⁻¹. If quantitatively close (within sampling noise), BMT fits Syracuse — this would be a **Tier 2 disconfirmation of the D2 disposition**. Effort: 30 min using D1's numerical data.

2. **Direct moment-formula check of bigraph Defn 1.4 against Syracuse rows (b), (d), (f).** Bigraph is more flexible than BMT (supports all 5 natural independences as pairwise relations + edge sets E_1, E_2). With careful choice of bigraph G, Syracuse may fit row (f) non-zero while keeping row (b) zero. Effort: 1-2 hours focused.

3. **Re-frame D2 disposition** to distinguish:
   - **Tier 1 (anti-mono, bi-mono, α-free) cumulant theories** — definitively fail Syracuse via iid-copies obstruction at the cumulant level.
   - **Tier 2/3 BMT + bigraph moment relations** — potentially fit Syracuse as moment relations without invoking iid copies. NOT YET FALSIFIED.
   - **Tier 4+ spreadability** — also requires iid copies per Hasebe 2010 §4.3, NOT a clean iid-free escape.

4. **Audit the D2 disposition's "outside universal-product framework" claim.** The correct claim is "outside the universal-product CUMULANT framework," not "outside universal-product altogether." This matters for downstream writeups.

---

## Files

- This audit: `C:/Collatz/D2_TIER1_ADVERSARIAL_AUDIT.md`
- Verified D2 sources:
  - `C:/Collatz/D2_TIER1_DISPOSITION.md`
  - `C:/Collatz/D2_TIER1_TAXONOMY.md`
  - `C:/Collatz/D2_TIER1_FIT_CHECK.md`
- PDFs verified verbatim:
  - `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_2010_three_state_independence.pdf` (p2, p5, p7 fn1, §4.3 p21-27)
  - `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_2011_conditionally_monotone.pdf` (p23 Defn 8.2)
  - `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_monograph_p28` (p28 antimono "just the flip")
  - `C:/Users/Nate/OneDrive/Documents/closure hunt/bi-monotonic_gu_hasebe_skoufranis_2017.pdf` (p1, p5 Defn 2.1, p6 Defn 2.4, p13 §3.1 dot operation)
  - `C:/Users/Nate/OneDrive/Documents/closure hunt/bmt_independence_2023.pdf` (p9 Defn 3.4 moment relation; p20 CLT iid)
  - `C:/Users/Nate/OneDrive/Documents/closure hunt/bigraph_independence_mixture_2026.pdf` (p4 Defn 1.4 moment relation; p3 Hilbert-space iid construction)
- Extraction artifacts in `c:/tmp/` (hasebe2010_p21_28.txt, ghs2017_p1_18.txt, hasebe2011_p21_28.txt, bmt_p1_10.txt, bmt_p16_22.txt, bigraph_p1_10.txt, hasebe_monograph_p26_32.txt, hasebe2010_p7_9.txt)
