# D2 — Verbatim BMT Defn 3.4 and Bigraph Defn 1.4 (Mode E)

**Date:** 2026-05-15
**Mode:** E — verbatim PDF quotes with page citations.
**Method:** `pypdf` direct extraction.
**Sources:**
- `C:/Users/Nate/OneDrive/Documents/closure hunt/bmt_independence_2023.pdf` (Arizmendi–Mendoza–Vázquez-Becerra, arXiv:2309.04123)
- `C:/Users/Nate/OneDrive/Documents/closure hunt/bigraph_independence_mixture_2026.pdf` (Gilliers–Jekel, arXiv:2601.15215)
- Extraction artifacts: `C:/tmp/bmt_p4_p6.txt`, `C:/tmp/bmt_p7_p13.txt`, `C:/tmp/bigraph_p1_p7.txt`

---

## 1. BMT (Arizmendi–Mendoza–Vázquez-Becerra 2023, arXiv:2309.04123)

### 1.1 Defn 2.8 — kerG[i] (p6)

> "**Definition 2.8.** The kernel of a function i : S → V subordinated to a digraph G = (V, E) is the partition of S determined by the equivalence relation k ∼ k′ if and only if i_k = i_{k′} and (i_ℓ, i_k) is an edge of G whenever i_k ≠ i_ℓ and either k < ℓ < k′ or k′ < ℓ < k. We denote this partition by ker_G[i]."

(BMT p6.)

### 1.2 Defn 3.3 — Boolean / Monotone / Tensor independence (p8–9)

> "**Definition 3.3.** Let (A, ϕ) be a non-commutative probability space. A family (A_i)_{i∈I} of subalgebras of A is said to be
>
> (B) **boolean independent** if for any alternating product a_1 ··· a_m of elements of (A_i)_{i∈I} we have ϕ(a_1 ··· a_m) = ϕ(a_1) ··· ϕ(a_m)
>
> (M) **monotone independent** if I has a linear order < and for any alternating product a_1 ··· a_m of elements of (A_i)_{i∈I} with a_j ∈ A_{i_j} we have
>   - **M.1** ϕ(a_1 ··· a_m) = ϕ(a_k) ϕ(a_1 ··· a_{k−1} a_{k+1} ··· a_n) if i_{k−1} < i_k > i_{k+1} for some k ∈ [2, m − 1]
>   - **M.2** ϕ(a_1 ··· a_m) = ∏_{ℓ=1}^{m} ϕ(a_ℓ) if i_1 > ··· > i_{k−1} > i_k < i_{k+1} < ··· < i_m for some k ∈ [m].
>
> (T) **tensor independent** if for any (not necessarily alternating) product a_1 ··· a_m of elements of (A_i)_{i∈I} with a_j ∈ A_{i_j} we have ϕ(a_1 ··· a_m) = ∏_{B ∈ ker(i)} ϕ(→∏_{k∈B} a_k)."

(BMT p8–9.)

### 1.3 Defn 3.4 — BMT independence (p9) — LOAD-BEARING

> "**Definition 3.4.** Let (A, ϕ) be a non-commutative probability space. Suppose (A_i)_{i∈I} is a family of sub-algebras of A and G = (I, E) is a digraph on the set of indices I. The family (A_i)_{i∈I} is said to be **BMT independent with respect to the pair (ϕ, G)** if for every integer m ≥ 1 and variables a_1 ∈ A_{i_1}, a_2 ∈ A_{i_2}, ..., a_m ∈ A_{i_m} we have
>
>   ϕ(a_1 a_2 ··· a_m) = ∏_{B ∈ ker_G[i]} ϕ[(a_k)|_B]."

(BMT p9 Defn 3.4.)

### 1.4 Prop 3.5(iii) — monotone digraph (p9)

> "(iii) the algebras (A_i)_{i∈I} are monotone independent if I has a total order < and G is the digraph of <, i.e., (i, j) is an edge of G if and only if j < i."

(BMT p9 Prop 3.5(iii).)

**Notational consequence.** Under the "monotone digraph" with linearly ordered I = N and the convention "(i, j) ∈ E iff j < i": the edge points from the LARGER index to the SMALLER. For Defn 2.8's edge-test "(i_ℓ, i_k) ∈ E whenever k < ℓ < k′", we need i_k < i_ℓ.

**No iid-copies anywhere.** Defn 3.4 is a moment-formula directly on a family (A_i)_{i∈I} of subalgebras with a digraph; no requirement that distinct A_i be iid copies. This is the load-bearing point auditing Claim 5.

---

## 2. Bigraph independence (Gilliers–Jekel 2026, arXiv:2601.15215)

### 2.1 Defn 1.1 — Bigraph (p3)

> "**Definition 1.1.** A bigraph is a triple G = (V, E_1, E_2) where
> • V is a (finite or countable) set of vertices.
> • E_1, E_2 ⊆ V × V are edge sets of two types.
> • E_1 is reflexive: ∆ := {(v, v) : v ∈ V} ⊆ E_1.
> • E_2 is irreflexive: E_2 ∩ ∆ = ∅.
> • E_2 is symmetric: (v, w) ∈ E_2 ⇔ (w, v) ∈ E_2.
>
> We also write E_1^op := {(v, w) ∈ V^2 : (w, v) ∈ E_1}."

(Bigraph p3 Defn 1.1.)

### 2.2 Defn 1.2 — Admissible partitions P(c, G) (p3)

> "**Definition 1.2.** Given a bigraph G, k ∈ N, and a labeling c : [k] → V, define P(c, G) ⊆ P(k) to be the set of partitions π satisfying:
>   (1) If i ∼_π j, then c(i) = c(j).
>   (2) If i_1 < j < i_2 and i_1 ∼_π i_2, then (c(i_1), c(j)) ∈ E_1.
>   (3) If i_1 < j_1 < i_2 < j_2 and i_1 ∼_π i_2 and j_1 ∼_π j_2 and i_1 ≁_π j_1, then (c(i_2), c(j_1)) ∈ E_2."

(Bigraph p3 Defn 1.2.)

### 2.3 Defn 1.4 — Bigraph independence (p4) — LOAD-BEARING

> "**Definition 1.4 (Bigraph independence).** Let G = (V, E_1, E_2) be a bigraph. In a non-commutative probability space (A, φ), a family of subalgebras (A_v)_{v∈V} is **G-independent** if for every k ∈ N, labeling c : [k] → V, and elements a_j ∈ A_{c(j)},
>
>   φ(a_1 ··· a_k) = ∑_{π ∈ P(c, G)} K^free_π (a_1, ..., a_k),                         (1.1)
>
> where K^free_π denotes the partitioned free cumulant (see [NS06])."

(Bigraph p4 Defn 1.4.)

### 2.4 Prop 1.5 — pairwise relations encoded by (E_1, E_2) (p4)

> "**Proposition 1.5.** Consider the situation of Definition 1.4. Let v, w ∈ V with v ≠ w. Then A_v and A_w are:
> (1) Boolean independent if (v, w) ∈ E_Bool := V × V \ (E_1 ∪ E_1^op);
> (2) monotone independent if (v, w) ∈ E_mono := E_1 \ E_1^op;
> (3) anti-monotone independent if (v, w) ∈ E_{mono}^op = E_1^op \ E_1;
> (4) freely independent if (v, w) ∈ E_free := (E_1 ∩ E_1^op) \ E_2;
> (5) tensor independent if (v, w) ∈ E_ten := E_1 ∩ E_1^op ∩ E_2;
> (6) if E_1 = V × V, then (A_v)_{v∈V} are graph-independent (ε-independent) in the sense of [Mło04] with respect to (V, E_2);
> (7) if E_2 ⊇ (E_1 ∩ E_1^op) \ ∆, then (A_v)_{v∈V} are BMT-independent in the sense of [AMVB25] with respect to the digraph (V, E_1 \ ∆)."

(Bigraph p4 Prop 1.5.)

### 2.5 Hilbert space iid copies — only the MODEL, not the relation (p4)

> "This definition is not vacuous: in Section 3 we constructG-independent copies of arbitrary C*-probability spaces."

(Bigraph p4, sentence immediately after Defn 1.4.)

**Critical reading.** The iid-copies construction is the **Hilbert space model** (Prop 1.9, p4) which CONSTRUCTS G-independent copies of arbitrary C*-probability spaces — i.e., a constructive existence proof. The DEFINITION 1.4 itself is a moment-relation between SUBALGEBRAS (A_v) of a fixed (A, φ); it does NOT require them to be iid copies of anything.

---

## 3. Cross-reference

| Source | Defn / page | Moment formula | iid-copies in defn? |
|---|---|---|---|
| BMT p9 Defn 3.4 | direct | ϕ(a_1 ... a_m) = ∏_{B ∈ ker_G[i]} ϕ[(a_k)\|_B] | NO |
| Bigraph p4 Defn 1.4 | direct | φ(a_1 ... a_k) = ∑_{π ∈ P(c,G)} K^free_π(a_1,...,a_k) | NO |

Both are direct moment-relation definitions on a family of subalgebras. The auditor's claim that D2 misread these as iid-copy-dependent is confirmed by the verbatim text.

---

## 4. Files

- This file: `C:/Collatz/D2_BMT_BIGRAPH_VERBATIM.md`
- Extraction:
  - `C:/tmp/bmt_p4_p6.txt` (BMT pages 4–6, Defns 2.1–2.10, kerG)
  - `C:/tmp/bmt_p7_p13.txt` (BMT pages 7–13, Defns 3.1–3.7, Prop 3.5, Thm 3.9)
  - `C:/tmp/bigraph_p1_p7.txt` (Bigraph pages 1–7, Defns 1.1, 1.2, 1.4, Prop 1.5)
