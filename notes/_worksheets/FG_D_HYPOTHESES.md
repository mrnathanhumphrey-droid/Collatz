# FG Candidate D — Benoist-Quint (stationary measures on homogeneous spaces)

**PDF:** Benoist_Quint_2016_Stationary_Measures_translation.pdf + Benoist_Quint_Random_Walks_Reductive_Groups_book.pdf.
**Extracted text:** `C:/tmp/fg/benoist_quint_2016.txt`, `C:/tmp/fg/benoist_quint_book.txt`.

---

## B-Q Theorem 1.1 (VERBATIM, BQ 2016 translation, p. 1)

> "Let G be a connected almost simple real Lie group, Λ a lattice in G, X = G/Λ and μ a probability measure on G with compact support, such that supp(μ) generates a Zariski-dense subgroup of G. Then any non-atomic μ-stationary Borel probability measure on X is the Haar measure on X."

### Hypotheses (typed):

- h_D.1.1.group: G **connected almost simple real Lie group**. [TYPE (i)]
- h_D.1.1.lattice: Λ lattice in G; X = G/Λ. [TYPE (i)]
- h_D.1.1.walk: μ probability on G, **compactly supported**. [TYPE (ii)]
- h_D.1.1.gen: supp(μ) generates Zariski-dense subgroup of G. [TYPE (ii)]
- h_D.1.1.stat: ν μ-stationary Borel probability on X. [TYPE (iii)]
- h_D.1.1.nonatomic: ν has no atom. [TYPE (iii)]

### Conclusion C_D.1.1:

- ν = Haar measure on X = G/Λ. [TYPE (iv) — **qualitative**, no rate]

---

## B-Q Theorem 1.3 (VERBATIM, BQ 2016 translation, p. 2)

> "Let Γ be a sub-semigroup of SL_d(Z) acting on R^d strongly irreducibly. Let μ be a measure on SL_d(Z) whose finite support generates Γ. Then any non-atomic μ-stationary probability measure on X = T^d is the Haar measure of X."

### Hypotheses (typed):

- h_D.1.3.group: Γ ⊂ **SL_d(Z)**, acting on R^d strongly irreducibly. [TYPE (i)]
- h_D.1.3.walk: μ on SL_d(Z) **finite support** generating Γ. [TYPE (ii)]
- h_D.1.3.stat: ν μ-stationary on T^d, non-atomic. [TYPE (iii)]

### Conclusion C_D.1.3:

- ν = Haar on T^d. [TYPE (iv) — qualitative]

---

## Note from BQ 2016, p. 2:

> "...this approach has the advantage of being more general; in particular, Benoist and Quint have been able to prove Theorem D without making the assumption (Γ-2). However their ergodic theoretic argument is not quantitative, certainly not in the sense of Theorem A. It also does not give equidistribution of ν^{*n} * δ_x as in Corollary B."

Also: p. 3 mentions "Theorem 1.1 and Corollary 1.2 can be extended, with no significant change to the proof, to **p-adic Lie groups G**" — this is the p-adic extension hook.

---

## Phase 1 — hypothesis × input matrix

| Hypothesis | (1)-(4) Disposition |
|---|---|
| h_D.group: G connected almost simple real Lie group | **FAILED**: Syracuse on (Z/3^n)* is abelian profinite, not simple Lie. The p-adic extension hook would replace G with a p-adic Lie group like SL_2(Q_3), but Syracuse's ambient (Z/3^n)* is *abelian* (= Z_3^* in the inverse limit, ≅ Z_3 × Z/2), NOT a simple p-adic Lie group. FAILED in both real and p-adic Lie versions. |
| h_D.lattice: Λ lattice in G | **N/A**: abelian profinite groups have no analog of "lattice in semisimple Lie group". |
| h_D.walk_compact_support | The Syracuse step distribution is Geom(2) on N+1 — *not finitely supported, not compactly supported in any ambient*. The 2-adic exponentials 2^{-v} ∈ Z_3 for v ∈ N+1 form an infinite countable subset; in the closure topology of Z_3, the closure is the orbit ⟨2⟩ ⊂ Z_3^*, but the measure itself has Geom(2) tails, **not compactly supported on (Z_3)\* in the sense BQ requires**. NEEDS_PROOF, but **likely FAILED.** |
| h_D.gen: support generates Zariski-dense subgroup | (Z/3^n)* is abelian → Zariski-closure of ⟨2⟩ is itself, which is *not Zariski-dense* in the natural ambient (a torus / 1-dim algebraic group); Zariski-density is vacuous or trivial in abelian case. FAILED. |
| h_D.stat: non-atomic stationary | π_n is Tao's stationary, which is non-atomic on Z/3^n (uniform-like with class-mass biases per R64.B). SATISFIED for each n; in the inverse limit, the stationary measure μ_∞ on Z_3^* is non-atomic. SATISFIED. |

**Phase 1: NO_FIT.** h_D.group + h_D.gen fail categorically (abelian → not simple Lie, Zariski-density vacuous).

---

## Phase 2 — conclusion shape

Even if hypotheses were satisfied, B-Q's conclusion is **qualitative**: "ν = Haar." This does not deliver a polynomial-in-A Fourier decay rate. The BFLM approach (Theorem C) does deliver such a rate (effective version), but BFLM's hypotheses fail on Syracuse for the same reasons as Theorem 1.3 here.

**CONCLUSION_SHAPE_MISMATCH**: qualitative classification result, not quantitative rate.

---

## Phase 3 — profinite extension

The "p-adic Lie group G" hook in BQ 2016 p. 3 mentions extension. But the load-bearing assumption "G connected almost simple" means G must be **non-abelian semisimple** even in the p-adic version. Concretely: BQ extends to SL_2(Q_p), SL_3(Q_p), Sp_2n(Q_p), etc. — none of these is abelian.

(Z/3^n)* (and its inverse limit Z_3^*) is a 1-dim **abelian** p-adic Lie group (≅ Z_3 × Z/2). It's not "almost simple"; in BQ terminology it's a torus / abelian group, the trivial case for which BQ's results give no content (the unique non-atomic Haar measure is the obvious one, and there's no orbit-equidistribution question).

**Phase 3 disposition: STRUCTURALLY_BLOCKED** by the abelian-vs-semisimple category gap.

---

## Disposition D: **NO_FIT** (categorical).

- Phase 1: group and walk-support hypotheses both fail; Zariski-density vacuous in abelian setting.
- Phase 2: conclusion shape mismatch (qualitative not quantitative).
- Phase 3: p-adic extension hook doesn't apply (semisimple-required).
