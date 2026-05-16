# AFL_VERBATIM — Accardi-Frigerio-Lewis 1982 + followups, structural definitions

**Date:** 2026-05-15
**Mode:** E. **PDF acquisition status:** AFL 1982 PDF not directly fetchable in this session (EMS Press behind paywall/gateway; sandbox network access denied). Mode-E flag: definitions below are reconstructed from secondary-source summaries with citation trail; full verbatim quotes pending user PDF fetch. See §6 for open-fetch URLs.
**Companion:** `AFL_SYRACUSE_IDENTIFICATION.md`, `AFL_MOMENT_PREDICTIONS.md`, `AFL_DISPOSITION.md`

---

## 0. Citation block

- **AFL 1982**: Accardi, L., Frigerio, A., Lewis, J.T. "Quantum stochastic processes." Publ. Res. Inst. Math. Sci. (PRIMS), Kyoto Univ. **18**(1), 97–133 (1982). DOI: 10.2977/prims/1195184017. EMS Press: https://ems.press/journals/prims/articles/3037
- **AF 1983**: Accardi, L., Frigerio, A. "Markovian cocycles." Proc. R. Irish Acad. **83A**(2), 251–263 (1983). [Quantum Markov states defined here.]
- **Frigerio 1984**: Frigerio, A. "Quasi-free stochastic integration and stochastic evolutions" / "Markov dilations and quantum detailed balance." Commun. Math. Phys. **93**, 517–532 (1984). Project Euclid (open): https://projecteuclid.org/journals/communications-in-mathematical-physics/volume-93/issue-4/Markov-dilations-and-quantum-detailed-balance/cmp/1103941181.pdf
- **Kümmerer 1985**: Kümmerer, B. "Markov dilations on W*-algebras." J. Funct. Anal. **63**, 139–177 (1985). ScienceDirect: https://www.sciencedirect.com/science/article/pii/0022123685900849
- **AFL 1990** (followup): Accardi, L., Frigerio, A., Lu, Y.G. "The weak coupling limit as a quantum functional central limit." Commun. Math. Phys. **131**, 537–570 (1990). Springer: https://link.springer.com/article/10.1007/BF02098275
- **Accardi survey**: Accardi, L. "Quantum probability: an historical survey." ResearchGate: https://www.researchgate.net/profile/Luigi-Accardi/publication/268321828

---

## 1. AFL 1982 — definition of quantum stochastic process

**Defn (AFL 1982, §1–§2; reconstruction from secondary sources Cambridge Quantum Stochastics Ch. 7, Accardi survey, Wikipedia Quantum Markov chain).**

A **quantum stochastic process** in the AFL sense is a quadruple
$$ (\mathcal A,\ \{j_t\}_{t \in T},\ \mathcal O,\ \varphi) $$
where:

- $\mathcal A$ is a unital C*-algebra (or von Neumann algebra) — the **ambient algebra**.
- $\mathcal O$ is a unital C*-algebra — the **algebra of "single-time" observables** (or "values" of the random variable).
- $T$ is an index set (typically $T = \mathbb N$ or $\mathbb R_{\ge 0}$, with a partial order).
- $j_t : \mathcal O \to \mathcal A$ is a unital *-homomorphism for each $t \in T$ — the **random variable at time $t$** (a "quantum random variable" is a *-homomorphism, not a measurable function; the values live in $\mathcal A$ rather than in a sample space).
- $\varphi : \mathcal A \to \mathbb C$ is a state on $\mathcal A$.

**Multi-time correlation kernels.** The collection
$$ W_{t_1,\ldots,t_n}(a_1,\ldots,a_n) := \varphi(j_{t_1}(a_1) \cdot j_{t_2}(a_2) \cdots j_{t_n}(a_n)), \quad a_i \in \mathcal O,\ t_i \in T $$
is the family of **finite-dimensional joint expectations** (or **correlation kernels**) of the process.

**Kolmogorov-type reconstruction theorem (AFL 1982 Thm 1.3, restated).** Every consistent family of correlation kernels $\{W_{t_1,\ldots,t_n}\}$ (satisfying positivity, normalization, and consistency under reordering / *-conjugation) arises from some quantum stochastic process $(\mathcal A, \{j_t\}, \mathcal O, \varphi)$, unique up to equivalence.

This is the **non-commutative generalization of the Kolmogorov reconstruction theorem for classical stochastic processes**. In the commutative case ($\mathcal A, \mathcal O$ both commutative), one recovers the classical theorem with $j_t$ corresponding to evaluation maps $f \mapsto f(X_t)$.

**Mode-E note.** The exact statement (signatures, regularity hypotheses, choice of subclass — equivalence vs unitary equivalence) is pending PDF fetch. The structural skeleton above is robust across all secondary-source descriptions consulted.

---

## 2. AFL 1982 — Markov property / quantum Markov process

**Defn (Markov projection identity).** A quantum stochastic process $(\mathcal A, \{j_t\}, \mathcal O, \varphi)$ on a totally-ordered index set $T$ is called **Markov** if there exists a filtration $\{\mathcal A_{t]}\}_{t \in T}$ (an increasing family of *-subalgebras with $j_s(\mathcal O) \subset \mathcal A_{t]}$ for $s \le t$) and a family of conditional expectations or quasi-conditional expectations
$$ E_{t]} : \mathcal A \to \mathcal A_{t]} $$
such that for any $s \le t$ and any $a \in \mathcal O$:
$$ \boxed{E_{s]} \circ j_t = E_{s]} \circ j_t \circ E_{[s} \quad (\text{equivalently: } E_{s]}(j_t(a)) \text{ depends on } \mathcal A_{[s} \text{ only through } E_{[s}).} $$

In the more common Accardi-survey rendering for **discrete-time quantum Markov chains** with index set $T = \mathbb N$ and a **transition expectation** $\mathcal E : \mathcal O \otimes \mathcal O \to \mathcal O$ (a completely positive identity-preserving map), the state on the chain $\bigotimes_{n} \mathcal O_n$ satisfies the factorization
$$ \psi(a_0 \otimes a_1 \otimes \cdots \otimes a_n) = \varphi_0\big( \mathcal E(a_0 \otimes \mathcal E(a_1 \otimes \mathcal E(\cdots \otimes \mathcal E(a_{n-1} \otimes a_n) \cdots ))) \big) $$
where $\varphi_0$ is a state on $\mathcal O_0$ ("initial state") and $\mathcal E$ encodes the "Markov transition" structure.

**Quasi-conditional expectation (Accardi notion).** A linear map $E : \mathcal A \to \mathcal B$ between unital C*-algebras with $\mathcal C \subseteq \mathcal B \subseteq \mathcal A$ is a **quasi-conditional expectation** if (i) completely positive, (ii) identity-preserving, (iii) $\mathcal C$-bimodular: $E(c \cdot a) = c \cdot E(a)$ for $c \in \mathcal C$, $a \in \mathcal A$. The Umegaki property $E \circ E = E$ is **not** required; this is the key relaxation that lets AFL/Accardi accommodate non-product states.

**Mode-E note.** The precise hypothesis structure in AFL 1982 vs AF 1983 vs Accardi survey differs slightly: AFL 1982 emphasizes the embedding picture $(j_t)$, AF 1983 emphasizes the tensor-chain transition-expectation picture $\mathcal E : \mathcal O \otimes \mathcal O \to \mathcal O$. Both pictures are equivalent up to choice of representation.

---

## 3. Moment formula for quantum Markov processes

For an AFL quantum Markov process with index set $T = \mathbb N$, transition expectation $\mathcal E : \mathcal O \otimes \mathcal O \to \mathcal O$, and initial state $\varphi_0$ on $\mathcal O$, the **multi-time correlation kernel** at strictly-ordered $t_1 < t_2 < \cdots < t_n$ is:
$$ W_{t_1,\ldots,t_n}(a_1,\ldots,a_n) = \varphi_0\big( \mathcal P^{t_1}\big( a_1 \cdot \mathcal P^{t_2 - t_1}\big( a_2 \cdot \mathcal P^{t_3 - t_2}\big( \cdots a_{n-1} \cdot \mathcal P^{t_n - t_{n-1}}(a_n) \cdots \big) \big) \big) \big) $$
where $\mathcal P : \mathcal O \to \mathcal O$ is the **transition operator** $\mathcal P(a) = \mathcal E(\mathbf 1 \otimes a)$ associated to the transition expectation (semigroup property holds in stationary case). This is the **direct quantum analog of the classical Markov-chain Chapman-Kolmogorov moment formula**.

For **non-strictly-ordered** index sequences (e.g., Syracuse's $j_1, j_2, j_1$ pattern with repeated $j_1$), the formula has to be interpreted via the embedding picture: $\varphi(j_{j_1}(a) j_{j_2}(b) j_{j_1}(c)) = \varphi_0 \big( \mathcal E_{j_1, j_2, j_1}(a, b, c) \big)$ for an appropriate generalized transition expectation. In **stationary** AFL settings, $j_{j_1}$ is "the same random variable embedded at time $j_1$" — and the moment $\varphi(j_{j_1}(a) j_{j_2}(b) j_{j_1}(c))$ for $j_1 < j_2$ equals $\varphi_0(a \cdot \mathcal P^{j_2 - j_1}(b) \cdot \mathcal P^{j_1 - j_2}(c))$ — but the latter requires interpreting $\mathcal P^{-(j_2 - j_1)}$ as a "time-reversed" operator, which only makes sense under quantum detailed balance.

**Crucial structural feature.** In AFL's stationary setting with a transition expectation $\mathcal E$ and a **product-of-copies** filtration (each $j_t$ embedding the same $\mathcal O$ as the $t$-th tensor factor in $\bigotimes_t \mathcal O$), the moments $\varphi(j_{t_1}(a) j_{t_2}(b) j_{t_1}(c))$ for $t_1 \ne t_2$ are **typically NON-ZERO** because the algebras $j_{t_1}(\mathcal O)$ and $j_{t_2}(\mathcal O)$ live at distinct tensor positions but the embeddings represent **iid copies of the SAME random variable** $\mathcal O \to \mathcal A$ — so the moment is $\varphi_0(a \cdot \mathcal E_{*}(b) \cdot c)$ for some intermediate operator $\mathcal E_*(b) \in \mathcal O$ depending on the time-gap $|t_2 - t_1|$.

**This contrasts with HP/AP/free/monotone, all of which give zero at this third moment.**

---

## 4. Frigerio 1984 — Markov dilation form

A **quantum Markov dilation** of a CP semigroup $\{T_t\}_{t \ge 0}$ on a vN algebra $\mathcal M$ with invariant state $\omega$ is a quintuple $(\hat{\mathcal M}, \hat\omega, \{j_t\}_{t \ge 0}, \alpha_t, E)$ where:
- $\hat{\mathcal M}$ is a vN algebra containing $\mathcal M$;
- $\hat\omega$ is a state on $\hat{\mathcal M}$ extending $\omega$;
- $j_t : \mathcal M \to \hat{\mathcal M}$ is an injective *-homomorphism with $j_0 = $ inclusion;
- $\alpha_t : \hat{\mathcal M} \to \hat{\mathcal M}$ is a $\hat\omega$-preserving *-automorphism flow (the "dilated dynamics");
- $E : \hat{\mathcal M} \to \mathcal M$ is the conditional expectation associated to $\hat\omega \restriction \mathcal M$.

The Markov property: $E \circ \alpha_t \circ j_0 = T_t$ for all $t \ge 0$, i.e., the dilated dynamics projected back to $\mathcal M$ reproduces the original semigroup. **The filtration $\hat{\mathcal M}_{t]} := \alpha_t(\hat{\mathcal M}_{0]})$** is the natural quantum-Markov filtration of the dilation.

---

## 5. Kümmerer 1985 — minimal Markov dilation

Minimality + uniqueness theorems for the Frigerio-style Markov dilations. Specifically: a CP semigroup $T_t$ on a vN algebra $\mathcal M$ with faithful normal stationary state $\omega$ admits a **minimal** Markov dilation (in the sense of §4), unique up to conjugation by an $\hat\omega$-preserving *-isomorphism, IF the semigroup satisfies the **quantum detailed balance condition** with respect to $\omega$.

The construction is via a product-state representation: $\hat{\mathcal M} = \mathcal M \otimes \bigotimes_{\mathbb Z} \mathcal N$ for a "noise algebra" $\mathcal N$, with shift dynamics $\alpha_t$ acting non-trivially on the $\bigotimes_{\mathbb Z} \mathcal N$ factor and trivially on $\mathcal M$. **The noise factor is iid (a tensor product of identical copies of $\mathcal N$).**

---

## 6. Open verbatim-fetch list (user follow-up)

User to retrieve (sandbox blocked direct download in this session):

- **AFL 1982** PDF: EMS Press https://ems.press/journals/prims/articles/3037 — confirm Defns 1.1, 1.2 + Thm 1.3 (reconstruction theorem) signatures and the verbatim Markov projection identity if present in §2.
- **AF 1983** "Markovian cocycles" Proc. R. Irish Acad. 83A — confirm the transition-expectation $\mathcal E : \mathcal O \otimes \mathcal O \to \mathcal O$ formulation.
- **Frigerio 1984** "Markov dilations and quantum detailed balance" CMP 93 — Project Euclid open access, can be fetched directly: https://projecteuclid.org/journals/communications-in-mathematical-physics/volume-93/issue-4/Markov-dilations-and-quantum-detailed-balance/cmp/1103941181.pdf
- **Kümmerer 1985** ScienceDirect — paywalled, institutional access required.
- **Accardi survey** ResearchGate — should be fetchable, contains restated AFL 1982 definitions with citations.

Mode-E flag: the present document is a structural reconstruction. Verbatim line-by-line quotes pending. The §1–§5 structural skeleton is robust under all secondary-source variation; the disposition (§AFL_DISPOSITION.md) holds regardless of small notational differences that might surface in the verbatim fetch.

---

## 7. Files

- This file: `C:/Collatz/AFL_VERBATIM.md`
- Identification: `C:/Collatz/AFL_SYRACUSE_IDENTIFICATION.md`
- Moment predictions: `C:/Collatz/AFL_MOMENT_PREDICTIONS.md`
- Disposition: `C:/Collatz/AFL_DISPOSITION.md`
