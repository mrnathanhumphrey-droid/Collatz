# AFL_MOMENT_PREDICTIONS — AFL 1982 moment predictions vs Syracuse numerics

**Date:** 2026-05-15
**Mode:** E. Numerical targets from `D1_DISPOSITION.md`; AFL predictions derived from `AFL_VERBATIM.md` §3 + `AFL_SYRACUSE_IDENTIFICATION.md` §3-§6.

---

## 0. The targets (Syracuse Reading B, sum_entries, n=4 alternating)

| Row | Moment | Syracuse target |
|---|---|---|
| (b) | $\varphi(\tilde X_{j_1} \cdot \tilde X_{j_2})$, $j_1 \ne j_2$ | $\approx 0$ (noise floor $1.08 \times 10^{-7}$) |
| (d) | $\varphi(\tilde X_{j_1} \cdot \tilde X_{j_2} \cdot \tilde X_{j_1})$ | **$0.108$** (4 orders above noise) |
| (f) | $\varphi(\tilde X_{j_1} \cdot \tilde X_{j_2} \cdot \tilde X_{j_1} \cdot \tilde X_{j_2})$ | **$0.609$** |
| Fubini | $F(v_1, v_1') = \mathbb E_{(v_2)}[\tilde X_{j_2} \cdot \tilde X_{j_1} \cdot \tilde X_{j_2}]$ | **constant $6.347 \times 10^{-2}$** across 12 grid points |

---

## 1. AFL moment formula for the ideal embedding picture

If Syracuse genuinely fit AFL with embeddings $j_j : \mathcal O \to \mathcal A$ and a transition expectation $\mathcal E$ (stationary case, see `AFL_VERBATIM.md` §3), the moments would be:
$$ \varphi(j_{j_1}(a) j_{j_2}(b) j_{j_3}(c) \ldots) = \varphi_0\big( a \cdot \mathcal P^{j_2 - j_1}(b \cdot \mathcal P^{j_3 - j_2}(c \cdot \ldots)) \big) $$
for $j_1 < j_2 < j_3 < \ldots$, where $\mathcal P(a) = \mathcal E(\mathbf 1 \otimes a)$ is the transition operator.

For the **alternating pattern** $j_1, j_2, j_1$ with $j_1 < j_2$ (repeated $j_1$), the formula requires interpreting the second appearance of $j_1$ at a time-position **after** $j_2$. In AFL's stationary embedding picture with time-ordered products, this is not directly defined — one must use the operator-product evaluation in $\mathcal A$, $\varphi(j_{j_1}(a) \cdot j_{j_2}(b) \cdot j_{j_1}(c))$, which depends on the commutation relations between $j_{j_1}(\mathcal O)$ and $j_{j_2}(\mathcal O)$.

### Row (b): two distinct times

AFL prediction (centered increments): if $\hat X_j := j_j(a_*)$ for $a_* \in \mathcal O$ centered (i.e., $\varphi_0(a_*) = 0$), then
$$ \varphi(\hat X_{j_1} \hat X_{j_2}) = \varphi_0(a_* \cdot \mathcal P^{j_2 - j_1}(a_*)). $$
This is the standard **autocovariance** $C(j_2 - j_1)$ of the AFL Markov chain. For stationary processes with mean-zero, $C(k) \ne 0$ generically.

✗ **AFL predicts $\varphi(\hat X_{j_1} \hat X_{j_2}) \ne 0$ generically.** Syracuse measures $\approx 0$ (Row b).

This is the **first numerical mismatch**: AFL's stationary Markov chain has non-zero autocovariance at distinct times, while Syracuse has cross-step second moment vanishing.

#### Caveat

If the AFL transition expectation $\mathcal P$ satisfies $\mathcal P(a_*) = 0$ (i.e., $a_*$ is in the **null space** of $\mathcal P$ on $\mathcal O$), then $C(k) = 0$ for $k \ge 1$. This is the **"orthogonal increment" / martingale condition** in AFL terms.

Under this condition: ✓ Row (b) = 0. AFL accommodates Syracuse's row (b) **conditionally on $a_*$ being in the kernel of $\mathcal P$**.

This is structurally feasible — AFL allows degenerate transition expectations. So row (b) is accommodated.

### Row (d): repeated $j_1$ with intervening $j_2$

AFL prediction for $\varphi(j_{j_1}(a) j_{j_2}(b) j_{j_1}(c))$ with $j_1 < j_2$ in the **stationary tensor-chain picture**:
The third factor $j_{j_1}(c)$ at "time $j_1$" — but we've already passed time $j_2$ in the product. This requires choosing a convention.

**Convention A (operator-product, no time-reorder).** Just evaluate the product of three operators in $\mathcal A$. The operators $j_{j_1}(a)$ and $j_{j_1}(c)$ live in the same algebra $j_{j_1}(\mathcal O)$ at time $j_1$, but the middle operator $j_{j_2}(b)$ at time $j_2$ may not commute with them in general.

If $j_{j_1}(\mathcal O)$ and $j_{j_2}(\mathcal O)$ commute (e.g., abelian $\mathcal O$ with classical chain), then:
$$ \varphi(j_{j_1}(a) j_{j_2}(b) j_{j_1}(c)) = \varphi(j_{j_1}(a c) \cdot j_{j_2}(b)) = \varphi_0(a c \cdot \mathcal P^{j_2 - j_1}(b)) $$
where we used commutativity to bring the two $j_{j_1}$ factors together as $j_{j_1}(a c)$.

If additionally $\mathcal P(b) = 0$ (consistent with the row-(b) condition), then $\mathcal P^{j_2 - j_1}(b) = 0$ for $j_2 > j_1$, and **the row (d) moment vanishes**: $\varphi_0(a c \cdot 0) = 0$.

✗ **Under the row-(b)-consistent kernel condition $\mathcal P(b) = 0$, AFL predicts row (d) = 0.** Syracuse measures $0.108$.

This is the **decisive numerical mismatch**. The kernel condition that explains row (b) FORCES row (d) = 0 in AFL — but Syracuse has row (d) non-zero.

**Convention B (operator-product, non-commuting $j_{j_1}(\mathcal O), j_{j_2}(\mathcal O)$).** Drop commutativity. Then the moment $\varphi(j_{j_1}(a) j_{j_2}(b) j_{j_1}(c))$ does not factorize. In general it's non-zero, and the actual value depends on how the $j_{j_1}$ and $j_{j_2}$ embeddings sit relative to each other in $\mathcal A$.

This **could** match Syracuse's $0.108$ — but only by **making the filtration non-abelian** (since the $j_{j_1}(\mathcal O), j_{j_2}(\mathcal O)$ generate the filtration). This contradicts the structural match in `AFL_SYRACUSE_IDENTIFICATION.md` §1 that AFL's value lay precisely in supporting an abelian filtration.

✗ **No AFL identification simultaneously gives (i) abelian filtration matching Syracuse's $\mathbb B_j$, and (ii) row (d) = 0.108.**

### Row (f): alternating $j_1, j_2, j_1, j_2$

Same analysis as row (d), one step further. Under Convention A (commuting embeddings, row-(b) kernel condition), row (f) factors as
$$ \varphi(j_{j_1}(a)^2 \cdot j_{j_2}(b)^2 \cdot (\text{intermediate transitions})) $$
which under $\mathcal P(b) = 0$ collapses to $\varphi_0(a^2) \cdot \varphi_0(b^2) \cdot 0 = 0$ or similar — generically **zero** under the kernel condition.

Syracuse measures $0.609$. ✗

### Fubini inner factor

Syracuse: $F(v_1, v_1') = 6.347 \times 10^{-2}$, **CONSTANT** across all 12 grid points.

Under an AFL identification with commuting embeddings $j_{j_1}(\mathcal O), j_{j_2}(\mathcal O)$, the inner factor reduces to $\mathbb E_{(v_2)}[\mathcal P^{j_2 - j_1}(b) \cdot j_{j_1}(c)]$, which generically **depends** on the $j_1$-fiber data $(v_1, v_1')$.

For $F$ to be **constant** in $(v_1, v_1')$, AFL would need a specific structural property: the inner integration $\mathbb E_{(v_2)}[\mathcal P^{j_2 - j_1}(b)]$ produces a scalar in $\mathcal O$'s center (i.e., a multiple of identity in $\mathcal O$), so the outer multiplication by $j_{j_1}(c)$ contributes only via $\varphi_0(c)$, independent of $(v_1, v_1')$.

This **could** be encoded in AFL by a specific choice of transition expectation $\mathcal E$. The structural feature would be: $\mathcal E$ has 1-dimensional fixed-point space spanned by identity. **Stationary AFL transition expectations** with this property exist — they correspond to "ergodic" classical Markov chains with unique invariant measure.

✓ **The Fubini constancy is structurally compatible with an AFL ergodic transition expectation.** This is consistent with Syracuse's $T_{\text{diag}}$ having a 1-d eigenspace on $(1, 4)$ at eigenvalue 1 (R77).

### Composite

| Row | Syracuse | AFL prediction | Status |
|---|---|---|---|
| (b) $\varphi(\tilde X_{j_1} \tilde X_{j_2})$ | $\approx 0$ | $C(j_2-j_1) = \varphi_0(a_* \mathcal P^{j_2-j_1}(a_*))$, generically $\ne 0$; $= 0$ if $\mathcal P(a_*) = 0$ | ✓ if kernel condition |
| (d) $\varphi(\tilde X_{j_1} \tilde X_{j_2} \tilde X_{j_1})$ | $0.108$ | $= 0$ under (i) abelian embedding + (ii) row-b kernel condition; **non-zero only under non-abelian embedding** | ✗ FAILS under abelian filtration |
| (f) $\varphi(\tilde X_{j_1} \tilde X_{j_2} \tilde X_{j_1} \tilde X_{j_2})$ | $0.609$ | $= 0$ under same conditions | ✗ FAILS |
| Fubini constant | $6.347 \times 10^{-2}$ | constant if $\mathcal P$ ergodic with 1-d fixed-point | ✓ compatible |

---

## 2. The structural reason AFL fails rows (d) and (f)

Both AFL and Syracuse have an abelian-side filtration. **In AFL, when the filtration is abelian (i.e., the embeddings $j_t(\mathcal O)$ pairwise commute), the embeddings $j_t$ effectively factor through a classical (commutative) Markov chain.** All multi-time moments reduce to classical Markov-chain moments via Chapman-Kolmogorov, with the **time-local kernel condition** $\mathcal P(a_*) = 0$ propagating through to **kill all moments containing more than one transition step**.

Syracuse violates this. **Syracuse's $\tilde X_j$ has non-zero higher moments at non-adjacent repeated indices, even though row (b) vanishes.** This is a fundamental **non-Markovian** property of $\tilde X_j$, in the classical sense:
- Row (b) = 0 says $\mathbb E[\tilde X_{j_1} \tilde X_{j_2} | \mathbb B_{j_1}] = 0$ for $j_1 \ne j_2$ (martingale-like)
- Row (d) = 0.108 says $\mathbb E[\tilde X_{j_1} \tilde X_{j_2} \tilde X_{j_1} | \mathbb B] \ne 0$ (third-order coupling)

In a Markov chain, the third-order coupling at $(j_1, j_2, j_1)$ would be ZERO if the second-order coupling at $(j_1, j_2)$ is zero — because Markov-chain moments propagate through the transition kernel multiplicatively, and a kernel that annihilates the test function at one step annihilates it at all subsequent steps.

**Syracuse's $\tilde X_j$ is not the embedding of a classical Markov chain** — it has a third-order coupling that survives despite second-order vanishing. This is **fundamentally inconsistent with the AFL framework when the filtration is abelian**.

### What could rescue AFL?

(i) **Non-abelian filtration.** Allow $j_{j_1}(\mathcal O)$ and $j_{j_2}(\mathcal O)$ to not commute. Then row (d) can be non-zero. But this gives up Syracuse's match with the abelian $\mathbb B_j$. ✗

(ii) **Non-stationary transition expectations.** Allow $\mathcal E$ to depend on time, with level-graded structure. Then rows (b) and (d) decouple. **This is the Accardi-Frigerio 1983 "generalized quantum Markov chains" generalization**, and it might be the right home for Syracuse — but it requires explicit construction of the level-graded $\mathcal E_j$, which would not be a standard AFL Markov process.

(iii) **Non-Markov AFL process.** AFL's quantum-stochastic-process definition (`AFL_VERBATIM.md` §1) does **not** require the Markov property. A general AFL process with arbitrary correlation kernels can match any moment pattern — but at the cost of carrying no predictive structure (it's just a relabeling of Syracuse's moment table).

### Verdict on rescue paths

(i) violates the structural match Syracuse had with AFL's abelian-filtration option.
(ii) is a **generalization of AFL**, not AFL proper.
(iii) is a **tautology** — any process is an AFL "quantum stochastic process" in this loose sense.

**No rescue path keeps AFL's predictive content while matching Syracuse's row (d) and (f).**

---

## 3. The "operator-valued" / dilation upgrade

The next-up framework — Frigerio 1984 "Markov dilations" — would let $\mathcal A = \mathcal M \otimes \bigotimes_{\mathbb Z} \mathcal N$ with shift dynamics. The "noise factor" $\bigotimes_{\mathbb Z} \mathcal N$ is then **iid**, with each $\mathcal N_k$ representing the noise contribution at step $k$.

For Syracuse, this requires $T_j$ to be (e.g.) the conjugation $\alpha_j \circ T_1 \circ \alpha_{-j}$ of a fixed $T_1$ by the shift $\alpha_j$. But again, $T_j$ is not the shift-translate of $T_1$ — the phase factor $\chi_j(b_{[1,j-1]})$ depends on **the cumulative shift accumulator**, not just on the shift index $j$.

✗ **Frigerio-Kümmerer Markov dilations don't capture Syracuse either**, for the same iid-noise / level-graded mismatch.

---

## 4. Summary

| Aspect | AFL accommodates | AFL fails |
|---|---|---|
| Abelian past filtration $\mathbb B_j$ | ✓ | |
| Single state $\varphi$, multi-time correlation kernels | ✓ (definition admits any kernel) | |
| Embeddings $j_t : \mathcal O \to \mathcal A$ as transports of a fixed $\mathcal O$ | | ✗ (Syracuse $T_j$ are not transports of a fixed algebra) |
| Row (b) cross-step second moment $= 0$ | ✓ (kernel condition on $\mathcal P$) | |
| Row (d) third moment $= 0.108$ at $(j_1, j_2, j_1)$ | | ✗ (kernel condition forces row (d) $= 0$) |
| Row (f) $= 0.609$ at $(j_1, j_2, j_1, j_2)$ | | ✗ same |
| Fubini constancy | ✓ (ergodic $\mathcal P$) | |
| Adapted operator $\tilde X_j = F_j \hat X_j$ rescue | | ✗ (Syracuse $\mathbb B$-content is inside the noise) |
| Markov dilation $\mathcal M \otimes \bigotimes \mathcal N$ rescue | | ✗ (iid noise vs level-graded $T_j$) |

**Net.** AFL is a structural advance over HP/AP (it admits abelian filtration), but fails Syracuse's distinctive row (d) at the same place HP/AP failed. The underlying issue is the same: **frameworks that posit "the same random variable transported by *-homomorphism at each time" cannot reproduce Syracuse's level-graded accumulator-coupled $\tilde X_j$.**

---

## 5. Mode-E gaps

- Verbatim AFL 1982 moment formula (e.g., Thm 3.x for general $n$-time products) — pending PDF fetch.
- The exact AF 1983 transition-expectation formulation and whether it admits **time-dependent (non-stationary) $\mathcal E_j$** with the required level-graded structure — not closed in this session.
- Whether the "generalized quantum Markov chain" of Accardi-Souissi-Soueidy 2020s admits Syracuse — flagged for future check; not closed here.

---

## 6. Files

- This file: `C:/Collatz/AFL_MOMENT_PREDICTIONS.md`
- Verbatim AFL: `C:/Collatz/AFL_VERBATIM.md`
- Identification: `C:/Collatz/AFL_SYRACUSE_IDENTIFICATION.md`
- Verdict: `C:/Collatz/AFL_DISPOSITION.md`
