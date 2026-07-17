# AFL_DISPOSITION — Accardi-Frigerio-Lewis 1982 vs Syracuse

**Date:** 2026-05-15
**Mode:** E. Self-adversarial. Honest verdict per brief's A/B/C menu.
**Companion:** `AFL_VERBATIM.md`, `AFL_SYRACUSE_IDENTIFICATION.md`, `AFL_MOMENT_PREDICTIONS.md`

---

## 0. One-sentence verdict

**Outcome B (near fit, structurally bounded) — AFL 1982 captures Syracuse's abelian-past filtration cleanly (a real structural advance over HP/AP), but FAILS at the row-(d) moment $0.108$ at distinct repeated indices for the same underlying reason HP/AP failed: AFL's embeddings $j_t : \mathcal O \to \mathcal A$ are *-homomorphism transports of a FIXED observable algebra, while Syracuse's per-step transfer operators $T_j$ are LEVEL-GRADED with accumulator coupling baked into the operator itself.**

The closest published AFL configuration that COULD match Syracuse requires (i) non-stationary, level-graded transition expectations $\{\mathcal E_j\}$, OR (ii) embedding into a "generalized AFL process" with arbitrary correlation kernels (tautological). Neither is AFL 1982 proper.

This is **structurally better than HP/AP (Outcome C)** — AFL's abelian-filtration option is a genuine match. But it's **not a clean fit (Outcome A)** — the row (d) failure is at $4 \times 10^4$ separation above noise, not fixable inside AFL's stated axioms.

---

## 1. Detailed verdict

### What AFL gets RIGHT about Syracuse (positive structural match)

1. **Abelian past filtration.** AFL admits $\mathcal A_{t]}$ abelian when $\mathcal O$ is abelian and embeddings $j_s(\mathcal O), j_t(\mathcal O)$ pairwise commute. Syracuse's $\mathbb B_j$ is abelian. **AFL is the first probe in the closure-hunt arc to natively support this**, as opposed to HP/AP which baked in non-commutative filtration via Fock or qubit-chain construction.

2. **Single state $\varphi$ + multi-time correlation kernels** with reconstruction theorem. AFL's framework explicitly takes the moment table $W_{t_1,\ldots,t_n}(a_1,\ldots,a_n) = \varphi(j_{t_1}(a_1) \cdots j_{t_n}(a_n))$ as the primary object, with consistency conditions corresponding to associativity / *-conjugation. Syracuse's moment table is the same kind of object.

3. **Quasi-conditional expectations** as the generalized "conditioning" operator. Syracuse's $\varphi = E_{\mathbb B}$ is a genuine Umegaki conditional expectation (a special case of quasi-conditional). AFL's relaxation accommodates the more general case if needed — though Syracuse doesn't require it.

4. **Row (b) cross-step second-moment vanishing** is compatible with AFL via the kernel condition $\mathcal P(a_*) = 0$ on the transition operator. Syracuse's row (b) $\approx 0$ matches AFL when this condition holds.

5. **Fubini constancy** $F(v_1, v_1') = 6.347 \times 10^{-2}$ is compatible with an **ergodic** AFL transition expectation (1-dim fixed-point space). This matches Syracuse's R77 $T_{\text{diag}}$ 1-dim invariant eigenspace at $(1, 4)$ — structurally clean.

**Net positive.** AFL is a **strictly better structural fit** than HP/AP for items 1, 2, 4, 5. Item 3 doesn't differentiate. This is real progress.

### Where AFL FAILS (the row (d) / (f) structural breakage)

1. **Row (d) $= 0.108$ at $(j_1, j_2, j_1)$ with $j_1 < j_2$.** Under any AFL identification matching items 1-5 above (abelian filtration, kernel condition on $\mathcal P$), the moment $\varphi(j_{j_1}(a_*) j_{j_2}(a_*) j_{j_1}(a_*))$ is **forced to zero**. The reason: with abelian embeddings and $\mathcal P(a_*) = 0$, the two $j_{j_1}$ factors commute past $j_{j_2}$ and contract via $j_{j_1}(a_*)^2$; the remaining $\mathcal P^{j_2 - j_1}(a_*) = 0$ by the kernel condition. So row (d) = 0 in AFL.

   Syracuse measures $0.108$, $\sim 2.5 \times 10^4$ above noise.

2. **Row (f) $= 0.609$ at $(j_1, j_2, j_1, j_2)$.** Same structural reason — the alternating pattern factors via abelian commutation and the kernel condition forces vanishing.

3. **The deeper structural reason.** AFL's *-homomorphism embeddings $j_t : \mathcal O \to \mathcal A$ transport a **fixed observable algebra** $\mathcal O$ to time $t$. The "random variable at time $t$" is the SAME random variable, just at a different time. In stationary processes, time-translation is an automorphism. In Syracuse, $T_j$ is NOT the time-translate of $T_1$: the phase factor $\chi_j(b_{[1, j-1]})$ involves the accumulator $b_{[1, j-1]} = v_1 + \cdots + v_{j-1}$, which is **the cumulative past**, not the time-index $j$ alone. There is no *-homomorphism $\mathcal O \to \mathbb A$ sending $T_1 \to T_j$ for $j > 1$.

   **AFL requires "the same random variable transported" — Syracuse has "the operator at each step depends on the cumulative past."**

4. **The adapted-process rescue fails.** Following HP/AP's idea: $\tilde X_j = F_j \cdot \hat X_j$ where $F_j \in \mathcal A_{j-1]}$ is past-measurable and $\hat X_j$ is the "pure AFL increment." Syracuse's $\tilde X_j$ has accumulator content **inside** the $\Sigma_{v,v'}$ integral, mixed with the shift content $\sigma_{-(v+v')}$. **No such factorization exists.** (Same finding as `QSC_SYRACUSE_IDENTIFICATION.md` §3.)

### Why this is Outcome B (not Outcome C)

AFL captures **half** of what Syracuse needs structurally: the abelian-filtration + single-state + multi-time-kernels skeleton. HP/AP captured **none** of this (their filtrations were forced non-commutative). The structural advance from HP/AP → AFL is genuine.

But AFL captures none of the **level-graded accumulator-coupled increment** structure that Syracuse exhibits in rows (d) and (f). This is the same gap HP/AP had, just in a different language: AFL's framework presumes the increments are time-translated copies of a fixed observable algebra, not level-graded objects with past-cumulative coupling baked in.

**Outcome B reflects this duality**: AFL is the right CATEGORY (algebraic probability spaces with filtration and conditional expectation), but the wrong INSTANCE (the embedding-of-fixed-$\mathcal O$ structure is too rigid for Syracuse's level-graded $T_j$).

---

## 2. What Syracuse needs beyond AFL

From `QSC_DISPOSITION.md` §2 (P1-P7) and the AFL analysis here:

| Property | AFL | Syracuse status |
|---|---|---|
| (P1) Abelian filtration | ✓ admits | ✓ achieves |
| (P2) Single fixed operator per step | weak (AFL embeddings are *-hom transports of fixed $\mathcal O$, not "single fixed" in Syracuse's sense) | ✓ achieves (one $\tilde X_j$ per step) |
| (P3) Cross-step second moment = 0 | ✓ admits via $\mathcal P(a_*) = 0$ | ✓ achieves |
| (P4) Distinct-index $n$-moments = 0 | ✓ admits via repeated application of $\mathcal P(a_*) = 0$ | ✓ achieves |
| (P5) Non-adjacent repeat moments NON-ZERO (e.g., $0.108$ at $(j_1, j_2, j_1)$) | ✗ FORCES zero under abelian + kernel condition | ✓ achieves $0.108$ |
| (P6) $\mathbb B$-content INSIDE the noise integral (not factored) | ✗ AFL presumes $\hat X_j$ is independent of prior accumulators; modifications go through adapted-process factorization, which Syracuse rejects | ✓ achieves (level-graded) |
| (P7) Constant Fubini inner factor | ✓ admits via ergodic $\mathcal P$ | ✓ achieves $6.347 \times 10^{-2}$ |

**AFL achieves P1, P3, P4, P7. Fails P5, P6.** Two structural feature failures (out of seven).

By contrast (per `QSC_DISPOSITION.md` §2):
- HP/AP fail P1 (filtration non-commutative), P5, P6, P7 ambiguous. **Four failures.**
- AFL fails **two**. Strictly better.

---

## 3. Mode-E gaps remaining

| Gap | Description | Effort |
|---|---|---|
| AFL-G1 | Verbatim AFL 1982 PDF (EMS Press paywalled, sandbox network blocked). Definitions in `AFL_VERBATIM.md` are reconstructed from secondary sources; verbatim line-by-line still pending. | 10 min user fetch |
| AFL-G2 | Verbatim Kümmerer 1985 (ScienceDirect paywall). Conclusion about Frigerio-Kümmerer Markov dilations is structural — iid noise factor incompatible with level-graded $T_j$ — but the verbatim form is unfetched. | 10 min user fetch (institutional) |
| AFL-G3 | "Generalized quantum Markov chains" (Accardi-Souissi-Soueidy 2020s) with non-stationary $\mathcal E_j$ — could those level-grade the way Syracuse needs? Not closed in this session. **Possible Outcome A pathway under a less restrictive AFL family.** | 1-2 day lit pull |
| AFL-G4 | The Belavkin 1989 quantum filtering framework (queued as primary recommendation per `QSC_DISPOSITION.md` §5(c)): explicitly handles classical observation filtration + non-commutative system dynamics. Brief deferred Belavkin to next probe; not closed here. | 1-2 day probe |
| AFL-G5 | The "level-graded refinement of monotone-B" hypothesis (per `QSC_DISPOSITION.md` §2(iii)) — not addressed here; queued. | 1-3 day |

---

## 4. Implications for the c = 7/45 derivation

**No change.** Per `THEOREM_C_745.md` and `D3_DERIVATION_AUDIT.md`, the leading c = 7/45 coefficient is rigorous unconditional via R75 + R76 + R77 + R64.B + HR74. AFL fit status is irrelevant to c = 7/45.

The AFL probe addresses the **higher-order structure** of Syracuse — specifically, identifying the framework that explains rows (d), (f), Fubini constancy — not the leading-order computation.

---

## 5. Recommended next step

Per the brief and `QSC_DISPOSITION.md` §5(c), the next pull is **Belavkin 1989 quantum filtering** (paywall TBD, arXiv variants may exist). Belavkin filtering is explicitly designed for the **classical observation filtration + non-commutative system algebra** setup — which is **Syracuse's exact category** (P1 + non-commutative $T_j$ separation).

The AFL probe rules out a clean AFL 1982 fit but **does not preclude Belavkin from being the right home**. AFL's "embeddings $j_t$ of fixed $\mathcal O$" rigidity is exactly what Belavkin loosens by separating observation algebra (classical, abelian) from system algebra (quantum, non-commutative) into TWO different objects, coupled by a filtering equation. Syracuse's $\mathbb B$ (abelian observation accumulator) vs $T_j$ (non-commutative system transfer operator) would map naturally to Belavkin's observation/system bipartition.

Belavkin is queued as the next probe.

---

## 6. Verdict in one paragraph

**AFL 1982 is a structural advance over HP/AP (it natively supports Syracuse's abelian past filtration $\mathbb B_j$), but it does not provide a clean fit. The decisive failure is at row (d): under the moment-table-consistent identification with abelian embeddings and the kernel condition on the transition operator $\mathcal P$, AFL FORCES $\varphi(j_{j_1}(a_*) j_{j_2}(a_*) j_{j_1}(a_*)) = 0$, while Syracuse measures $0.108$. The deeper reason is structural: AFL's $j_t : \mathcal O \to \mathcal A$ embeddings are *-homomorphism transports of a FIXED observable algebra to time $t$, but Syracuse's $T_j$ are LEVEL-GRADED operators with phase content $\chi_j(b_{[1, j-1]})$ depending on the cumulative past — NOT time-translates of a fixed operator $T_1$. No adapted-process factorization rescues this (the accumulator content is inside the noise integral, not multiplicatively factored out — same finding as the HP/AP probe). AFL captures Syracuse's items P1, P3, P4, P7 cleanly; fails P5 and P6. Outcome B — closer fit than HP/AP (Outcome C), but not closed. Belavkin queued next.**

---

## 7. Files

- This file: `C:/Collatz/AFL_DISPOSITION.md`
- Verbatim AFL: `C:/Collatz/AFL_VERBATIM.md`
- Identification: `C:/Collatz/AFL_SYRACUSE_IDENTIFICATION.md`
- Moment predictions: `C:/Collatz/AFL_MOMENT_PREDICTIONS.md`
- Project context: `QSC_DISPOSITION.md`, `H1_PRIME_DISPOSITION.md`, `D1_DISPOSITION.md`, `THEOREM_C_745.md`, `AMALG_FREENESS_SETUP.md`, `C1_TAO_RECURSION_FORM.md`
