# BELAVKIN_MOMENT_PREDICTIONS — Belavkin predictions vs Syracuse rows (b)/(d)/(f) + Fubini

**Date:** 2026-05-15
**Mode:** E. Numerical targets from `D1_DISPOSITION.md`; Belavkin predictions derived from `BELAVKIN_VERBATIM.md` §1.3–§2 + `BELAVKIN_SYRACUSE_IDENTIFICATION.md` §2–§3.

---

## 0. Targets (Syracuse Reading B, sum_entries, n=4 alternating)

| Row | Moment | Syracuse target |
|---|---|---|
| (b) | $\varphi(\tilde X_{j_1} \cdot \tilde X_{j_2})$, $j_1 \ne j_2$ | $\approx 0$ (noise floor $1.08 \times 10^{-7}$) |
| (d) | $\varphi(\tilde X_{j_1} \cdot \tilde X_{j_2} \cdot \tilde X_{j_1})$ | **$0.108$** |
| (f) | $\varphi(\tilde X_{j_1} \cdot \tilde X_{j_2} \cdot \tilde X_{j_1} \cdot \tilde X_{j_2})$ | **$0.609$** |
| Fubini | $F(v_1, v_1') = \mathbb E_{(v_2)}[\tilde X_{j_2} \cdot \tilde X_{j_1} \cdot \tilde X_{j_2}]$ | **constant $6.347 \times 10^{-2}$** across 12 grid points |

---

## 1. Belavkin's multi-time expectation formula applied to Syracuse

Under the identification of `BELAVKIN_SYRACUSE_IDENTIFICATION.md`:
- System Hilbert space: $\mathcal H_n = L^2((\mathbb Z/3^n)^*, \pi_n)$
- Observation algebra: $\mathbb B_j = vN(\{M_{b_{[1, k]}} : k \le j\})$ (abelian, coarsened to running sums)
- Kraus operators: $M_v^{(j, b_{[1, j-1]})} f(\xi) = 2^{-v/2} \cdot e^{-2\pi i \cdot 3^{2j-2} \cdot 2^{-b_{[1, j-1]} - v} \cdot (\text{phase}) / 3^n} \cdot f(\xi \cdot 2^{-v} \mod 3^n)$
- Centered system observable: $\tilde X_j = T_j - \mathbb E_{\mathbb B}[T_j]$

The Belavkin multi-time expectation for system observables along the trajectory is:
$$ \varphi(\tilde X_{j_1} \cdots \tilde X_{j_n}) = \sum_{v_{1:N}} \mu_{\mathrm{Geom}}^{\otimes N}(v_{1:N}) \cdot \mathrm{Tr}_{\mathcal H_n}\left( \tilde X_{j_n} M_{v_N}^{(N)} \cdots \tilde X_{j_1} M_{v_{j_1}}^{(j_1)} \cdots M_{v_1}^{(1)} \rho_0 (M_{v_1}^{(1)})^* \cdots (M_{v_{j_1}}^{(j_1)})^* \tilde X_{j_1} \cdots (M_{v_N}^{(N)})^* \right) $$

with $N \ge \max(j_1, \ldots, j_n)$ and the sum running over all Geom(2) outcome sequences $v_{1:N}$.

---

## 2. Row (b) — distinct-index second moment

**Syracuse:** $\varphi(\tilde X_{j_1} \tilde X_{j_2}) \approx 0$ for $j_1 \ne j_2$ (noise floor $1.08 \times 10^{-7}$).

**Belavkin prediction.** Take WLOG $j_1 < j_2$.
$$ \varphi(\tilde X_{j_1} \tilde X_{j_2}) = \sum_{v_{1:j_2}} \mu_{\mathrm{Geom}}(v_{1:j_2}) \, \mathrm{Tr}\left( \tilde X_{j_2} M_{v_{j_2}}^{(j_2)} \cdots \tilde X_{j_1} M_{v_{j_1}}^{(j_1)} \cdots \rho_0 \cdots (M_{v_{j_1}}^{(j_1)})^* \tilde X_{j_1}^* \cdots \right) $$

Centering: $\tilde X_j = T_j - \mathbb E_{\mathbb B}[T_j]$, so $\varphi(\tilde X_j) = 0$. The fact that **Syracuse measures $\varphi(\tilde X_{j_1} \tilde X_{j_2}) \approx 0$** is equivalent to the **martingale-difference property**:
$$ \mathbb E[\tilde X_{j_2} | \mathbb B_{j_2 - 1}] = 0 . $$

In the Belavkin Kraus form, this property follows from the centering: $\tilde X_j$ projects out the $\mathbb B$-measurable component of $T_j$. The conditional expectation $\mathbb E[\tilde X_j | \mathbb B_{j-1}]$ is, by definition of $\tilde X_j$, equal to zero (the $\mathbb B$-measurable part has been subtracted).

**By the tower property of conditional expectations:**
$$ \varphi(\tilde X_{j_1} \tilde X_{j_2}) = \varphi\left( \tilde X_{j_1} \cdot \mathbb E[\tilde X_{j_2} | \mathbb B_{j_2 - 1}] \right) = \varphi(\tilde X_{j_1} \cdot 0) = 0 . $$

✓ **Belavkin predicts row (b) = 0.** Matches Syracuse (and consistent with the "centered system observable on a Belavkin filter" being a $\mathbb B$-martingale-difference sequence.)

### Subtlety: the strict tower-property check

The tower property requires that $\tilde X_{j_1}$ commutes "enough" with $\tilde X_{j_2}$'s conditional expectation. In the Belavkin framework with the **non-demolition condition** between system algebra and observation algebra (verified in `BELAVKIN_SYRACUSE_IDENTIFICATION.md` §2.2), the tower property holds: $\tilde X_{j_1}$ as a system observable at step $j_1$ commutes with the observation algebra $\mathbb B_{j_2 - 1}$ which includes $\mathbb B_{j_1}$. So:
$$ \varphi(\tilde X_{j_1} \tilde X_{j_2}) = \varphi(\tilde X_{j_1} \mathbb E[\tilde X_{j_2} | \mathbb B_{j_2 - 1}]) = 0 . $$
✓ Rigorous under non-demolition.

---

## 3. Row (d) — repeated-index third moment

**Syracuse:** $\varphi(\tilde X_{j_1} \tilde X_{j_2} \tilde X_{j_1}) = 0.108$ for $j_1 < j_2$.

**Belavkin prediction.** This is the **decisive test**. The structural question: does Belavkin's framework predict a non-zero value for the repeated-index third moment $(j_1, j_2, j_1)$?

The Heisenberg-picture multi-time formula gives:
$$ \varphi(\tilde X_{j_1} \tilde X_{j_2} \tilde X_{j_1}) = \sum_{v_{1:j_2}} \mu_{\mathrm{Geom}}(v_{1:j_2}) \, \mathrm{Tr}\left( \tilde X_{j_1}^{\text{op}} M_{v_{j_2}}^{(j_2)} \cdots \tilde X_{j_2}^{\text{op}} \cdots \tilde X_{j_1}^{\text{op}} \cdots M_{v_1}^{(1)} \rho_0 (M_{v_1}^{(1)})^* \cdots (M_{v_{j_2}}^{(j_2)})^* \right) $$

**Why this is generically non-zero.** Unlike AFL (where the third moment vanished under abelian filtration via the kernel-condition contraction), Belavkin's framework has **no kernel condition** on the Kraus operators. The operators $M^{(j_2)}$ at step $j_2$ are arbitrary unitary-conjugation channels on $\mathcal H_S$, NOT shift-or-translate operators of a fixed object. They do NOT in general preserve the eigenspace structure of $\tilde X_{j_1}$.

Concretely, the operator chain
$$ \tilde X_{j_1}^{\text{op}} \cdot M^{(j_2)} \cdot \tilde X_{j_1}^{\text{op}} $$
sandwiches $\tilde X_{j_1}^{\text{op}}$ around the step-$j_2$ Kraus operator. Even though $\varphi(\tilde X_{j_1}) = 0$ (centering), the **double application** of $\tilde X_{j_1}^{\text{op}}$ with $M^{(j_2)}$ between them does NOT contract to zero — it produces a non-trivial system operator on $\mathcal H_S$, whose trace against the trajectory state is generically non-zero.

**Heuristic prediction.**
$$ \varphi(\tilde X_{j_1} \tilde X_{j_2} \tilde X_{j_1}) \sim \mathbb E\left[ \tilde X_{j_1}^2 \cdot M^{(j_2)} \right] $$
where the operator $M^{(j_2)}$ at the intervening step encodes the Geom(2)-averaged transfer over the gap $j_2 - j_1$. This is the **inner Fubini factor** appearing in Syracuse's row (d) structure.

**Belavkin specifically PREDICTS non-vanishing row (d)** because:
1. The two appearances of $\tilde X_{j_1}^{\text{op}}$ at different positions in the Heisenberg product are NOT a "kernel-condition contraction" — they are an actual operator product $\tilde X_{j_1}^{\text{op}} \cdot M^{(j_2)} \cdot \tilde X_{j_1}^{\text{op}}$.
2. Even though row (b) vanishes (first moment of $\tilde X_{j_2}$ on the inner Belavkin filter is zero), this does NOT imply that $\tilde X_{j_1}^{\text{op}} \cdot M^{(j_2)} \cdot \tilde X_{j_1}^{\text{op}}$ has zero trace — the squared $\tilde X_{j_1}^2$ is positive (non-zero second moment), and $M^{(j_2)}$ is a Kraus channel that preserves positivity.

✓ **Belavkin predicts row (d) NON-ZERO.** Matches Syracuse's $0.108$ qualitatively.

### Comparison with AFL's failure

In AFL, the operator at step $j_2$ was constrained to be the *-homomorphism image $j_{j_2}(b) \in \mathcal A$ of a fixed observable $b \in \mathcal O$, with the kernel condition $\mathcal P(b) = 0$ forcing $\varphi_0(a c \cdot \mathcal P^{j_2 - j_1}(b)) = 0$. This contraction depended on:
- (i) abelian embedding: $j_{j_1}(\mathcal O)$ and $j_{j_2}(\mathcal O)$ commute, allowing $j_{j_1}(a) j_{j_1}(c) = j_{j_1}(a c)$ via commutativity.
- (ii) kernel condition: $\mathcal P(b) = 0$ kills the intervening factor.

In Belavkin, neither (i) nor (ii) applies. The operators at different times are NOT *-homomorphism images of a fixed object; they are Kraus operators of independent unitary interactions $U_j$. Commutativity between $\tilde X_{j_1}^{\text{op}}$ at the two positions and $M^{(j_2)}$ is NOT assumed. The kernel-condition contraction does NOT occur.

---

## 4. Row (f) — alternating-index fourth moment

**Syracuse:** $\varphi(\tilde X_{j_1} \tilde X_{j_2} \tilde X_{j_1} \tilde X_{j_2}) = 0.609$ for $j_1 < j_2$.

**Belavkin prediction.** The same structural argument as row (d), one step further. The operator chain
$$ \tilde X_{j_2}^{\text{op}} \cdot M^{(j_1)} \cdot \tilde X_{j_1}^{\text{op}} \cdot M^{(j_2)} \cdot \tilde X_{j_2}^{\text{op}} \cdot M^{(j_1)} \cdot \tilde X_{j_1}^{\text{op}} $$
involves alternating $\tilde X_{j_1}, \tilde X_{j_2}$ with intervening Kraus operators. No kernel-condition contraction. Non-zero in general.

Heuristically:
$$ \varphi(\tilde X_{j_1} \tilde X_{j_2} \tilde X_{j_1} \tilde X_{j_2}) \sim \mathbb E\left[ \tilde X_{j_1}^2 \cdot M^{(j_2)} \cdot \tilde X_{j_2}^2 \right] $$
again positive in the natural Reading-B configuration where $\tilde X_j^2$ has non-zero second moment.

✓ **Belavkin predicts row (f) NON-ZERO.** Matches Syracuse's $0.609$ qualitatively.

### Magnitude consistency check

The brief notes row (f) $= 0.609$ is much larger than row (d) $= 0.108$. This ratio (~5.6×) should be reproduced by the Belavkin model if the framework fits. Heuristically:
$$ \frac{\text{row (f)}}{\text{row (d)}} \sim \frac{\mathbb E[\tilde X_{j_2}^2]}{\mathbb E[\text{Kraus channel from } j_1 \text{ to } j_2]} . $$
The denominator is the "transition kernel" at step $j_2 - j_1$, contracting toward the invariant state; the numerator is the second moment of $\tilde X_{j_2}$. Without specific calculation, the ratio is a free parameter in Belavkin; nothing forces it to a particular value. **This is a Mode-E gap** — explicit verification requires computing the Kraus channel structure of $T_j$ in the (1, 4)-eigenvector basis of $T_{\text{diag}}$ (R77).

The structural verdict "row (f) is non-zero of comparable magnitude to row (d)" holds, but the specific ratio $5.6×$ is an open numerical check.

---

## 5. Fubini constant

**Syracuse:** $F(v_1, v_1') = 6.347 \times 10^{-2}$ — **constant** across all 12 grid points in the $(v_1, v_1')$ fiber.

**Belavkin prediction.** The Fubini inner factor is:
$$ F(v_1, v_1') = \mathbb E_{(v_2)}\left[ \tilde X_{j_2}^{\text{op}}(v_2, v_2') \cdot \tilde X_{j_1}^{\text{op}}(v_1, v_1') \cdot \tilde X_{j_2}^{\text{op}}(v_2, v_2') \right] $$
where the expectation is over the Geom(2) measure at step $j_2$.

In Belavkin terms, the inner integration produces the **Kraus channel** at step $j_2$ applied around $\tilde X_{j_1}^{\text{op}}$:
$$ F(v_1, v_1') = \sum_{v_2, v_2'} \mu_{\mathrm{Geom}}(v_2) \mu_{\mathrm{Geom}}(v_2') \cdot (M_{v_2}^{(j_2)})^* \tilde X_{j_1}^{\text{op}}(v_1, v_1') M_{v_2}^{(j_2)} . $$

This is an operator on $\mathcal H_n$, and Syracuse takes its sum_entries (scalar reduction).

**For $F$ to be CONSTANT in $(v_1, v_1')$**, the inner channel must produce an operator whose sum_entries depends on $(v_1, v_1')$ only through the **invariant 1-d eigenspace** of the iterated Kraus channel.

This is precisely Syracuse's R77 result: $T_{\text{diag}}$ has a 1-dimensional eigenspace at eigenvalue 1, spanned by $(1, 4)$. The asymptotic operator content collapses onto this 1-d invariant, and the resulting sum_entries scalar is INDEPENDENT of the fiber $(v_1, v_1')$.

In Belavkin's framework, this corresponds to the Kraus channel $\mathcal M^{(j_2)}(\rho) = \sum_v M_v^{(j_2)} \rho (M_v^{(j_2)})^*$ being **ergodic** — having a unique invariant state (or 1-d invariant subspace under the appropriate normalization).

✓ **Belavkin predicts Fubini constancy** under the assumption of ergodic Kraus channels at level $\ge j_2 - j_1 = 1$. Syracuse's R77 (T_diag has 1-d invariant eigenspace at $(1, 4)$) is exactly the ergodicity condition. Match.

---

## 6. Composite scoring table

| Row | Syracuse target | Belavkin prediction | Status |
|---|---|---|---|
| (b) $\varphi(\tilde X_{j_1} \tilde X_{j_2})$ | $\approx 0$ | $= 0$ via tower property (Belavkin martingale-difference property of centered system observables) | ✓ |
| (d) $\varphi(\tilde X_{j_1} \tilde X_{j_2} \tilde X_{j_1})$ | $0.108$ | NON-ZERO (no kernel-condition contraction; operator product $\tilde X_{j_1} M^{(j_2)} \tilde X_{j_1}$ does not vanish) | ✓ qualitative; numerical value Mode-E gap |
| (f) $\varphi(\tilde X_{j_1} \tilde X_{j_2} \tilde X_{j_1} \tilde X_{j_2})$ | $0.609$ | NON-ZERO; magnitude ratio to row (d) is a free parameter | ✓ qualitative; ratio Mode-E gap |
| Fubini $F(v_1, v_1')$ | const $6.347 \times 10^{-2}$ | const under ergodic Kraus channel at step $\ge 1$ (Syracuse's R77 1-d invariant on $(1, 4)$) | ✓ qualitative + structural anchor via R77 |

**All four targets accommodated by Belavkin.** ✓

---

## 7. Comparison to AFL prediction table

| Row | AFL | Belavkin |
|---|---|---|
| (b) | ✓ via $\mathcal P(a_*) = 0$ kernel condition | ✓ via tower property (no kernel condition needed) |
| (d) | ✗ FORCED ZERO under abelian filtration | ✓ non-zero (no kernel-condition contraction) |
| (f) | ✗ FORCED ZERO under abelian filtration | ✓ non-zero |
| Fubini | ✓ ergodic $\mathcal P$ | ✓ ergodic Kraus channel (same property, different language) |

**Belavkin strictly dominates AFL on rows (d) and (f), tied on (b) and Fubini.**

---

## 8. The structural reason Belavkin succeeds

In one sentence: **Belavkin's level-graded adaptive Kraus operators $M_v^{(j, b_{[1, j-1]})}$ are not *-homomorphism transports of a fixed object, so the kernel-condition contraction that killed AFL's rows (d) and (f) does not occur in Belavkin.**

More precisely, Belavkin allows:
1. **Abelian observation filtration** (matches Syracuse's $\mathbb B$).
2. **Non-commutative system algebra** (matches Syracuse's $\mathcal A$ generated by $T_j$).
3. **Level-graded operators that are NOT time-translates** (matches Syracuse's $T_j$ depending on $j$ via Tao's recursion).
4. **Adapted phase coupling INSIDE the operator** (matches Syracuse's $\chi_j(b_{[1, j-1]})$ phase factor).

The combination (1)+(2)+(3)+(4) is precisely what AFL's *-homomorphism-of-fixed-$\mathcal O$ embedding could NOT supply. Belavkin's filtering framework was designed for exactly this setup — classical (abelian) observation feeding back adaptively into non-commutative quantum dynamics.

---

## 9. Mode-E gaps in moment predictions

| Gap | Description |
|---|---|
| MP-G1 | Numerical verification of $0.108$ specifically (vs just "non-zero") would require explicit construction of the Kraus channel $\mathcal M^{(j_2)}$ at step $j_2$ in the $(1, 4)$ basis of R77's $T_{\text{diag}}$ and computing the trace $\mathrm{Tr}(\tilde X_{j_1}^2 \cdot \mathcal M^{(j_2)})$. **This is the natural next computational check.** |
| MP-G2 | Same for row (f) = $0.609$ — derive from a closed-form trace computation. |
| MP-G3 | Verify the Fubini constancy is exactly $6.347 \times 10^{-2}$ from R77's eigenstructure + Kraus channel structure. The structural prediction is "constant"; the specific value $6.347 \times 10^{-2}$ is a numerical check involving R77 + R74 + R75 prefactors. |
| MP-G4 | Verbatim citation of BvHJ 2009 Section 4 filter equation in the equation-number-and-page format. The structure I derived is canonical, but Mode-E discipline asks for verbatim sourcing. |

**Mode-E status.** The qualitative verdict (Belavkin predicts non-zero rows (d), (f), zero row (b), constant Fubini) is structurally clean and matches Syracuse. The quantitative match (specific numerical values $0.108$, $0.609$, $6.347 \times 10^{-2}$) requires explicit computation of Kraus channel actions, deferred. This is the natural next computational pass beyond this disposition.

---

## 10. References

- Targets: `D1_DISPOSITION.md`
- Belavkin framework: `BELAVKIN_VERBATIM.md`
- Identification: `BELAVKIN_SYRACUSE_IDENTIFICATION.md`
- Prior probes: `AFL_MOMENT_PREDICTIONS.md`, `AFL_DISPOSITION.md`, `QSC_DISPOSITION.md`
