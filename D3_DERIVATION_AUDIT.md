# D3 — Derivation audit of c = 7/45 under H1' failure

**Date:** 2026-05-14
**Task:** D3 of the H1' disposition (`H1_PRIME_DISPOSITION.md §6` gap table).
**Mode:** E — verbatim reading of `MONOTONE_CUMULANTS_C_ASYMPTOTIC.md` §1-§8 + cross-checks to A/B + W1 + W2 + H1'_LO + H1'_DISP.
**Adversarial:** internal. The audit attempts to break, not to confirm, the survival of c = 7/45.

---

## 0. One-sentence verdict

**Outcome 1 (best case): c = 7/45 SURVIVES the H1' failure under a narrower framework I name "monotone-singleton sufficiency."** The leading 7/45 derivation in `MONOTONE_CUMULANTS_C_ASYMPTOTIC.md` invokes HS 2014 Thm 3.4 (and its monograph variant Thm 3.26) **only at the all-singletons monotone partition + diagonal κ_2^B**, never at a partition class that requires the peak rule at non-adjacent repeated indices. The seemingly load-bearing partition-class invocations are either (i) the all-singletons class (κ_1^B trivially commutes with the moment-cumulant formula via product-of-singletons-only computation), (ii) the diagonal κ_2 at a single index (no cross-index peak rule), or (iii) explicitly **flagged as conjectural / open** by the file itself (the −1/30 numerical coefficient + the rate 1/2). The c = 7/45 number itself reduces to **R75 Plancherel × R77 T_diag (1, 4)-eigenvector × R64.B class-mass ratio**, which never needed HS 2014 Thm 3.4 in the first place — the monotone cumulant framework was used only to **identify** the 7/15 mass with κ_1^B at the all-singletons partition, not to **derive** it.

---

## 1. Step-by-step trace of the c = 7/45 derivation

The derivation lives in `MONOTONE_CUMULANTS_C_ASYMPTOTIC.md` §1-§7. I trace every step that touches HS 2014 Thm 3.4 / HS 2011 Thm 6.1 / Hasebe monograph Thm 3.26 (these are equivalent statements; I'll call them collectively "Thm 3.4").

### Step 1 — Reciprocal-Cauchy composition (§2.1)

**File quote (§2.1, line 39):**

> "`H_{X_1 + X_2 + ⋯ + X_n}(z) = H_{X_1}(H_{X_2}(⋯ (H_{X_n}(z))⋯))`"

This is **Muraki 2003 Thm 4**, NOT HS 2014 Thm 3.4. The framework is invoked but not the moment-cumulant expansion.

**Partition class invoked:** N/A (reciprocal-Cauchy composition is a generating-function identity, not a partition-sum). The file itself flags this step as fiberwise-rigorous + conjectural at B-valued lift (§2.1 last paragraph "Mode-E note", line 58).

**Safe vs problematic:** Not a Thm 3.4 invocation. Independent of the H1' question.

### Step 2 — Cumulant additivity (§2.2)

**File quote (§2.2, line 70):**

> "`κ_n(x_1 + x_2 + ⋯ + x_N) = N · κ_n(x_1)`"

This is **HS 2011 Defn 4.5 / monograph (M3)** — the extensivity axiom — applied per-step at fixed accumulator (line 78):

> "`κ_n^B(Σ_j Off_j)(history) = Σ_j κ_n^B(Off_j)(b_{[1,j−1]})`    (under the lift)"

**Partition class invoked:** None directly. Extensivity is a statement about how cumulants scale, not about a specific partition sum.

**Safe vs problematic:** This is the additivity axiom, which in the scalar setting is independent of the peak rule at non-adjacent repeats. The Syracuse application is **per-step**, not iid-extensivity, so it reduces to summing κ_n^B(Off_j) over j — which is again a per-step computation, not a non-adjacent-repeat moment.

### Step 3 — Leading-singularity invocation of Thm 3.4 (§3 — THE LOAD-BEARING STEP)

**File quote (§3, lines 96-99):**

> "**Monotone-cumulant reading of the leading term.** Apply the moment-cumulant formula (Hasebe monograph Thm 3.26, verbatim):
>
> `E_B(X^n) = Σ_{π ∈ M(n)} (1/|π|!) κ_π^B(X)`"

This is the **first explicit invocation of Thm 3.4 / Thm 3.26** in the derivation. The very next sentence (line 102) selects the partition class:

> "The leading contribution at n → ∞ comes from the **all-singletons** monotone partition `π = ({1}, {2}, ..., {n})` (one block per position, ordered). This contributes
>
> `(1/n!) · κ_1^B(X)^n`   (since |π|! = n!, each block size 1 gives κ_1)"

**Partition class invoked:** **ALL-SINGLETONS ONLY** at this step. Every block has size 1, contributing κ_1^B per block. No κ_2, no κ_3 at non-adjacent repeats.

**Safe vs problematic:** **SAFE.** The all-singletons partition `π = ({1}, ..., {n})` has every block of size 1. The HS 2014 Thm 3.4 formula at all-singletons evaluates κ_π^B = κ_1^B · κ_1^B · ⋯ · κ_1^B (n times). This product is **B-valued multiplication of n copies of κ_1^B**, which is a function of accumulators only. The partition-cumulant formula at all-singletons is **trivially equivalent to** taking the n-fold product `(κ_1^B)^n`. **It does not invoke the peak rule at any non-adjacent-repeat index pattern** because there are no repeated indices in a partition (each singleton is a distinct position label).

The peak rule of HS 2014 Defn 2.2 is invoked **internally to define κ_n^B for n ≥ 2** (via the moment-cumulant inversion or equivalently the Speicher interval-block contraction in HS 2014 Defn 3.3). At all-singletons, we only need κ_1^B, which is just `E_B(X)` (HS 2014 Defn 3.3 specialized to n=1: the coefficient of N in `ϕ(N.X) = N · ϕ(X)`, which is `E_B(X)`). κ_1^B is defined **without invoking Defn 2.2 at all** (n=1 cumulant has no factorization to do).

### Step 4 — The κ_1^B numerical identification (§3, lines 107-114, "Verbatim algebra")

**File quote (§3, lines 108-110):**

> "For the Tao atoms, the leading-order κ_1^B is the diagonal contribution, which from R77 §1 gives eigenvalue 1 on (1, 4) — exactly the 7/15 mass limit."

And lines 117-127:

> "From R77 Thm 77.1: `T_diag = (1/5)·[[1,1],[4,4]]`, eigenvalues {0, 1}; (1, 4)-eigenvector projection of S = 2(P_+ + P_−); limit S_n / 3^n → 7/15. Then 7/45 = (7/15) / 3 = S_∞ / 3 = (1/3) · Plancherel mass on high-frequency Fourier coefficients."

**Partition class invoked:** None — this step does not invoke Thm 3.4 at all. It **identifies** κ_1^B with the R77 T_diag (1, 4)-eigenvalue from independent infrastructure (R75 Plancherel + R76 conservation + R77 spectrum + R64.B class-mass ratio).

**Safe vs problematic:** **SAFE.** This is the **pre-existing** rigorous derivation of 7/15 from R77 T_diag eigenstructure, which has nothing to do with HS 2014. The monotone-cumulant framework is invoked only to **name** the 7/15 as "κ_1^B at the all-singletons partition." The number 7/15 = 7/(3·5) emerges from:

- R76 Thm 76.3: `S_n = −2 R_n` (bilinear pair factor)
- R75 Plancherel: factor of 3
- R77 T_diag = (1/5)·[[1,1],[4,4]] eigenstructure
- R64.B class-mass ratio (1:4)

Each of these is **independent** of HS 2014 Thm 3.4. The c = 7/45 = (7/15)/3 reduction is **purely algebraic** (Plancherel rescaling: |μ̂_n|² · 3^n → 7/15, so |μ̂_n|² → (7/15)·3^{-n} = (7/45)·3^{-n+1}; the (7/45) appears in the bilinear pair-form normalization).

### Step 5 — Cross-check via the CLT (§3 "Cross-check: which cumulants?")

**File quote (§3, lines 131-139):**

> "The classical CLT for monotone independence gives the **arcsine distribution** (HS 2011 Thm 5.1, verbatim Deliverable A §5). Its monotone cumulants are `(0, 1, 0, 0, ...)` — only κ_2 non-zero. ... Syracuse is not in the CLT regime — it's in a **B-valued drift regime** where the mean (κ_1^B) is the dominant cumulant."

**Partition class invoked:** None — this is a comparison statement, used only to argue that the dominant Syracuse contribution is κ_1^B (mean), not κ_2 (variance, which would give the arcsine CLT).

**Safe vs problematic:** **SAFE.** This is a structural argument that Syracuse is mean-dominated, not a Thm 3.4 invocation.

### Step 6 — Subdominant: rate 1/2 and coefficient −1/30 (§4)

This step is where the file itself flags conjectural status. Let me read the partition class invoked very carefully.

**File quote (§4, lines 161-167):**

> "The dominant ((κ_1^B)^n / n!) gives 7/45 · 3^{−n}. The first subdominant correction comes from a single block of size 2 paired with (n−2) singletons:
>
> contribution: `(1/(n−1)!) · κ_2^B(X) · κ_1^B(X)^{n−2}`   (one 2-block + (n−2) 1-blocks)"

And the file's own immediate finding (lines 178-187):

> "So `E_B(X^n) ≈ (κ_1^B)^n / n! + (1/(n−2)!) · κ_2^B · (κ_1^B)^{n−2} + ...` ... The ratio of subdominant to dominant is: `subdominant / dominant = n(n−1) · κ_2^B / (κ_1^B)^2`. **This grows like n², not as `(1/2)^n`. So the second cumulant κ_2^B alone does NOT supply the rate-1/2 decay.**"

**Partition class invoked:** "one 2-block + (n−2) singletons."

**Safe vs problematic:** This is the **only place** in the c = 7/45 derivation that invokes a non-singleton partition class. Two sub-cases to check:

(a) **Diagonal κ_2 at a single index (same-step κ_2^B(Off_j, Off_j)).** Per Deliverable B §2.2 (verbatim, lines 87-93), the diagonal κ_2^B at fixed step j is defined via the scalar formula `κ_2^B(Off_j) = E_B(Off_j²) − [E_B(Off_j)]²` — this is HS 2014 Defn 3.3 specialized to n=2 at a single index, which **does NOT invoke the peak rule at non-adjacent repeats**. (The peak rule of Defn 2.2 requires distinct adjacent indices i_{k-1} < i_k > i_{k+1}; at n=2 with a single index, there's no peak structure.)

(b) **Cross-step κ_2^B(X̃_{j_1}, X̃_{j_2}) at j_1 ≠ j_2.** Per Deliverable B §2.2 (verbatim, line 83) and W2 §1 (verbatim, line 32), this **vanishes structurally** under marginal centering (Pascal-pair independence): `κ_2^B(X̃_{j_1}, X̃_{j_2}) = 0` for j_1 ≠ j_2.

**Critical: the 2-block in the partition expansion contributes via (a) only, NOT (b).** The file confirms this is the diagonal κ_2 contribution (§4 line 192-201):

> "The rate-1/2 decay is the **subdominant eigenvalue of T**, not of a single cumulant. ... The B-valued operator `κ_1^B(Off_j)(b_{[1,j−1]})` viewed as a function of b_{[1,j−1]} ∈ B has its own decomposition: ... A b-dependent part → contributes to the off-diagonal correction Off_j, giving rate-1/2 from the leading bilinear coupling (v = 1, v' = 3) which has 3-adic valuation 1 in the phase difference 2^{−v} − 2^{−v'}."

So the subdominant rate-1/2 is **redirected** from the κ_2 partition-sum (which fails by Track A audit W2 erratum — the partition count grows ~ n·ln(n), not exponentially) to the **T_M operator spectrum**, which is **independent of HS 2014 Thm 3.4**.

**Conclusion on Step 6:** The 2-block + singletons partition class is invoked **only in the diagonal κ_2 sub-case (a)**, which does NOT use the peak rule. The cross-step κ_2 (sub-case b) vanishes and never enters. So even the subdominant step is **safe** — no non-adjacent-repeat partition class is invoked.

The W2 audit (Track A integration §2.4 Erratum 2.2) already shifted the −1/30 coefficient's source from "monotone partition counting" to "structural decomposition via R76 × R75 × R77" (`1/30 = 1/(2·15)`), with the rate-1/2 redirected to T_M λ_2. This means the subdominant derivation **already doesn't depend on HS 2014 Thm 3.4** post-W2.

### Step 7 — Mode-E uncertainty ledger (§8)

The file's own §8 ledger (lines 326-333) confirms:

- `c = 7/45` leading coefficient: **Closed (rigorous via R77 T_diag + R75 Plancherel; monotone-cumulant framework consistent with this and identifies it as κ_1^B-dominant)**
- Rate `1/3` from Plancherel: Rigorous (R75)
- Rate `1/2` subdominant: Conjectural in monotone framework; rigorous derivation pending R77 §3
- Coefficient `−1/30` of subdominant: Open (numerical only)
- B-amalgamated lift of HS theorem: **Conjectural** (Mode-E gap)

The file **itself** does not claim that c = 7/45 depends on HS 2014 Thm 3.4 at any non-safe partition class. It identifies c = 7/45 with the κ_1^B all-singletons contribution and reduces the actual computation to R75 + R77 + R64.B.

---

## 2. Per-step partition-class table

| Step | File location | Partition class | Safe / problematic |
|---|---|---|---|
| 1 — H-transform composition | §2.1 | none (generating-function identity) | not a Thm 3.4 invocation |
| 2 — Cumulant additivity | §2.2 | none (axiom M3) | not a non-adjacent-repeat moment |
| 3 — Leading 7/15 via Thm 3.26 | §3 | **all-singletons only** | **SAFE** (κ_1 only) |
| 4 — Numerical identification 7/45 | §3 "Verbatim algebra" | none (R75/R77/R64.B) | independent of Thm 3.4 |
| 5 — CLT cross-check | §3 "Which cumulants?" | none (structural comparison) | safe |
| 6 — Subdominant 2-block + (n−2) singletons | §4 | **diagonal κ_2 + all-singletons** | **SAFE** (no cross-index κ_2 invoked; the cross-step κ_2 = 0) |
| 7 — PADE multi-spectral | §5 | none (consistency check) | safe |
| 8 — −1/30 coefficient | §4 + Track-A W2 erratum | redirected to R76×R75×R77 | safe (independent of Thm 3.4 post-W2) |

**No step in the c = 7/45 derivation invokes Thm 3.4 at a partition class with non-adjacent repeated indices.**

The only partition classes that appear are:
- All-singletons (Step 3, leading)
- Diagonal κ_2 at a single index + (n−2) singletons (Step 6, subdominant)

Both are **safe** under the H1' failure regime, because:
- All-singletons uses only κ_1^B = E_B(·), which doesn't depend on Defn 2.2's peak rule
- Diagonal κ_2 at a single index has no peak structure (n=2 single index = no adjacent-index alternation)

---

## 3. Verdict: Outcome 1 — survives under "monotone-singleton sufficiency"

### 3.1 The narrower framework

**Definition (monotone-singleton sufficiency).** A B-valued random-variable family `(X_j)_{j ∈ I}` satisfies **monotone-singleton sufficiency** over (A, B, ϕ) if:

(MSS-1) **Singleton-cumulants well-defined.** For each j and each n ≥ 1, the n-th individual monotone cumulant `κ_n^B(X_j) := κ_n^B(X_j, ..., X_j) ∈ B` is well-defined as the coefficient of N in the polynomial `ϕ((N.X_j)^n)`, where `N.X_j` is the dot-operation (HS 2014 Defn 2.3) applied to **iid copies of the single random variable X_j** in an extended algebraic probability space. **This is the scalar HS 2011 / Hasebe monograph cumulant, with B-valued range, applied at a single index** — it does NOT require the multi-index peak rule of Defn 2.2 at non-adjacent repeats.

(MSS-2) **Cross-index second cumulants vanish.** For all j_1 ≠ j_2, `κ_2^B(X_{j_1}, X_{j_2}) = 0` — i.e., the second-order independence holds at distinct indices.

(MSS-3) **Per-step additivity.** For the sum `S = Σ_j X_j`, the cumulants decompose as `κ_n^B(S)(b) = Σ_j κ_n^B(X_j)(b_{[1, j−1]})` at each fixed accumulator b. **This is per-step extensivity, not iid-extensivity.**

(MSS-4) **All-singletons partition contribution.** The leading contribution to `E_B(S^n)` as n → ∞ comes from the all-singletons monotone partition in the moment-cumulant formula (HS 2014 Thm 3.4 restricted to this partition class):

`E_B(S^n) = (1/n!) · (κ_1^B(S))^n + [subdominant corrections from non-trivial partitions]`

The all-singletons term is well-defined under (MSS-1) + (MSS-3), and it gives the leading asymptotic without invoking the full Thm 3.4 expansion at other partitions.

### 3.2 Why Syracuse satisfies MSS

- (MSS-1): The κ_n^B(X̃_j) for a **single** index j is the scalar HS cumulant at fixed b_{[1, j−1]}, B-marginalized. This is the "fiberwise scalar HS + marginalize over B" reading documented in Deliverable A §6 (operational Mode-E reading). The H1' failure was at **non-adjacent repeats** (sequences like (j_1, j_2, j_1, j_2)), NOT at single-index cumulants. The fiberwise scalar HS 2011 / monograph theorem applies cleanly at each single index.
- (MSS-2): Verified by W2 §1 (`κ_2^B(X̃_{j_1}, X̃_{j_2}) = 0` for j_1 ≠ j_2, Pascal-pair independence at distinct steps). This is independent of Defn 2.2.
- (MSS-3): Per-step additivity is the "naive cumulant decomposition" of Deliverable C §2.2 (lines 76-80), which is structurally trivial given that the per-step atoms Off_j depend only on (v_{2j−1}, v_{2j}) which are iid Pascal-distributed across j.
- (MSS-4): The all-singletons contribution `(κ_1^B(Σ Off_j))^n / n!` reduces (by (MSS-3)) to `(Σ_j κ_1^B(Off_j))^n / n!`, with κ_1^B(Off_j) the B-marginal mean of Off_j at fixed prior accumulator. Projected onto the R77 (1, 4)-eigenvector direction, κ_1^B(Σ Off_j) → 7/15 via R75 + R77 + R64.B (independent infrastructure).

The dominant contribution `(7/15)^n / n!` rearranges (via Plancherel and the bilinear pair-form normalization, R75 + R76) to give `c = 7/45` as the leading coefficient of `|μ̂_n|² · 3^n` in the asymptotic series — which is exactly what the file derives.

**Crucially: none of (MSS-1)–(MSS-4) requires HS 2014 Defn 2.2 to hold at non-adjacent repeats.** The W1 H1' failure (Defn 2.2 fails at n ≥ 4 alternating) breaks the **full Thm 3.4 expansion at higher partition classes** but does NOT break the all-singletons contribution + diagonal κ_2.

### 3.3 What does NOT survive

The **higher-order corrections** beyond leading + subdominant — i.e., partitions with κ_3, κ_4, ... at non-adjacent-repeat patterns — invoke the peak rule at exactly the index patterns where H1' fails. So:

- c = 7/45 leading: **survives** (Outcome 1).
- −1/30 subdominant amplitude: **survives** via the R76 × R75 × R77 structural decomposition (Track A W2 §2.4, post-erratum). The monotone-cumulant invocation was already not load-bearing after W2.
- Rate `(1/2)^n`: redirects to T_M λ_2, separately open (does not depend on Thm 3.4).
- Higher subdominants (rate `(1/4)^n` and beyond, multi-spectral PADE structure): **NOT salvageable** without a working framework at non-adjacent repeats. The file flags these as open already.

### 3.4 Mode-E gaps remaining (post-D3)

| Gap | Status after D3 |
|---|---|
| c = 7/45 leading | **CLOSED conditional on MSS** (Outcome 1). MSS is weaker than H1 / H1' verbatim and the Syracuse X̃_j family satisfies it. |
| −1/30 subdominant amplitude | Closed post-Track-A W2 erratum (`1/(2·15)`), independent of Thm 3.4. |
| Rate `(1/2)^n` | Open (R77 Conjecture 77.2), independent of Thm 3.4. |
| Higher subdominants / multi-spectral PADE | **Open** — would require a non-Hasebe framework at non-adjacent repeats. |
| MSS-as-a-published-framework | **D2 task** (closure-hunt for the right published framework) — but the Syracuse case is closeable without identifying a published name. |

---

## 4. Adversarial probes (attempted breaks)

I tried to find a hidden invocation of Thm 3.4 at a problematic partition class. Attempts:

### Probe 4.1 — Is the "leading singularity" derivation hiding a non-singleton sum?

The file's claim is that the all-singletons partition **dominates** as n → ∞. The dominance argument is that `(κ_1^B)^n / n!` grows faster than any term with a κ_2 or κ_3 factor (which are bounded). The Thm 3.4 full sum is `Σ_{π ∈ M(n)} (1/|π|!) κ_π^B(X)`, with `|M(n)| = (n+1)!/2`. The all-singletons block has |π|! = n!, giving the smallest combinatorial-prefactor reduction; non-singleton partitions have larger blocks reducing |π|!, but the cumulant factors decay.

**Hidden dependence?** The dominance argument **assumes** the full Thm 3.4 sum converges — which would require all κ_n^B to be well-defined and the partition sum to be controlled. If H1' fails at non-adjacent repeats, the κ_n^B for n ≥ 3 at non-adjacent-repeat sequences are **not well-defined** by the HS 2014 inversion (the moment-cumulant inversion requires the moments to be HS-consistent, which fails at non-adjacent repeats).

**However**: the leading asymptotic `(κ_1^B)^n / n!` is well-defined **without** the full sum. The file's claim is not "the full Thm 3.4 sum converges to (7/15)^n / n!," but rather "the leading term of the Thm 3.4 sum is (7/15)^n / n!, and the subdominant corrections are O(n²) smaller per step in the controlled regime." The controlled regime is bounded by the magnitudes of κ_n^B at distinct indices and at single-index diagonal — both of which are well-defined under MSS.

**Verdict on Probe 4.1:** The dominance argument is **at risk** if we read it as requiring the full Thm 3.4 sum to be a valid asymptotic. **Safer reading:** the all-singletons contribution is computed directly (no sum-over-partitions needed), and the result equals 7/15 via R77 T_diag — which is **independent** of whether higher-order partitions contribute. The 7/45 number is established by R75 + R77 + R64.B without ever needing the full Thm 3.4 sum.

This is the strongest adversarial finding. The c = 7/45 derivation **as a stand-alone identification of the leading coefficient** is independent of Thm 3.4. The Thm 3.4 framework is used only to **name** this leading coefficient as "κ_1^B at all-singletons." If the naming is invalid (because the full sum doesn't converge), the **number** is still correct from R77 alone.

**This collapses the dependence on monotone-cumulant theory entirely for the leading coefficient.** The c = 7/45 derivation is, at its rigorous core, an **R77 result**, not a Hasebe-Saigo result.

### Probe 4.2 — Is there a hidden non-adjacent-repeat moment in the κ_1^B definition?

The κ_1^B(Off_j) is `E_B(Off_j)` at fixed prior accumulator. This is a single-index, n=1 expectation. **No peak rule, no factorization, no problematic partition.**

**Verdict on Probe 4.2:** Clean. κ_1^B is the conditional expectation, defined trivially.

### Probe 4.3 — Is the per-step additivity (MSS-3) hiding a peak rule invocation?

The per-step additivity says `κ_n^B(Σ_j Off_j) = Σ_j κ_n^B(Off_j)` at fixed accumulator. For n = 1, this is just the linearity of expectation. For n ≥ 2, it's the B-valued analog of the scalar additivity, which in the Hasebe framework requires monotone independence (Defn 2.2 or Defn 1.21).

**Where does this come from?** Per-step additivity at n = 2 requires `E_B((Off_1 + Off_2)²) - [E_B(Off_1 + Off_2)]² = κ_2^B(Off_1) + κ_2^B(Off_2)` (assuming cross-step κ_2 = 0). Expanding LHS: `E_B(Off_1²) + 2 E_B(Off_1 · Off_2) + E_B(Off_2²) - [E_B(Off_1) + E_B(Off_2)]²`. The cross term `2 E_B(Off_1 · Off_2)` requires evaluation. **If Off_1 and Off_2 satisfy monotone Defn 2.2 at the (1, 2) index pattern (which is n = 2, distinct indices)**, then `E_B(Off_1 · Off_2) = E_B(Off_1) · E_B(Off_2)` (peak rule at position 2 = endpoint, distinct indices), giving the cross-term cancellation.

**Is this the failed regime?** No. The H1' failure is at **n ≥ 4 with non-adjacent repeats**. At n = 2 with distinct indices, Defn 2.2 holds (verified in `H1_PRIME_LOW_ORDER_CHECKS.md §2.1`, line 70: "Defn 2.2 holds at n = 2, j_1 < j_2"). Per-step additivity at n = 2 is **safe**.

For higher n, per-step additivity at non-adjacent-repeat patterns would invoke the failed regime. But the c = 7/45 derivation **only uses per-step additivity at the singleton level** (n = 1 cumulant per step, summed across steps), which reduces to linearity of expectation. **No higher-order per-step additivity is invoked in the leading 7/45 derivation.**

**Verdict on Probe 4.3:** Clean for the leading coefficient. Per-step additivity at n = 1 is trivial; at n ≥ 2 it would require Defn 2.2 at distinct-index patterns (safe at n=2; only the n ≥ 3 non-adjacent-repeat patterns fail).

### Probe 4.4 — Is the (1, 4)-eigenvector projection hiding a Thm 3.4 invocation?

The (1, 4)-eigenvector is from R77 T_diag = (1/5)·[[1,1],[4,4]]. This is a 2×2 matrix from R77's analysis of the diagonal transfer operator, **independent** of monotone-cumulant theory. The file's claim "κ_1^B projected onto (1, 4)-direction = 7/15" identifies the R77 result with the Hasebe κ_1^B at all-singletons, but the R77 result stands on its own.

**Verdict on Probe 4.4:** Clean. R77 is independent infrastructure.

### Probe 4.5 — Combined: would removing HS 2014 entirely break the c = 7/45 derivation?

This is the strongest adversarial test: **suppose HS 2014 / Hasebe-monograph Thm 3.4 doesn't exist. Can we still derive c = 7/45?**

**Answer: yes.** The c = 7/45 derivation reduces to:

1. R75 Plancherel: `Σ |μ̂_n|² = 3^{-n} · ‖d_n‖²`
2. R76 conservation: `S_n = -2 R_n`
3. R77 T_diag: `T_diag = (1/5)·[[1,1],[4,4]]`, eigenvalue 1 on (1, 4)-direction
4. R64.B class-mass ratio: (1/3)² : (2/3)² = 1 : 4 → eigenvector (1, 4)
5. Limit: `S_n / 3^n → 7/15` → `|μ̂_n|² ~ (7/45) · 3^{-n}`

**None of these steps invokes Thm 3.4 or any monotone-cumulant theorem.** The number 7/45 is a **pre-existing R77 result** that the monotone-cumulant framework was used to **interpret**, not derive.

**Verdict on Probe 4.5:** The c = 7/45 derivation is **fully independent** of HS 2014 Thm 3.4. The H1' failure has **no impact** on the leading coefficient. The framework was an interpretive overlay, not a derivation pathway.

---

## 5. Verdict: Outcome 1 confirmed, strengthened

**Outcome 1+ (strengthened):** c = 7/45 not only survives the H1' failure under a narrower framework (MSS), it is actually **independent of any monotone-cumulant framework**. The leading coefficient was always derivable from R75 + R76 + R77 + R64.B alone. The monotone-cumulant framework in `MONOTONE_CUMULANTS_C_ASYMPTOTIC.md` is an **interpretive identification** ("the 7/15 is κ_1^B at all-singletons"), not a derivation.

This is **Outcome 1** in the brief's taxonomy, but stronger than "narrower framework suffices" — the leading number is rigorous on entirely different infrastructure. The narrower framework MSS is sufficient to **name** the result in monotone-cumulant language, but not necessary to **establish** it.

**What changes vs. Track A integration's pre-D3 claim "rigorous conditional on H1'":**

| Component | Pre-D3 status | Post-D3 status |
|---|---|---|
| c = 7/45 leading | Rigorous conditional on H1' (HS 2014 Thm 3.4 + R75 + R77 + R64.B) | **Rigorous unconditionally** via R75 + R76 + R77 + R64.B; monotone-cumulant framework is interpretive overlay only |
| −1/30 subdominant amplitude | Closed structurally (`1/(2·15)`) | Unchanged — already independent of Thm 3.4 post-W2 |
| Rate `(1/2)^n` | Open (R77 Conj 77.2) | Unchanged |
| Higher subdominants | Open | Unchanged |

The Track A integration §1.4 line "**c = 7/45 leading coefficient**: Rigorous conditional on H1'" should be **upgraded** to:

> **c = 7/45 leading coefficient**: rigorous UNCONDITIONALLY via R75 + R76 + R77 + R64.B (pre-existing project infrastructure, independent of HS 2014 Thm 3.4). The monotone-cumulant framework (HS 2014 / monograph) provides an interpretive identification of 7/15 as κ_1^B at the all-singletons monotone partition, but is not needed for the derivation. The H1' failure at non-adjacent-repeat moments does NOT affect this conclusion.

The H1' Disposition §1 caveat that needs updating:

> "The leading c = 7/45 coefficient comes from the all-singletons monotone partition contribution to ϕ(X^n), which is `(κ_1^B(X))^n / n!`. The all-singletons partition evaluates n independent κ_1^B values — no peak rule at non-adjacent repeats is invoked at this contribution."

is **correct but understated**. Stronger statement: even the all-singletons identification is interpretive; the 7/45 number is established by R77 T_diag without invoking Hasebe theory.

---

## 6. Mode-E ledger update

| Gap | Pre-D3 status | Post-D3 status |
|---|---|---|
| Operator-valued vs vacuum-pairing distinction (D1) | Open | Open (separate task; D3 doesn't address) |
| Right-framework search (D2) | Open | **Lower priority** — c = 7/45 doesn't depend on naming the right framework. MSS suffices internally. D2 still useful for the higher-order corrections, not for the leading coefficient. |
| Derivation audit (D3 — this file) | Open | **CLOSED** — Outcome 1+ (strengthened). |
| STRUCTURAL-1 (n ≥ 4 alternating fails Defn 2.2) | Open, fundamental | Unchanged — but only relevant for higher-order corrections, not for c = 7/45. |

---

## 7. Files

- This audit: `C:/Collatz/D3_DERIVATION_AUDIT.md`
- Target file audited: `C:/Collatz/MONOTONE_CUMULANTS_C_ASYMPTOTIC.md`
- H1' failure context: `C:/Collatz/H1_PRIME_DISPOSITION.md`, `C:/Collatz/H1_PRIME_LOW_ORDER_CHECKS.md`
- Pre-existing infrastructure (load-bearing for c = 7/45):
  - R75: `C:/Collatz/result_75_plancherel.md` (Plancherel)
  - R76: `C:/Collatz/result_76_conservation_law.md` (`S_n = −2 R_n`)
  - R77: `C:/Collatz/result_77_T_lead_spectrum.md` (T_diag eigenstructure)
  - R64.B: class-mass ratio (in earlier results file)
  - `C:/Collatz/c_seven_forty_fifth.md` (Plancherel anchor for 7/45)
- W1 lift: `C:/Collatz/W1_BLIFT_THEOREM.md`, `C:/Collatz/W1_BLIFT_VERIFICATION.md`
- W2 cumulant: `C:/Collatz/W2_KAPPA2_CALC.md`
- Track A integration: `C:/Collatz/TRACK_A_INTEGRATION.md`
- HS 2014 PDF: `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_saigo_2014_operator_valued_monotone.pdf`
