# Result 12 (qx+1 paper) — literature: Tao 2019 does NOT give our result, and the reason NAMES our criticality. `H₁` vs `H₂`, and `D₂ = 1 ⟺ q = 3`.

**Date:** 2026-07-15. **Type:** literature read + framing (no probe). **Verdicts: Tao Prop 1.17 does NOT apply (different quantity) / Tao Remark 1.15 IS our address-count architecture, with a DIFFERENT entropy / ★ the contrast gives the paper its front-matter statement.**

**Headline: Tao counts addresses by SHANNON entropy (`H₁ = log 4` → 4ⁿ tuples); our L² object is governed by COLLISION entropy (`H₂ = log 3` → 3ⁿ). `H₂(Geom(1/2)) = log 3 = log q` exactly at q=3 ⇒ `D₂ = 1 ⟺ q = 3`. Collatz sits exactly at the critical point of its own family.**

Source: `archemedian_eigenvalues/pdfs/arxiv_1909.03562_Tao_Collatz.pdf` (arXiv:1909.03562v5, *Almost all orbits of the Collatz map attain almost bounded values*). Read on the R11 routing "READ TAO before deriving anything else."

## 1. Tao's Prop 1.17 is NOT our result — do not cite it as if it were

> **Prop 1.17 (Decay of characteristic function).** Let `n ≥ 1` and `ξ ∈ Z/3ⁿZ` not divisible by 3. Then `|E e^{−2πiξ·Syrac(Z/3ⁿZ)/3ⁿ}| ≪_A n^{−A}` for any fixed `A > 0`, **uniform in n and ξ**.

This is a **pointwise** superpolynomial bound on *individual* high-frequency coefficients. Our object is the **ℓ² sum**:

&nbsp;&nbsp;&nbsp;&nbsp;`M_n = Σ_{3∤ξ} |μ̂_n(ξ)|² → 7/15`

**Summing Prop 1.17 over the ~2·3^{n−1} admissible frequencies gives `≪ 3ⁿ·n^{−2A}`** — an upper bound that *grows* like 3ⁿ, and is therefore worthless against a sum converging to a constant. **Prop 1.17 cannot deliver our result at any A.** The quantities are genuinely different: his is pointwise high-frequency decay (equivalently, by his Remark 1.18, a TV/oscillation statement `Osc_{m,n} ≪_A m^{−A}`, Prop 1.14); ours is the ℓ² mass of the finest-scale deviation.

His method (Section 7: a two-dimensional renewal process interacting with a union of well-separated "triangles", after conditioning on the Pascal-distributed pair sums `b_j = a_{2j−1}+a_{2j}`) is built for the pointwise bound. **It is not an overlap count** and does not obviously adapt.

## 2. ★ Tao's Remark 1.15 IS our address-counting architecture — with a different entropy

Tao's heuristic justification of fine-scale mixing, verbatim:

> *"The geometric random variable Geom(2) can be computed to have a **Shannon entropy of log 4**; thus, by asymptotic equipartition, the random variable Geom(2)ⁿ is expected to behave like a uniform distribution on **4^{n+o(n)} separate tuples** in (N+1)ⁿ. On the other hand, the range Z/3ⁿZ of the map `a ↦ F_n(a) mod 3ⁿ` only has cardinality 3ⁿ. ... if one models this map by a random map from 4^{n+o(n)} elements to Z/3ⁿZ one is led to the estimate (1.23) (in fact this argument predicts a stronger bound of `exp(−cm)` for some c > 0, which we do not attempt to establish here)."*

**This is exactly R8's address picture** — addresses `(v_1..v_n)` mapping into `Z/3ⁿ`, and the question is how they collide. Tao counts them by **Shannon entropy**. Our L² object counts them by **collision (Rényi-2) entropy**. Both are exact:

| | definition | value for Geom(1/2) | address count |
|---|---|---|---|
| **Shannon** `H₁` | `−Σ_v p_v log p_v = E[v]·log2` | `2·log2 = **log 4**` | `4ⁿ` |
| **Collision** `H₂` | `−log Σ_v p_v² = −log Σ_v 4^{−v}` | `−log(1/3) = **log 3**` | `3ⁿ` |

(`H₂ ≤ H₁` always — Rényi entropy is non-increasing in its order. Here the gap is `log(4/3)`.)

**And that gap is precisely the gap between his problem and ours:**

| | governing entropy | addresses | vs `3ⁿ` slots | consequence |
|---|---|---|---|---|
| **Tao** (TV mixing) | `H₁ = log 4` | `4ⁿ` | `4ⁿ ≫ 3ⁿ` — **comfortable surplus** | mixing; heuristically `exp(−cm)` |
| **Ours** (ℓ² mass) | `H₂ = log 3` | `3ⁿ` | `3ⁿ = 3ⁿ` — **exactly critical** | `‖π_k‖² ~ (7/15)·k·3^{−k}`, the factor of k |

Tao's surplus is why mixing is *comfortable* for him. Our exact saturation is why q=3 is *critical* for us — and why R8's domination fails there by exactly a factor of k.

## 3. ★★ The front-matter statement

> **`H₂(Geom(1/2)) = −log Σ_{v≥1} 4^{−v} = log 3 = log q` — exactly, at q = 3.**
> Correlation dimension `D₂ = H₂ / log q = log 3 / log q`, and **`D₂ = 1 ⟺ q = 3`.**
> **Collatz sits exactly at the critical point of its own family.**

The `3` is `1/Σ_v p_v² = 1/E[2^{−v}] = 2² − 1` — a property of the **halving law alone** (base 2), carrying no information about the multiplier. The Collatz multiplier just *happens* to equal that number. That is the entire content of "3 is special", stated in one line of standard vocabulary.

**This retroactively explains R85 rung 2** (`result_85_bridge.md`): varying the multiplier `3 → p` gave `S_k(p) ~ (p/3)^k` with **the 3 fixed** across p=3,5,7,11,13. Of course it did — the 3 is the halving law's collision entropy and is multiplier-blind. Rung 2's refutation of `7 = Φ_p(2)` and its "rate 1/3 independent of multiplier" byproduct are the same fact seen twice.

**Fifth independent derivation of the same phase boundary**, now the cleanest: R6 (`(q−1)/2 = 1 ⟺ q=3`) · R7 (geometric series `Σ(q/3)^j` diverges iff q≤3) · R7 (`X_k` linear growth at q=3) · R8 (domination fails by a factor of k at q=3) · **R12 (`D₂ = 1 ⟺ q = 3`)**.

## 4. What to do with this

- **Cite Tao for the architecture** (Remark 1.15's address count) and **state explicitly that Prop 1.17 does not apply to the ℓ² sum** — a reader will otherwise assume it does, and the assumption is wrong by a factor of 3ⁿ.
- **`D₂ = 1 ⟺ q = 3` belongs in the abstract**, not buried. It is the paper's cleanest sentence and it lands the Bernoulli-convolution template (Erdős/Solomyak/Hochman: one-parameter family, phase transition identified by L² methods) as *the same kind of statement*, not an analogy.
- **`H₁` vs `H₂` is the honest positioning against Tao**: same object, adjacent question, different Rényi order. Not a duplicate, not a competitor.
- **The `exp(−cm)` remark** (Tao's own, unproven) is the regime he declined to enter and the one our exponential results live in. Worth a sentence.

## 5. Not at stake / still open

Prop 1.17 is not a route to Result 1's domination bound; the open step remains the **family-(b) collision count** (R11). Nothing here changes R5's rate, R6, R7, R8, or R10's law. THEOREM_C_745, Th 78.1–78.3, R81b, ε_k untouched.

**Still to read before publishing:** Siegel, *(p,q)-adic Analysis and Collatz* (USC 2022 / Springer 2024–25) — closest sibling. And the **p-adic Bernoulli convolution** literature, never covered by the 2026-05-04 dive (which looked only at archimedean BCs and Siegel).

_Reporting discipline: the read was routed BEFORE further derivation precisely to avoid re-deriving Tao, and the honest outcome is "he did the adjacent problem, not ours" — reported as such rather than inflated into either a duplicate-scare or a free lemma. The claim that Prop 1.17 is insufficient is quantified (summing gives 3ⁿ·n^{−2A}), not asserted. The `H₁`/`H₂` values are exact, not fitted._
