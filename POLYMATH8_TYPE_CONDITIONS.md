# POLYMATH8 TYPE I/II/III — Verbatim Conditions and Essentiality

**Purpose.** Catalog the exact hypotheses of the Polymath8a Type I, Type II, Type III estimates. Classify each as ESSENTIAL (the proof uses it directly and breaks without it) or WORKAROUND (the proof goes through under a weaker variant).

**Primary source.** Polymath8a paper, arXiv:1402.0811, Algebra & Number Theory 8-9 (2014). For verbatim text, two informants in the chain: (i) Tao blog "Estimation of the Type I and Type II sums" (2013-06-12), (ii) Tao blog "Estimation of the Type III sums" (2013-06-14), (iii) Tao blog "An improved Type I estimate" (2013-07-27). All three closely mirror the paper's definitions.

---

## 1. Verbatim definitions (Type I/II post, restated unchanged in the paper)

### Definition 1 (coefficient sequence)

> "A coefficient sequence is a finitely supported sequence α: ℕ → ℝ that obeys the bounds |α(n)| ≪ τ^O(1)(n) log^O(1)(x) for all n, where τ is the divisor function."

**Subconditions:**

- **(i) Scale.** α is at scale N if it is supported on the interval `[(1 − O(log^{-A₀} x))N, (1 + O(log^{-A₀} x))N]`.
- **(ii) Discrepancy.** `Δ(α; a (q))` measures deviation from equidistribution across primitive residue classes.
- **(iii) Siegel-Walfisz property.** α obeys Siegel-Walfisz if `|Δ(α 1_{(·,q)=1}; a (r))| ≪ τ(qr)^O(1) N log^{-A} x` for any `q, r ≥ 1`, any fixed `A`, and any primitive residue class `a (r)`.
- **(iv) Smoothness.** α at scale N is **smooth** if `α(n) = ψ(n/N)` for some smooth `ψ: ℝ → ℂ` supported on `[1 − O(log^{-A₀} x), 1 + O(log^{-A₀} x)]` with `ψ^{(j)}(t) = O(log^{jA₀} x)` for all fixed `j ≥ 0`.

> **Smoothness and Siegel-Walfisz are NOT both required of all sequences.** They are alternative auxiliary properties any individual sequence may satisfy. Different theorems use one or the other.

### Definition 2 (singleton congruence class system on `S_I`)

> "**S_I denotes the square-free numbers whose prime factors lie in I.**"
> "A singleton congruence class system on I is a collection `C = ({a_q})_{q ∈ S_I}` of primitive residue classes `a_q ∈ (ℤ/qℤ)×` for each `q ∈ S_I`, obeying the Chinese remainder theorem property `a_{qr} (qr) = (a_q (q)) ∩ (a_r (r))` whenever `q, r ∈ S_I` are coprime."

> **The class is square-free by definition.** Prime-power moduli `p^r` with `r ≥ 2` are excluded a priori.

### Definition (dense divisibility — Improved Type I post, post 8)

> "A natural number n is k-tuply y-densely divisible if for any i, j ≥ 0 with `i + j = k − 1` and `1 ≤ R ≤ n`, one can factor `n = qr` with `y^{-1} R ≤ r ≤ R` such that q is i-tuply y-densely divisible and r is j-tuply y-densely divisible."

> Dense divisibility GENERALIZES Zhang's y-smoothness (which required every prime factor ≤ y) to a factorization-at-every-scale condition. But the class remains within `S_I`, hence square-free.

---

## 2. Type I/II estimate (Theorem 3 of the 2013-06-12 post; matches Polymath8a paper Theorem)

### Verbatim hypothesis

> "Let ϖ, δ, σ > 0 be fixed quantities such that
> `11ϖ + 3δ + 2σ < 1/4`  and  `29ϖ + 5δ < 1/4`
> and let α, β be coefficient sequences at scales M, N respectively with `x ≪ MN ≪ x` and `x^{1/2−σ} ≪ N ≪ M ≪ x^{1/2+σ}` with β obeying a Siegel-Walfisz theorem. Then for any `I ⊂ [1, x^δ]` and any singleton congruence class system with controlled multiplicity we have
> `Σ_{q ∈ S_I: q < x^{1/2+2ϖ}} |Δ(α ∗ β; a_q)| ≪ x log^{-A} x`."

### Improved Type I (Theorem 4 of the 2013-07-27 post)

> "We have `Type^{(4)}_I[ϖ, δ, σ]` whenever `(160/3)ϖ + 16δ + (34/9)σ < 1` and `64ϖ + 18δ + 2σ < 1`."

> Hypothesis: α, β coefficient sequences at scales M, N with β obeying Siegel-Walfisz; `q ∈ S_I ∩ D^{(k)}_{x^δ}` (k-tuply x^δ-densely-divisible AND in `S_I`).

### Hypothesis catalog and essentiality

| # | Hypothesis | Where in proof it is used | Verdict |
|---|------------|---------------------------|---------|
| H1 | α is a coefficient sequence (polynomial divisor-bounded) | Standard archimedean control; used everywhere in completion-of-sums and divisor-function moment bounds | ESSENTIAL |
| H2 | β obeys Siegel-Walfisz | Used to handle the "main term" / Bombieri-Vinogradov-type input for arithmetic progressions to small moduli, after the bilinear Cauchy-Schwarz reduction | ESSENTIAL |
| H3 | α need NOT be smooth, NOR Siegel-Walfisz | α enters the bilinear sum on the un-completed side; the Cauchy-Schwarz step squares it out into a divisor-function moment | WORKAROUND (rough α is fine) |
| H4 | β need NOT be smooth (only Siegel-Walfisz) | The β-completion step uses Siegel-Walfisz directly to control AP residues | WORKAROUND for smoothness; but Siegel-Walfisz itself is ESSENTIAL |
| H5 | `q ∈ S_I` (square-free) | Used in: (a) the Chinese remainder factorization of `(ℤ/qℤ)×`, (b) the multiplicativity of Ramanujan/Kloosterman sums, (c) the trace function machinery (Fouvry-Kowalski-Michel), (d) avoiding ramification in the Weil-bound algebraic geometry | **ESSENTIAL** — every downstream tool in 8a's proof assumes it |
| H6 | `q` is x^δ-densely-divisible (or k-tuply) | Used for the multi-scale factorization `q = q_1 q_2 ... q_l` that enables van der Corput / q-van der Corput iteration | ESSENTIAL **as factorization tool** — but the WORKAROUND of dense divisibility GENERALIZES Zhang's y-smoothness, showing it was originally too restrictive |
| H7 | Scale balance `x^{1/2−σ} ≪ N ≪ M ≪ x^{1/2+σ}` | Required for bilinear Cauchy-Schwarz to give nontrivial gain | ESSENTIAL |
| H8 | Range `q < x^{1/2+2ϖ}` | The exponent of distribution we're proving | OUTPUT, not hypothesis |
| H9 | Singleton congruence class system with controlled multiplicity (CRT-compatible) | Required for the Bombieri-style averaging over residue classes | ESSENTIAL |

---

## 3. Type III estimate (Theorem 2 of 2013-06-14 post)

### Verbatim hypothesis (paraphrased from the post — the post lays out the same conditions used in Polymath8a)

> Let M, N₁, N₂, N₃ ≫ 1 with `x ≪ M N₁ N₂ N₃ ≪ x`, `N₁ ≫ N₂, N₃`, and `N₁^4 N₂^4 N₃^5 ≫ x^{4 + 16ϖ + δ + c}`. Let α be a coefficient sequence at scale M; let ψ₁, ψ₂, ψ₃ be **smooth** coefficient sequences at scales N₁, N₂, N₃. Let `q ∈ S_I` densely-divisible. Then the appropriate trilinear sum is bounded by `x log^{-A} x`.

### Hypothesis catalog and essentiality

| # | Hypothesis | Use | Verdict |
|---|------------|-----|---------|
| T1 | α at scale M, polynomial-bounded only | enters as the "outer" sequence; gets squared via Cauchy-Schwarz | ESSENTIAL |
| T2 | **ψ₁, ψ₂ smooth** | Completion-of-sums via Poisson on these two variables — requires Schwartz-class control of derivatives | **ESSENTIAL — smoothness is load-bearing for these two** |
| T3 | ψ₃ smooth | The third smoothness was also assumed in the paper; reading at the time noted "only two of three need smoothness" because the third is eliminated by Cauchy-Schwarz — but as STATED the theorem assumes all three smooth | ESSENTIAL as stated; **WORKAROUND** in principle if one Cauchy-Schwarzes ψ₃ out first |
| T4 | `q ∈ S_I` square-free, densely-divisible | Same Weil / Chinese-remainder machinery as Type I/II; here driving Birch-Bombieri 3-variable exponential-sum bounds on a smooth variety | **ESSENTIAL** — Deligne's Weil II for the surface only works for proper smooth varieties; ramification (i.e., `p^r` with `r ≥ 2`) breaks the algebraic geometry |
| T5 | Bombieri-Birch / Deligne | Source of √-savings on 3-d exponential sums | ESSENTIAL — the proof IS a Weil-bound proof |
| T6 | Singleton congruence class system | CRT-decompose | ESSENTIAL |

---

## 4. Where each ESSENTIAL condition fails for the Collatz target

The user's target is a bilinear character sum with:

- Modulus `q = p^{r+1}`, prime-power. Fails H5/T4.
- Rough amplitude `1̂(ξ) = Σ_{u=0}^{N-1} e_q(ξ u)` (Dirichlet kernel, no smooth cutoff). Fails T2/T3 if used in a Type III role; fails H4's spirit if smoothness is the only completion tool — but Siegel-Walfisz is not the right substitute either, because the target is a fixed-prime-power Fourier sum, not an AP residue.
- Empirical saving is `√q` from a coset cardinality. Replaces neither the Weil bound (which would also give `√q` here but is unavailable due to ramification) nor an amplifier.

The three failures align with the three relaxations the user is looking for.

---

## 5. Summary table

| Condition | Type I/II/III treatment | Essentiality | Relaxable in chain? |
|-----------|-------------------------|--------------|---------------------|
| Polynomial-bounded coefficient (divisor-fcn growth) | Required for all sequences | ESSENTIAL | No |
| Smoothness | Required for: Type III amplitudes ψᵢ; the *β* of Type I/II only needs Siegel-Walfisz | ESSENTIAL where used | **Partial WORKAROUND: Siegel-Walfisz substitutes for smoothness in Type I/II β-role only** |
| Siegel-Walfisz | Required for β of Type I/II | ESSENTIAL | No |
| Square-free modulus (`S_I`) | Required everywhere | ESSENTIAL | **No.** Universal across the chain. |
| Densely-divisible modulus | Required where iterated factorization is used | ESSENTIAL **as factorization tool**, but is itself a relaxation of Zhang's y-smoothness | WORKAROUND for y-smoothness; not for square-free |
| Weil bound / Deligne | Source of √-savings in every Type I/II/III proof | ESSENTIAL | **No** parallel construction in the chain |
| Cauchy-Schwarz halving | Bilinear-to-quadratic reduction at every step | ESSENTIAL | No |
| Singleton CRT congruence system | Required for the averaging | ESSENTIAL | No |
