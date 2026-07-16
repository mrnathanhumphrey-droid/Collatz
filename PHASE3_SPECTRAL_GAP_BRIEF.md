# Phase 3 — Proof Brief: Result 1 as a spectral gap of a transfer operator

**Purpose.** A working brief for the rigorous proof of Result 1 (the domination bound / universal rate). It states the theorem, its reduction to a spectral gap, the operator, why the non-archimedean setting should make the gap *elementary* to prove, and a concrete lemma skeleton. Everything marked **[NUM]** is confirmed numerically (R22–R27, gate-validated); everything marked **[TP]** (to prove) is the pen-and-paper work. Written 2026-07-16 after R22–R28.

---

## 1. The theorem

**Result 1.** For odd `q ≥ 5`, `‖π_k‖² ~ C_q · 3^{-k}` as `k → ∞`, with `C_q ≥ 1`. Equivalently `X_k := q^k‖π_k‖² ~ C_q·(q/3)^k`. At the critical `q = 3` the normalized object grows **linearly** (`~ (7/15)·k`) rather than converging — the phase boundary.

Here `π_k` is the stationary law of the qx+1 Syracuse chain on `(Z/q^k)^*`, and (R8) `π_k` is the mod-`q^k` image of the **q-adic self-similar measure** of the IFS
```
    T_v(x) = (q·x + 1)·2^{-v},   v ≥ 1,   weight p_v = 2^{-v}     (Σ p_v = 1).
```
Every `T_v` contracts by **exactly** `1/q` in the q-adic metric (the `q·x` factor; `2^{-v}` is a q-adic unit).

## 2. Reduction to a spectral gap

By the address/overlap decomposition (R8),
```
    ‖π_k‖² = DIAG_k + OVERLAP_k,     DIAG_k = (Σ_v p_v²)^k = (1/3)^k,
```
since `Σ_{v≥1} 2^{-2v} = 1/3`. Write the normalized excess
```
    cross(k) := ‖π_k‖² / (1/3)^k − 1 − within(k),      P2 := Σ p_v² = 1/3  (limit).
```
`within(k)` is the within-cell overlap (closed form, k-flat for k≥3, R15). Result 1 ⟺ **`cross(k)` is bounded in k** for `q ≥ 5` (and linear for `q = 3`).

**The operator identity.** There is a transfer operator `L` (§3) with `cross(k) = Σ_i A_i μ_i^k`, where `μ_i = λ_i / λ_1` are the eigenvalue ratios of `L`, `λ_1 = 1/3` the Perron eigenvalue **[NUM: R25, λ₁=1/3 at every q to 5 digits]**, and `A_i` the modal amplitudes. Therefore:

> **`cross(k)` bounded ⟺ `|μ_i| < 1` for every subdominant mode with `A_i ≠ 0` ⟺ `L` has a SPECTRAL GAP above `λ_1`.**
>
> Define **`r_q := max{ |λ_i|/λ_1 : i ≥ 2, A_i ≠ 0 }`** = the spectral-gap ratio. Then:
> - `r_q < 1` ⟺ `cross` converges ⟺ Result 1 holds (q ≥ 5).
> - `r_q = 1` ⟺ gap closes ⟺ linear growth (q = 3).

**★ This kills R21's "wrong axis" objection.** R21 worried Konyagin's subgroup-sum bounds are in `q` while we need uniformity in `k`. But once `cross(k) = Σ A_i μ_i^k`, **a fixed-q gap gives k-uniform geometric decay for free** — the operator carries `k` via its powers. We only ever need a **fixed-q strict-contraction** statement, which is exactly what subgroup exponential-sum machinery delivers.

**Numerical status of the gap [NUM]:** `r_3 = 1` (degenerate Perron, R25), `r_5 ≈ 0.62`, `r_7 ≈ 0.38` (R26 amplitude + R27 direct high-k, gate-validated). The gap is *observed*; the theorem is to *prove* it.

## 3. The operator (two equivalent presentations)

**(a) Arithmetic / cascade — `M` (R25, the one that is gate-validated).** State `s = (a, b, γ)` with `a = 2^{-S} mod q^L`, `b = 2^{-S'} mod q^L` in `H = ⟨2⟩`, carry `γ ∈ Z/q^L`. One step prepends a coordinate to a pair of independent addresses, keeps the branch iff the q-adic digit condition `(γ + (a'−b')) ≡ 0 (mod q)` holds, updates the carry. Then `sum(M^k v_0) = ‖π_k‖²` **exactly for k ≤ L [NUM: gate to machine precision]**. `M` is nonnegative, sub-stochastic; Perron `λ_1 = 1/3`. Good for computation and validation; depth-`L` truncated (→ genuinely infinite-dimensional as `L → ∞`, because collisions are tower-generated, R24).

**(b) Analytic / Ruelle–Perron–Frobenius — `L`.** The transfer operator of the IFS acting on the **pair correlation**. On Fourier data `F_k(j) = μ̂(j/q^k)`, self-similarity gives `μ̂(ξ) = Σ_v 2^{-v} e(ξ 2^{-v}) μ̂(q 2^{-v} ξ)`, and `‖π_k‖² = q^{-k} Σ_{j mod q^k} |F_k(j)|²`. Forming `|F|²` couples two arguments `2^{-v}j, 2^{-v'}j` — the **same `(a,b)` pair structure as `M`**. `L` acts on locally-constant functions on `(Z_q^*)²`; it is a **finite matrix on each q-adic level** and preserves the level filtration. Good for the proof.

`M` and `L` are the same operator in two bases. Use `L` for the gap argument, `M` to check every claim numerically.

## 4. Why the non-archimedean setting should make the gap ELEMENTARY

This is the central optimism of the brief. Real self-similar-measure spectral gaps need heavy thermodynamic formalism (bounded distortion, transversality, cone conditions). **None of that is needed here:**

- **Exact contraction `1/q`** — no distortion, no Hölder cocycle. The transfer operator has *constant* weights on cylinders.
- **Ultrametric: balls are nested-or-disjoint** — the cylinder partition is exact and clopen; `L` is literally block-diagonalizable along the level filtration. Function-space truncation (§3b) is **exact and finite-dimensional per level**, not an approximation.
- **Overlaps are exact algebraic coincidences** (`2^{-S} ≡ 2^{-S'} mod q^k`), not geometric near-collisions — no transversality/Erdős machinery.
- **Perron datum is explicit** — the diagonal eigenvector is known (R8), `λ_1 = 1/3 = Σ p_v²` in closed form.

So the gap should reduce to a **finite linear-algebra + one uniform subgroup-sum bound**, per fixed `q`. That is the whole content of Result 1's remaining rigor.

## 5. The phase boundary — why q = 3 is exactly the non-gap

`d := ord_q(2)`, `H = ⟨2⟩ ⊆ F_q^*`, `|H| = d`. At `q = 3`: `d = 2`, `H = {1,2} = F_3^*` is the **whole** multiplicative group, and `(q−1)/2 = 1` (R6). The off-diagonal operator then **resonates**: its leading mode is not strictly contracted, `λ_2 = λ_1 = 1/3`, Perron degenerate (Jordan block) → `cross(k)` linear. For `q ≥ 5`, `H` is a **proper** subgroup (`d ≤ q−1` with genuine dispersion), the off-diagonal modes disperse, and `λ_2 < λ_1` strictly. **The theorem is: `H ⊊ F_q^*` (equivalently `q ≥ 5`) ⟹ strict contraction of `L` off the diagonal.** This is the register in which the Konyagin subgroup sums live (R21).

## 6. Lemma skeleton (the pen-and-paper program)

- **L1 (well-posedness) [TP, elementary].** `L` is bounded on locally-constant functions on `(Z_q^*)²` and preserves the finite q-adic level filtration; on each level it is a finite nonnegative matrix. *(From exact self-similarity + clopen cylinders.)*

- **L2 (Perron) [TP, ~R8 + Perron–Frobenius].** `λ_1 = 1/3` is a simple, strictly dominant eigenvalue of `L` on the **diagonal** subspace, with the explicit diagonal eigenvector. `C_q ≥ 1` is forced (Cauchy–Schwarz, R8). *(Already essentially in hand.)*

- **L3 (THE GAP) [TP, the crux].** On the **off-diagonal** (`a ≠ b`) subspace, `spec-radius(L_off) < 1/3` for `q ≥ 5`. Route: bound `‖L_off^n‖` by the collision / subgroup-incidence sum (R18/R20: condition `2^{-A'} ≡ 2^{-A} + j·s·2^{-A''} mod q`, `a,b ∈ H`) = a **Konyagin subgroup exponential sum**; show it is **strictly** `< 1` when `H ⊊ F_q^*`. *Note:* we need **strictness**, not the sharp value — so even the "vacuous-at-small-d" Konyagin bounds (R21's worry) may suffice, because strictness is a fixed-q qualitative statement.

- **L4 (assembly) [TP, routine once L2+L3].** `cross(k) = ⟨𝟙, L^k v_0⟩`-type `= A_1 + O(r_q^k)` with `r_q = 3·spec-radius(L_off) < 1` for `q ≥ 5` → `cross` bounded → Result 1. At `q = 3`, L3 fails (`H` full) → `r_3 = 1` → linear. **k-uniformity is automatic** (operator powers).

**The whole theorem is L3.** L1, L2, L4 are bookkeeping. L3 is "a fixed multiplicative subgroup `H ⊊ F_q^*` cannot additively self-align enough to preserve the leading off-diagonal mode" — a subgroup-sum non-degeneracy statement.

## 7. What `r_q` is (and the open sub-question)

`r_q = 3·λ_2(L_off)` is an **algebraic number** (root of a fixed characteristic polynomial per `q`). `r_3 = 1`, `r_5 ≈ 0.62`, `r_7 ≈ 0.38` **[NUM]**. **No elementary closed form** (R28: `3/q` refuted — `r_5 > 3/5` but `r_7 < 3/7`). The closed form is a *separate, compute-limited* question (pin `r_11, r_13` via Lambda-scale high-k or higher-`L` modal); it is **not needed for Result 1**, which only asserts `r_q < 1`.

## 8. Literature to plug in (→ the "gather" task, option a)

For **L3** specifically:
- Konyagin, *Exponential sums over multiplicative groups in fields of prime order* (`Bourgain-Konyagin/Konyagin_Lectures.pdf`, on disk) — our object is his §1.2 verbatim (R21). Need the **strict-positivity / non-degeneracy** form, not the asymptotic bound.
- Bourgain–Chang, `Bourgain-Konyagin/122 NewExp.pdf` (on disk) — composite moduli `q^k` + Heilbronn sums (Stepanov method), if the `F_q` reduction needs lifting.
For the **RPF framework** (L1, L2, L4) in the ultrametric setting:
- p-adic / non-archimedean self-similar (Hutchinson) measures and their transfer operators — **flagged as never covered** in project memory; the target of the gather.
- Siegel 2024 p,q-adic Collatz (`references/Q-sweep/Siegel2024_pq_adic_Collatz_consolidated.pdf`, on disk) — read before publication.

---

*Discipline note carried from the arc: structural claims here are gate-validated [NUM]; the value `r_q` has no committed closed form (quantitative priors 0-for-8 this arc). The brief asserts only the STRUCTURE (gap ⟺ Result 1, `H ⊊ F_q^*` ⟹ gap), which is what the proof needs.*
