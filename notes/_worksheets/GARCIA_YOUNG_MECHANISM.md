# Garcia–Young 2023: where the √q secondary main term `A` actually arises

**Date:** 2026-05-11. Companion to `GARCIA_YOUNG_TRANSLATION.md`. Mechanism trace from the paper's Theorem 1.2 / Theorem 1.1 proofs.

The brief asks: their secondary main term `A` (resp. `A'`) is of rough size `q^{1/2}` and is the most suggestive structural parallel to our empirical `|K| / √N ≈ 2.0` saturation. Does it come from a bilinear-sum mechanism that translates to our setting?

## 1. The chain of sums

In §3 the paper sets up (eq 2.13, Lemma 2.18):
> `M = Σ_{χ even mod d} |L(1/2, χ·ψ)|² = ϕ(d) · Σ_± Σ_{m≡±n mod d, (mn,q)=1} ψ(m)ψ̄(n)/√(mn) · V(mn/q)`.

This is the AFE (approximate functional equation) opening: the second moment along the coset is converted to a divisor-style sum over pairs `(m,n)` with `m ≡ ±n (mod d)`.

§3.1 isolates the diagonal `m = n` which gives `D` (eq 1.7) via Lemma 2.19. **This is not the source of `A`.**

§3.2 (eq 3.6) writes the off-diagonal sum with `m = ±n + dl`, `l ≥ 1`:
> `B^±_{m>n}(M, N) = Σ_{l ≥ 1} Σ_{n ≥ 1} ψ(±n + dl) ψ̄(n) W_{M,N}(±n+dl, n)`.

The `n = 1` term of this bilinear sum (i.e., the `n=1, m=1+dl` slice) is the most-unbalanced range identified in the heuristic eq (1.9) and the **rigorous** entry point for `A`.

## 2. Poisson on `l` (§3.3) — this is where Lemma 2.15 enters

§3.3 applies **Poisson summation in the variable `l`** to `B^±_{m>n}(M, N)`. The argument is in the proof of Lemma 3.3 (paper p. 11):

After unfolding the Poisson dual and reorganizing, the paper arrives at (just before eq 3.13):
> `B^±_{m>n}(M, N) = (d/q) · Σ_{n≥1} Σ_{j∈Z} Ŵ^±_{n,d}(±dj/q) · S_{q,d}(ψ, jn)`,
where `S_{q,d}(ψ, jn)` is exactly the GY auxiliary sum defined in eq (2.10).

**Lemma 2.15 is applied here**, evaluating `S_{q,d}(ψ, jn)`. The result is eq (3.13):
> `B^±_{m>n}(M, N) = (-2a_ψ̄/(q/d²)) · ε_q · (d/√q) · Σ_{n≥1} Σ_{j≠0, jn≡-a_ψ mod d} Ŵ^±_{n,d}(±dj/q) · e_q(2a_ψ̄ (jn + a_ψ)²)`.

So the **`√q` enters via Lemma 2.15's evaluation of `S_{q,d}(ψ, jn)`**: GY's `S_{q,d}` has magnitude exactly `√q` on the (single) supporting residue class `k ≡ -a_ψ (mod d)`, and the Poisson chain leaves a `d/√q` prefactor outside the `(j,n)` sum. The combined size of the `j,n`-sum is `d/√q × O(1)` (because only one residue class contributes, plus the `Ŵ`-rapid-decay restricts `n | b_ψ` to a divisor sum), giving the **net size `(d/√q) · q/d² · √q · ... = ... ϕ(d)/d · √q × σ_0(|b_ψ|)/√|b_ψ|`** after one more reduction (§3.5 Lemma 3.6, paper p. 13–14).

The end-state (Lemma 3.7) is:
> `A_{m>n}(ψ) = (-2a_ψ̄/(q/d²)) · ε_q · ϕ(d)/2 · e_q(2a_ψ̄ (a_ψ − b_ψ)²) · σ_0(|b_ψ|) / √|b_ψ|`.

After symmetric addition with `A_{m<n}` (§3.6), this becomes the full `A'` in eq (1.11), of size roughly `ϕ(d) · σ_0(|b_ψ|)/√|b_ψ|` — but with the `1/d · √q` factor from Lemma 2.15's evaluation absorbed, leaving the apparent `q^{1/2}` only implicit in the heuristic (1.9) form `ϕ(d) · q^{1/2}/(d√|a_ψ|)` which we recover by reading `σ_0(|b_ψ|)/√|b_ψ| · ϕ(d)` and noting `|b_ψ| ≤ d/2 ≤ √q/2`.

## 3. So where is the √q **really** coming from?

There are **two** distinct √q contributions that ultimately compose:

**(a) Gauss-sum magnitude of one `S_{q,d}(ψ, k)`.** Lemma 2.15 says `|S_{q,d}(ψ, k)| = √q` on its single non-vanishing residue class. This is **not** a bilinear cancellation: it is the magnitude of a single complete sum of length `q/d` whose phase is quadratic in `u` (after Lemma 2.13(2) truncation). Standard quadratic Gauss-sum identity (Lemma 2.4) gives this. **Size: √q. Mechanism: Gauss sum.**

**(b) The Plancherel-side prefactor `d/√q`.** This comes from the Poisson summation step (paper line just before eq 3.13: `(d/q)·...·S_{q,d}`, and `|S_{q,d}| = √q` gives the `d/√q` net). **Origin: Plancherel duality on Z/q-frequencies.**

Combined: `(d/√q) · √q = d`. The √q's **cancel exactly** in the prefactor; the residual size of `A` is then determined by the *number of frequencies contributing* (the divisor count `σ_0(|b_ψ|)`) and the *amplitude* `1/√|b_ψ|` from the Fourier transform of the V-weight evaluated at the stationary frequency (Lemma 2.5 / 3.6).

**So in the GY proof, the √q is NOT a bilinear-cancellation √q.** It is:
- a Gauss-sum magnitude `√q` (from one Postnikov sum), and
- a Plancherel-duality prefactor `1/√q` (from Poisson),
- which combine to a *deterministic* factor that scales the size of the secondary term `A`.

The **size** of `A` (as `≪ √q`) appears because:
- the diagonal main `D` is `ϕ(d)² ϕ(q)/q · (log q)` (size `ϕ(d)² ϕ(q)/q · log q ≈ d · log q` since `ϕ(d) ≈ d`),
- the residual main `A` after the √q cancellation is `ϕ(d) · σ_0(|b_ψ|)/√|b_ψ|`, with `|b_ψ| ≤ d/2`,
- so `|A| ≤ ϕ(d) · σ_0(|b_ψ|) · √(2/d)` ≪ `ϕ(d) · q^ε / √d` ≤ `d · q^ε / √d` = `√d · q^ε` ≤ `√q · q^ε`.

The bound `|A| ≪ √q` is just `|b_ψ| ≥ 1` plus the trivial divisor bound on `σ_0`. **The √q is a side-effect of the parameter regime `q ≼ d²` (and `q ≼ d³`), not the result of any bilinear cancellation argument.**

## 4. Lemma 2.15 input/output as standalone

| Aspect | Lemma 2.15 |
|---|---|
| **Input** | A pair `(ψ, k)` where `ψ` is primitive mod `q` and `k ∈ Z` |
| **Hypothesis** | `(q, 3) = 1`, `d² | q`, `q | d³` |
| **Object evaluated** | `S_{q,d}(ψ, k) = Σ_{u (mod q/d)} ψ(1+du) e_q(dku)` |
| **Output (case `k ≡ -a_ψ mod d`)** | `ε_q · √q · e_q(2a_ψ̄(k+a_ψ)²) · (Jacobi(-2a_ψ̄ / (q/d²)))` |
| **Output (case `k ≢ -a_ψ mod d`)** | `0` |
| **Mechanism inside the proof** | Use Postnikov (Lemma 2.14) to rewrite `ψ(1+du) = e_q(a_ψ · L_q(1+du))`; use Lemma 2.13(2) to truncate `L_q(1+du) ≡ du - 2̄(du)² mod q`; **change variables `u ↦ u + q/d²`** to detect the periodicity, reducing to a sum mod `q/d²`; apply the quadratic Gauss-sum evaluation Lemma 2.4. |

**Key features for our purposes:**
1. The output is **size √q** uniformly on its single supporting residue class.
2. The phase `e_q(2a_ψ̄(k+a_ψ)²)` is **quadratic in `k`**, with a Gauss-sum-style normalization `ε_q · (Jacobi)`.
3. The vanishing condition `k ≢ -a_ψ (mod d)` is a **support condition** on the bilinear `(j,n)` lattice: only `(j,n)` with `jn ≡ -a_ψ (mod d)` contribute, restricting the bilinear range from `(d/√q) · (q^{1+ε}/M) × N` to a divisor-sum-style restriction.
4. **Lemma 2.15 is NOT itself a bilinear bound.** It is a closed-form evaluation of a complete sum of length `q/d`. The bilinear input is what is summed *outside* of Lemma 2.15 (the `(j,n)` indices), and the bilinear bound on that residual sum is supplied by the Ŵ-rapid-decay (eq 3.11) plus the divisor restriction.

## 5. Where the "Burgess-shape" actually lives in the paper

The closest object in GY to a standalone bilinear bound is Lemma 2.2 (= Heath-Brown 1978 [HB1] Lemma 9), eq (2.2):
> `Σ_{|h|≤A} Σ_{|n|≤B} |S(q; χ, hq_0, n)| ≪ q^{1/2+ε} · (AB · q_0^{-1/2} + (qAq_0)^{1/4})`
where `S(q; χ, h, n) = Σ_m χ(m+h)χ̄(m) e(mn/q)` is a Kloosterman-shape character sum (eq 2.1).

**This** is a bilinear bound (sum over `(h,n)` of a hyper-Kloosterman-style inner sum). It is used in §3.4 (Lemma 3.4) to handle the **balanced range** `M ≈ N` of the off-diagonal bilinear, via Poisson in `n` rather than in `l`.

So GY's proof has two bilinear inputs:
- **Unbalanced range (M ≫ N):** Poisson in `l`, evaluated via Lemma 2.15 (the Postnikov-truncated-log Gauss-sum lemma) — this yields the secondary main term `A`.
- **Balanced range (M ≈ N):** Poisson in `n`, evaluated via Lemma 2.2 (Heath-Brown's hybrid bilinear bound for `S(q;χ,h,n)`) — this is treated purely as an error term and produces no main term.

The √q in `A` emerges from the **unbalanced** range. The "Burgess-shape" bilinear bound from Heath-Brown lives in the **balanced** range and provides only an upper-bound error term.

## 6. Verdict on the mechanism

The √q in GY's secondary main term `A` comes from:
- **Gauss-sum magnitude of one Postnikov-coset complete sum** (Lemma 2.15's `√q`), arising at depth `d² | q | d³` where the Postnikov log truncates to quadratic in `Z/q[u]`.

It does **NOT** come from:
- a standalone bilinear bound of the form `|Σ x_a · y_a| ≤ C · √(Σ|x|²) · √(Σ|y|²)` with `Σ|x|² ≪ N`.
- a Burgess-shifted-pairs argument.
- a stationary-phase cancellation in the bilinear range.

The bilinear range `(j, n)` is handled by **rapid decay of the V-weight Fourier transform** (Ŵ falls off like `(|x|M/d)^{-C}` for any C, eq 3.11), which restricts the effective range of `(j, n)` to a divisor-of-`b_ψ` set — **not** by cancellation among bilinear factors. There is no √N saving on the principal-unit coset analogous to what we want.

**Pre-registration outcome:** This rules out **H_EXTRACTION_WORKS**. Lemma 2.15 by itself is not "the bilinear bound translated into our notation"; it is a Gauss-sum evaluation that GY use **inside** a larger argument to extract a main term, not a sub-Weyl bilinear estimate on the principal-unit coset.

The question now becomes whether **the Gauss-sum evaluation itself** (Lemma 2.15) is what we need (just for a different purpose than GY use it for) — i.e., whether it translates to a closed-form evaluation of our **`G[a] = F̂(p·a)/p`** that, plugged into our bilinear via Plancherel, gives the bound we want. This is the H_PARTIAL scenario, examined next in `GARCIA_YOUNG_TRANSLATION_ATTEMPT.md`.
