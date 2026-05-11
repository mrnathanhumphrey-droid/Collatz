# Garcia–Young 2023 → Burgess bundle: notation translation

**Date:** 2026-05-11. Reading + translation only. No new mathematics.

**Paper:** Bradford Garcia, Matthew P. Young, "Asymptotic second moment of Dirichlet L-functions along a thin coset." arXiv:2312.08482v1 (13 Dec 2023). Forum Math Sigma vol 13:e83 (DOI 10.1017/fms.2025.44). Lemma numbering checked against v1; published Forum Math Sigma version may differ — flagged where relevant.

## 1. Object correspondence

The brief asks whether Garcia–Young's Lemma 2.15 (or its analogue Lemma 2.16, the q ≼ d² regime relevant to Theorem 1.1) gives our bilinear bound after notation conversion. The first task is to lock down the dictionary between their setting and ours.

### Their setting (Garcia–Young, §§1–2)

- `q` (their) — a positive integer modulus.
- `d` — positive divisor of `q` with `q | d²` (Theorem 1.1) or `d² | q | d³` (Theorem 1.2). Throughout `d` and `q` share all prime factors (`a ≼ b` notation, paper p. 2 line 41).
- `ψ` — primitive even Dirichlet character mod `q`.
- `χ` — auxiliary character mod `d` summed over in the second moment.
- `L_q(1+dx)` — Postnikov logarithm, formal power series defined in Definition 2.6 eq (2.4):
  `L_q(1+dx) = Σ_{k≥1} (-1)^{k+1} (d^k / k) x^k`,
  viewed as an element of `(Z/qZ)[x]` after Lemma 2.10.
- `a_ψ ∈ Z/(q/d)Z` — the **Postnikov exponent** of `ψ`, defined via Lemma 2.14 eq (2.9): for all `x ∈ Z`, `ψ(1+dx) = e_q(a_ψ · L_q(1+dx))`.
- **The auxiliary sum (eq 2.10):**
  `S_{q,d}(ψ, k) := Σ_{u (mod q/d)} ψ(1+du) e_q(dku)`,    `k ∈ Z`.

### Our setting (Burgess bundle, locked in `PRECISE_ASK.md`)

- `p` — odd prime (≥ 3).
- `q = p^{r+1}` (their `q` ↔ our `q`).
- `period = p^r` (multiplicative order of `1+p` in `(Z/q)*`).
- `N = p^{r-1} = q/p² = period/p` (length of the short window).
- `c ∈ (Z/q)*` — the c-parameter of f; canonical case `c = 1`.
- `f(u) = e_q(c · (1+p)^u)`, periodic with period `p^r`.
- **Object F̂** (eq, brief §1): `F̂(ξ) = Σ_{u=0}^{q-1} f_periodic(u) · e_q(-ξu)`, supported on `{p·a : a ≡ 1 (mod p) in Z/p^r}`, magnitude `p·√q = p^{(r+3)/2}` (verified theorem, `FHAT_THEOREM_VERIFICATION_RESULTS.md`).
- **Object 1̂** (eq, brief §1): `1̂(ξ) = Σ_{u=0}^{N-1} e_q(ξu)` (length-N Dirichlet kernel).

### The translation

The mapping is forced once we identify what plays the role of `1+du` in our setting.

| Their object | Our analogue | Note |
|---|---|---|
| `q` | `q = p^{r+1}` | exact match — both are the L-function modulus |
| `d` | `p` | depth 1: their `d` divides their `q`, here `p` divides `p^{r+1}` |
| `q/d` | `p^r = period` | length of the coset / period |
| `1+du` (their parametrization of the principal-unit coset `K = ker((Z/qZ)* → (Z/dZ)*)`) | `(1+p)^u` (our parametrization of the **same** subgroup — for odd `p`, `K = ⟨1+p⟩`, see `FHAT_THEOREM_VERIFICATION_RESULTS.md` §3 Step 2) | Both parametrize the depth-`d` principal-unit coset, but **via different bijections**: GY use `u ↦ 1+du` (linear), we use `u ↦ (1+p)^u` (exponential). These differ — see §3 below. |
| `ψ(1+du)` | `e_q(c · (1+p)^u) = f(u)` | After Postnikov, with our `c` ↔ their `a_ψ`. **See §2 below for the precise statement.** |
| `a_ψ` | `c` | our c-parameter is the Postnikov exponent of the character `u ↦ e_q(c·(1+p)^u)` if/once we straighten the parametrization (§3 caveat) |
| `e_q(dku)` | (no analogue in our F̂) | their `dku` is a linear twist (from the off-diagonal sum); we have **no** corresponding `k`-parameter on the F̂ side |
| `S_{q,d}(ψ, k)` (their eq 2.10) | a length-`q/d` = length-`p^r` sum over the principal-unit coset of `ψ(1+du)·e_q(dku)` | this **IS** structurally our `G(a) = F̂(p·a)/p` (length-period DFT of f), with `k` ↔ `a` after the parametrization swap |

### Their `S_{q,d}(ψ, k)` is structurally our `G[a]`

Side-by-side after dictionary substitution:

```
GY 2.10:    S_{q,d}(ψ, k) = Σ_{u (mod q/d)} ψ(1+du) · e_q(dku)

Our G[a]:   G[a]          = Σ_{u (mod p^r)} f(u)    · e_{p^r}(-au)
                          = Σ_{u (mod p^r)} e_q(c·(1+p)^u) · e_{p^r}(-au)
                          = Σ_{u (mod p^r)} e_q(c·(1+p)^u) · e_q(-p·au)
                                                              (using q = p·p^r, so e_{p^r}(-au) = e_q(-p·a·u))
```

The length of summation matches (`q/d = p^r`). The principal-unit argument of `ψ` / `f` matches (`1+du` ↔ `(1+p)^u`, both parametrize `K = ⟨1+p⟩`). The linear twist matches in form: their `e_q(dku)` ↔ our `e_q(-p·a·u)` (with `d` ↔ `p` and `k` ↔ `-a`).

**Conclusion of §1.** Their `S_{q,d}(ψ, k)` and our `G[a] = F̂(p·a)/p` (or, equivalently, `F̂(p·a)` up to the factor of `p`) are **structurally the same object**, up to (i) the parametrization mismatch between `u ↦ 1+du` vs `u ↦ (1+p)^u`, and (ii) the role of `k`/`a` (their `k` is summed over after Poisson; our `a` is the index of the F̂ support point and is summed over in the bilinear sum). Lemma 2.15 evaluates this object explicitly under the hypotheses `d² | q | d³`. Lemma 2.16 evaluates it under `q | d²`. Mapping the hypotheses to our setting is the work of §4.

---

## 2. Postnikov correspondence

GY's Lemma 2.14 gives, for primitive `ψ` mod `q`:
> `ψ(1+dx) = e_q(a_ψ · L_q(1+dx))`,
where `L_q(1+dx) ≡ dx - 2̄(dx)² (mod q)` if `q | d³` (Lemma 2.13(2)), and `L_q(1+dx) ≡ dx (mod (q,d²))` if `q | d^∞` (Lemma 2.13(1)).

Specialization to `d = p`, `q = p^{r+1}`:
- Lemma 2.13(1) applies iff `q | d²`, i.e., `p^{r+1} | p²`, i.e., `r ≤ 1`. So at `r ≤ 1` (and only then), `L_q(1+px) ≡ px (mod q)`.
- Lemma 2.13(2) applies iff `q | d³`, i.e., `r ≤ 2` (and `(q,3) = 1` — which **fails** at `p = 3`).
- For `r ≥ 3`, **neither truncation lemma applies directly**. The L_q expansion at depth ≥ 4 contributes nontrivial higher-degree terms in `Z/qZ[x]`.

Compare to our `R78.4` factorization (`result_78_extended.md` Theorem 78.4):
> `F̂(3·a) = 3 · e_q(1) · G(a)`, `G(a) = Σ_{s=0}^{period-1} e_q(P_a(s))`, `P_a(s) = 3s − C_a · L(1+3s)`,
where `L` is the **Cochrane truncated 3-adic log** (truncated at some explicit r-dependent depth).

GY's `L_q` and Cochrane's truncated p-adic log are **literally the same power series** (eq 2.4 vs the Cochrane convention used in `result_78_extended.md`), reduced mod `q` in the same way (via Lemma 2.10). The cubic phase `P_a(s) = 3s − C_a · L(1+3s)` is therefore — by direct identification — the GY phase `dku + a_ψ · L_q(1+du)` with:
- `k` ↔ a parameter `k_a` to be determined from the change of variables,
- `a_ψ` ↔ `-C_a` (sign convention),
- `du` ↔ `3s`.

**The phase classes match.** This was already flagged in `milicevic_banks_verification.md` §"Milićević 2014" for the primal Milićević F-class. Garcia–Young use the same Postnikov framework, so the phase-class match carries over.

---

## 3. Parametrization mismatch: `1+pu` vs `(1+p)^u`

This is the load-bearing subtlety.

GY parametrize the kernel `K = ⟨1+p⟩` of `(Z/qZ)* → (Z/pZ)*` **additively**: every element of `K` is `1+pu (mod q)` for some `u (mod q/d) = u (mod p^r)`. This is a bijection because for odd `p`, `K` has order `p^r` and consists exactly of residues `≡ 1 (mod p)` in `(Z/qZ)*`.

We parametrize the same `K` **multiplicatively**: every element is `(1+p)^u (mod q)` for `u (mod p^r)`. This is also a bijection (`1+p` generates `K` cyclically).

These are two parametrizations of the same finite cyclic group of order `p^r`. They are related by the discrete logarithm: there is a bijection
> `φ : Z/p^r Z → Z/p^r Z`,    `1+pφ(u) ≡ (1+p)^u (mod q)`.

Equivalently, `φ(u) = (((1+p)^u − 1)/p) mod p^r`. By the binomial expansion mod p^{r+1}:
> `(1+p)^u = 1 + pu + C(u,2)·p² + ... = 1 + p·u + p²·(u(u-1)/2) + p³·(...)`,
so `φ(u) ≡ u + p·(u(u-1)/2) + p²·(...) (mod p^r)`. **The bijection `φ` is itself a Cochrane-log-style power series in `u`** — specifically, `φ(u) = u · (some unit power series in `pu`)`.

### Consequence

Substituting `u ↦ φ(u)` in GY's `S_{q,d}(ψ, k)` (which has `ψ(1+pu)` and `e_q(pku)`) converts:
- `ψ(1+p·φ(u)) → e_q(c · (1+p)^u)` (with `c = a_ψ`) — this **is** our `f(u)`.
- `e_q(p·k·φ(u)) → e_q(p·k·φ(u))` — but `φ(u)` is **nonlinear in `u`**, so this is **not** the same as `e_q(-p·a·u)` (our linear twist).

In words: the substitution `u ↦ φ(u)` turns GY's `ψ(1+du)`-side into our `f(u)`-side correctly (they're literally the same object), but the linear twist `e_q(dku)` becomes a **nonlinear** twist `e_q(d·k·φ(u))` in the multiplicative coordinate `u`.

This is the **first** point of friction. We address its size in §3 of `GARCIA_YOUNG_TRANSLATION_ATTEMPT.md`. Two scenarios:
- **Scenario A.** The nonlinearity of `φ` produces only "Hensel correction" perturbations to the saddle, leaving the leading evaluation invariant (modulo an explicit unit multiplier). Result: GY's Lemma 2.15 closed form carries over up to a determined factor.
- **Scenario B.** The nonlinearity of `φ` is essential and the saddle integral evaluates differently in the multiplicative coordinate. Result: their lemma evaluates a closely related but distinct object, and the closed form does not transfer.

Both scenarios are consistent with what's *visible* in the paper: GY never need to leave their `1+du` coordinate, so the multiplicative-coordinate evaluation is **not** done in their text.

---

## 4. Side-by-side heuristic eq (1.9) vs our bilinear

GY's heuristic for the size of `A` (paper §1, eq 1.9):
> `φ(d) · Σ_{l ≥ 1} e_{q/d}(a_ψ l) / √(dl) ≈ φ(d)/√d · Σ_{l < q/(d|a_ψ|)} 1/√l ≈ φ(d) · q^{1/2} / (d · √|a_ψ|)`.

Translation to our setting (`d ↔ p`, `q ↔ p^{r+1}`, `a_ψ ↔ c`, `φ(d) = p−1 ↔ ϕ(p) = p−1`):
- `q/d = p^r`,
- `q/(d|a_ψ|)` = `p^r / |c|` (depends on the unit `c`),
- the sum becomes `(p-1) · Σ_{l < p^r/|c|} e_{p^r}(c·l) / √(pl)`.

Our bilinear sum (`PRECISE_ASK.md` §3 Form B):
> `Σ_{a ≡ 1 (mod p) in Z/p^r} 1̂(p·a) · F̂(p·a)`.

These are **structurally different**:
- GY (1.9) sums **`e_{q/d}(a_ψ l) / √(dl)`** over `l` — a Dirichlet-series-like sum where the modulus is `q/d` and the amplitude `1/√(dl)` comes from the V-weight in the AFE (Lemma 2.17).
- Our bilinear sum has **`1̂(p·a) · F̂(p·a)`** — Plancherel-conjugate of a partial sum, weighted by the Dirichlet kernel `1̂(p·a) = Σ_{u=0}^{N-1} e_q(p·a·u)`.

The heuristic (1.9) is **already on the Plancherel-dual side** in a specific sense: `Σ_l e_{q/d}(a_ψ l) / √(dl)` is the Mellin transform along `Re(s) = 1/2` of `Σ_l e_{q/d}(a_ψ l) l^{-s}`, which is the Lerch transcendent `Φ(e_{q/d}(a_ψ), 1/2, ...)`. The `√q / √(d|a_ψ|)` size **emerges from the cancellation in the exponential character sum truncated at `q/(d|a_ψ|)`**, by partial summation against `1/√l`.

**This is not our bilinear sum.** Our `1̂(p·a)` and `F̂(p·a)` are both already evaluated dual-side objects; GY's heuristic is a primal-side L-function moment heuristic. The √q size in (1.9) does **not** come from a bilinear-sum bound; it comes from **truncating an oscillating geometric series at its first stationary point** `l ≈ q/(d|a_ψ|)`. Different mechanism.

What about the **rigorous** appearance of √q in GY? See `GARCIA_YOUNG_MECHANISM.md`. The short answer: it appears in (3.13) via Lemma 2.15, but Lemma 2.15's √q is `ε_q · √q` — a **Gauss-sum normalization factor**, attached to a **single saddle point**, not the result of cancellation across the bilinear range. Mechanism described in detail next.

---

## 5. Side-by-side: Lemma 2.15 vs our F̂ theorem

**Lemma 2.15 (GY, p. 8, verbatim).** Let `(q,3) = 1`, `ψ` have conductor `q`, and suppose `d² | q` and `q | d³`. Also let `k ∈ Z`. If `k ≢ −a_ψ (mod d)` then `S_{q,d}(ψ, k) = 0`. If `k ≡ −a_ψ (mod d)` then
> `S_{q,d}(ψ, k) = ε_q · √q · e_q(2a_ψ̄ (k + a_ψ)²) · (−2a_ψ̄ / (q/d²))`
where `2a_ψ̄` is the multiplicative inverse of `2a_ψ` mod `q`, the bracketed object is the Jacobi symbol, and `ε_q = 1` if `q ≡ 1 mod 4`, `ε_q = i` if `q ≡ 3 mod 4`.

**Our F̂ Theorem (verified, `FHAT_THEOREM_VERIFICATION_RESULTS.md` §8).** For every prime `p ≥ 3` and `r ≥ 1`, the full-period DFT `F̂_p^full` satisfies:
- support `{p·a : a ≡ c (mod p), a ∈ Z/p^r}`, size `p^{r-1}`,
- `|F̂_p^full(p·a)| = p^{(r+3)/2}` uniformly on support.

**Structural comparison.**

| Feature | Lemma 2.15 | Our F̂ theorem |
|---|---|---|
| Object | `S_{q,d}(ψ, k)`, length `q/d` | `G[a] = F̂(p·a)/p`, length `q/p = p^r = q/d` (matches) |
| Hypothesis | `(q,3) = 1, d² | q | d³` | `p ≥ 3, r ≥ 1` (no q‑mod-3 restriction; depth controlled by r) |
| Support condition (vanishing) | `k ≢ −a_ψ (mod d)` ⇒ vanish | `a ≢ c (mod p)` ⇒ vanish |
| Support cardinality | one residue class mod d, of size `q/d²` | `{a ≡ c mod p} in Z/p^r`, size `p^{r-1} = q/d²` (matches) |
| Magnitude | `\|S_{q,d}\| = √q` (since `\|ε_q\| = \|Jacobi\| = 1`) | `\|F̂(p·a)\| = p · √q`; equivalently `\|G[a]\| = √q · √(q/p)`... let's check: `p^{(r+3)/2} / p = p^{(r+1)/2}`, and `√(q/d²) = √(p^{r-1})` is `p^{(r-1)/2}`. So `\|G[a]\| = p^{(r+1)/2}` = `p · p^{(r-1)/2}` = `p · √(q/d²)`. Compare GY: `\|S_{q,d}\| = √q = √(p^{r+1}) = p^{(r+1)/2}`. **Magnitudes match exactly: `\|S_{q,d}\| = \|G[a]\| = p^{(r+1)/2}`.** |
| Explicit form on support | Closed form: Gauss-sum `ε_q √q × (Jacobi)(q/d²) × e_q(2a_ψ̄(k+a_ψ)²)` | Theorem 78.6 (`result_78_extended.md`): `G(a)/√q = e_q(P_a(s*(C_a)))` with saddle `s*(C_a) = (C_a − 1)/3`, modulo Hensel corrections at `r ≥ 4` |

**The structural match is exact for magnitude** and **partial for closed form**: GY have a clean closed-form Gauss-sum-style evaluation under their hypothesis `d² | q | d³`; we have an explicit saddle-point closed form (Theorem 78.6) under our setting `d = p, q = p^{r+1}` which has Hensel corrections at `r ≥ 4`. The hypothesis mismatch is the crux — see §4 of `GARCIA_YOUNG_MECHANISM.md` and adversarial check A2 of `GARCIA_YOUNG_TRANSLATION_ATTEMPT.md`.

Importantly: **GY's Lemma 2.15 gives the evaluation of a single F̂(p·a) value, NOT a bilinear sum.** The √q in Lemma 2.15 is the magnitude of one Gauss-sum/Postnikov sum, not the bilinear cancellation we need. See `GARCIA_YOUNG_MECHANISM.md` for where the √q **bilinear** secondary main term `A` arises in the proof.
