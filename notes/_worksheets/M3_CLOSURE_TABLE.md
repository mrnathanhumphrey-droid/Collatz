# M3_CLOSURE_TABLE — parameterized Nisoli closure inequality evaluation

**Date:** 2026-05-11. Phase 3 of M3 probe. Evaluates `|K| · K^{−A} · M_3 < 1` over candidate (K, A) values, with M_3 set at the project's quoted range and |K| set by the bilinear bound disposition (PATH2 strict 2√N at r ≤ 3, HENSEL polylog-free 2√p · √N at r ≥ 4).

## 1. Reading the closure inequality

From `result_77_2_nisoli_certification.md` §3.1 + PRECISE_ASK.md §4, the Nisoli Lemma 2.9 closure requires:

> **η  =  ε_K · M_3  <  1**

with `ε_K = |K| / √q` (PRECISE_ASK §4) and `q = 3^{K}` for K-th level data (q = p^{r+1} with r = K−1, p = 3, so q = 3^K).

Substituting:

> **(|K| / √q) · M_3  <  1**
> ⟺ `|K| < √q / M_3`
> ⟺ `|K| · q^{−1/2} · M_3 < 1`

The user's task statement parameterizes this as `|K| · K^{−A} · M_3 < 1`. The mapping is:

- The literal factor is `q^{−1/2} = 3^{−K/2}` — a fixed exponential decay in K.
- Tao Prop 1.17 provides a tighter decay `n^{−A}` for the (unrelated) full-period Fourier coefficient `|μ̂_n(ξ)|`, where the "n" is the level index and A is Tao's parameter (uniform in n, ξ).
- The user's parameterization `K^{−A}` is **Tao's spectral decay**, which can be inserted **in addition to** (or replacing) the trivial `q^{−1/2}`.

We tabulate over **both** the trivial decay and Tao-with-A parameterization.

## 2. Bilinear bound |K| at level K

From `PATH2_DISPOSITION.md` + `HENSEL_DISPOSITION.md`:

- **r ≤ 3 (i.e., K ≤ 4):** `|K(r)| ≤ 2√N` strict rigorous (PATH2 family-level).
- **r ≥ 4 (i.e., K ≥ 5):** `|K(r)| ≤ 2√p · √N` polylog-free (HENSEL Hensel-digit-extraction), with `p = 3`, so `|K| ≤ 2√3 · √N ≈ 3.46 · √N`.
- The TIGHTEN_* agent (running in parallel) may upgrade r ≥ 4 to strict 2√N. We use the safer polylog-free bound here per A1.

`N = p^{r−1} = 3^{r−1} = 3^{K−2}`. So:

| K | r = K−1 | N = 3^{K−2} | √N | |K| bound (r ≤ 3: 2√N, r ≥ 4: 2√3·√N) |
|---|---------|-------------|------|---------------------------------------|
| 6 | 5 | 81 | 9 | 2√3·9 ≈ 31.18 |
| 10 | 9 | 6561 | 81 | 2√3·81 ≈ 280.6 |
| 15 | 14 | 4,782,969 | 2187 | 2√3·2187 ≈ 7575 |
| 20 | 19 | 3.49 × 10⁹ | ≈ 59049 | 2√3·59049 ≈ 204519 |
| 30 | 29 | 2.06 × 10¹⁴ | ≈ 1.43 × 10⁷ | 2√3·(1.43 × 10⁷) ≈ 4.96 × 10⁷ |

## 3. M_3 value

From R77.2 §3.3 (with caveat: against falsified spectrum):

- **R77.2 quoted range:** `M_3 ≈ 800 – 1000` (sharper estimate)
- **R77.2 crude upper:** `M_3 ≤ 11320`
- **Approach A (this probe):** `M_3 ≤ 944`, `M_3 ≥ 8`
- **Anticipated numerical (Approach C):** `M_3 ∈ [50, 200]`

We tabulate with **M_3 = 1000** as the canonical "loose-but-explicit" project value, and **M_3 = 100** as an optimistic anticipated numerical value.

## 4. Closure table — trivial decay `q^{-1/2}` only (no Tao)

Inequality: `|K| · q^{−1/2} · M_3 < 1`, i.e., `|K| · M_3 < √q`.

| K | q = 3^K | √q | |K| bound | M_3 = 1000 product `|K|·M_3` | <√q? | M_3 = 100 product | <√q? |
|---|---------|-----|----------|------------------------------|------|-------------------|------|
| 6  | 729       | 27.0     | 31.2     | 31,200       | NO (1156× too large)  | 3,120        | NO (115× too large) |
| 10 | 59,049    | 243.0    | 280.6    | 280,600      | NO (1154× too large)  | 28,060       | NO (115× too large) |
| 15 | 14.3 M    | 3,788    | 7,575    | 7.58 M       | NO (2001× too large)  | 757,500      | NO (200× too large) |
| 20 | 3.49 G    | 59,049   | 204,519  | 2.05 × 10⁸   | NO (3464× too large)  | 2.05 × 10⁷   | NO (346× too large) |
| 30 | 2.06 × 10¹⁴ | 1.43 × 10⁷ | 4.96 × 10⁷ | 4.96 × 10¹⁰ | NO (3464× too large) | 4.96 × 10⁹ | NO (346× too large) |

**Verdict (no Tao):** Closure **fails at every tabulated K**, with no improvement as K grows because `|K|` and `√q` both scale `√q · 2√3/√q · √q = O(√q)` essentially — i.e., the ratio `|K| / √q` is asymptotically `2√3 ≈ 3.46`, **constant in K**. M_3 needs to be `< 1/3.46 ≈ 0.29` for trivial-decay closure, impossible since M_3 ≥ spectral-radius lower bound 8.

**Structural fact:** With polylog-free bilinear `|K| ≤ 2√3·√N` and `√q = √3·√N`, the ratio `|K|/√q = 2`. With strict 2√N (PATH2 r ≤ 3 or TIGHTEN_* upgrade), the ratio is `2/√3 ≈ 1.155`. **Trivial decay never closes** without an additional decay factor.

## 5. Closure table — with Tao Prop 1.17 `K^{−A}` decay

Inequality: `|K| · q^{−1/2} · K^{−A} · M_3 < 1`, i.e., `K^{−A} < √q / (|K| · M_3)`.

From §4, `|K| · M_3 / √q` is `2 · M_3 / √3` for r ≥ 4 (polylog-free) and `2 · M_3` for r ≤ 3 (strict). With M_3 = 1000: required `K^{−A} < 1/(2 · 1000 / √3) ≈ 1/1155` for K ≥ 5, or `K^{−A} < 1/2000` for K ≤ 4. With M_3 = 100: required `K^{−A} < 1/115` and `1/200` respectively.

| K | required `K^{−A}` (M_3=1000) | A needed (M_3=1000) | required `K^{−A}` (M_3=100) | A needed (M_3=100) |
|---|------------------------------|---------------------|------------------------------|---------------------|
| 6 (r=5) | 1/1155 = 8.66 × 10⁻⁴ | log(1155)/log(6) ≈ **3.94** | 1/115 = 8.7 × 10⁻³ | log(115)/log(6) ≈ **2.65** |
| 10 (r=9) | 1/1155 | log(1155)/log(10) ≈ **3.06** | 1/115 | log(115)/log(10) ≈ **2.06** |
| 15 (r=14) | 1/1155 | log(1155)/log(15) ≈ **2.60** | 1/115 | log(115)/log(15) ≈ **1.75** |
| 20 (r=19) | 1/1155 | log(1155)/log(20) ≈ **2.36** | 1/115 | log(115)/log(20) ≈ **1.58** |
| 30 (r=29) | 1/1155 | log(1155)/log(30) ≈ **2.07** | 1/115 | log(115)/log(30) ≈ **1.39** |

**Verdict (with Tao K^{−A}):** Closure **can hold** if Tao Prop 1.17 delivers effective constant `C_A = 1` at any of these A values:

- **M_3 = 1000:** Need `A ≥ 4` at K = 6, dropping to `A ≥ 2.1` at K = 30.
- **M_3 = 100:** Need `A ≥ 2.7` at K = 6, dropping to `A ≥ 1.4` at K = 30.

User's task says **Tao Prop 1.17 plausibly delivers A ∈ {2..10}**. So:

- At M_3 = 100 (optimistic numerical), closure holds at A ≥ 3 for **all** K ≥ 6, well inside Tao's plausible range.
- At M_3 = 1000 (R77.2 loose estimate), closure holds at A ≥ 4 for K = 6, at A ≥ 3 for K ≥ 10; **still inside** Tao's range.

## 6. Closure table — explicit (K, A) grid

`|K| · K^{−A} · M_3 < 1` with `|K| = 2√3 · 3^{(K−2)/2}` for K ≥ 5 (and trivial `q^{−1/2}` collected into K^{−A}; we use the formulation that |K|/√q is constant and K^{−A} carries the decay):

Define `Q(K, A, M) = (|K|/√q) · M · K^{−A} = (2 if K≤4 else 2/√3 · 2 = 3.46/1.73 = 2) · M · K^{−A}`.

Wait — more carefully: from §4 the ratio `|K|/√q` is exactly `2/√3 · √3 = 2` for r ≥ 4 polylog-free, and `2/√3` for r ≤ 3 strict. The "2" comes from `2√3 · √N / (√3 · √N) = 2`. **Constant.** So Q simplifies to `2 · M_3 · K^{−A}` for K ≥ 5 (and `(2/√3) · M_3 · K^{−A} ≈ 1.155 · M_3 · K^{−A}` for K ≤ 4).

| K | A | K^{−A} | Q (M_3=1000) | <1? | Q (M_3=100) | <1? |
|---|---|--------|---------------|------|--------------|------|
| 6 | 2 | 0.0278 | 55.6 | NO | 5.56 | NO |
| 6 | 3 | 0.00463 | 9.26 | NO | 0.926 | YES |
| 6 | 5 | 1.29 × 10⁻⁴ | 0.257 | YES | 0.0257 | YES |
| 6 | 10 | 1.65 × 10⁻⁸ | 3.3 × 10⁻⁵ | YES | 3.3 × 10⁻⁶ | YES |
| 6 | 20 | 2.74 × 10⁻¹⁶ | 5.5 × 10⁻¹³ | YES | 5.5 × 10⁻¹⁴ | YES |
| 10 | 2 | 0.01 | 20 | NO | 2.0 | NO |
| 10 | 3 | 0.001 | 2.0 | NO | 0.20 | YES |
| 10 | 5 | 1 × 10⁻⁵ | 0.020 | YES | 0.002 | YES |
| 10 | 10 | 1 × 10⁻¹⁰ | 2 × 10⁻⁷ | YES | 2 × 10⁻⁸ | YES |
| 10 | 20 | 1 × 10⁻²⁰ | 2 × 10⁻¹⁷ | YES | 2 × 10⁻¹⁸ | YES |
| 15 | 2 | 0.00444 | 8.89 | NO | 0.889 | YES |
| 15 | 3 | 2.96 × 10⁻⁴ | 0.593 | YES | 0.0593 | YES |
| 15 | 5 | 1.32 × 10⁻⁶ | 0.00263 | YES | 2.63 × 10⁻⁴ | YES |
| 15 | 10 | 1.73 × 10⁻¹² | 3.5 × 10⁻⁹ | YES | 3.5 × 10⁻¹⁰ | YES |
| 15 | 20 | 3.0 × 10⁻²⁴ | 6.0 × 10⁻²¹ | YES | 6.0 × 10⁻²² | YES |
| 20 | 2 | 0.0025 | 5.0 | NO | 0.50 | YES |
| 20 | 3 | 1.25 × 10⁻⁴ | 0.25 | YES | 0.025 | YES |
| 20 | 5 | 3.125 × 10⁻⁷ | 6.25 × 10⁻⁴ | YES | 6.25 × 10⁻⁵ | YES |
| 20 | 10 | 9.77 × 10⁻¹⁴ | 1.95 × 10⁻¹⁰ | YES | 1.95 × 10⁻¹¹ | YES |
| 30 | 2 | 1.11 × 10⁻³ | 2.22 | NO | 0.222 | YES |
| 30 | 3 | 3.7 × 10⁻⁵ | 0.074 | YES | 0.0074 | YES |
| 30 | 5 | 4.1 × 10⁻⁸ | 8.2 × 10⁻⁵ | YES | 8.2 × 10⁻⁶ | YES |
| 30 | 10 | 1.69 × 10⁻¹⁵ | 3.4 × 10⁻¹² | YES | 3.4 × 10⁻¹³ | YES |

## 7. Reading the table

**With M_3 = 1000 (R77.2 loose estimate):**
- A = 2 closes **nowhere** in the table.
- A = 3 closes at **K ≥ 15**, NOT at K = 6, 10.
- A = 5 closes everywhere from K = 6 up.
- A = 10 closes everywhere with massive margin.

**With M_3 = 100 (anticipated numerical):**
- A = 2 closes at **K ≥ 20**, NOT at K = 6, 10, 15.
- A = 3 closes everywhere K ≥ 6.
- A = 5+ closes everywhere with margin.

**Boundary case:** At Tao's likely-lowest A = 2 and our middle-of-the-road M_3 = 100, closure requires K = 20 — i.e., level-19 data with q = 3^20 ≈ 3.5 × 10⁹ modulus. The bilinear sum at K=20 has |K| ≤ 2√3 · √(3^18) ≈ 2.05 × 10⁵; not computationally bad but well above the current empirical range (R79b goes to K = 17 = r = 16, with |K_max| ≈ 11022).

## 8. Caveats (Phase 4 A1, A2, A3, A4)

(A1) **Bilinear bound at r ≥ 4 honest scope:** uses 2√p·√N polylog-free per HENSEL_DISPOSITION.md. If TIGHTEN_* succeeds upgrading to strict 2√N, the constant `|K|/√q` improves from `2` to `2/√3 ≈ 1.155`, shaving ~0.24 from required `log(K^A)` — small effect on closure requirements.

(A2) **K^{−A} placeholder honesty:** Tao Prop 1.17's effective `C_A` is **INFEASIBLE this session** (R77.2 §3.4; Route 1 result). Parameterizing over A assumes `C_A = 1`; **if `C_A` turns out large, the entire table shifts**. R77.2 §3.4: "C_A could be e^{200}" via Banks-Shparlinski-style buried constants. The table is therefore **conditional** on `C_A` being O(1).

(A3) **R77.3 spectrum falsification:** M_3 ≈ 1000 or 100 is computed against the falsified T_3. The "correct" T (which would govern actual ε_n) is **uncharacterized** (R77.4 erratum §parked). Without a corrected M_3 on a corrected T, the entire closure table is **structurally hollow** — it parameterizes a quantity whose operator-theoretic anchor doesn't exist.

(A4) **Operator-norm vs spectral-radius:** T_3 is non-normal (companion form). Approach A's M_3 ≤ 944 uses the κ(V) condition number; this is the operator norm upper bound. Approach C's anticipated 50–200 is the likely actual operator norm. Pre-R77.3, M_3 = 100 (vs 1000) shifts required A by ~1 across the table. Post-R77.3, the operator-norm precision is moot.

## 9. What closure would actually require, end-to-end

Putting Phase 1 + Phase 2 + Phase 3 together, **Nisoli closure for c = 7/45** requires:

1. **Effective Tao C_A:** This session marks Route 1 INFEASIBLE for this. Independent project required to extract effective constants from Tao §7.2–7.3.
2. **Characterized operator T with spectrum captured by a finite T_K:** R77.3 falsifies the 3-mode candidate; R77.4 erratum rules out the natural K_k; no characterized alternative exists in the project today.
3. **M_T < ~100 on a contour isolating 1/2 (or the rate eigenvalue, whatever it is):** depends on (2). If T exists and is "nicely conditioned" near 1/2, M can be O(100). If T turns out to have many eigenvalues clustering near 1/2 (R77.4 erratum's "no single shape dominates"), M can grow with truncation.
4. **K and A in the (table-feasible) range:** if items 1–3 deliver, table §6 shows closure at A ≥ 3 (M=100) or A ≥ 5 (M=1000) for K ≥ 6.

Items 1, 2, 3 are all currently open. The closure table is well-defined as a **conditional/parameterized object**, not a rigorous closure path.

## 10. Files

- `M3_DEFINITION.md` — Phase 1 definition.
- `M3_APPROACH_A.md` — Phase 2A spectral lower/upper bounds.
- `M3_APPROACH_B.md` — Phase 2B perturbation series diverges.
- `M3_APPROACH_C.md` — Phase 2C numerical specification + anticipated 50–200.
- `M3_CLOSURE_TABLE.md` (this file) — Phase 3 parameterized table.
- `M3_DISPOSITION.md` — top-level disposition.
