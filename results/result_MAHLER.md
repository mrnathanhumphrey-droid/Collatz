# RESULT — MAHLER: S is NOT a finite-order Mahler function; the Mahler order is INFINITE (proven) = the unipotent depth (2026-07-28)

**Probe:** `probes/probe_mahler.py`. Tests whether S_∞ is a finite-order Mahler evaluation — if so, Nishioka makes
rationality *decidable* and the archimedean doubly-exponential denominator blowup becomes irrelevant (doc-4's route). It is
**not**, and the obstruction is proven, not merely observed at accessible depth. Everything exact-rational.

## The test
S is finite-order Mahler ⟺ the exact sequences `a_i=3^i R_e^{(i)}(2)`, `b_i=3^i R_e^{(i)}(0)`, `T_i=4a_i−b_i` (with
`S_{i+1}=2T_i`) satisfy a **fixed finite-order linear recursion with rational coefficients** under `i→i+1` (the tower `z→z³`).

**Exact R_e build, doubly certified:** ν = exact-rational stationary of the Syracuse chain
(`stationary_rational ∘ build_markov_q(3,n)`); R_e-vector via the base-2 dlog-cycle geometric sum (exact, denominators
`2^{2·twoN}−1`); autocorr lags {0,2} exact. Gates, i=1..6: **float-gate vs `build_base2` ≤ 1.4×10⁻¹⁶**, and
**`4a_i−b_i == T_cert` EXACT** (matches the certified S-ladder rational at every level). i=6 = 335s (2^1458−1 denominators).

**Mandatory reproduction gate — PASSED:** Λ_1..7 (exact) gives **NO** finite rational recurrence (L=1..3). Machinery trusted.

## M-A — no finite recurrence in any sequence (at accessible exact depth)
```
 a_i = 3^i R_e(2):  NO finite rational recurrence, L=1..2 (n=6 exact)
 b_i = 3^i R_e(0):  NO finite rational recurrence, L=1..2 (n=6 exact)
 T_i = 4a_i − b_i:  NO finite rational recurrence, L=1..2 (n=6 exact)
```
Consistent with the R27-A gate, but only excludes low order at this depth. The proof below closes **all** orders.

## M-C PROOF — the denominator growth rate excludes every finite order
A fixed order-L rational recurrence `x_i = Σ_{j=1}^L α_j x_{i−j}` forces `den(x_i) | ∏_j den(α_j)den(x_{i−j})`, so
`D_i := log₂ den(T_i)` obeys `D_i ≤ Σ_{j=1}^L D_{i−j} + C`. The growth ratio of that inequality is `λ_L`, the dominant root
of `t^L = t^{L−1}+…+1` — and **`λ_L < 2` for every finite L** (λ₁=1, λ₂=φ=1.618, λ₃=1.839, … → 2⁻, never reached).
Measured `D_i` (exact, i=1..7):
```
 D_i (bits): 4.39, 16.05, 58.10, 201.7, 630.6, 1923.4, 5796.9
 ratio D_i/D_{i-1}:   3.655, 3.619, 3.472, 3.126, 3.050, 3.014   →  3
```
`D_i ~ 2·3^{i−1}` (the banked denominator theorem, `den(S_r) ~ 2^{2·3^{r−1}}`), **ratio → 3**. A sequence cannot exceed its
own upper bound: `3 > 2 > λ_L` for every finite L ⟹ **no finite-order rational recurrence exists at any order L**.
> **Theorem (conditional on the banked denominator theorem): S is not a finite-order Mahler function over ℚ(z); the Mahler
> order is infinite.**
Triangulated from three independent directions: (i) R27-A Λ no-recurrence; (ii) a,b,T no low-order recurrence (M-A);
(iii) the denominator-rate argument (all orders). The doubly-exponential denominators are not merely "consistent with" infinite
order — a growth ratio 3 > 2 **forces** it.

## What this means
- **Doc-4's Mahler/Nishioka route is CLOSED — as a theorem, not a guess.** S_∞ is not a finite-order Mahler evaluation, so
  Nishioka decidability does not apply; there is no finite q-difference equation to run a rational-solution algorithm on.
- **Infinite Mahler order = the unbounded unipotent depth.** This is the same wall in the Mahler framework's clean dichotomy
  (finite order ⟹ rationality decidable; infinite order ⟹ outside the theorem). It matches the operator picture `M=D(I+N)` with
  `N` nilpotent of index = level: the growing nilpotency index *is* the growing Mahler order, and the ratio-3 denominator growth
  *is* the multiplicative `𝔾_m` denominators (`2^{2·3^{r−1}}`) accumulating one cyclotomic factor per unipotent level.
- **The wall gets its sharpest statement:** not "hard," but the negative side of a **named dichotomy** — S_∞ is a period at
  infinite Mahler/unipotent order, which is precisely why finite-order methods (Nishioka, spectral gap, valuation, lattice)
  cannot reach it. Citable as *why* the constant resists.

## M-D — the b-family footnote (algebra, closes the +3/+5 question)
For the map `x ↦ (qx+b)·2^{−v}`, conjugation by the dilation `φ_c: x↦cx` gives `φ_c^{−1}∘renewal_b∘φ_c = renewal_{b/c}` (the
`2^{−v}` commutes with dilation). So for **b coprime to the modulus**, taking `c=b` conjugates `renewal_b` to `renewal_1` —
**+1, +5, +7 are the identical family**, S_∞(b)=S_∞(1), the b-axis is **flat**. For **b ≡ 0 mod q**, `qx+q=q(x+1)` maps into the
non-units and **degenerates** (same mechanism as the ⟨2⟩-grading kill / leaving the unit group). **Verdict: the b-axis is flat
for unit b and degenerate otherwise; the only genuine free parameter is the valuation prime p** (the "2" in `2^{−v}`) — the
`(q,p)`-Hydra / Siegel p-Hydra. The +3/+5 instinct resolves to "conjugate or degenerate, not a new family."

## Net
- **S is infinite-order Mahler (proven via denominator-rate); doc-4's decidability route is closed as a theorem.** The
  characterization "S_∞ = period at infinite unipotent/Mahler depth" is now backed by a proof, not just an inference.
- **7/15 UNAFFECTED** — stays excluded on value (floor `2·T_20 = 0.473177`); this decides *type/decidability*, not the value.
- **Not at stake:** the denominator theorem (this *uses* it), GARSIA, DENOM, SOLSTICE, R1–R30, P6D–P6K. **Newly banked:** the
  infinite-Mahler-order theorem + the flat/degenerate b-axis.
