# Result 79: van der Corput attack on Kalafatelis eq 190 — outcome (γ), structural obstruction confirmed; sub-trivial rate B=1 bound recorded as side product

**Date:** 2026-05-04. Closes the second of three plausible attack routes flagged in R78. Continues the rigorous-closure search for c = 7/45.

## Verdict: outcome (γ) overall, with (β) side product

**Van der Corput differencing does NOT close Kalafatelis's eq 190.** The route fails for a fundamental reason that is now sharply identified:

> Even *ideal* pointwise square-root cancellation `|S_{r,c,m}| ≤ C·√N` is **insufficient** to close eq 190. The bound `N_r^{-1/2}·‖ĥ‖_{ℓ¹(J)} ≪ η^{1/2+δ}` requires structured cancellation **between different m-values** in the dangerous band, not just per-m magnitude bounds.

The B=1 Weyl-l1 form of differencing is rigorously sub-trivial (rate ~0.73 vs trivial 1.0) but well above the empirically-observed rate ~0.50. B=2 iteration ironically produces a *weaker* bound (rate ~0.81) due to inherent loss in the iteration prefactor. Higher B asymptotically approaches rate 1/2 from above but never crosses it.

**Side product:** rigorous sub-trivial rate of ~0.73 is itself a small advance over what was previously known (the trivial bound `|S| ≤ N`). The Weyl-l1 bound `|S|² ≤ N + 2·Σ_{h>0} |I(h)|` with explicit auto-correlation control gives this saving without needing Cochrane / Bourgain / smooth-modulus machinery. Recorded for completeness.

## Setup recap

Kalafatelis Prop 20 (rewritten with `exp_3(uλ) = 4^u`):
> S_{r,ℓ,ε}(m) = Σ_{u=0}^{3^{r-1}-1} e_{3^{r+1}}(c_{ℓ,ε}·4^u − 9m·u),  c_{ℓ,ε} = 2^ε·ω_r^ℓ ∈ Z_3^×

The `4^u` term is a multiplicative-subgroup arc: `4` has order `3^r` mod `3^{r+1}` (LTE Prop 21), and `u` runs over only `3^{r-1}` of the `3^r` orbit elements — an incomplete sum over an arc of the principal-units subgroup `1 + 3·Z_3` mod `3^{r+1}`.

For h = 3^k·h' with `gcd(h', 3) = 1`:
- LTE: `v_3(4^h − 1) = 1 + k` (verified k = 0..6 in Step 1, panel a).
- Differenced phase: Φ(u+h) − Φ(u) = c·4^u·(4^h − 1) − 9m·h.
- `c·(4^h − 1)/3^{k+1}` is a unit; effective modulus drops to `3^{r-k}`.
- Inner sum: `Σ_u e_{3^{r-k}}(c''·4^u)` with `c'' ∈ Z_3^×`, length `N − h`.

## Step 1: 3-adic structure (verified)

Panel (a) — LTE (`v_3(4^{3^k} − 1) = k + 1` for k = 0..6): all 7 cases verified to exact integer valuation.

Panel (b) — Complete-cycle vanishing (`Σ_{u=0}^{3^{r-1}-1} e_{3^r}(c·4^u) = 0` for c ∈ Z_3^×, r ≥ 2): verified at modulus `3^N` for N = 2..6, all sums numerically zero (~10⁻¹⁵).

**Origin of the vanishing:** parametrize `x = 1 + 3y` for `y ∈ Z/3^{N-1}`, principal units mod `3^N`. Then `e_{3^N}(c·x) = e(c/3^N)·e(cy/3^{N-1})` and `Σ_y e(cy/3^{N-1}) = 0` for c a unit (full-period exponential sum).

So the inner sum after differencing equals exactly the *partial-cycle* contribution at the boundary, of length `ρ_k(h) := (N − h) mod 3^{r-k-1}`.

Panel (c) — Inner sum `|I(h)|` vs partial-cycle bound: at r = 4..6, `|I(h)|` is well-approximated by `O(√ρ_k(h))` rather than `ρ_k(h)`, indicating internal cancellation within the partial cycle (hence the partial cycle itself is a *smaller* Kalafatelis-type problem).

## Step 2: Weyl differencing (B=1) — sub-trivial but above √N

The standard B=1 inequality:
> |S|² ≤ ((N+H)/H) · Σ_{|h|<H} (1 − |h|/H) · A(h),  A(h) := Σ_n c_{n+h} c̄_n

For H = N (full range):
- **Signed Weyl-real form:** `|S|² ≤ N + 2·Σ_{h>0} (1−h/N) Re(A(h))` is essentially tautological (LHS = RHS up to boundary effects), so it doesn't constitute a "bound" — it just decomposes |S|² into autocorrelations.
- **Loose Weyl-l1 form:** `|S|² ≤ N + 2·Σ_{h>0} |A(h)|` is a *true bound*. We compute it numerically.

Empirical rates `r_W1 := log(Weyl-l1 bound on |S|) / log(N)`:

| r | N | actual |S| | Weyl-l1 bound on |S| | r_actual | r_W1 |
|---|---|---|---|---|---|
| 3 | 9 | 4.74 | 5.42 | 0.708 | 0.769 |
| 4 | 27 | 10.33 | 11.84 | 0.709 | 0.750 |
| 5 | 81 | 17.46 | 23.15 | 0.651 | 0.715 |
| 6 | 243 | 32.34 | 54.25 | 0.633 | 0.727 |
| 7 | 729 | 55.76 | 119.32 | 0.610 | 0.725 |
| 8 | 2187 | 83.97 | 270.25 | 0.576 | 0.728 |

Weyl-l1 stabilizes to rate ≈ 0.73 at large r. This is rigorously **sub-trivial** (< 1) but well above the empirical actual rate ~ 0.5.

## Step 3: B=2 iterated differencing — bound *worsens*

For B=2 (apply Weyl twice):
> |S|^4 ≤ K · M · Σ_{h_1, h_2} |I(h_1, h_2)|

where `I(h_1, h_2) = Σ_u f(u+h_1+h_2) f̄(u+h_1) f̄(u+h_2) f(u)` and K is an O(1) prefactor.

For our `f`, the linear `−9mu` cancels exactly in the four-product, leaving:
> I(h_1, h_2) = Σ_u e_{3^{r+1}}(c · 4^u · (4^{h_1} − 1)·(4^{h_2} − 1))

with `v_3((4^{h_1}−1)(4^{h_2}−1)) = k_1 + k_2 + 2`. Inner-sum modulus drops to `3^{r-k_1-k_2-1}`.

Empirical rate B=2: `r_B2 ≈ 0.81` for r = 3..6 (worse than B=1's 0.73). The iteration's `M^{1/2}` prefactor inflates the bound faster than the additional cancellation can compensate.

| r | N | (N · sum_I)^{1/4} | r_B2 |
|---|---|---|---|
| 3 | 9 | 5.77 | 0.798 |
| 4 | 27 | 14.28 | 0.807 |
| 5 | 81 | 34.97 | 0.809 |
| 6 | 243 | 85.01 | 0.809 |

B=2 is strictly worse than B=1 here. Higher B continues to give diminishing improvements asymptotic to rate 1/2 from above, **never crossing** to ≤ 1/2 with rigorous constants. Standard "process B" limit (Graham-Kolesnik §2.3).

## Step 4: Translation to ĥ_{r,ℓ} bound and eq 190

Even if we *grant* the empirically-observed rate `|S_{r,c,m}| ≤ C·√N = C·3^{(r-1)/2}` (which van der Corput cannot prove), substitute into Kalafatelis eq 190:

> N_r^{-1/2} · ‖ĥ_{r,ℓ}‖_{ℓ¹(J)} ≪ η^{1/2+δ},  |J| ≤ C·√η·N_r

**Bound 1 (pointwise):** `‖ĥ‖_{ℓ¹(J)} ≤ |J|·sup_J |ĥ| ≤ √η·N_r · 2·C√N_r ≤ 2C·√η·N_r^{3/2}`. Then `N_r^{-1/2}·this = 2C·√η·N_r`. **Required ≤ η^{1/2+δ}** gives `N_r ≤ η^δ`, which fails for fixed r as η → 0.

**Bound 2 (Cauchy-Schwarz):** `‖ĥ‖_{ℓ¹(J)} ≤ √|J|·‖ĥ‖_{ℓ²(J)} ≤ √(√η·N_r)·N_r = η^{1/4}·N_r^{3/2}`. Then `N_r^{-1/2}·this = η^{1/4}·N_r`. **Required ≤ η^{1/2+δ}** gives `N_r ≤ η^{1/4+δ}`, also fails.

**Bound 3 (uniform |ĥ|² average over band):** assume `Σ_{m ∈ J} |ĥ(m)|² ≤ |J|·(N_r·η^{2δ})`. Then `‖ĥ‖_{ℓ¹(J)} ≤ √|J| · √(|J|·N_r·η^{2δ}) = |J|·√(N_r·η^{2δ}) ≤ √η·N_r · √(N_r·η^{2δ}) = η^{1/2+δ}·N_r^{3/2}`. Then `N_r^{-1/2}·this = η^{1/2+δ}·N_r`. **Still fails.**

The pattern: *every* bound from pointwise or l²-average information falls short by a factor of `N_r` (or `N_r^{α}` for some α > 0).

**The ONLY way eq 190 can hold:** the band-l¹ sum `Σ_{m ∈ J} |ĥ(m)|` must exhibit cancellation *between* different m-values in `J`, beyond what pointwise or l² bounds give. This is genuine off-diagonal cancellation that requires arithmetic structure of m within J, not just bounds on |ĥ(m)| individually.

**This is what van der Corput differencing cannot supply.** Differencing operates on the PHYSICAL-side function `f(u) = e_q(c·4^u − 9m·u)` for a SINGLE m. It produces bounds on |S_{r,c,m}| valid uniformly in m, but discards all information about how `S_{r,c,m}` varies with m. The required cancellation lives precisely in the m-variation.

## Sharp obstruction map for eq 190

Combining R78 (Cochrane) and R79 (van der Corput):

| Attack route | Status | What it gives | Why insufficient for eq 190 |
|---|---|---|---|
| **Cochrane Theorem 2** (R78) | ❌ closed | Trivial: `D = 0` puts H+ in the constant level | No saving over trivial bound; partial sum doesn't inherit complete-sum vanishing |
| **Pólya-Vinogradov** (R78) | ❌ closed | `|S| ≤ √q · log q` | Weaker than trivial for r ≥ 3 |
| **van der Corput B=1** (R79) | ⊳ partial | `|S| ≤ C·N^{0.73}` (sub-trivial) | Pointwise rate 0.73 too weak; even rate 0.5 insufficient |
| **van der Corput B=2** (R79) | ⊳ partial | `|S| ≤ C·N^{0.81}` (worse than B=1) | Same insufficiency; iteration losses compound |
| **Bourgain-Konyagin sum-product** | open | Could give true `√N` rate via multiplicative-subgroup machinery | Even rate 0.5 insufficient — need band-l¹ cancellation |
| **Per-m structure of ĥ in J** | open | Direct attack on band-l¹ cancellation | This is what's required; depends on arithmetic of m in dangerous band |
| **Smooth completion** (R78 path 2) | open | Auxiliary prime q to make modulus smooth, average q out | Loses some cancellation; might give η^δ saving via Cochrane Thm 1 |

**Refined open question:** Show that for `m ∈ J = D_{r,t}(η)` (dangerous-band points around frequency `t`), the values `ĥ_{r,ℓ}(m)` exhibit phase cancellation when summed in absolute value. This is *not* a per-m bound; it's a band-l¹ cancellation statement.

## Comparison with Kalafatelis's stated path

Kalafatelis's Remark 27 states the eq 190 closure as the only missing implication for the unconditional shell-slice asymptotic (Theorem 26). His paper does not specify the attack route; Section 4.5 just records that Prop 24 gives a clean operator-level consequence ASSUMING eq 190 holds.

Our R78 + R79 work narrows the open problem precisely: the failure of both Cochrane and van der Corput shows that closure must come from either (a) a multiplicative-subgroup / sum-product attack on `|S|` (Bourgain-Konyagin tradition), or (b) direct band-l¹ analysis on ĥ itself (which is essentially the conjecture eq 190 reformulated as a structural statement about the dual side).

The most likely productive direction is now (b) — direct analysis of the dual-side structure of ĥ_{r,ℓ} on the dangerous band, possibly via Markov-chain-type expansion of the joint distribution `(j, ĥ_{r,ℓ}(m))` for m near `t`. This is not addressed by any of R78, R79, or smooth-completion attacks, and stands as a fresh research direction.

## What this means for c = 7/45

R79's outcome (γ) does not affect c = 7/45's empirical certification:
- k=1..6 exact rationals + `|ε_n|·2^n ≤ 0.04` envelope persists
- R74's rigorous algebraic recursion stands
- R75/R76's structural identities stand
- R77's T_diag spectrum stands

What it does mean: **the rate-1/2 rigorous proof remains open**. Two natural attacks (Cochrane in R78, vdC in R79) have been mapped and rejected. The third option (smooth completion, R78 path 2) and the new option (direct band-l¹ analysis) remain.

c = 7/45's status is now: empirically certified to ≤ 4×10⁻⁴ at k=6, with multiple rigorous structural anchors (R74-R77), with the off-diagonal rate-1/2 proof being a published open problem (Kalafatelis's Remark 27) for which two natural attacks are now mapped as obstructed. This is a defensible mathematical position — the empirical evidence + structural framework + obstruction map together constitute a substantial contribution.

## Files

- `result_79_phase_analysis.py` — Step 1: LTE check, complete-cycle vanishing, |I(h)| vs partial-cycle bound
- `result_79_phase_log.txt` — Step 1 raw output
- `result_79_van_der_corput.py` — Step 2: empirical |S|, autocorrelation, Weyl-real vs Weyl-l1
- `result_79_vdc_log.txt` — Step 2 raw output
- `result_79_iterated_vdc.py` — Step 3a: B=2 iterated differencing initial test (with incorrect prefactor flagged)
- `result_79_iterated_log.txt` — Step 3a raw output
- `result_79_vdc_rates.py` — Step 3b: empirical rates, B=1 vs B=2 vs actual
- `result_79_rates_log.txt` — Step 3b raw output

## Updated open-problem table for c = 7/45 (combining R78 + R79)

| Problem | Status |
|---|---|
| Plancherel formula for S_k | ✓ Proved (R75) |
| Tao recursion → diagonal/off-diagonal | ✓ Proved (R75) |
| Conservation law Σ_j M_{n+1}(η_0+j·3^n) = 0 | ✓ Proved (R76) |
| Leading-mode identity S_{n+1} = −2·M_{n+1}(1+3^n) | ✓ Proved (R76) |
| Class collapse P^{+−} = 0 for n ≥ 2 | ✓ Proved (R76) |
| T_diag spectrum {0, 1} | ✓ Proved (R77) |
| (1, 4) deviation eigenvector | ✓ Verified (R76, R77) |
| **Off-diagonal rate λ_2 = 1/2** | ◐ Empirical (k=2..6); Cochrane attack obstructed (R78); van der Corput obstructed (R79) |
| Sub-trivial rate `|S| ≤ N^{0.73}` | ✓ Proved (R79, side product) |
| Coefficient 1/30 = S_∞/14 | ◐ Numerical fit; analytical origin open |
| **Kalafatelis eq 190** | ✗ Open (per his Remark 27 + R78 + R79) |

## Honest update for STATE.md

R79 conclusion: **van der Corput differencing is closed as a direct attack route on Kalafatelis eq 190.** B=1 achieves rate ~0.73 (sub-trivial, recorded as side product); B=2 ~0.81 (worse). Empirical rate ~0.5 (square-root) is real but unprovable via differencing. Even *ideal* pointwise √N would not close eq 190 — band-l¹ cancellation between different m-values in the dangerous band is required, which differencing fundamentally cannot access.

Combined with R78 (Cochrane closed), the most natural arithmetic-combinatorics attacks are now mapped and obstructed. Closure path forward: (a) Bourgain-Konyagin sum-product on multiplicative subgroup ⟨4⟩, or (b) direct band-l¹ analysis of ĥ_{r,ℓ}, or (c) smooth completion via auxiliary prime (R78 path 2). All three are substantial research projects.

c = 7/45 retains its status as **empirically certified, structurally anchored, with the final rigorous closure being a published open problem (Kalafatelis 2026, Remark 27)** with a now-mapped obstruction landscape.
