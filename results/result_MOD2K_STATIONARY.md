# RESULT — MOD2K-STATIONARY: the inverse-tree node measure has 2-adic valuation ~ Geometric(1−1/λ), λ=(3+√21)/6 (2026-07-27)

**Probe:** `probes/probe_mod2k_stationary.py`. Closes inverse-tree Phase-5 open question #2 (`inverse_tree_findings.md`):
*"closed form for the residue stationary distribution mod 2^k"* + *"explain the increasing chi² (residues do NOT
equidistribute)."* Both answered. Builds directly on `result_BRANCH_BIAS.md` (same λ).

## Operator (exact, and why)
Node `n`: children `left=2n` (always), `right=(n−1)/3` (only n≡4 mod 6). **3 is a unit mod 2^k** (`3^{-1} mod 64 = 43`),
so the `/3` right-child map is **exact** on the 2-adic part — `(n−1)/3 mod 2^k = (r−1)·43 mod 2^k`, no 3-adic lift needed.
Only the branching *test* (n≡4 mod 6) needs the mod-3 component. So the offspring operator on states
`(u = n mod 2^k, w = n mod 3^J)` is exact; its Perron eigenvector (dominant λ) is the stationary node measure. The mod-2^k
marginal is `J`-stable to 5e-16 (J=1 vs 2), λ = (3+√21)/6 = 1.263763 at every k=3..6, and the marginal **matches the empirical
`tree_d50.parquet` deep-layer distribution to |err|max = 2×10⁻⁴** at k=3..6.

## CLOSED FORM — the 2-adic valuation is Geometric
The stationary mass at 2-adic valuation `v₂(n)=v` is **exactly geometric** (machine-precision, err ≤ 5×10⁻¹⁶):
> **`P(v₂(n)=v) = (1 − 1/λ)·(1/λ)^v`**,  ratio `1/λ = (√21−3)/2 = 0.791288`,  success `1−1/λ = (5−√21)/2 = 0.208712 = P(n odd)`.
```
 v   mass(v) empirical   (1-1/λ)(1/λ)^v      err
 0     0.208712          0.208712          2.5e-16
 1     0.165151          0.165151          1.1e-16
 2     0.130682          0.130682          8.3e-17
 3     0.103407          0.103407          4.2e-17
 4     0.081825          0.081825          0.0
 5+    0.310222          0.310222(tail)    5.0e-16
```

## Mechanism (renewal) + why it never equidistributes
The valuation does a **+1/reset renewal**: the always-on doubling child `2n` pushes `v₂ → v₂+1`; the branch child `(n−1)/3`
is **odd**, resetting `v₂ → 0`. The stationary distribution of a +1/reset renewal is geometric — hence the exact
`Geometric(1−1/λ)`. This is precisely why **residues do NOT equidistribute**: mass concentrates on high-valuation (highly
even) residues as `(1/λ)^v`, never uniform. The reported "chi² increasing" is increasing **in the modulus resolution k**
(1.40, 2.22, 3.31, 4.78 at k=3,4,5,6) — the measure is non-uniform at every 2-adic scale — while chi² **at fixed k is stable
across depth** (~4.78 mod 64, depths 20→55), confirming convergence to this fixed point, not divergence.

## What remains (finer sub-structure, honest)
The valuation *marginal* is closed. The distribution **within** a valuation level is non-uniform (odd residues span
8×10⁻⁵…6.5×10⁻²): odd residues are populated only by the branch pushforward `r ↦ (r−1)·3^{-1} mod 2^k`, so the within-level
law is the self-referential pushforward of the ≡4-mod-6 sub-measure — a recursive sub-structure, not yet in closed form. The
one-line headline (valuation ~ Geometric(1−1/λ)) is exact and fully explains the equidistribution failure.

## Net
Inverse-tree stationary node measure mod 2^k: **valuation `v₂ ~ Geometric` with ratio `1/λ=(√21−3)/2`**, matching the tree to
2×10⁻⁴, mechanism = doubling(+1)/branch(reset). Closes Phase-5 OQ#2 at the marginal level and explains the chi² non-uniformity
as geometric valuation concentration. Same `λ=(3+√21)/6` as `result_BRANCH_BIAS.md` — the growth rate, the branching bias, and
now the valuation law are one algebraic object. 2-adic prefix arc; not at stake for the 7/15 / SOLSTICE / GARSIA program.
