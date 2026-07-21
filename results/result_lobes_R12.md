# Probe R12 — support law / closed loop / lobes — **A/B gates PASS, F certified, C/D/E delivered**

**Date:** 2026-07-21  Exact gates + labeled measurement. Probe `probes/probe_lobes_R12.py` (reuses R7 + R10).
Follows walk-back #31 (R11-B: U is a *sparse* unitary, not dense-flat). This session gates Wilson's corrected
kernel geometry and closes the reformulation arc.

## R12-A — SUPPORT LAW (forced): **GATE PASS**
U's primitive block has zero pattern **exactly {k ≢ ξ mod 3}**, r = 2…5: off-class entries machine-zero (≤1e-14),
in-class entries |U| = 3^{−(r−1)/2} exact, conjugate-block U(−k,−ξ) = conj U(k,ξ) exact.

| r | off-class max\|U\| | in-class \|U\|−3^{−(r−1)/2} | conj-block dev |
|---|---|---|---|
| 2 | 6.6e−16 | 1.1e−16 | 6.2e−16 |
| 3 | 1.6e−15 | 1.5e−15 | 1.9e−15 |
| 4 | 4.6e−15 | 4.0e−15 | 6.5e−15 |
| 5 | 9.9e−15 | 1.1e−14 | 1.4e−14 |

Confirms the fiber derivation (z=z₀+3^{r−1}s, β(z₀+3^{r−1}s)=β(z₀)+3^{r−1}s·β′(z₀), β′≡1 mod 3, remainder v₃≥r ⟹
inner s-sum = 3·δ[k≡ξ mod 3]). **U is (mod-3 class-preserving) ⊕ of two conjugate dense-flat blocks** — R11-B's
half-support is the *unique* structure forced by flatness + the unitarity budget, not an anomaly. Walk-back #32 not
incurred; the corrected geometry is exact. The invariant U preserves is the frequency's ± class (R64.B v-parity,
the (1,4)/(P₊,P₋) split): the bridge relocates onto the corpus's class-resolved machinery.

## R12-B — CLOSED-LOOP WELD (forced): **GATE PASS**
**Σ_{k prim mod 3^r} |θ̂(k)|² e(km/3^r) = C_{r+1}(m)/3** exact, r = 2,3,4, for m = 1, 2, 3 and DC (m=3^r) — every
tested moment byte-equal (sign convention e(+km), locked). **Every angular moment of the layer profile is a banked
C-table entry.** The four representations — additive γ, strata C̄, character layers Λ, chirp kernel — have closed on
each other with no residue: profile ⟺ C-tables ⟺ engine ⟺ ε-increments. There is no fifth coordinate; the object is
irreducibly **one scalar sequence plus one angular profile whose moments are that sequence's own raw material.**
Reformulation ends here.

## R12-F — CONVENTION CERTIFICATION (required before r≥7 use): **CERTIFIED**
The historical ε-table (`result_77_7_eps_exact_through_k8_v2`) is **byte-identical** to the exact d_k = S_k − 7/15:

| k | d_k = S_k−7/15 | historical ε_k | match |
|---|---|---|---|
| 1 | 1/5 | 1/5 | ✅ |
| 2 | 1/105 | 1/105 | ✅ |
| 3 | −5191/1019445 | −5191/1019445 | ✅ |
| 4 | −11346676448406637/4627031617157687115 | = | ✅ |
| 5 | −(60-digit)/(62-digit) | = | ✅ |
| 6 | −(189-digit)/(191-digit) | = | ✅ |

**ε_k = d_k = S_k − 7/15, no scale/sign conversion — the two campaigns' sequences are literally the same rational
table.** The historical table is **exact through k=8** (not merely float), so it hands exact ε₇, ε₈. **Cross-check:**
S₇ from the ε-table = S₇ from R11's independent renewal μ₇ build (True) — two independent routes agree. The wall is
now at k≥9 (float-only in the old ε-CSVs); the merger is certified exact for all r ≤ 7.

## R12-E — LEDGER EXTENSION (exact, via certified ε): **Λ₁…Λ₇ exact**
| r | Λ_r | float | sign | source |
|---|---|---|---|---|
| 1 | −2/21 | −9.52e−02 | − | char/R10 |
| 2 | −1490/203889 | −7.31e−03 | − | char/R10 |
| 3 | +2849957897648150/… | +1.32e−03 | + | char/R10 |
| 4 | +… | +6.50e−04 | + | char/R10 |
| 5 | +… | +3.27e−04 | + | char/R10 |
| **6** | −(…)/… | **−3.39e−04** | **−** | ε-table |
| **7** | +(…)/… | **+2.15e−04** | **+** | ε-table |

**Λ signs −,−,+,+,+,−,+** — the period-9-in-k oscillation of the ε-sequence, now read directly off the exact layer
ledger through r=7. (The μ₈ renewal build — supp 4374, autocorr ~19M Fractions — is the exact wall for the
character-side route; the certified ε-table bypasses it for r=6,7.)

## R12-C — LOBE LEDGER (measurement, NO fit): **a new stable observable**
Λ_r = Σ_{k prim}|θ̂(k)|² Re w(k/3^r), Re w(x) = (4cos2πx−1)/(17−8cos2πx); split by sign into L_r (Re w>0) and M_r
(−Σ over Re w<0):

| r | L_r (Re w>0) | M_r (Re w<0) | L_r−M_r | Λ_r (exact) | L_r+M_r |
|---|---|---|---|---|---|
| 2 | 1.7437e−2 | 2.4745e−2 | −7.31e−3 | −7.31e−3 | 4.218e−2 |
| 3 | 3.6826e−2 | 3.5506e−2 | +1.32e−3 | +1.32e−3 | 7.233e−2 |
| 4 | 3.7573e−2 | 3.6923e−2 | +6.50e−4 | +6.50e−4 | 7.450e−2 |
| 5 | 3.6499e−2 | 3.6172e−2 | +3.27e−4 | +3.27e−4 | 7.267e−2 |
| 6 | 3.5606e−2 | 3.5945e−2 | −3.39e−4 | −3.39e−4 | 7.155e−2 |

**Both lobe masses stabilize to a common value ≈ 0.0358** (L_r+M_r → ~0.072 for r≥3), while Λ_r = L_r−M_r is the
tiny oscillating residue between them — exactly Wilson's prediction (both lobes → common L_∞, Λ_r the transient
approach). **The lobe-dominance flips track the Λ sign** (r=2 M>L → Λ<0; r=3,4,5 L>M → Λ>0; r=6 M>L → Λ<0): the
period-9 phase alignment lives in *which lobe leads*, on top of a stable common mass. L_r, M_r are a new observable
the corpus never tracked — strong evidence a limiting angular mass ψ exists with ⟨ψ, Re w⟩ = 0.

## R12-D — CLASS-RESOLVED PROFILE (forced check + measurement): mirror EXACT, classes balanced
The conjugation mirror |θ̂(k)|² = |θ̂(−k)|² holds **exact** (≤1e-14) at every r, and the two mod-3 classes (k≡1,
k≡2) contribute **equally**: class-lobe sums c₁ = c₂ exactly, each half of Λ_r (r=2: −3.654e−3 each; r=3: +6.599e−4
each; …). The two mod-3 frequency classes are the balanced (P₊,P₋) ± pair — R64.B's v-parity classes as literal
chirp geometry, connecting the layer ledger to the corpus's (1,4) class machinery.

## Status
**R12: A/B gates PASS** (support law {k≡ξ mod 3} forced and exact; the reformulation loop closed — angular moments
ARE C-table entries, no fifth coordinate), **F certified** (ε_k = d_k byte-equal, exact through k=8, S₇ cross-checked
against the renewal build — the merger is exact for r≤7), **C/D/E delivered** (exact ledger to Λ₇, signs −,−,+,+,+,−,+;
lobe masses L_r,M_r stabilizing to a common ≈0.0358 = a new ψ-existence observable; class mirror exact, ± pair
balanced). The object is now irreducibly one scalar sequence + one angular profile obeying two derived constraints
(⟨ψ,Re w⟩=0 orthogonality; the half-step class-flip walk's mixing against Re w). **Still owed (pen):** Σ_{r≥1}Λ_r =
−1/10 / the decorrelation rate ρ, now attackable as the lobe-difference limit L_∞−M_∞ = 0 with the −1/10 the
integrated transient, or via the Tao-side quadratic form ⟨μ̂, K_rμ̂⟩. No fitting; exact gates, labeled numeric
lobes; walk-back #31 corrected geometry verified, #32 not incurred.
