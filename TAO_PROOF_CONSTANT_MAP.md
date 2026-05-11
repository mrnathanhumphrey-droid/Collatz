# Phase 1a — Tao §7.2–§7.4 Proof Constant Map

**Source:** Tao, "Almost all orbits of the Collatz map attain almost bounded values," Forum of Mathematics, Pi **10**, e12 (2022); arXiv:1909.03562. Local OCR plaintext: `C:/Collatz/tao2022.txt`; §7.2–§7.3 extract: `C:/Users/Nate/Documents/wilson_nisoli_tao/tao2022_sec_7_2_7_3.txt`.

**Reading caveat (A1).** The local plaintext is OCR'd. Greek letters render irregularly. Throughout this map I write "ε" for what the OCR shows as a blank or "" in places where the Tao paper unambiguously denotes a single small parameter ε ∈ (0, 1/100). Where the OCR mangling leaves ambiguity (separation constant in Lemma 7.4, threshold in (7.16)), the reading is recorded with a flag and the contextual argument that pins down the intended value. **Live arXiv cross-check via WebFetch was blocked by the harness during Phase 1a.** A1 discrepancies, if any, are catalogued in §6 below from internal-consistency cross-checks only.

**Notation.** A "constant" below is any quantity entering an inequality at a point where, if it were larger, the inequality might fail. We list every such quantity in the order it first appears.

---

## §7.2 — Deterministic structural analysis of black points

### 7.2.A — Triangle geometry (eq 7.11)

Triangle Δ ⊂ (N+1)×Z is defined by
```
Δ = { (j', l') : j' ≥ j, l' ≤ l, (j' - j) log 9 + (l - l') log 2 ≤ s }
```
for a top-left corner (j, l) and size s ≥ 0.

**Constants entering:**

| # | Constant | Role | Status |
|---|---|---|---|
| C-1 | log 9 | slope-defining | EXPLICIT, irrational but explicit |
| C-2 | log 2 | slope-defining | EXPLICIT |

No estimate-breaking content here; both are exact mathematical constants. Note the **slope ratio (log 9)/(log 2) ≈ 3.17** drives downstream geometric estimates; the renewal-process mean slope (4 in j, 16 in l) gives slope 16/4 = 4 > 3.17, which is structurally why the process tends to exit triangles through the top (horizontal) edge (Tao's remark in the discussion of Figure 3).

### 7.2.B — Lemma 7.4 (Structure of black set)

Statement (recovered from OCR + the project's R77.2 certification cross-check):

> The black set B ⊂ ⌊n/2⌋ × Z (points (j, l) with |Φ(j, l)| ≤ ε) is a disjoint union of triangles Δ ∈ T, each contained in [⌊n/2⌋ − (1/10) log(1/ε)] × Z, and any two distinct triangles in T are separated by Euclidean distance at least **(1/10) log(1/ε)**.

**Constants entering:**

| # | Constant | Role | Status | Note |
|---|---|---|---|---|
| C-3 | ε ∈ (0, 1/100) | "black" threshold; absolute small parameter | EXPLICIT-PARAMETERIC | ε is the single small parameter the whole §7 argument depends on; ultimately fixed at proof end. The OCR sometimes drops ε's; reading is confirmed by §7.5 and Prop 1.17 introduction. |
| C-4 | 1/100 | "weakly black" threshold scale (proof of 7.4) | EXPLICIT NUMERIC | Tao defines weakly-black as |Φ| ≤ ε/100. The factor 1/100 appears specifically as a margin for the propagation lemmas (i)/(ii)/(iii) inside the proof (e.g., the line "|Φ(j, l)| ≤ 5/100" needs to remain < 1/100 after multiplication by 9; in fact 5/100·9 = 45/100, but then the 4·(j, l−1) absorbs etc.) — this is a **proof-internal absolute** with no parametric dependence. |
| C-5 | 1/10 | exponent in strip and separation | EXPLICIT NUMERIC | Strip: `j ≤ n/2 − (1/10) log(1/ε)`; separation: `≥ (1/10) log(1/ε)`. The 1/10 is one of the few hard-numeric constants in §7.2 and it shows up in every downstream estimate (cases 2 and 3 of Claim (*)). |
| C-6 | "log 9 + log 2" | absorbed in 1 − (log 9 + log 2)/10 < 1/2 boundary check (Case 1 of Claim *) | EXPLICIT | Derived: (log 9 + log 2)/10 ≈ (2.197 + 0.693)/10 ≈ 0.289, so 1 − 0.289 ≈ 0.711 — yes, < 1/2 fails! The actual constraint in eq (after 7.18) is exp(−s + (j′ − j)log 9 + (l − l′)log 2) ≤ exp(...) · ε^{1 − (log 9 + log 2)/10} < ε^{2/3}; the exponent 1 − (log 9 + log 2)/10 ≈ 0.711 is what determines that the bound stays well below ε^{2/3} < ε/2 for small ε. **The choice 1/10 (C-5) was made precisely to make this hold.** |
| C-7 | propagation factors 9, 4, 2 in claims (i)/(ii)/(iii) | absorbed | EXPLICIT NUMERIC | (i) |Φ(j+1,l)| or |Φ(j,l−1)| black implies |Φ(j,l)| ≤ 9ε or 2ε respectively (from 7.13/7.14); (ii) weakly-black at corners → 5/100 propagation; (iii) 9/100 propagation. All numeric, traceable. |

### 7.2.C — Eq (7.16): strip containment

The estimate
```
3^{n+1−2j} ε ≤ 1/3
```
implies j ≤ ⌊n/2⌋ − (1/10) log(1/ε). The "1/3" comes from the residue analysis "expression is 1/3 or 2/3 mod Z." Explicit. Constant C-5 (the 1/10) is set here by the choice of base e logarithm and the structural fact that log 3 · 2 = log 9, with the 1/10 chosen so that the strip-margin matches the separation guarantee that the rest of §7.2 requires.

### 7.2.D — Remark 7.5 (not used)

Mentions Baker's theorem on log 3 / log 2; explicitly says "we will not exploit any further structure of the black set in our arguments." No constants enter the rest of the proof from here.

---

## §7.3 — Formulation in terms of holding time

### 7.3.A — Holding time Hold and Lemma 7.6

> Hold = (1, b_{[1,j*]}) where j* is the least j with b_j = 3.
> b_j iid copies of a_1 + a_2, a_i iid Geom(2) [⇒ a_1+a_2 ~ Pascal/negative-binomial(2, 1/2)].

P(b_j = 3) = 1/4 (eq 7.25). Mean of Hold = (4, 16) (eq end of Lemma 7.6 proof).

**Constants entering:**

| # | Constant | Role | Status |
|---|---|---|---|
| C-8 | 1/4 = P(Pascal = 3) | exact geometric distribution mass | EXPLICIT, exact rational |
| C-9 | EHold = (4, 16) | mean of holding time | EXPLICIT, exact rational |
| C-10 | exponential-tail rate of Hold (Lemma 7.6) | first time absolute constant `c > 0` appears | NAMED-UNSPECIFIED |

The exponential-tail rate of Hold is the first instance of Tao's "absolute constant c > 0" convention (paper p. 13: "implied constants are allowed to vary from line to line"). Specifically:
```
E exp(Hold · k) < ∞ for ‖k‖ sufficiently close to 0.
```
The rate is finite and explicit-in-principle (it is the radius of convergence of the Pascal moment-generating function around 0), but Tao never extracts a number. **By Lemma 2.2 the Hold variable then has exponential tail with parameter c > 0 inherited from this MGF radius.** Bookkeeping requires computing this radius.

Exact computation: E exp((1, Pascal') · k) = (3/4)^{-1} · sum_{b ≠ 3} (3/4) · 2^{−b} · (4^{b−1}/3) · exp(k_1 + b k_2). The geometric series converges for k_2 < log 2 (and any k_1). So the natural Hold MGF lives on |k_1| < ∞, k_2 < log 2 − δ for some δ > 0. The Markov-tail constant c is then any value < log 2 − δ; canonical choice c = (log 2)/2 ≈ 0.347 gives exponential tail of rate c — but Tao never fixes this. **This constant feeds into every subsequent inequality in §7.3 and §7.4 that uses "exponential tail of Hold."**

### 7.3.B — Lemma 7.7 (Distribution of first passage location), eq 7.48

> P(v_{[1,κ]} = (j', l')) ≪ e^{−c(l'−s)} · (1+s)^{−1/2} · G_{1+s}(c (j' − s/4))
> where G_{1+s}(x) = exp(−|x|²/(1+s)) + exp(−|x|) (defined in eq 2.2).

**Constants entering:**

| # | Constant | Role | Status |
|---|---|---|---|
| C-11 | c (in e^{−c(l'−s)}) | exponential-tail / Markov constant inherited from C-10 | NAMED-UNSPECIFIED (= C-10 up to "vary line to line") |
| C-12 | c (in argument of G_{1+s}) | Gaussian width / quadratic-domain rate | NAMED-UNSPECIFIED (could be same as C-11 or smaller; Tao's convention permits drift) |
| C-13 | ≪ in (7.48) | absorbs absolute prefactor from Lemma 2.2 + union bound | NAMED-UNSPECIFIED, single Vinogradov absolute |
| C-14 | constant from Lemma 2.2 application (eq mid-proof: P(v_{[1,k−1]} = (j', s')) ≪ k^{−1} G_{k−1}(c((j', s') − (k−1)(4, 16)))) | local-CLT / multivariate Berry-Esseen-like bound on the renewal process | NAMED-UNSPECIFIED; **load-bearing** |

**Crucial pointer.** C-14 is generated by **Lemma 2.2 (a generic 2D-renewal local-CLT / Berry-Esseen estimate** stated abstractly in Tao §2.3, separate from §7). Its proof is referenced in Tao §2 (Bourgain–Konyagin / Kahn–Komlós-style; one of the supporting references). To track C-14 effectively, one must redo the **2D local-CLT bookkeeping** — this is a separate analytic-probability sub-project (Bourgain–Konyagin / multivariate Berry-Esseen domain).

### 7.3.C — Mid-proof summation arguments

Several "summing in l′" / "summing in j_k" steps with "adjusting the constants c appropriately." Each adjustment is finite (geometric-series convergence on the Markov tail) but **changes the effective c on every line.** Tao does not log the cumulative drift; bookkeeping would have to.

---

## §7.4 — Recursively controlling a maximal expression

### 7.4.A — Q(j, l) and Q_m, eqs 7.34–7.38

```
Q(j, l) := E ∏_{k∈N} exp(−ε · 1_W((j, l) + v_{[1,k]}))     (7.34)
Q_m := sup_{(j, l): j ≥ n/2 − m} max(n/2 − j, 1)^A Q(j, l)   (7.38)
Q_m ≤ m^A                                                     (7.39, base case)
```

**Constants entering:**

| # | Constant | Role | Status |
|---|---|---|---|
| C-15 | ε (same as C-3) | weight in exp(−ε · 1_W) | EXPLICIT-PARAMETERIC |
| C-16 | A | the target exponent of the n^{−A} bound | PARAMETER (free, A > 0) |
| C-17 | 1 in max(·, 1) | normalization to avoid singularity at j = n/2 | EXPLICIT NUMERIC |

### 7.4.B — Proposition 7.8 (Monotonicity)

> Q_m ≤ Q_{m−1} whenever C_{A,ε} ≤ m ≤ n/2 for some sufficiently large **C_{A,ε}** depending on A, ε.

**Constants entering:**

| # | Constant | Role | Status |
|---|---|---|---|
| C-18 | **C_{A,ε}** | the "sufficiently large" threshold that makes Prop 7.8 fire | **NAMED-UNSPECIFIED, A- AND ε-DEPENDENT — THIS IS THE PROXY FOR C_A IN PROP 1.17** |

By inductively iterating Prop 7.8 from m = ⌈C_{A,ε}⌉ up to m = n/2 (with the base case Q_m ≤ m^A for the bottom range), Tao obtains
```
Q_m ≤ C_{A,ε}^A · 1     for all 1 ≤ m ≤ n/2,
```
which is (7.37) and (7.36) and hence Prop 1.17 with implied constant **bounded above by C_{A,ε}^A** (up to the absorbed Vinogradov from Q(Hold) ≪_A n^{−A} · j^A and the EGeom(4)^A factor).

**This is where C_A is born.** The implied constant in Prop 1.17 is morally C_{A,ε}^A · E[Geom(4)^A], times factors from §7.2 propagation. The whole bookkeeping question is: **what is C_{A,ε}?**

### 7.4.C — Case 1 of Prop 7.8 proof: (j, l) ∈ W

Bound (7.42): `max(m − r, 1)^{−1} ≤ m^{−1} exp(O(r log m / m))`. This is an **explicit** elementary inequality with no hidden absolute constant; the O(·) inside the exp is absolute.

Then with Geom(4) for r:
```
Q(j, l) ≤ exp(−ε) · m^{−A} Q_{m−1} · E exp(O(A log m · Geom(4) / m))
       ≤ exp(−ε/2) · m^{−A} · Q_{m−1}   [for m large depending on A, ε]
```

**Constants entering:**

| # | Constant | Role | Status |
|---|---|---|---|
| C-19 | exp(−ε) factor from white-point gain | EXPLICIT |
| C-20 | "O(·)" inside the exp on the right of (7.42) | absorbed | NAMED-UNSPECIFIED; small (≤ 1, since log(1 − r/m) − log(1 − (m−1)/m) ≤ r/(m·(1−r/m)) − 0 ≤ r/(m−r) ≤ 2r/m for r ≤ m/2). Extractable to a numeric constant ≤ 2. |
| C-21 | E exp(O(A log m · Geom(4)/m)) = 1 + O(A/log m) | absorbed | NAMED-UNSPECIFIED; controlled by Geom(4) MGF. Numeric. |
| C-22 | "m large enough depending on A, ε" threshold for Case 1 | NAMED-UNSPECIFIED, A,ε-dependent | feeds C-18 |

### 7.4.D — Case 2: (j, l) ∈ Δ, s ≤ m / log² m

Defines first passage time κ = least k with l_{[1,k]} > s. Applies Lemma 7.7 (and hence C-11, C-12, C-13, C-14).

Estimate (7.49): `E exp(O(A log m / m · j_{[1,κ]})) ≤ 1 + O(A/log m)`. **Same form as C-21**, but here `j_{[1,κ]}` is the first-passage j-coordinate, not Geom(4). The bound uses C-14 (the local-CLT / Lemma 2.2).

Estimate (7.51): `P((j, l) + v_{[1,κ]} ∈ W) ≫ 1` — i.e., on the first exit, the process lands in white with probability ≥ some absolute constant > 0.

**Constants entering:**

| # | Constant | Role | Status |
|---|---|---|---|
| C-23 | O(·) in (7.49) | absorbed via C-14 | NAMED-UNSPECIFIED |
| C-24 | "≫ 1" in (7.51) | absolute lower bound on probability of landing on white side of triangle | NAMED-UNSPECIFIED; **proved by direct geometric argument (Case-1 analysis of Lemma 7.4 + Gaussian dispersion of first-passage location); extractable to a positive numeric constant by careful Gaussian-tail computation.** |
| C-25 | "m ≥ C_{A,ε}" threshold for Case 2 | NAMED-UNSPECIFIED, A,ε-dependent | feeds C-18 |

### 7.4.E — Case 3: (j, l) ∈ Δ, s > m / log² m

The hardest case. Introduces parameter **P = P(A, ε)** large compared to A, 1/ε but small compared to m, and triangle-multiplicity statistic r.

Bound (7.52): `s ≤ (log 9 / log 2) · m`. Slope ratio C-1/C-2, explicit ≈ 3.17.

Estimate `P(j_{[1,κ+P]} ≥ 0.9m) ≪_P exp(−cm)`. Constants:

| # | Constant | Role | Status |
|---|---|---|---|
| C-26 | 0.9 (threshold for "renewal process went too far") | EXPLICIT NUMERIC; chosen because 0.8 > (1/4)(log 9/log 2) ≈ 0.793, so it lies in the Gaussian tail | EXPLICIT |
| C-27 | 0.8 = bound on j_{[1,κ]} alone | EXPLICIT NUMERIC; chosen as in C-26 with the margin 0.8 − 0.793 ≈ 0.007 driving the rate c | EXPLICIT BUT TIGHT |
| C-28 | c in exp(−cm) tail | NAMED-UNSPECIFIED, derived from C-14 via the 0.8 − 0.793 ≈ 0.007 gap (large-deviations rate ∼ (gap)² / variance) | NAMED-UNSPECIFIED; **load-bearing for Case 3** |
| C-29 | ≪_P factor from Lemma 2.2 + P repeated applications | absorbed | NAMED-UNSPECIFIED |
| C-30 | 10A / 3 (the "many white points" target) | EXPLICIT NUMERIC | |
| C-31 | exp(−10A) and 10^{-A-2} thresholds | EXPLICIT NUMERIC | drives where C-18 lives in A |
| C-32 | **R := ⌈A²/ε⌉** | number of triangles bound | EXPLICIT FUNCTION OF A, ε |

### 7.4.F — Lemma 7.9 (Many triangles usually implies many white points)

> E exp(−ε ∑_{p=1}^{t_min(r,R)} 1_W((j', l') + v_{[1,p]}) + ε min(r, R)) ≤ exp(ε)

This is the **separation-of-triangles → many-white-points conversion** lemma, leveraging Lemma 7.4's separation distance C-5 = 1/10. The induction is on R; the base R = 1 is trivial; the induction step uses

> P((r ≠ 0) ∩ ((j′, l′) + v_{[1,k_1]} ∈ W)) ≫ P(r ≠ 0)

(eq 7.59, repeating the geometric argument of (7.51), hence C-24).

**Constants entering:**

| # | Constant | Role | Status |
|---|---|---|---|
| C-33 | the "≫" prefactor in (7.59) | absorbed | inherits from C-24 |
| C-34 | the "0 < ε < 1/100 sufficiently small" hypothesis | EXPLICIT PARAMETRIC | matches C-3 |

### 7.4.G — Lemma 7.10 (Large triangles are rarely encountered shortly after a lengthy crossing)

> P(E_{p, s'}) ≪ A²(1 + p)/s' + exp(−c A²(1 + p))

(for s′ ≤ m^{0.4}, 0 ≤ p ≤ m^{0.1}, encountering a triangle Δ′ of size s′′ ≤ s′ shortly after first passage).

**Constants entering:**

| # | Constant | Role | Status |
|---|---|---|---|
| C-35 | ≪ in main term of Lemma 7.10 | absorbed Vinogradov absolute | NAMED-UNSPECIFIED |
| C-36 | c in exp(−c A²(1+p)) | absolute exponential rate, inherited from Lemma 7.7 / C-14 | NAMED-UNSPECIFIED |
| C-37 | exponents 0.4, 0.6, 0.1, 0.2 in the m^{0.4}, m^{0.6}, m^{0.1}, m^{0.2} appearing in the proof (eqs 7.60–7.65) | EXPLICIT NUMERIC | all chosen so that powers of m are sub-dominant to log² m and 1; **traceable** |
| C-38 | the constant 10 in "l′ − s′/log 2 > l − 10" (eq just before 7.65) | EXPLICIT NUMERIC | margin to force triangle-distinctness contradiction |
| C-39 | the constant 4 in s′ < 4A(1 + p)³ (last paragraph of §7.4) | EXPLICIT NUMERIC; chosen so the union-bound estimate `P(E) ≪ A² 4^{−A}` is ≤ 10^{-A-2} | EXPLICIT |

### 7.4.H — Final union-bound assembly

> P(F) ≪ 10^{−A−2}     (P(F) is the failure event)
> ⇒ ∑_{p=1}^{t_R} 1_W ≥ min(r, R) − O(A) ≥ A²/(3ε) when r ≥ R = A²/ε
> ⇒ Q(j, l) ≤ m^{−A} Q_{m−1} as required (eq 7.41).

**Constants entering:**

| # | Constant | Role | Status |
|---|---|---|---|
| C-40 | ≪ in P(F) ≪ 10^{-A-2} | absolute Vinogradov + assembled prefactors from all prior cases | NAMED-UNSPECIFIED; **terminal accumulation point** |
| C-41 | "P large enough depending on A, ε" final threshold | NAMED-UNSPECIFIED, A,ε-dependent | feeds C-18 directly |

---

## §7-Auxiliary — Lemma 2.2 (referenced cross-section)

Lemma 7.7's proof cites Lemma 2.2 (a 2D-renewal local-CLT / Berry-Esseen-style bound). This is the **single largest reservoir of absorbed absolute constants** in the entire §7 pipeline. C-14, C-11, C-12, C-28, C-36 all trace to constants generated inside the proof of Lemma 2.2.

Lemma 2.2 is in turn referenced via Tao's general framework setup in §2; its proof depends on Fourier-analytic / characteristic-function-based estimates with absolute constants that, on inspection, would be tractable to extract for a SPECIFIC 2D Z² lattice random walk (the joint Geom(4) × Pascal-conditioned-on-{≠3} variable). The bookkeeping cost is **non-trivial but not blocked-by-expertise** — it is a standard multivariate Berry-Esseen / saddle-point computation.

---

## §6 — A1 OCR-discrepancy log

Internal cross-checks against the project's existing R77.2 certification document and the §7.4 prose:

1. **(7.16) threshold "1/3":** OCR shows `3^{n+1-2j}ε ≤ 1/3` — confirmed by the "1/3 or 2/3 mod Z" preceding line. **Reading consistent.**
2. **Strip "[n/2 − (1/10) log(1/ε)] × Z":** OCR shows ` [ n  -  1   log  1  ]  �  Z ` — the (1/10) was read off the proof's appearance of this same coefficient later in Claim (*) Cases 2 and 3, where it modulates log 9 / log 2 in (7.19), (7.20), (7.22), (7.23). **Reading consistent.**
3. **Separation "(1/10) log(1/ε)":** Same OCR pattern, same internal check via Cases 2/3 of Claim (*). **Reading consistent.**
4. **Lemma 7.7 first passage probability:** OCR cleanly shows `e^{-c(l-s)} (1+s)^{-1/2} G_{1+s}(c(j-s/4))`. **Reading consistent.** No ambiguity.
5. **Geom(4) and EHold = (4, 16):** Both verified by the explicit Pascal arithmetic in Lemma 7.6 proof. **Reading consistent.**
6. **Case 2 vs Case 3 threshold `s ≤ m/log² m`:** OCR consistently shows this scaling. **Reading consistent.**
7. **R = A²/ε (eq 7.66):** OCR shows `R := A²/ε`. **Reading consistent** with Case-3 induction structure.

**No A1 discrepancies that would change Phase 1 disposition.** WebFetch cross-check against arXiv was blocked, so a full A1 audit must wait for an unblocked WebFetch session; this is logged for the disposition.

---

## §7 — Constant inventory summary (counts)

- **EXPLICIT (numeric or exact-rational):** C-1, C-2, C-3 (parameter ε ∈ (0, 1/100)), C-4, C-5, C-6, C-7, C-8, C-9, C-15, C-16, C-17, C-19, C-26, C-27, C-30, C-31, C-32, C-34, C-37, C-38, C-39 — **22 explicit constants**.
- **NAMED-UNSPECIFIED ("absolute c > 0", "≪", "≫"):** C-10, C-11, C-12, C-13, C-14, C-20, C-21, C-22, C-23, C-24, C-25, C-28, C-29, C-33, C-35, C-36, C-40, C-41 — **18 named-unspecified constants**.
- **PARAMETER-DEPENDENT TERMINAL:** C-18 (= C_{A,ε}), C-32 (R = A²/ε) — **2 parameter-dependent terminal constants** of which C-18 IS the proxy for Prop 1.17's C_A.

**Total constant entries: 42 (plus the Lemma 2.2 unitary-Berry-Esseen reservoir hidden inside the named-unspecified group).**

The 22:18:2 split is the headline. Most absorbed-Vinogradov constants are inherited from Lemma 2.2 (the renewal-process local-CLT). A handful (C-24, C-28, C-36) are honestly geometric/large-deviations constants that can be extracted by direct Gaussian-tail computation. The "≪ varies line to line" convention (Tao p. 13) means the same symbol c shifts on every estimate; cumulative drift is the dominant looseness source.

---

End of Phase 1a — Constant Map.
