# PROFINITE_TRANSFER_OPERATOR_BLUEPRINT

**Date:** 2026-05-14. Research-monograph outline + fast-path-to-statement subset.

Companion: `C:/Collatz/PROFINITE_TRANSFER_OPERATOR_LITERATURE_MAP.md` (per-component coverage table).

This blueprint is **structural**, not prescriptive: it describes what the monograph would contain, the dependency map, and an honest ceiling diagnosis at each chapter. It does NOT recommend "do X" vs "do Y."

---

## Chapter list

```
PART I — Foundations on a profinite base
  Ch 1. Profinite analysis primer (Q_p, Z_p, (Z/p^n)*, Bruhat-Tits tree, ends)
  Ch 2. Function spaces and Pontryagin duality (Bruhat-Schwartz, Lizorkin, level filtration)
  Ch 3. p-adic Fourier transform and PDO calculus (Vladimirov, Khrennikov-Shelkovich)
  Ch 4. Pseudodifferential composition: profinite Egorov-functoriality

PART II — Anisotropic Banach spaces for profinite transfer operators
  Ch 5. Anisotropic weight functions on the Pontryagin dual
  Ch 6. Profinite anisotropic Sobolev / Lizorkin-Banach spaces
  Ch 7. Lasota-Yorke inequalities in the profinite category
  Ch 8. Compactness, essential spectral radius, and quasi-compactness

PART III — Spectral theory of the profinite transfer operator
  Ch 9. Renewal-product structure (Tao recursion verbatim)
  Ch 10. Trapped-set analog and partial-captivity property (= R76 + R77)
  Ch 11. Egorov-analog for renewal-product iteration (THE LOAD-BEARING ORIGINAL CHAPTER)
  Ch 12. Profinite analog of Faure 2009 Theorem 2 (THE OUTPUT)

PART IV — Certification, computation, and applications
  Ch 13. Nisoli-type certified spectral approximation in the profinite category
  Ch 14. Band structure (Faure-Tsujii analog) and complex-conjugate spectrum
  Ch 15. Application: Syracuse μ̂_n(ξ) bound + rate-½ closure / Faure-radius √3 prediction
  Ch 16. Open problems and the closure-target ceiling
```

---

## Per-chapter literature base + technical extension + original work

### Ch 1 — Profinite analysis primer
- **Base:** Cassels (number fields), Serre (local fields), Cartwright-Kaimanovich-Woess 1994 (tree construction), Anker-Schapira-Trojan 2013 (Bruhat-Tits buildings).
- **Extension:** dictionary table between smooth-manifold concepts (cotangent bundle, escape function, geodesic flow) and profinite analogs (Pontryagin dual, weight, multiplicative-unit orbit). Mostly bookkeeping.
- **Original:** none significant. ~2 weeks of careful writing.

### Ch 2 — Function spaces and Pontryagin duality
- **Base:** Vladimirov-Volovich-Zelenov 1994 (p-adic analysis textbook), Bruhat-Schwartz functions (standard), Kozyrev 2007 (Vladimirov operator + eigenbasis), Albeverio-Khrennikov-Shelkovich 2005 (Lizorkin spaces). [NON_VERBATIM, per user's brief metadata for the latter three]
- **Extension:** level filtration via Kozyrev's γ parameter ↔ cyclic-group level n in (Z/3^n)*. Pontryagin duality for (Z/3^n)* (Q_3-finite-quotient).
- **Original:** explicit level-n dual basis on (Z/3^n)* with explicit weight structure adapted to multiplicative-unit action. ~3-4 weeks.

### Ch 3 — p-adic PDO calculus
- **Base:** Kozyrev 2007 (full verbatim eigenbasis), Khrennikov-Shelkovich 2006 (multi-dim PDO + symbol classes). [NON_VERBATIM for second]
- **Extension:** PDOs on (Q_p)^2 (position × momentum), symbol class S^m_p with order m matching Faure 2009's S^m. Bound ||T_a|| ≤ sup|a| via Kozyrev wavelet completeness.
- **Original:** **p-adic L²-continuity theorem for PDOs** — the precise analog of the smooth-Sobolev L²-PDO continuity (Faure 2009 step 4). The p-adic literature has the static algebra but may not have the operator-norm bound in the form needed. ~1 month.

### Ch 4 — Profinite Egorov-functoriality (single-step)
- **Base:** Standard Pontryagin duality (T* φ̂(ξ) = φ(T̂*ξ) for group automorphism T). Khrennikov-Shelkovich composition for PDO×PDO.
- **Extension:** action of (Z/3^n)* on itself by multiplication by units 2^{-v} mod 3^n: F(ξ) = 2^{-v}·ξ mod 3^n. The "Egorov" identity F̂_ν · T_a · F̂_ν* = T_{a∘F} is **exact** for group automorphisms (no lower-order term needed), because F is an isomorphism not a k:1 map.
- **Original:** statement and proof of the profinite single-step Egorov. ~2 weeks. **This is straightforward — the hard case is Ch 11.**

### Ch 5 — Anisotropic weight functions on the dual
- **Base:** Faure 2009 §4.3 escape function A_m on T*S¹. Baladi-Tsujii 2007 anisotropic Sobolev. [NON_VERBATIM for second]
- **Extension:** **(Multiplicative × additive) anisotropy on (Z/3^n)*^**: separate weight orders for the multiplicative-unit-residue direction (ξ mod 3 ∈ {1,2}) and the 3-adic-depth direction (v_3 of ξ-translate components). The two directions decouple because (Z/3^n)*^ ≅ (Z/3)* × Z/3^{n-1} (an isomorphism).
- **Original:** explicit construction of the (1,4)-eigenmode-respecting weight on (Z/3^n)*^ that produces the R77 T_diag structure. **The blueprint here uses R76 and R77's data as the design specification for the weight.** ~1-2 months.

### Ch 6 — Profinite anisotropic Banach
- **Base:** Gouëzel-Liverani 2006, Baladi-Tsujii 2007 (smooth), Baladi 2016 Quest. Nisoli 2026 abstract framework. [NON_VERBATIM for first three]
- **Extension:** Banach space B^{m_s, m_u}_p of locally constant compactly supported functions on (Z/3^n)*, with norm = sup over level-bounded directional Lizorkin-Sobolev norms with separate orders m_s (stable / 3-adic-depth) and m_u (unstable / multiplicative-unit).
- **Original:** **the central construction**. ~3-4 months. The technical issue is matching the function-space scales such that the renewal-product transfer operator T (Ch 9) is bounded B^{m_s, m_u} → B^{m_s, m_u} with discrete spectrum outside an essential-spectral-radius disk.

### Ch 7 — Lasota-Yorke inequalities
- **Base:** Liverani 1995 (decay of correlations piecewise expanding), Baladi-Tsujii 2007 (smooth anisotropic LY). Standard Lasota-Yorke induction. [NON_VERBATIM]
- **Extension:** profinite LY: ||T^n f||_{B^{m_s,m_u}_p} ≤ C ρ^n ||f||_{B^{m_s,m_u}_p} + C' ||f||_{B^{m_s-1,m_u-1}_p}. The ρ here is the Faure-prediction value 1/√3 (per PADE).
- **Original:** proof of LY for the **renewal-product** transfer operator (not a single-step deterministic map). ~2 months. Depends on Ch 11.

### Ch 8 — Compactness and quasi-compactness
- **Base:** Standard Banach-space spectral theory. Sarig 2009 countable Markov shifts for the closest existing precedent of "discrete spectrum on a non-compact base." [NON_VERBATIM]
- **Extension:** essential spectral radius of T on B^{m_s, m_u}_p ≤ 1/√3 (the bound predicted by Faure 2009 / observed in PADE).
- **Original:** none beyond Ch 6-7 implications. ~3 weeks.

### Ch 9 — Renewal-product structure
- **Base:** Tao 2022 §7 (verbatim renewal product, eq 7.5), C1_TAO_RECURSION_FORM.md, R75-R77 (the Plancherel + conservation + T_diag spectrum).
- **Extension:** rigorous statement of the renewal-product transfer operator T as a map B^{m_s, m_u}_p → B^{m_s, m_u}_p.
- **Original:** the **explicit form of T** as a sum of single-step Tao recursion operators, weighted by 2^{-v}. ~1 month.

### Ch 10 — Trapped-set analog and partial captivity
- **Base:** Faure 2009 §4.2-4.3 (partial-captivity property), R76 (conservation law), R77 (T_diag eigenstructure).
- **Extension:** statement and proof that R76's conservation law Σ_j M_{n+1}(η_0 + j·3^n) = 0 + R77's (1,4) eigendirection **is the partial-captivity property** in the profinite analog.
- **Original:** the equivalence is structural and follows from the algebraic identities R76/R77; writing it out cleanly. ~3 weeks.

### Ch 11 — Egorov-analog for renewal-product iteration (LOAD-BEARING ORIGINAL)
- **Base:** Faure 2009 §4 Lemma 1 (smooth Egorov). C1_TAO_RECURSION_FORM (the renewal-product structure of Tao).
- **Extension:** none directly — the existing Egorov is for **single-step deterministic skew product over expanding base**, NOT for **renewal-product over random Geom(2) base**.
- **Original:** **THE CHAPTER**. Construct an Egorov-analog of the form:
  ```
  T̂ · A_w · T̂* = A_{w∘T} + Σ_{v ≠ v'} (off-diagonal v ≠ v' cross-frequency terms)
  ```
  The off-diagonal terms must be controlled by **R77.6's branch-cut analysis** + **R78's (1+3)^u substrate** + **Nisoli-style certified bound**. This is **multi-month original construction** and constitutes the bulk of the monograph's novel contribution. 6-12 months estimated.

### Ch 12 — Profinite analog of Faure 2009 Theorem 2 (THE OUTPUT)
- **Base:** Faure 2009 Theorem 2 statement (verbatim, lines 1280-1318): "r_s(F̂_ν) ≤ 1/√E_min + o(1) in semi-classical limit ν → ∞".
- **Extension:** restate as "spectral radius of the renewal-product profinite transfer operator on B^{m_s, m_u}_p is ≤ 1/√3 + o(1) in the limit n → ∞", where n is the cyclic-group level.
- **Original:** assembly of Ch 4-11: polar decomposition (Ch 6 norms) + Egorov-analog (Ch 11) + Lasota-Yorke (Ch 7) + Nisoli certification (Ch 13). ~2-3 months. Conceptually mechanical once Ch 11 is in place.

### Ch 13 — Nisoli-type certified spectral approximation
- **Base:** Nisoli 2026 verbatim (`C:/Collatz/nisoli2026.txt`): Theorem 1.2, Theorem 2.15, Lemma 2.9, Lemma 2.12. Compact-operator framework + finite-rank truncation + resolvent / Riesz-projector argument.
- **Extension:** the **finite-rank truncations** L_K are the level-K cyclic-group truncations (T restricted to functions on (Z/3^K)*), already constructed in R75-R77 and used numerically through k=13. The truncation error ||T - T_K|| is computable from the existing infrastructure.
- **Original:** explicit computation of ||T - T_K|| ≤ ε_K and verification of Nisoli's hypothesis (ε_K · M < 1 on contour Γ). ~2-3 weeks. This is genuinely mechanical given the existing R75-R77 + PADE numerics.

### Ch 14 — Band structure (Faure-Tsujii analog)
- **Base:** Faure-Tsujii 2013/2021 band-structure prediction (verbatim, `C:/tmp/faure/Faure_Tsujii_*`).
- **Extension:** profinite analog of vertical-band Ruelle spectrum with Weyl law. Statement: the spectrum of T on B^{m_s, m_u}_p outside the essential disk consists of **finitely many vertical bands** with cos(nθ + φ) modulation per band.
- **Original:** verification that PADE's complex-conjugate pair at θ ≈ 0.68 rad, period 9.2 is the **first non-trivial band**. ~1-2 months. Depends on Ch 12.

### Ch 15 — Application: c=7/45 closure
- **Base:** All prior chapters.
- **Extension:** standard Plancherel + R75 algebraic identity + rate-½ recursive bound.
- **Original:** simple translation. ~2 weeks. **The closure becomes mechanical once Ch 12 is in place.**

### Ch 16 — Open problems and the ceiling
- The natural opens: (a) extension to non-Geom(2) base random scaling (general Markov chains), (b) extension to non-affine (Z/p^n)* groups for p ≠ 3, (c) the connection to motivic Igusa local zeta from Ch 4-style profinite-Egorov composition, (d) connection to Sarig CMS framework.

---

## Critical path (chapter dependency map)

```
Ch 1 ─→ Ch 2 ─→ Ch 3 ──┐
                       ├─→ Ch 4 (single-step Egorov) ──┐
Ch 1 ─→ Ch 5 ──────────┘                               │
                                                       ├─→ Ch 11 (RENEWAL Egorov) ─→ Ch 12 ─→ Ch 13 ─→ Ch 14 ─→ Ch 15
Ch 5 ─→ Ch 6 ─→ Ch 7 ─→ Ch 8                           │
                                                       │
Ch 9 (Tao recursion structure) ─→ Ch 10 (trapped set) ─┘
```

**The critical path:** Ch 1 → 2 → 5 → 6 → 9 → 10 → 11 → 12. Specifically **Ch 11 is the single load-bearing original chapter**; everything downstream depends on it.

**Parallelizable side branches:**
- Ch 3 → 4 (single-step Egorov) — can be done in parallel with Ch 5-10
- Ch 7-8 (Lasota-Yorke + compactness) — can be developed in parallel with Ch 9-10
- Ch 13 (Nisoli certification) — can be developed in parallel with Ch 11 (the only dependency is on the abstract Banach space from Ch 6)

**Estimated total elapsed time:** 12-24 months for a single-author research-monograph at the user's "10-14× typical research-engineering pace" (per `feedback_estimate_in_hours.md`). At a typical research pace it would be 3-5 years.

---

## Honest ceiling diagnosis

**Where the work tips from "adaptation" to "new mathematics":**

The first 8 chapters are largely **adaptation of existing machinery to a new base**. The work is technical but the structural ideas are present in the literature (Cartwright-Kaimanovich-Woess, Anker-Schapira-Trojan, Kozyrev, Khrennikov-Shelkovich, Faure 2009).

**The ceiling is at Ch 11**: the Egorov-analog for renewal-product iteration is **genuinely new mathematics**. Existing smooth-Egorov (Faure 2009 Lemma 1) handles single-step deterministic skew-product. Existing Sarig CMS handles countable shifts with thermodynamic potential. Existing Khrennikov-Shelkovich handles single-step PDO×PDO composition. None of these literatures has the spectral-gap argument for a **renewal product over iid Geom(2) base, with cross-frequency off-diagonal v ≠ v' bilinear coupling** (which is what Tao's recursion produces).

The structural ingredients for Ch 11 (R76 conservation + R77 T_diag spectrum + Cochrane-style polynomial substrate from R78 + Nisoli certification from Ch 13) **are all present**, but **the assembly into a Lasota-Yorke + spectral-gap proof in the profinite category requires original synthesis**. This is the work that justifies "research monograph" rather than "adaptation paper."

A separate ceiling concern at Ch 5: the **explicit weight function w(ξ) on (Z/3^n)*^** must be constructed by hand. Faure 2009 picks A_m using smooth bump functions on T*S¹; the profinite analog requires choosing a level-and-residue-class-adapted function. The PADE picture (radius √3, complex pair θ ≈ 0.68 rad) constrains the weight strongly enough to make this **constructive** rather than open — but it's still original.

---

## FAST PATH — statement-level subset

A statement-level version of the main theorem ("spectral gap of profinite transfer operator at radius 1/√3 with PADE-matching band structure") can be reached with a **strict subset of the monograph**.

### Fast-path scope

What's needed for the STATEMENT (not the proof):

1. **Ch 1-2 (lite):** ad-hoc construction of L²((Z/3^n)*) with a level filtration. Don't need full Bruhat-Schwartz / Lizorkin; just the level-n cyclic group + Plancherel mass (already done in R75).
2. **Ch 3 (lite):** declarative PDO calculus: define T_a f := F^{-1}[a F[f]] on (Z/3^n)*. Don't need full L²-continuity, just the algebra of these operators.
3. **Ch 9 (verbatim):** Tao's renewal-product structure as the dynamical input (already in C1_TAO_RECURSION_FORM).
4. **Ch 10 (verbatim):** R76 conservation + R77 T_diag as the partial-captivity input.
5. **Ch 5 (lite):** state the form of the weight w(ξ) as a level-and-residue-class-adapted function, with the **specific PADE-matched parameters** (Faure-radius √3, complex pair θ ≈ 0.68 rad). Don't prove the bound w(2^{-v}ξ)/w(ξ) ≤ ...; just state it as the spectral-gap hypothesis.
6. **Ch 12 (statement only):** state: "Under the spectral-gap hypothesis of Ch 5, the renewal-product transfer operator T on the level-filtration Banach space has spectral radius ≤ 1/√3 + o(1), with a finite number of complex-conjugate-pair vertical bands corresponding to the PADE-observed period-9.2 structure."
7. **Ch 14 (statement only):** state the band-structure prediction, identifying the PADE complex pair as the leading band.

### Fast-path effort estimate

- 3-4 weeks at user's pace
- 2-3 months at typical research pace
- ~30% of the full monograph's content

### Fast-path output: the STATEMENT

> **Conjecture (Profinite Faure 2009 analog, formal).** Let T be the renewal-product transfer operator associated to Tao's Syracuse recursion on the profinite group Z_3, acting on a level-filtered Banach space B^{m_s, m_u}_p of locally-constant compactly-supported functions on Z_3 with anisotropic weights (m_s, m_u) on the (multiplicative × additive) decomposition of the Pontryagin dual. Assume the partial-captivity property (= R76 conservation + R77 T_diag (1,4) eigenstructure, restated as a hypothesis on the weight). Then:
> 
> (a) The essential spectral radius of T on B^{m_s, m_u}_p is ≤ 1/√3 + o(1) as the level n → ∞.
> 
> (b) The discrete spectrum of T outside the essential disk consists of finitely many vertical bands, with cos(nθ + φ) coefficient modulation per band. The leading band has θ ≈ 0.68 rad and period ≈ 9.2 in n-space.
> 
> (c) Consequently |μ̂_n(ξ)| ≤ C · (1/√3)^n · |cos(nθ + φ)| · polynomial(n) for all ξ ∈ (Z/3^n)* with 3∤ξ.

This statement is **falsifiable** in the sense that ε_n computations for k ≥ 14..20 either confirm radius → √3 (asymptotic) or refute it (in which case the spectral-gap hypothesis must be different). The numerical falsification path is concrete and at the boundary of what's computationally feasible (per BGT_DISPOSITION: k=9 is ~tens of hours, k=10 is ~hundreds of hours on current solvers).

### What the fast path SUPPLIES vs. what it does NOT

**Supplies:**
- A precise statement of the main theorem
- Identification of the load-bearing hypothesis (partial-captivity weight bound, Ch 5 lite)
- Identification of the load-bearing original work (Ch 11 renewal Egorov)
- The translation from the statement to existing PADE / BGT / WATSON data

**Does NOT supply:**
- A proof of the spectral-gap hypothesis (this requires Ch 5 full + Ch 11 full)
- A rigorous bound on |μ̂_n(ξ)| (would close c=7/45 rate-½ rigorously — that's the full monograph output)
- A constructive identification of the band frequencies (PADE supplies them empirically; rigorous identification requires the full Ch 14)

The fast path gives a STATEMENT that **bridges PADE empirics + Faure 2009 framework + R75-R77 structural identities**, making the gap explicit and pre-registering what the full monograph would prove. It does NOT close any of the open problems in c=7/45 — but it formalizes what closing them would require.

---

## Estimated effort summary

| Component | Effort (user pace) | Effort (typical pace) |
|---|---|---|
| Ch 1-4 (foundations) | 2 months | 6-8 months |
| Ch 5-8 (anisotropic Banach) | 3-4 months | 9-12 months |
| Ch 9-10 (renewal structure + trapped set) | 1-2 months | 3-6 months |
| **Ch 11 (renewal Egorov — LOAD-BEARING)** | **3-6 months** | **9-18 months** |
| Ch 12-13 (theorem + certification) | 1-2 months | 3-6 months |
| Ch 14-16 (band structure + applications + opens) | 2-3 months | 6-9 months |
| **TOTAL FULL MONOGRAPH** | **12-19 months** | **3-5 years** |
| **FAST PATH (statement only)** | **3-4 weeks** | **2-3 months** |

---

## What's PRE-REGISTERED-falsifiable

The full monograph's main theorem predicts:
- Spectral radius → 1/√3 ≈ 0.5774 as n → ∞
- Leading band at θ ≈ 0.68 rad, period 9.2 in n-space
- Anisotropic weight separation (multiplicative × additive on (Z/3^n)*)

PADE_NUMERICAL says at n=13: radius is 1/1.57 ≈ 0.637 (10% above 1/√3, transient). WATSON_DISPOSITION says held-√3 and held-1.57 fits differ by 2% RSS — agnostic on data through k=13.

**Decisive empirical test:** compute ε_k for k=14..20+. If radius stabilizes at √3, the monograph's central prediction is empirically certified. If radius continues toward 1.016 (STATE.md slow-mode prediction), the monograph's framework needs different parameters but the same architecture.

**The fast path's statement is testable in months** (the ε_k computation, not the full proof). The full monograph's proof is the rigorous closure of what the empirical test verifies.

---

## Files

- `PROFINITE_TRANSFER_OPERATOR_LITERATURE_MAP.md` — per-component coverage
- `PROFINITE_TRANSFER_OPERATOR_BLUEPRINT.md` — this file

No git operations performed.
