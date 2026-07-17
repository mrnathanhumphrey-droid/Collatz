# IGUSA_DISPOSITION — Igusa local zeta on R78 (1+3)^u substrate, headline

**Date:** 2026-05-14. Probe IGUSA. Working dir: C:/Collatz/. Mode E (verbatim from 10 PDFs at C:/Users/Nate/OneDrive/Documents/igusa_local_zeta/pdfs/, extracted via pypdf to C:/tmp/igusa/).

Pre-registration: `C:/Collatz/IGUSA_PRE_REGISTRATION.md`.

---

## Headline

**NO_SELECTED. Net disposition: 10 NO_FIT (categorically dominated by a single structural barrier).**

A **THIRD INDEPENDENT, CATEGORICAL STRUCTURAL BARRIER** to log_3(2) as an Igusa pole emerged at Phase 2:

> **All real-parts of Igusa local zeta poles are NEGATIVE RATIONAL** (Igusa rationality Thm 1.3.2 + Kashiwara theorem on b-function roots + Monodromy Conjecture requiring roots of unity).
>
> The probe target **s = log_3(2) ≈ +0.6309 is POSITIVE and IRRATIONAL.** It cannot be an Igusa pole real-part for any polynomial f ∈ Q[x_1,…,x_n] over any local field Q_p.

This barrier is independent of substrate choice — it is a property of the *category of objects* "Igusa local zeta functions of algebraic polynomials". The R77.6 numerical anchor s = log_3(2) lives in a fundamentally different category (irrational positive real, characteristic of branch-cut singularities or transcendental analytic structures), not in the Igusa-pole category (negative rational, characteristic of algebraic-polynomial Mellin transforms).

The **ninth category-of-object barrier** in the systematic obstruction map for c = 7/45 closure.

---

## Pre-registered probability vs. realized outcome

| Outcome | Pre-Phase-1 | Realized |
|---|---|---|
| SELECTED | 20-25% | NOT |
| PARTIAL | 30% | NOT |
| NO_FIT | 30% | **REALIZED — dominantly** |
| BLOCKER | 10% | NOT |
| MODE_H_CIRCULAR | 5-10% | NOT |

The pre-registration explicitly flagged a "structural pre-warning" (page-bottom note) about R78's D=0 implying g(u) is a 3-adic unit, hence trivial Igusa zeta. That warning was confirmed at Phase 2A. A *second* structural barrier (positive/irrational target vs negative/rational Igusa poles) emerged separately at Phase 2B and was confirmed independently at Phase 2D (Monodromy) and Phase 2H (b-function).

---

## Summary table

| Code | Theorem | Phase 0 | Phase 1 hyp check | Phase 2 conclusion | Phase 3 substrate match | Disposition |
|---|---|---|---|---|---|---|
| A | Igusa rationality (Denef Bourbaki Thm 1.3.2) | EXTRACTED | SATISFIED (any nonzero polynomial) | Substrate 1: trivial Z=1 (no poles) by R78 D=0. Substrate 2: pole at s=-1 | POLE_LIST_MISMATCH / POLYNOMIAL_FORM_TRIVIAL | **NO_FIT** |
| B | Denef-Hoornaert Newton polyhedron (Thm 0.27 / Thm 2.9) | EXTRACTED | SATISFIED for univariate (with non-degeneracy vacuous) | Single facet v=1, m(v)=1, σ(v)=1 → s=-1 only | POLE_LIST_MISMATCH — log_3(2) is positive irrational, Igusa poles negative rational | **NO_FIT (categorical)** |
| C | Bories-Veys non-degenerated surfaces (Thm 0.12) | EXTRACTED | FAILED h_DIMENSION (n=3 required, substrate n=1) | N/A | N/A | **NO_FIT (dimension)** |
| D | Monodromy Conjecture (Veys Conj 2.12) | EXTRACTED | Vacuous (downstream of pole identification) | Forces e^{2πi Re(s_0)} ∈ roots of unity → Re(s_0) ∈ Q | log_3(2) irrational, RULED OUT BY MONODROMY | **NO_FIT (categorical, 2nd barrier)** |
| E | Motivic Igusa (Denef-Loeser Thm 5.5) | EXTRACTED | Inherits A's hypotheses | Same pole real-parts as A: {-ν_i/N_i} | Same as A | **NO_FIT** |
| F | Bories-Cluckers polynomial mapping (Thm 2.5) | EXTRACTED | FAILED h_POLY_MAPPING_FORM (substrate is single poly, not mapping/ideal pair) | 1D reduction = candidate B | Same as B | **NO_FIT** |
| G | Stationary phase / Watson on Igusa (Cor 1.4.5) | EXTRACTED | SATISFIED vacuously (C_g empty for substrate 1) | Asymptotic E_Φ(z) depends on Igusa poles — no poles for substrate 1; positive pole would imply growing asymptotic violating boundedness | POLE_INCOMPATIBLE_WITH_BOUNDEDNESS | **NO_FIT (reinforces categorical barrier)** |
| H | Bernstein-Sato b-function (Veys Thm 2.2-2.6) | EXTRACTED | SATISFIED | b_g(s) = s+1 for generic univariate g. Igusa poles ⊂ {root(b_f) - k : k ≥ 0} ⊂ Q_{<0} (Kashiwara: roots of b_f are negative rational) | log_3(2) irrational AND positive — RULED OUT BY b-FUNCTION | **NO_FIT (categorical, 3rd barrier)** |
| I | Monomial / linear-form Igusa (Denef Bourbaki §1.3) | EXTRACTED | g(u) does not reduce to monomial/linear after change of variables (constant mod p) | Monomial poles: -1/a_j ∈ negative rationals; linear: -n | log_3(2) ∉ {-1/a, -n} | **NO_FIT (categorical)** |
| J | (1+p)^u-specific Igusa | OPEN_SEARCH | None exists in corpus. Hypothetical reformulation collapses to candidate I (linear monomial) | Pole at s=-1 again | log_3(2) still not in pole list | **NO_FIT (corpus negative)** |

Total candidates with extractable Phase-0 statements: **10/10**. No BLOCKER.
- 0 SELECTED
- 0 PARTIAL
- 10 NO_FIT (all dominated by the three categorical barriers: trivial substrate, positive vs negative real part, irrational vs rational real part)

No surfaced candidate K beyond the search-for-(1+p)^u-specific theorem in candidate J (none exists).

---

## Final disposition: **NO_FIT (categorically dominant)**

The Igusa local zeta substrate is closed as a closure route for c = 7/45.

**Three independent categorical barriers, each by itself sufficient:**

### Barrier 1: R78's D=0 makes the natural univariate substrate trivial

The polynomial g(u) = c · Σ_{k=0}^r C(u,k) p^k − p²·m·u, viewed as element of Z_p[u], has g(u) ≡ c (mod 3) where c ∈ (Z/3^{r+1})* is a unit. Therefore |g(u)|_3 = 1 uniformly on Z_3, and

> Z(s, g, 3) = ∫_{Z_3} |g(u)|_3^s du = ∫_{Z_3} 1 · du = 1.

The Igusa local zeta is the constant function 1, with no poles. **The natural Syracuse polynomialization gives no Igusa-pole information at all.**

Cochrane-Pinner Postnikov-style alternative (substrate 2, R78 Feature 4/5 cubic phase P_a(s)) gives Igusa zeta with a single pole at s = -1 — the standard p-adic-log-canonical-threshold pole — not at log_3(2).

### Barrier 2: log_3(2) is POSITIVE; Igusa poles have NEGATIVE real part

By Igusa's rationality theorem (Denef Bourbaki Thm 1.3.2), poles are of form s = -ν_i/N_i + 2πik/(N_i log q) with (N_i, ν_i) numerical data from an embedded resolution. Both ν_i ≥ 1 and N_i ≥ 1 (positive integers from divisor multiplicities), so **Re(s_0) = -ν_i/N_i ≤ -1/N_max < 0**.

By convergence: ∫|f|^s converges for Re(s) > 0, so any pole must have Re(s) ≤ 0.

**log_3(2) ≈ +0.631 > 0 is OUTSIDE the Igusa pole region.**

### Barrier 3: log_3(2) is IRRATIONAL; Igusa pole real-parts are RATIONAL

By the Monodromy Conjecture (Veys Conj 2.12, established for surfaces by Bories-Veys, broadly conjectural elsewhere): if s_0 is an Igusa pole, then exp(2πi Re(s_0)) is a monodromy eigenvalue of f, hence a **root of unity** (Veys Prop 2.5(1)). Therefore Re(s_0) ∈ Q.

Even unconditionally (without invoking Monodromy): Igusa rationality gives Z(s, f, p) as rational in p^{-s}. Poles of rational functions in p^{-s} (over the appropriate algebraic closure) have **rational real parts** (since the polynomial in p^{-s} has algebraic-number roots whose magnitudes log/log p are rational over Q).

The Bernstein-Sato confirmation (Veys Thm 2.6 + Kashiwara theorem): roots of b_f are negative rational, and Igusa poles ⊂ {root(b_f) − k : k ∈ Z_{≥0}} ⊂ Q_{<0}.

**log_3(2), log_3(45/43), log_3(1/0.984) are all irrational** (since their arguments 2, 45/43, 1/0.984 are not integer powers of 3). **CATEGORICAL BAR.**

---

## What category of theorem is missing

For ANY theorem to produce a pole at s = log_3(2) (or similar irrational positive location), the theorem must operate on an object whose pole locations are not constrained to (negative rational ⊂ real). Candidate object categories include:

1. **Branch-cut singularities** of generating functions f̃(z) = Σ ε_k z^k where the singularity at z = z_0 corresponds to s = log_p(z_0) for any positive real z_0. This is exactly the R77.6 reading and the Padé-numerical reading. Branch cuts are NOT poles in the Tate/Igusa Mellin sense; they're a different analytic category.

2. **Spectral gap of a transfer operator on a partially expanding map** — Faure-style semiclassical spectral theory (R77.5/R77.6 already pointed at this). The spectral radius can be any positive real, including irrational; the corresponding "Mellin pole" is at s = -log_p(spectral_radius), which can be irrational positive (with sign flip relative to Igusa).

3. **Asymptotic expansion of incomplete exponential sums via Erdős-Ko-Rado / Burgess-style bounds** — gives non-polynomial-in-A asymptotics with irrational exponents tied to character-sum critical exponents.

4. **Stochastic-iteration Mellin transform on a profinite multiplicative group** with controlled deviation envelope — the hybrid object identified in ADELIC_DISPOSITION's "missing category".

**None of these are Igusa.** The Igusa category is fundamentally constrained to rational pole real-parts.

---

## Strategic position

Pre-IGUSA: 8 categorical barriers mapped (5-probe Fourier-decay, Cluster 1 Cochrane, Cluster 2 BMP, Bruhat-Tits/BKL, Tauberian, Furstenberg-Guivarc'h, BGT, ADELIC).

Post-IGUSA: **9 categorical barriers mapped.** The Igusa local zeta substrate is closed.

The IGUSA probe is structurally different from the prior 8 in one important way: **the categorical barrier here is identified MORE SHARPLY** than in any prior probe. Specifically:

- Prior probes (BGT, ADELIC, Tauberian, Cluster 1/2, etc.) found "no theorem in this corpus fires", which is corpus-specific.
- IGUSA finds: **no theorem in this CATEGORY can fire**, period — by the rationality / monodromy / b-function structural argument. The barrier is intrinsic to the Igusa category of objects, not to the specific corpus.

This makes IGUSA's NO_FIT verdict **strong and final** for the Igusa-school approach. There's no "missing PDF" that could change it.

### Sharpened picture of the open question

The structural identification of log_3(2) as the R77.6 anchor was a STRONG signal — "irrational positive real, characteristic of branch-cut not pole, characteristic of transfer-operator spectral radius not algebraic-polynomial Mellin pole". The IGUSA probe sharpens this:

> **The Syracuse μ_n asymptotic rate is encoded in an object whose natural singularity is a BRANCH CUT (or transfer-operator spectral gap), NOT an algebraic-polynomial Mellin POLE.**

This is consistent with:
- R77.5/R77.6 branch-cut reading (now sub-leading per Padé-numerical, but structurally branch-cut not pole).
- PADE_NUMERICAL_DISPOSITION's leading singularity at |z| ≈ 1.5..1.7 with complex-conjugate-pair structure (consistent with transfer-operator-spectral-gap origin, not Igusa pole).
- T_LEAD_CORRECTED's discrete eigenvalue 43/45 over Q (an algebraic eigenvalue of an operator — also NOT an Igusa pole; eigenvalues are different objects).

**The closure target lives at the analytic operator-theoretic side, not the algebraic-polynomial Mellin side.** Faure 2009 semiclassical spectral gap is now the **only remaining categorically-distinct route** in the prior probe map.

---

## SECONDARY ROUTING

Per pre-registration (priority-ordered):

1. **Faure 2009 semiclassical spectral gap on partially expanding maps.** Categorically distinct from Igusa (operator-theoretic, not algebraic-polynomial Mellin). Operates on the transfer operator of Tao recursion as a partially-expanding-map iteration. Spectral radius can be irrational, real, or complex; produces exactly the kind of branch-cut singularity R77.6 was diagnosing. **PRIORITY: HIGH — now the sole remaining categorically-distinct probe.**

2. **Watson lemma / saddle-point on R78/R79 bilinear off-diagonal sum.** Operates closer to chain-side. Independent of Igusa. May give saddle-point asymptotic for the bilinear off-diagonal sum where the k=7 third-mode contribution lives. Doesn't directly close c=7/45 but tightens structural picture. **PRIORITY: MODERATE.**

3. **Multi-singularity Tauberian extension (Flajolet-Sedgewick Ch. VI §VI.4-VI.5).** Pre-flagged in BGT/Tauberian dispositions. Operates on the *sequence* ε_k directly via generating-function transforms with multi-singularity structure. Categorically distinct from Igusa AND from Faure. **PRIORITY: MODERATE.**

4. **Stochastic-iteration Mellin transform on profinite multiplicative group** (the "missing category" identified in ADELIC_DISPOSITION). Open construction problem; not yet a theorem in any corpus. **PRIORITY: LOW (theoretical construction, no existing theorem).**

5. **Direct computation of Z(s, g, 3) via Mainfile thesis rationality proof** (last resort, now MOOT given the three categorical barriers): even if executed exactly, would produce only negative-rational poles, none of which match the target. **PRIORITY: NONE (moot).**

### Recommended top-priority secondary route: **Faure 2009 semiclassical spectral gap**

Rationale (sharpened by IGUSA findings):

- **Faure operates on the right CATEGORY OF OBJECTS for log_3(2):** spectral radii of transfer operators can be any positive real, including irrational. The Igusa category (negative rational pole real-parts) is now CATEGORICALLY EXCLUDED.
- **Faure produces a branch-cut not a pole:** matches R77.6's diagnostic exactly.
- **Faure addresses the Tao recursion transfer operator directly:** chain-side, the operator μ̂_n → μ̂_{n+1} is partially-expanding (the 2-adic Geom(2) factor expands; the 3-adic level shift is neutral). This is exactly Faure's setup.
- **Faure was already pre-flagged in BOTH BGT_DISPOSITION and ADELIC_DISPOSITION as priority-2 secondary routing.** IGUSA's NO_FIT elevates it to priority-1.

The Faure semiclassical spectral gap probe is the natural next step.

---

## Surprises in the inputs

### Surprise 1: R78's D=0 finding is the SAME structural obstruction that closes IGUSA

The R78 Cochrane-Theorem-2 negative result ("g(u) mod 3 is constant; D=0; complete-sum trivially vanishes; partial-sum doesn't inherit the trivial vanishing"). This is identically the structural fact that closes IGUSA candidate A:

g(u) mod 3 constant ⇒ |g(u)|_3 ≡ 1 ⇒ Z(s, g, 3) = 1.

Same algebraic fact (g(u) ≡ c mod 3) reads as a NEGATIVE result for both Cochrane Theorem 2 (no character-sum saving) AND Igusa local zeta (trivial integral). The R78 obstruction propagates **all the way through Igusa** because the underlying issue is the same: 4 = 1+3 in Z_3 makes (1+3)^u behave 3-adically as 1 + O(3), which is a unit, not a vanishing-or-singular polynomial.

### Surprise 2: The categorical barrier is SHARPER than expected

Pre-registration estimated 25-30% NO_FIT, 30% PARTIAL. Realized: NO_FIT dominant due to three INDEPENDENT structural barriers, each sufficient by itself. The expected "PARTIAL — pole exists but at wrong location" outcome did NOT materialize: there's no Igusa theorem that even *could* produce log_3(2) as a pole, because positive irrational ∉ Igusa pole locus categorically.

This is a sharper barrier than the typical hypothesis-mismatch barrier of prior probes.

### Surprise 3: R77.6 branch-cut at z=2 is REFUTED at n=13 (per PADE_NUMERICAL_DISPOSITION)

The original numerical anchor s = log_3(2) ≈ 0.631 was derived from R77.6's z=2 branch-cut reading. PADE_NUMERICAL_DISPOSITION (n ≤ 13) shows Hadamard radius is **1.57** at n=13, trending inward toward 1.046 or 1.016. The "z = 2" anchor was a transient fingerprint at n=2..6; the actual leading singularity sits closer to z = 1.

This **further weakens** the case for log_3(2) being the relevant pole location, but the categorical IGUSA barriers (positive/irrational/Bernstein-Sato/Monodromy) close the Igusa route regardless of where the actual singularity sits, as long as it's positive and irrational (which 1.57, 1.046, 1.016 all are: log_3(1.57)≈0.41, log_3(1.046)≈0.041, log_3(1.016)≈0.014 — all positive irrationals).

### Surprise 4: IGUSA's NO_FIT is FINAL in a way prior probes' NO_FIT was not

Prior probes (BGT, ADELIC, etc.) left open "maybe a paper not in corpus would fire". IGUSA's NO_FIT is unconditional: no theorem of any kind in the Igusa-local-zeta-of-algebraic-polynomial category can produce an irrational positive pole. This is a property of the category, not a property of the corpus. No future corpus extension can change it.

The only way to get an irrational positive Mellin pole is to **leave the Igusa category** (use transfer-operator spectra, or branch cuts of non-algebraic generating functions, or transcendental-Mellin objects).

---

## What category of theorem WOULD fit

The closure target requires:

1. A theorem operating on a sequence ε_k or a generating function f̃(z) = Σ ε_k z^k or an iteration operator T whose spectrum encodes the rate.
2. Whose natural singularity / pole / spectral-radius locations are POSITIVE IRRATIONAL REALS (or complex with positive irrational modulus), not constrained to negative rationals.
3. That accepts non-symmetric, forward-only, partially-expanding dynamics (Tao recursion is exactly this).
4. That produces explicit asymptotic expansion convertible to |μ̂_n(ξ)| envelope via R75 Plancherel + R76 conservation + R77 spectrum.

**Existing theorem fitting (1)-(4): Faure 2009 semiclassical spectral gap** (priority-1 secondary route).

Possibly also: Dolgopyat-style spectral gap on hyperbolic dynamics (operates on transfer operator of a hyperbolic flow); Pollicott-Sharp style anisotropic Banach space spectral theory for transfer operators. These all live in the operator-theoretic / dynamical-systems category, NOT the Igusa-polynomial-Mellin category.

---

## Deliverables

In C:/Collatz/:

- `IGUSA_PRE_REGISTRATION.md` — pre-reg (locked 2026-05-14 before Phase 1)
- `IGUSA_A_HYPOTHESES.md` — Igusa rationality theorem
- `IGUSA_B_HYPOTHESES.md` — Denef-Hoornaert Newton-polyhedron
- `IGUSA_C_HYPOTHESES.md` — Bories-Veys non-degenerated surface
- `IGUSA_D_HYPOTHESES.md` — Monodromy Conjecture
- `IGUSA_E_HYPOTHESES.md` — Motivic Igusa
- `IGUSA_F_HYPOTHESES.md` — Bories-Cluckers polynomial mapping
- `IGUSA_G_HYPOTHESES.md` — Watson / stationary phase on Igusa
- `IGUSA_H_HYPOTHESES.md` — Bernstein-Sato b-function
- `IGUSA_I_HYPOTHESES.md` — monomial / linear-form Igusa
- `IGUSA_J_HYPOTHESES.md` — open search for (1+p)^u-specific Igusa
- `IGUSA_DISPOSITION.md` (this file) — headline

PDF extractions (UTF-8 from pypdf): C:/tmp/igusa/*.txt (10 PDFs, all extracted cleanly).

No git operations performed (per discipline).

---

End disposition.
