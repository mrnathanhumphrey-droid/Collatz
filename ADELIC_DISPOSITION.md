# ADELIC_DISPOSITION — adelic Mellin / Tate-style probe headline

**Date:** 2026-05-14. Probe ADELIC. Working dir: C:/Collatz/. Mode E (verbatim from 11 PDFs at C:/Users/Nate/OneDrive/Documents/adelic_mellin/pdfs, extracted via pypdf to C:/tmp/adelic/).

Pre-registration: `C:/Collatz/ADELIC_PRE_REGISTRATION.md`.

---

## Headline

**NO_SELECTED. Net disposition: 10 NO_FIT + 0 PARTIAL + 1 STRUCTURAL_NEAR_MISS_KEPT_AS_NO_FIT (candidate E).**

The Tate adelic-Mellin framework and all closely related modern adelic / affine-sieve / heat-kernel-on-buildings theorems FAIL the Syracuse μ_n hypothesis class categorically. The failure modes cluster into three categorical barriers:

1. **Additive vs multiplicative.** Tao's μ̂_n(ξ) is the *additive* Fourier coefficient on Z/3^n; Tate's local zeta is a *multiplicative* Mellin against an *idele class character*. The two transforms live on dual groups in different categories (Z/3^n's Pontryagin dual is itself Z/3^n; Tate operates on the multiplicative idele group). Candidates A, C, D fail at this barrier.

2. **Exceptional (parabolic) vs non-exceptional (CKW non-exceptional / BGS Zariski-dense).** Tao's family {M_v = [[3, 1], [0, 2^v]]} sits inside the Borel B ⊂ PGL_2(ℚ_3) = stabilizer of ω ∈ ∂T_3 — which is the *exceptional / unimodular* case in CKW's terminology, the case CKW's main theorems explicitly *exclude*. BGS / Sarnak / Kontorovich require non-trivial limit set (Hausdorff dimension > 1/2 or Zariski density); Tao's limit set is a single boundary point. Candidates E, H, I fail at this barrier.

3. **Symmetric / K-bi-invariant / reversible vs forward-only.** ASTrojan heat-kernel requires K-bi-invariant (spherical) isotropic random walks; Saloff-Coste requires symmetric reversible μ. Tao's walk is neither (asymmetric forward-only, parabolic-fixing). Candidates F, J fail at this barrier.

The combined failure is the **eighth category-of-object barrier** in the systematic obstruction map for c = 7/45 closure (following 5-probe Fourier-decay, Cluster 1 Cochrane, Cluster 2 BMP, Bruhat-Tits/BKL, Tauberian, Furstenberg-Guivarc'h, BGT).

---

## Pre-registered probability vs. realized outcome

| Outcome | Pre-Phase-0 | Realized |
|---|---|---|
| SELECTED | 12% | NOT |
| PARTIAL | 30% | NOT (E was the closest, but FAILED at h2 non-exceptional) |
| NO_FIT | 30% | **REALIZED** |
| MODE_H_CIRCULAR | 15% | partial (candidates A, C, G show Mode H fingerprints but the primary failure is hypothesis mismatch, not pure circularity) |
| BLOCKER | 13% | NOT |

The prior of 30% NO_FIT was slightly conservative. Realized NO_FIT is the dominant outcome.

---

## Summary table

| Code | Theorem | Phase 0 | Phase 1 hyp check | Phase 2 conclusion | Phase 3 factorization | Disposition |
|---|---|---|---|---|---|---|
| A | Tate adelic Mellin FE (Binder Thm 5.18) | EXTRACTED | FAILED h2 (no L¹(𝔸) function), h3, h4 (additive vs mult), h5 | ADELIC FE delivered IF hyps held; doesn't here | ADELIC_FACTORIZATION_INHERENT | **NO_FIT (Mode H fingerprint)** |
| B | Tate archimedean local factor (Binder p. 13-14) | EXTRACTED | FAILED h1 (Syracuse has no archimedean component), h2, h3 | Archimedean Γ-factor; not applicable | ARCHIMEDEAN_VISIBLE (vacuously) | **NO_FIT** |
| C | Tate p-adic local factor (Binder §3.7) | EXTRACTED | FAILED h3 (additive vs mult), h4 | Local L; resonance with R77 spectrum but no formal identification | NON_ARCH_ONLY | **NO_FIT (structural resonance kept)** |
| D | Adelic Poisson summation (Binder Thm 5.14) | EXTRACTED | FAILED h3 (μ_n not L¹(𝔸)) | Adelic Poisson; not applicable. R75 finite-group Plancherel is the correct analog | ADELIC_FACTORIZATION_INHERENT | **NO_FIT** |
| E | Cartwright-Kaimanovich-Woess (Thms 2, 9, §4) | EXTRACTED | **PARTIAL FIRE on h1 (Tao's family DOES live in AFF(ℚ_3) Borel)** but FAILED h2 (Tao's group is exceptional, fixes ω = ∞) and h5 (not spread-out) | If h2/h5 held, would give Poisson boundary ≡ harmonic measure on ∂T_3 — *single-place* closure only | NON_ARCH_ONLY | **NO_FIT (structural near-miss)** |
| F | Anker-Schapira-Trojan 2013 heat kernel (Thm 3) | EXTRACTED | FAILED h2 (Tao not K-bi-invariant / isotropic), h3 | Sharp asymptotic p(n; v_n); not applicable | NON_ARCH_ONLY | **NO_FIT (categorical)** |
| G | Chambert-Loir / Tschinkel 2009 Igusa (Thms 1.2.1, 1.3.1) | EXTRACTED | FAILED h1 (no variety), h2 (no height), h3 | Volume asymptotic via meromorphic continuation of Z(s); not applicable. Closer relative = plain Igusa local zeta | ADELIC_FACTORIZATION_INHERENT | **NO_FIT (Mode H fingerprint)** |
| H | Bourgain-Gamburd-Sarnak (Thms 1.1, 1.2, 1.3) | EXTRACTED | FAILED h1 (Γ in PGL_2(ℚ_3) not SL_2(ℤ)), h2 (limit set = single point) | Spectral gap on congruence covers; not applicable | GLOBAL_BUT_PLACE_BLIND | **NO_FIT (categorical)** |
| I | Kontorovich 2014 levels (BV-style equidist) | EXTRACTED | FAILED h1 (Tao not group action), h2 | Almost-primes in orbit; not applicable | GLOBAL_BUT_PLACE_BLIND | **NO_FIT (categorical)** |
| J | Saloff-Coste 2001 (Thms 2, 3, 4) | EXTRACTED | FAILED h2/h3 (Tao not symmetric / reversible) | φ(n) decay via volume growth; not applicable | GLOBAL_BUT_PLACE_BLIND | **NO_FIT (categorical, FG-confirmed)** |

Total extractable: 10/10 (no BLOCKER).
- 0 SELECTED
- 0 PARTIAL
- 10 NO_FIT
- 2 with Mode H fingerprint embedded (A, G — but primary failure is hyp mismatch, not pure circularity)

No surfaced candidate K (Sarnak survey extends BGS, same disposition; Affine_Sieve_Beyond_Expansion extends Kontorovich, same disposition; ASTrojan 2007 is earlier ASTrojan 2013, same disposition).

---

## Final disposition: **NO_FIT (dominant)**

The adelic / Tate-style category does not contain a theorem that fires on Syracuse μ_n as currently formalized. This **closes the adelic-Mellin substrate as a closure route for c = 7/45**.

Of particular note are the **structural near-misses**:

### Near-miss 1 (candidate E): Tao's family DOES live in AFF(ℚ_3)

The Cartwright-Kaimanovich-Woess setup for random walks on Aff(local field) matches Tao's setting *up to one critical hypothesis*. The Tao-generated subgroup of AFF(ℚ_3) is precisely the *exceptional case* in CKW's classification — the case where their main theorems (2-9) explicitly *fail to apply* because Γ fixes a boundary point.

This is structurally identical to BT_DISPOSITION's Negative-case Q2 finding: "Tao's algebraic content is 'a discrete random walk in the Borel B ⊂ PGL_2(ℚ_3)'." The Borel B = AFF(ℚ_3) is the parabolic / exceptional subgroup. CKW's theorems do not apply *because* Γ is in the parabolic case — and this is structural, not coincidental.

**Conclusion: the categorical reason CKW fails is the same reason BT fails: Tao's group fixes ω ∈ ∂T_3. The exceptional / parabolic case is exactly the categorical barrier.**

### Near-miss 2 (candidate C): R77 T_diag eigenvalues {0, 1} mirror Tate local poles at s ∈ {0, 1}

R77 T_diag = (1/5) [[1, 1], [4, 4]] has characteristic polynomial λ² − λ, hence eigenvalues 0 and 1. Tate's unramified ℚ_3 local L-factor 1/(1 − 3^{-s}) has simple poles at s = 0 and (via FE) s = 1. The two pieces have the *same pole structure* (in different variables), suggestive of an underlying Tate-style identification.

But no formal identification was established — they live in different operator categories (R77 T_diag acts on (P_+, P_−) class-resolved deviation space; Tate's local L acts on Schwartz-Bruhat functions on ℚ_3*). The resonance is *suggestive but not formal*. Treating it as a closure path would require building the bridge from scratch, which is Mode H circular (the bridge would need to derive F_3 from R77 data, and "F_3 from chain-side data" is the closure target).

### Near-miss 3 (candidate G/A): the meromorphic-continuation-past-first-pole template

Both CLT (Igusa) and Tate have theorems whose CONCLUSION is "Z(s) admits meromorphic continuation past the first pole, encoding the rate". This is *exactly* the closure-target shape. The hypotheses for both theorems involve geometric/adelic ingredients Syracuse doesn't supply.

If a future probe COULD supply those ingredients (e.g., by lifting μ_n to a measure on the principal-unit coset of ℚ_3* and applying Igusa local zeta to the (1+3)^u polynomial), then the meromorphic-continuation conclusion would directly close c = 7/45. This is the BGT-flagged Igusa local zeta route — categorically distinct from CLT's adelic-Mellin form, but using the same Mellin-of-polynomial-on-local-field machinery.

---

## What category of theorem is missing

The closure target requires a theorem that:

1. Operates on a *finite-group measure* (Syracuse μ_n on (Z/3^n)*), or its inverse limit (measure on ℤ_3*).
2. Produces a *Mellin-style transform* in the multiplicative direction (since R77.6's branch-cut at z = 2 corresponds, via z = 1/(1/2)= q^s with q = 3, s = log_3(2), to a non-Tate Mellin location — i.e., needs a non-standard Mellin variable).
3. Handles the *non-symmetric forward-only* dynamics of Tao without requiring reversibility / K-bi-invariance / spread-out / non-exceptional.
4. Produces an explicit *polynomial-in-A bound* on |μ̂_n(ξ)| or a closed-form for F_∞(s) or F_3(s).

**No theorem in the 11-PDF adelic-Mellin corpus matches all four criteria.** The adelic-Mellin family is built on (1') Schwartz-Bruhat functions on full adele rings 𝔸; (2') multiplicative idele class characters; (3') symmetric / K-bi-invariant / non-exceptional / spread-out hypotheses; (4') analytic continuation conclusions whose hypotheses include the closure target.

**Missing category:** something like a "stochastic-iteration Mellin transform on a profinite multiplicative group, with controlled deviation envelope via additive Fourier coefficients" — a hybrid object that the adelic-Mellin tradition doesn't address. The closest existing relatives are:

- **Igusa local zeta** (operates on a single local field, single polynomial f; no random walk; could potentially fire on R78 (1+3)^u). NOT in adelic-Mellin corpus.
- **Faure 2009 semiclassical spectral gap** (operates on partially expanding maps; could potentially fire on Tao recursion's transfer operator). NOT in adelic-Mellin corpus.
- **Bingham-Ostaszewski sequential RV** (operates on the ε_k sequence; PARTIAL fired in BGT probe). Outside adelic-Mellin.

---

## SECONDARY ROUTING

Per pre-registration locked priorities:

1. **Igusa local zeta** (top priority). Operates on R78 (1+3)^u polynomial directly. Categorically distinct from adelic Mellin (it's a *single-variable* meromorphic-continuation statement on a *single* local field, no adele product). Would explicitly compute Z_3(s) = ∫_{ℚ_3} |f(x)|^s dx for f related to (1+3)^u − c, and identify the singularity controlling the rate. **PRIORITY: HIGH** — same priority assignment as BGT_DISPOSITION.

   Why it might work where the present probe failed:
   - Bypasses the additive-vs-multiplicative barrier (Igusa Z(s) is a multiplicative Mellin of |f|^s, while Syracuse's additive Fourier μ̂_n(ξ) can be related to Igusa Z via Mellin-of-Fourier-coefficient identities — a known dual);
   - Bypasses the exceptional-parabolic barrier (Igusa doesn't care about group structure of the integration variable, just the polynomial);
   - Bypasses the symmetric-reversible barrier (Igusa is a static integral, no random walk);
   - Delivers explicit meromorphic-continuation poles with computable residues.

2. **Faure 2009 semiclassical spectral gap** (next priority). Operates on the transfer operator of partially expanding maps. Tao recursion's transfer operator is exactly this kind of object. Would directly address R77 off-diagonal rate-½.

3. **Heat-kernel narrowing** (deprioritized): F was NO_FIT not PARTIAL, so no narrowing-from-PARTIAL is possible.

4. **Watson lemma / saddle-point on R78/R79 bilinear** (parked).

### Recommended top-priority secondary route: **Igusa local zeta**.

This is the same recommendation as BGT_DISPOSITION made on 2026-05-13. The present ADELIC probe *confirms and strengthens* the Igusa recommendation:
- Adelic-Mellin corpus is structurally close to Igusa (CLT generalizes Igusa to varieties), but CLT requires geometric ingredients Syracuse doesn't supply. Pure Igusa local zeta on R78's (1+3)^u polynomial is the cleanest match.
- The Igusa route is **NON_ARCH_ONLY** (delivers F_3 only, not F_∞). This is consistent with the BT_DISPOSITION finding that "the attractor lives at the archimedean place" — Igusa would give a partial closure (F_3 piece) but leave the F_∞ piece open. That is acceptable as a *next* probe outcome: a PARTIAL on F_3 would tighten the structural picture even if it doesn't close F_∞.

---

## Surprises in the inputs

### Surprise 1: Tao's family is precisely the exceptional case in CKW

The CKW probe surfaced that Tao's family {M_v} not only lives in AFF(ℚ_3) (already noted by BT_DISPOSITION as the Borel B picture), but lives there in CKW's *EXCEPTIONAL* / *unimodular* case — which is exactly the case CKW's main theorems exclude. This is a sharper structural fact than BT's "lives in the Borel": CKW *names* the exceptional case and proves all its main theorems require non-exceptional.

So the obstruction is not just "Tao fixes ∞" (BT) but "Tao is in the CKW unimodular category, where the random-walk geometry is fundamentally different (HOR(T) acts on R simply transitively, harmonic measure is trivial or 1-dimensional)."

This unifies BT and ADELIC obstructions into a single structural statement: **Tao's algebraic content sits in the CKW-exceptional / parabolic-stabilizer case at every local-field place, including ℚ_3.**

### Surprise 2: R77 T_diag spectrum {0, 1} = Tate ℚ_3 pole structure {s = 0, s = 1}

Both objects produce the pair (0, 1) as their analytic data: R77 T_diag eigenvalues, Tate ℚ_3 local L poles. The coincidence is structurally interesting and might be more than coincidence: the (1, 4) deviation eigenvector of T_diag corresponds to "preserving total Plancherel mass S = P_+ + P_−"; the s = 1 pole of Tate's ζ_3 corresponds to the "trivial / Tamagawa" character residue. Both are "mass-preserving" in their respective senses.

However, no formal identification can be made within the adelic-Mellin corpus. This would require a new theorem bridging deviation operators on class-resolved Plancherel spaces to Tate-local pole structure — outside the corpus.

### Surprise 3: R77.6 branch-cut at z = 2 corresponds to s = log_3(2) in Tate-Mellin variable

If we identify the generating-function variable z with the multiplicative-Mellin variable via z = q^s, q = 3, then R77.6's branch-cut at z = 2 sits at s = log_3(2) ≈ 0.631. This is *not* a standard Tate pole location (which would be s = 0 or s = 1). So R77.6's branch-cut is NOT a Tate-local-factor pole.

It could be a *non-Tate* pole — e.g., an Igusa local zeta pole for a polynomial f with v_3-valuation properties that produce s = log_3(2) as a pole. This would be the kind of identification Igusa local zeta could probe directly. So Igusa is the right secondary route specifically for the R77.6 branch-cut diagnosis.

### Surprise 4: BT_DISPOSITION's "archimedean attractor" finding restated

BT said "the 1-attractor is archimedean" — viz., r_n stops being > 1 in the archimedean norm. The ADELIC probe sharpens this: the trajectory r_n's archimedean behavior is **measure-level invisible** because μ_n lives entirely on the non-archimedean side (profinite (Z/3^n)*). The archimedean information is in the *deterministic* trajectory r_n, not in the *measure* μ_n.

This means the closure of c = 7/45 (which is a measure-level statement: Plancherel mass) cannot directly invoke the archimedean attractor. The archimedean attractor is a *consequence* of the measure-level closure, not its proof input.

This reframes the BT-flagged adelic substrate requirement: not "the closure proof must use archimedean information" but "the closure proof must explain why the non-archimedean measure-level behavior (Plancherel mass on (Z/3^n)*) reflects the archimedean attractor (r_n → 1 in ℝ)." The two are linked structurally but not symmetrically.

**This is a genuinely new finding from the ADELIC probe** — sharpens BT_DISPOSITION's strategic note into a precise question.

---

## Strategic position

Pre-ADELIC: c = 7/45 closure had 7 categorical barriers mapped, 1 PARTIAL pending (BGT sequential RV in plateau k=2..6).

Post-ADELIC: c = 7/45 closure has 8 categorical barriers mapped. The adelic-Mellin substrate is closed as a direct closure route. The structural picture sharpens:

- **Tao's natural algebraic home** is AFF(ℚ_3) Borel = parabolic stabilizer of ω ∈ ∂T_3 = CKW-exceptional. This is *categorically incompatible* with the modern random-walks-on-affine-groups / Bruhat-Tits-tree / heat-kernel-on-buildings / affine-sieve / Tate-adelic family of theorems, all of which require *non-exceptional* / *Zariski-dense* / *K-bi-invariant* / *symmetric* hypotheses.

- **The archimedean attractor lives at the trajectory level**, not the measure level. The closure (which is a measure-level Plancherel-mass statement) doesn't directly need archimedean information — it needs a *non-archimedean / 3-adic local zeta-style analytic continuation* that explains why S_n = Σ|μ̂_n(ξ)|² → 7/15 at rate (1/2)^n.

- **The Igusa local zeta route remains open**, recommended as top secondary priority. It's the only one of the four pre-registered secondary routes that (a) avoids the categorical barriers identified in the present probe, (b) directly operates on the R78 (1+3)^u algebraic structure, (c) delivers explicit meromorphic-continuation conclusions with computable singularities, and (d) is fundamentally a non-archimedean single-place statement (consistent with the post-ADELIC understanding that the closure is non-archimedean).

---

## Deliverables

In C:/Collatz/:

- `ADELIC_PRE_REGISTRATION.md` — pre-reg (locked before Phase 0)
- `ADELIC_A_HYPOTHESES.md` — Tate adelic Mellin FE
- `ADELIC_B_HYPOTHESES.md` — Tate archimedean local factor
- `ADELIC_C_HYPOTHESES.md` — Tate non-archimedean local factor
- `ADELIC_D_HYPOTHESES.md` — Adelic Poisson summation
- `ADELIC_E_HYPOTHESES.md` — Cartwright-Kaimanovich-Woess 1994
- `ADELIC_F_HYPOTHESES.md` — Anker-Schapira-Trojan 2013 heat kernel
- `ADELIC_G_HYPOTHESES.md` — Chambert-Loir / Tschinkel 2009 Igusa
- `ADELIC_H_HYPOTHESES.md` — Bourgain-Gamburd-Sarnak
- `ADELIC_I_HYPOTHESES.md` — Kontorovich 2014 levels
- `ADELIC_J_HYPOTHESES.md` — Saloff-Coste 2001
- `ADELIC_DISPOSITION.md` (this file) — headline

PDF extractions (UTF-8 from pypdf): C:/tmp/adelic/*.txt (11 PDFs).

No git operations performed (per discipline).

---

End disposition.
