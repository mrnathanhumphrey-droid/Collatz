# FAURE_PRE_REGISTRATION — Faure-school semiclassical spectral gap on Syracuse / Tao recursion transfer operator

**Date:** 2026-05-14. Locked BEFORE Phase 0 reading.
**Mode:** E (verbatim from PDFs in C:/Users/Nate/OneDrive/Documents/faure_semiclassical/pdfs/).
**Probe context:** TRIPLE convergent secondary routing from BGT, ADELIC, IGUSA dispositions. IGUSA's NO_FIT elevated Faure 2009 to priority-1 (the only remaining categorically-distinct route in the prior probe map). The closure target requires POSITIVE IRRATIONAL spectral radii (negative-and-rational Igusa poles excluded), a category natively supported by transfer-operator spectral theory.

## Candidate list

- **A. Faure 2009** — Semiclassical origin of spectral gap for partially expanding maps (FOUNDATIONAL)
- **B. Faure-Sjostrand 2010** — Upper bound on density of Ruelle resonances for Anosov flows
- **C. Faure-Tsujii 2013** — Band structure of Ruelle spectrum of contact Anosov flows (CRAS)
- **D. Faure-Tsujii 2021** — Micro-local analysis + band structure refinement
- **E. Faure 2008** — Semiclassical Anosov diffeomorphisms and Ruelle resonances
- **F. Dyatlov-Guillarmou 2014** — Pollicott-Ruelle resonances for open systems
- **G. Datchev-Dyatlov-Zworski** — Sharp polynomial bounds on number of Pollicott-Ruelle resonances
- **H. Faure 2025** — "Can quantum dynamics emerge from classical chaos?" (overview)
- **I. Liverani 1995 / Gouëzel-Liverani 2006** — anisotropic Banach space framework (cross-ref, supporting)
- **J. Pollicott-Sharp 2018 / Parry-Pollicott 1990** — dynamical zeta function framework (cross-ref)
- **K. Any candidate surfaced during Phase 0 reading**

## Pre-registered priors

| Candidate | Prior | Rationale |
|---|---|---|
| A. FAURE_2009 | HIGHEST (0.40) | Foundational, partially-expanding maps, positive irrational spectral radii natively supported, PADE multi-spectral picture matches resonance band structure prediction |
| B. FAURE_SJOSTRAND_2010 | moderate-high (0.20) | Density of resonances, useful for converting band structure to polynomial-in-A bound |
| C. FAURE_TSUJII_2013 | moderate (0.15) | Band-structure prediction matches PADE multi-spectral; could SELECT if hypotheses match |
| D. FAURE_TSUJII_2021 | moderate (0.10) | Micro-local refinement; relevant if Syracuse needs micro-local treatment |
| E. FAURE_2008 | moderate (0.10) | Diffeomorphism case; relevant if Tao recursion is diffeomorphism rather than partially expanding |
| F. DYATLOV_GUILLARMOU | moderate (0.10) | Open systems; check if Syracuse has natural escape behavior |
| G. DATCHEV_DYATLOV_ZWORSKI | moderate-high (0.15) | Sharp polynomial bound is exactly the conclusion shape we want |
| H. FAURE_2025 | uncertain (0.05) | Overview, check for direct applicability |
| I. LIVERANI_GOUEZEL | supporting (0.05) | Anisotropic Banach framework; not standalone |
| J. POLLICOTT_SHARP / PARRY_POLLICOTT | low (0.05) | Zeta-function alternative; less direct match |

## Decision rules

- **SELECTED**: a candidate's hypotheses are SATISFIED by Syracuse / Tao recursion as a partially-expanding map (or close enough) + theorem produces a spectral-gap statement that converts via R75/R76/R77 to polynomial-in-A bound, AND accommodates the multi-spectral PADE picture (transient + asymptotic + complex pair at period ~9.2)
- **PARTIAL**: theorem applies in spirit but Syracuse doesn't literally satisfy the partial-expansion hypothesis; specify what technical work (e.g., adapting Faure 2009 proof for a profinite-space partially-expanding map) would close the gap
- **NO_FIT**: at least one hypothesis FAILED for every candidate; Faure-school category cannot accommodate Syracuse's specific class (TENTH category-of-object barrier)
- **BLOCKER**: theorem statement UNVERIFIABLE; need more literature
- **MODE_H_CIRCULAR**: theorem requires |μ̂_n(ξ)| or stationary-measure decay as INPUT (rare for Faure-school, which derives spectrum from map data)

## Locked thresholds

- Hypothesis SATISFIED requires verbatim match between Syracuse object and theorem's required object class
- Partial-expansion verification: must specify expansion direction(s), contraction direction(s), and uniform-vs-measure-theoretic nature of expansion
- Smoothness gap (Faure works in C^∞ setting; Syracuse is profinite Markov chain) is a LOAD-BEARING discriminator. If smoothness is required at the proof technique level (not just statement), this is NO_FIT or PARTIAL
- PADE picture match: complex pair period ≈9.2 with θ ≈ 0.68 rad, asymptotic at z ≈ 1.016, leading at z ≈ 1.57 (transient at n=13). Theorem's spectral-gap location must predict at least one of these
- Conversion check: theorem's conclusion must be convertible via R75 Plancherel + R76 conservation + R77 spectrum to |μ̂_n(ξ)| envelope. If conversion is non-mechanical, PARTIAL

## Pre-registered SECONDARY ROUTING if NO_FIT or PARTIAL

1. **Watson lemma / saddle-point on R78/R79 bilinear off-diagonal sum** — recommended by IGUSA as priority-2 alternative if Faure fails. Doesn't directly close c=7/45 but tightens structural picture at chain side
2. **Multi-singularity Flajolet-Sedgewick VI.4-VI.5** — recommended by PADE_NUMERICAL_DISPOSITION. The relevant FS section once z=2 is no longer the dominant singularity and complex-conjugate pair structure is visible
3. **Direct construction of an adapted spectral theorem for Syracuse's specific profinite setting** — research-grade work, last resort. Wilson's PADE picture (complex pair period 9.2 + asymptotic z=1.016) provides the roadmap for what such a theorem must produce

## Discipline notes

- Mode E: hypotheses VERBATIM from PDF, no inheritance
- pypdf only for extraction; UTF-8 file write before read
- Pre-registration locked before Phase 0
- Don't git commit (Nathan commits manually)
- Compute is FREE

## Expected outcomes (honest priors)

- SELECTED (~25-35%): Faure 2009 or Faure-Tsujii band structure directly applies
- PARTIAL (~30%): Theorem fires for idealized smooth Syracuse, profinite extension needed
- NO_FIT (~20%): Smoothness gap load-bearing (tenth category-of-object barrier)
- BLOCKER (~10%): Specific paper missing from corpus
- MODE_H_CIRCULAR (~5%): unlikely

End pre-registration. Locked at this state.
