# Collatz analytical-phase reading set

Generated 2026-05-06. Three reading clusters mapping to the Hank brief.

The empirical phase has characterized ρ_slow ≈ 0.83 as the genuine inverse-limit convergence rate of finite-k truncations to π_∞ on Z_3 (in L¹/TV/L²/L∞), but this rate doesn't appear as an eigenvalue of K_k or singular value of R_k on natural L² function spaces. Reading goal: find the operator + function space where ρ_slow IS an isolated dominant eigenvalue.

---

## CLUSTER A — Anisotropic Banach space methods for transfer operators

**Folder:** `C:/Collatz/fluid_dynamics/`

**Technique in one line:** construct a function space anisotropically tuned to the dynamics so the transfer operator has discrete spectrum instead of band-supported / continuous spectrum on natural L².

### Pedagogical entry point (start here)
- **`demers_2019_gentle_intro_anisotropic_banach_arxiv.pdf`** ★★★ — Demers, "A gentle introduction to anisotropic Banach spaces" (arXiv:1901.00131). Exactly the survey-level treatment of the technique you described. **Start here.**

### Foundational
- **`baladi_tsujii_2007_anisotropic_arxiv.pdf`** ★★★ — Baladi-Tsujii (Ann Inst Fourier 2007 — your high priority Cluster A pick). arXiv math/0505015.
- **`faure_sjostrand_2011_anosov_resonances_arxiv.pdf`** ★★★ — Faure-Sjöstrand (Comm Math Phys 2011 — your high priority). arXiv 1003.0513.

### Reviews / orientation (your second priority)
- **`baladi_2003_dynamical_zeta_arxiv.pdf`** — Baladi 2003 "Dynamical zeta functions" survey. Earliest free precursor to the 2018 book.
- **`baladi_2014_quest_anisotropic_banach_arxiv.pdf`** — Baladi 2014 "The quest for the ultimate anisotropic Banach space" (arXiv 1408.2937). Direct precursor to the 2018 book.
- **`baladi_2019_kneading_determinants_arxiv.pdf`** — Baladi 2019 (arXiv 1907.03453). Most-recent in this lineage.

### Demers-Liverani technique line
- **`demers_liverani_pene_2013_martingale_anisotropic_arxiv.pdf`** ★★ — "Martingale approximations and anisotropic Banach spaces with an application to the time-one map of a Lorentz gas" (arXiv 1301.0168). Concrete application of the 2007 framework.
- **`demers_liverani_pene_2021_stability_piecewise_hyperbolic_arxiv.pdf`** ★★ — "Stability of Statistical Properties in Two-Dimensional Piecewise Hyperbolic Maps" (arXiv 2104.06947). Same authors, more recent.
- **`demers_liverani_2008_stability_open_systems_arxiv.pdf`** ★ — Demers-Liverani "Stability of statistical properties in two-dimensional open systems" (arXiv 0710.2456). Predecessor.
- **`liverani_saussol_2006_clt_sequential_arxiv.pdf`** — "Central Limit Theorem for Sequential Dynamical Systems" (arXiv math/0602067). Technique reference for non-stationary settings.
- **`demers_2016_projective_cones_sequential_billiards_arxiv.pdf`** — "Projective Cones for Sequential Dispersing Billiards" (arXiv 1506.02836). Cone-method variant.

### Adjacent (already on disk, related topic)
- `liverani_2004_contact_anosov_annals.pdf` — Liverani 2004 "On contact Anosov flows" (Annals)
- `giulietti_liverani_pollicott_2013_anosov_zeta_annals.pdf` — Capstone application of the program

### Liverani lecture notes (substitutes for the missing 2004 survey)
- `liverani_notes_ravello.pdf`, `liverani_notes_barcellona.pdf`, `liverani_stonybrook_notes.pdf`, `liverani_lezioni_sigh.pdf` — All 2016 pedagogical notes covering similar material to the requested 2003-2005 survey.

### ❌ Not obtained
- **Liverani 1995 Annals "Decay of correlations"** — paywalled at JSTOR (https://www.jstor.org/stable/2118636); pre-2000 Annals doesn't host free PDFs. Annals stable: https://annals.math.princeton.edu/1995/142-2/p02
- **Liverani 2004 survey "Invariant measures and their properties: a functional analytic point of view"** — Birkhauser/CIME book chapter, not on arXiv or Liverani's homepage. The 2016 lecture notes substitute decently.

---

## CLUSTER B — Spectral analysis of random walks on profinite groups

**Folder:** `C:/Collatz/varju_followups/`

### Direct fits — abelian / profinite / multiplicative-mod-q (load-bearing)
- **`varju_2013_random_walks_compact_arxiv.pdf`** ★★★ — Varjú 2013 (arXiv:1209.1745). The starting point.
- **`eberhard_varju_2020_sharp_mixing_arxiv.pdf`** ★★★ — Eberhard-Varjú "Sharp mixing time of x_{n+1}=ax_n+b_n on Z/qZ" (arXiv 2003.08117). Mixing rate governed by entropy of an associated **Bernoulli convolution** — your "specific functional" question, abelian, multiplicative mod q. Bullseye for the framework.
- **`hermon_olesker_taylor_2024_cutoff_abelian_arxiv.pdf`** ★★★ — "Cutoff for Almost All Random Walks on Abelian Groups" (arXiv 2403.12355). Directly abelian, all groups including Z/p^k.
- **`hermon_olesker_taylor_2021_cutoff_abelian_arxiv.pdf`** ★★ — Earlier paper in the series (arXiv 2102.02809).
- **`hussain_lamzouri_2023_functional_clt_arxiv.pdf`** ★★★ — Functional CLT for normalized partial-character-sum paths on Z/pZ (arXiv 2304.13025). Skorokhod convergence in C[0,1]. Directly the "convergence of specific functionals" axis.
- **`pierce_weisbart_2024_padic_donsker_arxiv.pdf`** ★★★ — First weak/Skorokhod functional convergence of a discrete random walk on Z_p to p-adic Brownian motion (arXiv 2407.05561). All three boxes: functional + abelian + profinite.
- **`ayyer_singla_2019_random_walks_commrings_arxiv.pdf`** ★★ — Character-theoretic eigenvalue analysis of additive+multiplicative random walks on commutative rings, with Z/p^k Z explicitly worked out (arXiv 1605.05089).
- **`weisbart_2024_padic_brownian_arxiv.pdf`** ★★ — p-adic Brownian motion technical foundation (arXiv 2010.05492).

### Technique reference (non-abelian but methodological)
- **`bourgain_gamburd_2016_spectral_gap_sud_arxiv.pdf`** — Bourgain-Gamburd-style "Spectral Gap Theorem in SU(d)" (arXiv 1607.01530). The closest free Bourgain-Gamburd paper to the 2008 SU(2) original.
- **`bekka_guivarch_2011_spectral_gap_unitary_reps_arxiv.pdf`** — Bekka-Guivarc'h "A spectral gap property for random walks under unitary representations" (arXiv 1108.3146). Closest available substitute for the requested Bekka-de la Harpe Property (T) abelian-groups chapter.
- **`guivarch_2009_ergodicity_spectral_gap_arxiv.pdf`** — Guivarc'h "Ergodicity of group actions and spectral gap, applications to random walks and Markov shifts" (arXiv 0908.0637).
- **`guivarch_lepage_2017_spectral_gap_linear_walks_arxiv.pdf`** — Guivarc'h-Le Page "Spectral gap properties for linear random walks and Pareto's asymptotics for affine stochastic recursions" (arXiv 1705.09593).
- **`benoist_quint_2012_recurrence_linear_groups_arxiv.pdf`** — Benoist-Quint "Recurrence and ergodicity of random walks on linear groups and on homogeneous spaces" (arXiv 1204.6004).

### Cutoff lineage (foundational technique)
- **`nestoridi_olesker_taylor_2020_cutoff_nilpotent_arxiv.pdf`** — "Cutoff for random Cayley graphs of nilpotent groups" (arXiv 2008.08564).
- **`lubetzky_peres_2016_cutoff_ramanujan_arxiv.pdf`** — "Cutoff for Ramanujan graphs via degree inflation" (arXiv 1610.04357).
- **`nestoridi_2015_refined_random_walks_symmetric_arxiv.pdf`** — Nestoridi "Refined estimates for some basic random walks on the symmetric and alternating groups" (arXiv 1512.02361).
- **`breuillard_varju_2022_arxiv.pdf`** — Breuillard-Varjú joint (arXiv 1909.09053).

### Background technique
- `lamzouri_zaharescu_2011_character_sums_arxiv.pdf`, `hussain_2022_limiting_dist_character_sums_bristol.pdf` — random-walk-on-Z/mZ models for character sums.
- `chatterjee_diaconis_2020_arxiv.pdf`, `bate_connor_2018_arxiv.pdf`, `hildebrand_2008_arxiv.pdf`, `shkredov_2021_arxiv.pdf` — adjacent technique.

### ❌ Not obtained
- **Bourgain-Gamburd 2008 SU(2)** Inventiones — paywalled, IAS/Brooklyn homepages 404. The 2016 SU(d) paper is the closest free substitute.
- **Bekka-de la Harpe "Kazhdan's Property (T)"** book — Cambridge, paywalled. No arXived chapter.
- **Diaconis-Shahshahani 1981** "Random transpositions" Z. Wahrsch. — predates arXiv, paywalled at SpringerLink.
- **Saloff-Coste survey on tower-of-quotients** — most of his survey work is in Springer/Birkhauser book chapters, paywalled. No arXived survey.

---

## CLUSTER C — Inverse-limit operator construction on Z_p

**Folder:** `C:/Collatz/inverse_limit_operators/`

### Bendikov-Pittet line (your starting reference, 2011 already in references/)
- **`bendikov_pittet_2023_spectral_locally_finite_arxiv.pdf`** ★★★ — "Spectral properties of a class of random walks on locally finite groups" (arXiv 2307.04538). **2023 follow-up** to the 2011 paper you have. This is the explicit non-locally-finite extension target.
- **`bendikov_pittet_2013_l2_isoperimetric_arxiv.pdf`** ★★ — Bendikov-Pittet "Spectral distribution and L²-isoperimetric profile of Laplace operators on groups" (arXiv 1304.6271).
- **`bendikov_pittet_2009_braiding_schur_arxiv.pdf`** — Bendikov-Pittet "Braiding and asymptotic Schur's orthogonality" (arXiv 0901.0271). Earlier in the series.

### Cross-reference (already in `varju_followups/`)
The Pierce-Weisbart 2024 + Weisbart 2024 papers are the most directly relevant Cluster C content — they construct operators on Z_p inheriting structure from finite quotients.

### ❌ Not obtained / didn't find
- **Recent (2020+) Ruelle-Pollicott resonances for finite-state Markov chains, profinite limits** — arXiv search returned no direct hits. The corner is sparse / underdeveloped per the prior Varjú-followup hunt: "no one has assembled all three legs (functional + abelian + profinite)."
- **Faure or Sjöstrand specifically on non-Anosov / discrete settings** — their work is overwhelmingly in the Anosov flow / hyperbolic dynamics literature. No discrete/Markov-chain pivot found in 2020+ arXiv.

---

## Summary

| Cluster | Folder | Papers on disk | Strongest hits |
|---|---|---:|---|
| A — anisotropic Banach | `fluid_dynamics/` | 16 | Demers 2019 gentle intro, Baladi-Tsujii 2007, Faure-Sjöstrand 2011 |
| B — profinite random walks | `varju_followups/` | 23 | Eberhard-Varjú 2020, Hussain-Lamzouri 2023, Pierce-Weisbart 2024, Hermon-Olesker-Taylor 2024 |
| C — inverse-limit operators | `inverse_limit_operators/` | 3 | Bendikov-Pittet 2023 |
| **Total** | | **42 PDFs** | |

## Reading order recommendation

1. **`demers_2019_gentle_intro_anisotropic_banach_arxiv.pdf`** (Cluster A, pedagogical) — get the technique
2. **`baladi_2014_quest_anisotropic_banach_arxiv.pdf`** (Cluster A, meta-survey) — understand "what function space?" as a research question, the closest meta-match for your (Z/3^k)* question
3. **`eberhard_varju_2020_sharp_mixing_arxiv.pdf`** (Cluster B, bullseye) — they've done abelian + multiplicative-mod-q + functional convergence, just for the additive-multiplicative affine recursion not the qx+1 map
4. **`hussain_lamzouri_2023_functional_clt_arxiv.pdf`** (Cluster B, functional CLT) — directly functional convergence on Z/pZ
5. **`pierce_weisbart_2024_padic_donsker_arxiv.pdf`** (Cluster B/C bridge) — first profinite + functional + abelian convergence theorem
6. **`bendikov_pittet_2023_spectral_locally_finite_arxiv.pdf`** (Cluster C) — operators inheriting structure from finite quotients
7. **`baladi_tsujii_2007_anisotropic_arxiv.pdf`** + **`faure_sjostrand_2011_anosov_resonances_arxiv.pdf`** (Cluster A foundations) — once you have technique, dig into the original constructions

## Honest gap report

The user's exact question — *functional convergence of inverse-limit transfer operators on (Z/p^k)* with isolated rate ρ_slow* — has no direct hit in the literature. Three building-block lineages each cover ~⅔ of the combination but no one has fused all three. The framework would be filling this corner.
