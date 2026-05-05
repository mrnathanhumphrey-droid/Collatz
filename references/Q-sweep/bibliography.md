## Critical-q phenomenon in qx+1 — literature characterization (2026-05-04)

### Context
Q-Sweep Test 1 (May 4, 2026): for odd primes q ∈ {3, 5, 7, 11, 13}, the partial Plancherel sums S_n^(q) of the qx+1 Syracuse measure satisfy S_n^(q) / S_{n−1}^(q) → q/3 (rel. error 9.6e-5 at q=13). This makes q=3 the unique critical odd prime where the ratio equals 1, S_n converges only at q=3 (constant c = 7/45), and S_n diverges geometrically for q ≥ 5 with explicit common-ratio q/3. The literature search below evaluates whether the *partial-sum-ratio law* and the *q=3 phase-boundary* framing have been stated previously.

### Findings

**(FRAMEWORK) Matthews (2010), "Generalized 3x+1 Mappings: Markov Chains and Ergodic Theory"** — http://www.numbertheory.org/PDFS/matthews-final-revised.pdf
The standard reference for Markov-chain analysis of generalized qx+1 maps. Provides the apparatus that *predicts* the convergent/divergent dichotomy through expected log-growth per step (≈ ½ log(q/2) under uniform-binary-valuation heuristic), but does not state the q/3 ratio law on Plancherel partial sums and does not name q=3 as a critical/marginal point. Frames the divide via density-1 statements about iterates of T_q, not via Fourier analysis. Closest framework match for the user's finding but does not articulate the partial-sum ratio.

**(FRAMEWORK / NEAR-DIRECT) Kontorovich & Lagarias (2009), "Stochastic Models for the 3x+1 and 5x+1 Problems"** — arXiv:0910.1944, https://arxiv.org/pdf/0910.1944
Provides parallel-but-distinct stochastic models for q=3 and q=5 (sections 2–6 vs 7–8) under section heading "5x+1 Function: Stochastic Models and Results." Establishes the rigorous heuristic that q=5 trajectories diverge for almost all starts. Does **not** contain the words "critical", "phase", "boundary", or "marginal" in the body, does not state the q/3 partial-sum ratio, and does not parametrize q continuously (treats q=3 and q=5 as separate case studies). This is the closest pre-existing comparison but the q-as-parameter / Plancherel-S_n framing is absent.

**(FRAMEWORK / DIRECTLY ADJACENT) Siegel, M. C. (2022 PhD; arXiv:1909.09733; Springer p-Adic Numbers Part I 2024 = doi:10.1134/S2070046624020055; Part II 2025 = doi:10.1134/S2070046625020062; arXiv:2208.11082; arXiv:2412.02902; arXiv:2601.17030)** — "(p,q)-Adic Analysis and the Collatz Conjecture", "Conservation of Singularities… Dreamcatchers for Hydra Maps", "Hydra Map and Numen Formalisms"
The single most-on-target framework in the literature. Defines, for each odd prime q, a "Numen" χ_q : Z_2 → Z_q and an associated Fourier-Stieltjes transform; develops a (p,q)-adic Plancherel/Wiener-Tauberian theory; and proves a Correspondence Principle linking q-divergent dynamics to non-vanishing of |χ̂_q|, q=3 (Collatz) to χ̂_3 having Fourier-Stieltjes mass arbitrarily q-adically close to zero infinitely often. **However:** Siegel works in *q-adic* absolute value (non-archimedean), not the archimedean partial-sum S_n on the unit interval; he does not state the q/3 ratio law, does not use the words "critical/phase/marginal" comparing q's, and the framework analyzes each q separately rather than scanning q. The user's archimedean-Plancherel q-sweep is a parallel observable, not the same object as Siegel's q-adic Fourier-Stieltjes magnitude. Strong structural cousin; not a duplicate.

**(NEAR-FRAMEWORK) Tao (2019), "Almost all Collatz orbits attain almost bounded values"** — arXiv:1909.03562, https://arxiv.org/pdf/1909.03562; companion blog https://terrytao.wordpress.com/2020/01/25/equidistribution-of-syracuse-random-variables-and-density-of-collatz-preimages/
Uses Plancherel's theorem on Z/3^n Z to control characteristic functions of "Syracuse random variables" and prove superpolynomial Fourier-coefficient decay. **q=3 only** — neither paper nor the equidistribution blog post extends the Plancherel argument to q ≥ 5 nor explains why the framework requires q=3. So: same archimedean tool (Plancherel) the user is invoking, but never run as a q-sweep.

**(FRAMEWORK / SECONDARY) Aicardi, "A simple Markov chain for the extended Collatz problem"** — IJRDO Comp Sci Eng, https://www.ijrdo.org/index.php/cse/article/download/804/752/ ; and Santos (2021), "On the Collatz general problem qn+1" — arXiv:2005.00346
Both reproduce the q=3-convergent / q≥5-divergent dichotomy (q=3 = "critical threshold" / Crandall conjecture) by elementary 3-state Markov-chain or probabilistic argument. **Neither states the q/3 ratio law on Plancherel sums.** They confirm the qualitative dichotomy is folkloric but the user's quantitative scaling law is not in either source.

**(NEGATIVE) Sinai (2003), "Statistical (3x+1) Problem"** — arXiv:math/0201102. q=3 only. Does not generalize to qx+1. No Plancherel-S_n analysis.

**(NEGATIVE) Lagarias annotated bibliography II 2000-2009** — arXiv:math/0608208. Catalog of ~150 entries; the "phase transition" / "critical-q" framing for the qx+1 family does not surface as a topic header. Crandall-style q=3-special results appear, but always as "q=3 is conjecturally the only convergent qx+1," never with a quantitative ratio law on Fourier/Plancherel partial sums.

**(PARALLEL) Erdős (1939); Solomyak (1995), "On the random series Σ ±λ^n (an Erdős problem)" Annals 142, 611–625; Peres-Schlag-Solomyak, "Sixty Years of Bernoulli Convolutions" — https://gauss.math.yale.edu/~ws442/papers/sixty.pdf**
The cleanest structural parallel. Bernoulli convolution ν_λ on R is **singular for λ < 1/2, uniform on [-1,1] at λ = 1/2 (the critical value), and absolutely continuous for almost every λ ∈ (1/2, 1)** with a known countable singular exceptional set (Pisot λ). The phase transition is at a single critical parameter λ=1/2; Plancherel/L² methods (Garsia, Solomyak, Peres-Solomyak, Hochman) detect it via Fourier-coefficient decay of the convolution. **This is the structural template the user's q=3 finding fits into:** a single critical parameter where a measure-theoretic Plancherel/Fourier object transitions from convergent to divergent behavior. Quote-worthy: "the distribution ν_λ is singular for λ < 1/2, ν_{1/2} is uniform on [-1,1]" — directly analogous to "S_n^(q) converges with constant 7/45 only at q=3, diverges for q ≥ 5."

**(PARALLEL) Anashin & Khrennikov (2009), "Applied Algebraic Dynamics"; Anashin-Khrennikov-Yurova ergodicity-via-van-der-Put-basis line.**
General framework for p-adic dynamical systems of arithmetic maps including Collatz-type. Provides ergodicity / measure-preservation criteria but does not run the q-sweep over qx+1 nor produce a q/3 ratio law. Useful as the broader p-adic-dynamics context Siegel sits inside.

**(PARALLEL) Random matrix β ensembles (Forrester; Dumitriu-Edelman 2002, "Matrix models for beta ensembles", J. Math. Phys. 43, 5830).**
β = 1, 2, 4 are special points (orthogonal/unitary/symplectic) and β as continuous parameter has been studied. Loose parallel: a discrete family of "natural" parameter values plus a continuous interpolation in which the natural points are special. **Weak parallel** — the criticality framing is different; flagged for completeness.

### Gap analysis

The user's *quantitative* finding — that the Plancherel partial-sum ratio S_n^(q)/S_{n−1}^(q) tends exactly to q/3 — does **not appear in the literature reviewed**. The qualitative dichotomy (q=3 convergent, q≥5 divergent) is folkloric (Crandall conjecture) and reproduced by multiple stochastic / Markov-chain / p-adic frameworks (Matthews, Kontorovich-Lagarias, Aicardi, Santos, Siegel). No source surveyed states the closed-form ratio q/3, the constant c = 7/45 at q=3, or the phase-boundary framing in archimedean-Plancherel terms. Tao 2019 uses the archimedean Plancherel apparatus on a Syracuse-3 measure but never sweeps q. Siegel develops the closest analog (a q-indexed Fourier-Stieltjes object) but in q-adic absolute value with a non-vanishing-vs-vanishing dichotomy, not a partial-sum ratio.

**Cleanest existing language to adopt for the user's framing:** the **Bernoulli-convolution template (Erdős–Solomyak–Hochman)** — "a one-parameter family of measures with a critical parameter at which the L²/Plancherel object transitions from convergent (absolutely continuous / finite ratio) to divergent (singular / blow-up ratio); the critical value is the unique fixed point of the natural ratio." That is the direct dynamical-systems analogue and is the language the qx+1 result should be framed in. **Closest prior work in qx+1 itself is Siegel (FRAMEWORK / DIRECTLY ADJACENT)** — the user's finding is potentially statable as an archimedean-Plancherel companion to Siegel's q-adic Fourier-Stieltjes Correspondence Principle, with q=3 as the unique parameter where archimedean S_n converges. The combined package would be a novel quantitative phase-boundary statement for the qx+1 family.

The framing **q=3 as the unique critical / phase-boundary odd prime in qx+1, detected by archimedean Plancherel partial-sum ratio q/3** appears to be novel.

### Local PDFs (in this folder)

| File | Source | Size |
|---|---|---|
| `KontorovichLagarias2009_3x1_5x1_stochastic.pdf` | [arXiv:0910.1944](https://arxiv.org/abs/0910.1944) | 824 KB |
| `Santos2021_qn1_general.pdf` | [arXiv:2005.00346](https://arxiv.org/abs/2005.00346) | 310 KB |
| `Matthews_generalized_3x1_markov.pdf` | [numbertheory.org](http://www.numbertheory.org/PDFS/matthews-final-revised.pdf) | 426 KB |
| `Siegel_pq_adic_Fourier_part2.pdf` | [arXiv:2208.11082](https://arxiv.org/abs/2208.11082) (OA mirror of Springer p-Adic Numbers, paywalled at https://link.springer.com/content/pdf/10.1134/S2070046625020062.pdf) | 746 KB |
| `Siegel2026_HydraMap_Numen.pdf` | [arXiv:2601.17030](https://arxiv.org/abs/2601.17030) | 670 KB |
| `Siegel2024_pq_adic_Collatz_consolidated.pdf` | [arXiv:2412.02902](https://arxiv.org/abs/2412.02902) | 3.8 MB |
| `Siegel2019_Conservation_Singularities_Dreamcatchers.pdf` | [arXiv:1909.09733](https://arxiv.org/abs/1909.09733) | 814 KB |
| `Tao2019_almost_all_collatz.pdf` | [arXiv:1909.03562](https://arxiv.org/abs/1909.03562) | 2.0 MB |
| `PeresSchlagSolomyak_sixty_years_bernoulli.pdf` | [u.math.biu.ac.il](https://u.math.biu.ac.il/~solomyb/RESEARCH/sixty.pdf) (Solomyak's Bar-Ilan mirror; Yale URL was dead) | 413 KB |
| `Solomyak_Bernoulli_notes.pdf` | [u.math.biu.ac.il](https://u.math.biu.ac.il/~solomyb/RESEARCH/Bernotes.pdf) | 542 KB |
| `Ruelle_dynamical_zeta_transfer_operators.pdf` | [ams.org Notices](https://www.ams.org/notices/200208/fea-ruelle.pdf) | 155 KB |
| `KontorovichSinai2002_dgh_maps_structure.pdf` | [arXiv:math/0601622](https://arxiv.org/abs/math/0601622) (Bull. Braz. Math. Soc. 33(2), 213-224, 2002) | 155 KB |
| `Volkov2006_5k1_probabilistic.pdf` | [Volkov's Lund preprint](https://www.maths.lth.se/matstat/staff/s.volkov/PAPERS/5k+1.pdf) (Stoch. Proc. Appl. 2006) | 270 KB |
| `Sinai2003_statistical_3x1.pdf` | [arXiv:math/0201102](https://arxiv.org/abs/math/0201102) (Comm. Pure Appl. Math. 56, 1016-1028, 2003) | 151 KB |
| `Lagarias1987_3x1_and_generalizations.pdf` | Lagarias, "The 3x+1 problem and its generalizations" (Amer. Math. Monthly 92, 3-23, 1985; widely cited as the 1987 reprint in Conway et al. *Organic Mathematics* / *Selected Papers* anthology) — **local copy provided by user** | 2.1 MB |

### Reading focus (per user note 2026-05-04)
Read **Kontorovich-Sinai (d,g,h)-Structure Theorem** + **Volkov 5k+1 probabilistic model** before deciding whether the K-S transfer-operator framework actually has the machinery to attack the rate-1/2 problem in the Q-sweep finding.


