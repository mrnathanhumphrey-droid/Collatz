# Varjú Follow-ups: Annotated Bibliography

User question: papers that extend Péter Varjú's "Random walks in compact groups" program (arXiv:1209.1745) to **abelian profinite groups** with **convergence of specific functionals** (not just total-variation distance to Haar). Multiplicative-dynamics-mod-p^k flavor preferred.

Each entry: title, authors, venue, year, arXiv URL (free, legal), one-paragraph summary, and a verdict — **STRONG MATCH** (extends Varjú directly into abelian profinite + functional layer), **TANGENT** (similar technique, different problem), or **BACKGROUND** (foundational reference assumed).

---

## Anchor paper (the one we're branching from)

### Varjú 2012 — Random walks in compact groups
- **Authors:** Péter Pál Varjú
- **Venue:** Documenta Mathematica 18 (2013), 1137–1175
- **arXiv:** https://arxiv.org/abs/1209.1745
- **Summary:** Poly-logarithmic convergence rates of products of i.i.d. random group elements to Haar measure on compact groups; the result is uniform across a family parametrized by quasirandomness. Improves Solovay/Kitaev/Gamburd/Shahshahani/Dinai. Stated for compact groups generally (finite + compact semi-simple Lie); the abelian case is implicit but is not the paper's focus.
- **Verdict:** ANCHOR (this is the paper the user wants follow-ups for, not itself a follow-up).

---

## Strong matches — direct extensions into abelian + multiplicative-mod-p territory

### Eberhard & Varjú 2020 — Mixing time of the Chung–Diaconis–Graham random process
- **Authors:** Sean Eberhard, Péter P. Varjú
- **Venue:** Probability Theory and Related Fields 181 (2021), 1361–1411
- **arXiv:** https://arxiv.org/abs/2003.08117
- **Summary:** Sharp mixing time of X_{n+1} = a·X_n + b_n mod q (a = 2 by default, b_n iid in {-1, 0, 1}) on the abelian group Z/qZ. Mixing time = (c + o(1))·log_2 q with explicit constant c ≈ 1.01136 for almost all odd q. Generalizes to any fixed integer a; mixing rate is governed by entropy of an associated self-similar Bernoulli convolution, tying CDG mixing into Varjú's Bernoulli-convolution program (arXiv:1610.09154 etc.). Convergence is in total variation to Haar (uniform); functional / character-level pieces appear in the proof via Fourier analysis.
- **Why it matters for Collatz:** The map x → ax+b mod q is the canonical abelian-profinite affine random walk; this is the closest existing analog of Syracuse-style multiplicative dynamics on Z/qZ. The constant c being entropy-of-a-Bernoulli-convolution is exactly the kind of arithmetic-combinatorial functional the user is asking about.
- **Verdict:** **STRONG MATCH** (Varjú himself extending his program into abelian Z/qZ, with a functional — the Bernoulli-convolution entropy — as the convergence rate).

### Breuillard & Varjú 2019/2022 — Cut-off phenomenon for the ax+b Markov chain over a finite field
- **Authors:** Emmanuel Breuillard, Péter P. Varjú
- **Venue:** arXiv preprint, last revised 2022
- **arXiv:** https://arxiv.org/abs/1909.09053
- **Summary:** Cut-off (sharp transition) for x_{n+1} = a·x_n + b_n on F_p, for "most" primes p and "most" a ∈ F_p, conditional on GRH for Dedekind zeta of relevant number fields; unconditional upper bounds on mixing time. Companion to Eberhard-Varjú; both push CDG analysis into the abelian Z/pZ regime with sharp constants.
- **Verdict:** **STRONG MATCH** (Varjú extending his program; abelian; affine multiplicative; sharp result, conditional on a number-theoretic hypothesis — feels like a hint the deep machinery is arithmetic-functional).

### Hermon & Olesker-Taylor 2021 — Cutoff for almost all random walks on abelian groups
- **Authors:** Jonathan Hermon, Sam Olesker-Taylor
- **Venue:** Journal of the European Mathematical Society (accepted Sept 2025); arXiv 2021
- **arXiv:** https://arxiv.org/abs/2102.02809
- **Summary:** Cutoff for random walks on random Cayley graphs of finite abelian groups (including Z/p^k Z and arbitrary products of cyclic groups), under k - d(G) ≫ 1 where k is the number of generators and d(G) the minimal generating size. Cutoff time is characterized via the entropy of an associated Z^k random walk. Verifies a 1997 Wilson conjecture for nilpotent groups; extends an Aldous-Diaconis conjecture. Functional layer: the entropy-rate of a Z^k walk is the "functional" governing the convergence rate.
- **Why it matters:** This is the cleanest known statement of "Varjú-style sharp convergence in the abelian case at full generality", and Z/p^k is in scope. Companion technique to Eberhard-Varjú.
- **Verdict:** **STRONG MATCH** (most general abelian extension of the spectral-gap/cutoff thread; covers Z/p^k explicitly).

### Hussain & Lamzouri 2023 — The limiting distribution of Legendre paths
- **Authors:** Ayesha Hussain, Youness Lamzouri
- **Venue:** Journal de l'École polytechnique — Mathématiques, 2024 (vol. 11)
- **arXiv:** https://arxiv.org/abs/2304.13025
- **Summary:** Functional convergence (in the space C[0,1] under Skorokhod-style topology) of the polygonal path of normalized partial character sums Σ_{n≤t·p} χ_p(n)/√p to a random Fourier series, as p ranges over primes in [Q, 2Q] with Q → ∞. Removes the GRH assumption from Hussain's earlier (thesis/2022) version. The random Fourier series limit is built from Rademacher random completely-multiplicative functions.
- **Why it matters:** This is precisely "functional convergence of a specific functional (partial character sum) of an abelian-profinite random walk." It IS the kind of next-layer-down theorem the user is asking about — only the underlying group is Z/pZ rather than Z/p^k Z, and the "randomness" averages over primes rather than over walk steps. Hussain's earlier solo paper (IMRN 2022) is the conditional precursor.
- **Verdict:** **STRONG MATCH** (functional convergence + abelian Z/pZ + multiplicative structure; closest single hit for "functional CLT for multiplicative-character random walk").

### Hussain 2022 — The limiting distribution of character sums (Hussain's solo, GRH-conditional)
- **Author:** Ayesha Hussain
- **Venue:** International Mathematics Research Notices, 2022(20), 16292–16322
- **arXiv:** (preprint at https://people.maths.bris.ac.uk/~zj18371/papers/limitingdistributionsofcharactersums.pdf)
- **Summary:** Same Legendre-path functional convergence theorem as Hussain-Lamzouri 2023, but assuming GRH. Predecessor / sets up the random Fourier series limit object.
- **Verdict:** **STRONG MATCH** (precursor to 2304.13025; cite both).

### Lamzouri & Zaharescu 2011 — Randomness of character sums modulo m
- **Authors:** Youness Lamzouri, Alexandru Zaharescu
- **Venue:** Journal of Number Theory 132 (2012), 1037–1067
- **arXiv:** https://arxiv.org/abs/1104.4957
- **Summary:** Models real character sums probabilistically as a random walk on the additive abelian group Z/mZ; proves the values of certain real character sums are uniformly distributed in residue classes mod m. The setup — character sums = partial sums of an abelian-group-valued process — is the historical entry point for the Hussain-Lamzouri functional convergence story.
- **Verdict:** **STRONG MATCH** (foundational for the Hussain functional-convergence theorem; abelian Z/mZ + functional question).

### Ayyer & Singla 2019 — Random Motion on Finite Rings, I: Commutative Rings
- **Authors:** Arvind Ayyer, Pooja Singla
- **Venue:** Algebras and Representation Theory 23 (2020), 2473–2503
- **arXiv:** https://arxiv.org/abs/1605.05089
- **Summary:** Markov chain on a finite commutative ring R using uniformly random addition and arbitrary (possibly random) multiplication. Eigenvalues/multiplicities of the transition operator computed via character theory of the underlying additive abelian group. **Z/p^k Z (the finite chain ring of length k) gets a dedicated section** with explicit stationary measure and constant mixing time (independent of k!).
- **Why it matters:** This is the rare paper that names Z/p^k Z explicitly and treats both the multiplicative and additive layer simultaneously, with full character-theoretic / functional analysis. Direct relevance to qx+1 / Syracuse on Z/p^k.
- **Verdict:** **STRONG MATCH** (explicit Z/p^k Z, character-theoretic / functional analysis, both layers of structure).

---

## Tangents — same machinery, adjacent question

### Bate & Connor 2018 — Mixing time and cutoff for a random walk on the ring of integers mod n
- **Authors:** Michael E. Bate, Stephen B. Connor
- **Venue:** Bernoulli 24(2) (2018), 993–1009
- **arXiv:** https://arxiv.org/abs/1407.3580
- **Summary:** Random walk on Z/nZ that at each step does either an additive step or a multiplicative jump (probability of jump → 0 as n → ∞). Establishes total-variation pre-cutoff for the full walk, and a true cutoff for the subsampled-at-jump-times process. n is general (not restricted to prime); does not isolate the prime-power case but the result applies.
- **Verdict:** TANGENT (Z/nZ generally, not the (Z/p^k)* multiplicative subgroup specifically; functional layer is implicit in TV-cutoff machinery).

### Shkredov 2021/2023 — On the multiplicative Chung–Diaconis–Graham process
- **Author:** Ilya D. Shkredov
- **Venue:** Sbornik: Mathematics 214(6) (2023), 878–895
- **arXiv:** https://arxiv.org/abs/2106.09615
- **Summary:** Multiplicative version of CDG on F_p: X_{n+1} = X_n with prob 1/2, else X_{n+1} = f(X_n)·ε_n. Mixing time exp(O(log p / log log p)). Uses additive-combinatorial / Sidon-set machinery.
- **Verdict:** TANGENT (multiplicative CDG, but only on F_p, no p^k extension; convergence is TV not functional).

### Chatterjee & Diaconis 2020 — Speeding up Markov chains with deterministic jumps
- **Authors:** Sourav Chatterjee, Persi Diaconis
- **Venue:** Probability Theory and Related Fields 178 (2020), 1193–1214
- **arXiv:** https://arxiv.org/abs/2004.11491
- **Summary:** General framework: alternating a Markov chain with a deterministic move can speed mixing dramatically. Examples include CDG-type processes on abelian groups. Provides one of the abstract frameworks under which Eberhard-Varjú lives.
- **Verdict:** TANGENT (general framework, abelian as one example, no functional convergence per se).

### Hildebrand 1993 — Random Processes of the Form X_{n+1} = a_n·X_n + b_n (mod p)
- **Author:** Martin Hildebrand
- **Venue:** Annals of Probability 21(2) (1993), 710–720
- **JSTOR:** https://www.jstor.org/stable/2244691 (paywalled — see hunt log)
- **Project Euclid:** https://projecteuclid.org/euclid.aop/1176989264 (paywalled)
- **Summary:** Generalizes CDG: both a_n and b_n are random, iid sequences. Gets O((log p)^2) mixing on Z/pZ. Discrete Fourier transform machinery. Foundational reference for the whole thread; cited by Eberhard-Varjú.
- **Verdict:** BACKGROUND (foundational; user likely already has it via Eberhard-Varjú citations).

### Hildebrand 2008 — A lower bound for the Chung-Diaconis-Graham random process
- **Author:** Martin Hildebrand
- **arXiv:** https://arxiv.org/abs/0801.3094
- **Summary:** Matching log_2 p lower bound for the CDG process on Z/pZ.
- **Verdict:** BACKGROUND.

---

## P-adic / continuous-time profinite (orthogonal but listed because user mentioned profinite)

### Pierce & Weisbart 2024 — Brownian Motion in the p-Adic Integers is a Limit of Discrete Time Random Walks
- **Authors:** Tyler Pierce, David Weisbart
- **Venue:** Journal of Statistical Physics 192 (2025), article 41
- **arXiv:** https://arxiv.org/abs/2407.05561
- **Summary:** Discrete-time random walks on Z_p (the additive profinite abelian group) converge weakly, as processes, to a real-time diffusion (p-adic Brownian motion) generated by the Vladimirov-Kochubei operator. **Explicit functional convergence in a profinite abelian group** — Skorokhod-style weak convergence of measures on path space. Authors describe it as "the first example of discrete approximation of Brownian motion in the setting of a profinite group."
- **Verdict:** **STRONG MATCH** (functional convergence + profinite abelian group; orthogonal flavor — additive p-adic diffusion, not multiplicative — but the framework is exactly what the user asked for).

### Weisbart 2024 — p-Adic Brownian Motion is a Scaling Limit
- **Author:** David Weisbart
- **Venue:** Journal of Physics A 57 (2024)
- **arXiv:** https://arxiv.org/abs/2010.05492
- **Summary:** Earlier in the same line; scaling-limit construction of p-adic Brownian motion from discrete walks.
- **Verdict:** TANGENT (same theme, prior to 2407.05561).

### Bakken & Weisbart 2019 — p-Adic Brownian Motion as a Limit of Discrete Time Random Walks
- **Authors:** Erik Makino Bakken, David Weisbart
- **Venue:** Communications in Mathematical Physics 369 (2019), 371–402
- **DOI:** https://doi.org/10.1007/s00220-019-03447-y (publisher paywall; see hunt log for free preprint search)
- **Summary:** Original paper of the Weisbart line: the first construction of p-adic BM as a scaling limit of discrete-time walks.
- **Verdict:** TANGENT / BACKGROUND.

### Cruz-López, Estala-Arias & Murillo-Salas 2016 — A random walk on the profinite completion of Z
- **Authors:** Manuel Cruz-López, Samuel Estala-Arias, Antonio Murillo-Salas
- **Venue:** Statistics & Probability Letters 109 (2016), 130–138
- **DOI:** https://www.sciencedirect.com/science/article/abs/pii/S0167715215003831 (publisher paywall)
- **Summary:** Continuous-time random walk on Ẑ = profinite completion of Z (the canonical abelian profinite group, isomorphic to ∏_p Z_p). Computes the infinitesimal generator and gives recurrence/transience results.
- **Verdict:** TANGENT (additive profinite-Z dynamics, not multiplicative; no functional CLT).

### Albeverio & Karwowski 1994 — A random walk on p-adics: the generator and its spectrum
- **Authors:** Sergio Albeverio, Witold Karwowski
- **Venue:** Stochastic Processes and their Applications 53 (1994), 1–22
- **DOI:** https://www.sciencedirect.com/science/article/abs/pii/0304414994900329 (publisher paywall)
- **Summary:** Foundational Lévy-process-on-Q_p / Z_p construction via Chapman–Kolmogorov; spectrum of the generator. Whole p-adic-diffusion school descends from this.
- **Verdict:** BACKGROUND (foundational reference; Pierce-Weisbart and Cruz-López both cite it).

---

## Quick scoreboard

- 6 STRONG MATCHES: Eberhard-Varjú, Breuillard-Varjú, Hermon-Olesker-Taylor, Hussain-Lamzouri, Hussain (solo), Lamzouri-Zaharescu, Ayyer-Singla, Pierce-Weisbart (8 if you count the two Hussain papers separately and the two Weisbart papers separately).
- 4 TANGENTS: Bate-Connor, Shkredov, Chatterjee-Diaconis, Weisbart 2024, Cruz-López et al.
- 3 BACKGROUND: Hildebrand 1993, Hildebrand 2008, Bakken-Weisbart 2019, Albeverio-Karwowski 1994.

## Honest assessment of the literature gap

The user's *exact* question — "functional convergence (in the Donsker sense) of a multiplicative random walk on (Z/p^k Z)*, with k → ∞" — does not appear to have a direct hit. The closest existing pieces are:

1. **Multiplicative + (Z/p^k)\* + character analysis**: Ayyer-Singla 2019 has the Z/p^k case explicitly, but their convergence is to a stationary measure in TV, not a functional CLT.
2. **Functional convergence + abelian + character sums**: Hussain-Lamzouri 2023 has functional convergence (to a random Fourier series) of a partial-character-sum process, but on Z/pZ (no p^k lift).
3. **Profinite + functional convergence**: Pierce-Weisbart 2024 has Skorokhod-functional convergence on Z_p (profinite abelian), but additive (not multiplicative; so misses the qx+1 angle).

So the gap the user is sniffing — multiplicative dynamics on (Z/p^k)\* with a functional CLT as k → ∞ — is **a relatively underdeveloped corner**. The three building blocks above exist; the user's project would be combining them, not replicating an existing theorem.
