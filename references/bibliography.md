# Wiener-Hopf for W_j — annotated bibliography

Fetch run: 2026-05-02. Goal: close the third constant in a Bayesian-Tao-bridge analysis of the Collatz conjecture via Wiener-Hopf factorization for W_j.

Legend: ✓ downloaded | ✗ paywalled (note method)

---

## Tier 1 — Core framework (read first, in order)

### ✓ Denisov & Wachtel (2026), "Fluctuations of Discrete-Time Random Walks" — `DenisovWachtel2026_fluctuations.pdf` (732 KB)
- arXiv:2602.18081 — https://arxiv.org/pdf/2602.18081
- **Why:** Most directly applicable. Discrete-time, one-dimensional, covers Wiener-Hopf factorization AND the newer "universality approach" which is more robust. Exactly the setting (unbounded log-walk, first passage). **Primary reference**, not Doney.

### ✗ Doney (2007), "Fluctuation Theory for Lévy Processes" (Saint-Flour XXXV)
- Springer LNM 1897. ~150 pages.
- Springer paywall: https://link.springer.com/book/10.1007/978-3-540-48511-7
- **Status:** No open-access preprint located (author homepage URLs 404'd; ResearchGate has metadata-only "Request PDF" entry; Manchester staff pages returned 404). Only pirate mirrors found in search — skipped per legal posture.
- **How to obtain:** Institutional library access via Springer Link, ResearchGate "Request PDF" from Doney directly, or interlibrary loan.
- **Why:** The named reference. Needed for the rigorous Wiener-Hopf identities and ladder variable theory.

### ✓ Kyprianou (2014), "Introductory Lectures on Fluctuations of Lévy Processes with Applications" (2nd ed.) — `Kyprianou2014_levy_fluctuations.pdf` (3.5 MB)
- Springer textbook. ~450 pages. **Chapter 6 is the Wiener-Hopf chapter.**
- Source: https://public.econ.duke.edu/~get/browse/courses/883/Spr16/COURSE-MATERIALS/Z_Papers/Ky.pdf (Duke course materials mirror)
- **Why:** Pedagogically clearer than Doney. Best entry point if Denisov-Wachtel feels too terse.

---

## Tier 2 — Foundational identities (the actual machinery)

### ✓ Spitzer (1956), "A combinatorial lemma and its application to probability theory" — `Spitzer1956_combinatorial_lemma.pdf` (1.1 MB)
- Trans. AMS 82, 323-339.
- Source: https://www.ams.org/journals/tran/1956-082-02/S0002-9947-1956-0079851-X/S0002-9947-1956-0079851-X.pdf (AMS open access)
- **Why:** The original Spitzer identity — **the engine of W_j**. Cite for the Spitzer-Baxter identity for the joint distribution of (τ, S_τ).

### ✗ Alili & Doney (1999), "Wiener-Hopf factorisation revisited and some applications"
- Stochastics and Stochastics Reports 66, 87-102.
- Warwick WRAP record: http://wrap.warwick.ac.uk/63692/ — metadata only, **no full-text attached**.
- **Status:** No open-access PDF located.
- **How to obtain:** Institutional access to Stochastics and Stochastics Reports (Taylor & Francis), or request from Doney.
- **Why:** Modern reformulation of Spitzer. Often easier to apply than the original.

### ✗ Doney (1982), "On the exact asymptotic behavior of the distribution of ladder epochs"
- Stoch. Proc. Appl. 12, 203-214.
- Elsevier paywall: https://www.sciencedirect.com/science/article/pii/0304414982900424
- **Status:** Pre-arXiv era; no preprint mirror located.
- **How to obtain:** Institutional ScienceDirect access.
- **Why:** Ladder epoch tails — directly relevant to W_j since first-passage decomposes through ladder structure. Result: for E(X)=0, E(X²)=∞, right tail asymptotically larger than left, P{T+ ≥ n} ~ n^(1/β-1) · L+(n) iff P{X ≥ x} ~ 1/(x^β · L(x)) with 1<β<2.

---

## Tier 3 — Small-target / discrete-set first passage

### ✓ Caravenna & Doney (2019), "Local large deviations and the strong renewal theorem" — `CaravennaDoney2019_local_large_deviations.pdf` (647 KB)
- arXiv:1612.07635 — https://arxiv.org/abs/1612.07635
- **Why:** "Discrete small-target set" requires local (not just integral) limit theorems. This handles the regime.

### ✓ Vidmar (2015), "Fluctuation theory for upwards skip-free Lévy chains" — `Vidmar2015_skip_free.pdf` (579 KB)
- arXiv:1309.5328 — https://arxiv.org/abs/1309.5328
- **Why:** If the m_j structure has any skip-free property (worth checking), this gives explicit scale functions and may make W_j **closed-form** rather than just asymptotic.

### ✓ Kwaśnicki (2019), "Random walks are determined by their trace on the positive half-line" — `Kwasnicki2019_random_walks_trace.pdf` (168 KB)
- arXiv:1902.08481 — https://arxiv.org/abs/1902.08481
- **Why:** Complex-analytic methods for Wiener-Hopf factors. May give cleaner derivation depending on the walk's symbol structure.

---

## Tier 4 — Markov-modulated case (only if needed)

### ✓ Alsmeyer & Buckmann (2018), "Fluctuation theory for Markov random walks" — `AlsmeyerBuckmann2018_markov_random_walks.pdf` (690 KB)
- arXiv:1608.08377 — https://arxiv.org/abs/1608.08377
- **Why:** Extends Wiener-Hopf to Markov-modulated walks. Read only if the iid framework breaks (when the absorbing-chain structure for P(j) makes W_j genuinely Markov-driven rather than iid).

---

## Added during fetch run (user-specified)

### ✓ Kuznetsov, Kyprianou, Pardo (2010), "Meromorphic Lévy processes and their fluctuation identities" — `KuznetsovKyprianouPardo2010_meromorphic_levy.pdf` (1.5 MB)
- arXiv:1004.4671 — https://arxiv.org/abs/1004.4671
- **Note:** Family of Lévy processes with **explicitly computable Wiener-Hopf factorizations** for first-passage problems. Likely directly useful for closed-form constant-closing if the Collatz-derived walk fits the meromorphic class.

---

## Summary

- **Downloaded:** 8 of 11 (8 of the 10 in the original brief, plus 1 user-added)
- **Paywalled / no open-access located:** 3 (Doney 2007 St-Flour, Alili-Doney 1999, Doney 1982)
- **All 4 priority items:** ✓ Denisov-Wachtel ✓ Kyprianou ✓ Caravenna-Doney ✗ Doney 2007 St-Flour (paywalled)

The three paywalled items are all **named foundational references** but not directly required if Kyprianou (2014) Ch. 6 is read first — Kyprianou cites and reproduces the key Spitzer/Alili-Doney/Doney identities.

---

## Postnikov 1956 — character formula for L-series mod p^k

Fetch run: 2026-05-03. One-time legal-research add-on; unrelated to Wiener-Hopf.

### Citation
- A. G. Postnikov (1956), "On Dirichlet L-series with the character modulus equal to the power of a prime number." J. Indian Math. Soc. (N.S.) 20, 217–226.
- Companion (Russian): A. G. Postnikov (1955), "On the sum of characters with respect to a modulus equal to a power of a prime number." Izv. Akad. Nauk SSSR. Ser. Mat. 19, 11–16.

### Primary status: ✗ NOT located
- **J. Indian Math. Soc. Vol. 20 (1956)** is a Silver Jubilee Commemoration Volume. The IMS digitization effort hosts back issues at https://informaticsjournals.co.in/index.php/jims/issue/archive but Vol. 20 (1956) was not directly retrievable in this fetch (WebFetch was permission-blocked on the candidate issue URLs `informaticsjournals.com/index.php/jims/issue/view/14XX`; only Vol. 16 = 1476 was confirmed via web search snippet).
- **HathiTrust** holds JIMS Vols. 1–18 freely readable; **Vol. 20 was not confirmed** in this fetch — would need direct catalog navigation at https://catalog.hathitrust.org/Record/008320850.
- **Internet Archive** (`dli.ernet`) hosts scattered JIMS volumes (2, 8, 9, 11, 14, 15) but no Vol. 20 surfaced.
- **mathnet.ru** Postnikov author page (personid=9568) blocked by WebFetch permission this session; the Russian companion (Izv. SSSR 19, 1955) is the natural mirror but was not retrieved.
- **For institutional access:** any university with HathiTrust full-text rights, or the IMS journal archive direct download. Not pirate-mirrored per legal posture.

### Secondary status: ✓ formula captured via modern reproduction

### ✓ Banks & Shparlinski (2016), "Bounds on short character sums and L-functions for characters with a smooth modulus" — `BanksShparlinski2016_Postnikov_reproduction.pdf` (267.8 KB, **cached but not yet copied to references/**)
- arXiv:1605.07553v2 — https://arxiv.org/pdf/1605.07553
- **Cached at:** `C:\Users\Nate\.claude\projects\c--As-Above-So-Below-Master\16dd9aa3-ce85-4e15-8d3d-5e463abf1666\tool-results\webfetch-1777873234448-riasc3.pdf`  *(Bash/PowerShell copy was permission-blocked this session — please move it manually with `cp <cache> C:\Collatz\references\BanksShparlinski2016_Postnikov_reproduction.pdf`)*
- **Cites Postnikov 1956 directly** as [12] (alongside the 1955 Russian companion as [11]) — the paper's entire approach extends "the method of Postnikov [11,12]."
- **Postnikov's character formula appears here as Lemma 4.1**, attributed via a chain Postnikov 1956 → Gallagher [4, Lemma 2] (1972) → Iwaniec [7, Lemma 2] (1974). The truncated-log polynomial is given just before the lemma in equation (4.1).

### Foundational character formula — verbatim block-quote (Banks–Shparlinski 2016, §4.2):

> **(4.1)** For an integer d ≥ 1 we use F_d to denote the polynomial approximation to log(1 + x) given by
>
>     F_d(x) = Σ_{r=1}^{d} (−1)^{r−1} · x^r / r .
>
> **Lemma 4.1.** Let χ be a primitive character modulo q. Let d be an integer such that q² | q_♯^d, and put
>
>     τ = 2  if 4 | q,
>     τ = 1  otherwise.
>
> Then  **χ(1 + τ q_♯ x) = e(f(x))**,  where f is a polynomial of the form
>
>     f(x) = q^{−1} m · F_d( τ q_♯ x )
>
> with an integer m for which gcd(m, q) = 1, and r | m for every integer r ∈ [1, d] coprime to q.

Notation:
- `q_♯ = ∏_{p|q} p` is the **core/kernel** of q (i.e., the squarefree radical).
- `e(t) = e^{2πit}` (additive character on R/Z).
- `χ` is a primitive Dirichlet character mod q.
- For the prime-power case `q = p^γ` (the original Postnikov setting), `q_♯ = p` and the lemma becomes the direct statement: **for x with vp(x) = 0 small enough that 1+px is a unit, χ(1+px) is the additive character e(f(x)) of a truncated-logarithm polynomial whose coefficients live in (1/q)Z.** That is the foundational Postnikov character formula — multiplicative characters mod p^γ are *exponentials of a truncated p-adic-logarithm-style polynomial*, exactly as the user described.

### Equivalent classical statement (the form Postnikov 1956 actually proved):
For odd prime p, γ ≥ 2, and any primitive Dirichlet character χ mod p^γ, there exists an integer a (depending on χ, with gcd(a, p) = 1) such that for all integers t,

    χ(1 + p t) = e(  a · L_γ(t) / p^γ  ),

where L_γ(t) = p t − (p t)²/2 + (p t)³/3 − ⋯ + (−1)^{γ−1} (p t)^γ / γ is the truncated logarithm — equivalent to F_d above with τ = 1, q_♯ = p, q = p^γ, and m = a. (This is the form quoted in Iwaniec–Kowalski, *Analytic Number Theory*, AMS Colloq. 53 (2004), §12.6, citing Postnikov via Gallagher; Banks–Shparlinski Lemma 4.1 is the smooth-modulus generalization.)

### ✓ Banks & Shparlinski (2015), "Estimates for character sums and Dirichlet L-functions to smooth moduli" — `BanksShparlinski2015_smooth_modulus_Postnikov_reproduction.pdf` (455.8 KB, **cached but not copied**)
- arXiv:1503.07156 — https://arxiv.org/pdf/1503.07156
- **Cached at:** `C:\Users\Nate\.claude\projects\c--As-Above-So-Below-Master\16dd9aa3-ce85-4e15-8d3d-5e463abf1666\tool-results\webfetch-1777873235910-8vkcd4.pdf`
- The earlier paper in the same Banks–Shparlinski sequence; uses the same Postnikov-style reduction. Extracted formula content was lower-quality from the binary fetch — Banks–Shparlinski 2016 above is the cleaner reproduction.

### Other modern sources that explicitly reproduce/use Postnikov's formula (not fetched this run):
- Iwaniec & Kowalski, *Analytic Number Theory*, AMS Colloq. 53 (2004), Theorem 12.16 and §12.6 — textbook reproduction.
- Konyagin & Shparlinski, *Character Sums with Exponential Functions and Their Applications*, Cambridge Tracts 136 (1999), §1.2.
- Gallagher (1972), "Primes in progressions to prime-power modulus." Invent. Math. 16, 191–201 — Lemma 2 is the first English reproof.
- Iwaniec (1974), "On zeros of Dirichlet's L series." Invent. Math. 23, 97–104 — Lemma 2 extends to smooth moduli.
- Milićević (2016), "Sub-Weyl subconvexity for Dirichlet L-functions to powerful moduli." Compos. Math. 152, 825–875.
- Chang (2014), "Short character sums for composite moduli." J. d'Analyse Math. (2014), 1–33.

### Summary
- **Primary (Postnikov 1956 PDF):** ✗ not located in open-access channels this run; institutional Vol. 20 JIMS access required.
- **Secondary reproduction:** ✓ Banks–Shparlinski 2016 (arXiv:1605.07553) PDF cached locally — clean modern restatement of the formula as Lemma 4.1, with explicit truncated-log polynomial F_d, citing Postnikov 1956 directly as ref [12].
- **Formula captured:** ✓ verbatim above.

