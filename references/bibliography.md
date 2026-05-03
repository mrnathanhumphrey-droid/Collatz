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
