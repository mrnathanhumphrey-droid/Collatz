# Literature defensibility note — prefix decomposition

A 30-minute due-diligence pass on the obvious places someone looking for prior
art on the prefix-decomposition observation (residue class mod 2^k → symbolic
state a_final · m + c_final with a_final ∈ {3^j : 1 ≤ j ≤ k}) would check.

Goal: be able to say honestly to Bonacorsi "I checked X, Y, Z and didn't find
this specific decomposition; if it's there I missed it."

## What's already in the literature

**Terras (1976), "A stopping time problem on the positive integers", Acta
Arithmetica 30, pp. 241-252.** The canonical predecessor. Two relevant lemmas
(via the standard exposition at https://www.ericr.nl/wondrous/terras.html):

- **Lemma 3:** The parity vector P_k(n) of the first k Collatz iterates depends
  only on n mod 2^k. Equivalently, P_k(M) = P_k(N) iff M ≡ N (mod 2^k).
- **Lemma 4 (the Terras approximation):** After k iterations,
  S_k ≈ S_0 · 3^d(k) · 2^(-k), where d(k) is the number of odd steps in the
  first k iterations of S_0.

Lemma 4 is the *asymptotic-approximation* version of our prefix decomposition.
The factor 3^d(k) · 2^(-k) is, up to relabeling, the same a_final / 2^k that
appears in our writeup.

**Differences from what we use:**

1. **Stopping criterion.** Terras's Lemma 4 fixes exactly k steps. Our prefix
   stops at the first step where the symbolic multiplier a becomes odd, which
   is generally < k and varies by residue class (so prefix length itself is a
   per-class deterministic function of r, not a constant). This is what gives
   a_final ∈ {3^1, ..., 3^k} rather than a single 3^d(k) for fixed k.
2. **The additive term.** Terras's approximation S_k ≈ S_0 · 3^d(k)/2^k drops
   the additive correction. We retain it explicitly as c_final, because
   empirically it modulates the conditional distribution σ(n) | n mod 2^k at
   second order (~3× smaller effect than the cross-cluster a_final effect, but
   real and detectable — see experiments/05_cfinal_ks_analysis.py and
   experiments/07_anderson_darling.py).
3. **What we do with it.** Terras uses Lemma 4 to prove almost-sure finite
   stopping (CLT on the parity vector). We use the *exact* symbolic state
   (a_final, c_final) as a *covariate* parameterizing the conditional
   distribution σ | n mod 2^k, which collapses 2^(k−1) odd residue classes
   onto k primary clusters indexed by a_final.

## What I could not find anywhere

No reference that:
- writes the exact symbolic state (a_final · m + c_final) with variable prefix
  length (stop at first odd a) explicitly, or
- uses (a_final, c_final) as a covariate / hierarchy level for modeling the
  conditional distribution σ(n) | n mod 2^k, or
- cites the σ | a_final clustering in any form.

## Where I checked

Time-boxed search across the standard places:

- **Lagarias, "The 3x+1 problem: An annotated bibliography (1963–1999)",
  arXiv:math/0309224.** Abstract/metadata only via WebFetch — could not
  full-text-grep the annotated entries from the public landing page. Cross-
  checked the canonical results indexed elsewhere.
- **Lagarias, "The 3x+1 problem: An annotated bibliography II (2000–2009)",
  arXiv:math/0608208.** Same caveat as above. Worth a manual grep of the PDFs
  on Bonacorsi's side if he wants belt-and-suspenders.
- **Sinai (2003), "Statistical (3x+1) problem", Comm. Pure Appl. Math 56,
  1016–1028, arXiv:math/0201102.** Probabilistic / CLT framework on the
  Syracuse map T : odd → odd via T(x) = (3x+1)/2^k. Treats the geometric
  distribution of the 2-adic valuation k(x), not the residue-mod-2^k symbolic
  state. Not the same decomposition.
- **Korec (1994), "A density estimate for the 3x+1 problem", Mathematica
  Slovaca 44(1), 85–89.** Asymptotic density result (almost all orbits dip
  below n^c for c > log_4 3). Density bound, not symbolic state analysis.
- **Tao (2019/2022), "Almost all Collatz orbits attain almost bounded values"
  + accompanying blog post.** Syracuse-map framework with 3-adic Fourier
  analysis on Z/3^n Z. Not the residue-mod-2^k symbolic decomposition; in fact
  Tao explicitly works in the 3-adic setting that's complementary to ours.
- **Stérin, "The Collatz process embeds a base conversion algorithm",
  arXiv:2007.06979.** Quasi-cellular-automaton on base-2/base-3 digit grids.
  Different symbolic framework (digit-wise, finite-state), not residue-mod-2^k
  affine evolution.
- **Stérin / Khan / others, recent finite-state automaton models
  (e.g., arXiv:2506.21728, arXiv:2104.12135, etc.).** Various base-digit and
  tree-topology frameworks. None expose a_final / c_final as an explicit
  covariate for stopping-time distributions.
- **Wikipedia "Collatz conjecture" article.** Quotes Terras Lemma 3 ("parity
  sequences agree in first k terms iff M ≡ N mod 2^k") but does not present
  the exact (a_final, c_final) decomposition or its application to
  conditional distributions.
- **Bonacorsi & Bordoni (2026), arXiv:2603.04479.** Cites Terras 1976,
  Lagarias 1985, Sinai 1992, Tao 2022. Their hierarchical model uses 32 mod-64
  random effects without referencing a structural decomposition for those
  classes.

## Honest verdict

The **asymptotic identity** S_k ≈ S_0 · 3^d(k) / 2^k is **Terras 1976,
Lemma 4**. We should — and do — cite it as the predecessor.

The **exact symbolic identity** state = a_final · m + c_final with a_final ∈
{3^j : 1 ≤ j ≤ k} (variable prefix length, c_final retained), and the
**modeling consequence** that this collapses 2^(k−1) odd residue classes onto
k primary distributional clusters parameterized by a_final, are what we
believe is novel — at least in the form presented. Both are *elementary*
refinements of Terras: anyone could derive them. But "elementary" and
"explicitly written down somewhere I can find" are different things, and
across a focused search of the canonical 3x+1 literature I did not find this
combination.

If a colleague (Bonacorsi, Lagarias, anyone with deeper survey memory of the
50-year corpus) recognizes this from somewhere I missed, that is the natural
correction to make. I'd rather be told the reference than miss it.
