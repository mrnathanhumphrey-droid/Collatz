# LIT DISPOSITION — is the denominator-growth Mahler-ORDER criterion (λ_L → 2) novel? **General principle = COROLLARY; sharp λ_L→2 = SHARP-VERSION-NOVEL but UNVERIFIED (likely folklore/companion-matrix corollary).** (2026-09-06)

Places the criterion extracted from `result_MAHLER.md` in the literature. Spec: Wilson (kill-first lit-check).
Two *separate* dispositions, as the guardrails demand — the general principle vs the sharp order-indexed ceiling.

## The claim being placed
CRITERION: a rational Mahler equation of order `L`, `Σ_{i=0}^{L} A_i(z)·f(z^{q^i}) = B(z)` with `A_i,B ∈ ℚ(z)`,
forces the log-denominator growth ratio of `f`'s coefficients `≤ λ_L`, where `λ_L` is the dominant root of
`x^L = x^{L−1}+…+1` (the "(L+1)-bonacci" constant): `λ_1 = φ ≈ 1.618`, `λ_2 ≈ 1.839` (tribonacci),
`λ_L → 2⁻`. Hence **log-denominator ratio `> 2` ⟹ no finite-order rational Mahler equation** (infinite order),
with NO zero-estimate. Our `S`: diagonal denominators `2^{2·3^{r−1}}−1`, log-denom ratio = `q = 3 > 2` ⟹ infinite order.

---

## DISPOSITION 1 — the GENERAL PRINCIPLE ("big coefficient denominators ⟹ no finite-order rational Mahler"): **COROLLARY.**

This direction is classical and its Mahler-specific form is a **theorem in the current literature**. The prior was
correctly COROLLARY, and it holds:

- **Eisenstein (1852) lineage — controlled denominators of algebraic series.** An algebraic `f ∈ ℚ[[z]]` has
  `f ∈ ℤ[1/S][[z]]` for a single `S`: denominators are built from *finitely many primes* and grow at most
  **geometrically** (single common denominator `D` with `D^n a_n ∈ ℤ`). Explicit quantitative versions exist
  (Dwork–Robba, Schmidt, van der Poorten; recent explicit form arXiv:2412.11738). Contrapositive = the principle.
- **D-finite analogue.** A D-finite power series has denominators divisible by `≤ K·ℓ·ln ℓ` primes (Eisenstein-type
  bound). Same shape: "nice equation ⟹ controlled denominators," so uncontrolled denominators ⟹ not nice.
- **Mahler-specific — the direct references.**
  - **Mahler's own denominator bound `Δ_n`** (induction on `n` from the functional equation): in the *linear*
    case (`m = 1`) the exponent of each factor is **independent of `n`**; only the nonlinear case (`m ≥ 2`) grows
    with `n`. This is exactly "finite-order linear Mahler ⟹ controlled denominators" (see arXiv:2502.09259 for the
    modern write-up of `Δ_n`).
  - **Bell–Chyzak–Coons–Dumas, "Becker's conjecture on Mahler functions"** (arXiv:1802.08653; Trans. AMS 372
    (2019) 3405–3423): `k`-regular `F` ⟹ `F/R` satisfies a Mahler equation with `a_0(z)=1` and `R = z^γ Q(z)`
    explicit — the structural controlled-denominator statement for Mahler functions.
  - **Dumas** (thesis, *Récurrences mahlériennes, suites automatiques*): the classical structure/asymptotics of
    Mahler recurrences underpinning all of the above.

⟹ **We do NOT claim novelty for the general principle.** Its load-bearing use for `S` (denominators too big ⟹ no
finite-order rational Mahler eq) is a **corollary** of the Eisenstein/D-finite/Mahler controlled-denominator
lineage; cite Adamczewski–Bell–Smertnig + Bell–Chyzak–Coons–Dumas and drop any novelty flag on the principle.

---

## DISPOSITION 2 — the SHARP `λ_L → 2` CEILING (order-indexed, L-bonacci, sup exactly 2): **SHARP-VERSION-NOVEL, but UNVERIFIED — treat as likely folklore / routine companion-matrix corollary.**

**What the literature's nearest quantitative refinement actually is — and why it is a DIFFERENT shape.**

The closest existing object is the **height gap theorem** (the state of the art for "how fast can Mahler-function
coefficients' arithmetic complexity grow"):

> **Adamczewski, Bell & Smertnig, "A height gap theorem for coefficients of Mahler functions"** (arXiv:2003.03429,
> J. Eur. Math. Soc.), abstract verbatim: *"We study the asymptotic growth of coefficients of Mahler power series
> with algebraic coefficients, as measured by their logarithmic Weil height. We show that there are five different
> growth behaviors, all of which being reached. Thus, there are gaps in the possible growths."*
> Two of the five anchors are identified in the follow-up literature: `k`-automatic ⟺ coefficients in a **finite
> set** (bounded height); `k`-regular ⟺ height **`O(log n)`**. (Extended to Hahn series by Faverjon–Roques,
> "Hahn series and Mahler equations: Height Gap Theorem.")

The critical structural facts for our verdict:
1. **It is a FIVE-CLASS GAP theorem, not a continuum.** The achievable growths are five discrete regimes with
   *forbidden gaps between them* — qualitatively unlike our `λ_L` **continuum** `φ, 1.839, 1.928, … → 2`.
2. **It is order-INDEPENDENT.** The five classes are intrinsic growth regimes of the *function*; they are **not
   indexed by the order `L`** of the Mahler equation. Our statement is precisely an *order-indexing* (`λ_L` rises
   with `L`) — a refinement along an axis the height-gap theorem does not use.
3. **No bonacci constant and no "2" as a threshold appear** in the abstract or the searchable statements. The named
   anchors are `O(1)` and `O(log n)`, not `φ`/tribonacci/`2`.

So the sharp `λ_L → 2` statement, **as stated**, was **not found** in the literature. Under the spec's taxonomy that
is SHARP-VERSION-NOVEL. **But three skeptic's caveats keep this well short of a novelty claim:**

- **(a) I could not read the two most relevant papers' internals.** Both `2003.03429` and the Faverjon–Roques
  height-gap paper returned as undecodable binary PDFs; I have only their abstracts + secondary snippets. One of
  the five height-gap classes could coincide with a rate-2 boundary, or their internal proof could contain the
  companion-matrix bound below. **"Not found" here = "not found in abstracts + search snippets," not "read cover to
  cover."**
- **(b) `λ_L` = spectral radius of the order-`L` companion matrix — a standard object.** An order-`L` Mahler
  recurrence expresses `a_n`'s denominator from `a_{n/q}, …, a_{n/q^L}`, so `D(n) := log(den a_n)` obeys
  `D(n) ⪅ D(n/q)+…+D(n/q^L)`, whose growth exponent is the dominant root of `x^L = x^{L−1}+…+1` — literally the
  `(L+1)`-bonacci constant / spectral radius of the all-ones companion matrix. This is a **routine linear-recurrence
  computation**; an expert may well regard the `λ_L → 2` ceiling as *immediate* from the recurrence's shape, i.e.
  folklore, even though it is not written in these words. (The searches surfaced exactly this companion-matrix
  characteristic-polynomial framing `p_c(λ)=λ^d − c_1λ^{d−1} − …` in the adjacent Mahler literature.)
- **(c) The load-bearing use does not need the sharp version.** For `S` we only need "ratio `= 3 > 2 = sup λ_L` ⟹
  infinite order," and DISPOSITION 1 already delivers that as a corollary of controlled denominators. The sharp
  `λ_L` sequence is a *decoration* on a conclusion the general principle already gives; nothing downstream (MAHLER,
  the capstone) rests on `λ_L` being new.

⟹ **Verdict for the sharp version: SHARP-VERSION-NOVEL as literally stated, but UNVERIFIED and most likely
folklore/COROLLARY.** Recommended posture: **do NOT assert public novelty for `λ_L → 2`** without an expert check
(Lean Zulip / Adamczewski-circle), because (b) makes it the kind of statement an expert calls "obvious." Use it
internally as a clean quantitative *illustration* of the general principle; if it is ever put in front of experts,
frame it as "the order-`L` companion-matrix / `(L+1)`-bonacci refinement of the controlled-denominator bound —
is this written down anywhere?" and let them place it.

---

## VERDICT SUMMARY
| Claim | Verdict | Action |
|---|---|---|
| General principle: big denominators ⟹ no finite-order rational Mahler | **COROLLARY** | cite Eisenstein + Adamczewski–Bell–Smertnig + Bell–Chyzak–Coons–Dumas; no novelty flag |
| Sharp `λ_L = (L+1)`-bonacci `→ 2` order-indexed ceiling | **SHARP-VERSION-NOVEL (as stated) — but UNVERIFIED, likely folklore/companion-matrix corollary** | don't claim novelty; expert-check before any public framing; not load-bearing anyway |

**Skeptic's note on "not found":** searched (i) Eisenstein / controlled-denominator lineage, (ii) Mahler-order
denominator/height growth rate, (iii) the height gap theorem (Adamczewski–Bell–Smertnig; Faverjon–Roques Hahn-series
extension), (iv) Bell–Chyzak–Coons–Dumas / Dumas structure theory, (v) explicit "golden ratio / tribonacci /
multinacci / spectral radius / companion matrix" + Mahler order. The general principle is *everywhere*; the exact
`λ_L → 2` order-indexed sequencing is *nowhere I could read* — but the two definitive papers' internals were
un-extractable (binary PDFs), so this is "absent from abstracts + snippets," a **high-but-not-exhaustive** bar. The
default remains COROLLARY-flavored.

**Sources**
- Adamczewski, Bell & Smertnig, *A height gap theorem for coefficients of Mahler functions* — https://arxiv.org/abs/2003.03429 · https://ems.press/journals/jems/articles/5898523
- Bell, Chyzak, Coons & Dumas, *Becker's conjecture on Mahler functions* — https://arxiv.org/abs/1802.08653
- Faverjon & Roques, *Hahn series and Mahler equations: Height Gap Theorem* — https://math.univ-lyon1.fr/~roques/HSMEHGT.pdf
- *Arithmetic properties of the Taylor coefficients of differentially algebraic power series* (modern write-up of Mahler's `Δ_n` denominator bound; linear case exponent independent of `n`) — https://arxiv.org/abs/2502.09259
- *About Eisenstein's Theorem* (explicit quantitative controlled-denominator bound for algebraic series) — https://arxiv.org/abs/2412.11738
