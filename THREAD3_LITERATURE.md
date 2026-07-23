# Thread-3 analytic step — literature map (novel / robust / derivable)

**Date:** 2026-07-22. Consolidates four independent literature hunts on the single remaining analytic step of the
`S_∞ = 7/15` program, integrated with the **W2/R81 reframe** (the asymptotic rate is the period-9 ~0.984 mode, not
½ — see `results/result_W2_tailor.md`). The step, in its rate-independent form: **`S_∞ = 7/15 ⟺ Σ_r |A_r(m)| < ∞`
on the window `m ≲ r`** — a geometric-rate second-moment (Plancherel) Fourier-decay / conditional-digit
equidistribution statement on `1+3ℤ₃`, in the **non-resonant** regime (2 a primitive root mod 3^r ⟹ ⟨2⟩ = full
unit group). Summability needs only rate `< 1`; the exact rate (½ subcritically, ~0.984 asymptotically) does not
change *whether* it's summable, only *which mode* controls it.

---

## 1. NOVEL? — Yes, decisively (all four hunts converge)

- `7/15`, the primitive-shell Plancherel mass `S_r`, the layer/`Λ_r` theory, and the spectral picture are **not in
  the literature** — checked Tao (2019/2022 + the 2020 preimages blog), Kontorovich–Lagarias, Wirsching (LNM 1681),
  Krasikov–Lagarias, Assani–Ebbighausen (2023), and Siegel's non-archimedean Collatz program.
- Sharper: the **geometric-rate decay is an *unproven conjecture of Tao's*** — his `c_n = 3^{−n+o(n)}` (the "β=1"
  hypothesis, 2020 blog) — for *our exact object* `Syrac(ℤ/3ⁿ)`. He proved only **superpolynomial** decay
  `E[e^{−2πiξ·Syrac/3ⁿ}] ≪_A n^{−A}` (uniform, 3∤ξ). A geometric-rate result **sharpens the current frontier**, it
  is not a corollary of it.
- Closest adjacent program: **M. Siegel, "The Collatz Conjecture & Non-Archimedean Spectral Theory"** (arXiv
  2007.15936, 2208.11082; + "Hydra Map / Numen" 2601.17030) — the only other 3-adic/(p,q)-adic Fourier program on
  Collatz; a **novelty-boundary citation and a p-adic-Fourier toolbox**, but it has no 7/15 and no `|λ₂|` gap.

## 2. ROBUST? — Yes; the one real red flag is provably neutralized

- **The Pisot red flag is genuine.** The *archimedean* sibling of our measure (contraction ⅓ on ℝ) is the
  **middle-third Cantor measure, which has no Fourier decay** (non-Rajchman, because 3 is Pisot — Varjú–Yu
  arXiv:2004.09358, Brémont arXiv:1910.03463 make it an iff). Naïvely this predicts our decay *fails*.
- **It provably does not transfer.** Our contraction-⅓ comes from the 3× lift, but the *multiplier* `2^{−v}` is a
  **3-adic unit** generating the **full** unit group (Gauss: 2 primitive root mod 3^r). That is the
  **decay-forcing (dense, non-arithmetic)** side of the dichotomy — opposite the Cantor measure. Exact overlap is
  benign here: branches **tile** the shells, they don't drop dimension. Two hunts confirm from different directions
  (Pisot dichotomy; and Konyagin small-subgroup non-cancellation is a *small*-|H| pathology, inapplicable to the
  full group). **No obstruction; 7/15 and gap-survival are the expected outcome.**
- Our **"no finite transfer operator" proof (R29) is consistent with the state of the art**, not a dead end: the
  archimedean theory *also* abandoned finite transfer matrices (Garsia, Alexander–Zagier) for flattening /
  continuous-parameter renewal, for exactly the growing-spaces reason.
- **W2/R81 caveat folded in:** the four numerical lines put a mode at ½, but R81 showed **½ is the
  subcritical/transient mode; the asymptotic rate is the slower period-9 ~0.984** (critical-only). Both are `< 1`,
  so **summability — hence the theorem — is unaffected**; the reframe only re-points the *value*-derivation.

## 3. DERIVABLE? — Yes, with leads; several survive the no-finite-operator obstruction

Ranked, each naming the theorem/method and what **we** must supply. Re-pointed at the actual target (summability
of `A_r(m)`; the dominant mode is the period-9 ~0.984, but summability needs only rate `< 1`).

1. **Bourgain–Gamburd L²-flattening / sum-product in ℤ/3ⁿ (top).** Read `S_r` as an L²-mass; each halving step
   flattens the primitive-shell L² mass by a fixed factor ⟹ geometric decay. It is an **inequality, not a spectral
   recurrence — bypasses the no-finite-operator obstruction entirely.** **Baker–Khalil–Sahlsten, "Fourier decay
   from L²-flattening" (arXiv:2407.16699, 2024)** proves `L²-flattening ⇒ sup-norm decay`, so **our L² target is
   the *easy* half** (their input). Supply: non-concentration of the digit measure (Geom(½) tail gives it cheaply),
   full-unit-group driver (have, Gauss), a sum-product/expansion input in ℤ/3ⁿ (p-adic sum-product,
   arXiv:1602.00400).
2. **Sharpen Tao's characteristic-function estimate (arXiv:1909.03562).** Tao already built the 3-adic Fourier
   scaffold (2-D renewal × triangles) for our exact object; the geometric-rate improvement is his own open
   conjecture. **Crux (Hunter 2):** certify the **frequency-multiplicity of `A_r(m)`** — if `O(poly r)` low modes
   carry it (which the mode structure suggests), Tao's `n^{−A}` already **over-satisfies** `O(r^{−1−δ})` and the
   step is closed modulo bookkeeping; if the full ~3^r shell is needed with no cancellation, it is not.
3. **Chung–Diaconis–Graham non-concentration (abelian character sums).** Our linear part is abelian; the R(d)
   operator is a multiplier in the multiplicative-character basis (infinite but structured — the object R29 proved
   non-finite). Supply the single uniform estimate **`sup_{χ≠1} |E_v[χ(2^{−v})]| ≤ ρ < 1`** (a lacunary/geometric
   character-sum bound; primitive-root gives the qualitative backbone) ⟹ `|λ₂| ≤ ρ`. Cites: CDG (Ann. Probab. 15,
   1987); **Eberhard–Varjú (arXiv:2003.08117)** tie the `a·x+b mod q` gap to Bernoulli-convolution entropy +
   Diophantine data.
4. **Chains with complete connections / g-measures (matches R30-D exactly).** R30-D measured "prefix-influence
   decays with depth" (asymptotically depth-Markov). **Fernández–Maillard (arXiv:math/0305026)**: geometric decay
   of the variation coefficients `var_k` ⇒ geometric loss-of-memory ⇒ spectral gap (Doeblin–Fortet). Supply
   `var_k ≤ Cγ^k`. The cleanest formalization of the R30-D obstruction.
5. **Li's continuous-parameter twisted renewal operator.** Fourier decay of stationary measures via a renewal
   theorem with a *continuous* spectral parameter — **explicitly not a finite recurrence**, matching our
   growing-spaces finding. Cites: Li, "Decrease of Fourier coefficients of stationary measures" (Math. Ann.);
   Li–Sahlsten–Stevens.
6. **Named character-sum instruments (band-ℓ¹, the R79-flagged missing tool).** **Postnikov–Korobov** character sums
   mod 3^r (p-adic-log linearization of χ over the smooth ⟨2⟩ orbit; arXiv:1605.07553) for a citable geometric
   exponent; **Cochrane–Zheng** p-adic stationary phase (`p^{m/2}` at non-degenerate critical points) as an
   independent confirmation that cancellation is expected.

**Wrong tools (explicitly excluded, Hunter 3).** Dolgopyat cocycle method (Algom–Hertz–Wang arXiv:2306.01275) and
BFLM/He–Lakrec–Lindenstrauss require **non-linearity / non-abelian** linear parts; our system is **affine +
abelian**, so it mixes by *arithmetic non-resonance*, not oscillatory cancellation. First-moment equidistribution
of ν is meanwhile **trivially solved** (Diaconis–Freedman contraction, Hennion–Hervé quasi-compactness; rate ⅓ =
the metric contraction) — but that is the first moment; the theorem is second-moment.

## 4. The value puzzles, placed

- **cos(π/q) and the ½.** The single-digit Bernoulli factor `|½(1+e^{2πi/q})| = cos(π/q)` gives `q=2→0`, `q=3→½`,
  `q≥4→1`, colliding with the algebraic `λ_c=(q−1)/(q+1)=½` only at q=3. **But W1/W2 showed this ½ is the
  fixed-λ/first-moment/transient value** — the actual random-`v` (Geom ½) renewal has the ½ only as a subcritical
  sub-mode; the asymptotic rate is the period-9 ~0.984. So `cos(π/3)=½` *explains the transient*, not the theorem's
  rate. (Hunter 3's "|λ₂| = renewal tail rate = 1−p under Geom(p)" is the same first-moment story, and its falsifier
  — vary p, see if the *transient* ½ tracks 1−p — is worth one exact run.)
- **The asymptotic 0.984 (period-9).** Critical-only, absent subcritically, so the archimedean self-similar-measure
  literature (which is about the stationary/subcritical object) does not name it directly. Its home is the corpus's
  own **inter-level operator** (T_lead within-level spectrum `{43/45, 0}`; the period-9.2 sign oscillation) — a
  Pisot–Vijayaraghavan/log₃-over-log₂ Diophantine-surface flavor, still the pen's to pin.

## 5. One-line answers

- **Novel?** Yes — 7/15, the S_r theory, and the geometric-rate decay are new; the sharp decay is an unproven Tao
  conjecture for our exact object.
- **Robust?** Yes — the Pisot red flag is real but provably neutralized by the unit-multiplier / full-orbit
  non-resonance; no obstruction, and the ½→0.984 reframe leaves summability (the theorem) intact.
- **Derivable once and for all?** Not yet a theorem, but **inside reach of one machine** — Bourgain–Gamburd
  L²-flattening in ℤ/3ⁿ (our L² target is the easy half, per Baker–Khalil–Sahlsten) — with the decisive shortcut
  being the **frequency-multiplicity of `A_r(m)`** (if `O(poly r)`, Tao's existing bound already closes it).

## 6. Citations (author, year, id, relevance)

- T. Tao 2022, *Almost all orbits of the Collatz map…*, Forum Math. Pi 10:e12, **arXiv:1909.03562** — Syrac(ℤ/3ⁿ)
  = our ν; superpolynomial char-fn decay `n^{−A}`; sharp rate left open.
- T. Tao 2020 blog, *Equidistribution of Syracuse random variables…* — `c_n`, conjecture `c_n=3^{−n+o(n)}` (= our
  target), unproven.
- P. Varjú & H. Yu 2022, *Fourier decay of self-similar measures…*, **arXiv:2004.09358**; J. Brémont 2021,
  *Self-similar measures and the Rajchman property*, **arXiv:1910.03463** — Pisot ⇔ non-Rajchman (the red flag).
- S. Baker, O. Khalil, T. Sahlsten 2024, *Fourier decay from L²-flattening*, **arXiv:2407.16699** — `L²-flattening ⇒
  sup-norm decay` (our L² is the easy half). p-adic sum-product: **arXiv:1602.00400**.
- Chung–Diaconis–Graham 1987, Ann. Probab. 15; Eberhard–Varjú 2021, **arXiv:2003.08117** — abelian `a·x+b` gap ⇔
  non-resonance/Diophantine.
- Fernández–Maillard 2005, **arXiv:math/0305026** — chains with complete connections; decaying past-influence ⇒
  geometric mixing (the R30-D route).
- J. Li (Math. Ann.), Li–Sahlsten–Stevens — continuous-parameter renewal Fourier decay (no finite recurrence).
- Mauduit–Rivat 2010 (Ann. Math. 171), Drmota–Mauduit–Rivat 2011, Bassily–Kátai 1995 — digit-equidistribution with
  geometric per-block savings (the classical template).
- Postnikov–Korobov char sums mod 3^r, **arXiv:1605.07553**; Cochrane–Zheng (Illinois J. Math. 2000) — named
  geometric/`p^{m/2}` character-sum instruments (band-ℓ¹).
- Diaconis–Freedman 1999 (SIAM Rev. 41); Hennion–Hervé LNM 1766 (2001) — first-moment equidistribution solved by
  contraction (rate ⅓), no finite reduction.
- M. Siegel, **arXiv:2007.15936 / 2208.11082 / 2601.17030** — non-archimedean Collatz Fourier program
  (novelty-boundary + p-adic toolbox).
- Excluded (wrong tool): Algom–Hertz–Wang **arXiv:2306.01275** (Dolgopyat needs non-linearity); BFLM (Invent. Math.
  187, 2012) + He–Lakrec–Lindenstrauss **arXiv:2003.03743** (need non-abelian) — framing only.
