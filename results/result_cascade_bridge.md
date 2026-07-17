# Result — the Syracuse multiplier as a cascade generator, and the turbulence bridge (2026-07-15)

**Probe:** `probe_cascade_bridge_2026_07_15.py` → `experiments_output/cascade_bridge_2026_07_15.json`
**Cost:** $0, local, seconds. **Zero JHU** — turbulence side is banked (`D:/Turbulence/data/processed/probe11A_analysis.json`, CLM-246).
**Prior, stated to lose:** FALSIFY the bridge. **Outcome: FALSIFIED — but not for the reason predicted (see §4, the prior lost twice).**

## Origin
User idea: "build a Collatz snowglobe universe — was the universe built by Collatz dynamics?" Sharpened into a falsifiable
question by noting that **both Collatz and the turbulent cascade are random multiplicative processes**, so the claim reduces to:
*does the turbulent cascade's multiplier law match the one Collatz induces?*

## 1. ★ THE DURABLE POSITIVE (bank this — it is a real, citable fact about Collatz)

**The Syracuse multiplier `W = q/2^v` (v ~ Geom(1/2)) is a legitimate log-geometric compound-Poisson cascade generator.**
All verified in-probe (closed form checked against brute force at q ∈ {−0.5, 2/3, 1, 4/3, 2, 3}):

- **Closed form:** `E[W^s] = q^s / (2^(s+1) − 1)`, converging iff `s > −1`.
- **Energy-conserving:** `E[W] = 3·E[2^−v] = 3·(1/3) = 1` **exactly** — the condition a physical cascade must satisfy.
  ⚠️ **But this is CIRCULAR**: it is a restatement of R5/R8's banked identity `3 = 1/E_Geom[2^−v]`. True by construction.
  It cannot be evidence for a bridge. *(This killed the author's first "intriguing" reading — see §4.)*
- **Kahane–Peyrière non-degeneracy:** `E[W ln W] = ln3 − (4/3)ln2 = 0.1744 < ln2 = 0.6931` ✓ (and < ln3). **Valid generator.**
- **Mean-critical but a.s. decaying:** `E[W] = 1` while `E[ln W] = ln(3/4) = −0.2877 < 0`. A real fact about Collatz.
- **Log-infinitely-divisible:** `ln W = ln3 − v·ln2` with v ~ Geom(1/2) = NegBinom(1,1/2), which **is ID** (Feller: ID on ℕ ⟹
  compound Poisson). An affine map of an ID law is ID. **⇒ satisfies Novikov's scale-similarity constraint; lands inside the
  Barral–Mandelbrot compound-Poisson class.**
- **Structural sibling of She–Levêque:** Collatz is `β^N` with **N ~ Geometric**; She–Levêque is `β^N` with **N ~ Poisson**.
  Collatz is the *log-geometric* cascade next door to the field's most celebrated *log-Poisson* one.

**This fits the register R8 already established** (π_k as a q-adic self-similar measure; Bernoulli-convolution framing per the
2026-05-04 lit dive). It is a free structural remark for the qx+1 paper. **It is NOT evidence the universe runs on Collatz.**

## 2. ★ THE "3" HAS TWO READINGS THAT COINCIDE ONLY BY ACCIDENT (bank — Collatz-side observation)

- `E[2^−v] = Σ_v 2^−v·2^−v = Σ4^−v = 1/3` — the **annealed multiplier expectation / breakdown coefficient** (cascade reading).
- `Σ_v p_v² = Σ_v (2^−v)² = Σ4^−v = 1/3` — R8's **participation ratio** ⟹ `D₂ = log3/log q`, the **correlation dimension**.

These are the same number **only because `p_v = 2^−v`**: the geometric weights *are* the values being averaged. That coincidence
is why the cascade reading and the multifractal reading both find a "3", and it is special to Geom(1/2) — it would not survive a
different halving statistic. Worth a line in the paper: it explains why R8's D₂ and R5's rate name the same constant.

## 3. THE BRIDGE IS FALSIFIED — two independent kills

### 3a. THE KILL: negative skewness (λ-independent, therefore decisive)
Jouault, Schmiegel & Greiner (arXiv:chao-dyn/9909033): only **positively skewed** weight distributions reproduce the observed
unconditional *and* conditional multiplier distributions; this is what **rules out log-Poisson and log-stable**. Collatz:
- `skew(W) = −0.3742` · `skew(ln W) = −2.1213` (= −skew(Geom(1/2))). **Negative in both conventions.**
- **λ never touches W's law** (it enters ζ(p) only as a log base) ⇒ **this kill is λ-independent and cannot be tuned away.**
- Not bad luck: Collatz **inherits its log-Poisson sibling's known failure mode** (§1).
⚠️ *Caveat: JSG full text unretrievable (PDF extraction failed on every host); could not verify whether "positively skewed" is
stated on W or ln W. Conclusion is robust either way — Collatz is negative in both.*

### 3b. THE λ TENSION (option 2 — Collatz-specific, NOT degenerate; this is the interesting one)
Asymptotics: as p→∞ the `v=1` atom dominates (`E[W^s] → (1/2)(3/2)^s`), giving
`ζ(p) → (p/3)(1 − log_λ(3/2)) + log_λ2`, so **h_min = (1/3)(1 − log_λ(3/2))** and **C_∞ = log_λ2 = codimension of the most
intense structures**. Analytic and numeric agree exactly (log-space evaluation; the naive direct form overflows at large p).

Fit to **our banked ESS** ζ(p)/ζ(3) at R_λ=1280 (ESS is Re-invariant and cancels the fit-window bias in raw ζ(3)=1.0421):

| λ | max\|dev\| on ζ(2..6) | C_∞ |
|---|---|---|
| 2 — the 2-adic choice | 0.1082 | 1.000 (sheets) |
| 3 — R8's q-adic contraction | 0.0256 | 0.631 |
| **2.835 — best fit** | **0.0132** | **0.665** |
| √2 — forced to C_∞=2 | 0.4708 | 2.000 (filaments) |
| She–Levêque / measured | — | **2.000** |

**⇒ No single λ satisfies both.** The λ fitting ζ(2..6) to 1.3% implies codim ≈ 0.665 (near space-filling); turbulence measures
codim 2 (vortex filaments — the physical input She–Levêque gets right). Forcing codim 2 blows the exponents up 36×.
**Collatz cannot be simultaneously quantitatively and geometrically right.** This tension is Collatz-specific and survives the
degeneracy objection of §3c — it is the one non-vacuous structural finding here.
⚠️ *Medium confidence: C_∞=2 rests partly on She–Levêque's ansatz, and Iyer/Sreenivasan/Yeung (PRF 5, 054605, 2020) find
transverse ζ(p) SATURATES at ≈2.1 for p≥10 — fatal to Collatz's linear growth but equally fatal to She–Levêque and every log-ID
cascade. **Do not claim high-p saturation as a Collatz-specific kill.***

### 3c. AND THE TEST WAS NEARLY VACUOUS ANYWAY — ζ(p) cannot discriminate
**Eggers & Greiner, Nucl. Phys. Proc. Suppl. 92, 179 (2001)**, verbatim: cumulant ratios *"can successfully distinguish between
splitting functions **while multifractal scaling exponents and multiplier distributions cannot**."* That rules out **both**
discriminators originally contemplated. Corroborated by: Jouault/Schmiegel/Greiner (discrete and continuous cascades have
*"indistinguishable multiplier statistics"*); the **log-stable precedent** (fits ζ(p), ruled out by multiplier systematics —
our exact scenario, already run); **Molchan, Phys. Fluids 9, 2387 (1997)** (even K62 lognormal agrees over p ∈ (1,18)).
**Any one-parameter log-ID family fits ζ(p) over the measurable window — and Collatz's W is log-ID (§1), i.e. it sits in the
family where degeneracy is MAXIMAL.** The 0.013 fit at λ=2.835 measures the degeneracy, not the universe.

**Sharpest available discriminator if ever revisited:** cumulant ratios of ln ε (Eggers & Greiner 2001); second, the two-point
correlation of ε (Cleve, Greiner & Sreenivasan; arXiv:physics/0312113). **Not worth spending on** — §3a already kills it free.

## 4. THE PRIOR LOST TWICE (protocol note — both author priors were wrong, the conclusion survived anyway)
The handoff pre-registered FALSIFY, for two stated reasons. **Both reasons were wrong:**
1. *"E[W] ≠ 1, not energy-conserving"* → **E[W] = 1 exactly.** And the fact is circular (§1), so it is not evidence either way.
2. *"discreteness is falsifying"* → **it is not.** `τ(q) = log_b Σ p_i w_i^q` is real-analytic and strictly convex for ANY
   non-degenerate W, so f(α) is smooth for a 2-atom, countable, or continuous law alike (Halsey et al., PRA 33, 1141 (1986),
   built the formalism on a two-scale Cantor set). **She–Levêque's own log-Poisson multiplier is atomic.** The p-model
   (Meneveau & Sreenivasan, PRL 59, 1424 (1987)) has a two-valued non-ID generator and fits the whole exponent spectrum.
   Discreteness cannot be the falsifier.

The bridge dies on **skewness** (§3a) and the **λ tension** (§3b), plus vacuity (§3c) — none of which were predicted.
**Verdict correct, mechanism wrong.** Logged because this arc's protocol scores mechanisms, not just verdicts.

## 5. A METHOD NOTE ON THE NEAR-MISS (why λ=3 was seductive and still wrong)
R8 derives that the IFS `T_v(x) = (qx+1)/2^v` contracts by **exactly 1/q** — so λ=q=3 looked *derived*, not fitted, and at λ=3
the match to banked ESS is ≤0.026 zero-parameter. **The leap is illegitimate:** R8's 1/q contraction is **q-adic**
(non-archimedean — balls nested-or-disjoint, no transversality, collisions are exact algebraic coincidences). Turbulence's λ is
a **spatial, archimedean** scale ratio. Equating them is a **pun on "contraction," not a derivation.** The vacuity did not die at
λ=3; **it relocated into the metric.** Note also the best fit is λ*=2.835, not 3 — the "3-matches-3x+1" reading is not even what
the data picks.

## 6. LIT STATUS (subagent hunt, 2026-07-15)
- **Collatz ↔ turbulence: territory EMPTY.** ~10 query formulations across {Collatz, Syracuse, 3x+1, 2-adic} × {turbulence,
  cascade, multiplier, breakdown coefficient, intermittency, multifractal}. Zero hits — not real work, not crank work.
  ⚠️ **"Found no paper," NOT "verified empty"** — no MathSciNet/Zentralblatt/Scopus access. Confidence ~75%.
- **Nearest neighbours (both serious):** **V. I. Arnold (2004)**, "Number-theoretic turbulence…" — coined the *analogy*
  (statistics of Euler's function as "very intermittent"); an analogy, not a cascade claim. **A. Migdal**, "Decaying Turbulence
  and the Riemann Hypothesis" (arXiv:2604.12207) + Euler Ensemble — the one real occupant, **contested** (arXiv:2509.18992
  argues the momentum-loop measure cannot be nonnegative unless trivial).
  **Crucial asymmetry: Migdal DERIVES number theory OUT of Navier–Stokes. This proposal ran the opposite direction — importing a
  multiplier law from an unrelated arithmetic map and hoping it matched. That direction is exactly where the vacuity lives, and
  Migdal's work does not license it.**
- **Measured multiplier law:** continuous, bounded [0,1], symmetric about 1/2, Beta(β,β) (Chhabra & Sreenivasan, PRL 68, 2762
  (1992); Sreenivasan & Stolovitzky, J. Stat. Phys. 78, 311 (1995)). **But weaker than its reputation:** β unsettled (≈3
  atmospheric vs ≈8 in DNS, Hartlep/Cuzzi/Weston PRE 2017); **the symmetry P(m)=P(1−m) is an ESTIMATOR ARTIFACT** (m is a ratio
  to its container ⟹ m₁+m₂=1 identically) — *so the measured symmetry does NOT falsify Collatz's asymmetric multiplier; the data
  cannot see it*; and **cross-level independence is FALSE by S&S's own admission**, so every ζ(p)-from-multipliers formula
  (including this probe's) is known wrong in detail.

## 7. DISPOSITION
**Bridge: DEAD.** Killed λ-independently by negative skewness (§3a); structurally by the λ/codimension tension (§3b); and the
test that was supposed to adjudicate it is degenerate anyway (§3c). No further spend. **Do not revive without a positively
skewed variant AND a non-ζ(p) discriminator — and note that a "variant" chosen to fix skewness would be fitted, i.e. vacuous.**
**Kept:** §1 (log-geometric compound-Poisson generator, Kahane–Peyrière-valid, She–Levêque sibling) and §2 (the two readings of
the "3") — both are Collatz-side facts, independent of turbulence, usable in the qx+1 paper's Bernoulli-convolution framing.
**Not at stake:** R5 rate, R7 object identification, R8 mechanism, pillars 1–3, THEOREM_C_745. This probe touched none of them.
