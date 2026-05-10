# §3.5 Conjecture Statement — Link-Function Translation Clarification

## Section 1 — Current §3.5 text (verbatim, per brief quotation)

> "The mod-8 random effect u_{n mod 8} in the Bonacorsi-Bordoni NB2-GLM
> is, up to a global additive constant, log(a_final(n mod 8)) where
> a_final is the terminal value of the deterministic Collatz prefix at
> modular resolution k=3."

## Section 2 — Identified ambiguity (the unit-conversion error)

The conjecture as written equates two quantities living on different scales:

- **u_r** is a random effect inside an NB2-GLM with **log link**, i.e.,
  log(μ_i) = α + β·log(n_i) + u_r. So u_r is an additive shift on the
  **log-mean** of σ.
- **log(a_final(r))** is the closed-form prefix-decomposition quantity
  whose role in §4 is as a per-class additive offset on the **σ-additive**
  scale (via α_det(r) = prefix_steps(r) + K_h·log(a_final(r)/2^k), which
  is the per-class intercept of σ vs log(n) on the σ-additive scale).

Setting u_r equal to log(a_final(r)) directly is dimensionally inconsistent.
A reader who tests the literal claim by extracting B1's u_r posterior means
and comparing them to log(a_final(r)) will find a ~13× magnitude
discrepancy that is purely a unit-conversion artifact. The HMC validation
(B2 coefficient check, 2026-05-09) confirmed that once the link-function
translation is performed correctly, the empirical GLM-scale slope on
log(a_final) (0.0774, 95% CI [0.0676, 0.0869]) is consistent with the
predicted K_h/E[σ] = 0.0675 within linearization-approximation noise.

The §3.5 statement must therefore make the link-function translation
explicit so a reader cannot make the same unit-conversion error.

## Section 3 — Three rewrite options with commentary

### Option A — function-form first, scaling second

> "The mod-8 random effect u_{n mod 8} in the Bonacorsi-Bordoni NB2-GLM
> is monotone-increasing in log(a_final(n mod 8)) where a_final is the
> terminal value of the deterministic Collatz prefix at modular resolution
> k=3. Up to a global additive constant, the per-residue contribution to
> σ on the σ-additive scale is log(a_final(r)); the corresponding
> contribution on the GLM log-mean scale is approximately
> log(a_final(r))/E[σ]."

**Non-Bayesian-stats reader**: lead claim ("monotone-increasing in
log(a_final)") is comprehensible without inference machinery. The two-scale
follow-up may parse as background hedging rather than central content;
"GLM log-mean scale" is opaque without context.

**Lagarias (careful empiricist)**: would notice the lead claim weakens
the conjecture from a literal-magnitude identification to a directional
statement. Given the HMC validation now supports the literal form once
translated, the directional softening is conservative beyond what the
data requires. Would respect the explicit two-scale follow-up but ask why
the lead is weaker than what the evidence supports.

**Most natural empirical reference at this point in the paper**: §4 (the
σ-additive scale connects to α_det predictions); §5 (the GLM scale)
appears later. Option A's two-scale statement keeps both connections live,
which is appropriate at the end of §3.

### Option B — literal claim with explicit scale

> "The mod-8 random effect u_{n mod 8} in the Bonacorsi-Bordoni NB2-GLM
> (log link) is, up to a global additive constant,
> log(a_final(n mod 8))/E[σ] at modular resolution k=3. Equivalently, the
> per-residue effect on σ on the σ-additive scale is
> log(a_final(n mod 8)) up to a constant."

**Non-Bayesian-stats reader**: concrete formula on first read. The
"equivalently" clause helps. The 1/E[σ] coefficient is unfamiliar
notation; a reader without exposure to GLM linearization will not
immediately see why E[σ] is the right denominator.

**Lagarias**: would recognize 1/E[σ] as the leading-order link-function
rescaling. Would object that E[σ] is a data-dependent quantity — it
shifts with the n-distribution of the dataset — and that hardcoding it
into the conjecture statement makes the conjecture itself
distribution-dependent. Would prefer the conjecture stated on the
distribution-invariant σ-additive scale, with the GLM rescaling as a
derivation rather than the primary form.

**Most natural empirical reference**: §5 (the GLM dummies directly test
the log(a_final)/E[σ] form). §4's α_det evidence appears via the
"equivalently" clause.

### Option C — structural statement separating scales

> "At modular resolution k=3 of the Collatz prefix, the per-residue
> effect on σ is log(a_final(n mod 8)) up to a global additive constant.
> The corresponding random effect u_{n mod 8} in the Bonacorsi-Bordoni
> NB2-GLM is therefore log(a_final(n mod 8)) rescaled by the GLM's link
> function, with leading-order coefficient 1/E[σ] ≈ 1/154.4 at the scale
> of the data tested."

**Non-Bayesian-stats reader**: leads with a σ-scale claim ("per-residue
effect on σ is log(a_final)") that is mathematically clean and reads as
a structural identity. The link-function rescaling is then a derived
consequence; the reader who doesn't follow the rescaling can still take
the σ-additive form away. The numerical 1/154.4 is concrete and helpful.

**Lagarias**: would prefer the σ-additive lead — it is the cleanest
mathematical content of the conjecture, distribution-invariant, and
directly testable via §4's per-class regressions of σ on log(n).
"Leading-order coefficient" is honest acknowledgment of the
linearization, and "at the scale of the data tested" appropriately
hedges the 1/E[σ] empirical fit. Would accept this as the intended
conjecture.

**Most natural empirical reference**: §4 first (σ-additive structural
identity), §5 second (the rescaling check). This matches the paper's
order: §3.5 sits between §3.4 (the bridge to Tao on σ-additive scale)
and §4 (more σ-additive empirical evidence), with §5's GLM check arriving
later as the validation.

## Section 4 — Recommended option

**Recommendation: Option C.**

Justification:

1. **Structural primacy on σ-additive scale**: log(a_final(r)) is a
   structural quantity from the deterministic Collatz prefix. Stating
   the conjecture on σ-additive scale puts the structural identity at
   the center, where it is distribution-invariant and directly testable
   via §4. The GLM rescaling becomes a derived consequence rather than
   the primary form.

2. **Honest about linearization**: "leading-order coefficient" makes
   transparent that 1/E[σ] is the first-order Jacobian of the log link
   evaluated at the class-mean σ. This is the precise content that the
   B2 coefficient check (2026-05-09) verified to ~5% under a
   prefix_steps refinement. A reader cannot make the unit-conversion
   error because the rescaling is named.

3. **Distribution-invariant primary claim**: the σ-additive form does
   not reference E[σ]; only the GLM-rescaling derivative does.
   Lagarias-style readers will prefer this because the conjecture's
   structural content does not change with the n-sampling distribution.

4. **Scope-narrowing language**: "at the scale of the data tested" is
   exactly the kind of empirical-fit hedge appropriate for a conjecture
   whose validation is at N=10⁴ HMC, with full-N HMC pending at the
   Bonacorsi-Bordoni group's production scale.

5. **Position in the paper**: §3.5 sits between §3.4's σ-additive bridge
   to Tao (5.15) and §4's σ-additive empirical content, so the σ-additive
   lead is the most natural continuation. The GLM rescaling forward-points
   to §5's empirical check.

Option A's monotone-increasing lead is conservative beyond what the
data supports. Option B's lead form embeds E[σ] in the primary
conjecture statement, making the conjecture distribution-dependent.
Option C avoids both.

## Section 5 — Recommended §3.5 LaTeX block (ready to drop in)

```latex
\subsection{The conjecture}
\citet{bonacorsi2026} introduce a $\bmod 8$ residue-class random effect
$u_{n \bmod 8}$ in their NB2-GLM. We conjecture:

\begin{conjecture}
At modular resolution $k=3$ of the Collatz prefix, the per-residue
effect on the total stopping time $\sigma$ is $\log(a_{\text{final}}(n \bmod 8))$
up to a global additive constant. The corresponding random effect
$u_{n \bmod 8}$ in the Bonacorsi-Bordoni NB2-GLM with log link is
therefore $\log(a_{\text{final}}(n \bmod 8))$ rescaled by the GLM's
link function, with leading-order coefficient $1 / \mathbb{E}[\sigma]$
on the log-mean scale, which equals approximately $1/154.4$ at the
scale of the data tested in \S5.
\end{conjecture}
```

---

*Inputs: brief-quoted §3.5 text and HMC validation results from
docs/paper1_b2_coefficient_check.md and
docs/paper1_prefix_steps_confound_check.md (2026-05-09). No edits to
main.tex; this document is a recommendation for editorial review.*
