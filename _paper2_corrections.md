**Paper 2 corrections --- single apply-in-order checklist**

A Symbolic Prefix Decomposition of the Collatz Map. Cross-referenced
against closed_form_findings.md (Result 1, lines 105--141) and the repo
experiment index.

**The structural shift this document makes:** Your headline result
⟨α_det⟩ = log(6)/log(4/3) is currently presented as a heuristic identity
verified empirically. It is actually a theorem provable for all k ≥ 1,
with one named heuristic input (the K_h random-walk constant). The proof
reduces to Terras 1976 Lemma 1 (the parity-vector bijection). After this
change, the prefix-algebra portion of the paper becomes fully
theorem-grade. Empirical Observation 1 collapses to a one-line
corollary. Two of the three \"open questions\" in §6 and §7 disappear.
The empirical content of §4 remains exactly what it is --- validation of
the K_h embedding\'s predictive power on per-class σ statistics --- but
it\'s no longer carrying the closed-form result on its back. This is the
rare review where the fix makes the paper substantially stronger, not
weaker.

How to use this document

Apply the steps in order, top to bottom. Each step is self-contained: it
tells you where to look, what\'s wrong, the exact current text, the
exact replacement text, and why. Later steps assume earlier ones are
applied --- for example, Step 6 references Lemma\~\\ref{lem:j-binomial}
introduced in Step 2.

There are 18 steps. Time estimate: a long afternoon. The first six are
the heavy ones (the proof and its immediate downstream renames). Steps
7--14 are the textual cleanup that ripples through main.tex once the
proof is in. Steps 15--18 are non-proof-related fixes.

Line numbers reference the main.tex you sent (32,734 bytes, May 17
2026). After you apply Step 2, line numbers in later steps will shift
--- search by quoted text instead of line numbers if that becomes a
problem.

Priority levels:

- CRITICAL --- referee will spot this in minutes.

- MAJOR --- real defect, will draw serious referee comment.

- DOWNSTREAM --- a mechanical rename or rephrase forced by an earlier
  step. Fast to apply, but if you skip you\'ll have internal
  contradictions.

- MINOR --- line-level cleanup.

- NOTE --- style/structure suggestion.

Summary --- all 18 steps in apply order

  -------------------------------------------------------------------------------------------
  **\#**   **Priority**     **Where**        **What changes**
  -------- ---------------- ---------------- ------------------------------------------------
  **1**    **MAJOR**        §3.2 Theorem 1   Add missing sentence: odd-steps cannot occur
                            proof (lines     consecutively because 3c+1 is even when c is
                            86--94)          odd.

  **2**    **CRITICAL**     §3.5 (lines      Replace Heuristic Identity 1 + Derivation with:
                            113--125) ---    Lemma (block decomposition), Lemma (Syracuse
                            replace whole    equivalence), Lemma (binomial j-distribution,
                            subsection       proven via Terras 1976 Lemma 1), Theorem (closed
                                             form ⟨α_det⟩ = log 6/log(4/3) under K_h
                                             heuristic only).

  **3**    **DOWNSTREAM**   §3.3 (lines      Empirical Observation 1 → Corollary of
                            100--106)        Lemma\~\\ref{lem:j-binomial}. One-line proof:
                                             C(k−1, j−1) ≥ 1 for j ∈ {1,...,k}.

  **4**    **DOWNSTREAM**   Preamble (lines  Remove the \\newtheorem for {heuristic} and
                            17, 19)          {observation} for §3.5/§3.3 results; keep
                                             {observation} only for §4 empirical claims. (Or
                                             repurpose, see step.)

  **5**    **DOWNSTREAM**   §1 abstract      Update self-description: ⟨α_det⟩ is a theorem
                            paragraph (line  under K_h heuristic, not a heuristic identity.
                            35)              Realization is a corollary, not an empirical
                                             observation.

  **6**    **DOWNSTREAM**   §1 organization  Update the §3 outline to match new structure:
                            (line 43)        theorem instead of heuristic identity, corollary
                                             instead of empirical observation.

  **7**    **DOWNSTREAM**   §1 (line 37) --- Reframe finding (iii): the closed form is now a
                            empirical        theorem; §4.3 is numerical verification of a
                            findings (iii)   theorem, not empirical verification of a
                                             heuristic.

  **8**    **DOWNSTREAM**   §3.4 ahead of    Update the paragraph introducing α_det: remove
                            the new §3.5     the framing \"prediction holds up to a global
                            (line 109--111)  additive constant\" since the closed form is now
                                             exact, and tighten the \"derivation is
                                             heuristic\" paragraph.

  **9**    **DOWNSTREAM**   §4.3 (lines      Observation 3 reframed: numerical verification
                            184--188)        of Theorem\~\\ref{thm:alpha-det-closed-form},
                                             not empirical verification of a heuristic.

  **10**   **DOWNSTREAM**   §5 (lines        Substantially trim. Frame the prefix
                            194--203)        decomposition as actually using Terras 1976
                                             Lemma 1 (parity-vector bijection) to give a
                                             constructive sharpening of Terras 1976 Lemma 4.

  **11**   **DOWNSTREAM**   §6.1(a) (line    Remove this open question. It is now a theorem
                            213)             (Lemma\~\\ref{lem:j-binomial}, proven). Renumber
                                             §6.1(b) and §6.1(c) to (a) and (b).

  **12**   **DOWNSTREAM**   §7 paragraph 1   Reword: closed form is a theorem, not a
                            (line 228)       heuristic identity. The only heuristic input is
                                             K_h.

  **13**   **DOWNSTREAM**   §7 paragraph 2   Reword \"empirically verified predictive power\"
                            (line 230)       to be clearer that §4 validates the K_h
                                             embedding, not the closed form itself.

  **14**   **DOWNSTREAM**   §7 paragraph 3   Two open questions remain (post-prefix
                            (line 232)       conditional distribution, and Tao 2022
                                             connection), not three. Remove the realization
                                             item.

  **15**   **MAJOR**        Appendix A (line Wrong scripts cited. Replace
                            258)             alpha_beta_gamma_decay.py +
                                             analytical_abc_derivation.py with
                                             experiments/01_alpha_decomposition.py +
                                             experiments/24_k_sweep_alpha_decomposition.py.

  **16**   **MAJOR**        §4.1 (line 155)  Promote the matched-n controlled comparison out
                                             of parentheses into its own paragraph + small
                                             table.

  **17**   **MAJOR**        §4.2 (lines      Add a figure. The k-distinct-clusters result
                            158--182)        needs to be shown, not just stated.

  **18**   **MAJOR**        references.bib   Bibliography too thin. Add Sinai 2003, Chang
                                             2026 (both arXiv IDs), Bonacorsi-Bordoni 2026.
                                             Cite humphrey_collatz_2026 in §1 or remove it.
  -------------------------------------------------------------------------------------------

Minor and Note items (style cleanup) follow after Step 18.

Dependency map

Which steps depend on which. Lines indicate \"this earlier step must be
done first\":

Step 1 (Theorem 1 proof fix) ─── independent

Step 2 (§3.5 new proof) ───┬─► Step 3 (Obs. 1 → Corollary)

├─► Step 4 (remove \\newtheorem)

├─► Step 5 (§1 abstract reframe)

├─► Step 6 (§1 organization reframe)

├─► Step 7 (§1 finding (iii) reframe)

├─► Step 8 (§3.4 intro reframe)

├─► Step 9 (§4.3 reframe)

├─► Step 10 (§5 trim + new Terras framing)

├─► Step 11 (§6.1 remove open question)

├─► Step 12 (§7 ¶1 reword)

├─► Step 13 (§7 ¶2 reword)

└─► Step 14 (§7 ¶3 update open questions)

Steps 15--18 ─── independent (Appendix A, §4.1, §4.2, references.bib)

Practical reading of this map: if you apply Step 2 cleanly, Steps 3--14
are mostly mechanical find-and-replace. Don\'t skip them --- every one
is a place where the paper would otherwise contradict itself internally.

The 18 steps

**STEP 1 MAJOR · §3.2, proof of Theorem 1 --- last paragraph (lines
89--91)**

**What\'s wrong:** The j ≤ k bound is asserted via \"the number of
odd-steps interleaved among them is at most k\", but this doesn\'t
actually rule out arbitrarily many odd-steps before the first even-step
(case (ii) doesn\'t change v_2(a)). The real reason j ≤ k is that after
any case (ii) step, c\_{t+1} = 3c_t + 1 is even, so the next step must
be case (i). Odd-steps cannot occur consecutively. That sentence is
missing.

**Current text in main.tex:**

The lower bound \$j \\geq 1\$ follows from the initial state: \$a_0 =
2\^k\$ is even and \$c_0 = r\$ is odd (since \$r\$ is odd by
hypothesis), so the first step of the prefix iteration is case (ii), an
odd-step application. Therefore the prefix consumes at least one
odd-step application, giving \$j \\geq 1\$. The upper bound \$j \\leq
k\$ follows from the 2-adic valuation argument above: each even-step
decrements \$v_2(a)\$ by exactly \$1\$, the prefix terminates when \$v_2
= 0\$, so the prefix length is bounded above by \$k\$ even-steps and the
number of odd-steps interleaved among them is at most \$k\$.

**Replace with:**

The lower bound \$j \\geq 1\$ follows from the initial state: \$a_0 =
2\^k\$ is even and \$c_0 = r\$ is odd, so the first step is case (ii),
giving \$j \\geq 1\$.

For the upper bound \$j \\leq k\$: after any case-(ii) step, \$c\_{t+1}
= 3 c_t + 1\$ is even (since \$c_t\$ is odd, \$3 c_t\$ is odd, and \$3
c_t + 1\$ is even), so the next step is necessarily case (i). Odd-steps
cannot occur consecutively. The prefix consumes exactly \$k\$ case-(i)
steps before termination (by the \$2\$-adic valuation argument). Each
odd-step is immediately followed by an even-step, so the number of
odd-steps is at most the number of even-steps, giving \$j \\leq k\$.

**Why:** *One sentence patch. Without it, a referee reading the j ≤ k
argument stops. Important: this sentence is also a load-bearing input to
Step 2\'s block decomposition lemma --- the proof of the binomial
j-distribution uses this exact parity argument. Do this step first so
the §3.5 proof has its prerequisite.*

**STEP 2 CRITICAL · §3.5, replace lines 113--125 entirely (from
\\begin{heuristic} through \\end{derivation} and the discussion
paragraph that follows)**

**What\'s wrong:** Heuristic Identity 1 currently asserts ⟨α_det⟩ = log
6/log(4/3) using an independent-Bernoulli heuristic that gives ⟨j⟩ = 2
and produces a k-dependent ⟨α_det⟩. With ⟨j⟩ = 2: ⟨α_det⟩ = k(1 − K_h
log 2) + 2(1 + K_h log 3), coefficient of k is (1 − K_h log 2) ≈ −6.23,
not zero. So the heuristic as stated doesn\'t give the result it claims
to. The actual derivation in closed_form_findings.md Result 1 uses a
binomial j-distribution: the count of odd r mod 2\^k with j odd-steps is
exactly C(k−1, j−1), giving ⟨j⟩ = (k+1)/2 and exact cross-cancellation.
This binomial j-distribution is a theorem provable via Terras 1976 Lemma
1 --- not a heuristic at all. This step replaces the §3.5 content with a
full proof.

**Current text in main.tex:**

\\begin{heuristic}\[Closed-form \$\\langle \\alpha\_{\\text{det}}
\\rangle\$\]

\\label{heur:alpha_det}

Under the heuristic that the prefix terminates at \$a\_{\\text{final}} =
3\^j\$ with probability \$2\^{(-j)}\$ (modeling each odd-step decision
as an independent Bernoulli event), the class-averaged
\$\\alpha\_{\\text{det}}\$ satisfies \$\\langle \\alpha\_{\\text{det}}
\\rangle = \\log(6)/ \\log(4/3)\$.

\\end{heuristic}

\\begin{derivation}

\\label{der:bern}

Under the independent-Bernoulli heuristic, the expected number of
odd-step applications during the prefix is \$\\mathbb{E}\[j\] =
\\sum\_{j \\geq 1} j \\cdot 2\^{(-j)}=2\$ (truncated by the bound \$j
\\leq k\$, but the truncation is negligible at \$k \\geq 6\$).
Substituting into \$\\alpha\_{\\text{det}}\$ and computing the class
average across the heuristic distribution gives \$\\langle
\\alpha\_{\\text{det}} \\rangle = \\langle \\text{prefix\\\_steps}
\\rangle + K_h \\cdot \\langle \\log(a\_{\\text{final}} / 2\^k)
\\rangle\$.

The two terms combine via the random-walk heuristic to yield \$\\log(6)/
\\log(4/3)\$. The full derivation is given in the accompanying
repository \\verb\|(closed_form_findings.md)\|.

\\end{derivation}

We name the result a heuristic identity rather than a theorem because
the underlying combinatorial assumption \-- that each odd-step decision
is an independent Bernoulli event \-- is a heuristic. The actual
distribution of \$a\_{\\text{final}}\$ across residue classes is
determined by the prefix iteration\'s parity dynamics, not by
independent coin flips. The heuristic is an approximation. The empirical
verification in \\S4.3 confirms that \$\\langle \\alpha\_{\\text{det}}
\\rangle\$ matches \$\\log(6)/ \\log(4/3) \\approx 6.226\$ to within
sampling precision at modular resolutions \$k \\in \\{6, 7, 8, 9\\}\$.

**Replace with:**

\\begin{lemma}\[Block decomposition\]

\\label{lem:block-decomposition}

The prefix trace of any odd \$r \\bmod 2\^k\$ decomposes uniquely as a
sequence of \$k\$ blocks, each of one of two types: Type O is a
case-(ii) step immediately followed by a case-(i) step (the substring
\`oe\'); Type E is a single case-(i) step (the substring \`e\'). The
first block of every prefix trace is of type O, and the number of type-O
blocks equals \$j\$, the number of odd-step applications in the prefix.

\\end{lemma}

\\begin{proof}

By the argument in the proof of Theorem\~\\ref{thm:termination}, after
any case-(ii) step \$c\_{t+1} = 3 c_t + 1\$ is even, so the next step is
case (i): odd-steps cannot occur consecutively. Grouping each odd-step
with its successor partitions the trace into type-O blocks (\`oe\') and
type-E blocks (\`e\'). The first block is type O because \$c_0 = r\$ is
odd. By Theorem\~\\ref{thm:termination} the prefix consumes exactly
\$k\$ case-(i) steps in total; each block contains one case-(i) step, so
the total number of blocks is \$k\$. The number of type-O blocks equals
the number of odd-steps, \$j\$.

\\end{proof}

\\begin{lemma}\[Block sequence equals Syracuse parity vector\]

\\label{lem:syracuse-equivalence}

Let \$S \\colon \\mathbb{N} \\to \\mathbb{N}\$ denote the Syracuse
compressed map \$S(c) = (3c+1)/2\$ if \$c\$ is odd and \$S(c) = c/2\$ if
\$c\$ is even. Define \$c_0 := r\$ and \$c\_{t+1} := S(c_t)\$ for \$t
\\geq 0\$. Each block of the prefix trace at level \$k\$ corresponds to
one application of \$S\$: the \$t\$-th block is type O if \$c_t\$ is odd
and type E if \$c_t\$ is even. The block sequence of the prefix is the
length-\$k\$ Syracuse parity vector \$\\Pi_k(r) := (c_0 \\bmod 2, c_1
\\bmod 2, \\ldots, c\_{k-1} \\bmod 2)\$.

\\end{lemma}

\\begin{proof}

Type O at step \$t\$: \$c_t\$ is odd, so the odd-step gives an
intermediate value \$3 c_t + 1\$ (even), followed by the forced
even-step yielding \$c\_{t+1} = (3 c_t + 1)/2 = S(c_t)\$. Type E at step
\$t\$: \$c_t\$ is even, so the single even-step gives \$c\_{t+1} = c_t /
2 = S(c_t)\$. In both cases the block applies \$S\$ once to \$c\$ and
records \$c_t \\bmod 2\$.

\\end{proof}

\\begin{lemma}\[Binomial \$j\$-distribution\]

\\label{lem:j-binomial}

Among the \$2\^{k-1}\$ odd residue classes \$r \\bmod 2\^k\$, the number
whose prefix iteration terminates with exactly \$j\$ odd-step
applications is \$\\binom{k-1}{j-1}\$, for \$j \\in \\{1, 2, \\ldots,
k\\}\$. Consequently \$\\langle j \\rangle = (k+1)/2\$.

\\end{lemma}

\\begin{proof}

By Lemma\~\\ref{lem:syracuse-equivalence}, the block sequence of the
prefix at level \$k\$ equals the Syracuse parity vector \$\\Pi_k(r)\$,
and the number of type-O blocks equals the Hamming weight of
\$\\Pi_k(r)\$. The classical Syracuse parity-vector bijection
(\\citet{terras1976}, Lemma 1; \\citet{lagarias1985}, \\S 2) states that
\$\\Pi_k \\colon \\mathbb{Z}/2\^k\\mathbb{Z} \\to \\{0,1\\}\^k\$ is a
bijection. Restricting to odd \$r\$ corresponds to fixing the first
coordinate \$c_0 \\bmod 2 = 1\$, giving a bijection between odd residues
mod \$2\^k\$ and \$\\{1\\} \\times \\{0,1\\}\^{k-1}\$. The number of
vectors in this restricted set with Hamming weight \$j\$ (i.e., with
\$j-1\$ ones among the \$k-1\$ free coordinates) is
\$\\binom{k-1}{j-1}\$.

The mean is \$\\langle j \\rangle = 2\^{-(k-1)} \\sum\_{j=1}\^{k} j
\\binom{k-1}{j-1} = 2\^{-(k-1)} \\sum\_{i=0}\^{k-1} (i+1)
\\binom{k-1}{i} = 2\^{-(k-1)} \\left\[ (k-1) 2\^{k-2} + 2\^{k-1}
\\right\] = (k+1)/2\$.

\\end{proof}

\\begin{theorem}\[Closed form for \$\\langle \\alpha\_{\\text{det}}
\\rangle\$\]

\\label{thm:alpha-det-closed-form}

Under the random-walk heuristic (\$K_h = 3/\\log(4/3)\$) applied to the
post-prefix state and Lemma\~\\ref{lem:j-binomial},

\$\$\\langle \\alpha\_{\\text{det}} \\rangle = \\frac{\\log
6}{\\log(4/3)}, \\quad \\text{independent of } k.\$\$

\\end{theorem}

\\begin{proof}

From \$\\alpha\_{\\text{det}}(r) = \\text{prefix\\\_steps}(r) + K_h
\\log(a\_{\\text{final}}(r)/2\^k)\$, Theorem\~\\ref{thm:termination}\'s
\$\\text{prefix\\\_steps}(r) = k + j(r)\$ and \$a\_{\\text{final}}(r) =
3\^{j(r)}\$:

\$\$\\alpha\_{\\text{det}}(r) = (k + j) + K_h (j \\log 3 - k \\log
2).\$\$

Averaging over odd \$r \\bmod 2\^k\$ and using \$\\langle j \\rangle =
(k+1)/2\$ from Lemma\~\\ref{lem:j-binomial}:

\$\$\\langle \\alpha\_{\\text{det}} \\rangle = k(1 - K_h \\log 2) +
\\frac{k+1}{2}(1 + K_h \\log 3).\$\$

Writing \$L = \\log(4/3) = 2 \\log 2 - \\log 3\$ and using \$K_h =
3/L\$:

\$\$1 - K_h \\log 2 = \\frac{L - 3 \\log 2}{L} = -\\frac{\\log 6}{L},
\\qquad 1 + K_h \\log 3 = \\frac{L + 3 \\log 3}{L} = \\frac{2 \\log
6}{L},\$\$

so

\$\$\\langle \\alpha\_{\\text{det}} \\rangle = -\\frac{k \\log 6}{L} +
\\frac{k+1}{2} \\cdot \\frac{2 \\log 6}{L} = \\frac{\\log
6}{L}\\left\[-k + (k+1)\\right\] = \\frac{\\log 6}{\\log(4/3)}.\$\$

\\end{proof}

Theorem\~\\ref{thm:alpha-det-closed-form} is exact in the prefix algebra
and Lemma\~\\ref{lem:j-binomial}; the only non-deductive input is the
application of \$K_h\$ from the random-walk argument to the post-prefix
state. Direct enumeration of \$\\alpha\_{\\text{det}}(r)\$ across all
\$2\^{k-1}\$ odd classes at \$k \\in \\{6, 8, 10, 12, 14\\}\$ recovers
\$\\log 6 / \\log(4/3) \\approx 6.228263\$ to machine precision
(relative error \$\\leq 10\^{-14}\$), as verified in
\\S\\ref{sec:closed-form-verification}.

**Why:** *This is the structural shift. Three lemmas + one theorem
replace one heuristic identity + one shaky derivation. The proof reduces
to Terras 1976 Lemma 1, so the paper\'s framing as a Terras sharpening
becomes literal --- you use Terras\'s parity-vector bijection (Lemma 1)
to give a constructive version of Terras\'s iterate-ratio identity
(Lemma 4). The dependency on Terras is structural, not just narrative.
Verified computationally through k=16; the proof works for all k.*

**STEP 3 DOWNSTREAM · §3.3, lines 100--106 (Empirical Observation 1)**

**What\'s wrong:** After Step 2 introduces Lemma\~\\ref{lem:j-binomial},
the realization of every j ∈ {1,...,k} is a one-line corollary (C(k−1,
j−1) ≥ 1). Keeping it as an Empirical Observation is wrong --- it\'s
deductive.

**Current text in main.tex:**

\\begin{observation}\[Realization\]

\\label{obs:realization}

For modular resolutions \$k \\in \\{4, 5, 6, 7, 8, 9, 10, 11, 12\\}\$,
every \$j \\in \\{1, \\ldots, k\\}\$ is realized as \$a\_{\\text{final}}
= 3\^j\$ by some odd residue class \$r \\bmod 2\^k\$.

The observation has been verified by direct computation of the prefix
iteration for all \$2\^{(k-1)}\$ odd residue classes at each cited
\$k\$, with the realized set of \$a\_{\\text{final}}\$ values enumerated
in each case. We name the observation as empirical rather than as a
theorem because we do not have a proof that all \$k\$ values are
realized at every \$k\$ for arbitrary \$k\$. The observation is
consistent with the available evidence and the structure of the prefix
iteration suggests (but does not prove) that the realization extends to
arbitrary \$k\$. We present the result as what we have verified
directly, without extrapolation.

\\end{observation}

**Replace with:**

\\begin{corollary}\[Realization\]

\\label{cor:realization}

For every \$k \\geq 1\$ and every \$j \\in \\{1, \\ldots, k\\}\$, there
exists an odd residue class \$r \\bmod 2\^k\$ with
\$a\_{\\text{final}}(r) = 3\^j\$.

\\end{corollary}

\\begin{proof}

By Lemma\~\\ref{lem:j-binomial}, the number of such classes is
\$\\binom{k-1}{j-1} \\geq 1\$.

\\end{proof}

**Why:** *Mechanical. The corollary is shorter and stronger (all k, not
just k ≤ 12) than the observation it replaces.*

**STEP 4 DOWNSTREAM · Preamble, lines 17 and 19**

**What\'s wrong:** After Steps 2 and 3, the {heuristic} and {derivation}
environments are no longer used for §3.5, and {observation} is no longer
used for §3.3. The {observation} environment is still used in §4 for
genuine empirical claims. Keep {observation}; remove or repurpose
{heuristic} and {derivation}.

**Current text in main.tex:**

\\newtheorem{observation}\[theorem\]{Empirical Observation}

\\theoremstyle{plain}

\\newtheorem{heuristic}\[theorem\]{Heuristic Identity}

\\newtheorem{derivation}\[theorem\]{Derivation Sketch}

**Replace with:**

\\newtheorem{observation}\[theorem\]{Empirical Observation}

\% {heuristic} and {derivation} environments are no longer needed:

\% \\S3.5\'s results are now Theorem\~\\ref{thm:alpha-det-closed-form}
(proper theorem)

\% and \\S3.3\'s realization is now Corollary\~\\ref{cor:realization}.

\% \\newtheorem{heuristic}\[theorem\]{Heuristic Identity} % removed

\% \\newtheorem{derivation}\[theorem\]{Derivation Sketch} % removed

**Why:** *Cleanup. Unused theorem environments are harmless but
distracting. Once Steps 2 and 3 are in, these have no remaining
references in the body.*

**STEP 5 DOWNSTREAM · §1, second paragraph (line 35) --- the
contribution self-description**

**What\'s wrong:** Currently says \"A closed-form quantity ⟨α_det⟩ = log
6/log(4/3) is derived from the prefix iteration under a heuristic
combinatorial argument; we name this a heuristic identity rather than a
theorem because the underlying combinatorial assumption is itself
heuristic.\" After Step 2, this is now a theorem under K_h only.

**Current text in main.tex:**

This paper establishes properties of the prefix decomposition and
presents empirical evidence that the decomposition\'s terminal value
\$a\_{\\text{final}}\$ determines per-class statistics of the total
stopping time \$\\sigma\$ on odd integers in residue class \$r \\bmod
2\^k\$. The contribution has two parts. The algebraic part, presented in
\\S3, establishes the properties of the prefix iteration that are
deductive consequences of the symbolic state\'s parity dynamics. The
principal algebraic statement is the termination theorem (\\S3.2), which
establishes that the prefix iteration terminates within at most \$k\$
odd-step applications and produces a terminal \$a\_{\\text{final}}\$ of
the form \$3\^j\$ for some \$1 \\leq j \\leq k\$. A closed-form quantity
\$\\langle \\alpha\_{\\text{det}} \\rangle = \\log(6)/ \\log(4/3)\$ is
derived (\\S3.5) from the prefix iteration under a heuristic
combinatorial argument; we name this a heuristic identity rather than a
theorem because the underlying combinatorial assumption is itself
heuristic.

**Replace with:**

This paper establishes properties of the prefix decomposition and
presents empirical evidence that the decomposition\'s terminal value
\$a\_{\\text{final}}\$ determines per-class statistics of the total
stopping time \$\\sigma\$ on odd integers in residue class \$r \\bmod
2\^k\$. The contribution has two parts. The algebraic part, presented in
\\S3, establishes the properties of the prefix iteration that are
deductive consequences of the symbolic state\'s parity dynamics. The
principal algebraic statements are the termination theorem (\\S3.2; the
prefix iteration consumes exactly \$k\$ even-step applications and at
most \$k\$ odd-step applications, with terminal \$a\_{\\text{final}}
\\in \\{3\^j : 1 \\leq j \\leq k\\}\$), the binomial \$j\$-distribution
lemma (\\S3.5; the count of odd classes mod \$2\^k\$ with \$j\$ odd-step
applications is \$\\binom{k-1}{j-1}\$, proven via the classical Syracuse
parity-vector bijection), and the closed-form theorem \$\\langle
\\alpha\_{\\text{det}} \\rangle = \\log(6)/\\log(4/3)\$ (\\S3.5; exact
in the prefix algebra, with the random-walk constant \$K_h\$ as the only
heuristic input).

**Why:** *Self-description must match the new structure. This is the
first paragraph a referee reads; getting it wrong puts the whole paper
on the wrong footing.*

**STEP 6 DOWNSTREAM · §1, last paragraph (line 43) --- organization
outline**

**What\'s wrong:** The outline of §3 still describes the old structure:
\"empirical observation about realization\", \"heuristic closed form\".
Update to match.

**Current text in main.tex:**

The paper is organized as follows. Section 2 places the decomposition in
its lineage within the standard Collatz random-walk heuristic,
\\citet{terras1976} Lemma 4, and symbolic approaches to Collatz residue
dynamics. Section 3 defines the prefix decomposition, states and proves
the termination theorem, presents the empirical observation about
realization of all \$k\$ \$a\_{\\text{final}}\$ values, defines the
closed-form predictor \$\\alpha\_{\\text{det}}\$, and derives the
heuristic closed form \$\\langle \\alpha\_{\\text{det}} \\rangle =
\\log(6) / \\log(4/3)\$. Section 4 presents empirical regularities at
modular resolutions \$k \\in \\{6, 7, 8, 9\\}\$ and \$N\$ up to
\$2\^{27}\$. Section 5 places the decomposition relative to
\\citet{terras1976} explicitly. Section 6 discusses what the empirical
regularities suggest, what is required for a stronger theorem-grade
result, and the \$qx+1\$ extension. Section 7 concludes.

**Replace with:**

The paper is organized as follows. Section 2 places the decomposition in
its lineage within the standard Collatz random-walk heuristic,
\\citet{terras1976} Lemma 4, and symbolic approaches to Collatz residue
dynamics. Section 3 defines the prefix decomposition, states and proves
the termination theorem (\\S\\ref{sec:termination}), establishes
realization of all \$k\$ values of \$a\_{\\text{final}}\$ as a corollary
(\\S\\ref{sec:realization}), defines the closed-form predictor
\$\\alpha\_{\\text{det}}\$, and proves \$\\langle \\alpha\_{\\text{det}}
\\rangle = \\log(6)/\\log(4/3)\$ via the binomial \$j\$-distribution
lemma (\\S\\ref{sec:closed-form}). Section 4 presents empirical
regularities at modular resolutions \$k \\in \\{6, 7, 8, 9\\}\$ and
\$N\$ up to \$2\^{27}\$ that validate the random-walk-heuristic
embedding \$K_h\$ underlying \$\\alpha\_{\\text{det}}\$. Section 5
places the decomposition relative to \\citet{terras1976}. Section 6
discusses what the empirical regularities suggest, what remains open,
and the \$qx+1\$ extension. Section 7 concludes.

**Why:** *Pulls the organization paragraph into alignment with the new
§3 structure. Also reframes §4 honestly: §4 validates the K_h embedding,
not the closed form itself (the closed form is proven in §3.5). Note:
this assumes you add \\label{sec:termination}, \\label{sec:realization},
\\label{sec:closed-form} to the corresponding subsections in §3 --- do
that as part of this step.*

**STEP 7 DOWNSTREAM · §1, third paragraph (line 37) --- empirical
findings list, item (iii)**

**What\'s wrong:** Finding (iii) currently says \"the closed form
⟨α_det⟩ = log 6/log(4/3) is verified empirically at each tested k to
within sampling precision.\" After Step 2, the closed form is a theorem
and §4.3 is a numerical sanity check, not empirical verification.

**Current text in main.tex:**

The empirical part, presented in \\S4, presents regularities verified at
modular resolutions \$k \\in \\{6, 7, 8, 9\\}\$ and \$N\$ up to
\$2\^{27}\$. The principal empirical findings are: (i) per-class
intercept \$\\alpha(r)\$ of \$\\sigma\$ vs \$\\log(n)\$ is predicted by
\$\\alpha\_{\\text{det}}(r)\$ with regression \$\\text{R}\^2 \\geq
0.985\$ and signal-to-noise ratio \$\\text{SD(residual)} /
\\text{mean(SE)}\$ in \$\[0.91, 0.99\]\$ across the four tested
resolutions, with per-class residuals scaling with sampling noise and no
detectable structure above sampling noise; (ii) per-class central
moments of \$\\sigma\$ residuals cluster by \$a\_{\\text{final}}\$, with
the number of distinct clusters equal to \$k\$ at every tested \$k\$;
and (iii) the closed form \$\\langle \\alpha\_{\\text{det}} \\rangle =
\\log(6)/ \\log(4/3)\$ is verified empirically at each tested \$k\$ to
within sampling precision. The empirical findings are presented as such;
we do not claim the prefix decomposition\'s predictive content is
theorem-proven for arbitrary \$k\$ or \$N\$.

**Replace with:**

The empirical part, presented in \\S4, presents regularities verified at
modular resolutions \$k \\in \\{6, 7, 8, 9\\}\$ and \$N\$ up to
\$2\^{27}\$. The principal empirical findings are: (i) per-class
intercept \$\\alpha(r)\$ of \$\\sigma\$ vs \$\\log(n)\$ is predicted by
\$\\alpha\_{\\text{det}}(r)\$ with regression \$R\^2 \\geq 0.985\$ and
signal-to-noise ratio \$\\text{SD(residual)} / \\text{mean(SE)}\$ in
\$\[0.91, 0.99\]\$ across the four tested resolutions, with per-class
residuals scaling with sampling noise and no detectable structure above
sampling noise; and (ii) per-class central moments of \$\\sigma\$
residuals cluster by \$a\_{\\text{final}}\$, with the number of distinct
clusters equal to \$k\$ at every tested \$k\$. These findings validate
the random-walk embedding \$K_h\$ underlying \$\\alpha\_{\\text{det}}\$
on per-class \$\\sigma\$ statistics; we do not claim that the per-class
prediction is theorem-proven for arbitrary \$k\$ or \$N\$. The closed
form \$\\langle \\alpha\_{\\text{det}} \\rangle = \\log(6)/\\log(4/3)\$
itself is Theorem\~\\ref{thm:alpha-det-closed-form}, not an empirical
claim; \\S\\ref{sec:closed-form-verification} reports the corresponding
numerical sanity check.

**Why:** *Removes finding (iii) from the empirical list --- it doesn\'t
belong there anymore --- and replaces it with a sentence clarifying that
the closed form is a theorem. The empirical content is now exactly two
findings, both about validation of the K_h embedding.*

**STEP 8 DOWNSTREAM · §3.4, the paragraph defining α_det and the one
discussing the heuristic (lines 109--111)**

**What\'s wrong:** These paragraphs frame α_det\'s derivation as
heuristic in a way that, after Step 2, conflates two separate heuristic
inputs: (a) the K_h random-walk constant on the post-prefix state
(genuinely heuristic), and (b) the per-class structural offset (now
exact in the prefix algebra). Be clear about which is which.

**Current text in main.tex:**

By the standard random-walk heuristic on the post-prefix state, the
per-class intercept of \$\\sigma\$ vs \$\\log(n)\$ is heuristically:
\$\$\\alpha\_{\\text{det}}(r) := \\text{prefix\\\_steps}(r) + K_h \\cdot
\\log (a\_{\\text{final}}(r)/2\^k)\$\$ where
\$\\text{prefix\\\_steps}(r)\$ is the total number of Collatz steps in
the prefix from \$(a_0, c_0) = (2\^k, r)\$ until \$a\$ becomes odd (the
sum of even-step and odd-step applications during the prefix), and \$K_h
= 3/ \\log(4/3)\$. The prediction holds up to a single global additive
constant absorbed into the global intercept \$\\mu_a\$.

The derivation of \$\\alpha\_{\\text{det}}\$ is heuristic: it applies
the standard Collatz random-walk argument \\citep{lagarias1985} to the
post-prefix state \$a\_{\\text{final}} \\cdot m + c\_{\\text{final}}\$,
assuming the post-prefix dynamics behave as a random walk with i.i.d.
\$\\text{Geometric}(\\frac{1}{2})\$ 2-adic valuations. The heuristic is
the same one that produces the leading-order prediction
\$\\mathbb{E}\[\\sigma \\vert n\] \\approx K_h \\cdot \\log(n)\$ for
unconditioned \$\\sigma\$; \$\\alpha\_{\\text{det}}\$ adapts the
heuristic to a specific post-prefix initial configuration. We present
\$\\alpha\_{\\text{det}}\$ as a closed-form quantity computable from
\$r\$ and verify its empirical predictive content in \\S4. We do not
claim that \$\\alpha\_{\\text{det}}\$ is the exact per-class intercept;
we claim that it predicts the per-class intercept empirically.

**Replace with:**

We define the per-class structural offset \$\\alpha\_{\\text{det}}(r)\$
by

\$\$\\alpha\_{\\text{det}}(r) := \\text{prefix\\\_steps}(r) + K_h \\cdot
\\log(a\_{\\text{final}}(r)/2\^k) = (k + j(r)) + K_h \\cdot (j(r) \\log
3 - k \\log 2),\$\$

where \$j(r)\$ is the number of odd-step applications consumed during
the prefix iteration starting from \$r\$, \$\\text{prefix\\\_steps}(r) =
k + j(r)\$, \$a\_{\\text{final}}(r) = 3\^{j(r)}\$, and \$K_h =
3/\\log(4/3)\$. The quantity \$\\alpha\_{\\text{det}}(r)\$ is a
closed-form function of \$r\$ alone: both \$j(r)\$ and
\$\\text{prefix\\\_steps}(r)\$ are deterministic outputs of the prefix
iteration (Lemma\~\\ref{lem:prefix-determinism},
Theorem\~\\ref{thm:termination}).

The form of \$\\alpha\_{\\text{det}}\$ is motivated by the standard
Collatz random-walk heuristic \\citep{lagarias1985}: applying \$K_h
\\log(\\,\\cdot\\,)\$ to the post-prefix state \$a\_{\\text{final}}
\\cdot m + c\_{\\text{final}}\$ recovers the leading-order behavior of
\$\\sigma\$ on the post-prefix orbit, separating the deterministic
contribution of the prefix (which depends only on \$r\$) from the
stochastic contribution of the post-prefix dynamics (which depends on
\$m\$). This is the only heuristic input to the construction; we present
\\S\\ref{sec:empirics} as empirical validation of the per-class
predictive content of this heuristic.

**Why:** *Cleaner separation. Two paragraphs now do two things: (a)
define α_det as a closed-form function of r (deductive), (b) name K_h as
the single heuristic input. The old \"prediction holds up to a global
additive constant\" framing carried the implicit claim that the
per-class part was approximate --- after Step 2 the per-class part of
α_det is exact and that framing is misleading.*

**STEP 9 DOWNSTREAM · §4.3 (lines 184--188) --- closed-form verification
observation**

**What\'s wrong:** Currently framed as \"empirical verification of the
heuristic identity\". After Step 2 it\'s a numerical sanity check of a
theorem.

**Current text in main.tex:**

\\subsection{Verification of the closed form \$\\langle
\\alpha\_{\\text{det}} \\rangle = \\log(6)/ \\log(4/3)\$}

\\begin{observation}\[Closed-form \$\\langle \\alpha\_{\\text{det}}
\\rangle\$ verification\]

\\label{obs:closed-form-verification}

Computed empirically at each tested \$k\$ from per-class
\$\\alpha\_{\\text{det}}\$ values, the class-averaged \$\\langle
\\alpha\_{\\text{det}} \\rangle\$ matches \$\\log(6)/ \\log(4/3)
\\approx 6.226\$ to within sampling precision at \$k \\in \\{6, 7, 8,
9\\}\$. The empirical verification confirms the heuristic identity from
\\S3.5 at the cited resolutions.

\\end{observation}

**Replace with:**

\\subsection{Numerical verification of
Theorem\~\\ref{thm:alpha-det-closed-form}}

\\label{sec:closed-form-verification}

Theorem\~\\ref{thm:alpha-det-closed-form} predicts \$\\langle
\\alpha\_{\\text{det}} \\rangle = \\log 6 / \\log(4/3) \\approx
6.228263\$ exactly, independent of \$k\$. Direct computation of
\$\\alpha\_{\\text{det}}(r)\$ across all \$2\^{k-1}\$ odd residue
classes (no sampling) recovers this value to machine precision at every
tested resolution:

\\begin{center}

\\begin{tabular}{rrr}

\\hline

\$k\$ & \$2\^{k-1}\$ classes & \$\\langle \\alpha\_{\\text{det}}
\\rangle - \\log 6 / \\log(4/3)\$ \\\\

\\hline

6 & 32 & \$-5.3 \\times 10\^{-15}\$ \\\\

8 & 128 & \$-7.1 \\times 10\^{-15}\$ \\\\

10 & 512 & \$-7.1 \\times 10\^{-15}\$ \\\\

12 & 2048 & \$-8.9 \\times 10\^{-15}\$ \\\\

14 & 8192 & \$-1.1 \\times 10\^{-14}\$ \\\\

\\hline

\\end{tabular}

\\end{center}

This is a sanity check on Theorem\~\\ref{thm:alpha-det-closed-form}
rather than empirical evidence for it: the closed form is a deductive
consequence of Lemma\~\\ref{lem:j-binomial} and the definition of
\$\\alpha\_{\\text{det}}\$, both of which are non-empirical inputs. We
report the verification at \$k \\leq 14\$ because beyond that the
enumeration over \$2\^{k-1}\$ classes becomes expensive, not because the
result is range-bounded; the theorem holds for all \$k\$.

**Why:** *Honest framing. The numbers come from closed_form_findings.md
Result 1 (your own data). The new framing also expands the k-range to
{6,8,10,12,14} which matches that data --- broader than §4\'s
σ-regression range of k ∈ {6,7,8,9} because the closed-form verification
is enumeration, not sampling.*

**STEP 10 DOWNSTREAM · §5 \"Relation to Terras 1976\" (lines 194--203)
--- substantially rewrite**

**What\'s wrong:** After Step 2, the relationship to Terras is no longer
just narrative framing --- the proof of Lemma\~\\ref{lem:j-binomial}
literally uses Terras 1976 Lemma 1 (parity-vector bijection). Rewrite §5
to make this structural relationship explicit, and drop the duplicated
framing from §1/§2.2 that the current §5 mostly repeats.

**Current text in main.tex:**

\\section{Relation to \\citep{terras1976}}

\\citet{terras1976} Lemma 4 establishes the asymptotic identity for the
\$k\$th Collatz iterate\'s residue-class structure: \$S_k \\approx S_0
\\cdot 3\^{d(k)}/ 2\^k\$, where \$d(k)\$ is the number of odd-step
applications among the first \$k\$ Collatz steps starting from the
residue class. The prefix decomposition retains the symbolic state
\$(a,c)\$ through the deterministic prefix and terminates at
\$a\_{\\text{final}} = 3\^j\$ explicitly, where \$j\$ is the number of
odd-step applications consumed during the prefix.

There are two ways to view the relationship:

\\begin{enumerate}

\\item The prefix decomposition is a constructive version of Terras
Lemma 4 at finite \$k\$. Where Terras gives the asymptotic ratio \$S_k /
S_0 \\approx 3\^{d(k)} / 2\^k\$ for the iterate\'s residue-class
structure, the prefix decomposition gives the explicit terminal symbolic
state \$(a\_{\\text{final}}, c\_{\\text{final}},
\\text{prefix\\\_steps})\$ at finite \$k\$. The two results address
related but distinct objects: Terras\'s lemma describes the \$k\$th
iterate\'s residue distribution; the prefix decomposition describes the
symbolic state evolution from initial residue \$r\$.

\\item The prefix decomposition strengthens Terras Lemma 4 by retaining
the constructive content \$(a\_{\\text{final}}\$ and
\$\\text{prefix\\\_steps}\$) at finite \$k\$ as a closed-form covariate.
This constructive content predicts per-class statistics of \$\\sigma\$
at the level of the empirical regularities documented in \\S4. The
asymptotic content of Terras\'s lemma, by itself, does not produce these
per-class predictions.

We do not claim the prefix decomposition is novel relative to
\\citep{terras1976}. The decomposition is a sharpening of Terras Lemma 4
that retains constructive information Terras\'s lemma discards. The
acknowledgment is explicit and located here so that no reader
interpreting the decomposition as \"new\" relative to the broader
literature goes away with that impression incorrectly.

\\end{enumerate}

**Replace with:**

\\section{Relation to \\citet{terras1976}}

The prefix decomposition relates to \\citet{terras1976} in two places.
First, the proof of Lemma\~\\ref{lem:j-binomial} relies on
\\citet{terras1976} Lemma 1, the Syracuse parity-vector bijection
\$\\Pi_k \\colon \\mathbb{Z}/2\^k\\mathbb{Z} \\to \\{0,1\\}\^k\$. The
block decomposition of the prefix trace
(Lemma\~\\ref{lem:block-decomposition}) and its identification with the
Syracuse parity vector (Lemma\~\\ref{lem:syracuse-equivalence}) reduce
the count of odd classes with \$j\$ odd-step applications to counting
parity vectors of given Hamming weight under Terras\'s bijection. The
prefix-algebra content of this paper rests on Terras 1976 Lemma 1.

Second, the prefix decomposition is a constructive sharpening of
\\citet{terras1976} Lemma 4. Terras\'s Lemma 4 gives the asymptotic
ratio \$S_k \\approx S_0 \\cdot 3\^{d(k)}/2\^k\$ for the \$k\$th
iterate\'s residue-class structure but does not retain the symbolic
state past the asymptotic identity. The prefix decomposition retains the
symbolic state \$(a, c)\$ through the deterministic prefix, terminates
at \$a\_{\\text{final}} = 3\^j\$ explicitly, and exposes
\$a\_{\\text{final}}\$ and \$\\text{prefix\\\_steps}\$ as closed-form
covariates of per-class \$\\sigma\$ statistics at finite \$k\$. The
constructive content is what enables the per-class predictions of
\\S\\ref{sec:empirics}; Terras\'s asymptotic content, by itself, does
not produce them. We do not claim novelty relative to
\\citet{terras1976}; the prefix decomposition is a constructive use of
Terras 1976 Lemma 1 to give a finite-\$k\$ refinement of Terras 1976
Lemma 4.

**Why:** *§5 was repeating §1 and §2.2 framing. Now it does something
only §5 can do: explicitly trace the Terras 1976 → prefix-decomposition
dependency on both Lemma 1 (used in the proof) and Lemma 4 (sharpened by
the result). The new §5 is about half the length of the old one and
carries more information. A specialist who skips to §5 to check
positioning gets the actual story in two paragraphs.*

**STEP 11 DOWNSTREAM · §6.1(a) (line 213) and the renumbering of (b),
(c)**

**What\'s wrong:** §6.1(a) currently lists \"proof of all k values of
a_final realized for arbitrary k\" as an open question. After Step 3
it\'s Corollary\~\\ref{cor:realization}, proven. Remove the item;
renumber (b) and (c) to (a) and (b).

**Current text in main.tex:**

Three pieces of analytical work would strengthen the empirical
regularities of \\S4 toward theorem-grade results:

\\begin{enumerate}\[\\alph\*)\]

\\item Proof of \"all \$k\$ values of \$a\_{\\text{final}}\$ are
realized at every \$k\$\" for arbitrary \$k\$ (Empirical Observation 1
strengthened to a theorem). The realization at small \$k\$ can be
enumerated; whether the realization extends to arbitrary \$k\$ is a
question about the structure of the prefix iteration that we have not
resolved.

\\item Analytical characterization of the post-prefix dynamics from
state \$(a\_{\\text{final}}, c\_{\\text{final}}, m)\$. The standard
Collatz random-walk heuristic gives leading-order behavior on the
post-prefix state but does not capture the higher-moment structure that
the prefix-collapse predicts empirically. A characterization that
recovers the empirical Pearson correlations of \\S4.2 (variances,
kurtosis, skewness clustering by \$a\_{\\text{final}}\$) from the
post-prefix dynamics analytically would be a substantial strengthening.

\\item Connection to \\citep{tao2022}\'s almost-all-\$N\$ theorem. The
prefix decomposition\'s \$\\alpha\_{\\text{det}}\$ extends to predict
per-class mean first-passage time at arbitrary thresholds \$f(N)\$, with
the relationships \$s\_{\\text{mean}}(r; f) \\approx
\\alpha\_{\\text{det}}(r) + K_h \\cdot \\log(N/f(N))\$ verified
empirically across forty cells in the accompanying repository (and
reported separately). A theorem-grade version would derive the per-class
structural offset \$\\alpha\_{\\text{det}}\$ from Tao\'s framework
directly, connecting the prefix decomposition to the leading-term
content of Tao\'s (5.15).

\\end{enumerate}

**Replace with:**

Two pieces of analytical work would strengthen the empirical
regularities of \\S\\ref{sec:empirics} toward theorem-grade results:

\\begin{enumerate}\[\\alph\*)\]

\\item Analytical characterization of the post-prefix dynamics from
state \$(a\_{\\text{final}}, c\_{\\text{final}}, m)\$. The standard
Collatz random-walk heuristic gives leading-order behavior on the
post-prefix state but does not capture the higher-moment structure that
the prefix-collapse predicts empirically. A characterization that
recovers the empirical Pearson correlations of
\\S\\ref{sec:moment-clustering} (per-class variance, kurtosis, skewness
clustering by \$a\_{\\text{final}}\$) from the post-prefix dynamics
analytically would be a substantial strengthening.

\\item Connection to \\citep{tao2022}\'s almost-all-\$N\$ theorem.
Preliminary empirical work in the accompanying repository extends
\$\\alpha\_{\\text{det}}\$ to predict per-class mean first-passage time
at thresholds \$f(N)\$, of the form \$s\_{\\text{mean}}(r; f) \\approx
\\alpha\_{\\text{det}}(r) + K_h \\cdot \\log(N/f(N))\$; this extension
is reported separately and not claimed here. A theorem-grade version
would derive the per-class structural offset \$\\alpha\_{\\text{det}}\$
from Tao\'s framework directly, connecting the prefix decomposition to
the leading-term content of Tao\'s (5.15).

\\end{enumerate}

**Why:** *Removes the now-resolved open question and reframes the
Tao-bridge item more cautiously (\"preliminary empirical work\...
reported separately and not claimed here\") so the 40-cell mention
doesn\'t half-smuggle a claim. The first open question (post-prefix
dynamics) is genuinely open and important.*

**STEP 12 DOWNSTREAM · §7, first paragraph (line 228)**

**What\'s wrong:** Currently says \"The closed form ⟨α_det⟩ = log
6/log(4/3) is a heuristic identity derived under independent-Bernoulli
combinatorial assumption (Heuristic Identity 1) and is verified
empirically.\" After Step 2, all wrong.

**Current text in main.tex:**

The algebraic content: the prefix iteration is well-defined (Lemma 1),
terminated within at most \$k\$ odd-step applications (Theorem 1), and
produces a terminal \$a\_{\\text{final}} \\in \\{3\^1, 3\^2, \\ldots,
3\^k\\}\$. Theorem 1 is the deductive content of the symbolic state\'s
parity dynamics and is proven directly. The closed form \$\\langle
\\alpha\_{\\text{det}} \\rangle = \\log(6)/ \\log(4/3)\$ is a heuristic
identity derived under independent-Bernoulli combinatorial assumption
(Heuristic Identity 1) and is verified empirically.

**Replace with:**

The algebraic content: the prefix iteration is well-defined
(Lemma\~\\ref{lem:prefix-determinism}), terminates with exactly \$k\$
even-step applications and at most \$k\$ odd-step applications
(Theorem\~\\ref{thm:termination}), and produces a terminal
\$a\_{\\text{final}} \\in \\{3\^j : 1 \\leq j \\leq k\\}\$. The
\$j\$-distribution across the \$2\^{k-1}\$ odd classes mod \$2\^k\$ is
exactly binomial: \$\\binom{k-1}{j-1}\$ classes terminate at
\$a\_{\\text{final}} = 3\^j\$ (Lemma\~\\ref{lem:j-binomial}), with
realization of every \$j \\in \\{1, \\ldots, k\\}\$ as a corollary
(Corollary\~\\ref{cor:realization}).
Theorem\~\\ref{thm:alpha-det-closed-form} establishes \$\\langle
\\alpha\_{\\text{det}} \\rangle = \\log(6)/\\log(4/3)\$ exactly in the
prefix algebra, with the random-walk constant \$K_h\$ as the only
heuristic input.

**Why:** *Conclusion paragraph 1 (\"algebraic content\") must reflect
the new structure. This is the version a reader skimming to §7 sees as
the headline summary.*

**STEP 13 DOWNSTREAM · §7, second paragraph (line 230) --- the empirical
content paragraph**

**What\'s wrong:** Currently says α_det \"predicts per-class intercept
α(r) with R² ≥ 0.985 \... empirically verified predictive power\".
Should clarify: §4 validates the K_h embedding, not the closed form
(which is now Theorem).

**Current text in main.tex:**

The empirical content: at modular resolutions \$k \\in \\{6, 7, 8,
9\\}\$ and \$N\$ up to \$2\^{27}\$, the prefix decomposition\'s
\$\\alpha\_{\\text{det}}\$ predicts per-class intercept \$\\alpha(r)\$
of \$\\sigma\$ vs \$\\log(n)\$ with regressions \$R\^2 \\geq 0.985\$ and
signal-to-noise ratio \$\\text{SD(residual)} / \\text{mean(SE)}\$ in
\$\[0.91, 0.99\]\$; per-class central moments of \$\\sigma\$ residuals
cluster by \$a\_{\\text{final}}\$, with the number of distinct clusters
equal to \$k\$ at every tested \$k\$. The regularities are verified at
the cited resolutions and are not extrapolated. The relationship to
existing literature: the prefix decomposition is a sharpening of
\\citep{terras1976} Lemma 4. Where Terras\'s lemma gives the asymptotic
ratio for the \$k\$th iterate\'s residue-class structure, the prefix
decomposition retains the symbolic state through the deterministic
prefix and exposes \$a\_{\\text{final}}\$ and
\$\\text{prefix\\\_steps}\$ as closed-form covariates for per-class
\$\\sigma\$ statistics. The contribution of the present paper is the
constructive content (the explicit symbolic prefix and its empirically
verified predictive power), not novelty relative to Terras\'s asymptotic
identity.

**Replace with:**

The empirical content: at modular resolutions \$k \\in \\{6, 7, 8,
9\\}\$ and \$N\$ up to \$2\^{27}\$, the random-walk-heuristic embedding
\$K_h = 3/\\log(4/3)\$ underlying \$\\alpha\_{\\text{det}}\$ is
validated at the per-class level. Per-class intercept \$\\alpha(r)\$ of
\$\\sigma\$ vs \$\\log(n)\$ is predicted by
\$\\alpha\_{\\text{det}}(r)\$ with \$R\^2 \\geq 0.985\$ and
signal-to-noise ratio \$\\text{SD(residual)}/\\text{mean(SE)} \\in
\[0.91, 0.99\]\$; per-class central moments of \$\\sigma\$ residuals
cluster by \$a\_{\\text{final}}\$ into exactly \$k\$ distinct groups at
every tested \$k\$. These regularities are verified at the cited
resolutions and are not extrapolated. The contribution of this paper
combines the prefix-algebra results of \\S\\ref{sec:prefix}
(Theorem\~\\ref{thm:termination}, Lemma\~\\ref{lem:j-binomial},
Corollary\~\\ref{cor:realization},
Theorem\~\\ref{thm:alpha-det-closed-form}) with the empirical validation
of \\S\\ref{sec:empirics}; the prefix-algebra results constitute a
constructive sharpening of \\citet{terras1976} Lemma 4 using Terras\'s
own Lemma 1 (the Syracuse parity-vector bijection).

**Why:** *Distinguishes the now-deductive §3 content from the empirical
§4 content cleanly, and locates the Terras 1976 dependency explicitly.
The line \"empirically verified predictive power\" in the old text
conflated the closed form (now Theorem) with per-class predictions
(genuinely empirical) --- separating them removes the ambiguity.*

**STEP 14 DOWNSTREAM · §7, third paragraph (line 232) --- open questions
list**

**What\'s wrong:** Lists three open questions, but after Step 11 only
two remain.

**Current text in main.tex:**

Three open questions remain: rigorous extension of Empirical Observation
1 (realization of all \$k\$ \$a\_{\\text{final}}\$ values) to arbitrary
\$k\$; analytical characterization of the post-prefix conditional
distribution that recovers the empirical clustering of per-class central
moments by \$a\_{\\text{final}}\$; and the formal connection of
\$\\alpha\_{\\text{det}}\$ to \\citep{tao2022}\'s almost-all-\$N\$
theorem at the level of per-class structural offset.

**Replace with:**

Two open questions remain: analytical characterization of the
post-prefix conditional distribution that recovers the empirical
clustering of per-class central moments by \$a\_{\\text{final}}\$, and
the formal connection of \$\\alpha\_{\\text{det}}\$ to
\\citet{tao2022}\'s almost-all-\$N\$ theorem at the level of per-class
structural offset.

**Why:** *Mechanical update. Matches §6.1 after Step 11.*

**STEP 15 MAJOR · Appendix A, last sentence (line 258)**

**What\'s wrong:** The cited scripts (alpha_beta_gamma_decay.py and
analytical_abc_derivation.py) are NOT prefix iteration implementations.
alpha_beta_gamma_decay.py computes Markov chain sub-cell mass shares on
(Z/3\^k Z)\*; analytical_abc_derivation.py does inverse-tree
path-counting for \|μ̂(1/3)\|². Neither has anything to do with the
prefix iteration. The actual implementations, per your repo experiment
index, are in experiments/01_alpha_decomposition.py and
experiments/24_k_sweep_alpha_decomposition.py --- which Appendix B
already correctly cites. So this is an internal inconsistency between
Appendix A and Appendix B as well as factually wrong.

**Current text in main.tex:**

The iteration terminates when \$a\$ becomes odd (the loop\'s exit
condition). Theorem\~1 guarantees termination within at most \$2k\$
iterations. Reference implementations in the accompanying repository:
\\verb\|alpha\\\_beta\\\_gamma\\\_decay.py\| and
\\verb\|analytical\\\_abc\\\_derivation.py\|

**Replace with:**

The iteration terminates when \$a\$ becomes odd. By
Theorem\~\\ref{thm:termination}, exactly \$k\$ case-(i) iterations are
consumed and at most \$k\$ case-(ii) iterations are interleaved, so the
total iteration count is at most \$2k\$. Reference implementations of
the prefix iteration and per-class \$\\alpha\_{\\text{det}}\$
computation are in the accompanying repository at
\\verb\|experiments/01_alpha_decomposition.py\| (hierarchical Bayesian
fit at \$k=6\$, \$N=2\^{25}\$) and
\\verb\|experiments/24_k_sweep_alpha_decomposition.py\| (\$k\$-sweep at
\$N=2\^{27}\$ for \$k \\in \\{4,\\ldots,12\\}\$).

**Why:** *Independent of Step 2. Just wrong scripts cited. A referee or
repo-checker would notice immediately. Also tightens the \"at most 2k
iterations\" wording.*

**STEP 16 MAJOR · §4.1, currently the last sentence of Observation 2
(line 155)**

**What\'s wrong:** The matched-n controlled comparison is the strongest
empirical control in §4 --- it isolates \"α_det predicts well\" from
\"sample size grows, R² grows\". Currently buried as parenthetical
inside an Observation block. Promote it.

**Current text in main.tex:**

A controlled comparison: at \$k=8\$ with \$N = 2\^{27}\$ (524K per
class, matching the original \$k=6, N=2\^{25}\$ data scale of 524K per
class), \$R\^2 = 0.9918\$, essentially identical to the original \$k=6\$
\$R\^2\$ of \$0.9967\$ \-- same data per class produces same \$R\^2\$
regardless of modular resolution.

\\end{observation}

**Replace with:**

\\end{observation}

\\paragraph{Controlled comparison: matched per-class sample size.}

The decline in \$R\^2\$ across \$k \\in \\{6, 7, 8, 9\\}\$ in
Table\~\\ref{tab:residual-calibration} reflects shrinking per-class
sample sizes, not weakening of \$\\alpha\_{\\text{det}}\$ as a
predictor. At fixed per-class \$n\$, the \$R\^2\$ is essentially
independent of \$k\$:

\\begin{center}

\\begin{tabular}{rcrr}

\\hline

\$k\$ & \$N\$ & per-class \$n\$ & \$R\^2\$ \\\\

\\hline

6 & \$2\^{25}\$ & 524K & 0.9967 \\\\

8 & \$2\^{27}\$ & 524K & 0.9918 \\\\

\\hline

\\end{tabular}

\\end{center}

Same per-class sample size at different modular resolutions produces
essentially the same \$R\^2\$. The \$k\$-gradient in
Table\~\\ref{tab:residual-calibration} is sampling noise, not signal
degradation: as \$k\$ grows, the \$2\^{k-1}\$ classes share the same
total \$N\$, so per-class \$n\$ shrinks by half at each \$k\$, and
per-class \$\\alpha\$ posterior SE grows correspondingly. The
signal-to-noise band \$\[0.91, 0.99\]\$ in
Table\~\\ref{tab:residual-calibration} confirms residuals are at the
noise floor across all four resolutions.

**Why:** *Independent of Step 2. This is the single most load-bearing
empirical sentence in §4 and it should not be parenthetical. A
specialist who reads the controlled comparison sentence dismisses the
alternative interpretation \"maybe α_det just looks predictive because
there\'s more data\" --- without it, that\'s an open suspicion.*

**STEP 17 MAJOR · §4.2, between the Pearson table at line 175 and the
paragraph at line 177**

**What\'s wrong:** The \'k distinct clusters\' result is described
verbally, supported by Pearson correlations that actually understate the
result (the clustering is discrete, not linear). The paper currently has
no figures.

**Current text in main.tex:**

(no figure exists; add one)

**Replace with:**

Generate from your existing alpha_beta_gamma_values.csv plus
experiments/09_multi_stat_decomposition.py output. Suggested figure: 2×2
panel at k=8 (gives 8 clusters):

Panel (a): per-class Var(σ residual) vs α_det(r), colored by a_final.

Panel (b): per-class kurtosis vs α_det(r), same coloring.

Panel (c): per-class skewness vs α_det(r), same coloring.

Panel (d): number of distinct a_final realized vs k for k ∈ {6,\...,16}.

LaTeX inclusion stub:

\\begin{figure}\[ht\]

\\centering

\\includegraphics\[width=0.95\\textwidth\]{figures/moment_clustering_k8.pdf}

\\caption{Per-class central moments of \$\\sigma\$ residuals at \$k=8\$,
\$N=2\^{27}\$, plotted against \$\\alpha\_{\\text{det}}(r)\$ and colored
by \$a\_{\\text{final}}(r)\$. Variance (a), kurtosis (b), and skewness
(c) cluster into 8 discrete groups indexed by \$a\_{\\text{final}} =
3\^j\$ for \$j \\in \\{1, \\ldots, 8\\}\$, as predicted by
Lemma\~\\ref{lem:j-binomial}. Panel (d): the number of distinct
\$a\_{\\text{final}}\$ values realized as a function of \$k\$,
illustrating Corollary\~\\ref{cor:realization}.}

\\label{fig:moment-clustering}

\\end{figure}

**Why:** *Independent of Step 2 but the caption can now reference
Lemma\~\\ref{lem:j-binomial} and Corollary\~\\ref{cor:realization} from
the new §3.5. A specialist looking at panel (a) sees 8 horizontal bands
and is convinced in 5 seconds, where reading \"Pearson 0.99989 with
discrete substructure\" requires the reader to imagine that plot. Single
highest-leverage figure to add.*

**STEP 18 MAJOR · references.bib**

**What\'s wrong:** Five references is too thin for a Lagarias-targeted
paper. Three substantive gaps: Sinai 2003 (trajectory measure framing in
§2.1); Chang 2026 (parallel mod-8 work --- your STATE notes flag this
for citation); humphrey_collatz_2026 is in the bib but never cited (cite
in §1 as Paper 1 establishing the mod-8 case at k=3, or remove).

**Current text in main.tex:**

(see references.bib in the zip)

**Replace with:**

Add to references.bib:

\@article{sinai2003,

author = {Sinai, Yakov G.},

title = {Statistical (3x+1) problem},

journal = {Communications on Pure and Applied Mathematics},

volume = {56}, number = {7}, pages = {1016\--1028}, year = {2003}

}

\@misc{chang2026a,

author = {Chang, \[fill\]},

title = {\[fill\]}, year = {2026},

eprint = {2603.25753}, archivePrefix = {arXiv}

}

\@misc{chang2026b,

author = {Chang, \[fill\]},

title = {\[fill\]}, year = {2026},

eprint = {2603.11066}, archivePrefix = {arXiv}

}

\@misc{bonacorsi_bordoni_2026,

author = {Bonacorsi, Stefano and Bordoni, \[fill\]},

title = {\[fill\]}, year = {2026},

note = {fill from literature\\\_check.md}

}

Citations to add in the body:

§1 (after Step 5 update), add at end of contribution paragraph:

\"\...generalizing the \$k=3\$ (mod 8) case established in
\\citep{humphrey_collatz_2026}.\"

§2.3 (line 53), one sentence after the existing literature paragraph:

\"Independent parallel work by Chang \\citep{chang2026a, chang2026b}
obtains a related decomposition of the burst-gap indicator at \$K \\in
\\{3, 4, 5\\}\$ via block total-variation; Chang\'s framework addresses
a different statistical question than the prefix decomposition
introduced here.\"

§2.1 (line 47), after the existing K_h derivation:
\"\\citep{sinai2003}\" added to the Lagarias citation.

**Why:** *Pre-empts three referee comments: (1) \"aren\'t you aware of
Sinai?\", (2) \"how does this compare to Chang?\", (3) \"you cite your
own paper, where do you use it?\". Each cheap to fix now, expensive to
defend in revisions.*

Minor cleanup (apply after Step 18 or in any order)

**STEP M1 MINOR · §1, opening paragraph (line 33)**

**What\'s wrong:** Forward-reference to a_final before it\'s defined.

**Current text in main.tex:**

For odd integer \$n\$ with residue \$r = n \\bmod 2\^k\$, the Collatz
iteration starting from \$n\$ can be tracked symbolically. Writing
\$n=2\^k \\cdot m+r\$ and following the state \$a \\cdot m + c\$ through
Collatz iteration, the state\'s parity is determined by \$r\$ alone
(i.e., by \$a\\bmod 2\$ and \$c \\bmod 2\$ jointly) for the prefix
region while \$a\$ remains even, terminating when \$a\$ becomes odd. The
terminal \$a\_{\\text{final}}\$ lies in \$\\{3\^1, 3\^2, \\ldots,
3\^k\\}\$; the number of odd-step applications consumed is at most
\$k\$.

**Replace with:**

For odd integer \$n\$ with residue \$r = n \\bmod 2\^k\$, the Collatz
iteration starting from \$n\$ can be tracked symbolically. Writing \$n =
2\^k \\cdot m + r\$ and following the symbolic state \$a_t \\cdot m +
c_t\$ (initial state \$(a_0, c_0) = (2\^k, r)\$) through Collatz
iteration, the state\'s parity is determined by \$r\$ alone while
\$a_t\$ remains even. We call this the \\emph{prefix iteration} and
write \$a\_{\\text{final}}\$ for the value of \$a_t\$ at the first step
where it becomes odd. The terminal \$a\_{\\text{final}}\$ lies in
\$\\{3\^j : 1 \\leq j \\leq k\\}\$, and the number of odd-step
applications consumed is at most \$k\$.

**Why:** *Self-contained intro paragraph.*

**STEP M2 MINOR · §3.1 (line 57)**

**What\'s wrong:** Clarify m=0 edge case.

**Current text in main.tex:**

Fix an odd integer \$n\$ with residue \$r = n \\bmod 2\^k\$, and write
\$n = 2\^k \\cdot m + r\$, where \$m\$ is a non-negative integer.

**Replace with:**

Fix an odd integer \$n\$ with residue \$r = n \\bmod 2\^k\$, and write
\$n = 2\^k \\cdot m + r\$, where \$m \\geq 0\$ is an integer (with \$m =
0\$ corresponding to \$n = r\$, the smallest representative of the
class).

**Why:** *Defuses a small ambiguity.*

**STEP M3 MINOR · Appendix B (lines 260--270) --- reproducibility**

**What\'s wrong:** Doesn\'t specify σ definition or whether vanilla
Collatz vs Syracuse.

**Current text in main.tex:**

Detailed experimental records are in \\verb\|findings.md\| and
\\verb\|closed\\\_form\\\_findings.md\| within the repository.

**Replace with:**

\\paragraph{Conventions.} The total stopping time \$\\sigma(n)\$ is
defined as the number of Collatz steps required to reach \$n = 1\$ (not
the steps-to-first-descent-below-\$n\$ variant). Iteration uses the
unaccelerated Collatz map \$T(n) = 3n+1\$ if \$n\$ odd, \$n/2\$ if \$n\$
even, with no Syracuse compression. All computations of \$\\sigma\$ and
\$\\alpha\_{\\text{det}}\$ are deterministic; no random sampling is
used. Per-class \$\\alpha(r)\$ is estimated by OLS regression of
\$\\sigma(n)\$ on \$\\log(n)\$ within each residue class \$r \\bmod
2\^k\$ over \$n \\leq N\$.

Detailed experimental records are in \\verb\|findings.md\| and
\\verb\|closed_form_findings.md\| (the latter contains the derivation of
Theorem\~\\ref{thm:alpha-det-closed-form} as Result 1) within the
repository.

**Why:** *A specialist who wants to re-run a single class needs the
conventions.*

Closing note

After applying these 18 steps (plus the three minor items), the paper
has:

- One heuristic input (K_h from the random-walk argument).

- Two theorems (termination, closed-form ⟨α_det⟩).

- Three lemmas (prefix determinism, block decomposition, Syracuse
  equivalence, binomial j-distribution --- four if you count prefix
  determinism separately, which you should).

- One corollary (realization).

- Two genuine open questions (post-prefix conditional distribution; Tao
  2022 connection).

- An explicit structural dependency on Terras 1976 (both Lemma 1 used in
  the proof and Lemma 4 sharpened by the result).

This is what Paper 2 should look like for Lagarias. The empirical
content of §4 is now cleanly separable from the deductive content of §3
--- §3 stands on its own as theorem-grade prefix algebra, and §4
validates the K_h heuristic embedding empirically on per-class σ
statistics. A specialist can read either section first and the work
makes sense.

Slow is smooth.
