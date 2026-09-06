# RECON — lifting the channel engine (R7–R10) to general multiplier-equivariant Markov chains. **ALTITUDE: ADD-A-PEAK.** (2026-09-06)

Lit-recon + shopping for the generalization of the R7–R10 "channel engine" from the Syracuse chain to a general
theorem on Markov chains on finite abelian groups equivariant under a multiplier `x↦qx`. Spec: Wilson (kill-first
recon). Three background agents swept three distinct massifs (Diaconis rep-theory; lumping/symmetry-reduction;
abelian Fourier + Ramanujan/supercharacter bridge); ~24 web searches total; the two definitive PDFs read in full
via local `pdftotext` (WebFetch decode-failures flagged below). **Do NOT re-derive R7–R10 — pure recon.**

## The target being placed
> Markov chain on a finite abelian group `G`, EQUIVARIANT under `x↦qx`, stationary `μ` ⟹ the off-diagonal
> L²/Plancherel correlation mass decomposes as a weighted sum of **twisted ⟨q⟩-orbit character sums** of `μ`,
> **vanishing identically iff `μ` is Haar**. (Syracuse instance: `OffDiag_k = (2/3)Σ_m 4^{−m} C_k(m)`, `C_k(m)` a
> twisted `⟨4^m⟩`-orbit character sum of `μ_{k−1}`, period `3^{k−1}`; Mersenne `4^P−1` = character-group trace of
> the weight `1/(4ω−1)`.)

---

## ALTITUDE VERDICT: **ADD-A-PEAK** (known range; the assembled theorem is a real, unpublished peak, but a short climb on well-cited legs).

Not **ALREADY-DONE**: none of the three massifs contains the assembled statement (a *stationary, possibly non-Haar*
`μ`; *off-diagonal Plancherel correlation mass*; expanded in *twisted `⟨q⟩`-orbit character sums on the dual*; with
a *vanishing-iff-Haar* dichotomy). Not clean **FIRST-ASCENT** either — and the guardrail's high bar is why: this is
a heavily-worked range, every *component* is standard, and two competitors reach the same orbit-on-the-dual
structure by a different route. The honest call is **ADD-A-PEAK**: the range (equivariant harmonic analysis of
Markov chains + supercharacter/Ramanujan orbit sums) is settled; our contribution is the specific **wiring** — the
join of two well-cited legs — plus the off-diagonal-correlation **vanishing dichotomy** under a *single-multiplier*
`⟨q⟩` sub-orbit twist. Whether that peak merits a standalone theorem vs a remark is an expert call (Zulip / a
Diaconis-circle reader), same posture as the denominator-criterion disposition.

### Sub-verdict (1) — UNIFORM-KILL ("correlation vanishes iff μ = Haar"): **KNOWN / base camp.**
Textbook character orthogonality, then packaged as the Diaconis–Shahshahani upper-bound lemma. This is
acclimatization, not the summit.
- Character orthogonality (Tao, 245C Notes 2, Lemma 17, **verbatim**): *"for a multiplicative character χ, the
  integral ∫_G χ dμ equals 1 when χ is trivial and 0 when χ is nontrivial."* Haar is the unique measure with zero
  Fourier coefficient at every nontrivial χ (Fourier uniqueness) ⟹ uniform-kill.
- Markov packaging — Diaconis–Shahshahani **upper-bound lemma**, abelian form (*paraphrase*, as reproduced across
  the mixing literature): `‖P − U‖²_{TV} ≤ (1/4) Σ_{ρ≠id} dim(ρ)·tr(P̂(ρ)P̂(ρ)*)` → on an abelian group,
  `4‖μ − U‖²_{TV} ≤ Σ_{g≠0}|μ̂(g)|²`. The L² gap to uniform = sum of squared *nontrivial* Fourier mass, zero iff
  `P = U`. **This is the coarse flat-Fourier ancestor of our correlation identity.** Cite Diaconis 1988 + Terras 1999.

### Sub-verdict (2) — the ⟨q⟩-ORBIT REFINEMENT (the altitude decider): **the pieces exist SEPARATELY; the wiring is UNPUBLISHED.**
Two published half-legs, standing apart:
- **Leg A (uniform-kill):** Diaconis/Terras (above).
- **Leg B (Ramanujan sum = orbit sum of additive characters under the units action):** KNOWN and NAMED. The
  Ramanujan sum is *by definition* an orbit sum (Encyclopedia of Math, **verbatim**): `c_k(n) = Σ_h exp(2πi·n·h/k)`,
  *"where h ranges over all non-negative integers less than k that are relatively prime to k (the units modulo k)."*
  Elevated to the **supercharacter** framework by Fowler–Garcia–Karaali, *Ramanujan sums as supercharacters*
  (abstract **verbatim**): *"The theory of supercharacters, recently developed by Diaconis–Isaacs and Andre, can be
  used to derive the fundamental algebraic properties of Ramanujan sums."* — `G = ℤ/nℤ`, `Γ = (ℤ/nℤ)×` acts, and the
  supercharacters are sums of additive characters over `Γ`-orbits = the Ramanujan sums. **KEY GAP:** this uses the
  FULL unit group `Γ`, not the cyclic subgroup `⟨q⟩ = ⟨4⟩`. **Our single-multiplier `⟨4⟩` twist (period `3^{k−1}`)
  is a *sub*-orbit refinement of the supercharacter/Ramanujan construction** — that specific `⟨q⟩`-subgroup twist
  is not published. And the Mersenne `4^P−1` trace is downstream of this same orbit-sum machinery.

What is **absent everywhere searched**: the *join* — a `q`-equivariant Markov correlation expanded in twisted
`⟨q⟩`-orbit character sums of a general stationary `μ`, vanishing iff Haar. Leg A and Leg B are two unconnected
published results; the wiring is the contribution.

---

## LOAD-BEARING HYPOTHESIS — **CONFIRMED: multiplier-equivariance is ESSENTIAL (not arbitrary chains).**
Both competitors build on it and neither survives without it:
- Boyd–Diaconis–Parrilo–Xiao: the whole construction is on chains that **commute with the group action** (their
  object is literally *"the set of all Markov chains on V that commute with the action of G"*).
- Dubail–Massoulié: the Fourier factorization is derived *"by the linearity of the dynamics… its Fourier transform
  factorizes."*
⟹ our target's hypothesis (multiplier-equivariant kernel, not arbitrary chains) is the **right** load-bearing
assumption. Drop equivariance and the orbit structure collapses to flat per-character.

---

## NEAREST PUBLISHED STATEMENTS TO THE TARGET (verbatim) + the exact gap

**1. Dubail–Massoulié, "Accelerating Abelian Random Walks with Hyperbolic Dynamics" (arXiv:2106.10079) — the nearest MECHANISM.**
> **Lemma 1.** For all `t ≥ 1`, for all `ξ ∈ T_n`, `P̂ᵗ(ξ) = ∏_{j=0}^{t−1} μ̂((Aᵀ)ʲ ξ)`. *(Proof: "= P̂ᵗ⁻¹(Aᵀξ)·μ̂(ξ) … which gives the result by induction.")* — fed into the upper-bound lemma `4‖μ − U‖²_{TV} ≤ Σ_{g≠0}|μ̂(g)|²`.

This is the exact multiplier-twisted-orbit structure on characters (the `⟨Aᵀ⟩`-orbit of `ξ`), driven by an
equivariant/linear kernel — the mechanistic sibling of `C_k(m)`. **Exact gap:** (a) their orbit object is a
**product along a forward time-trajectory** governing a **mixing rate**, not a **weighted sum decomposing stationary
correlation mass**; (b) their stationary measure is **uniform by construction**, so there is **no non-Haar μ** and
therefore **no vanishing-iff-Haar dichotomy**. Our off-diagonal-mass estimand + the dichotomy is the delta.

**2. Boyd–Diaconis–Parrilo–Xiao, "Symmetry Analysis of Reversible Markov Chains" (Internet Math. 2005) — the nearest ORBIT-CHAIN.**
> **Theorem 3.1.** *"Let G be the automorphism group of the reversible Markov chain (K,π). Suppose V = O₁ ⊔ … ⊔ O_k as a disjoint union of G-orbits… Then all eigenvalues of K occur among the eigenvalues of {K_{H_i}}… every eigenfunction of K occurs by translating a lift of an eigenfunction of some K_{H_i}."*
> (+ block-diagonalization: `ρ(g)P = Pρ(g) ∀g ⟹ T⁻¹PT` block diagonal, one block per isotypic component; abelian ⟹ 1-D character blocks.)

**Exact gap:** these orbits are of the automorphism group on the **state space V** (lumped quotient chains,
eigenvalue bookkeeping), **not** orbits of a multiplier `x↦qx` on the **dual/character group**; and it is a
spectrum-preserving block-diagonalization, with **no twisted character sum, no 4⁻ᵐ weighting, no Mersenne-trace
denominators, no vanishing criterion.** It block-diagonalizes and halts.

**3. Ayyer–Singla, "Random Motion on Finite Rings I" (arXiv:1605.05089) — the nearest CHARACTER-SUM eigenvalues.**
Character-theoretic eigenvalues of an additive-multiplicative chain on a finite abelian group — but the eigenvalue
formula is a **character sum over a coset `F_a`, not over `⟨q⟩`-orbits on the dual**, and carries **no twist**
("orbit" appears nowhere in it).

**4. Tao, "Almost all orbits of the Collatz map…" (Forum Math. Pi 10, 2022; arXiv:1909.03562) — the nearest ANCESTOR (our exact setting, coarse form).**
Defines the Syracuse random variable (**verbatim**): `Σ_{m=1}^n 3^{n−m} 2^{−a_m−⋯−a_n}`, the `a_i` i.i.d. geometric,
and *"estimates the characteristic function of a skew random walk on a 3-adic cyclic group ℤ/3ⁿℤ at high
frequencies."* **Exact gap:** this is the **flat-Fourier** version — Fourier coefficients of the Syracuse law to
prove approach-to-uniform — with **no `⟨q⟩`-orbit refinement, no Ramanujan-sum closed form, no Markov-correlation
identity.** This is the right thing to cite as "the coarse version of what we sharpen."

---

## SHOPPING HAUL (climbing gear, tagged)
1. **[base-camp]** P. Diaconis, *Group Representations in Probability and Statistics*, IMS Lecture Notes 11 (1988) — the upper-bound lemma; deviation-from-uniform = nontrivial Fourier mass.
2. **[base-camp]** A. Terras, *Fourier Analysis on Finite Groups and Applications*, Cambridge (1999) — companion; character theory + Ramanujan sums on finite abelian groups.
3. **[base-camp / climbing tool]** Encyclopedia of Mathematics, *Ramanujan sums* — https://encyclopediaofmath.org/wiki/Ramanujan_sums — verbatim additive-character-over-units (orbit-sum) definition.
4. **[climbing tool — the orbit-sum bridge]** C. Fowler, S. R. Garcia, G. Karaali, *Ramanujan sums as supercharacters*, Ramanujan J. 35 (2014) — https://arxiv.org/abs/1201.1060 — Ramanujan sum = sum of additive characters over `(ℤ/nℤ)×`-orbits (our `⟨q⟩` twist is a sub-orbit refinement). Foundational: Diaconis–Isaacs, *Supercharacters and superclasses for algebra groups*, Trans. AMS 360 (2008).
5. **[competitor]** Boyd, Diaconis, Parrilo, Xiao, *Symmetry Analysis of Reversible Markov Chains*, Internet Math. 2(1):31–71 (2005) — https://web.stanford.edu/~boyd/papers/pdf/symmetry.pdf — orbit chains + isotypic block-diagonalization (state-space orbits).
6. **[competitor]** Boyd, Diaconis, Parrilo, Xiao, *Fastest Mixing Markov Chain on Graphs with Symmetries*, SIAM J. Opt. 20(2):792–819 (2009) — https://arxiv.org/abs/0709.0955 — orbit-lumping (Eq. 16) + Thm 4.1/5.1 block structure.
7. **[competitor — nearest mechanism]** J. Dubail, L. Massoulié, *Accelerating Abelian Random Walks with Hyperbolic Dynamics* (2022) — https://arxiv.org/abs/2106.10079 — Fourier factorization over the `⟨Aᵀ⟩`-orbit of a character (Lemma 1), equivariant kernel; forward-time product, uniform μ.
8. **[competitor — character-sum eigenvalues]** A. Ayyer, P. Singla, *Random Motion on Finite Rings I: Commutative Rings* (2020) — https://arxiv.org/abs/1605.05089 — coset-indexed (not `⟨q⟩`-orbit), no twist.
9. **[climbing tool]** *Diaconis–Shahshahani upper bound lemma*, abelian ℤ/n form as used in cutoff work — e.g. https://arxiv.org/abs/1407.3580 (*Mixing Time and Cutoff for a Random Walk on the Ring of Integers mod n*).
10. **[motivating example — nearest ancestor]** T. Tao, *Almost all orbits of the Collatz map attain almost bounded values*, Forum Math. Pi 10 (2022) — https://arxiv.org/abs/1909.03562 — Fourier coefficients of the Syracuse distribution on ℤ/3ⁿ (coarse form of our target).
   - *Leads (not read verbatim — treat as pointers):* Asci, affine walk `X_{n+1}=aX_n+b_n` on `(ℤ/nℤ)^d` (2001, the source Dubail–Massoulié attribute Lemma 1 to); Hildebrand, *Enumeration and random random walks on finite groups*, Ann. Probab. 24 (1996); Saloff-Coste, *Random Walks on Finite Groups* (2004); *Schur–Weyl duality for diagonalizing a Markov chain on the hypercube*, arXiv:2512.23285.

## Read-failures (flagged, NOT nulls)
- WebFetch could not decode the binary/FlateDecode PDFs for arXiv:0709.0955, arXiv:1605.05089, arXiv:2106.10079,
  Garcia's `RSS.pdf`, Garrett's `fin_ab_fourier.pdf`, arXiv:1407.3580, and a Berkeley lecture PDF (302→404).
  Recovered via: local `pdftotext` (0709.0955 Symmetry-paper sibling `symmetry.pdf`, and 2106.10079 — read in full,
  the Dubail–Massoulié Lemma 1 and Boyd–Diaconis theorem quotes are from these); ar5iv/HTML renderings (0709.0955,
  1605.05089); and abstract pages (1201.1060). Springer typeset Ayyer–Singla is auth-gated (used ar5iv preprint).
  Diaconis 1988 itself: **only Google-Books metadata readable** — its content is inferred from competitors' citations
  of "[Dia88] Chapter 1," not read verbatim.

## Skeptic's note — why the wiring is judged UNPUBLISHED
~24 searches across the three massifs (equivariant/commuting-action isotypic + orbit chains; automorphism/multiplier-
twisted walks; abelian Fourier + Dirichlet-character multiplicative walks; supercharacters/Ramanujan orbit sums;
"off-diagonal Plancherel correlation vanishes iff uniform"; Syracuse/Collatz Fourier). Every hit lands in one of
three separate buckets — flat per-irrep (Diaconis–Shahshahani), state-space-orbit chains (Boyd–Diaconis), or the
forward-time mixing-rate orbit product (Asci/Dubail–Massoulié) — plus the standalone Ramanujan-sum-as-orbit-sum fact
(supercharacters). None joins "Markov L² off-diagonal correlation mass" to "twisted `⟨q⟩`-orbit character sum of a
general stationary μ" with a "vanishing-iff-Haar" iff. **The range is known; the specific peak (the join + the
dichotomy under a single-multiplier `⟨q⟩` sub-orbit twist) is the contribution.** Bar met for ADD-A-PEAK, not for a
clean FIRST-ASCENT.
