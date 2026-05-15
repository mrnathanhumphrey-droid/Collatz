# W4 Adversarial Audit — Faure √3 partial identification

**Date:** 2026-05-14
**Auditor:** Adversarial verification agent (independent of W4)
**Source files inspected (Mode E, verbatim where load-bearing):**
- `W4_DISPOSITION.md`, `W4_OPERATOR_SPECTRUM.md`, `W4_FAURE_VERBATIM.md`
- `w4_cumulant_spectrum.py`
- `result_77_T_lead_spectrum.md`, `bilinear_pair_operator.py`, `C1_TAO_RECURSION_FORM.md`
- `experiments_output/M_n_bilinear_moments.csv`
- `PADE_NUMERICAL_DISPOSITION.md`, `PADE_NUMERICAL_TRAJECTORY.md`
- `C:/tmp/faure/Faure_2009_Semiclassical_Spectral_Gap_Partially_Expanding_pdfminer.txt`
  (pdfminer-decoded Faure 2009 arXiv:0903.2747v1)

---

## Claim 1 audit (operator identification)

### Is T_dev pre-existing or new?

**FINDING: T_dev is a new object introduced by W4.** Grep of `C:/Collatz/`
shows the string "T_dev" appears ONLY in three files, all W4-generated:
`W4_DISPOSITION.md`, `W4_OPERATOR_SPECTRUM.md`, `w4_cumulant_spectrum.py`.

The pre-existing infrastructure (R75–R78) defines:
- `T_diag = (1/5)·[[1,1],[4,4]]` with spectrum `{0, 1}` (R77 Thm 77.1, **rigorous**)
- Full `T = T_diag + Off_n` with conjectural eigenvalues `{1, 1/2}` (R77 Conj 77.2,
  empirical at k=2..6)
- The bilinear pair-form `M_n(η) = Σ_ξ μ̂_n(ξ) μ̂_n*(ξ·η)` from
  `bilinear_pair_operator.py` (R76)

Nowhere in R75–R78 is there an object called "the bilinear deviation propagator
T_dev on {M_n(η) : η ≠ 1}". The W4 agent introduces this object and then
identifies its spectral radius with 1/√3. **The object is not derived from the
Tao recursion in any of the existing infrastructure files** — it is named and
its spectrum is asserted, but the level-to-level map d_n(η) → d_{n+1}(η) for
deviation vectors over (Z/3^n)* of changing dimension is never explicitly
written down. The W4 disposition itself flags this: "T_dev (bilinear deviation) | ≤ 1/√3 per step
| ≤ 1/√3 | Faure bound (analog); spectral radius ≤ 0.577" — the entry is
asserted, not derived.

**Verdict on this sub-claim:** the operator T_dev is a NEW posited construct,
not anchored in R75–R78. Identifying its spectrum requires either explicit
matrix entries (not given) or a proof that an asymptotic spectral radius
exists (not given).

### k=3 branching derivation soundness

**FINDING: The "k=3 branching → probability 1/3 → L² amplitude 1/√3" derivation
contains a category error.**

The Tao recursion (verbatim from `C1_TAO_RECURSION_FORM.md` §1 and
`bilinear_pair_operator.py` docstring lines 9):

> "μ̂_{n+1}(ξ) = Σ_v 2^{-v} e^{-2πi ξ 2^{-v}/3^{n+1}} μ̂_n(ξ·2^{-v} mod 3^n)"

The sum index `v` ranges over `{1, 2, 3, ...}` (truncated to `M = 2·3^{n-1}` in
the script) with iid **Geom(2)** weights `2^{-v}`. This is the **2-adic** part
of the recursion. From `C1_TAO_RECURSION_FORM.md` §2 — Tao's pair-grouping uses
`b_j := a_{2j-1} + a_{2j}` distributed as **Pascal(2, 1/2)** on ℕ+2; weights
remain 2-adic geometric.

There is **no per-step probability of 1/3** in the Tao recursion. The number 1/3
that the W4 agent invokes is `Σ_v 4^{-v} = 1/(4-1) = 1/3`, which is the
"diagonal v=v' aggregate weight" in the bilinear pair form. In R77 §1, the
relevant aggregate is:

> "Σ_{v even} 4^{-v} = 1/15 ... times 3 (level-n+1 to level-n cover factor) ..."

The "1/3" emerges arithmetically as the bilinear diagonal trace, not as a
trapping probability. Faure's "1/k probability of remaining trapped" is a
property of a SMOOTH expanding map `E(x) = kg(x) mod 1` on `S¹` where each
preimage of a point has cardinality `k`. The Syracuse/Tao recursion has
**countably infinite** v-domain with geometric weight, not k=3 equal-weight
branches. The two `k`'s are not the same: in Faure k is the preimage
cardinality of E; in W4's "3-adic branching" k is the cover index
|(Z/3^{n+1})*/(Z/3^n)*| = 3, which is a number-theoretic structural index, not
a dynamical fan-out.

`W4_FAURE_VERBATIM.md` §5 itself admits the hypothesis check:

> "Overall applicability of Faure 2009 Theorem 2 to Syracuse: **NOT APPLICABLE**
> (hypotheses fail)."

After admitting non-applicability the W4 agent re-imports the conclusion via the
"L² amplitude" mechanism — but the mechanism IS the Faure theorem (it's the
intuition behind the proof), so the value-transfer is logically illegitimate.

### L² vs bilinear pair moment interpretation

**FINDING: M_n(η) is a bilinear PAIR moment, not an L² amplitude in Faure's
sense.**

`bilinear_pair_operator.py` line 6:

> "M_n(η) := Σ_{ξ ∈ Z/3^n, 3∤ξ} μ̂_n(ξ) · μ̂_n*(ξ·η)"

This is a frequency-correlation moment (it equals the Fourier transform at η of
`|μ_n|²`, evaluated on (Z/3^n)*). It is real-valued in the sense of being the
expectation of an additive character composed with a sum, and `M_n(1) = S_n`
(Plancherel mass). It is NOT the L² norm of an operator iterate. The "L²
amplitude" framing in W4 conflates: (a) the L² norm of the Faure wave packet
projected onto trapped set K under iteration, and (b) the bilinear moment
M_n(η). These have different meaning and different scaling.

---

## Claim 2 audit (numerics)

### Independent re-computation from M_n_bilinear_moments.csv

Manual re-computation from the CSV (rows verbatim):

- Level 1, η=2: |M_1(2)| = 0.333333. Deviation L² = 0.333333. n_modes=1.
  norm_avg = 0.333333. **(matches W4)**

- Level 2, η∈{2,4,5,7,8}: |M_2(·)| = {0.095238, 0.238095, 0.095238, 0.238095,
  0.190476}. Sum of squares = 0.009070+0.056689+0.009070+0.056689+0.036281 =
  0.167800. L² = √0.167800 = 0.40964. n_modes=5. norm_avg = 0.40964/√5 =
  0.18325. **(matches W4)**

- Level 3, η∈{2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26}:
  17 modes, sum of squares = 0.383106, L² = 0.61895, n_modes=17. norm_avg =
  0.61895/√17 = 0.15013. **(matches W4)**

Ratios:
- Normalized L1→L2: 0.18325/0.33333 = **0.5497** ≈ 0.550 (matches W4)
- Normalized L2→L3: 0.15013/0.18325 = **0.8194** ≈ 0.819 (matches W4)

**The numbers themselves are correctly reported.**

### Ratio formula precision check

`W4_OPERATOR_SPECTRUM.md` §6 defines normalization as `‖d_k‖_2 / √n_modes`.
This is the **per-mode L² average**, not an operator norm.

**Critical observation: the RAW L² norms (no √n_modes division) are**:
- L1: 0.33333
- L2: 0.40964
- L3: 0.61895

So `‖d_k‖_2` is **GROWING** with k, not decaying. The "0.550" ratio is the
ratio of `‖d_k‖_2/√n_dev_k`, which decreases ONLY because `n_dev_k` grows by
a factor ≈ 3 per level (1 → 5 → 17, asymptotically (Z/3^n)*\{1} has
2·3^n − 1 elements, growing 3× per level).

Algebra: if `‖d_k‖²` stays bounded as k grows (which the data is consistent
with), then `‖d_k‖_2/√n_dev_k ~ const/√(3^k)`, giving a ratio of ~ 1/√3 by
**arithmetic of mode-counting**, not by spectral decay. The W4 framing
attributes the 0.550 to "the spectral radius of T_dev"; the data is equally
consistent with "‖d_k‖² is roughly conserved while mode count grows like 3^k",
which would make 1/√3 a TAUTOLOGY from coset growth, not a spectral statement.

**There is no W4-internal control distinguishing "spectral 1/√3" from
"counting-tautology 1/√3".**

### T_2 matrix construction soundness

The 6×6 T_2 matrix construction (`w4_cumulant_spectrum.py` lines 425–444) uses
the Tao recursion form `μ̂_{n+1}(ξ_out) = Σ_v 2^{-v}/Z_v · e^{-2πi ξ_out
2^{-v}/27} μ̂_2(ξ_out·2^{-v} mod 9)`. The construction is internally
consistent: it expresses one application of Tao's recursion at n=2 → n=3 as a
matrix on `(Z/9)*`. The eigenvalues `{0.792, 0.697, 0.557, 0.462, 0.352,
0.349}` look plausibly computed.

**However**: this matrix is the level-2-to-level-3 stepper LIFTED to act on
(Z/9)*-valued vectors. It is **not** the asymptotic transfer operator. There
is no analysis given that its spectrum converges as n → ∞. The claim that
"0.7921 is transient toward 1/√3" rests on **one** level computation (n=2)
plus a single PADE Hadamard-radius extrapolation; neither is a convergence
argument.

### "Transient" claim at 2 data points

**FINDING: The "transient" narrative requires the ratio to be MOVING TOWARD
1/√3 = 0.5774 over levels, but the data moves AWAY from it.**

The two ratios are:
- L1→L2 normalized: 0.550 (4.7% **below** 1/√3 = 0.5774)
- L2→L3 normalized: 0.819 (41.7% **above** 1/√3)

Direction is monotone **upward** (0.550 → 0.819), not toward 0.5774. The W4
explanation:

> "transient: the max norm persists but spreads over more modes"

is a verbal patch with no quantitative model of how 0.819 should subsequently
drop to 0.577. Two points define a slope and the slope is +0.269 per level —
extrapolating linearly gives 1.088 at L3→L4, not 0.577. The W4 agent does
NOT compute the L3→L4 ratio (which would require M_4, available in the script
infrastructure at `pis[4]` but not extracted into M_n_bilinear_moments.csv —
this is a glaring gap).

**Two data points moving in the wrong direction do not establish
"convergence toward 1/√3 as transient".**

---

## Claim 3 audit (Faure hypothesis applicability)

### Verbatim Faure hypotheses from PDF

From `C:/tmp/faure/Faure_2009_Semiclassical_Spectral_Gap_Partially_Expanding_pdfminer.txt`
line 1280–1318:

> "Theorem 2. Spetral gap in the semilassial limit.
> if the map f is partial ly aptive (de(cid:28)nition given page 15) (and m smal l
> enough), then the [...] spetral radius of the operator [...] does not depend
> on m and satis(cid:28)es in the [...] semi-lassial limit ν → ∞:
>  r_s(F̂_ν) ≤ 1/√E_min + o(1)   (11)
> whih is strily smal ler than 1 from (3)."

Hypotheses (from the pdfminer extraction, page 2-4): the map E: S¹ → S¹ is
`E(x) = kg(x) mod 1` with `g` a C^∞ diffeomorphism of S¹, and the skew product
`f: T² → T²` is `(x, s) ↦ (E(x), s + τ(x)/(2π) mod 1)`. F̂_ν is the
semiclassical transfer operator on `H^m(S¹)` with frequency `ν`. The
"partially captive" condition refers to the trapped set on the cotangent
dynamics having sub-exponential growth N(n) = O(1).

`W4_FAURE_VERBATIM.md` §5 explicitly lists 5 hypotheses and marks all FAIL for
Syracuse:
1. C^∞ smooth compact T²: Syracuse uses profinite Z_3*
2. Uniformly expanding base E_min > 1: Syracuse has stochastic Geom(2),
   no manifold
3. C^∞ skew-product structure: Syracuse uses Tao renewal product
4. Pseudodifferential calculus on T*S¹: no cotangent on Z_3*
5. Anisotropic Sobolev H^m: no smooth structure

These 5 FAILs are correctly reported.

**However**: in Faure's own Figure 2 example (line 1359), he uses **k=2**, not
k=3. The k=3 reading is NOT in Faure's paper; it is W4's projected analog for
Syracuse via "3-adic cover factor". This is a paraphrase, not a citation.

### Legitimacy of value-transfer when hypotheses fail

**FINDING: The value transfer is logically unsupported.**

The W4 disposition says:

> "Overall applicability of Faure 2009 Theorem 2 to Syracuse: NOT APPLICABLE
> (hypotheses fail).
> The CONCLUSION-SIDE numerical value √3 ≈ 1.732 is empirically consistent
> with the PADE Hadamard radius trajectory ... but Theorem 2 cannot be formally
> invoked."

Yet `W4_DISPOSITION.md` §1 immediately re-imports it:

> "**√3 = (1/r_s(T_dev))** ... via the L²-amplitude mechanism of Faure 2009
> Theorem 2 with k=3."

The "L²-amplitude mechanism" IS the content of Faure's Theorem 2. There is no
separate mechanism that survives when the smooth-manifold and pseudo-
differential hypotheses fail. The W4 agent is invoking Faure's proof-intuition
"`1/k` probability → `1/√k` amplitude" as a self-standing principle, but this
intuition is itself proved under those failing hypotheses. **You cannot
transfer the value when the mechanism that produces the value requires
exactly the hypotheses that fail.**

The W4 disposition's own §4 admits this:

> "r_s(T_dev) ≤ 1/√3 (analogy) | **CONJECTURAL EXTENSION** (profinite analog not proved)"
> "√3 = 1/r_s(T_dev) identification | **PARTIAL** — mechanism clear, rigorous
>  proof requires profinite semiclassical analysis"

This is honest about the gap, but the headline claim ("√3 IS the inverse
spectral radius of T_dev") is asserted in §1 of the disposition without the
caveat. The headline overreaches the caveat-section.

### PADE-to-√3 extrapolation soundness

**FINDING: The PADE trajectory does NOT support √3 as the asymptotic.**

`PADE_NUMERICAL_TRAJECTORY.md` reports Hadamard radius from `|ε_n|^(1/n)`:

| n  | implied ρ |
|----|-----------|
| 10 | 2.06 |
| 11 | 1.81 |
| 12 | 1.66 |
| 13 | 1.57 |

The trajectory passes THROUGH √3 = 1.732 between n=11 (1.81 > √3) and
n=12 (1.66 < √3). It does NOT stop at √3 — it continues inward to 1.57 at
n=13.

`PADE_NUMERICAL_DISPOSITION.md` headline:

> "H_TWO_SINGULARITIES_VISIBLE ... the Hadamard radius estimate from |ε_n|^(1/n)
> at n=10..13 trajectory (2.06 → 1.81 → 1.66 → 1.57) places the leading
> singularity at |z| ≈ 1.57 at n=13, with **monotone-inward trend**"

and

> "**z=2 is empirically refuted as the leading singularity at large n**"

and the expected asymptotic per STATE.md (cited in PADE_NUMERICAL_DISPOSITION):
**z ≈ 1.016**, NOT √3.

W4_DISPOSITION.md §1 claims:

> "(PADE trajectory 1/r_s converging to √3 from below starting at 1.57 at n=13)"

This is **factually inverted**. 1.57 is BELOW √3, but the trajectory is
descending, not ascending — so "from below" (which would mean approaching √3
from values smaller than √3) is the wrong direction-word, AND the trajectory's
projected asymptote per the source file is 1.016, not √3.

**The PADE evidence cited by W4 actively refutes the W4 reading.**

---

## Claim 4 audit (PADE consistency)

### Convergence direction verification

Restating: the data is 2.06 → 1.81 → 1.66 → 1.57. This is:
- **monotone decreasing** ✓ (verified)
- passes through √3 ≈ 1.732 between n=11 and n=12 (verified)
- direction of motion at n=13: still decreasing toward smaller values
  (PADE_NUMERICAL_DISPOSITION confirms "monotone-inward trend")
- expected asymptote per cited STATE.md: 1.016 (from the slow-mode ρ ≈ 0.984)

The W4 agent treats √3 as the asymptote because √3 sits in the middle of the
n=10..13 window. But the trajectory does not slow at √3 — it accelerates
through it (steps: -0.25, -0.15, -0.09; decelerating but still moving inward).
**There is no signal in the PADE data that singles out √3** as the asymptotic
target rather than any other value in [1.0, 1.7].

`PADE_NUMERICAL_TRAJECTORY.md` explicitly says:

> "For the slow-mode prediction (z = 1/0.984 ≈ 1.016) to be the asymptotic,
> the Hadamard trajectory must continue inward by another factor of ~1.5"

i.e., the source-file reading is that 1.016 IS the asymptote and n=13 is not
yet near it. W4 cherry-picks √3 from the same data.

---

## Verdict: FAILED (with caveats)

The W4 partial identification of √3 with the inverse spectral radius of a
posited bilinear deviation propagator **fails adversarial verification on
multiple load-bearing steps**:

1. **T_dev is invented, not derived.** It does not appear in R75–R78
   infrastructure. Its matrix entries are never written down. Its
   "spectrum ≤ 1/√3" is asserted, not computed.

2. **The "k=3 branching → 1/3 probability → 1/√3 amplitude" chain conflates
   two distinct k's.** Faure's k is the preimage cardinality of a smooth
   expanding map; the Tao recursion has Geom(2) weights summed over countable
   `v`, with the "3" arising as the coset-cover index of (Z/3^{n+1})*/(Z/3^n)*.
   These are different mathematical objects with the same numeral.

3. **The numerical "0.550 close to 1/√3" is plausibly a coset-counting
   tautology.** Dividing raw L² by √(mode count) where mode count grows ~3×
   per level introduces a 1/√3 factor by construction. The raw L² norms grow
   (0.333 → 0.410 → 0.619); the only thing decaying is the per-mode average.

4. **The "transient" narrative is unsupported by the data.** The two computed
   ratios move AWAY from 1/√3 (0.550 → 0.819, direction +0.269 per level).
   No L3→L4 ratio is reported despite the infrastructure supporting it
   (`pis[4]` is built in the script).

5. **The PADE trajectory cited as supporting √3 actively refutes it.** The
   trajectory passes through √3 monotonically en route to a much smaller
   asymptote (per cited STATE.md: ρ ≈ 1.016). The W4 disposition's claim
   "converging to √3 from below" is directionally wrong (1.57 < √3, but
   trajectory is still descending) and asymptotically wrong (asymptote is
   not √3).

6. **The value transfer from Faure 2009 is logically invalid.** The
   "L²-amplitude mechanism" IS the proof of Faure's Theorem 2. The W4
   disposition admits all 5 Faure hypotheses fail for Syracuse, then
   re-imports the conclusion via the very mechanism that proves the theorem
   under those (failing) hypotheses.

**What survives:**
- The numerical values reported (0.550, 0.819, 0.7921, the M_n(η) CSV)
  reproduce correctly when recomputed.
- The honest acknowledgment that Faure's 5 hypotheses FAIL for Syracuse
  (§5 of W4_FAURE_VERBATIM and the §4 status table of W4_DISPOSITION).
- The disposition file does flag the identification as PARTIAL.

**What fails:**
- The headline-level claim that √3 = 1/r_s(T_dev) is a "mechanism-clear
  partial identification". The mechanism is not clear (it requires the failed
  hypotheses), the object T_dev is not defined, the numerics are ambiguous
  between spectral and tautological, and the PADE evidence is misread.

**Sharpest error (a single quotable inversion):**

W4_DISPOSITION.md §1 states the PADE evidence as
"(1/r_s converging to √3 from below starting at 1.57 at n=13)". The PADE
trajectory is 2.06 → 1.81 → 1.66 → 1.57; it passes through √3 ≈ 1.732
between n=11 and n=12 and **continues to descend**. 1.57 is not "approaching
√3 from below" — it is on the opposite side of √3, still moving away. The
cited source file (`PADE_NUMERICAL_TRAJECTORY.md`) states the trajectory's
expected asymptote is z ≈ 1.016, not √3.

---

## Compared to W1 + W2 audits

I do not have a W1 or W2 audit file in `C:/Collatz/` to compare against
(no files matching `W1_AUDIT*`, `W2_AUDIT*` under the prior probe sweep).
This audit was conducted independently of any prior adversarial frame.

**Rigor self-assessment:** This audit is more rigorous than the W4
agent's self-assessment because (a) it independently re-derives the CSV
numerics; (b) it tests the W4 PADE reading against the cited source
disposition file rather than W4's paraphrase of it; (c) it identifies a
load-bearing inversion between W4's text and the source PADE disposition;
(d) it tests the load-bearing derivation step ("1/3 per step") against the
verbatim Tao recursion and finds a category error (Geom(2) vs k=3); (e) it
identifies the coset-counting tautology that fits the "0.550 ≈ 1/√3" datum
equally well as the spectral hypothesis, and notes that W4 supplies no
control distinguishing the two.

**Caveats:**
- Independent computation was performed by hand from the CSV; PowerShell
  and Bash were both denied. A scripted recomputation would tighten the
  L³ sum of squares to higher decimal places (the hand arithmetic carries
  ~5 sig figs).
- The L3→L4 normalized ratio (which would decisively test "transient toward
  1/√3 = 0.577 vs. monotone divergence") was not computed; the
  infrastructure is in place (`pis[4]`, `compute_M_n` at k=4) but is not run.
  Recommend extracting M_4(η) and computing the L3→L4 normalized ratio as the
  most decisive single experiment for this question.
- The "T_dev is well-defined" question is partly a matter of philosophy: one
  could *posit* T_dev as the limiting operator in a yet-to-be-constructed
  injective system on profinite L²((Z/3^n)*). The audit's claim is the
  weaker one that W4 has not constructed it, not that it cannot be
  constructed.
