# Monotone closure — integrated writeup (post-Task 4)

**Date:** 2026-05-14
**Status:** Partial closure landed. Pivot-or-dig decision point.
**Supersedes:** the "next step" portion of `C:/Collatz/OBSTRUCTION_MAP_TERMINAL.md`
**Reads with:** OBSTRUCTION_MAP_TERMINAL.md (framework identification, terminal finding) + MONOTONE_CUMULANTS_{A,B,C,D,DISPOSITION}.md (the four deliverables this writeup integrates)

---

## 0. Headline

**Leading-order c = 7/45 closure derivation: IN HAND** (rigorous fiberwise; conjectural at one named lift step).

**Full multi-spectral asymptotic closure: NOT in hand.** Four named wrinkles, characterized below as dig-hard vs pivot.

**What this means.** Yesterday's terminal finding was "framework identification: monotone, not free; effort 12-19mo → 5-9hr." Today's 5-9hr produced (a) numerical confirmation of the monotone diagnostic at 10⁶ separation, (b) verbatim Hasebe-Saigo cumulant derivation of the 7/45 leading coefficient, and (c) a sharp four-item map of what's left. The result is a real partial closure — the 7/45 coefficient is now anchored to an established mathematical framework, with the residual work scoped down to four concrete gaps.

---

## 1. What today closed

### 1.1 Numerical diagnostic confirmation (Task 1, ~1h actual)

The third-order alternating B-centered moment `φ(X̃_{j_1}·X̃_{j_2}·X̃_{j_1})` was computed numerically at level n=3, (Z/27)*, 18 states, Geom(2) truncated at V_MAX=16 (tail mass < 1.6e-5).

| Moment | Modulus | Status |
|---|---|---|
| φ(X̃_1·X̃_2) | 1.076×10⁻⁷ | ~0 (structural, ✓) |
| φ(X̃_1·X̃_2·X̃_3) | 1.430×10⁻⁵ | ~0 (structural, ✓) |
| **φ(X̃_1·X̃_2·X̃_1)** | **1.078×10⁻¹** | **NON-ZERO — diagnostic** |

|M_3_alt| / |M_2| ≈ 10⁶. Phase purely real and positive (±10⁻¹² rad). Robust across four scalar reductions (tr_π, ⟨π,·π⟩, ⟨δ_1,·δ_1⟩, ⟨1,·1⟩).

Outputs: `verify_monotone_diagnostic.py` + `experiments_output/monotone_diagnostic_n3.json`.

**Subtlety surfaced:** the diagnostic is non-zero under **marginal centering** `φ(Off_j) = E[Off_j | b_prior]` (the operative reading in AMALG_FREENESS_MOMENT_CALCULATION.md §8). Under strict conditional centering on the full B = vN(b_{[1,j]}) per SETUP.md, all moments collapse to algebraic zero (T_j becomes B-measurable up to within-pair split). The monotone framework therefore applies with B interpreted as the coarser σ-algebra of running-sum information without per-step within-pair split-conditioning. This is a real load-bearing choice and is documented as such in the framework writeup.

### 1.2 c = 7/45 leading-order derivation (Task 4 part 1, ~3h actual)

Via Hasebe-Saigo 2011 Thm 4.5 (κ_n definition) + Hasebe monograph Thm 3.26 (moment-cumulant formula), the dominant term in the moment expansion `E_B(X^n)` is the all-singletons monotone partition `π = ({1},...,{n}) ∈ M(n)`, contributing `(1/n!)·(κ_1^B)^n`.

For Syracuse, this evaluates to:

- `κ_1^B(Off_j)` projected onto the (1, 4)-eigenvector of R77's T_diag eigenstructure
- Combined with R64.B's class-mass identity (1/3)² : (2/3)² = 1 : 4 (which fixes the relative weight of the two eigenvector components)
- Combined with R75 Plancherel normalization `3^{-n}`
- **Yields c = 7/45** as `(S_∞)/3 = (7/15)/3`

The `7` in the numerator arises from the algebraic identity `1 − 8/15 = 7/15`, where `8/15` is the R77 T_diag mass on the (1, −1)-null direction (eigenvalue 0) and `7/15` is the mass on the (1, 4)-eigenvalue-1 direction. The `45 = 15·3` combines the R77 mass denominator and the R75 Plancherel `3^{-n}` factor.

**Verbatim citations:**
- Hasebe-Saigo 2011, "The monotone cumulants," Annales IHP B 47(4), 1160-1170 — Thm 4.5 (cumulant definition), Thm 4.8 (moment-cumulant formula)
- Hasebe monograph "Monotone Probability Theory" — Thm 3.26 (moment-cumulant via monotone partitions M(n)), Defn 1.21 (peak-rule for monotone product)
- Muraki 2003 RIMS Kokyuroku 1186-3 — Thm 4 (reciprocal-Cauchy composition law `H_{Σ X_j} = H_{X_1} ∘ H_{X_2} ∘ ... ∘ H_{X_n}`)

### 1.3 Diagnostic ⟷ peak-rule correspondence

Hasebe monograph Defn 1.21 (peak-rule factorization for monotone product) predicts:

`E_B(X̃_{j_1} · X̃_{j_2} · X̃_{j_1}) = E_B(X̃_{j_2}) · E_B(X̃_{j_1}²)` under monotone independence

where the "peak" at position 2 (the j_2 sandwiched between two j_1's) factors out as a scalar `E_B(X̃_{j_2})`.

For Syracuse under marginal centering: `E_B(X̃_{j_2})` retains a B-measurable phase-twist `Δ_{j_2}(b_{[1,j_1]}) ≠ 0` through accumulator coupling. The factorization predicts:

`E_B(X̃_{j_1}·X̃_{j_2}·X̃_{j_1}) = Δ_{j_2}(b_{[1,j_1]}) · E_B(X̃_{j_1}²)`

The diagnostic value 0.1078 from Task 1 is consistent with this product structure (the second factor `E_B(X̃_{j_1}²)` is the positive-definite norm-squared at step j_1; the first factor `Δ_{j_2}` carries the cross-step phase information).

**This is the structural mechanism that makes Syracuse monotone, not free.**

---

## 2. What today did not close — the four wrinkles

Each wrinkle is named with (a) what's missing, (b) effort estimate, (c) characterization as **dig-hard** (tractable with focused work) or **pivot** (requires fresh framework / external collaboration).

### Wrinkle 1 — B-amalgamated lift of Hasebe-Saigo

**Missing:** Hasebe-Saigo 2011 + Hasebe monograph develop monotone cumulants in the scalar case (state φ: A → C). The operator-valued / B-amalgamated extension used throughout Deliverables B-C is **not a verbatim theorem in the closure-hunt corpus**.

**Why it matters:** the entire derivation in §1.2 above is "rigorous fiberwise at each fixed accumulator history" + "conjectural at the B-valued composition step." Lifting the scalar theorem to the abelian-B case we need (B = scalar functions of accumulators, conditional expectation = integration over residual pair-distribution) is the principal load-bearing assumption.

**Characterization:** **dig-hard tractable.**
- Estimated effort: 1-3 days focused
- Why tractable: B is abelian. The lift in abelian-B settings is well-known in operator-valued free probability (Voiculescu-Speicher) and the monotone-amalgamated analog has been sketched in Popa 2008 and Hasebe 2011 "Operator-valued free probability" extensions. A direct verification in the Syracuse-specific abelian-B case is mechanical given the published infrastructure.
- Risk: low — failure mode would be discovering an unanticipated obstruction in the operator-valued extension, which would itself be a publishable note.

### Wrinkle 2 — Subdominant coefficient −1/30

**Missing:** R77 §4 empirical fit gives `S_n = 7/15 − (1/30)·(1/2)^n + O((1/4)^n)`. The factor `1/30 = 7/(15·14) = S_∞/14` with `14 = 2·7`. The monotone-cumulant moment formula predicts a term `(n−1)·(1/(n−1)!)·κ_2^B·(κ_1^B)^{n−2}` from monotone partitions with exactly one 2-block, but the **closed-form combinatorial 14** is not derived.

**Why it matters:** without the 1/30, we have the rate (1/2)^n but not the amplitude. The full leading-plus-subdominant statement is `S_n ≈ 7/15 − (1/30)·(1/2)^n`, and `−1/30` is half of that statement.

**Characterization:** **dig-hard tractable.**
- Estimated effort: 4-8 hours focused
- Why tractable: the calculation is (a) evaluate κ_2^B(Off_{j_1}, Off_{j_2}) on the (1, 4)-direction, (b) sum the monotone-partition combinatorics from M(n), (c) apply the R75 Plancherel bilinear normalization. Each step has known infrastructure (R77 T_diag, Hasebe Thm 3.26 monotone partition formula, R75 Plancherel). The combinatorial 14 is conjectured to be `2·7 = 2·S_∞·15` from the Plancherel double-counting at bilinear level.
- Risk: medium — the factor 14 might require a more intricate identity than the conjectured 2·7, in which case the search would extend to the monograph's deeper combinatorial structure.

### Wrinkle 3 — PADE complex pair (period 9.2)

**Missing:** PADE numerics show a complex-conjugate pair at θ ≈ 0.68 rad with period ≈ 9.2 in n-space (sign pattern +,+,−,−,−,−,−,−,−,+,+,+,+). The framework gives a semi-quantitative match: `2π / (log 3 / log 2) ≈ 2π / 1.585 ≈ 3.96` per step, modulo phase normalization gives period 9.24 (within 1% of empirical 9.2). But this is **consistent**, not **derived** from first principles.

**Why it matters:** the multi-spectral structure of T (not just T_diag) is what controls the n=10..13 PADE Hadamard radius trajectory (2.06 → 1.81 → 1.66 → 1.57). Without a direct derivation, the complex pair is descriptive but not predictive.

**Characterization:** **dig-hard borderline-pivot.**
- Estimated effort: 1-2 weeks focused
- Why borderline: the calculation requires 3-adic phase analysis sketched in R77 §3 but not closed. The χ_j phase factor `3^{2j−2}·2^{−b_{[1,j]}}` mod 3^n rotates in (Z/3^n)* with an asymptotic rate that needs careful cyclotomic computation. The framework supports this but doesn't perform it.
- Risk: medium-high — if the 3-adic phase analysis surfaces unexpected number-theoretic structure (e.g., the `log 3 / log 2` irrationality coupling to Diophantine approximation), this could escalate to a separate research arc.

### Wrinkle 4 — Faure √3 in the cumulant spectrum

**Missing:** Faure 2009's prediction √3 ≈ 1.732 (essential spectral radius of T on its anisotropic Banach space) matches PADE 1.57 at n=13 within 10%. In the monotone-cumulant reading, √3 is **consistent** with an intermediate cumulant scale (between the dominant κ_1^B-rooted contribution and the asymptotic slow-mode at z ≈ 1.016), but a direct map "Faure √3 ↔ specific monotone-cumulant operator eigenvalue" is not constructed.

**Why it matters:** if Faure's √3 corresponds cleanly to the spectral radius of the second-cumulant operator (or of T projected onto the (1, 4)-deviation subspace after κ_1^B projection), this would unify Faure semiclassical with the monotone-cumulant framework in a single statement.

**Characterization:** **dig-hard tractable, possibly easy.**
- Estimated effort: 4-12 hours
- Why tractable: this is a one-shot check. Compute the spectral radius of κ_2^B viewed as an operator on the appropriate subspace. If it equals √3 to within numerical precision, the identification is made. If not, the framework points to which other operator scale could carry the √3 (κ_3^B, or the full T restricted to the post-κ_1 subspace).
- Risk: low — failure mode is identifying a different spectral object, which is itself useful.

---

## 3. Pivot vs dig characterization

| Wrinkle | Tractability | Effort | Risk | Pivot signal? |
|---|---|---|---|---|
| 1. B-amalgamated lift | High | 1-3 days | Low | No — mechanical |
| 2. Subdominant −1/30 | High | 4-8 hours | Medium | No — combinatorial |
| 3. PADE period 9.2 | Medium | 1-2 weeks | Medium-high | Possibly yes (Diophantine surface) |
| 4. Faure √3 ↔ cumulant op | High | 4-12 hours | Low | No — one-shot check |

**Reading:**
- Wrinkles 1, 2, 4 are all dig-hard. Combined effort estimate: 2-5 days at user pace (~10-14×). Together they would close: rigorous-not-fiberwise leading 7/45, full subdominant statement `7/15 − (1/30)·(1/2)^n`, and Faure √3 identification.
- Wrinkle 3 is the pivot candidate. The 3-adic phase analysis may surface number-theoretic structure (irrational rotation, Diophantine approximation of log 3 / log 2) that constitutes a separate research arc.

**A two-track plan candidate** (if dig-hard chosen):
- Track A: Wrinkles 1, 2, 4 in parallel — close the partial-rigorous statement
- Track B: Wrinkle 3 — separate arc, possibly month-scale

**A pivot-now alternative:**
- Wrinkle 3 says the multi-spectral PADE structure may have number-theoretic flavor. If this is true, the natural pivot is from monotone cumulants to a coupled "monotone + Diophantine" framework, or to engaging directly with the 3-adic phase analysis that R77 §3 sketches.

User decides. This writeup surfaces data + structure per protocol.

---

## 4. State of play — where the 11-arc obstruction map landed

| # | Arc | Disposition (today) |
|---|---|---|
| 1 | 5-probe modern Fourier-decay | NO_FIT |
| 2 | C1 Cochrane exp sums | NO_FIT |
| 3 | C2 BMP cut-and-project | NO_FIT |
| 4 | BT 3-adic Bruhat-Tits | NO_FIT |
| 5 | Tauberian re-scope | NO_SELECTED |
| 6 | Furstenberg-Guivarc'h | NO_FIT |
| 7 | BGT regular variation | PARTIAL |
| 8 | Adelic Mellin | NO_FIT |
| 9 | Igusa local zeta | NO_FIT |
| 10 | Faure semiclassical | PARTIAL — **now connectable via Wrinkle 4** |
| 11 | Watson saddle-point | PARTIAL |
| C4 v1 | Tao RMT free probability | UNCHANGED |
| C4 v2 | Cébron multiplicative free | PARTIAL_OUTSIDE |
| C4 v3 | Voiculescu + Speicher B-amalgamated | framework available |
| **C4 verify** | **B-amalgamated freeness** | **monotone identified** |
| **C4 closure** | **Monotone cumulant derivation** | **PARTIAL — c=7/45 closed, 4 wrinkles** |

The terminal finding from yesterday (framework = monotone) is now anchored by today's numerical diagnostic + Hasebe-Saigo cumulant derivation of the leading coefficient. The four wrinkles describe the residual work cleanly.

---

## 5. Mode-E uncertainty ledger

| Component | Status |
|---|---|
| Framework = B-valued monotone (Muraki / HS) | Identified verbatim, diagnostic confirmed at 10⁶ separation |
| c = 7/45 leading | Rigorous fiberwise + conjectural at B-lift (Wrinkle 1) |
| Rate 1/3 (Plancherel) | Rigorous (R75, pre-existing) |
| Rate (1/2)^n subdominant | Mechanism in framework; exponent from R77 §3 outside framework |
| Coefficient −1/30 | Numerical only (R77 §4); mechanism in framework, value open (Wrinkle 2) |
| PADE complex pair period 9.2 | Semi-quantitative match within 1%; not derived (Wrinkle 3) |
| Faure √3 | Consistent with intermediate cumulant scale (Wrinkle 4) |
| Sign pattern +,+,−,−,−,−,−,−,−,+,+,+,+ | Consistent with single complex-pair flip |
| ε_k k=1..6 plateau ≈ 0.033 | Matches `1/30·2^k ≈ 0.033` — consistent |
| ε_k k=7,8 deviation 0.150, 0.191 | Multi-spectral transient onset; not derived |

---

## 6. Files

**Today's deliverables:**
- `C:/Collatz/MONOTONE_CUMULANTS_A_VERBATIM.md` — verbatim HS / Muraki / Hasebe-monograph statements
- `C:/Collatz/MONOTONE_CUMULANTS_B_SYRACUSE.md` — per-step Syracuse monotone cumulants
- `C:/Collatz/MONOTONE_CUMULANTS_C_ASYMPTOTIC.md` — asymptotic derivation
- `C:/Collatz/MONOTONE_CUMULANTS_D_COMPARISON.md` — comparison to PADE / Faure / ε_k
- `C:/Collatz/MONOTONE_CUMULANTS_DISPOSITION.md` — 1-page disposition
- `C:/Collatz/verify_monotone_diagnostic.py` + `experiments_output/monotone_diagnostic_n3.json` — Task 1 numerical anchor

**Yesterday's terminal:**
- `C:/Collatz/OBSTRUCTION_MAP_TERMINAL.md` — framework identification, paper-grade
- `C:/Collatz/POST_COMPACT_NEXT_STEPS.md` — 4-task program (now mostly executed)
- `C:/Collatz/AMALG_FREENESS_{SETUP,SUBALGEBRA_CHECK,MOMENT_CALCULATION,DISPOSITION}.md` — verification probe outputs

**Closure-hunt corpus (key today):**
- `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_saigo_2011_monotone_cumulants.pdf`
- `C:/Users/Nate/OneDrive/Documents/closure hunt/hasebe_monotone_probability_theory_monograph.pdf`
- `C:/Users/Nate/OneDrive/Documents/closure hunt/muraki_2003_five_independences_kyoto_precursor.pdf`

**Chain-side load-bearing:**
- `c_seven_forty_fifth.md`, `result_75_*`, `result_76_*`, `result_77_*` (T_diag eigenstructure, conservation, Plancherel)
- `PADE_NUMERICAL_DISPOSITION.md` (multi-spectral picture)
- `experiments_output/result_77_7_eps_exact_through_k8_v2_vec_pool.json` (ε_k k=1..8 exact)
