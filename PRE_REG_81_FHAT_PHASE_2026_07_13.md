# PROBE 81: F̂ Phase Profile on the R78 Support (PRE_REG)

**Agent:** compute (Agent 1/2/3)
**Repo:** `C:\Collatz`
**Date locked:** 2026-07-13 (base pre-reg carried over from desktop draft; amendments dated inline).
**Status:** pre-registration — fire only after reading the source files in §1.
**Do not write to any file outside the new deliverables listed in §6.**

> **Amendment log (all PRE-fire, so still a pre-reg — not deviations):**
> - **§3′ (2026-07-13):** multi-denominator congruence sweep. Supersedes the single-3^r fit in §4.3.
> - **§4.6 exactness (2026-07-13):** test `F̂²/(9q)` is a root of unity (NOT `F̂/√q`). |F̂| = 3√q, not √q.
> - **DWM BRIDGE (2026-07-13):** §6 deliverable D5 + new §2d moment cross-check. This probe is now a candidate *proof of the Syracuse = DWM identification*, not only a route-selector. See §0.

---

> **▓▓ OUTCOME — FIRED 2026-07-13 (parallel agent). VERDICT: H_GROWING_DEGREE.**
> `arg F̂(3a)` is certified exactly (F̂(3a)²/(9q) a root of unity in ℤ[ζ_q]; σ=−1 for r even = the √−3 twist you predicted) and is **NOT a fixed-degree polynomial** — its degree in b GROWS with r (r=3→3, 4→4, 5→4, 6→5, ≈⌊r/2⌋+2), 3-adic-analytic. **H_QUAD, H_LIN, H_POLY_HIGHER(fixed), and H_PSEUDO(random) are ALL excluded** — the obstruction is growing degree, not randomness. Derived+verified identity: `F̂(3a) = 3·Σ_j e_q(c·4^j)·e_d(−aj) = 3·ĝ(a)`, the DFT of the **exponential chirp** e_q(c·4^j) (flat magnitude, non-quadratic phase — why the Gauss-sum picture can't hold). Closes the smooth-completion / stationary-phase route as a *uniform* √-saving mechanism (Weyl completion would need Θ(r) differencing). Thms 78.1–78.3 untouched. Full disposition: `result_81_fhat_phase.md`; STATE.md 2026-07-13 R81 entry.
>
> **⚠ RETRACTION (2026-07-13):** the two "H_PSEUDO ⇒ the DWM match is a coincidence" claims below (§3 H_PSEUDO bullet, §7) are **WRONG and retracted.** Probe 81 tested the phase in the *frequency* variable `a` (the R78 bilinear pairing `Σ 1̂(3a)F̂(3a)`). The DWM bridge is a **different pairing** — the Geom(½)-weighted sum over the *orbit step* `v`, `Σ_v 2^{−v} χ(2^{−v})`. H_GROWING_DEGREE does NOT settle it; it **sharpens** it: the chirp e_q(c·4^j) is now explicit, so Phase 2a is concrete — flat DFT (F̂) vs Geom-weighted sum (DWM) of the SAME explicit chirp. **1B (weighted PGF at z=−1, ω) is now the live central test.**

---

## 0. Why this probe matters more than "which of three routes" (2026-07-13)

R78's F̂ and the DWM adaptive Kraus operator are the **same character on the same 2-adic orbit, against two different measures.** This is the corrected bridge (2026-07-13) — the earlier "same object, phase = phase" framing was WRONG and would have read a true bridge as a refutation. See the CAUTION box below.

- **R78 F̂(ξ) = Σ_u e_{3^{r+1}}(c·4^u − ξu):** the Fourier transform of the **counting measure** on the orbit ⟨4⟩ (every orbit point weight 1 — a *flat* sum). Magnitude flat (|F̂| = 3√q on all q/9 support points, Thm 78.3), all structure in `arg F̂`. The generator is 4 = 2², ord(4 mod 3^{r+1}) = 3^r (index 2 — the principal units), which is exactly why supp = {a ≡ 1 mod 3}. That is forced, not coincidental.
- **DWM Kraus `M_v^{(j,b)} = 2^{-v/2} · exp(−2πi · 3^{2j-2} · 2^{−b−v} · phase / 3^n) · σ_{-v}`** (`FRAMEWORK_IDENTIFICATION.md` line 45): the **same character on the same orbit**, but evaluated along a **Geom(½)-weighted walk** — weight `2^{−v}` per step, not weight 1. Ergodic invariant is (1,4) (line 53) — the same squared-class-mass 9 that sits in 7/45 = 7/(9·5).

> **CAUTION — do NOT test phase = phase (2026-07-13).** Same group, same generator ⟨4⟩, same 3-power character, **different measure on the orbit** (flat vs geometric). So `arg F̂ = DWM phase` will NOT come out as an equality even when the identification is TRUE. The correct bridge (§2a/§2d): F̂ is the FT of the counting measure; the DWM moments are the same character **integrated against the Geom(½) weight**. They differ by exactly the weighting. Probe 81 supplies the character (`arg F̂`); Phase 2a tests whether the DWM moments are recoverable as the weighted orbit sum `Σ_v 2^{−v} · χ(2^{−v})` with χ read off F̂. That is still a proof of the DWM identification if it lands — it is just not phase = phase.

The DWM identification is currently NUMERICAL ONLY (6 sig digits, `DWM_MP_G1_RESULT.md`, never derived). The load-bearing new stake: **if the character read off `arg F̂`, integrated against the Geom(½) weight, reproduces the four DWM reductions, the 6-sig-digit coincidence upgrades to an algebraic identity.** See §2d, §CONJ, and deliverable D5.

### §CONJ — the weighting is where 7, 9, and 5 come from (CONJECTURE, 2026-07-13)

The measure that separates F̂ from DWM is the Geom(½) halving weight. Its probability generating function is the whole halving process in one object:

    G(z) = Σ_{v≥1} 2^{−v} z^v = z / (2 − z)

Evaluate it at roots of unity — precisely "2-adic orbit read through a 3-power character":

- **z = −1 (order 2):** `G(−1) = −1/3` ⇒ P(v even) = (1+G(−1))/2 = **1/3**, P(v odd) = 2/3 ⇒ class mass (0, 1/3, 2/3) ⇒ squared 1 : 4 ⇒ **the 9 in 7/45**.
- **z = ω (primitive cube root):** `G(ω) = ω/(2−ω)`, so `|G(ω)|² = 1/|2−ω|² = 1/N(2−ω)`. In the Eisenstein integers ℤ[ω], `N(2−ω) = 4 + 2 + 1 = 7` ⇒ `|G(ω)|² = 1/7` ⇒ **the 7 in 7/45**, as the norm of the halving weight against the cube-root character.
- **the 5:** `T_diag = (1/5)[[1,1],[4,4]]`; the 1/5 normalizes the (1,4) eigendirection, `1 + 4 = 5`.

**Conjecture (not result):**

    7/45  =  N(2 − ω) / ( 3² · (1 + 4) )
          =  7 / ( 9 · 5 )

numerator = norm of the halving weight against the cube-root character; denominator = squared class mass (3²) × the (1,4) eigendirection normalization (5).

**Where the teeth are.** The arithmetic `7/45 = 7/(9·5)` is just factoring. The content is whether the 7 that R75/R76 pull out of the Plancherel sum **is** `N(2−ω) = |G(ω)|^{−2}` mechanistically — i.e. whether the halving-weight-against-cube-root-character is the origin, not a numerical rhyme. That is a Phase-2 **derivation** target, NOT settled by Probe 81 or 1B. Probe 81 supplies the character; 1B (re-scoped) tests robustness of the inputs −1/3 and 1/7; the norm-mechanism is the separate prize. This subsection collapses the old Phase 2c ("pin the 5") — the 5 is identified.

---

## 1. Context (read first, do not re-derive)

- `result_78_FINAL.md` — Theorems 78.1 / 78.2 / 78.3.
- `result_78d_fourier_sparsity.py` — existing code that already builds F̂.
- `FRAMEWORK_IDENTIFICATION.md` — DWM Kraus operator (line 45), (1,4) invariant (line 53), moment table (§Moment-pattern fit).
- `DWM_MP_G1_RESULT.md` — the four verified scalar reductions (target numbers for §2d).
- `band_l1_analysis.md` — band-ℓ1 route, CLOSED. Do not re-run it.
- `bk_moments_analysis.md` — C2/BGK additive-energy route, random-like. Do not re-run it.

Established (rigorous, do not re-verify beyond a smoke check):

- **78.1** For all r ≥ 2, ℓ ∈ {0,1,2}, ε ∈ {0,1}, m ∈ Z: the complete sum
  `Σ_{u=0}^{3^{r+1}-1} e_{3^{r+1}}(c_{ℓ,ε}·4^u − 9mu) = 0`.
- **78.2** `F̂(ξ) = Σ_{u=0}^{3^{r+1}-1} e_q(c·4^u − ξu)` with `q = 3^{r+1}` is supported on
  `supp(F̂) = {3a : a ∈ Z/3^r, a ≡ 1 mod 3}`, `|supp| = 3^{r-1} = q/9`.
- **78.3** `|F̂(ξ)| = 3√q` for every ξ ∈ supp(F̂). **Constant magnitude.**

**The gap this probe attacks.** R78 §"Crucial observation" onward bounds `Σ 1̂(3a)·F̂(3a)` by Cauchy–Schwarz on the two factors separately and recovers only the trivial bound. Its own stated conclusion: the square-root saving must come from **phase cancellation in the product**, not from either factor. Since |F̂| is constant on the support, **all** of F̂'s structure lives in `arg F̂(3a)`. That phase profile has never been mapped. This probe maps it.

---

## 2. Object

For each level r, define on the support index `a ∈ A_r := {a ∈ Z/3^r : a ≡ 1 mod 3}`:

    θ_r(a) := arg F̂(3a)          (principal branch, unwrapped along a)

Normalize to units of the character modulus (primary denominator per §3′):

    φ_r(a) := θ_r(a) · D / (2π)   (real-valued, tested mod D for D ∈ {3^r, 3^{r+1}, 2·3^r})

### 2d. Character extraction for the weighted-sum bridge (2026-07-13, corrected)

**NOT phase = phase.** Probe 81 extracts the **character** χ from `arg F̂`, so that Phase 2a can test whether the DWM moments are recoverable as the **Geom(½)-weighted orbit sum**

    DWM moment  ?=  Σ_{v≥1} 2^{−v} · χ(2^{−v})

with χ supplied by F̂. Record φ_r(a) as a function of **the 2-adic orbit variable** (the exponent u, equiv. 4^u = 2^{2u}), NOT of the raw support index a — this exposes χ as a function of the orbit point, which is what the weighted sum needs. Tabulate (u, χ(orbit point u) := exp(iθ_r(a(u)))) under the identification character-modulus `3^{r+1} ↔ 3^n`. This is the table Phase 2a feeds into the weighted reconstruction. Do **not** tabulate against `Φ_DWM = 3^{2j-2}/3^n` for an equality check — that comparison is the retired phase=phase framing.

---

## 3. Hypotheses (pre-registered, mutually exclusive)

- **H_QUAD** — `φ_r(a) ≡ α·a² + β·a + γ (mod D)` with α ≢ 0 mod 3.
  F̂ is a genuine Gauss sum on the principal-unit group. **Two payoffs:** (i) opens a stationary-phase / Weyl-differencing route on `Σ 1̂(3a)F̂(3a)` not in the paper's current route list; (ii) a quadratic `arg F̂` fixes **the character** χ on the orbit — the DWM proof is then whether the Geom(½)-weighted sum `Σ_v 2^{−v} χ(2^{−v})` lands on the four DWM reductions (§2d, §2a), NOT coefficient equality with Φ_DWM.
- **H_LIN** — same with α ≡ 0 mod 3, β ≢ 0. The product sum collapses to a shifted short sum; direct evaluation possible. Prior: low.
- **H_POLY_HIGHER** — best fit needs degree ≥ 3 in a. Routes to Weyl-differencing at higher order.
- **H_PSEUDO** — no low-degree polynomial fit at ANY denominator in §3′; residuals consistent with equidistribution on Z/D. **The Burgess wall is genuine.** Certifying negative — retires the smooth-completion route cleanly and is publishable as such. ~~(Also implies the DWM match is a numerical coincidence, not a phase identity.)~~ **[RETRACTED post-fire — see OUTCOME banner: the DWM bridge is a different pairing (weighted sum over v), not settled by this probe. Actual verdict was H_GROWING_DEGREE, not H_PSEUDO.]**

**Most-likely outcome (state it, then test it honestly):** H_QUAD, on the grounds that constant magnitude on a coset of the principal units is the fingerprint of a Gauss sum, AND the DWM Kraus phase is manifestly quadratic-shaped in the orbit exponent. Do not let this prior bias the fit acceptance thresholds in §4.

### 3′. Denominator sweep (AMENDMENT 2026-07-13 — replaces single-3^r fit)

The natural Gauss-sum phase may sit at denominator `3^{r+1}` (the character modulus F̂ actually lives on) or `2·3^r` (full group order incl. the b_prior parity), not `3^r` (the support index set). Run the congruence fit of §4.3 at **all three** denominators:

    D ∈ {3^r, 3^{r+1}, 2·3^r}

plus a **fourth case**: allow (α, β, γ) rational with small denominator (catches a completing-the-square half-shift; 2 is a unit mod 3^r so it stays integral, but a character twist can force a 1/2 — the 2·3^r column should also catch this).

- **Primary denominator: 3^{r+1}** (highest weight — F̂ is defined on Z/3^{r+1}; the base draft's 3^r came from the support *index* a ∈ Z/3^r, which is the index set, not the character modulus).
- Report a pass/fail table across all four cases × all r.
- **H_PSEUDO fires only if all four fail at every r ≥ 3.**

---

## 4. Method

For r ∈ {2, 3, 4, 5} and (extend if cheap) r = 6, and for every `(c_{ℓ,ε})` combination used in `result_78d_fourier_sparsity.py`:

1. **Smoke check.** Recompute |F̂(3a)| across the support. Assert constancy to 1e-12 relative. If this fails, STOP and report — 78.3 is wrong and nothing downstream is valid.
2. Compute `θ_r(a)` for all a ∈ A_r. Use exact-arithmetic phase where feasible (`sympy`/`Fraction` on the exponent sums); float64 `np.angle` is acceptable as a first pass but record the precision used.
3. **Fit over Z/D, not over R** (per §3′). The right test is a *congruence* fit, not least-squares — the phase is an integer (or small rational) mod D. Solve for (α, β, γ) by exact linear algebra on three support points, then **verify on all remaining |A_r| − 3 points**. Pass/fail, not R².
   - Report: number of support points satisfying the congruence exactly / total, per denominator D.
   - **Decision rule:** H_QUAD FIRES at denominator D iff the fit from any 3 points verifies on **100%** of the remaining points, at **every** tested r ≥ 3. Anything less than 100% at any r ≥ 3 refutes the polynomial model at that D; go to §3′ next D, then step 4.
4. If all denominators in §3′ fail at degree ≤ 4: test H_PSEUDO. Compute the distribution of `φ_r(a) mod D` and the pair-correlation / discrepancy against uniform. Report a discrepancy statistic and whether it is consistent with equidistribution at each r, at each D.
5. **Cross-r coherence.** If a polynomial fires, report how (α, β, γ) depend on r. A clean r-recursion is a strong second signal; an incoherent one is a red flag that step 3 overfitted at small |A_r|.
6. **Exactness certification (AMENDMENT 2026-07-13 — corrects the base draft).** If a polynomial fires at any denominator, certify the Gauss-sum claim exactly:
   - Test whether **`F̂(3a)² / (9q)`** is a root of unity in `Q(ζ_{3^{r+1}})`. Since |F̂|² = 9q **exactly**, the denominator is a rational integer and the check is fully exact — no irrational normalization.
   - **Do NOT test `F̂/√q`** (fails by a factor of 3 since |F̂| = 3√q, not √q) and **do NOT test `F̂/(3√q)`** directly — squaring via F̂²/(9q) also sidesteps the √(−3) quadratic-Gauss twist a genuine char-3 quadratic Gauss sum carries. F̂ is a sum of roots of unity → an algebraic integer in `Z[ζ_{3^{r+1}}]` for free; the check is cheap and upgrades H_QUAD from "float fit at r ≤ 6" to a theorem.

**Bug guards (these have bitten this project before):**
- `|A_2| = 3` — at r=2 the quadratic fit has zero degrees of freedom and fits *by construction*. **r=2 carries no evidential weight for H_QUAD.** Require r ≥ 3 for any fire; state this explicitly in the write-up.
- Do not apply a magnitude threshold filter to the support (the R3 j-saturation bug of 2026-05-17 was exactly this).
- Branch/unwrapping choice is load-bearing. Do the fit as a congruence mod D (step 3), which is branch-free; use the unwrapped float phase only for plotting.

---

## 2a / 2d routing (Phase 2 consumers — for the conductor, not this run)

- **2a (weighted-sum bridge — NOT equality):** take the character χ extracted in D5 and test whether the DWM moments are recoverable as the Geom(½)-weighted orbit sum `Σ_{v≥1} 2^{−v} χ(2^{−v})`. F̂ supplies χ (counting measure, flat); the weight 2^{−v} supplies the measure DWM carries. Landing ⇒ DWM identification proved. Do NOT test phase = phase against Φ_DWM (retired 2026-07-13).
- **2d (moment cross-check, free triangulation):** feed the weighted reconstruction through the existing DWM channel construction (`dwm_kraus_match_syracuse.py` for 3-alternating, `dwm_kraus_match_g2.py` for 4-alternating) and check it reproduces the four already-verified reductions:
  - (d) 3-alt `ϕ(X̃_{j1}X̃_{j2}X̃_{j1})` = **+1.078308×10⁻¹** (sum_entries).
  - (f) 4-alt `ϕ(X̃_{j1}X̃_{j2}X̃_{j1}X̃_{j2})` = **+6.088793×10⁻¹** (sum_entries), **+5.357225×10⁻²** (tr_π), **+5.742026×10⁻²** (delta_1), **+4.775479×10⁻³** (vac_π).
  - Acceptance: match to the precision those numbers are already known (~6 sig digits). If the phase profile is the DWM phase, this MUST reproduce — no new heavy compute, the target numbers already exist.

*These two are Phase 2 work; do NOT run them inside Probe 81. This probe only characterizes the phase and emits the comparable table (D5). Routing comes after.*

---

## 5. What NOT to do

- Do not attempt the bilinear bound `Σ 1̂(3a)F̂(3a)` in this probe. This probe **only characterizes the phase.** Routing comes after.
- Do not run the Phase 2a/2d consumers here (§2a/2d above are for the conductor).
- Do not touch `lagarias_framework_synthesis.docx` or any external-facing file.
- Do not re-open band-ℓ1 or BGK.
- Do not refit anything on r=2 alone.

---

## 6. Deliverables

- `result_81_fhat_phase_profile.py` — the probe.
- `result_81_fhat_phase_data.csv` — one row per (r, ℓ, ε, a): a, orbit-exponent u, |F̂|, θ_r(a), φ_r(a) at each D ∈ {3^r, 3^{r+1}, 2·3^r}, congruence residual per D.
- **D5 (2026-07-13): `result_81_fhat_phase_character.csv`** — the §2d character table: (u, orbit point 4^u mod 3^{r+1}, χ(u) := exp(iθ_r(a(u)))) at denominator 3^{r+1}. This is the character F̂ supplies; Phase 2a feeds it into the Geom(½)-weighted sum `Σ_v 2^{−v} χ(2^{−v})`. NOT a Φ_DWM equality table.
- `result_81_fhat_phase.md` — disposition. Structure: hypothesis fired / refuted per §3 **at each denominator (§3′ table)**, the pass/fail table from step 3 with r=2 explicitly excluded from evidence, the (α, β, γ) cross-r table if applicable, the exactness certification result (§4.6) if a polynomial fired, the discrepancy statistic if H_PSEUDO, a one-paragraph routing statement (which of the paper routes this opens or closes), and a one-paragraph **DWM-bridge statement** (is the extracted character χ clean enough to feed the Geom(½)-weighted sum `Σ_v 2^{−v} χ(2^{−v})` — the actual bridge, pending Phase 2a — noting explicitly that phase = phase is NOT the test).
- `result_81_log.txt` — run log incl. precision used and denominator sweep.

Append the disposition to `STATE.md` under a new dated entry. **Do not rewrite existing STATE.md content.**

---

## 7. Reporting discipline

Report the outcome that fired, including a null. A clean H_PSEUDO is a *result*, not a failure — it certifies the Burgess wall and retires a route. ~~AND settles the DWM match as coincidence-not-identity.~~ **[RETRACTED — see OUTCOME banner; the DWM bridge is a separate weighted-sum pairing.]** Do not soften a refutation into a "partial." If the fit verifies on 97% of points at every denominator, that is a refutation of the polynomial model, not a near-miss. Do not let the H_QUAD prior (§3) relax the 100%-verification bar.
