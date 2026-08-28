# CORPUS LEDGER — the internal-relations map (Euler hunt, collect phase complete, 2026-08-28)

Built from a 12-agent read-only sweep of the whole `C:\Collatz` corpus (2,073 tracked files,
250 `result_*.md`, 555 probes, ~130 dispositions). Each agent produced a compact object-ledger +
already-stated relations for one thread. This file is the merge: the **map**, not every row.

⚠️ Values here are copied from the agent reports (which copied from the source files). Before the
disprove phase *gates* any figure, re-grep it from its own `result_*.md` ([[copy_dont_recall]]).

---

## 0. THE ONE-PARAGRAPH MAP

The corpus is saturated with **q=3 criticality** (the same phase boundary re-derived ~6 ways) and rests
on one **established cross-thread spine** (Mahler ↔ Galois-tower ↔ denominator, three derivations) that
*is* the **structure/value split** (unipotent/associated-graded/finite-place = law-governed & symmetric;
reductive/tail/archimedean = where the value lives, walled). Nearly every *beautiful* internal identity in
the corpus lands on the **leading-mode rational `7/15`** — which is **dead as the value** (true
`S_∞ ≈ 0.475`, floor `2·T_20 = 0.473177`). The true value = leading mode **+ the infinite unipotent tail**.
The most Euler-shaped moves (sum=product for the spectrum; the cyclotomic-norm reading of the "7"; order
reciprocity) **were already tried and DISPOSED NEGATIVE**. What is **never linked** — the genuine open
targets — are the *cross-thread* pairings: **inverse-tree `λ` ↔ Plancherel `S_∞`**, and the
**archimedean/entropy side ↔ the finite-place denominators**.

---

## 1. q=3 CRITICALITY — the same boundary, ~6 independent derivations (a convergence, mostly by-construction)

All of these are the assertion "Collatz sits at the critical point of its own (q,p) family," each from a
different framework. They agree because they are the *same* boundary, so this is a **convergence, not an
Euler identity** — but it is the corpus's strongest "many faces, one object" fact.

- `D₂ = H₂/log q = log 3 / log q`, and `D₂ = 1 ⟺ q = 3`   (correlation dimension) — result_12, RESEARCH_ARC
- `r_q = 1 ⟺ λ₁ degenerate`   (transfer operator RPF) — result_25
- Siegel `α_H(0) = (1+q)/4 = 1 ⟺ σ_H = 1 ⟺ q = 3`; exact set-equality `{q:α_H(0)=1} = {3} = {q:r_q=1}` — result_41 (**the one exact analytic-NT bridge**)
- `(q,p)`-Hydra boundary `q(p−1) = p+1`; Collatz `(3,2)` **doubly** critical (`q/3=1` AND `3(p−1)/(p+1)=1`) — result_PHYDRA_FAMILY
- collision entropy `H₂(Geom½) = −log Σ 4^{−v} = log 3 = log q` at q=3 — QX1 writeup
- `(q−1)/2 = 1 ⟺ q = 3`   (leading-mode uniqueness of the conservation law) — result_6

Hub constant: **the "3" `= 1/Σ_v 4^{−v} = 1/E[2^{−v}]`** (participation ratio of the halving law).

---

## 2. THE ESTABLISHED SPINE — Mahler ↔ Galois-tower ↔ denominator (THREE derivations, genuinely cross-thread)

The corpus's one **real, proven, cross-thread** identity chain. Three frameworks, one statement:

> **infinite Mahler order = unbounded unipotent depth = `𝔾_m` emerges only in the inverse limit
> = doubly-exponential Mersenne denominators `2^{2·3^{r−1}}−1` (ratio → 3).**

- denominator side: `den S_r | 3^a·(2^{2·3^{r−1}}−1)`, `log₂ D_r ~ 2·3^{r−1}`, ratio → 3 — result_DENOM_OBSTRUCTION, result_MAHLER
- spectral side: "infinite Mahler order **is** the unbounded unipotent depth… ratio-3 denominator growth **is** the `𝔾_m` denominators" — result_MAHLER §47-49
- Galois side: "the tower `∏_i 3^i` does not converge to a coboundary — its non-gaugeable growth **is** MAHLER's denominator rate; `𝔾_m` is a pure pro-object phenomenon" — result_QDIFF2 §29, result_QDIFF3

**This chain IS the structure/value split** (5+ confirmations across frameworks):
unipotent `𝔾ₐ` / associated-graded / finite place = **the real, unconditional, symmetric content**;
reductive `𝔾_m` / tail / archimedean = **the contested piece where the value lives, delivered by no framework.**
Probabilistic restatement (free-prob thread): the **diagonal is free** (→ `7/15`); the **off-diagonal
cross-frequency coupling is NON-free (monotone)** — the non-freeness *is* the correction toward `0.475`.

---

## 3. LEADING-MODE IDENTITIES — all real, all landing on the DEAD 7/15 (the "nice" relations)

These are proven or near-proven and beautiful, but they resolve the **leading mode** `S_∞ → 7/15`,
`c = 7/45`, which `THEOREM_C_745` itself labels *"leading-order convergence"* and `BOUNDARY_THEOREM_MAP`
marks *conditional on an unobserved log₃-turnover near r≈27 ⟹ real `S_∞ ≈ 0.477`*. They are the **leading
term of the object we actually want.**

- `S_{k+1} = 3^{k+1}·‖d_{k+1}‖²`  (R74) ; `c = 7/45 = (1/3)·S_∞`  — c_seven_forty_fifth, THEOREM_C_745
- `‖R_k‖² = ‖d_{k+1}‖²` (literally one vector) ; `‖R_k‖²·3^k = S_{k+1}/3 → 7/45` — result_77_5
- Conservation Law `Σ_{j=0,1,2} M_{n+1}(η₀+j·3^n) = 0` (proven, ports to all odd q incl. composite) — result_76, result_6
- Leading-mode `S_{n+1} = −2·M_{n+1}(1+3^n) = −2·M_{n+1}(1+2·3^n)` (proven) — result_76 Thm 76.3
- `D_k = M_k(1) = q^k‖d_k‖²` — "one object, three names" — result_7
- `c̃_q = M_k(1)/(q/3)^k` (the qx+1 pillar-2 constant IS R76's normalized S_k) — result_7
- Gowers `U²` identity: `mean_k|γ_k(m)−1|² = Σ_{a≠0}|μ̂_k(a)|⁴`, `U² = 0.29754` level-invariant, between-class part **exactly 2/9** — result_MAXMODE2
- `7/45 = N(2−ω)/(3²·(1+4))`, `N(2−ω)=7` — **q=3 FACT, but see §5: does NOT generalize** — result_1b
- "one object, three coordinates: additive `γ` on ℤ₃, stratum `C̄`, multiplicative `Λ`" — BOUNDARY_THEOREM_MAP
- inverse-tree: `λ=(3+√21)/6`, `3λ²−3λ−1=0`; conditional `1/λ`, branching `(√21−3)/6`, valuation law `(1−1/λ)(1/λ)^v` — **three empirical faces of one algebraic λ** — result_BRANCH_BIAS, result_MOD2K_STATIONARY

Subleading, partly-open: `S_n = 7/15 − (1/30)(1/2)^n + …` with `1/30 = 1/(2·15) = S_∞/14` ("Why 14? Open");
**the `(1/2)^n` subdominant rate is itself DISPOSED-NEGATIVE** (`|ε_7|·2^7 = 0.150`, THEOREM_C_745 §4.1).

---

## 4. LOAD-BEARING CONSTANTS (the hub nodes for pairing)

| constant | exact form | threads it appears in |
|---|---|---|
| the "3" | `1/Σ_v 4^{−v} = 1/E[2^{−v}]` | dimension, rate, family, entropy, denominator |
| `μ` (drift) | `log(4/3) = 0.287682` | stopping-time, entropy (`= H₁−H₂`), base-fit (`K_h=3/μ`), Tao |
| `H₁, H₂` | Shannon `log 4`, collision `log 3` (of Geom½) | entropy/dimension |
| `K_h` | `3/log(4/3) = 10.4282` | base-fit slope, stopping-time |
| `⟨α_det⟩` | `log 6/log(4/3) = 6.228263` | base-fit (matched to `10^{-15}`) |
| `λ` (inverse tree) | `(3+√21)/6 = 1.263763`, `3λ²−3λ−1=0` | inverse-tree ONLY (unlinked to S_∞) |
| `S_∞` | `≈ 0.475`, floor `2·T_20 = 0.473177`; leading `7/15` DEAD | ladder / the target |
| `α` | `log₂3 = log₄9 = 1.58496` | two-walls (archimedean face) |
| `ord₃(2)=2`, `2²−1=3` | | denominator, reciprocity, Cramér-root |
| subdominant rates | inverse-limit `ρ=0.834`; SOLSTICE tail `0.87`; contraction `0.80`; conjecture `1/2` | (`0.834` resolved = inverse-limit rate, ∉ other spectra) |
| Chang dim | `log(φ)/log 2 = 0.6942`, `φ=(1+√5)/2` | multifractal (distinct √5, vs inverse-tree √21) |
| Faure radius | `1/√3 = 0.577` | analytic-NT semiclassical |

---

## 5. THE EULER-MOVE GRAVEYARD — already tried, DISPOSED NEGATIVE (do NOT re-run)

- **sum = product for the spectrum:** `Ĉ(m) = Π_j w(2⁻ʲm/N)`, `w=1/(5−4cos2πu)` → corr ≈ 0 at all J — result_PRODFORM. *(The literal "our π²/6" product form, killed. Survivor: `A(K)=C(K)=autocorr(ρ)≥0` via `ρ≥0`.)*
- **cyclotomic-norm reading of the 7:** `7 = N(2−ω)` true at q=3, but as a family law `Φ_p(2)=2^p−1` it's a 2-point (p=3,5) coincidence — breaks ~0.65× at p=7,11,13 — result_85 rung 2
- **order reciprocity** `ord_8(3) ↔ ord_3(2)`: "one object" TRUE, **coupling law REFUTED** (power-of-2 coincidence; `d_q` generically not pow-2) — result_44
- **mirror value symmetry:** `S_∞(2,3) ≈ 0.459 ≠ S_∞(3,2) ≈ 0.475`; involution is a symmetry of *type* not *value* — result_MIRROR
- **inverse `D` = f(forward `S`):** all of `D·S, D+S, D/S`, diagonal fail; `D_n(1)/S_1 = 1/3` is a "trivial coincidence" — duality_S_vs_D_verdict
- **`c̃_q`'s `2^{ord}−1` = denominator's `2^M−1`:** "a rhyme, not a unifying identity" — result_CTILDE_EXTEND
- **Esscher duality** `E[L⁻]=E_P*[L⁺*]`: falsified 1480σ — closed_form_findings
- **imported frameworks:** freeness (alternating triple ≠ 0), free-RMT (coupling *is* non-freeness), Belavkin (A→B, actually Davies–Wiseman–Milburn), Garcia–Young, Polymath8, Pascadi, Bourgain–Konyagin/Korobov/Burgess — all DISPOSED off-object
- **"independent" H-dim methods:** five values in a 0.07 window, no algebraic identity (~10% agreement) — independence_audit
- **walk-back flagged:** `K_k` in Ayyer–Singla/Diaconis–Graham group-walk lineage (framework_cohesion, May 6) FALSIFIED later (RESEARCH_ARC, "+1 breaks group-walk") — later disposition wins

---

## 6. THE OPEN TARGETS — cross-thread pairings NEVER asserted-equal (the actual remaining hunt)

The seed-set agent's closing finding: the `EULER_RELATIONS_HUNT.md` endpoints are *"candidate endpoints not
yet asserted-equal anywhere."* The corpus is dense with *within-spine* identities (all → `7/15`) but has
**never linked across these threads.** These are the leads:

- **T1 — `λ = (3+√21)/6` (inverse-tree)  ↔  `S_∞` (Plancherel).** ❌ **DISPOSED-NEGATIVE** (result_EULER_T1):
  both natural candidates (`2logλ=0.4682`, `λ−1/λ=3−2λ=0.4725`) fall **below the exact floor `2·T_20=0.473177`**;
  the value region has no low-complexity λ-form; opposite sides of the structure/value split; forward↔inverse
  already disposed (duality_S_vs_D). Cannot prove nonexistence, but no support and no mechanism.
- **T2 — `μ = H₁ − H₂ = log 4 − log 3`  (drift = Shannon−collision entropy gap).** ✅ TRUE / ❌ **DISPOSED as a
  discovery** (result_EULER_T2): holds iff `q=p+1` (family test), and factors into `μ=H₁−log q` (definitional)
  + `H₂=log q` (known criticality). Clean restatement, leading-layer only, moves no new information.
- **T3 — `ε_S ≈ log 4` (archimedean renewal residue)  ↔  finite-place denominators / `H₁`.** `ε_S = log 4 = H₁`
  is *unconfirmed* (50M-orbit indecisive). If real, it ties the archimedean wall's residue to the Shannon
  entropy — an across-the-two-places statement, the opposite of the "each obstruction is local to its place"
  caveat. High-risk, high-reward.
- **T4 — product form on the *tail*** ❌ **SPACE CLOSED** (result_EULER_T4): A EXCLUDED (S_n is a Plancherel
  sum not a dynamical trace tr(Lⁿ) → no orbit product; leading residue = dead 7/15; subdominant pole → q=3
  EP → continuous spectrum, no discrete factor survives L→∞). B VACUOUS (Λ_i denominators factor into the
  cyclotomic cofactors trivially; numerators carry FOREIGN primes {5,149,883,396022747,…} with no cascade
  product). C RE-RUN⊕EXCLUDED (μ̂ FE is a weighted sum → iteration = renewal, not product; Plancherel |·|²
  off-diagonal = 60% of total, breaks Σ→Π). MAHLER backstop: any product's value-extraction needs a closing
  FE ⟺ finite Mahler order ⟺ contradicts MAHLER. **No SURVIVED-A-KILL.**
- **T3 — `ε_S ≈ log4` ↔ denominators** — NOT YET RUN (deferred; `ε_S=log4` itself unconfirmed).

---

## 7. STATUS

Collect phase: **COMPLETE** (12/12 threads). Disprove phase: **T2, T1, T4 disposed** (result_EULER_T2/T1/T4).
- **T2** ✅true/❌not-a-discovery (definitional ⊕ known criticality).
- **T1** ❌not-supported (natural candidates below the exact floor; opposite sides of the split).
- **T4** ❌space-closed (no structured product; A EXCLUDED, B VACUOUS, C RE-RUN⊕EXCLUDED, MAHLER backstop).
- **T3** deferred (`ε_S=log4` unconfirmed).

**Net Euler-hunt finding:** the corpus's beautiful internal identities all resolve the leading-mode rational
`7/15`; the true value `S_∞≈0.475` lives in the tail, and **the tail admits no non-trivial internal relation
we could find** — no cross-thread algebraic link (T1), no product form (T4). The Euler move (`π²/6`) needs a
closing functional equation; MAHLER is exactly the theorem that no such closure exists here. The wall is
re-confirmed from the internal-relations angle. Remaining live thread: T3, and whether a *genuinely new*
object (not in the current ledger) could bridge — Wilson's call.
