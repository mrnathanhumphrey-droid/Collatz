# TAUBERIAN_RESCOPE_DISPOSITION

**Date:** 2026-05-13.
**Probe:** Tauberian arc re-scope — single-theorem selection for c=7/45 closure.
**Mode:** E (verbatim theorem hypotheses from PDF, no inheritance from prior project files).

---

## Headline

**No SELECTED.** No theorem in the 20-PDF corpus admits inputs (1)+(2)+(3)+(4) as hypothesis-satisfying instances under Mode E discipline. The Tauberian arc closes **BLOCKER-dominant** with one categorical NO_FIT (F). The blocker pattern is uniform: **the analytic-continuation property used as Tauberian hypothesis is essentially the closure target itself** (Mode H circular). No candidate's load-bearing parameter (Chevalier M, FS α, Wiener-Ikehara A) is realizable from inputs (1)-(4) without first proving the closure result.

---

## Summary table

| Candidate | Theorem | Pre-reg priors (S/P/N/B) | Realized | File |
|---|---|---|---|---|
| A | Flajolet-Sedgewick VI.1-VI.5 (singularity analysis) | 20 / 30 / 35 / 15 | **BLOCKER** | A_HYPOTHESES, A_HYPOTHESIS_CHECK |
| B | Chevalier 2507.15394 Theorem 1.16 | 15 / 35 / 30 / 20 | **BLOCKER** (+ M unrealizable) | B_HYPOTHESES, B_HYPOTHESIS_CHECK, B_M_PARAMETER |
| C | Korevaar 2002 Wiener-Ikehara (Thm 4.2) + Newman-Korevaar (Thm 6.1, 8.1) | 10 / 25 / 50 / 15 | **BLOCKER** (Mode H) | C_HYPOTHESES, C_HYPOTHESIS_CHECK |
| D | Newman 1980 / Zagier 1997 Analytic Theorem | 10 / 20 / 55 / 15 | **BLOCKER** (Mode H) | D_HYPOTHESES, D_HYPOTHESIS_CHECK |
| E | Alberts 2508.20814 Theorem 1.1 (twisted moments) | 5 / 20 / 60 / 15 | **BLOCKER** (Mode H) | E_HYPOTHESES, E_HYPOTHESIS_CHECK |
| F | Singha Roy 2511.15928 LSD Theorem 1.1 | 5 / 20 / 60 / 15 | **NO_FIT** (categorical) | F_HYPOTHESES, F_HYPOTHESIS_CHECK |
| G | Tao 2020 Notorious Collatz | 0 / 5 / 30 / 65 | **BLOCKER** (pointer-only) | G_HYPOTHESES |
| H1-H12 | Borwein, Holland, Ingham-EM, SD-remarks, Korevaar CV, Häggström, Wiener-Mandrekar, Riemenschneider, GuideTauberian, Lagarias×3 | 5 / 30 / 50 / 15 (aggregate) | **NO_FIT or BLOCKER** | H_HYPOTHESES |

(S = SELECTED, P = PARTIAL, N = NO_FIT, B = BLOCKER.)

Total candidates examined: **8 primary (A-G + H aggregate) + 12 secondary (H1-H12 individually scanned) = 20 PDFs in corpus, all consulted**. Zero new candidates surfaced during reading beyond what was in the brief's list.

---

## Detailed dispositions

### B (Chevalier 1.16) — BLOCKER

The pre-registered favorite. Theorem 1.16's conclusion b_n = D n^{M - 3/2} (1 + d_1/n + …) requires:
- A generating function g(z) = Σ b_n z^n analytic on D, continuous on D̄.
- An explicit meromorphic profile h_p on a neighborhood of D(1,1)^{1/2} with a single pole at 0 of multiplicity M ≥ 1, satisfying g(z) = h_p(√(1-z)).

Input (1) gives 8 numerical coefficients with the empirical pattern:
- |ε_k|·2^k ≈ 0.04 (roughly flat) for k = 2..6 — consistent with |ε_k| ~ C · 2^{-k}.
- |ε_7|·2^7 = 0.150, |ε_8|·2^8 = 0.191 — a **sharp jump at k = 7**, breaking the geometric pattern.

Backing out M from log-log fit on b_n = |ε_n|·2^n for n=2..6 gives **M − 3/2 ≈ 0**, i.e. **M ≈ 3/2** — but Theorem 1.16 requires M to be a *positive integer* ≥ 1. Backing out M from b_n = |ε_n| for n=3..6 gives slope ≈ −3.35 i.e. **M ≈ −1.85** — also incompatible with M ≥ 1.

The k=7 jump itself is the most informative observation: a clean single-pole-at-0 asymptotic n^{M-3/2} predicts monotone behavior for large n. Empirical b_n is non-monotone at the transition n=6→7. **No integer M ≥ 1 fits.** See `TAUBERIAN_RESCOPE_B_M_PARAMETER.md`.

Even ignoring (1)'s 8-coefficient limitation, hypotheses h_4 (∃ meromorphic h_p) and h_6 (g = h_p ∘ √(1-·)) require knowing g globally — not deliverable from any of inputs (1)-(4) without first proving the closure.

### A (Flajolet-Sedgewick VI.1-VI.5) — BLOCKER

Singularity analysis applied to f(z) = Σ ε_n z^n requires:
- Identification of dominant singularities ζ_j on |z| = ρ (radius of convergence).
- Δ-analyticity: analytic continuation to ζ_j · Δ_0 for each j.
- Singular expansion in standard scale S = {(1-z)^{-α} λ(z)^β}.

Inputs (1)-(4) supply: 8 coefficients, structural renewal-walk form, BMP F_1 diffraction (a different object), archimedean-place obstruction. **None gives the analytic continuation.**

Empirically the k=7 jump is consistent with a *multi-singularity* setup (Theorem VI.5, two dominant singularities each contributing a power-law) — a structurally interesting observation, but verifying VI.5's hypotheses requires the same analytic continuation A doesn't have.

### C (Wiener-Ikehara + Newman-Korevaar Dirichlet) — BLOCKER (Mode H circular)

Both Theorem 4.2 (W-I Laplace-Stieltjes) and Theorem 6.1/8.1 (Newman-Korevaar Dirichlet) require: g(z) = f(z) − A/(z-1) extends continuously/holomorphically to the boundary line Re z = 1. This is **the polynomial-in-A target itself dressed as a Dirichlet-series analytic-continuation property**. Mode H circular.

The nondecreasing-S(t) hypothesis (h_3 of 4.2) is satisfiable by construction (cumulative |ε_k|·2^k). But the continuous-extension condition is not.

### D (Newman-Zagier Analytic Theorem) — BLOCKER (Mode H circular)

Strictly weaker than C: requires g(z) extending holomorphically to Re z ≥ 0; conclusion is convergence of ∫ f(t) dt (no asymptotic). Same Mode H circularity.

### E (Alberts 2508.20814) — BLOCKER (Mode H circular)

Replaces pointwise vertical bounds (C, D) with twisted-moment bounds for better error terms. Underlying analytic-continuation hypothesis (h_4: L(s, F) meromorphic continuation to Re s ≥ σ_a − δ) is unchanged. Same Mode H trap.

### F (Singha Roy LSD) — NO_FIT (categorical)

Property P(ν, {α_χ}_χ; c_0, Ω) requires Dirichlet-character decomposition over a *fixed-integer modulus* (ℤ/qℤ)^* with explicit L-function exponents (α_χ)_{χ mod q}. Syracuse's structure on (ℤ/3^n)^* with *growing modulus 3^n* is a Markov-chain stationary distribution on a profinite group — categorically different from a Dirichlet-character decomposition mod a fixed q. **h_1 categorically fails.**

This is the *only* candidate that fails categorically rather than via Mode H circularity. The failure is informative: it indicates LSD-class machinery (Sathe-Selberg-Delange for arithmetic progressions) is not the right tool category.

### G (Tao 2020) — BLOCKER (pointer-only)

Pop-math slide deck. No Tauberian theorem statement extracted. Pointer to Lagarias bibliography (also covered in H10-H12 below).

### H1-H12 (secondaries) — uniformly NO_FIT or BLOCKER

- **H3 (Ingham-EM Bringmann et al.):** requires γ > 0 (exponential e^{c√N} growth of partial sums) — Syracuse ε_k has decay, not partition-type growth. **NO_FIT (wrong category).**
- **H4 (SD Remarks):** subsumed by F. **NO_FIT.**
- **H6 (Häggström textbook):** subsumed by C, D. **NO_FIT.**
- **H9 (Guide to Tauberian arithmetic apps):** Theorem A requires Dirichlet-series analytic continuation — same Mode H trap. **BLOCKER.**
- **H10-H12 (Lagarias surveys):** pointer-only. **BLOCKER.**
- **H1 Borwein, H2 Holland, H7 Mandrekar, H8 Riemenschneider:** all either subsumed by primary candidates or pointer-only. **NO_FIT or BLOCKER.**

---

## Surprises in the inputs

### Surprise 1: the k=7 jump in |ε_k|·2^k

Pre-registration assumed |ε_k|·2^k might follow a clean asymptotic. Empirical pattern:
- k=2..6: nearly flat at 0.03-0.04 (consistent with |ε_k| ~ C · 2^{-k}).
- k=7,8: jumps to 0.15, 0.19 — a factor of ~4× larger than k=6.

This is the most informative datum in (1). It is **STRUCTURALLY INCOMPATIBLE WITH ANY SINGLE-POLE / SINGLE-SINGULARITY ASYMPTOTIC** (Chevalier 1.16, FS VI.4, Theorem VI.7 polylogarithms). The clean single-power-law n^{α-1} cannot reproduce a transition from one decay rate to another at finite n.

The pattern is consistent with **either**:
(a) A multi-singularity generating function (FS VI.5) where the dominant singularity changes at k ≈ 7 (the slower-decay term becomes dominant). This would mean ε_k has two distinct "modes" — sign-coherent decay in one regime and resurgent in the other.
(b) An *exact* algebraic transition at k = 7 driven by the 7/45 substrate threshold itself: the algebraic equation defining the threshold gains a new branch at k = 7. This is consistent with the c=7/45 being a *boundary* of some algebraic / number-theoretic regime.

The brief notes that the c=7/45 closure is the goal; the empirical k=7 transition in ε_k is suspiciously coincident. This is a **structural observation, not an answer**, but it strongly suggests that any Tauberian closure for c=7/45 must account for the k=7 transition — and the available candidate theorems do not, because they all assume single-singularity / clean-asymptotic behavior.

### Surprise 2: sign pattern + + - - - - - -

ε_1 = +1/5, ε_2 = +1/105 are positive; ε_k for k ≥ 3 are negative. This sign mixing at low k:
- Defeats Theorem 4.2 / Theorem 8.1's positivity / one-sided-bound hypothesis directly applied to ε_k. (Can be patched by |ε_k|, but then h_3 of 1.16 / etc. needs to be reanalyzed for the patched sequence.)
- Is consistent with the generating function f(z) = Σ ε_n z^n having a contribution from a singularity at z = -1 (alternating signs from k=3 onward suggests near-cancellation near a (1+z)^β-type contribution). FS Theorem VI.5 (multiple singularities) would in principle handle this, but the analytic-continuation step is unverified.

### Surprise 3: the C2 BMP F_1 finding is the unweighted-support diffraction, not Syracuse μ_n's weights

This was noted in C2_DISPOSITION but bears repeating in the Tauberian rescope context: the *only* clean BMP / cut-and-project deliverable for Syracuse is the diffraction of the 3-coprime-integer SUPPORT, not the weighted Markov measure. For Tauberian closure, we need the weights. So input (3) is informative about the support (a pure-point spectrum on rationals with cubefree-at-3 denominator) but does not feed a Tauberian theorem's hypotheses about the weighted Dirichlet/Laplace transform.

### Surprise 4: BT archimedean-place finding actively predicts non-fit

Input (4) tells us the c=7/45 closure lives at the *archimedean place* (visible only adelic-ly). All seven Tauberian candidates A-F operate on a *single complex variable* z (and at most one Dirichlet series); none uses adelic structure. So input (4) *predicts* that single-place Tauberian theorems cannot deliver the closure — the disposition matrix confirms this prediction.

---

## What category of theorem is missing

The closure target requires:
1. A theorem operating on a *generating series / Dirichlet series* (so categorically aligned with Tauberian).
2. That handles **multi-singularity / multi-regime asymptotics** (not just single n^{α} decay).
3. That accommodates **adelic / multi-place** input (since the c=7/45 lives at the archimedean place per BT).
4. Whose hypotheses can be verified from the Syracuse structural data (renewal-walk Tao §7, BMP F_1 support, Markov stationary on (ℤ/3^n)^*) — NOT just from a hypothesized analytic continuation that IS the closure target.

No theorem in the 20-PDF corpus satisfies all four desiderata. Specifically:
- A, B handle multi/single singularity but require single-z analytic continuation (Mode H).
- C, D, E handle Dirichlet series but require holomorphic extension past Re z = 1 (Mode H).
- F handles Dirichlet-character L-function decomposition but requires fixed-modulus Dirichlet structure (categorically NO_FIT for Syracuse).
- H3 (Ingham-EM) handles partition-type singularities at q = 1 but requires non-negative coefficients and exponential growth (categorically NO_FIT).
- The adelic / multi-place direction is not represented in any candidate.

**Missing-category description:** a Tauberian theorem for **adelic Dirichlet series / multi-place L-functions / Markov-chain stationary distributions on profinite groups** with **multi-regime asymptotic outputs**. None of the 20 candidates instantiates this.

---

## What input strengthening would unblock a candidate

Per the brief's PARTIAL category, here is what realizable strengthening of inputs would change the disposition (this is reporting, not recommending):

- **For B (Chevalier 1.16) → PARTIAL:** if ε_9..ε_K were available exact-rational for K ≥ 15-20, the asymptotic fit could be tested empirically with confidence. But the k=7 jump suggests fitting will not converge to a clean integer M; this is a structural finding, not a small-K artifact. Strengthening (1) alone is insufficient.
- **For C/D/E → PARTIAL:** an *independent* proof (not from the closure itself) that f(z) = Σ ε_n n^{-z} or Σ ε_n e^{-zt} has analytic continuation past Re z = 1 (resp. Re z = 0) would unblock the Mode H trap. Such a proof would essentially be a *substitute* closure proof; the Tauberian theorem would then deliver an *asymptotic strengthening* of an already-proven result, not the closure itself.
- **For F → still NO_FIT:** strengthening would not help; the failure is categorical (no Dirichlet-character mod-q structure on Syracuse).
- **For A → PARTIAL:** if the dominant-singularity structure of f(z) on |z| = 1 could be identified (e.g. proving two dominant singularities at z = 1 and z = -1 with specific exponents), VI.5 might fire. But identifying these singularities requires substrate-level analytic work beyond the 8-coefficient input.

---

## Final disposition

**NO SELECTED. Probe closes BLOCKER-DOMINANT.**

The Tauberian arc, the SOLE remaining live route after the four prior probe arcs closed NO-GO, fails to deliver a single-theorem selection. The failure mode is uniform and informative:

- **Mode H target-object circularity** (C, D, E, A, B): every analytic-continuation hypothesis is essentially the closure target.
- **Categorical mismatch** (F, H3): wrong object class altogether (Dirichlet-character mod q / partition-type generating series).
- **Pointer-only** (G, H10-H12, H1, H7): not theorem statements.
- **Subsumed by primary** (H4, H6, H9, H2, H8): no new content.

**Combined with the four prior NO-GO arcs (5-probe Fourier-decay, C1 Cochrane/BC/HB exp sums, C2 BMP cut-and-project, BT Bruhat-Tits / BKL billiards), the Tauberian-arc BLOCKER closes the modern-framework-transplant route comprehensively.** All five recognizable framework families (continuous-smooth-dynamical Fourier decay, discrete-arithmetic exponential sums, model-set diffraction, p-adic-tree dynamics, complex-Tauberian analysis) have been honestly probed and either ruled out (NO_FIT / NO-GO) or blocked at the load-bearing step (BLOCKER / Mode H circular).

**The c=7/45 closure does not have a known existing-framework path in the 20-PDF Tauberian corpus, nor in the cumulative literature scanned across all prior probes.**

The k=7 transition in the empirical (1) data and the archimedean-place finding in (4) together suggest the missing-category description above (adelic / multi-place / multi-regime Tauberian for Markov stationary on profinite groups). Whether such a theorem exists in the broader literature outside the 20-PDF corpus is a separate question — out of scope for this single-theorem selection probe.

---

## Files produced

- `TAUBERIAN_RESCOPE_PRE_REGISTRATION.md`
- `TAUBERIAN_RESCOPE_A_HYPOTHESES.md`, `TAUBERIAN_RESCOPE_A_HYPOTHESIS_CHECK.md`
- `TAUBERIAN_RESCOPE_B_HYPOTHESES.md`, `TAUBERIAN_RESCOPE_B_HYPOTHESIS_CHECK.md`, `TAUBERIAN_RESCOPE_B_M_PARAMETER.md`
- `TAUBERIAN_RESCOPE_C_HYPOTHESES.md`, `TAUBERIAN_RESCOPE_C_HYPOTHESIS_CHECK.md`
- `TAUBERIAN_RESCOPE_D_HYPOTHESES.md`, `TAUBERIAN_RESCOPE_D_HYPOTHESIS_CHECK.md`
- `TAUBERIAN_RESCOPE_E_HYPOTHESES.md`, `TAUBERIAN_RESCOPE_E_HYPOTHESIS_CHECK.md`
- `TAUBERIAN_RESCOPE_F_HYPOTHESES.md`, `TAUBERIAN_RESCOPE_F_HYPOTHESIS_CHECK.md`
- `TAUBERIAN_RESCOPE_G_HYPOTHESES.md`
- `TAUBERIAN_RESCOPE_H_HYPOTHESES.md`
- `TAUBERIAN_RESCOPE_DISPOSITION.md` (this file)

Also: PDF text extracts at `C:/Collatz/tauberian_extract/*.txt` for the 20 PDFs.

No git operations performed. Nathan commits manually.
