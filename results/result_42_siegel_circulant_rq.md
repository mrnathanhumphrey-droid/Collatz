# Result 42 (qx+1 paper) — LEAD 2: Siegel's circulant gives NO analytic shortcut to r_q, BUT the Fourier frequency-class decomposition has an exact SHIFT structure that recovers R7 (Pythagoras) and NAMES the primitive as the coprime-frequency second moment.

**Date:** 2026-07-16. **Verdicts: H_PARSEVAL ✓ (gate, F_k = Siegel circulant eigenvalues) / H_SINGLEVAL ✓ r_q is NOT a single Fourier value (q=7 refutes) / H_SECONDMOM ✓ r_q = second-moment rate = build_M / ★ H_DECOUPLE — exact SHIFT structure S_j(k)=S_0(k−j) ⇒ M_k = Σ_i S_0(i) = R7 recovered.**

**Headline: Lead 2 asked whether Siegel's Ch3 circulant theorem (single-copy operator eigenvalues = Fourier values χ̃_N(n)) hands us r_q. Answer: NO shortcut — r_q is the decay rate of the SECOND MOMENT `Σ_n|F_k(n)|²` (Parseval = our Lean), not a single eigenvalue (R28 stands). BUT decomposing that second moment by q-adic frequency class v_q(n) reveals an EXACT SHIFT: `S_j(k) = S_0(k−j)`, so `q^k‖π_k‖² = Σ_{i=0}^k S_0(i)` — a cumulative sum = R7's Pythagoras, recovered from the Fourier side. And it IDENTIFIES R7's primitive `M_i(1)` as the `v_q(n)=0` (coprime/primitive-frequency) part of Siegel's circulant spectrum. At q=3, S_0(i)→7/15 (linear divergence); q≥5 geometric. The circulant UNIFIES R7+R8 and names the primitive, but gives no closed form for r_q.**

Probe: `probe_42_siegel_circulant_rq.py`. Log: `result_42_siegel_circulant_rq_log.txt`. Runtime: ~1 min.

## The setup (Lead 2)

Siegel's dissertation Ch3 (Thm 3.39/Prop 3.61): the single-copy transition operator `M_N` is **circulant**, so its eigenvalues are the Fourier values `χ̃_N(n)`. For our self-similar measure `π` (= Siegel's `dμ_{H,ℓ}`, 2026 Prop 7.3), those Fourier values are `F_k(n) = μ̂(n/q^k) = DFT of π_k`. The bridge to our object is Parseval (our `Parseval.lean`): `‖π_k‖² = q^{−k}Σ_n|F_k(n)|²`. So the question is whether the second moment `Σ_n|F_k(n)|²` inherits the single-copy diagonalization — an analytic route to `r_q`.

## Method

`F_k = np.fft.fft(π_k)` with `π_k` from `stationary(q,k)` embedded into `Z/q^k` (measure lives on units). Verify Parseval gate. Decompose the second moment by q-adic frequency class: `S_j(k) = Σ_{v_q(n)=j}|F_k(n)|²`. q=3,5,7, k up to 5–6.

## Results

**H_PARSEVAL (gate):** `Σ_n|F_k(n)|² = q^k‖π_k‖²` to machine precision at every k, all q. So `F_k(n)` are genuinely Siegel's circulant eigenvalues, and the Parseval bridge holds.

**H_SINGLEVAL — r_q is NOT a single Fourier value.** Single-copy `|F_1(n)|` (Siegel symbol eigenvalues): q=7 → `{1.0, 0.723, 0.723, 0.510, 0.510, 0.467}`; `r_7 = 0.39` is **not present**. (q=5 has `|F_1|=0.609` near `r_5≈0.62`, but this is a non-robust coincidence — q=7 breaks it.) So no single-eigenvalue shortcut.

**H_SECONDMOM — r_q = second-moment rate = build_M.** The increment ratios of `X_k = M_k/(q/3)^k` → r_q (q=5: →0.63; q=7: →0.39; q=3: →1.0), matching build_M. `r_q` lives in the second moment, which is our pair-correlation operator — the circulant does not shrink it.

**★ H_DECOUPLE — the exact shift structure (the real find).** The per-class second moments satisfy
```
    S_j(k) = S_0(k − j)   (exactly)
```
Verified numerically: e.g. q=3 at k=6, classes `v_q = 0,1,2,3,4` = `0.4662, 0.4655, 0.4642, 0.4616, 0.4762` = `S_0(6), S_0(5), S_0(4), S_0(3), S_0(2)`. The classes are **shifted copies of one primitive sequence** `S_0`. Consequently:
```
    q^k‖π_k‖² = Σ_{j=0}^k S_j(k) = Σ_{i=0}^k S_0(i)     (cumulative sum)
```

## What the shift structure means — R7/R8 recovered and the primitive named

- **This IS R7 (Pythagoras) from the Fourier side.** R7 established `X_k = q^k‖π_k‖²` is the *cumulative sum* of the primitive `M_i(1) = q^i‖d_i‖²`. Here `M_k = Σ_{i=0}^k S_0(i)` is exactly that cumulative sum, and **`S_0(i)` (the `v_q(n)=0` coprime-frequency second moment) = R7's primitive `M_i(1)`.** So Siegel's circulant + Parseval + the frequency grading reproduce R7+R8 and *identify the primitive object as the primitive-frequency (coprime-`n`) part of the Fourier spectrum*.
- **The divergence at q=3 is visible in `S_0`.** `S_0(i) → 7/15 ≈ 0.4667` (constant) at q=3 ⇒ `M_k ~ (7/15)k` (linear) = R8's slope. For q≥5, `S_0(i) ~ (q/3)^i` (geometric) ⇒ `X_k` converges.
- **`r_q` is the subdominant rate of `S_0(i)/(q/3)^i → C_q`** — the same object as build_M and R23's increment recurrence. No closed form (R28 stands).

## Verdict on Lead 2

- **No analytic shortcut.** Siegel's circulant diagonalizes the single copy (eigenvalues = Fourier values), but `r_q` is the second-moment decay rate, which is the pair-correlation operator (build_M) — not a single eigenvalue, no closed form.
- **But a genuine unification + naming.** The Fourier/circulant side recovers R7 (cumulative sum) and R8 (Parseval + 7/15 slope), and names R7's primitive as the coprime-frequency (`v_q=0`) component of Siegel's spectrum. `r_q` is thereby named as the **L² Fourier-decay rate of a self-similar measure / Riesz product** — Kahane–Salem–Zygmund territory (the correct harmonic-analysis home for the gap, and where a bound should be sought).
- **For L3:** confirms the character-sum / transfer-operator route is the right one (the circulant doesn't bypass it), and places the object in the self-similar-measure Fourier-decay literature. The shift structure `S_j(k)=S_0(k−j)` is a renewal equation `M(z) = A(z)/(1−z)` with `A(z)=Σ S_0(i)z^i` — so `r_q` = subdominant singularity of the primitive generating function `A(z)`, a clean restatement (though still no closed form).

## Not at stake
R1–R41. This tests whether Siegel's circulant shortcuts r_q (it does not) and unifies the Fourier picture; it changes no value.

_Reporting discipline: the "r_q = a single Siegel Fourier value" shortcut hope was pre-registered AS predicted-to-lose (R28) and lost (q=7 refutes; q=5 near-match flagged as coincidence, not claimed). The shift structure `S_j(k)=S_0(k−j)` is an exact numerical identity (machine precision), reported as the payoff. An embedding bug (measure on units → zero-pad to Z/q^k) was caught via the Parseval gate failing, fixed, and the gate then passed at every k — disclosed, not silent._
