# FG Candidate F — Le Page CLT / local CLT + Li (Fourier decay of Furstenberg measure on flag variety)

**PDFs:** Le_Page_1982 + Li_Fourier_Decay_Renewal_Spectral_Gaps.pdf.
**Extracted text:** `C:/tmp/fg/le_page_1982.txt`, `C:/tmp/fg/li_renewal_spectral.txt`.

---

## Li 2018 Theorem 1.7 (Fourier decay of Furstenberg measure on flag variety) — VERBATIM (p. 4)

> "Let G be a connected R-split reductive R-group whose semisimple part is simply connected and let G = G(R) be its group of real points. Let μ be Zariski dense Borel probability measure on G with finite exponential moment. Let ν be the μ-stationary measure on the flag variety P.
> For every γ > 0, there exist ǫ_0 > 0, ǫ_1 > 0 depending on μ such that the following holds. For ξ > 0 large enough and any pair of real functions ϕ ∈ C^2(P), r ∈ C^γ(P) such that ϕ is (ξ^{ǫ_0}, r) good, ||r||_∞ ≤ 1 and c_γ(r) ≤ ξ^{ǫ_0}, then
>     | ∫_P e^{iξ ϕ(η)} r(η) dν(η) | ≤ ξ^{-ǫ_1}."

### Hypotheses (typed):

- h_F.1.7.group: G connected **R-split reductive R-group**, semisimple part simply connected. [TYPE (i)]
- h_F.1.7.walk: μ on G(R) Zariski-dense, **finite exponential moment**. [TYPE (ii)]
- h_F.1.7.stat: ν = unique μ-stationary measure on flag variety **P = G/B** (B Borel subgroup of G). [TYPE (iii)]
- h_F.1.7.test_fn: ϕ ∈ C^2(P), r ∈ C^γ(P), ϕ is (ξ^{ǫ_0}, r)-good. [TYPE (iv) — input regularity]

### Conclusion C_F.1.7:

- Polynomial Fourier decay: | ∫_P e^{iξ ϕ} r dν | ≤ ξ^{-ǫ_1} for ξ large. **Polynomial-in-ξ decay rate ǫ_1 > 0.**

### Theorem 1.5 (specialization to SL_2(R), VERBATIM p. 4):

> "Let μ be a Zariski dense Borel probability measure on SL_2(R) with a finite exponential moment. Let X = P(R^2) and let ν be the μ-stationary measure on X. For every γ > 0, there exist ǫ_0 > 0, ǫ_1 > 0 depending on μ such that … |∫ e^{iξϕ(x)} r(x) dν(x)| ≤ ξ^{-ǫ_1}."

### Important note from Li 2018, p. 4:

> "It would also be interesting to establish a similar Fourier decay for the group SL_2(Q_p) and the stationary measure on P^1_{Q_p}." [as of 2018, the p-adic version is **OPEN**]

---

## Le Page Théorème 6 (Local CLT) — already extracted in FG_B_HYPOTHESES.md

Same hypotheses on (P): G = SL(d,R), p with class-B_∞ (high integrability), non-arithmeticity (P_2^* or P_3^*).
Conclusion: Gaussian local-limit theorem for log||g_n…g_1 x||.

---

## Phase 1 — hypothesis × input matrix (Li Theorem 1.7)

| Hypothesis | (1)-(4) Disposition |
|---|---|
| h_F.1.7.group: G R-split reductive | **FAILED** — Syracuse on (Z/3^n)* is abelian profinite; reductive R-group requires non-abelian semisimple part. The p-adic version (SL_2(Q_p)) is explicitly noted as OPEN at the time of writing; the abelian (Z_3^*) version is *trivially* zero-content for the Li framework (no flag variety, no Furstenberg measure as defined). FAILED categorically. |
| h_F.1.7.walk: μ Zariski-dense, finite exp moment | The Syracuse step distribution Geom(2) on N+1 maps to ⟨2⟩ ⊂ Z_3^*; "Zariski dense in SL_d(R)" doesn't apply to abelian. Exponential moment **on the chain step**: ∫ ||g||^η dμ_step depends on what "g" means; if g = 2^{-v} ∈ Z_3 with v ~ Geom(2), then |g|_3 = 3^0 = 1 trivially (units), so exponential moment is trivially satisfied. But this is **vacuous** — the relevant moment in the Li framework is on the matrix norm in SL(d,R), which Syracuse has no analog of. FAILED categorically. |
| h_F.1.7.stat: ν Furstenberg measure on flag variety P | (Z/3^n)* has no flag variety (no non-abelian Lie structure → no Borel subgroup → no G/B). Tao's π_n on (Z/3^n)* is the closest analog, but it's a probability on the abelian profinite group itself, not on a flag variety. FAILED categorically. |

**Phase 1: NO_FIT** on group, walk, and stationary all three.

---

## Le Page Théorème 6 (local CLT):

Same category as Le Page Théorème 1 (SL(d,R), Lyapunov-cocycle structure). FAILED categorically on (Z/3^n)*. See FG_B_HYPOTHESES.md.

---

## Phase 2 — conclusion shape

Li 1.7's conclusion |∫ e^{iξϕ} r dν| ≤ ξ^{-ǫ_1} IS the polynomial-in-A Fourier decay shape we seek. With ǫ_1 dependent on (μ, regularity), this is polynomial decay, polynomial-in-A if we admit ǫ_1 as the A-parameter.

**Conclusion shape: STRONG MATCH** — exactly the kind of polynomial Fourier decay we want.

The obstruction is again purely at the hypotheses side.

---

## Phase 3 — profinite extension

The 2018 Li paper specifically flags "Fourier decay for SL_2(Q_p) and stationary measure on P^1_{Q_p}" as an open question. This is exactly the **p-adic analog** needed. Two routes:

**Route F.3.a:** The p-adic version of Li 1.7 — if proved, would apply to SL_2(Q_3) acting on P^1_{Q_3} (the projective line over Q_3 = the boundary of the Bruhat-Tits tree T_3). Does Syracuse fit?

Connection to Bruhat-Tits / Cluster BT_DISPOSITION:
- BT explicitly considered SL_2(Q_p) random walks on T_p. The disposition was H_BT_NONE_FIT — adelic visibility of c = 7/45 requires global / archimedean place, and the SL_2(Q_3) finite-place machinery alone doesn't see the 1-attractor (the constant 7/45 is an *archimedean* constant).
- Even if Li's p-adic extension were proved, it would give Fourier decay of the SL_2(Q_3)-stationary measure on P^1_{Q_3} — but **Syracuse's stationary is on Z_3^* (a maximal compact in Q_3^*), NOT on the boundary P^1_{Q_3} of the Bruhat-Tits tree**. These are different objects.

**Route F.3.b:** Abelian profinite version — does Li-style spectral-gap analysis apply to (Z/3^n)*?

The proof strategy of Li 1.7:
1. Spectral gap of transfer operator P_z (Pz f(x) = ∫ e^{z σ(g,x)} f(gx) dμ(g)) on Hölder functions on P. The spectral gap on the imaginary line (Re z = 0) is the key.
2. Convert spectral gap to Fourier decay via Dolgopyat-style oscillatory-integral argument.

On (Z/3^n)*: the transfer operator P_z f(x) = E_v[e^{-z v log 2} f(x · 2^{-v} mod 3^n)] (as derived in FG_B Phase 3). This operator on finite-dim functions on (Z/3^n)* has:
- A trivial eigenvalue 1 on constants.
- The character spectrum: P_z χ_ξ(x) = χ_ξ(x) · χ_ξ(2^{-v}) summed against the v-distribution = χ_ξ(x) · E[e^{-2πi ξ 2^{-v}/3^n} e^{-z v log 2}] — this is a *Mellin / Fourier-Laplace* combined transform.
- The spectral gap on Re z = 0: |P_{ib} χ_ξ(x) / χ_ξ(x)| = |E[e^{-ib v log 2 - 2πi ξ 2^{-v} / 3^n}]| ≤ 1 with equality iff phase is constant.

For non-trivial ξ ≠ 0 mod 3^n, this is **exactly the Syracuse Fourier coefficient |μ̂(...) up to factors|**, i.e., the closure target.

**Phase 3 disposition: REDUCES_TO_TARGET (Mode H circular).** The Li-style spectral-gap-on-imaginary-line analysis, when transferred to (Z/3^n)*, becomes the question of how |μ̂_n(ξ)| decays in ξ — the closure target itself.

---

## Disposition F: **MODE_H_CIRCULAR** (and NO_FIT at Phase 1).

- Phase 1: NO_FIT on group/walk/stationary categorically (R-split reductive vs abelian profinite).
- Phase 2: conclusion shape is strong match.
- Phase 3: extension reduces to closure target (the imaginary-line spectral gap of the transfer operator IS the Fourier decay of the chain stationary measure).

This is the load-bearing finding for the FG arc: **the FG framework's most-quantitative Fourier-decay theorem (Li 1.7) translates, in the profinite abelian setting, into the very statement we are trying to prove.** The mechanism that powers Li 1.7 in the SL_d(R) setting is the **non-abelian non-commutativity** giving spectral gap via Dolgopyat oscillation — abelian profinite groups have *commutative* transfer operators (diagonalizable by characters), so the Dolgopyat mechanism is absent.

This is the same finding as in Cluster BT (Bruhat-Tits): adelic visibility / archimedean structure needed; finite-place / abelian-place machinery doesn't deliver. Different route, same wall.
