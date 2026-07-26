# Probe MODREG — the modulation regression (Wilson's cheap check) — **⟨cos⟩ under the pure transport modulation w is EXACTLY 1/4 (ŵ(2)/ŵ(0), with the m=2m′ frequency doubling); but regressing the measured Ĉ(m) against a SINGLE w(2⁻¹m mod N) over 3∤m gives R²≈0 (0.0003 at r=12, 0.0000 at r=16) with a positive-but-vanishing correlation (+0.017→+0.006). So the single-w modulation is the RIGHT SIGN but the WRONG SHAPE: d1>0 does NOT reduce to "one modulation component positive." This is what's expected if the true accumulated modulation is the LACUNARY PRODUCT Π_j w(2⁻ʲm), not one factor — and R16-A does NOT certify that product (the transport is non-diagonal / a tower map between growing spaces, R16-crux + R17 caveat), so the product form is banked as Wilson's structural conjecture, unconfirmed. The next move is the lit hunt (dispatched), not compute.**

**Date:** 2026-07-25. Probe `probes/probe_modreg.py`, log `logs/modreg_run.log`. No build_nu (uses `scratchpad/rho_12,16.npy`). Tests Wilson's pen claim that d1's sign is carried by the transport-modulation-aligned component of the fluctuation spectrum.

## The exact number (Wilson's, verified analytically)
The Geom(½) transport kernel in the dilated coordinate `u=m′/N` is `w(u)=1/(5−4cos2πu) = 1/|1−2e^{iθ}|² = (1/3)Σ_n 2^{−|n|}e^{inθ}`, so `ŵ(n)=2^{−|n|}/3`. Since `dlog(−2)=2⁻¹` gives `m=2m′`, the weight `cos(2πm/N)=cos(4πu)` is at the **doubled** frequency n=2, so
$$\langle\cos\rangle_w=\frac{\hat w(2)}{\hat w(0)}=\frac{2^{-2}/3}{1/3}=\frac14.$$
**The white part of Ĉ contributes exactly 0 to ⟨cos⟩; the pure modulation contributes +1/4.** So d1 = ⟨cos⟩ ≈ 0.00194 (r=16) against a pure-modulation 0.25 means the fluctuation spectrum is white + a **~0.8% modulation-shaped ripple** — the entire sign lives in that ripple.

## The cheap check — regress Ĉ(m) ~ a + b·w(2⁻¹m mod N) over 3∤m
| r | slope b | R² | corr | ripple/white | ⟨cos⟩ (=banked d1) |
|---|---|---|---|---|---|
| 12 | +9.78e−8 | 0.0003 | +0.0170 | 2.01e−2 | +2.9636e−3 |
| 16 | +4.66e−10 | 0.0000 | +0.0059 | 7.72e−3 | +1.9392e−3 |

- **Sign: b > 0 at both levels** (correct direction — aligns with +1/4 and d1>0; consistent with the POINCARE dilated-band big-mod-positive finding).
- **Shape: R² ≈ 0** — a single `w(2⁻¹m)` explains essentially none of Ĉ's variance, and the weak positive correlation is **shrinking** (+0.017→+0.006).
- **Verdict (Wilson's dichotomy):** "if the fit is poor, the ripple isn't modulation-shaped and the 1/4 is a coincidence of the null." The fit IS poor ⟹ **the single-w 1/4 is a null coincidence; d1>0 does NOT cleanly reduce to a single provable modulation-amplitude sign.** The clean structural reduction fails the fit test.

## Why this doesn't kill the modulation idea — the product form
A single-w regression **cannot** capture a PRODUCT of r dilated copies. Wilson's accumulated modulation after r transport steps is the lacunary/Riesz-type product `Π_j w(2⁻ʲm mod 3^r)` — a product along the 2-adic dilation orbit of m. A poor single-w fit is exactly what that predicts (w is only the last factor). So MODREG redirects, not refutes: the operative object is `⟨cos⟩` of `Π_j w(2⁻ʲm)`, not the single w's 1/4.

## R16-A confirmation of the product form — NOT certified (Wilson's flag upheld)
R16-A (`result_transport_R16.md`): one transport step `θ_r(t)=E_v E_{X~μ_{r−1}}[1(dlog₄(1+3·2⁻ᵛX)=t)]` reproduces the frozen shell exactly (gate PASS r=2..6). Per level this IS the Geom(½) convolution (symbol w) in the dlog domain. **But the transport is NOT a clean per-level multiplication:** it is "a tower map between growing spaces" (R16 crux) and "the exact transport is **not diagonal**" (R17 caveat, `result_slowmode_R17.md`) — `ν̂_r(ξ)=e(ξ/3^{r+1})·E_v[ν̂_{r−1}(ξ2⁻ᵛ)]` is a coherent SUM over the whole Geom ladder of dilations, not a single-frequency factor. **So the clean product `Π_j w(2⁻ʲm)` is a heuristic from the shape of the recursion, NOT derivable from the certified transport.** Banked as Wilson's structural conjecture, explicitly unconfirmed against R16-A — one line of confirmation still owed by the pen before it's a fact.

## The named object for the lit hunt (dispatched, "Hank")
`⟨cos⟩` of a lacunary/Riesz-type product `Π_j w(2⁻ʲm mod 3^r)` = the cosine moment of a Riesz product built from the ×2 (×4) orbit on Z/3^r. Prior mapping on disk: R42 named `r_q` as the "L² Fourier-decay rate of a Riesz product" (Kahane–Salem–Zygmund); `L2_FLATTENING_DISPOSITION` recommended the Furstenberg-measure Rajchman route (Hochman–Solomyak) but flagged **all Fourier-decay machinery is sign-blind** (bounds |·|, not the sign — useless for a positivity question); R21 placed the per-level subgroup sum on the Konyagin shelf (bounds in q, not k). Lit hunt asks: is the **SIGN** of this cosine moment KNOWN / OPEN / rigidity-adjacent-HARD, and is any result sign-bearing rather than magnitude-only.

## Status
**MODREG (Wilson's cheap check):** ⭐**⟨cos⟩ under pure modulation w = EXACTLY 1/4** (ŵ(2)/ŵ(0), m=2m′ doubling; verified analytically). ⚠️**Single-w regression: sign RIGHT (b>0) but shape NULL (R²=0.0003→0.0000, corr +0.017→+0.006 shrinking)** ⟹ single-w 1/4 is a null coincidence, d1>0 does NOT reduce to one provable modulation-amplitude sign. ⭐**Redirects to the LACUNARY PRODUCT Π_j w(2⁻ʲm)** — the real accumulated modulation; a poor single-w fit is exactly what the product predicts. ⚠️**R16-A does NOT certify the product** (transport non-diagonal / tower map between growing spaces, R16-crux + R17); product = Wilson's structural conjecture, UNCONFIRMED — one line owed by the pen. ⭐**Next move = lit hunt (dispatched), not compute:** is the SIGN of ⟨cos⟩ of a Riesz product along the ×2 orbit on Z/3^r known/open/hard; is any result sign-bearing (Fourier-decay machinery is sign-blind). Not at stake: R1–R30, R80–R82, d1 ladder to r=16, MODES/RATIO-2/AC-LAGS/SPECTILT/POINCARE. rho_12..16 in scratchpad. commit pending; lit-hunt result pending.
