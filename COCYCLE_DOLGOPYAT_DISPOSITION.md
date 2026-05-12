# Cocycle Dolgopyat structural-compatibility probe — DISPOSITION

**Verdict: H_COCYCLE_DOLGOPYAT_LINEAR_EXCLUSION.**

The natural Syracuse-derived IFS — base-3 expansion Φ_3 = {x → (x+a)/3 : a ∈ {0,1,2}} on [0,1] — is **linear in ARHW's exact technical sense** (f''_a(x) ≡ 0 on K_Φ for every map). ARHW Theorem 1.1 requires the IFS to be **not conjugate to a linear IFS** (verbatim, p. 2); Corollary 1.2 requires "Φ contains a non-aﬃne map" (verbatim, p. 3). Φ_3 satisfies neither — it is the exact prototype of what the framework excludes.

The exclusion is structural, not technical. From ARHW p. 7 (Claim 2.2 proof): if all σ-periodic codings give equal cocycle-derivative limits, Φ is C² conjugate to linear. For Φ_3 the cocycle c(I,x) = −log|f'_I(x)| = |I| log 3 is **constant in x**, so d/dx(log f'_ξ − log f'_ζ) ≡ 0 identically for any codings ξ, ζ — UNI fails by definition. Non-linearity ⟺ UNI ⟺ spectral-gap mechanism fires. None of those fire for Φ_3.

## Adversarial checks resolved

- **(A1) Non-linearity exclusion load-bearing:** Confirmed. Just as Probe 2 fell on det = 0 (load-bearing geometry), this probe falls on f'' ≡ 0 (load-bearing dynamics). Both are vanishing-discriminant failures of identical structural type — discrete arithmetic objects landing in degenerate cases of frameworks designed for non-degenerate smooth-dynamical settings.
- **(A2) UNI verification:** UNI fails trivially for Φ_3 — derivative cocycle is constant in x.
- **(A3) Transfer from ν to μ_n:** Doesn't exist generically. ν lives on [0,1] with real-Fourier decay; μ_n lives on ℤ_3 with 3-adic Fourier. Real characters e^{2πiqx} don't restrict to 3-adic characters of ℤ_3. Same failure mode as Probe 2's T1 transfer. See `COCYCLE_DOLGOPYAT_TRANSFER.md` (b).
- **(A4) §5 inherited-claim discipline:** The corpus INDEX claim "T_lead transfer operator is a cocycle case" is verified to be **terminology overlap, not structural identification**. T_lead is a 2×2 rational rank-1 matrix on the finite-dim V_M space (Probe 2); ARHW's P_s is an infinite-dim operator on C¹([0,1]) built from a smooth-IFS derivative cocycle. Different categories. The INDEX claim does not survive verification — same pattern as the inherited claims falsified earlier in T_LEAD_CORRECTED and NISOLI_CLOSURE_CORRECTED.
- **(A5) PDF discipline:** All 53 pages extracted via pypdf (C:/tmp/arhw_full.txt). Phase 1 hypotheses include verbatim Theorem 1.1, Corollary 1.2, Claim 6.1 (analytic dichotomy), Definition 2.3 (derivative cocycle and transfer operator), Claim 2.2 (UNI), Theorem 2.4 part (5) (UNI in all parts), Theorem 2.8 (spectral gap). Fidelity discipline held.

## Routing

Three probes (L²-flattening / SL_2 embedding / cocycle Dolgopyat) all return structural negatives via the same meta-pattern: **discrete-arithmetic Markov-chain stationary measures don't fit modern Fourier-decay frameworks built for continuous/smooth-dynamical settings.** The failures are not technical hypothesis misfits; they are category-of-object mismatches:

- L²-flattening (Probe 1): BKS Step 3 separation requires the cocycle generator log|D_v| = −v log 2 to NOT lie on a single arithmetic progression. In the discrete Syracuse setting it does. Structural negative on Plancherel collapse.
- SL_2(ℝ) embedding (Probe 2): T_lead is rank-1 (det = 0), no SL_2(ℝ) action on P^1. Structural negative on action geometry.
- Cocycle Dolgopyat (Probe 3): Natural Syracuse IFS is linear (f'' ≡ 0); ARHW explicitly excludes linear/self-similar IFSs from polynomial-Fourier-decay conclusion. Structural negative on non-linearity hypothesis.

In each case the framework was built for a regime (smooth dynamics, continuous-group action, non-affine smooth IFS) that Syracuse's arithmetic-residue dynamics simply does not occupy.

**Next routing options:**

1. **Probe 4 — transfer-operator certified approximation (arxiv:2602.19435).** Rigorous-numerics approach. May be more accommodating to finite/discrete arithmetic transfer operators since it doesn't require smooth-IFS hypotheses — it works with the actual finite matrix data. T_lead's rank-1 structure remains a problem there too, BUT the certified-approximation framework may permit a different angle (e.g., perturbing T_lead within V_M, looking at higher-order V_{M^k}).

2. **Probe 5 — drift conditions (arxiv:2005.08145).** Markov-chain-native framework. Doesn't require IFS structure or smooth dynamics; works with the Syracuse Markov chain on its own terms. If Syracuse μ_n admits a Lyapunov function with the right drift, polynomial Fourier decay can sometimes be derived purely from chain mixing rates. **This is the most natural framework given that Syracuse IS a Markov chain.**

3. **Genuinely new technique.** Three negatives in a row through frameworks built for continuous-smooth settings is itself a signal: the polynomial-in-A bound for Syracuse may require a discrete-arithmetic Fourier-decay technique that does not yet exist in the literature in the right form. Closest existing fragments: Bourgain-Konyagin (in the bundle at C:/Collatz/Bourgain-Konyagin) on sum-product in ℤ/p^k ℤ; Burgess (in the bundle at C:/Collatz/Burgess) on character sums; the 73-PDF tauberian bundle (per MEMORY) explicitly opened for the Flajolet-Sedgewick Ch. VI / Chevalier 2507.15394 Thm 1.16 candidate.

**Recommended next step: Probe 5 (drift conditions, arxiv:2005.08145).** Markov-chain native, doesn't require IFS hypotheses, doesn't require continuous-group action. If Probe 5 also lands negative, that's the third structural-category negative (smooth-dynamical / continuous-group / Markov-chain-native frameworks all fail) and the routing becomes "Flajolet-Sedgewick / Chevalier tauberian arc OR Bourgain-Konyagin discrete sum-product OR new technique."

## Deliverables

- `COCYCLE_DOLGOPYAT_HYPOTHESES.md` — verbatim ARHW Theorem 1.1, Corollary 1.2, Definitions 2.3/Claim 2.2/Theorem 2.4(5)/Theorem 2.8, Claim 6.1 (analytic dichotomy)
- `COCYCLE_DOLGOPYAT_CANDIDATES.md` — four candidate IFS/cocycle structures, all fail
- `COCYCLE_DOLGOPYAT_TRANSFER.md` — even if framework fired, transfer ν → μ_n breaks
- `COCYCLE_DOLGOPYAT_DISPOSITION.md` — this file

## Pre-registration audit

Pre-registered favoring H_COCYCLE_DOLGOPYAT_LINEAR_EXCLUSION. Outcome: **H_COCYCLE_DOLGOPYAT_LINEAR_EXCLUSION confirmed.** Pre-registration was honest and correct.
