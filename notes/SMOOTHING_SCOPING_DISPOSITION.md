# Smoothing-route scoping probe — DISPOSITION

**Verdict: H_DELTA_EXISTS_TRANSFER_BROKEN.**

The pre-registered most-likely outcome is confirmed. Of five candidate Syracuse-derived δ_a constructions, three (C1, C2, C3) fall at the C² entry gate; one (C4) advances marginally with degenerate UNI strength; one (C5, the most structurally honest candidate using T_lead's (1, 4) eigenvector) clears the entry gate (C² + UC + UNI non-vacuous) but **fails the canonical-relation-to-μ_n gate via the structural pincer**:

> A C² conjugacy h: Φ' → Φ_3 exists **iff** Φ' is C²-conjugate-to-linear **iff** ARHW Theorem 1.1 does NOT fire. The conjugacy that gives the canonical ν → μ_n relation is precisely the conjugacy ARHW's hypothesis forbids.

This pincer is structural, not technical. It is the same failure mode as Probe 2 Candidate A (T1 transfer) and Probe 3 Candidate (b) (object-identity gap), now confirmed not to be fixable by upgrading the perturbation to a Syracuse-derived one.

## Candidate summary

| Candidate | Source data | C² | UC | UNI | Canonical-relation | Verdict |
|---|---|---|---|---|---|---|
| C1: v_2 perturbation | Syracuse-native | **FAIL** | — | — | — | No canonical C² interpolation of v_2 on ℝ |
| C2: T_lead (1,4) class-indicator | T_lead eigenstructure | **FAIL** | — | — | — | Class indicator τ is discontinuous |
| C3: Tao mod-3 stationary density | Syracuse Markov stationary | **FAIL** (pincer) | — | — | — | Non-constant ⟹ discontinuous; constant ⟹ affine |
| C4: parity-density Gaussian-smoothed | parity branches + ad-hoc σ | PASS | PASS (small ε) | marginal, m∝ε | **FAIL** | UNI degenerate as ε → 0, transfer broken |
| C5: class-broken two-branch IFS | T_lead (1,4) eigenvector | PASS | PASS | PASS | **FAIL** | Transfer broken via canonical-relation pincer |

## What closed each gate

- **Phase 1 (UNI recap).** Framework requires C² globally on [0,1], not piecewise. α exponent depends implicitly on UNI strength m; α → 0 as m → 0 (degenerate-limit risk for any smoothing route).
- **Phase 2 (candidates).** Three of five candidates fail at C² because canonical Syracuse data lives on discrete arithmetic objects (Z_3, parity sequences, v_2 valuations) with no canonical smooth interpolation to ℝ. Two pass the entry gate (C4 marginal, C5 honest).
- **Phase 3 (transfer).** Both surviving candidates fail at canonical-relation-to-μ_n via two independent obstructions: (a) the conjugacy pincer (the conjugacy that gives ν → μ_n is the same one ARHW forbids), (b) the real-vs-3-adic Fourier incompatibility (characters don't restrict canonically across κ: Z_3 ↪ [0,1]).
- **Phase 4 (limit behavior).** Moot — no candidate reaches Phase 4.

## Adversarial checks resolved

- **(A1) Pre-registration honesty.** Pre-registered favoring H_DELTA_EXISTS_TRANSFER_BROKEN; that is the outcome. C5 was the most optimistic candidate and was honestly evaluated; the conjugacy pincer is a real structural argument, not a manufactured negative.
- **(A2) Derivation fidelity.** C5 and C3 use Syracuse data verbatim (T_lead eigenvector / mod-3 stationary); C1 uses v_2 directly; C2 uses class indicator from v_2 parity. C4 partial — Gaussian width is ad hoc. The classification is strict: candidates were not labeled "Syracuse-derived" without an explicit derivation, and the ones that pass A2 still fail at later gates.
- **(A3) UNI verification by computation.** For C5: two distinct C² branch maps φ_+, φ_- with g_+ ≠ g_- give d/dx log φ'_+ ≠ d/dx log φ'_- generically, UNI fires non-vacuously. For C4: m(ε) = O(ε/σ²), UNI is non-vacuous for ε > 0 but degenerate in the ε → 0 limit. Verified.
- **(A4) Transfer-mechanism honesty.** Two independent obstructions (conjugacy pincer + Fourier-category mismatch) confirmed. The real-vs-3-adic Fourier incompatibility is the SAME structural issue that closed Probe 2 T1 transfer and Probe 3 Candidate b — not papered over.
- **(A5) Not repeating Probe 3 Candidate b.** C5's derivation is structurally Syracuse-encoding (T_lead's (1, 4) eigenvector is the asymptotic class-mass ratio, a genuine Syracuse invariant). The failure at the transfer gate is NOT due to ad-hoc-ness but due to the structural conjugacy pincer — the discovery is that even a properly Syracuse-encoding perturbation cannot bridge the real-vs-3-adic Fourier gap.

## What this means for the smoothing route

The ARHW + Syracuse-derived smoothing route is **structurally closed**, not merely lacking a candidate. The closure mechanism is dual:

1. **Inside-out:** the conjugacy that gives canonical ν → μ_n is exactly the conjugacy ARHW's non-linearity hypothesis forbids. Inside the IFS framework, the two conditions cannot both hold.

2. **Outside-in:** real Fourier (on [0,1]) and 3-adic Fourier (on Z_3 or Z/3^n Z) are different harmonic analyses on different topological groups. There is no general theorem converting polynomial decay of one to polynomial decay of the other.

This dual closure is the same pattern as L²-flattening (Probe 1) and SL_2 embedding (Probe 2): **discrete-arithmetic Markov-chain stationary measures don't fit modern Fourier-decay frameworks built for continuous/smooth-dynamical settings**, and smoothing does not bridge the gap because the gap is about the **Fourier-analytic category** of the target object, not about the smoothness of the source.

## Routing

Four routes through continuous/smooth-dynamical frameworks have now closed structurally:
- L²-flattening (Probe 1): discrete log|D_v| = −v log 2 lies on a single AP, no Plancherel collapse.
- SL_2 embedding (Probe 2): T_lead is rank-1 (det = 0), no SL_2(ℝ) action.
- Cocycle Dolgopyat (Probe 3): natural Syracuse IFS Φ_3 is affine (f'' ≡ 0), framework excludes.
- Smoothing through ARHW (this probe): Syracuse-derived δ_a clears entry gate but fails canonical-relation transfer.

**Three open routes remain, in order of structural promise:**

1. **Probe 5 — Markov-chain-native drift conditions (arxiv:2005.08145).** This is now the highest-priority next probe. Syracuse IS a Markov chain on Z_3; a framework that takes Markov-chain mixing as primary, rather than smooth-IFS structure or continuous-group action, is the natural next category to try. Doesn't require any of the obstructions that closed Probes 1-4.

2. **Tauberian arc.** Per MEMORY: the 73-PDF Tauberian bundle at C:/Collatz/burgess/literature/ was opened for the Flajolet-Sedgewick Ch. VI / Chevalier 2507.15394 Thm 1.16 candidate. Single-theorem selection is pending ε_7 exact-rational compute. This route doesn't require Fourier-decay at all — it works on generating-function singularity structure, which is natively combinatorial / arithmetic.

3. **Probe 4 — certified-approximation transfer operator (arxiv:2602.19435).** Rigorous-numerics route. T_lead's rank-1 problem persists but may be circumvent-able via V_{M^k} (higher-order moment subspaces). Lower priority than (1) and (2) because the rank-1 obstruction is the same Probe 2 obstacle in a different framework.

**Recommended next step: Probe 5 (Markov-chain drift conditions).** If Probe 5 also lands negative, the structural-category signal is unambiguous: smooth-dynamical, continuous-group, smooth-IFS, and Markov-chain-native frameworks have all closed. At that point routing collapses to the Tauberian arc OR Bourgain-Konyagin discrete sum-product OR new technique.

## Deliverables

- `SMOOTHING_SCOPING_UNI_RECAP.md` — Phase 1, ARHW UNI quantitative dependence
- `SMOOTHING_SCOPING_CANDIDATES.md` — Phase 2, five candidates with verdicts
- `SMOOTHING_SCOPING_TRANSFER.md` — Phase 3, canonical-relation gate closure
- `SMOOTHING_SCOPING_DISPOSITION.md` — this file
- Phase 4 (`SMOOTHING_SCOPING_LIMIT.md`) deliberately omitted — no candidate survives Phase 3

## Pre-registration audit

Pre-registered favoring **H_DELTA_EXISTS_TRANSFER_BROKEN** (most likely), with H_DELTA_EXISTS_BUT_UNI_DEGENERATE as fallback. Outcome: **H_DELTA_EXISTS_TRANSFER_BROKEN confirmed for C5** (the strongest candidate), **H_DELTA_EXISTS_BUT_UNI_DEGENERATE confirmed for C4** (the marginal candidate). The two pre-registered "exists past the entry gate but fails deeper" outcomes both materialized, with C5 being the more important confirmation because its derivation was most rigorously Syracuse-encoding (T_lead's (1, 4) eigenvector). Pre-registration was correct.
