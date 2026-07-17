# Result G0c′ (CARRIED GATE — blocks Phase 2): the amplitude-carrying modal PASSES at q=7 (0.38 vs bar 0.39±0.05) — the previously-failing case, now confirmed via ESPRIT on the exact ‖π_k‖² sequence. q=5 control is ESPRIT-UNRESOLVABLE (within 0.60 ≈ cross 0.62, 3% apart) — an estimator-resolution miss, NOT a framing miscalibration. Secondary: q=3 build_M is DEFECTIVE (R39 confirmed over R32). Adjudication to Nathan.

**Date:** 2026-07-16. Carried gate G0c′ (blocks Phase 2). Method: `build_M` at L=3 walls out (q7 dim 7.4M / ~1.6e11 nnz; q5 dim 1.25M / ~1e10 nnz), so — better than L=3 — ESPRIT/matrix-pencil on the exact `‖π_k‖²` sequence (`= Σ A_i λ_i^k` over build_M's amplitude-carrying eigenvalues, by the mass identity `Σ M^k v₀ = ‖π_k‖²`) recovers those eigenvalues at effectively L=∞. v-truncation vmax=64 all channels; sequence gate-validated (probe_27 references).

**Pre-registered bars (committed): q=7 = 0.39 ± 0.05 (primary, L=2 known bad 0.475); q=5 = 0.60–0.63 (control).**

## PRIMARY — q=7: PASS (0.38)

ESPRIT p=3 on `‖π_k‖²` (k=1..7, to n=705894):

| mode | \|λ\| | \|λ\|/λ₁ | amplitude A |
|---|---|---|---|
| Perron | 0.33334 | 1.0000 | +1.366 |
| **cross (r_q)** | **0.12661** | **0.3798** | **−0.2245** |
| sub | 0.05173 | 0.1552 | −0.033 |

**Amplitude-carrying subdominant = 0.3798 ∈ [0.34, 0.44] → PASS.** It sits **below** the within-mode ratio `3/7 = 0.4286`, matches probe_27's gate-validated cross-rho `0.39`, and beats the discredited L=2 raw modal `0.475`. The (M)-spectral framing gives r_q at the discriminating prime — **the first time q=7 modal has passed.**

## CONTROL — q=5: ESPRIT UNRESOLVABLE (explained miss, not a framing failure)

ESPRIT p=3 on `‖π_k‖²` (k=1..8, to n=312500):

| mode | \|λ\| | \|λ\|/λ₁ | amplitude A |
|---|---|---|---|
| Perron | 0.33392 | 1.0018 | +1.203 |
| (whisper) | 0.22092 | 0.6628 | −7.0e−4 |
| (dominant sub) | 0.16094 | 0.4828 | −0.150 |

**Dominant-amplitude subdominant = 0.48, a whisper at 0.66 — neither cleanly in [0.60, 0.63]. By the letter, a MISS (off by >0.05).** Cause, disclosed as a deviation not waved away: at q=5 the two subdominant eigenvalues are `1/q = 0.200` (within, ratio 3/5=0.60) and `r_q/3 = 0.207` (cross, ratio r_5≈0.62) — **3% apart**, unresolvable by ESPRIT from 8 points; the 0.48 estimate is a blend artifact, not a true mode. **This is an estimator-resolution failure BY DESIGN — q=5 is the concealing prime where within ≈ cross (the standing trap's own premise) — not a miscalibration of the (M) framing** (both true modes 0.60/0.62 lie inside the bar; ESPRIT simply cannot split them). The gate-validated cross-rho anchors q=5 at **0.62** (probe_27), which the successive-ratio estimator resolves where ESPRIT cannot.

## The honest tension (report, not patch)

- **q=7 (the discriminating prime, within 0.4286 ≠ cross 0.39): clean PASS at 0.38.** The (M)-spectral framing is confirmed where the modes are separable.
- **q=5 (the concealing prime, within 0.60 ≈ cross 0.62): independent modal UNAVAILABLE** — ESPRIT can't resolve 3%-apart modes from 8 points; only the (c)/cross-rho estimator gives 0.62, which is not an *independent* modal confirmation.
- **Standing trap (q=5 AND q=7 minimum):** q=7 satisfied independently; q=5 satisfiable only via (c), not (M). So the INDEPENDENT modal is load-bearing at q=7 but not at q=5.
- **ON-FAIL diagnosis does NOT fit:** the pre-registered fail-consequence was "(M)-spectral framing miscalibrated." It is not — q=7 passed. The q=5 shortfall is estimator resolution on a near-degenerate pair, not framing. So this is neither a clean pass (q=5 unmet by the letter) nor the pre-registered failure mode. **Reporting and stopping per protocol; Nathan adjudicates.**

## SECONDARY (recon, non-gating) — q=3 defectiveness: DEFECTIVE confirmed (R39 over R32)

q=3 build_M, L=2 (dim 324): `|λ₁|=0.34683, |λ₂|=0.34392, ratio 0.99161`; **top-pair eigenvector overlap `|⟨v₁,v₂⟩| = 0.9983` (→1 = merging), `cond(R) = inf` (eigenvector matrix SINGULAR — no eigenbasis).** So at finite L the q=3 top pair is **genuinely defective**, not distinct-and-merging: the R39 exceptional-point / Jordan reading is confirmed over the R32 distinct-modes reading. **The d=2 boundary is a real defective point** — carry to Phase 2 as an asset (the isometry/marginal argument breaks against a Jordan block, exactly as §9's boundary clause requires). (L=3 recon was cut by a background-wrapper issue; L=2 is conclusive — cond=inf.)

## Adjudication requested (Nathan)
Two readings of the trap:
1. **Accept:** q=7 clean-pass confirms the modal framing at the prime that can discriminate; q=5 is unresolvable *by nature* (within≈cross) and anchored at 0.62 by cross-rho. The framing is not miscalibrated. → Phase 2 proceeds on build_M.
2. **Hold:** demand an *independent* q=5 modal. Options: push ESPRIT to higher k (q=5 k=9 = 1.95M states, heavy) — but the 3% mode separation may stay unresolvable at any reachable k; or a degeneracy-aware estimator. If q=5 independent-modal is deemed structurally unavailable, restate Phase 2's q=5 checks in (c)-coordinates (cross-rho) while keeping (M) for q=7 and the q=3 boundary.

## Not at stake
R1–R46, Phase 0 closure. No r_q value changes (r_7=0.39, r_5=0.62 reconfirmed). This is the independent-modal check: PASS q=7, unresolvable q=5, q=3 defective.

_Reporting discipline: build_M L=3 infeasibility disclosed (not silently downgraded); ESPRIT used on the exact sequence = the R32 amplitude spectrum at L=∞. The q=5 deviation reported AS a deviation (0.48 vs bar 0.60–0.63) with mechanism (3% mode separation), NOT as "consistent with." No patching, no re-deriving around: reported and stopped for adjudication. q=3 defectiveness is exact (cond=inf), not a threshold call._
