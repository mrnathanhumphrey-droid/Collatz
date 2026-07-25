# Probe ILEDGER — the interference ledger (second-moment source decomposition) — **GATE PASSES; the two competing quantities are now separately computed and separately signed: the FRESH source ⟨s_j,Re w⟩ is uniformly NEGATIVE (toward 7/15), the PROPAGATION is positivity-definite (W⁽ᵏ⁾→ an all-positive fixed point, toward 0.477); the net is a knife's-edge magnitude competition the decomposition EXHIBITS but does not sign — outcome three, banked plainly.**

**Date:** 2026-07-25. Probe `probes/probe_ILEDGER.py`. Wilson's reformulation: δ_r is a second moment; the map on it
is not the pushforward of the measure map (that's why the measure-source was zero, checkpoint (c)). Diagonal
(|A|²+|B|²) = |D|² kernel; the cross term 2Re(A B̄) = the source (bilinear, phases). Recursion
`δ_r = T̃_diag δ_{r−1} + s_r`. Projection = iterated pullback of the fixed weight (no eigenvector — it doesn't
exist): `W⁽⁰⁾=Re w`, `W⁽ᵏ⁺¹⁾(x)` = |D|²-weighted average of W⁽ᵏ⁾ over the 3 preimages {x/3,(x+1)/3,(x+2)/3},
`wD(y)=1/(5−4cos2πy)`. All in the **dlog** coordinate (β linearizes ×4→+1; refinement→×3; R14 profile).

## P4 — GATE PASS (the #43 check)
`Σ_{j=2}^r ⟨s_j, W⁽ʳ⁻ʲ⁾⟩ = g_r = ⟨δ_r, Re w⟩` for r=2..7, exact to ~10⁻⁹ (`s_j` and `W⁽ᵏ⁾` computed by
**separate** code paths — s_j via forward `T̃_diag`, W⁽ᵏ⁾ via the pullback — built as an exact adjoint pair on the
dense primitive grid, so a disconnect would break the telescope). **Decomposition valid; convention consistent.**

## P1 — the pullback POSITIVIZES (this is the key structural fact)
`W⁽ᵏ⁾` from Re w (positive-arc Haar measure 0.4196):

| k | 0 | 1 | 2 | 3 | … | ∞ |
|---|---|---|---|---|---|---|
| Haar(W>0) | 0.4196 | **1.000** | 1.000 | 1.000 | | 1.000 |
| W⁽ᵏ⁾(0) | 0.3333 | 0.2275 | 0.1953 | 0.1891 | | **0.18758** |

**After a single pullback, W⁽ᵏ⁾ is positive everywhere**, and it converges to an all-positive fixed point
(W^∞(0)=0.18758, eigenvalue 1). The |D²|-weighted 3-fold pullback concentrates mass at DC (where wD=1, Re w=+⅓
maximal) and **destroys Re w's sign structure** — exactly your "fragmentation-toward-zero" turning into
positivization. So **any propagated source contribution (k≥1) is positive.**

## P2/P3 — the ledger splits by k: fresh source NEGATIVE, propagated POSITIVE
`A[j,k] = ⟨(T̃_diag)^k s_j, Re w⟩`:

| j\k | 0 (fresh) | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| 2 | **−7.42e−3** | +4.88e−3 | +1.08e−3 | +2.05e−4 | +3.9e−5 |
| 3 | **−2.02e−3** | +8.21e−4 | +1.69e−4 | +3.21e−5 | +6.1e−6 |
| 4 | **−4.96e−4** | +4.43e−4 | +9.17e−5 | +1.75e−5 | +3.3e−6 |
| 5 | **−1.14e−4** | +6.02e−4 | +1.20e−4 | +2.29e−5 | +4.4e−6 |
| 6 | **−1.49e−3** | +5.45e−4 | +1.17e−4 | +2.23e−5 | +4.2e−6 |
| 7 | **−2.35e−4** | +3.83e−4 | +8.12e−5 | +1.55e−5 | +2.9e−6 |

**Column k=0 (the freshly-manufactured coupling `⟨s_j, Re w⟩`) is NEGATIVE for every j=2..7.** Columns k≥1
(propagated through the all-positive W⁽ᵏ⁾) are POSITIVE for every j. So each source injects **one negative blip at
birth, then a positive decaying tail.**

## P5 — the two rates, and the competition
- **Propagation (row) rate → 0.190** (the pullback fixed-point contraction; clean, D-independent).
- **Source (column) rate ~0.66 in norm** but **erratic**: `‖s_j‖` = 0.23, 0.19, 0.13, 0.090, 0.059, 0.039, while
  `⟨s_j,Re w⟩` is non-monotone (the **s₆ anomaly**, −1.49e−3 ≫ its neighbors, is exactly the r=6 sign flip).
- **Net:** `g_r = Σ_j A[j,r−j] = (one fresh NEGATIVE, k=0) + (accumulating POSITIVE tail, k≥1)`. A single negative
  term against a sum of positives: generically the positive sum wins (g_r>0 → 0.477), and it does at r=3,4,5,7;
  the negative wins only when the fresh source is anomalously large (r=2, r=6). The effective margin
  `≈ 2.1·|c₁| − |c₀|` per source **straddles zero** (j=3 → net negative → 7/15; j=2,4,7 → net positive → 0.477).
  **Knife's edge, undecided.**

## Verdict — outcome three (exhibit, don't sign), banked plainly
The interference ledger **cleanly separates and signs the two competing quantities** you named:
**fragmentation-toward-zero → positivization (propagation POSITIVE, toward 0.477)** vs **the fresh source
(NEGATIVE, toward 7/15)**. It does **not** decide the net — that is a magnitude claim (`2.1|c₁|` vs `|c₀|`,
straddling 1), and reading it as a sign result would be #38. So: **sign-indefinite, on a knife's edge, undecided.**
Structurally it leans 0.477 (a single fresh negative vs an accumulating positive sum), but the margin is within the
source's own oscillation (the s₆ anomaly flips it), so the lean is not a verdict.

## ⚠️ Convention caveat (flagged, not patched)
`‖s_j‖/‖δ_j‖ ≈ 1.0–1.6` — the source is **O(1)×δ, not small**. With the **normalized-average** pullback (as you
specced — "weighted average"), the |D|² diagonal does **not** make the source a small interference correction. Your
other phrase — the diagonal is "**mass-contracting**" — points instead at the **un-normalized** |D|² kernel, which
would rescale W⁽ᵏ⁾ (no eigenvalue-1 fixed point) and shrink the source. The GATE passes for either (it telescopes),
so this is a genuine convention fork I did not resolve: **the sign results above are for the normalized-average
kernel as written.** If the certified diagonal is the mass-contracting (un-normalized) one, the source and its sign
should be re-read there before the knife's-edge lean is trusted.

## Status
**ILEDGER: gate PASS; the coupling competition is now exhibited with both sides separately signed** — propagation
positivity-definite (W⁽ᵏ⁾→all-positive fixed point, W^∞(0)=0.1876), fresh source `⟨s_j,Re w⟩` uniformly NEGATIVE
(j=2..7). Rates: propagation → 0.190, source ~0.66 (erratic; s₆ anomaly = the r=6 flip). **Net sign undecided —
outcome three, on a knife's edge (`2.1|c₁|` vs `|c₀|` straddles 1); structurally leans 0.477 (single fresh negative
vs accumulating positive) but within the source's own oscillation.** Convention fork flagged: normalized-average
(specced) vs mass-contracting (un-normalized) — sign results are for the former; re-read needed if the latter is the
certified diagonal. Not at stake: R1–R30, R80–R82. The 7/15-vs-0.477 decider (sign of the net coupling) remains
open — now reduced to a signed, computable competition rather than an eigenvalue question. **USER (pen):** confirm
the diagonal normalization, and — if signable — rank `2.1|c₁|` vs `|c₀|` (or its correct form) asymptotically.
