# Probe R5 — the deviation law — **candidate DEAD** (clean kill at k=3)

**Date:** 2026-07-21  Exact rationals, zero new computation. Source: frozen `CollatzVerify/Basic.lean` S₁..S₆
(= `s_infinity_exact.py` output). Candidate under test (flagged, D9 rules): **S_k = 7/15 + 1/(5·21^{k−1})**.
Verdict: **DEAD** — matched the two hand-derived welds (k=1,2), then wrong *sign* from k=3 on. No partial credit.

## R5-A — the judge table (rationals compared as rationals)
| k | d_k = S_k − 7/15 (exact) | float | candidate 1/(5·21^{k−1}) | verdict |
|---|---|---|---|---|
| 1 | **1/5** | +0.20000000 | +1/5 | **MATCH** |
| 2 | **1/105** | +0.00952381 | +1/105 | **MATCH** |
| 3 | **−5191/1019445** | **−0.00509199** | +1/2205 | **DEAD** (wrong sign) |
| 4 | −11346676448406637/4627031617157687115 | −0.00245226 | +1/46305 | DEAD |
| 5 | −(92434…4871)/(80256…3045) | −0.00115175 | +1/972405 | DEAD |
| 6 | −(49242…9259)/(98899…7805) | −0.00049791 | +1/20420505 | DEAD |

The candidate's clean geometric was a **two-point coincidence**: the hand-derived S₁ = 7/15+1/5 and
S₂ = 7/15+1/105 are exact and real (Judge One from the v-parity mod-3 law; Judge Two from the mod-9 stationary
measure), but 1/5, 1/105 do *not* continue as 1/(5·21^{k−1}). At k=3 the true deviation is **negative** — the
sequence has already **crossed 7/15 and overshot**, while the candidate stays positive. Kill is unambiguous.
(Cross-check: S₃ = 31370/67963 = R1's exact c₃ — the frozen sources agree.)

## R5-B — the true deviation sequence (the real fingerprint)
**Signs: +, +, −, −, −, −.** The deviation approaches 7/15 **from above** for k=1,2, **crosses zero between
k=2 and k=3, and approaches from below** thereafter — an **overshoot**, not a monotone one-sided decay.

**Tail contraction ratios** (k ≥ 3, where the sign is settled):
| ratio | value |
|---|---|
| \|d₄\|/\|d₃\| | 0.4816 |
| \|d₅\|/\|d₄\| | 0.4697 |
| \|d₆\|/\|d₅\| | 0.4323 |

The ratio is **≈ 0.43–0.48 and drifting down** — *not* the candidate's 1/21 = 0.0476, and not (yet) a constant.
The overshoot + a sub-geometric, still-varying ratio is the signature of **competing contributions of opposite
sign** (consistent with the two-mode / partner-braid structure: the c₀-adjacent term positive, the condensation
term negative, crossing near k=2.5). The denominators (1019445 = 3·5·7²·19·73, then exploding) carry no clean
power structure — there is no simple one-term geometric law.

## What this hands the derivation
- **7/15 is NOT read off a closed one-term deviation law.** S∞ = 7/15 remains a *limit*, reached by an
  overshooting, sign-flipping sequence whose contraction ratio is drifting (~0.43 at k=6, not settled).
- **The real target** is the two-sided law behind the +,+,−,−,−,− fingerprint — a difference of (at least) two
  geometric-ish terms whose crossover produces the overshoot. The exact d₃..d₆ above are its constraints.
- The two hand-derivations (Judge One 2/3, Judge Two 10/21) stand untouched — they are exact facts about S₁, S₂;
  only the *extrapolation* died. D9 protocol worked exactly as intended: candidate flagged, tested on frozen
  exact data, killed at the first term that discriminates.

## Status
Candidate S_k = 7/15 + 1/(5·21^{k−1}) **DEAD** (k=3, wrong sign, exact). True deviations d₁..d₆ frozen above
(Basic.lean). Fingerprint: overshoot at k≈2.5, sign +,+,−,−,−,−, tail ratio ~0.43–0.48 drifting. No fitting;
rationals compared as rationals; the honest sequence is the derivation's next input.
