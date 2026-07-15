# Result 17 — ADVERSARIAL AUDIT of the 2026-07-15 session (R6–R16). Four checks clean; **one real published error found and corrected.**

**Date:** 2026-07-15. **Verdicts: A1 PASS / A2 PASS / A3 PASS / A4 PASS / A5 PASS — but A5 exposed a MISLABELED published number in R15.**

**Headline: the machinery is sound (the power iteration converged everywhere at ≤2.2e−16; the grading holds at k=4, validating R15's extrapolation; dense == sparse; exact == float). But R15 published "off by 0.06%" against 7/15 when the true figure is 0.19% — the code compared to my own prediction while the label said 7/15. Corrected in `result_15`, STATE, and the probe.**

Probe: `probe_17_audit.py`. Log: `result_17_audit_log.txt`. Runtime: **7.6 s**. Motivated by the operator's observation that the session "has been easy so far" — an unusual run of confirmations warrants an attempt to break them.

## The four vulnerabilities, named before testing

**A1 — SINGLE POINT OF FAILURE.** Every ratio in R6–R16 traces to `probe_6.stationary()`, a float power iteration with an iteration cap and a 1e−15 L1 stopping rule. **I never checked whether it converged.** Would kill: *everything*.

> **PASS.** Worst stationarity residual `‖K^Tπ − π‖₁` = **2.166e−16** over 16 cases (q=3..59, k=2..8); simplex error ≤2.2e−16; row-sum error ≤2.2e−16.
> **★ Bonus finding:** it converges in **3–5 iterations**, never near the cap — because `iters ≈ k+1`. That is STATE's **K_k structural lemma** ("K_k mixes in exactly k Markov steps") appearing independently. **The convergence is structural, not lucky** — and R8's derivation predicts it (mod `q^k` the `r_0` term dies after k steps).

**A2 — UNTESTED EXTRAPOLATION.** R14 verified the triangular grading **only at k=3**; R15's `ratio_within(k)` *uses* it at k=4..8, and that produces the headline "within frozen at 0.71958983896 while cross grows at 7/15." Would kill: R15's within/cross split.

> **PASS — EXACT at k=4.** `value(v_1,v_2,v_3,v_4) = 2^{−v_4} + q·2^{−(v_3+v_4)} + q²·2^{−(v_2+v_3+v_4)} + q³·2^{−(v_1+…+v_4)} mod q⁴`:
>
> | q | d | M=ord_{q⁴}(2) | `v_1` mod `d` | `v_2` mod `dq` | `v_3` mod `dq²` |
> |---|---|---|---|---|---|
> | 3 | 2 | 54 ✓ | 26,624 checks, **0 bad** | 24,576, **0 bad** | 18,432, **0 bad** |
> | 5 | 4 | 500 ✓ | 253,952, **0 bad** | 245,760, **0 bad** | 204,800, **0 bad** |
> | 7 | 3 | 1029 ✓ | 351,918, **0 bad** | 345,744, **0 bad** | 302,526, **0 bad** |
>
> **R15's extrapolation from k=3 is validated.** This was the biggest structural worry and it survives.

**A3 — NEW CODE PATH.** R16 introduced a dense chain builder never cross-checked against the sparse one.

> **PASS.** dense == sparse at (3,4), (3,6), (5,3), (5,4), (7,3), (11,3); worst rel diff **1.06e−13**.

**A4 — THIN EXACT VERIFICATION.** Exact rational arithmetic ran at exactly **one** point all session (q=41, k=1, in R10).

> **PASS.** Exact `‖π_k‖²` vs float at three fresh points: q=5/k=2 (4.3e−16), q=7/k=2 (1.9e−16), q=3/k=3 (1.4e−16).

## ★ A5 — the headline survives, and it caught the error

**Design:** R15's "cross grows at 7/15" should be **immune** to A2, because `within` is *constant* in k for k≥4, so cross-differences ≡ total-differences. Test: recompute the slope from **raw totals only**, never touching `ratio_within`.

> **PASS.** `raw total ratios k=2..8` → differences `0.591354, 0.464274, 0.465515, 0.466169, 0.465491, 0.465921`; **mean of last 4 = 0.465774 vs 7/15 = 0.466667 — off 0.19%.** The 7/15 result uses no `within` formula and survives.

### ⚠️ THE ERROR

**R15 published "off by 0.06%" against 7/15. The true figure is 0.19%.**

The probe computed `off = abs(av - 0.4655) / 0.4655` — the deviation from **my own pre-committed prediction (0.4655)** — while the printed label read `vs 7/15=0.46667`. Both numbers were in the line; only one was what the label claimed.

**Propagated into:** `result_15_tower_k_count.md` (header + headline quote + body), the R15 STATE entry, and commit `0c5b4db`'s message.

**Impact:** the *verdict* is unaffected — 0.19% still passes the pre-registered 2% rule, and the agreement with 7/15 remains excellent. **But the published number was wrong, and wrong in the flattering direction**, which is precisely the class of error that does not get caught by accident.

**Corrections applied:**
- `result_15_tower_k_count.md` — header, headline quote, and body corrected to 0.19% with an explicit ⚠️ note stating what the original said and why it was wrong.
- `STATE.md` — R15 entry corrected inline with the same note.
- `probe_15_tower_k_count.py` — **the code fixed**, now reporting *both* figures, each labelled with what it is actually compared against:
  ```
  cross increments k=5..8: mean 0.46577, spread 0.15%
     vs my pre-committed prediction 0.4655 : off 0.06%
     vs 7/15 = 0.466667                    : off 0.19%
  ```
- Commit `0c5b4db`'s message is immutable; this result and the STATE note are the correction of record.

## Net

**Nothing is retracted.** All of R6–R16's structural claims survive: the power iteration is sound, the grading holds a level beyond where it was tested, the two code paths agree, exact matches float, and the 7/15 slope is independent of the formula it was reported alongside.

**One published figure was wrong** and is now corrected in three places plus the source. That is the audit paying for itself: the error was invisible to every check that existed — it lived in a *label*, not in arithmetic, and no consistency test would have found it.

**Lesson (durable):** when a probe compares a measurement to a target, the label and the arithmetic must name the *same* constant. Reporting both a prediction-deviation and a truth-deviation, each explicitly labelled, is the fix — and it is now in `probe_15`.

_Reporting discipline: the four vulnerabilities were named before any of them was tested, each with an explicit "KILLS IF FAILS" list, so a failure would have had a pre-committed blast radius. A5 was designed specifically to test the claim I most wanted to be true. The error found is reported as an error, in the direction it favoured, with its propagation enumerated — not softened by the fact that the verdict held._
