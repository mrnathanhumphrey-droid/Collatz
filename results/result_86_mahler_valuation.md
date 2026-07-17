# Result 86 — closed form for the Mahler valuation profile v₃(c_k)

**Date:** 2026-07-14. **Priority: LOW (puzzle-box).** **Verdict: NULL — no candidate reproduces the extended profile out of sample. v₃(c_k) remains deterministic, r-stable, ~1.3k, period-3 with digit-rollover defects, closed form OPEN.**

Probe `probe_86_mahler_valuation_form.py`; data `result_86_data.csv`; log `result_86_log.txt`.

**Not load-bearing.** R81b's analyticity certification, the degree law `degree(r)=max{k:v₃(c_k)≤r}`, the smooth-completion closure, and the log₃2 category obstruction all stand regardless of this outcome. This probe only asks whether the *valuation profile* has an obvious closed form.

## Feasibility (why not k=20)

`v₃(c_k) ≈ 1.3k`, and resolving c_k needs `3^{r+1} > 3^{v₃(c_k)}` ⇒ `r ≳ v₃(c_k)`. k=20 ⇒ v₃≈26 ⇒ r≈26 ⇒ `3^26 ≈ 2.5×10¹²`-term sums (and float loses the integer index past r≈20–22): infeasible. Per operator decision, extended to the feasible ceiling (confirmed to **k=13**). k=20 belongs to Task C's algebraic derivation, not a numerical extension.

## Task A — extended profile with per-k r-stability

| k | v₃(c_k) | status | Δ | ⌊4k/3⌋ | v₃−⌊4k/3⌋ | s₃(k) | carries(k+k) |
|---|---|---|---|---|---|---|---|
| 0 | 0 | confirmed |  | 0 | 0 | 0 | 0 |
| 1 | 2 | confirmed | 2 | 1 | 1 | 1 | 0 |
| 2 | 2 | confirmed | 0 | 2 | 0 | 2 | 1 |
| 3 | 3 | confirmed | 1 | 4 | -1 | 1 | 0 |
| 4 | 4 | confirmed | 1 | 5 | -1 | 2 | 0 |
| 5 | 6 | confirmed | 2 | 6 | 0 | 3 | 2 |
| 6 | 7 | confirmed | 1 | 8 | -1 | 2 | 1 |
| 7 | 8 | confirmed | 1 | 9 | -1 | 3 | 1 |
| 8 | 10 | confirmed | 2 | 10 | 0 | 4 | 2 |
| 9 | 11 | confirmed | 1 | 12 | -1 | 1 | 0 |
| 10 | 12 | confirmed | 1 | 13 | -1 | 2 | 0 |
| 11 | 15 | confirmed | 3 | 14 | 1 | 3 | 1 |
| 12 | 16 | confirmed | 1 | 16 | 0 | 2 | 0 |
| 13 | 17 | confirmed | 1 | 17 | 0 | 3 | 0 |
| 14 | 19 | single-r | 2 | 18 | 1 | 4 | 3 |
| 15 | — | unresolved |  | 20 |  | 3 | 2 |
| 16 | — | unresolved |  | 21 |  | 4 | 2 |

`confirmed` = valuation r-stable (agrees across ≥2 resolving r); `single-r` = resolved at only one r (tentative, not evidential); `unresolved` = modulus too small at every tested r. Only `confirmed` points carry weight.

## Task B — candidate models, out of sample

Fit window k≤11; out-of-sample k≥12 = **[12, 13, 14]**. Decision rule: integer coefficients, reproduce **every** point incl. k≥12, or REFUTED.

| model | integer coeff | in-window | out-of-sample | verdict |
|---|---|---|---|---|
| floor(4k/3) only | — | no | FAILS | refuted |
| + a*s3(k) | — | no | n/a | refuted |
| + a*carries(k+k) | — | no | n/a | refuted |
| + a*v3(k!) | — | no | n/a | refuted |

**All candidates refuted.** The Legendre/Kummer family (`⌊4k/3⌋` + integer·{s₃, carries, v₃(k!)}) does not reproduce the profile out of sample. Consistent with the pre-registered most-likely null and with the author's 0-for-4 structural-prior record in this arc. The period-3 `1,1,2` difference structure with digit-rollover defects is real and visible in the table, but no closed form in these features captures it.

## Task C — derivation

Not earned — no model survived Task B. The closed form for v₃(c_k) is reported **empirical and OPEN**. The natural derivation route (Mahler coefficients of `x↦4^x=(1+3)^x` from binomial structure) is where a real answer would come from; the numerical profile alone does not yield it.

## Scope

Untouched, explicitly: **R81b's certification**, the degree law, the smooth-completion closure, the log₃2 obstruction, THEOREM_C_745, and Thms 78.1–78.3. This is a low-priority puzzle-box and yields the machine to anything else that needs it.

_Reporting discipline: extended before fitting (Task A first); out-of-sample decision rule enforced; integer coefficients required; r-stability verified per point; the pre-registered candidate treated as a hypothesis to kill._
