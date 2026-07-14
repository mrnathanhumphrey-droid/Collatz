# R81 disposition — F-hat phase profile on the R78 support

**Date:** 2026-07-13. **Verdict: H_PSEUDO.**

Probe: `result_81_fhat_phase_profile.py`. Data: `result_81_fhat_phase_data.csv`. Log: `result_81_log.txt`.

Smoke (Theorem 78.3, |F_hat|=3√q constant on support, support {a≡1 mod 3}): **PASS**.

## Pass/fail table (per denominator D; r=2 EXCLUDED from evidence)

Cell = fraction of (ℓ,ε) c-family combos (6 total) at that r whose phase index J_a, rescaled to D-th roots, fits a congruence quadratic in a on 100% of support points. `q(α≠0)` marks that a genuine α≢0 mod 3 quadratic-in-a representative exists.

| r | 3^r | 3^(r+1) | 2*3^r | 2q | 4q |
|---|---|---|---|---|---|
| 2 *(excluded)* | 0/6 pass, 0/6 on-roots | 0/6 pass, 0/6 on-roots | 0/6 pass, 0/6 on-roots | 0/6 pass, 0/6 on-roots | 0/6 pass, 0/6 on-roots |
| 3 | 0/6 pass, 0/6 on-roots | 0/6 pass | 0/6 pass, 0/6 on-roots | 0/6 pass | 0/6 pass |
| 4 | 0/6 pass, 0/6 on-roots | 0/6 pass, 0/6 on-roots | 0/6 pass, 0/6 on-roots | 0/6 pass, 0/6 on-roots | 0/6 pass, 0/6 on-roots |
| 5 | 0/6 pass, 0/6 on-roots | 0/6 pass | 0/6 pass, 0/6 on-roots | 0/6 pass | 0/6 pass |
| 6 | 0/6 pass, 0/6 on-roots | 0/6 pass, 0/6 on-roots | 0/6 pass, 0/6 on-roots | 0/6 pass, 0/6 on-roots | 0/6 pass, 0/6 on-roots |

## Decision (§3′ rule)

- H_QUAD fires iff some D gives 100% with an α≢0 mod 3 rep at **every** r∈[3, 4, 5, 6]. Fired denominators: **NONE**.
- Any-polynomial-pass denominators (100% congruence, possibly linear/degenerate) at every r≥3: **NONE**.
- H_PSEUDO fires iff ALL denominators fail at every r≥3.

## Routing (which of the three paper routes this opens/closes)

No low-degree polynomial phase; the phase index is consistent with equidistribution on Z/3^r. **This is a certifying negative** — it **retires the smooth-completion / Gauss-sum route** cleanly: the square-root saving in Σ 1̂(3a)F̂(3a) cannot come from a completed quadratic phase, so R78’s residual bilinear bound genuinely needs Burgess-strength input (the Burgess wall is real). Publishable as a route-closing result; does not weaken Theorems 78.1–78.3.

_Reporting discipline: the fired outcome above is reported as-is, including a null. r=2 carries no evidential weight (|A_2|=3, zero dof). No within-support magnitude filter was applied._
