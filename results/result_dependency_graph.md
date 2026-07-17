# Result Dependency Graph

Directed graph showing which results consume which data sources / earlier
results. Reading: `A -> B` means B depends on A's data, methodology, or output.

## Data sources (roots)

```
Collatz forward map  m -> (3m+1)/2^v
        |
        +--> [forward orbits]                  random m_0 in [1, 2^N], walk to 1
        |       |
        |       +--> sigma cache               per-class mean sigma at modulus 2^k
        |       +--> chang_qsd_test.csv        D(r, t) at multiple t, 2M-10M orbits
        |       +--> qsd_late_t_avg.csv        D_avg averaged over late t
        |       +--> per-orbit v_t / log m_t   raw step traces (R61 temporal)
        |       +--> R60 K kernel                empirical (r, b) -> (r', b') counts
        |
        +--> [inverse tree]                    BFS from m=1 backward, value-truncated
        |       |
        |       +--> tree node set             {m <= max_value : 1 reachable forward}
        |       +--> subtree sizes (variant a)
        |       +--> branching counts (variant b)
        |       +--> sigma-weighted (variant e)
        |       +--> R23 BFS density           depth-asymptotic density per residue
        |
        +--> [transfer operators]              algebraic constructions on Z/2^k
                |
                +--> Chang P (forward, mod 64, depth-13 lift averaging)
                +--> M_closed (inverse, mod 2^k, doubling + inverse-3)
                +--> M_closed parametrized M(s) (Bowen-pressure form)
```

## Result dependencies

```
sigma cache -----+--> R23 lambda_max          eigenvalue of M_closed
                 |   - source: pure algebra (no sigma data needed)
                 |
                 +--> R24/R25 mean sigma per residue
                 +--> Tao bridge equation
                       - K_h = 3/log(4/3) is INPUT (Tao 2022)
                       - alpha_det = +6.23 is closed form (prefix algebra, no data)
                       - epsilon = -2.45 is empirical (sigma cache)

chang_qsd_test.csv (D at t=90)
                 |
                 +--> R50 QSD computation     fits 7 absorbing-set kernels
                 +--> R51 cylinder depth-extension test (5 absorption x 5 depths)
                 +--> R52 inverse tree (4 weighting schemes) -- depth-truncated
                 +--> R58 inverse tree (variant a/c/e) -- VALUE-truncated
                       - same source as R50/R51, different observable
                 +--> R57 R23 eigvec vs D_avg
                 +--> R59 zadic measure framework (variant b, sigma-band)

qsd_late_t_avg.csv (D_avg)
                 |
                 +--> R53 renewal kernel on visit-event space
                 +--> R60 size-stratified Markov (residue x log size)
                       - kernel from forward orbits (NEW data this session)

inverse tree -+
              |
              +--> R23 BFS density        depth-50 BFS per residue
              +--> R52 inverse tree weighting (depth-trunc)  REJECTED
              +--> R58 inverse tree weighting (value-trunc)  +0.86 Pearson
              +--> R59 zadic mass-dim_q2 sweep (k=5..15)
              +--> R61 spatial Z_q on Z_2 cylinders
                    - same inverse tree as R52/R58/R59
                    - extends q=2 sweep to q in [-5, 5]

R23 lambda_max --+
                 +--> R57 eigvec test against D_avg            REJECTED
                 +--> dim_h_validation: log(lam_max)/log(2)    walked back

Chang P ---------+
                 +--> Chang dim_H = log(phi)/log(2) = 0.694
                 +--> Chang invariant core I_2 = {7, 27, 31, 59, 63} mod 64
                 +--> Chang stationary pi (used by R50, R53, R60 normalization)

M_closed --------+
                 +--> R23 lambda_max
                 +--> Furstenberg branching dim 0.338
                 +--> Bowen-pressure parametrization (no closed match)
```

## Cross-validation graph

These pairs are sometimes cited as mutually confirming. The audit's verdict
on each:

```
R23 lam_max <----?----> Chang dim_H               (β) different operators, same dynamics
R58 prediction <--+0.92--> R60 prediction         (β) DECISIVE — same identification
R58 empirical <--+0.92--> R60 empirical           (β) different observables but correlated
R61 D_1 = 0.608 <-no algebra-> R23 0.675          (β) loose ~10% agreement
R61 D_1 = 0.608 <-no algebra-> Chang 0.694        (β) loose ~10% agreement
Chang heuristic 0.68 <-?-> R23 derived 0.675      (β) loose 0.005 agreement (heuristic)
TA1 epsilon (forward orbits) <-no agreement test->
                                bridge K_h (Tao)  (β) same data, K_h is INPUT
```

## Independent inputs

These are the GENUINELY independent inputs that the body of work draws on:

1. **Collatz map definition** — irreducible. All results derive from this.
2. **Tao 2022's K_h heuristic** — external citation, not re-derived.
3. **Chang 2603.11066v6** — external paper; we use Chang's pi and invariant
   core but compute everything else ourselves.
4. **Cramer-Lundberg rate function I(0) = 0.1465** — external classical
   random-walk theory; cited.

Everything else is OUR computation on Collatz dynamics. This means
independence claims must trace back to one of these four roots — and
in practice, almost everything traces to root #1 (Collatz map itself),
making most cross-result claims (β) consistent characterizations rather
than (α) independent confirmations.

## Two genuine cross-validations

The audit identifies just two pairs that ARE genuinely cross-validating:

1. **R51 cylinder consistency theorem (D_K invariant in K under fixed
   cylinder-21 absorption) confirms Chang's depth-13 cylinder kernel
   is the right depth.** Both depend on Collatz dynamics (same root) but
   the test is whether D_K is K-stable, which is an internal mathematical
   claim verified algorithmically. (β-strict — both depend on map.)

2. **Bridge equation's slope-1.000 verification of Tao's K_h.** K_h is
   from Tao (root #2), and our forward-orbit sigma data tests whether
   K_h is the correct empirical slope. This DOES use independent input
   (Tao's heuristic) and tests it against our data. **Genuine
   cross-validation: empirical sigma-data confirms Tao's K_h numerically.**

Outside these two, all "confirmations" trace back to single Collatz-map root
and should be framed as consistency, not independence.
