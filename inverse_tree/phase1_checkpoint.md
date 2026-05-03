# Phase 1 Checkpoint — Inverse Collatz Tree (d ≤ 50)

## Wall time
Tree construction: 0.072s

## Node counts
Total nodes (d=0..50): **379,600**

| depth | n_nodes |
|------:|--------:|
| 0 | 1 |
| 5 | 2 |
| 10 | 6 |
| 15 | 24 |
| 20 | 72 |
| 25 | 227 |
| 30 | 732 |
| 35 | 2,365 |
| 40 | 7,628 |
| 45 | 24,561 |
| 50 | 79,255 |

## Empirical growth rate (window d=10..50)
- slope = 0.2343
- intercept = -0.4309
- R² = 0.999865
- asymptote log(4/3) ≈ 0.2877 (finite-window slope expected lower)

## OEIS A005186 cross-check (d=0..20)
Result: **PASS**

## Forward-σ verification
Sample size: 10,000
Result: **PASS**

## Parquet output
Path: `C:/Collatz/inverse_tree/tree_d50.parquet`
Size: 4170.4 KiB (4,270,511 bytes)
Schema:
```
  n: Int64
  depth: Int64
  parent: Int64
  child_left: Int64
  child_right: Int64
  r_6: Int64
  ell_6: Int64
  a_final_6: Int64
  c_final_6: Int64
  alpha_det_6: Float64
  r_8: Int64
  ell_8: Int64
  a_final_8: Int64
  c_final_8: Int64
  alpha_det_8: Float64
  r_10: Int64
  ell_10: Int64
  a_final_10: Int64
  c_final_10: Int64
  alpha_det_10: Float64
```

## Surprises / anomalies
(none flagged at this checkpoint)