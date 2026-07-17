# Eigenvalue spectrum — diagnostic

## Solver, sparsity, timing

| k | n    | density | nz/row min | nz/row max | nz/row mean | build (s) | eig (s) | solver |
|---|------|---------|------------|------------|-------------|-----------|---------|--------|
| 5 | 162 | 0.3025 | 49 | 49 | 49.0 | 0.00 | 0.01 | dense (np.linalg.eig) |
| 6 | 486 | 0.1008 | 49 | 49 | 49.0 | 0.04 | 0.01 | sparse (scipy.sparse.linalg.eigs, LM) |
| 7 | 1458 | 0.0336 | 49 | 49 | 49.0 | 0.26 | 0.06 | sparse (scipy.sparse.linalg.eigs, LM) |

## Perron-Frobenius checks

| k | |lambda_1 - 1| | n_neg_components (eigvec_1) | PF OK |
|---|----------------|------------------------------|-------|
| 5 | 2.22e-16 | 0 | yes |
| 6 | 4.44e-16 | 0 | yes |
| 7 | 4.44e-16 | 0 | yes |

Kernel construction: float64 K_k from Tao-Syracuse build_markov_float (coprime states r in Z/3^k with r mod 3 != 0; transition r -> ((3r+1)*2^{-v}) mod 3^k with prob (1/2^v)/Z_M).
