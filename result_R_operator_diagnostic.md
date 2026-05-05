# R̃_k operator — diagnostic

## K_k cache (sparse build, sparse stationary)

| k | n | nnz(K_k) |
|---|---|----------|
| 4 | 54 | 2916 |
| 5 | 162 | 26244 |
| 6 | 486 | 236196 |
| 7 | 1458 | 1565892 |
| 8 | 4374 | 4697676 |

## Construction details

- Lift L_k: (n_k, n_{k+1}) sparse matrix, 3 nonzeros per row.
- Within-level kernel K_{k+1}: (n_{k+1}, n_{k+1}) sparse, approximately 49 nonzeros per row (after 2^{-v} float64 cutoff).
- Projection P: (n_{k+1}, n_k) sparse, 1 nonzero per row.
- R̃_k = L · K_{k+1}^m · P built via sparse multiplication, densified only for final eigvals call.

## Verification residuals (Step 1, all k)

- π_{k+1} @ P = π_k holds to ~1e-15 (machine precision). Confirms Syracuse dynamics descends cleanly through projection mod 3^k.
- L_A @ π_k = π_{k+1} fails by O(1) (uniform lift is wrong for the real stationary).
- L_B @ π_k = π_{k+1} holds to ~1e-15 by construction.

## Numerical conditioning notes

- All matrices in float64. R̃ is dense after sparse composition.
- Largest dense R̃ is at k=7: 1458 × 1458 ≈ 17 MB.
- Eigvals via np.linalg.eigvals (LAPACK GEEV).
