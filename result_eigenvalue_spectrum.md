# Eigenvalue spectrum of K_k (k = 5, 6, 7)

Top-10 eigenvalues by |lambda| of the Tao-Syracuse Markov kernel K_k on coprime states of Z/3^k. Pure data acquisition, no fitting.

## Summary table

| k  | n    | |lambda_2| | |lambda_3| | gap = 1-|lambda_2| | |lambda_2|-1/2 | |lambda_2|-|lambda_3| | PF OK |
|----|------|------------|------------|--------------------|------------------|-----------------------|-------|
| 5 | 162 | 0.00034186 | 0.00034183 | 0.99965814 | -0.49965814 | +0.00000003 | yes |
| 6 | 486 | 0.00044126 | 0.00044126 | 0.99955874 | -0.49955874 | +0.00000000 | yes |
| 7 | 1458 | 0.00181652 | 0.00180154 | 0.99818348 | -0.49818348 | +0.00001498 | yes |

## k = 5 (n = 162)

```
  i        Re(lambda)        Im(lambda)        |lambda|     arg(lambda)
  1   +1.000000000000   +0.000000000000  1.000000000000   +0.0000000000
  2   +0.000341861688   +0.000000000000  0.000341861688   +0.0000000000
  3   +0.000105595004   +0.000325114879  0.000341833277   +1.2567517561
  4   +0.000105595004   -0.000325114879  0.000341833277   -1.2567517561
  5   -0.000276525848   +0.000200877753  0.000341787092   +2.5133452171
  6   -0.000276525848   -0.000200877753  0.000341787092   -2.5133452171
  7   +0.000036474327   +0.000036479584  0.000051586205   +0.7854702348
  8   +0.000036474327   -0.000036479584  0.000051586205   -0.7854702348
  9   -0.000036474327   -0.000036469068  0.000051578769   -2.3562665774
 10   -0.000036474327   +0.000036469068  0.000051578769   +2.3562665774
```

## k = 6 (n = 486)

```
  i        Re(lambda)        Im(lambda)        |lambda|     arg(lambda)
  1   +1.000000000000   +0.000000000000  1.000000000000   +0.0000000000
  2   +0.000356960748   +0.000259395336  0.000441256066   +0.6284068500
  3   +0.000356960748   -0.000259395336  0.000441256066   -0.6284068500
  4   -0.000136392486   +0.000419569616  0.000441182018   +1.8850979842
  5   -0.000136392486   -0.000419569616  0.000441182018   -1.8850979842
  6   -0.000441136526   +0.000000000000  0.000441136526   +3.1415926536
  7   +0.000244434983   +0.000000000000  0.000244434983   +0.0000000000
  8   -0.000244094755   +0.000000000000  0.000244094755   +3.1415926536
  9   +0.000075447796   +0.000231797341  0.000243767055   +1.2561209926
 10   +0.000075447796   -0.000231797341  0.000243767055   -1.2561209926
```

## k = 7 (n = 1458)

```
  i        Re(lambda)        Im(lambda)        |lambda|     arg(lambda)
  1   +1.000000000000   +0.000000000000  1.000000000000   +0.0000000000
  2   +0.001816520615   +0.000000000000  0.001816520615   +0.0000000000
  3   +0.001095724010   +0.001430012756  0.001801540337   +0.9169863266
  4   +0.001095724010   -0.001430012756  0.001801540337   -0.9169863266
  5   -0.000438250528   +0.001708150651  0.001763474460   +1.8219435774
  6   -0.000438250528   -0.001708150651  0.001763474460   -1.8219435774
  7   -0.001565734041   +0.000728938304  0.001727099922   +2.7058771766
  8   -0.001565734041   -0.000728938304  0.001727099922   -2.7058771766
  9   -0.000716094806   +0.001349959211  0.001528130113   +2.0585113259
 10   -0.000716094806   -0.001349959211  0.001528130113   -2.0585113259
```

## Decision rubric (from brief)

- |lambda_2| ~ 0.5 with clean gap to |lambda_3| < 0.4: single-mode rate-1/2.
- |lambda_2| ~ 0.5 but |lambda_3| also ~ 0.5: multi-mode, comparable rates.
- |lambda_2| significantly off 1/2: rate-1/2 conjecture wrong at algebraic level.
- Spectrum character changes with k: K_k structure is k-dependent.

## Empirical verdict

- **k=5**: |lambda_2| = 0.000342 (diff to 1/2 = -0.499658), |lambda_3| = 0.000342 (gap |lambda_2|-|lambda_3| = +0.000000).
- **k=6**: |lambda_2| = 0.000441 (diff to 1/2 = -0.499559), |lambda_3| = 0.000441 (gap |lambda_2|-|lambda_3| = +0.000000).
- **k=7**: |lambda_2| = 0.001817 (diff to 1/2 = -0.498183), |lambda_3| = 0.001802 (gap |lambda_2|-|lambda_3| = +0.000015).

Numeric checks:

- max |lambda_2 - 1/2| over k=5..7: 0.499658
- min (|lambda_2| - |lambda_3|) over k=5..7: 0.000000
- |lambda_2| drift across k=5..7: 0.001475
- any |lambda_2| off 1/2 by >0.05?  yes
- any |lambda_2|-|lambda_3| < 0.05? yes

## Verdict

**Decision branch 3 fires at the tested k:** |λ_2| is far from 1/2 at k = 5,
6, 7. Branches concerning asymptotic spectrum behavior are not addressed —
three data points don't constrain the k → ∞ limit.

Three observations from the data (scope: k ∈ {5, 6, 7}):

1. **|λ_2| is three orders of magnitude below 1/2 at every tested k.**
   k=5: 3.42e-4. k=6: 4.41e-4. k=7: 1.82e-3. Spectral gap is essentially 1
   in this regime. K_k mixes in O(1) iterations at the tested k. The
   rate-1/2 envelope on epsilon_n is not located in K_k's within-level
   spectrum at these k. (This direct measurement extends
   `result_77_4_K_spectrum_erratum.md` to k=7 — K_k is not behaving as the
   rate operator at any k tested. STATE.md's R71 supersession ("R71
   conjecture λ_2(K_k) = 1/2 → R71.B + R73: convergence rate from
   level-lifting, NOT chain spectrum") is consistent with this.)

2. **The top 10 |λ_i| form a tight cluster, not a spectrum with one
   isolated sub-leading mode.** At k=7 the magnitudes run 0.001817,
   0.001802, 0.001802, 0.001763, 0.001763, 0.001727, 0.001727, 0.001528,
   0.001528 — all within ~16% of |λ_2|. No "second eigenvalue 1/2 with a
   clean gap to the rest" at any tested k.

3. **|λ_2| changes with k:** 3.42e-4 → 4.41e-4 → 1.82e-3, factors of
   ×1.29 (k=5→6) and ×4.12 (k=6→7). The trajectory is non-monotone in
   step ratio; three data points do not establish a fixed asymptotic
   limit.

**Complex eigenvalues with non-zero argument are present at every tested k.**
These produce structural oscillation under iteration of K_k. They live at
magnitude ~10⁻³, so their contribution to K_k^n dies in O(1) within-level
steps; they cannot drive long-range structure in iterates of K_k itself.
Whether they connect to inter-level structure is a separate question this
probe doesn't address.

**Implication.** The within-level operator K_k is not the rate operator at
the tested k. Whatever generates the empirical convergence behavior of
S_k → 7/15 (whose current best fit is the order-3 recurrence with ρ_slow
≈ 0.83 from STATE.md) is not located in K_k's spectrum at k=5,6,7. This
is consistent with prior R71 supersession; the new content is direct
spectrum measurement at k=5,6,7 confirming the same conclusion at higher k.
