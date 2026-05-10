# δ(r) Smoke Test — Closed-Form Probe of B1 Random Effect Residual at k=3

**Verdict (one line):** STRUCTURED, CLOSED-FORM CANDIDATE.

Probe '1: prefix_steps(r)' fires with |Spearman rho| = 0.951. Treat as smoke-test positive worthy of follow-up; 4 points cannot settle the closed form.


---


## Section 1: Inputs


### B1 HMC posterior random-effect parameters (n=10000, train=8000)


| r (mod 8) | u_r posterior mean | u_r posterior SD |
|---|---|---|
| 0 | -0.13371 | 0.03769 |
| 1 | +0.04266 | 0.03777 |
| 2 | -0.03516 | 0.03757 |
| 3 | +0.05402 | 0.03793 |
| 4 | -0.04536 | 0.03735 |
| 5 | -0.01423 | 0.03751 |
| 6 | +0.03502 | 0.03759 |
| 7 | +0.11877 | 0.03779 |

### Prefix predictions (k=3)


| r | a_final at k=3 | j (=log_3 a_final) | prefix_steps | log(a_final) |
|---|---|---|---|---|
| 0 | 1 = 3^0 | 0 | 3 | 0.00000 |
| 1 | 9 = 3^2 | 2 | 5 | 2.19722 |
| 2 | 3 = 3^1 | 1 | 4 | 1.09861 |
| 3 | 9 = 3^2 | 2 | 5 | 2.19722 |
| 4 | 3 = 3^1 | 1 | 4 | 1.09861 |
| 5 | 3 = 3^1 | 1 | 4 | 1.09861 |
| 6 | 9 = 3^2 | 2 | 5 | 2.19722 |
| 7 | 27 = 3^3 | 3 | 6 | 3.29584 |

*Collision note*: r=1 and r=3 both terminate at a_final=9=3^2 at k=3. 4 odd residues collapse onto 3 distinct a_final values (per the prefix-collapse structure noted in §3 of the paper).*


### k=6 aggregated values (8 residues mod 64 sharing each mod-8 class)


| r mod 8 | residues mod 64 | a_final values at k=6 | mean log(a_final) k=6 |
|---|---|---|---|
| 0 | 0, 8, 16, 24, 32, 40, 48, 56 | 1, 9, 3, 9, 3, 3, 9, 27 | 1.6479 |
| 1 | 1, 9, 17, 25, 33, 41, 49, 57 | 27, 81, 27, 27, 81, 243, 9, 81 | 3.8451 |
| 2 | 2, 10, 18, 26, 34, 42, 50, 58 | 27, 9, 81, 9, 9, 3, 27, 27 | 2.7465 |
| 3 | 3, 11, 19, 27, 35, 43, 51, 59 | 27, 27, 81, 243, 9, 81, 27, 81 | 3.8451 |
| 4 | 4, 12, 20, 28, 36, 44, 52, 60 | 9, 9, 3, 27, 27, 27, 9, 81 | 2.7465 |
| 5 | 5, 13, 21, 29, 37, 45, 53, 61 | 9, 9, 3, 27, 27, 27, 9, 81 | 2.7465 |
| 6 | 6, 14, 22, 30, 38, 46, 54, 62 | 9, 81, 27, 81, 27, 27, 81, 243 | 3.8451 |
| 7 | 7, 15, 23, 31, 39, 47, 55, 63 | 81, 81, 27, 243, 243, 243, 81, 729 | 4.9438 |

## Section 2: δ(r) table


| r | u_r (mean) | log(a_final) k=3 | δ_raw | δ (de-meaned) | u_r SD |
|---|---|---|---|---|---|
| 0 | -0.13371 | 0.00000 | -0.13371 | +1.51145 | 0.03769 |
| 1 | +0.04266 | 2.19722 | -2.15456 | -0.50940 | 0.03777 |
| 2 | -0.03516 | 1.09861 | -1.13377 | +0.51140 | 0.03757 |
| 3 | +0.05402 | 2.19722 | -2.14320 | -0.49804 | 0.03793 |
| 4 | -0.04536 | 1.09861 | -1.14397 | +0.50120 | 0.03735 |
| 5 | -0.01423 | 1.09861 | -1.11284 | +0.53232 | 0.03751 |
| 6 | +0.03502 | 2.19722 | -2.16220 | -0.51704 | 0.03759 |
| 7 | +0.11877 | 3.29584 | -3.17706 | -1.53190 | 0.03779 |

- max|δ(r)| = 1.53190
- δ_noise = sqrt(mean(σ_u_r²)) = 0.03765
- max|δ| / δ_noise = 40.687
- Threshold for informativeness: max|δ| > 1.5 × δ_noise = 0.05648


## Section 3: Probe correlations


| Probe | x values | Spearman ρ | p | Pearson r | p | Note |
|---|---|---|---|---|---|---|
| 1: prefix_steps(r) | 3, 5, 4, 5, 4, 4, 5, 6 | -0.9512 | 0.0003 | -0.9999 | 0.0000 |  |
| 2: j(r) | 0, 2, 1, 2, 1, 1, 2, 3 | -0.9512 | 0.0003 | -0.9999 | 0.0000 |  |
| 3: log a_final k6_agg - log a_final k3 | 1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 1.65 | nan | nan | nan | nan | constant — undefined correlation |
| 4: r mod 4 | 0, 1, 2, 3, 0, 1, 2, 3 | -0.6343 | 0.0912 | -0.7718 | 0.0249 |  |
| 5: r | 0, 1, 2, 3, 4, 5, 6, 7 | -0.5952 | 0.1195 | -0.6283 | 0.0953 |  |
| 6: parity of prefix_steps | 1, 1, 0, 1, 0, 0, 1, 0 | -0.1091 | 0.7970 | -0.0037 | 0.9931 |  |

*With only 4 data points, Pearson r values near ±1 are not surprising; Spearman ρ is the conservative test.*


## Section 4: Verdict


**Category: STRUCTURED, CLOSED-FORM CANDIDATE**


Probe '1: prefix_steps(r)' fires with |Spearman rho| = 0.951. Treat as smoke-test positive worthy of follow-up; 4 points cannot settle the closed form.


**Probe 3 (k=3 → k=6 delta) status: DOES NOT FIRE** (Spearman ρ = +nan, Pearson r = +nan). Probe 3 is the load-bearing test: it asks whether the residual at k=3 is *exactly* what k=6 captures. If yes (|ρ|>0.95), the random effect's content beyond k=3's a_final is the very thing the paper's B3 covariate adds.


## Section 5: Recommendation for Paper 1.5


Probe '1: prefix_steps(r)' is the seed of Paper 1.5. The next step is to verify the candidate at higher modular resolutions (k=6, 7, 8, 9) and at larger N where the posterior noise floor drops, then derive the structural mechanism if the empirical match continues to hold.


---


*Inputs: B1 HMC fit at N=10K (8000 train / 2000 test), seed 20260509. All chain CSVs at C:/Collatz/experiments_output/paper1_hmc_n10k/hmc_B1/. Smoke test only — 4-point verdict, not paper material.*
