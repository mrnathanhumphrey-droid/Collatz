# RESULT — PROBE CARRYCOV: the carry-covariance conjecture (2026-07-26)

**Probe:** `probes/probe_carrycov.py`. **Question (Wilson's pen):** does the fluctuation profile get
*carried along* by the carry — `dpi_{4x} ≈ T_{-c(x)} dpi_x` — so that `q_r(1)-1/3` becomes a **variance**
`E_mu ||dpi||^2 > 0` (the right SHAPE, the m=0 degeneration transported to m=1 by the carry twist)?

Machinery reused verbatim from `probe_carrylemma.gates_float`: `M=3^r`; `x`=class mod M; `c=floor(4x/M)∈{0,1,2,3}`
(DF Prop 5.1 carry); `L=4x mod M`; `dpi_x(d)=nu_hi(x+dM)/nu_low(x)-1/3` (3-vector); `mu(x) ∝ nu(x)nu(L)`.
Rotations mod 3, `(T_a f)(d)=f(d+a)`. Exact split `dpi_{4x} = T_{-c}dpi_x + rho_x`, `rho = dpL - roll(dpx, c%3)`.

## VERDICT: conjecture FAILS the sufficient bound; positivity is NOT a variance. Honest negative — plus one real gain.

### 1. The defect is NOT small — `D ≈ 1.39`, rock-stable, `> 1` at every level

Defect ratio `D = sqrt(E_mu||rho||^2)/sqrt(E_mu||dpi_x||^2)` (mu-weighted RMS):

| r | q-1/3 | VAR=E‖dpi‖² | RES=E⟨dpi,T_c rho⟩ | D (RMS) | ptw | VAR(1-D) | (q-1/3)/VAR |
|---|-------|-------|-------|--------|------|--------|-------|
| 4  | +1.902e-3 | +6.240e-2 | -6.050e-2 | 1.3925 | 1.44 | -2.45e-2 | 0.0305 |
| 8  | +1.011e-3 | +3.225e-2 | -3.124e-2 | 1.3873 | 1.60 | -1.25e-2 | 0.0313 |
| 12 | +6.396e-4 | +1.770e-2 | -1.706e-2 | 1.3902 | 1.82 | -6.91e-3 | 0.0361 |
| 16 | +4.179e-4 | +1.032e-2 | -9.903e-3 | 1.3921 | 2.08 | -4.05e-3 | 0.0405 |

`D` sequence r=4..16: `1.392 1.383 1.382 1.385 1.387 1.389 1.389 1.390 1.390 1.391 1.391 1.391 1.392` — dead
flat, `D[16]-D[12]=+0.0019`. **The carry does not carry the profile along: `dpL` and `T_{-c}dpx` are essentially
UNcorrelated** — `D≈1.39` sits just under `√2≈1.414`, the value for two independent equal-variance mean-zero
vectors. Only a razor of positive correlation pulls it below √2. The literal "carried along, defect small"
picture (which would need `D≪1`) is quantitatively false. `ptw` (pointwise mean defect/‖dpi‖) even RISES 1.44→2.08.

### 2. Cauchy–Schwarz sufficient condition (`D<1`) not met, never approached

`VAR(1-D) < 0` at every level (≈ `-4.0e-3` at r16), ~**10× the wrong way** relative to `q-1/3 = +4.2e-4`.
`slack/(q-1/3)` r=12..16 = `-10.8, -10.3, -10.0, -9.75, -9.68`. The C–S bound is sign-blind here (far too lossy);
positivity survives ONLY through the **signed** RES term, not a variance lower bound.

### 3. The split is EXACT but CIRCULAR — no new leverage (Wilson's reading iii, sharpened to an identity)

`q_r(1)-1/3 = VAR + RES` holds exactly (Fractions r=4,5: `VAR+RES == q-1/3` True). The norm expansion gives the
sharp threshold in closed form:
```
E_mu||rho||^2 = VAR + VAR_L - 2(q_r(1)-1/3),   VAR_L := E_mu||dpi_L||^2
  =>  q_r(1)-1/3 = ( VAR + VAR_L - E_mu||rho||^2 ) / 2
  =>  "defect small"  (E||rho||^2 < VAR + VAR_L)  <=>  q_r(1)-1/3 > 0   IDENTICALLY.
```
At r16: `E‖rho‖²=2.000e-2`, `VAR+VAR_L=2.084e-2` (VAR=1.032e-2, VAR_L≈1.052e-2), margin `8.4e-4 = 2(q-1/3)` ✓.
So the defect ratio `D` encodes **nothing beyond `LEM/VAR`**: "the defect is small" is not evidence for the sign,
it is an algebraic restatement OF the sign. The covariance framing does not reduce the problem. Positivity stays a
fine near-cancellation `+1.03e-2 + (-9.90e-3) = +4.18e-4` (residual = 4% of VAR), living in the same signed
cross-term `⟨dpi, T_c rho⟩` as the carry lemma's `T1+T2`. No escape to a variance.

*(Mild note against "low-r accident": `(q-1/3)/VAR` RISES 0.031→0.041 over r — the residual is a slowly growing
fraction of the variance, so it is not washing out relative to VAR. Weak, one decade.)*

### 4. THE GAIN — deparitied excess rate confirms the transitional-0.89 flag (independent of the conjecture)

Wilson's corrected rate check: deparity first (`s_r=(e_r+e_{r+1})/2` kills the period-2 term), read at the TOP.
- raw excess two-step `(|e_r/e_{r-2}|)^.5`: r12..16 = `0.907 0.916 0.907 0.900 0.891`
- **deparitied `s_{r+1}/s_r`: r11..14 = `0.916 0.907 0.901 0.891` — smooth MONOTONE DOWN, no leveling.**

This **confirms and sharpens FOURCELL's transitional-rate flag.** The honest asymptotic excess rate is `≤0.89 and
still falling`; if it steepens toward the ~0.80 component rate, the tail sum shrinks and `S_∞` drops toward ~0.473.
**Consequence for the value: the `S_∞≈0.477` figure rests on the transitional 0.89 rate and should be read as an
UPPER estimate; the deparitied trend puts `S_∞` below 0.475.** The SIGN is untouched (`e_r>0` at every r, VAR
structurally dominant) — this is a correction to the *value*, not the *sign*.

## Net
- Carry-covariance conjecture **refuted** as a proof mechanism: defect `D≈1.39` (not small), C–S bound fails by
  40%, and the exact "defect < variances" threshold is algebraically **identical** to `q-1/3>0` — circular, zero
  new leverage. Positivity remains a signed near-cancellation in the same cross-term we already had.
- **The sign continues to live in the signed carry cross-term `⟨dpi, T_c rho⟩` (= CARRYLEMMA `T1+T2`).** The
  covariance route does not remove the signed object; it repackages it.
- **Real gain: the deparitied top-of-ladder excess rate (0.916→0.891, monotone down) confirms the 0.89 rate is
  transitional ⟹ `S_∞` estimate should sit below 0.475, not at 0.477.**

**Gates:** VAR+RES == q-1/3 exact (Fractions r=4,5); float split rel < 1e-6 all r=4..16; `build_nu(0.5,16)` (349s).
**Not at stake:** CHANNEL_ID (`d1=A_r(1)/S_r`) + CARRYLEMMA identities, R1–R30, R80–R82, all Thread-3 probes.
