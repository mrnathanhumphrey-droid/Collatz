# Probe A — dominant singular vector support inspection

Inspecting u_1 (left singular vector in W_{k+1} embedded back into V_{k+1} natural basis) and v_1 (right singular vector in V_k natural basis) for R_k = P_{W_{k+1}} ∘ K_{k+1} ∘ L_k at k = 5, 6, 7.

## Per-k summary

| k | sigma_1 | u_1 nz/n_{k+1} | v_1 nz/n_k | top-10% E(u_1) | top-10% E(v_1) | cos(u_1, pi_{k+1}) | cos(v_1, pi_k) |
|---|---|---|---|---|---|---|---|
| 5 | 0.670554 | 486/486 | 161/162 | 0.4359 | 0.1960 | +0.0277 | +0.0253 |
| 6 | 0.670615 | 1458/1458 | 485/486 | 0.4384 | 0.1960 | -0.0162 | -0.0147 |
| 7 | 0.670622 | 4374/4374 | 1458/1458 | 0.4400 | 0.1973 | +0.0002 | +0.0002 |

## u_1 mod-3 residue structure across k

Energy fraction of u_1 in V_{k+1} on residues r mod 3 ∈ {1, 2}:

| k | r=1 mod 3 | r=2 mod 3 |
|---|---|---|
| 5 | 0.200000 | 0.800000 |
| 6 | 0.200000 | 0.800000 |
| 7 | 0.200000 | 0.800000 |

## u_1 mod-9 residue structure

Energy fraction of u_1 grouped by r mod 9 (only coprime residues r ∈ {1,2,4,5,7,8}):

| k | r=1 | r=2 | r=4 | r=5 | r=7 | r=8 |
|---|---|---|---|---|---|---|
| 5 | 0.112053 | 0.448213 | 0.080944 | 0.028013 | 0.007003 | 0.323774 |
| 6 | 0.112046 | 0.448183 | 0.080951 | 0.028011 | 0.007003 | 0.323806 |
| 7 | 0.112045 | 0.448180 | 0.080952 | 0.028011 | 0.007003 | 0.323809 |

## v_1 mod-3 residue structure (input direction in V_k)

| k | r=1 mod 3 | r=2 mod 3 |
|---|---|---|
| 5 | 0.500000 | 0.500000 |
| 6 | 0.500000 | 0.500000 |
| 7 | 0.500000 | 0.500000 |

## Stationary alignment

cos(u_1, pi_{k+1}) and cos(v_1, pi_k): if either is near ±1, the dominant SVD direction is the stationary itself; if near 0, it's orthogonal (a non-trivial deviation mode).

- k=5: cos(u_1, pi_6) = +0.027669, cos(v_1, pi_5) = +0.025346
- k=6: cos(u_1, pi_7) = -0.016162, cos(v_1, pi_6) = -0.014702
- k=7: cos(u_1, pi_8) = +0.000248, cos(v_1, pi_7) = +0.000224

The cosines decrease toward 0 as k grows (0.028 → 0.016 → 0.0002 for u_1 vs π).
The dominant SVD direction is therefore essentially **orthogonal to the
stationary** at large k — a non-trivial deviation mode.

## Structural findings

### 1. v_1 is mod-3-balanced and residue-class-uniform on coprime classes

Right singular vector v_1 ∈ V_k has identical energy fractions across all
6 coprime residue classes mod 9 (each = 1/6 to 4 decimal places), and across
all 18 coprime classes mod 27 (each = 1/18). v_1 is "residue-class energy
uniform" at the input side — element-wise non-uniform within each residue
class (||v||_∞ ≠ uniform amplitude), but per-class energy is constant.

Holds at all three k values tested.

### 2. u_1 mod-3 energy split is exactly (1/5, 4/5)

Energy ratio 1:4 between r=1 and r=2 mod 3 ⇒ amplitude ratio 1:2.

This matches the asymptotic Tao-chain stationary's mod-3 mass ratio
(b, c) = (1/3, 2/3) — amplitude 2× on r=2 — exactly. The dominant forcing
direction inherits the stationary's mod-3 asymmetry as a **squared** ratio.

### 3. u_1 mod-9 distribution is k-invariant

Six numbers (r=1: 0.1120, r=2: 0.4482, r=4: 0.0810, r=5: 0.0280, r=7: 0.0070,
r=8: 0.3238), stable to 4 decimal places across k = 5, 6, 7. Conditional on
r mod 3:
- r mod 3 = 1 (residues 1, 4, 7 mod 9): conditional dist (0.560, 0.405, 0.035)
- r mod 3 = 2 (residues 2, 5, 8 mod 9): conditional dist (0.560, 0.035, 0.405)

The two conditionals are equal up to a swap of the second and third entries
— a reflective structure between the two coprime-mod-3 cosets.

### 4. u_1 has uniform energy at moduli ≥ 27

Energy at mod 27 = (mod-9 energy) / 3 to machine precision; same at mod 81.
**u_1's residue-class structure exists only at modulus 9**; finer
subdivisions partition energy uniformly across children. The dominant
forcing direction lives in a 6-parameter family (one parameter per coprime
mod-9 class), uniformly extended at higher resolution.

### Interpretation

The dominant V → W forcing direction takes a **mod-3-balanced input** and
produces a **mod-9-structured output** whose mod-3 marginal matches the
asymptotic Tao stationary mass ratio. The 6-number mod-9 signature is
k-invariant. R_k's dominant action: "lift balanced input, project the
K_{k+1} image onto a universal mod-9 deviation pattern."

This is consistent with the rank-2n_k/3 structure of R_k: the kernel
(dim n_k/3) likely consists of inputs whose mod-3 marginal already
matches the stationary's (1/3, 2/3) — they lift+evolve into L_k(V_k)
and produce zero in W_{k+1}. The non-kernel directions are "deviations
from the stationary mod-3 ratio" and they all map into the same
6-parameter mod-9 family up to k-stability.

The σ_1 ≈ 0.6706 number doesn't yet have a closed form, but the dominant
input/output directions are now characterized in residue-class language.
