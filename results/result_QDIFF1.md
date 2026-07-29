# RESULT — QDIFF-1 (Move 1 foundation gate): the finite q-difference truncations EXIST, factor as clean D·(I+N), and NEST; the renewal is a pro-object of growing finite modules, not a fixed one. FOUNDATION HOLDS (pro sense). One spec correction (2026-07-28)

**Probe:** `probes/probe_qdiff1.py`, exact `Fraction` throughout. Foundation gate for Move 1 (realize the Syracuse renewal as a pro-object in the Tannakian category of q-difference modules). Data: exact `a_i=3^i R_e^(i)(2)`, `b_i=3^i R_e^(i)(0)` (i=1..5, MAHLER-certified `exact_Re`), `T_i=4a_i−b_i`, `Λ_i` from the certified S-ladder JSON. Substitution `z↦z³`, exponents `d_i=2·3^{i−1}`.

## Verdict
**Move 1's foundation HOLDS in the pro sense**, with the structure sharper than argued: the finite truncations exist, factor as a **clean `D·(I+N)`** with a *level-independent* unipotent integrator, and **NEST** (the pre-registered genuine uncertainty — resolved positive). But QD-B as written is vacuous — corrected below.

## The data (exact; `4a−b == T_cert` at every level)
```
   i     a_i                b_i           T_i=4a−b        Λ_i          log2 den(a_i)
   1     2/21               1/7           5/21           −2/21           4.4
   2     20614/203889       35401/203889  15685/67963    −1490/203889   17.6
   3     …/(62.5-bit den)   …             71597…/…       284995789764   62.5
   4     …                  …             …              347914009142   203.3
   5     …                  …             …              307926927625   629.3
```
`den(a_i)` bits `4.4→17.6→62.5→203→629`, ratios `4.0→3.55→3.25→3.10 → 3` — the doubly-exponential MAHLER rate.

## The structure (built in, then verified exactly)
- **The `+1` IS the unipotent integrator (N).** `T_i = T_{i−1} + Λ_i` holds **EXACTLY** for all i=1..7 (`T_0 = S_1/2 = 1/3` the base; `S_{i+1}=2T_i` banked). On `(T, Λ-source)` this is the Jordan block `N=[[1,1],[0,0]]`, `N²=0` — and it is **the same block at every level** (level-independent).
- **`D` is the `𝔾_m`/scaling.** `a_i,b_i,T_i` carry the `3^i` homogeneity = the `z↦z³` eigenvalue. So `A_1 = D·(I+N)`: `D` = multiply-by-3, `N` = the integrator — the argued `M=D(I+N)` shape is **present at finite level** (QD-D).
- **The growth source is Λ.** `Λ_i` has **no finite rational recurrence** (R27-A reproduced: none for L=1..3 on 7 exact terms). The module dimension grows *only* through the Λ-source coordinate — one new coordinate per level (dim increment = 1 = the unipotent-index increment).

## QD-B' — extrapolation (the real bite): no fixed finite module
A fixed affine map `V_{i+1}=M·V_i+c` on `V_i=(a_i,b_i)` (dim-2, 6 unknowns) fit *exactly* from transitions `1→2,2→3,3→4`, then asked to predict **held-out** `V_5`: **MISSES** (rel-err `+8.6e-4`). And rigorously **at any dimension**: `den(a_i)` is doubly-exponential, but a fixed constant-rational map gives at most single-exponential denominators `(fixed lcm)^i` — so **no fixed finite-dim constant-rational module exists** (the MAHLER denominator proof, applied to the vector). The truncation is a *good but imperfect* finite approximation (the near-miss), and the dimension **must grow** — this is the pro-object, not a fixed module.

## QD-C — nesting (the pre-registered genuine question): NESTS
Because the unipotent integrator `N` is **level-independent**, `A_{r+1}` restricts to `A_r` on the `(T, integrated-Λ)` block; the only growth is the Λ-source coordinate appended each level. So `lim← A_r` is well-defined, graded by Λ-level, dim increment constant = 1. **The tower nests cleanly.**

## ⚠️ Spec correction (flagged per "believe the gate, but the gate must bite")
**QD-B "does solving `A_r` reproduce the exact `a_i,b_i`" is VACUOUS for monomial-supported generating functions.** The `a_i` live only on exponents `d_i=2·3^{i−1}`; the naive realization `A_1(z)=[[z^{d_2−d_1},0],[0,1]]=[[z^4,0],[0,1]]` reproduces `a_1↦a_2` by pure `z`-power scaling and **encodes nothing** — yet "passes." So reproduction cannot be the gate. The content is **(B') extrapolation + (C) nesting + (D) the `D·(I+N)` shape**, which is what this probe tests. This is exactly the trap the spec warned of (the coincident 3's), and it bites the reproduction check itself — reproduction is necessary but not discriminating.

## Net — for Move 1
- **Foundation HOLDS (pro sense):** the renewal is the inverse limit of finite q-difference modules; the **unipotent part is finite, clean (`I+N`, integrator), and nesting**; the dimension growth is isolated in the Λ-source = the infinite Mahler depth. It is **not** a fixed finite q-difference module (MAHLER), so the Galois computation is on the **pro-object / associated graded**, not a single module — consistent with the whole arc (structure at the graded, value walled above).
- **Redirect, minor:** QD-B should be extrapolation+nesting, not reproduction (banked here so the skeleton uses the right gate).
- **Not at stake:** S_∞≈0.475 (floor 0.473177), MAHLER, MIRROR, PHYDRA, two-walls, R1–R30. Structure only; 7/15 excluded regardless.
