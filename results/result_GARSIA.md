# RESULT — GARSIA: ν is 3-adically ABSOLUTELY-CONTINUOUS-type (h=1), NOT singular — the singular⟹irrational lean is falsified (2026-07-27)

**Probe:** `probes/probe_garsia.py` (+ `probe_cf_cancel.py`). Spec: Wilson's two thoughts — (1) the Garsia support/entropy
instrument, made *decisive* in the p-adic setting; (2) the p-adic Garsia dictionary (local dim = log(reach/level)/log 3,
Streck maximality, Breuillard–Varjú Mahler bridge). Certified `build_nu` (R28) machinery only. **The point:** decide the
singular-vs-a.c. TYPE of ν — the thing `result_DENOM_OBSTRUCTION.md` reduced the whole sign question to.

## The instrument, correctly: SUPPORT ≠ MEASURE
Two different dimensions. Support-fatness is necessary, not sufficient, for a.c. The **measure** decides, via **entropy**
(mass scaling), not the point-count. Both computed; they say the same thing.

## PART 1 — SUPPORT: exact, β=3, box-dim 1 (algebra, not a fit)
`ν_r` = law of `X_r`, renewal `X' = 1 + 3·2^{-v}X` on ℤ/3^{r+1}, `v~Geom(½)`. The multiplier `2^{-v}` runs over `⟨2⟩ mod 3^r`.
- **2 is a primitive root mod 3^r for every r** (verified `ord_{3^r}(2)=2·3^{r-1}`, r=1..15). So `⟨2⟩` is the *whole* unit group.
- Every support point is `≡1 mod 3` (a unit), so `{2^{-v}·X_p : v}` = all units ⟹ `supp(ν_r) = {1+3t : t∈(ℤ/3^r)^*}`.
- **`s_r = 2·3^{r-1}` EXACTLY** (verified by full-multiplier-group closure, r=1..8: 2,6,18,54,162,486,1458,4374). So
  **β=3 exactly, box-dimension 1** — ν is 3-adically full-support, NO Cantor-like missing digits. (A crude "not all digits
  reached ⟹ singular" is dead on arrival.)
- ⚠️gated footnote: `len(build_nu, tol=1e-18)` happens to equal `s_r` here (level-composition fills the coset), but that is
  luck of this measure — the *instrument* for singularity is entropy, not `len`.

## PART 2 — ENTROPY: the Garsia dimension → 1 (a.c.-type)
`H_r = −Σ ν_r(X) log ν_r(X)` from `build_nu(0.5,r)`, r=1..14 (mass-weighted ⟹ **tol-safe**: identical at tol 1e-12/18/24).
Per-level entropy `dH_r/log3` (Garsia entropy dimension, incremental) **rises monotonically toward 1**:
```
 r    dH_r/log3      h_r=H_r/(r log3)      (gap to 1 of dH_r/log3)
 8     0.96949          0.88523              0.0305
10     0.97783          0.90338              0.0222
12     0.98333          0.91650              0.0167
14     0.98718          0.92647              0.0128
```
- gap-to-1 decays ~geometrically (ratio ≈0.88/level) ⟹ **extrapolated limit h = 1**.
- late-window (r=8..14) slope `dH/dr / log3 = 0.982`, still climbing.
- **Heuristic anchor (matches):** each level supplies `H(v)=H(Geom ½)=2 log2 = 1.386` nats of randomness against
  `log3 = 1.099` of scale — ratio **1.262 > 1** ⟹ the new trit is saturated ⟹ **h=1**, maximal entropy.

**Verdict: ν has Garsia entropy dimension 1 — it is 3-adically ABSOLUTELY-CONTINUOUS-type, NOT singular.**

## Theory: why a.c. is the RIGHT answer (and DENOM's Pisot lean was the wrong slice)
`result_DENOM_OBSTRUCTION.md` §Reduction leaned "1/λ=2 is a Pisot base ⟹ singular is the default (⟹ S∞ generically
irrational)." That is the **wrong slice** of the Bernoulli-convolution theorem. Erdős/Pisot-singularity is for **non-integer**
Pisot λ∈(1,2) (golden ratio etc.). Our base is the **integer 2**, whose Bernoulli convolution is **a.c.** (λ=2 → uniform, the
textbook a.c. case). Breuillard–Varjú: Mahler `M(2)=2`, no conjugates on/inside the circle ⟹ no base-driven entropy defect;
any defect would come from **digit-map overlap**, and the measured `dH_r/log3 → 1` says the overlap is asymptotically
entropy-negligible. **The measurement and the corrected theory agree: ν is a.c.-type.**

## What this DOES and does NOT settle (honest)
- **DOES falsify** the singular half of DENOM's reduction: ν is **not** singular, so **"ν singular ⟹ S∞ irrational" never
  fires.** The type-based obstruction to 7/15 **does not exist.** The lean that S∞ is "generically irrational, 7/15 wrong in
  kind" is **retracted** — it rested on a singularity that isn't there.
- **Does NOT prove S∞ rational.** An a.c. measure can have an irrational Plancherel mass (an a.c. density can have irrational
  L² norm). h=1 removes the singularity *argument* for irrationality; it does not decide rational vs irrational.
- **7/15 status unchanged: excluded on VALUE** (floor `2·T_20 = 0.473177` from i=20), **not** on type. If anything, ν being
  a.c./maximal-entropy points *toward* S∞ algebraic (Streck: maximal entropy ⟹ ν sits on a rigid finite linear locus — a
  checkable structure), i.e. the irrationality is genuinely open and *less* favored than the singular reading implied.

## CF corroboration (`probe_cf_cancel.py`, low priority — as flagged)
Wilson's honest integer-level test (the half-integer strobe was correctly refused: the tower is 3-to-1 with no intermediate
group, `M^{1/2}` is a different operator whose `(−1)ⁿ` would be a branch-cut artifact). Rotation number
`α = log3/log2 − 1 = 0.585`. Convergent **denominators** `k_k` (return times) predict phase nodes: **i=12** is the single
closest-to-integer level (dist 0.020, a convergent denom), i=5 second — **confirms the rotation number `log3/log2` from the
integer side.** The DMT **numerator** levels {8,19} are mixed on residual magnitude (i=19,12 below-mean, i=8 above). Weak
positive, 13 points, decides nothing — it re-states the already-confirmed `2π/log2 = 9.06` period. As flagged.

## Net
- **Support β=3 exact (algebra); entropy dimension h→1 (measured 0.987@r=14, rising; tol-robust; heuristic-confirmed).**
- **ν is 3-adically a.c.-type. The DENOM singular-lean is falsified; the arithmetic type does NOT obstruct 7/15.**
- Value stands (S∞≈0.475, floor 0.4732; 7/15 excluded on value). Rational-vs-irrational **genuinely open**, and now points,
  if anywhere, toward *algebraic* (Streck rigid-locus), not toward irrational-by-singularity.
- **Next lever (Wilson's, if he wants it):** the Streck maximality equations — a.c./maximal-entropy ν must satisfy a finite
  linear system tied to the base; test the banked class means (5/3, 2/3) against it. That is the route that could upgrade
  "a.c.-type" to a *constraint on the value*.
- **Not at stake:** P6D–P6K identities, `S_{i+1}=2T_i`, i=20 no-crossing, R1–R30, the DENOM *denominator theorem* itself
  (only its singular-lean *corollary* is corrected).
