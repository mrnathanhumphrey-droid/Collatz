# RESULT — GARSIA: ν is 3-adically DIMENSION-1 / NOT Garsia-singular (h=1) — the singular⟹irrational lean is falsified (2026-07-27)

> ⚠️**GUARDRAIL (read first):** the result is **dimension 1 / not Garsia-singular**, which is **NOT** the same as
> **absolutely continuous**. Dimension 1 is *necessary but not sufficient* for a.c. — there exist dimension-1 singular
> measures. This is exactly the gap the Bernoulli-convolution field lived in for decades (Solomyak: a.e. a.c.; Hochman:
> dimension 1 outside a zero-dimensional exceptional set — two different theorems). So everywhere below, read "not
> singular in the Garsia/dimension sense," **never** "absolutely continuous." The distinction is load-bearing: dimension 1
> **removes the structural (dimension-drop) reason S∞ would be irrational**; it does **not** hand you rationality, and a
> dimension-1 measure can still have a transcendental Plancherel mass.

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

**Verdict: ν has Garsia entropy dimension 1 — it is 3-adically DIMENSION-1 / NOT Garsia-singular. (Not "a.c." — see guardrail.)**

## Theory: why dimension-1 is the RIGHT answer (and DENOM's Pisot lean was the wrong slice)
`result_DENOM_OBSTRUCTION.md` §Reduction leaned "1/λ=2 is a Pisot base ⟹ singular is the default (⟹ S∞ generically
irrational)." That is the **wrong slice** of the Bernoulli-convolution theorem. Erdős/Pisot-singularity is for **non-integer**
Pisot λ∈(1,2) (golden ratio etc.). Base exactly **2** is the **a.c. edge**, not the Pisot interior (λ=2 → uniform). So the
singular-by-default expectation was a **misapplication**, and the support count then proves full dimension outright.
Breuillard–Varjú: Mahler `M(2)=2` ⟹ no base-driven entropy defect; any defect would come from **digit-map overlap**, and the
measured `dH_r/log3 → 1` says the overlap is asymptotically entropy-negligible. **Measurement and corrected theory agree:
ν is dimension-1 / not Garsia-singular.** (Still NOT a claim of absolute continuity — guardrail.)

## What this DOES and does NOT settle (honest, at the right size)
- **DOES kill the singular (dimension-drop) mechanism.** Full support (β=3, exact) + entropy dimension 1 ⟹ **no Garsia
  singularity.** So DENOM's "ν singular ⟹ S∞ irrational" **never fires — there is no singularity to fire it.** The "generically
  irrational, 7/15 wrong-in-kind" lean is **retracted**; it was a scaffold built on a misread theorem.
- **Does NOT establish absolute continuity, and does NOT prove S∞ rational.** Dimension 1 ≠ a.c.; and even a.c. can have an
  irrational L² / Plancherel mass. h=1 removes the *structural argument* for irrationality; it decides nothing about the value's
  arithmetic type.
- **The flip removes a wrong argument in BOTH directions** (the good kind of correction): the type of ν now gives **no
  obstruction to rationality either way** — 7/15 is neither excluded by singularity (there is none) nor supported by it. The
  rational-vs-irrational question comes **off the measure-type argument entirely** and goes back to what it was before the Pisot
  detour: an open question about a specific limit, value pinned at **≈0.475**, with **7/15 excluded only by the finite-depth /
  tail arithmetic** (the 47–150× exclusion, barring a hidden `ρ₃ > 0.999` real mode).
- If a direction must be named at all: ν being dimension-1/maximal-entropy is *mildly* consistent with S∞ algebraic (Streck's
  rigid-locus), but that is a lead to test, not a conclusion.

## CF corroboration (`probe_cf_cancel.py`, low priority — as flagged)
Wilson's honest integer-level test (the half-integer strobe was correctly refused: the tower is 3-to-1 with no intermediate
group, `M^{1/2}` is a different operator whose `(−1)ⁿ` would be a branch-cut artifact). Rotation number
`α = log3/log2 − 1 = 0.585`. Convergent **denominators** `k_k` (return times) predict phase nodes: **i=12** is the single
closest-to-integer level (dist 0.020, a convergent denom), i=5 second — **confirms the rotation number `log3/log2` from the
integer side.** The DMT **numerator** levels {8,19} are mixed on residual magnitude (i=19,12 below-mean, i=8 above). Weak
positive, 13 points, decides nothing — it re-states the already-confirmed `2π/log2 = 9.06` period. As flagged.

## Net
- **Support β=3 exact (algebra); entropy dimension h→1 (measured 0.987@r=14, rising; tol-robust; heuristic-anchored).**
- **ν is 3-adically dimension-1 / NOT Garsia-singular** (NOT a claim of a.c. — guardrail). The DENOM singular-lean is
  **falsified**; the arithmetic type now obstructs 7/15 in **neither** direction — it is off the table as an argument.
- Value stands (S∞≈0.475, floor 0.4732; **7/15 excluded only by the finite-depth/tail arithmetic**, barring `ρ₃>0.999`).
  Rational-vs-irrational is **genuinely open**, back to a plain question about a specific limit, no measure-type lean.
- **Two follow-ups (Wilson's pen; neither urgent):**
  1. **Upgrade h=1 from measured to proven.** Heuristic: per-level increment supplies `2log2=1.386` vs `log3=1.099`. ⚠️gate-catch:
     the *finite-r* increment is `dH_r/log3 = 0.987 < 1` (below log3), so a literal "increment ≥ log3 each level" is **false** at
     finite r — the proof needs the increment *→ log3* (a limiting/lower-bound-that-saturates argument), not a per-level bound.
     Given the primitive-root full-support fact, a clean asymptotic lower bound looks reachable and would make it a **lemma**.
  2. **Mahler cross-check (Breuillard–Varjú).** Dimension 1 ⟺ maximal Garsia entropy ⟺ trivial Mahler defect; base 2 has
     `M(2)=2`, entropy scale `log2`. The den(N) `(2^M−1)` structure is the arithmetic shadow — full entropy predicts a specific
     simple Mahler value. If it checks, the entropy result is corroborated **from the norm side already proven** (independent).
- **Not at stake:** P6D–P6K identities, `S_{i+1}=2T_i`, i=20 no-crossing, R1–R30, the DENOM *denominator theorem* itself
  (only its singular-lean *corollary* is corrected).
