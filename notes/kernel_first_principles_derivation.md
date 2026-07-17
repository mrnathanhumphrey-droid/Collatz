# Result 77 — Kernel first-principles derivation: outcome (γ); local dynamics match (per-state Pearson 0.96), but Perron eigvec under within-state-uniform-m loses the trajectory measure's structure

**Date:** 2026-05-03. Derives R60's 1024-state kernel K from first principles
(arithmetic-deterministic v rule + log_2 random walk under uniform-m-within-state)
and compares to R60's empirical K.

**Verdict (γ) with substantive insight:** the local transition structure
matches well (per-state row Pearson 0.96 median) but the global Perron
eigenvector under the uniform-within-state assumption is essentially
uniform — it does NOT recover D_avg.

**This reveals what R60 was actually identifying:** the trajectory measure
is NOT determined by Markov local dynamics on (r, b). The Markov dynamics
are too uniform; the trajectory measure's non-uniformity comes from
the empirical visit-weighting pattern, which is a separate structural fact.

## Construction (first-principles K_derived)

For each state s = (r mod 32, b) with r ∈ {1, 3, ..., 31}, b ∈ {0, ..., 63}:

1. **Conditional v distribution P(v | r mod 32):** from arithmetic determinism
   (R45/R68), v is determined by m mod 32 for 15 of 16 residues. r=21 is the
   singular boundary; v ≥ 5 with Geom(½) tail by bit-recursion.

2. **Conditional next-residue P(r' | r mod 32, v):** for r ≠ 21, depends on
   m mod 64 (or higher) — typically 2 values for v=1, 4 for v=2, 8 for v=3,
   16 for v=4, each with equal probability.

3. **Conditional Δb under uniform-m-within-state:**
   ```
     b' = b + Δb  with  Δb ∈ {−v + 1, −v + 2}, probs (1/3, 2/3)
   ```
   Derived from u = m / 2^b uniform in [1, 2) and the breakpoint u = 4/3 from
   log₂(3u) crossing integer values.

4. **K_derived[s, s'] = Σ_{(v, r')} P(v, r' | r mod 32) · P(b' − b | v)**

K_derived has 22,175 nonzero entries (2.1% sparsity), row sums ≈ 1.

## Per-state transition agreement (HIGH)

```
Per-state row Pearson(K_derived[s], K_emp[s]) over 626 common states:
  mean   = 0.890
  median = 0.958
  min    = 0.035
```

Sample state checks confirm structural agreement:
```
State (r=3, b=10):
  DERIVED: (5|11):0.333, (21|11):0.333, (5|10):0.167, (21|10):0.167
  EMPIRIC: (21|11):0.365, (5|11):0.330, (5|10):0.181, (21|10):0.124
  → SAME 4 transitions, comparable probabilities

State (r=15, b=20):
  DERIVED: (23|21):0.334, (7|21):0.333, (23|20):0.167, (7|20):0.167
  EMPIRIC: (23|21):0.300, (7|21):0.290, (23|20):0.210, (7|20):0.202
  → SAME 4 transitions, slight skew toward Δb=0 in empirical
```

The 4-transition pattern (2 next-residues × 2 size-bin shifts) matches
exactly. Empirical proportions deviate by 5-10% from the 1/3 ↔ 2/3 split,
consistent with finite-N sampling AND non-uniform-m-within-state effects.

## Failure case: state (r=5, b=8)

```
State (r=5, b=8):
  DERIVED: (23|6):0.042, (5|6):0.042, (3|6):0.042, ...  — 16 transitions ~uniform
  EMPIRIC: (29|5):0.799, (27|6):0.153, (9|6):0.033, ...  — concentrated on 2-3 states
```

For r=5 (v=4 deterministic), my derived K spreads over 16 next-residues
uniformly. The empirical K concentrates on (29|5) at 80%. This means the
m values that ARE visited at (r=5, b=8) on the trajectory measure are
NOT uniformly distributed mod 32^higher — they are heavily concentrated
at a specific higher-bit residue.

This is exactly the trajectory-measure structure that R60 was identifying
implicitly through visit-frequency weighting.

## Perron eigenvector recovery: FAILS

```
   r       rho_pred    D_pred     D_avg     diff
   1       0.062397    0.983      1.609     −0.626
   3       0.062553    1.002      1.236     −0.234
   5       0.062492    0.992      1.864     −0.872
  …       (all rho_pred ≈ 0.0625, all D_pred ≈ 1.0)
  31       0.062664    1.014      0.767     +0.247

  Total |D_pred − D_avg| = 5.40
  Pearson ρ = −0.57
  
  Reference: R60 empirical → total_dev = 3.40, ρ = +0.80
             trivial null  → total_dev = 4.72
```

K_derived's Perron eigenvector under uniform-within-state is essentially
uniform on (Z/32Z)*: π[r] ≈ 1/16 for every odd r. Marginalized to
D_pred(r) = π[r] / π_32[r] ≈ 1.0 for every r. **It does NOT recover the
trajectory measure D_avg.**

## What this means

R60's empirical kernel and its Perron eigenvector decompose as:

```
  K_emp ≈ K_dynamics × W_visit
```

where:
- K_dynamics encodes the local Syracuse transition structure (matches K_derived,
  per-state Pearson 0.96)
- W_visit encodes how often each (r, b) cell is visited by the actual
  trajectory measure — this is the absent piece in K_derived

The trajectory measure D_avg = Perron eigenvector of K_emp ≠ Perron
eigenvector of K_derived because the eigenvector calculation is sensitive
to W_visit weights, not just to K_dynamics.

**R60's "size-stratified Markov identifies D_avg" was technically correct
but obscured a key fact:** the identification works because R60's K is
empirically estimated FROM trajectory orbits, so the visit frequencies are
already baked in. A genuinely first-principles K (without visit-frequency
weighting) doesn't reproduce D_avg.

## Per brief outcomes

| Outcome | Status |
|---|---|
| (α) K_derived ≈ K_empirical structurally; rigorous identification | **REJECTED** for Perron eigvec |
| (β) Some structure matches, residuals identifiable | **APPLIES** for per-state rows; not for global Perron |
| (γ) K_derived differs systematically; assumption missing structure | **PRIMARY** for the eigvec recovery test |

## What's the missing structure?

The trajectory measure has within-state non-uniformity in m mod 2^k for
k > 5 (where 2^5 = 32 is the residue modulus). Specifically, R65 showed
the trajectory measure is mod-3 non-uniform (Z_3-Bohr concentration);
this propagates to non-uniform m mod 2^k structure within each (r, b)
state via the Syracuse map's intertwining of mod 2 and mod 3.

To derive K from first principles AND recover D_avg, would need:
1. The within-state distribution of m mod 2^k for high k (NOT uniform)
2. This distribution itself derives from the SAME trajectory measure
   we're trying to identify
3. → fixed-point problem: K_derived = K(D_avg) and D_avg = Perron(K_derived)

This circularity is informative: the trajectory measure is NOT a free
output of Markov dynamics on (r, b) — it requires self-consistency with
its own within-state profile.

## Implications for framework synthesis

**Walk back:** "R60 size-stratified Markov derives D_avg from first principles"
is not what the empirical kernel does. R60's empirical kernel encodes
visit frequencies, not pure dynamics.

**Reframe:** R60 is a **structural identification** of the form
  D_avg = Perron eigvec of (K_dynamics × W_visit)
where K_dynamics is derivable from arithmetic but W_visit is empirically
extracted from trajectory orbits. Both pieces are needed; neither alone
suffices.

**For external correspondence:** the framing "trajectory measure ↔ Perron
eigenvector of empirical Markov kernel" should explicitly note that the
kernel is not first-principles derivable without circularly using the
trajectory measure itself.

## Connection to other findings

- **R65 3-adic specificity:** the within-state non-uniformity that K_derived
  misses IS the 3-adic structure. Connecting these: K_derived would match
  R60 if the within-state distribution were the trajectory measure's
  3-adic Bohr-set decomposition.
- **R72/R73 distribution characterization:** the per-a Berry-Tabor
  exponential is an asymptotic property; here we see that the trajectory
  measure's STRUCTURE (not just distribution) requires self-consistency
  with its own profile.
- **R74 lifting recursion:** S_{k+1} = 3^(k+1)·||d||² is provable
  (algebraic). K_derived shows that the LOCAL Markov structure is also
  provable. But neither closes the trajectory measure's CONSTANT (the
  c = 7/45 that gives 7/15) without an additional self-consistency step.

## What would close outcome (α)

A first-principles derivation that recovers D_avg would need:
1. Show that the trajectory measure satisfies a fixed-point equation
   D_avg = Perron(K(D_avg))
2. Solve the fixed-point equation in closed form

Step 1 is straightforward (trajectory measure IS the visit measure of
Collatz orbits, and the Markov dynamics on (r, b) aggregate to its
identity). Step 2 is harder — it's the same difficulty as the open
c = 7/45 derivation in R74's analytical step.

Possibly the same algebraic identity that closes R74 also closes this.

## Files

- `kernel_first_principles_v2.py` — corrected derivation (proper joint
  P(v, r' | r mod 32))
- `kernel_first_principles_v2_log.txt` — full output
- `K_derived_v2.npz` — sparse 1024×1024 derived kernel
- `D_predicted_derived_v2.csv` — Perron marginal vs D_avg
- `K_derived_vs_empirical_v2.csv` — transition-by-transition comparison
- `kernel_first_principles_derivation.md` — this writeup

(`kernel_first_principles.py` v1 had a bug treating r' as deterministic
from (r, v); kept for reference.)

## Concrete next moves

1. **Within-state distribution measurement**: empirically compute the
   trajectory measure's distribution of m mod 2^k for k > 5 within each
   (r mod 32, b) cell. This gives W_visit explicitly.
2. **Fixed-point K_derived = K(W_visit)**: build K with W_visit (not
   uniform) and check whether Perron(K) recovers D_avg.
3. **Connection to R74's c = 7/45**: investigate whether the fixed-point
   structure here is the same algebraic identity that gives the leading
   ||d||² coefficient.
4. **Reframe R60 in framework synthesis**: "D_avg = Perron of the
   empirical kernel that already encodes visit weights" — not "first-
   principles derivation."
