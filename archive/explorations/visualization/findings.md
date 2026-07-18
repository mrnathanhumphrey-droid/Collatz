# qx+1 phase-space 3D visualization — findings

**Data:** 60,000 starting values (5,000 odd n per q) drawn uniformly from [3, 10⁶], for q ∈ {3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25}. Trajectory cap 10,000 steps or value > 2⁶² (int64 safe-cap). Generation: 5.5s.

**Convergence breakdown** (drives every other finding):

| q | converged | divergent | timeout |
|---|---|---|---|
| 3 | 5000 (100.0%) | 0 | 0 |
| 5 | 17 (0.3%) | 4945 (98.9%) | 38 (0.8%) |
| 7 | 0 | 5000 (100%) | 0 |
| 9 | 0 | 5000 (100%) | 0 |
| 11 | 0 | 5000 (100%) | 0 |
| 13 | 1 | 4999 | 0 |
| 15 | 1 | 4999 | 0 |
| 17–25 | 0 | 5000 (100%) | 0 |

Max σ observed (over all converged orbits): **381 steps**.

**Outputs:** [vis_a.csv](file:///C:/Collatz/visualization/vis_a.csv) (5.37M rows), [vis_b.csv](file:///C:/Collatz/visualization/vis_b.csv) (60K), [vis_b_div.csv](file:///C:/Collatz/visualization/vis_b_div.csv) (60K), [qx1_3d.blend](file:///C:/Collatz/visualization/qx1_3d.blend), eight PNGs in [snapshots/](file:///C:/Collatz/visualization/snapshots/).

---

## Q1. Does the trajectory cloud (A) show structure across q that's not visible at fixed q?

**Verdict: partially — the structure is q-stratified but the *across-q* signal is weak.**

The cloud forms 12 horizontal bands stacked along z (see [01_all_iso.png](file:///C:/Collatz/visualization/snapshots/01_all_iso.png), [02_all_front.png](file:///C:/Collatz/visualization/snapshots/02_all_front.png), [05_visA_only_iso.png](file:///C:/Collatz/visualization/snapshots/05_visA_only_iso.png)). One band per q. The geometric shape inside each band is what the existing experiments already tell us in 2D:

- **q = 3 band:** dense, compact, x ∈ ≈[1, 10] (log-values stay bounded — orbits converge); y dense up to ~200 (typical step counts). It is qualitatively the only stratum that "lives in a box."
- **q ≥ 5 bands:** all qualitatively similar to each other — log-value rises from log(n_start) to ≈42 (the int64 safe-cap log), then the orbit terminates. Each band is a wedge fanning right and slightly down as log-value grows. Bands for q = 7, 9, …, 25 are nearly translated copies of each other along the z-axis with slight shape differences in step count (high q reaches the cap in fewer steps).

What the 3D view *does* add over fixed-q slices: the visual confirmation that **q = 3 is geometrically isolated** from the q ≥ 5 stack — it doesn't smoothly connect to its neighbours. Everything ≥ 5 looks like a continuous family. That matches the known "q = 3 is the only converging value" result, but seeing it as a discontinuity in the stack is more vivid than reading it from a table.

What the 3D view does *not* add: any across-q geometric continuity beyond what the analytical Cramér-rate result already predicts. The bands for q = 7..25 differ only in slope and height; no qualitatively new feature emerges from the stacking.

## Q2. Does the stopping-time cloud (B) show structure across q?

**Verdict: structurally untestable for q ≥ 5 — there is no stopping time to plot.**

[06_visB_only_iso.png](file:///C:/Collatz/visualization/snapshots/06_visB_only_iso.png) shows the consequence: the only meaningful Vis B content is the 5,000-point cluster at z = q = 3 showing the well-known σ ≈ c·log(n) curve; every other q-stratum has at most a handful of points (the 17 + 1 + 1 lucky converged orbits) plus 55K sentinel points parked at y = -20.

The visualization confirms what the convergence table already says: at sample size 5K per q with cap = 10⁶ for n, the qx+1 system simply does not produce convergent orbits often enough at q ≥ 5 for a "stopping-time across q" hypothesis to be testable from this kind of brute-force sample. To get convergent orbits at high q you need either a much larger n-pool or to oversample residue classes that have small a★.

## Q3. When overlaid, do A and B show correlated geometric features?

**Verdict: only at q = 3, where both visualizations are populated.**

In [01_all_iso.png](file:///C:/Collatz/visualization/snapshots/01_all_iso.png) the q = 3 stratum shows both the dense Vis A trajectory band and the σ-vs-log(n) curve from Vis B occupying the same z-slice. They are visually disjoint within that slice (different y-ranges by construction — y is "step index" in A vs "σ" in B), and the y-axes aren't comparable, so this isn't an alignment finding so much as a co-residence finding.

For q ≥ 5, Vis B is empty or sentinel-only, so there is nothing to correlate. Vis B' (the divergence-aware variant) does give populated geometry at every q — see below.

## Q4. Is the a★ coloring spatially coherent (clusters/regions) or scattered (no structure)?

**Verdict: the cloud is a fairly uniform mix of all six a★ classes within each band — no obvious spatial clustering.**

The rainbow appearance in [02_all_front.png](file:///C:/Collatz/visualization/snapshots/02_all_front.png) and [08_low_angle_all.png](file:///C:/Collatz/visualization/snapshots/08_low_angle_all.png) initially looked like color stratification along x, but on inspection it is mostly the artifact of trajectories from different starting-n values overlapping in log-value space, not a true class-based clustering. Every trajectory is monochromatic in this coloring (a★ is a starting-class invariant), and adjacent trajectories with different a★ values occupy the same x-region.

Concretely: at a fixed (q, x), points from trajectories with a★ ∈ {1, q, q², q³, q⁴, q⁵, q⁶} all coexist with no apparent spatial separation. If there were class-based geometric structure — e.g., trajectories from a particular a★ class systematically occupying a particular log-value sub-band — it should show up as colored sheets within each q-band. It does not.

This is consistent with the analytical picture: a★ controls the deterministic *prefix* (the first few steps), not the long-range trajectory geometry, so once orbits leave their prefix region they should mix.

---

## Bonus: Vis B' (divergence-aware, log final value)

[07_visBdiv_only_iso.png](file:///C:/Collatz/visualization/snapshots/07_visBdiv_only_iso.png) shows what was added per the option-3 decision at the data-gen checkpoint. It is by far the cleanest geometric picture of the three views: 12 tight horizontal stripes at distinct y-values, one per q. The stripe height (y = log(final value at termination)) sits at log(2⁶²/q) for divergent orbits, which decreases by log(q) per stratum. The q = 3 stripe is at y ≈ 0 (orbits converge to value = 1).

This is essentially a triviality (all divergent orbits in the int64-bounded simulation hit the same overflow cap at the same log-value), but it does clearly visualize the q = 3 vs q ≥ 5 dichotomy as a discrete y-jump, which Vis B failed to deliver.

---

## One-line summary

**The 3D stacking shows that q = 3 is geometrically isolated from the q ≥ 5 family, but provides no across-q structure that the analytical Cramér-rate result and existing 2D plots don't already give. The a★ coloring shows no spatial clustering. A null result for the "geometry-emerges-from-stacking" hypothesis, with the q = 3 isolation as a small consolation prize.**
