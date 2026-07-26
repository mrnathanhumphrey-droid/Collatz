# LIT HUNT (Hank) — the object IS Tao's Syracuse char. function; superpolynomial-uniform is a THEOREM, C_A extractable not barred (2026-07-26)

**Dispatched:** effective/uniform power-saving for `Σ_{b<L} e(a·2^b/3^n)`. Redirected mid-run (Wilson): uniform power
saving is probably false → the real target is effectivity of Tao's superpolynomial constant `C_A`, and whether
entropy-decrement arguments admit effectivity at all. Hank read Tao §7 primary source directly and **corrected the
proposed reason**.

## Headline — a terminus WITH a reason, and better than "barred"
Our object is exactly the **Syracuse characteristic function** of Tao 1909.03562, and the decay is **superpolynomial,
not power-saving** — which is precisely why our measured δ shrinks with no floor. But Wilson's proposed *reason* (an
entropy-decrement / unlocated-scale ineffectivity barrier) is **NOT what Tao's proof uses.** The superpolynomial decay
is proved by a **Fourier/Plancherel reduction + a 2-D renewal process + elementary triangle geometry + downward
induction on scale** — fully quantitative, constructive, **no entropy decrement, no compactness, no
pigeonhole-over-unlocated-scales.** So the honest terminus is: **"uniform superpolynomial is a THEOREM; an explicit
C_A is not written down but is not blocked — it's a bookkeeping extraction from an effective-in-principle argument."**

## 1. The exact home
**Tao, "Almost all orbits of the Collatz map attain almost bounded values,"** arXiv:1909.03562 (Forum Math Pi 10
(2022) e12). Syracuse random variable on ℤ/3ⁿ (eq. 1.29): `Syrac = 2^{−a₁} + 3·2^{−a_{1,2}} + ⋯ + 3^{n−1}·2^{−a_{1,n}}`,
`(a₁..aₙ)~Geom(2)ⁿ` — **the Geom(2)-weighted (probabilistic) form of our Σ e(a·2^b/3ⁿ)**; the top orbit terms
3^{j−1}2^{−a} are our lattice, and `E[e(−ξ·/3ⁿ)]` is the weighted partial exponential sum whose flat/unweighted
ceiling is our tree bound ≤ 2^L.

**Prop 1.17:** for n≥1, 3∤ξ, every A: `|E e(−2πi ξ·Syrac/3ⁿ)| ≪_A n^{−A}`, with the implied constant **uniform in n
AND in ξ** (permitted to depend on A). Equivalent to fine-scale mixing Prop 1.14, `Osc_{m,n} ≪_A m^{−A}`.

## 2. Uniformity/shape verdict — matches our numerics precisely
- **Superpolynomial, not power-saving.** n^{−A} for every A is faster than any fixed L^{−δ}, but the constant grows
  with A. Converting to a per-level saving exponent, `δ_eff ~ (A log n − log C_A)/levels` → **appears to shrink toward 0
  as n grows at any fixed measurement budget, with no positive floor.** Our LAMBDA (A)-vs-(B) dichotomy resolves to
  **(B)-shaped**: no uniform fixed-δ power saving, and none expected. **δ ≈ 0.30 → 0.12 monotone-no-floor is the correct
  fingerprint of n^{−A} with A-dependent constant.** LAMBDA independently reproduced Tao's shape from the spectral side.
- **Uniformity IS achieved — in the surviving sense:** the constant is uniform in n and ξ. The right positive statement
  is "**uniform superpolynomial**," not "uniform power-saving." (Tao Remark 1.15: the *heuristic* predicts even stronger
  `exp(−cm)` fine-scale decay, "which we do not attempt to establish" — still no power-saving-in-n.)

## 3. Is C_A effective? (priority-1)
**As written: no explicit C_A anywhere** (≪_A throughout; Prop 7.8 "some sufficiently large C_{A,·}"). No follow-up makes
it explicit. **In principle: extractable, NO ineffective step.** Hank traced §7: 7.1–7.2 char. function = average of
products of cosines over a 2-D renewal walk, black set partitioned into disjoint separated triangles (Lemma 7.4,
constructive); 7.3–7.4 Prop 7.3 reduces to `E exp(−3·#{white points hit}) ≪_A n^{−A}`, a **non-negative large-deviation
estimate** ("does not require capturing cancellation"), closed by **downward induction on m≈n/2−j** (base Q_m≤m^A + a
monotonicity step, Prop 7.8). Every constant (renewal holding-time tails 7.6, first-passage bounds 7.7, the induction)
is quantitative. **Verdict: an explicit C_A is a tedious extraction from §7 — likely tower-/factorial-in-A but fully
computable — NOT a barrier.**

## 4. The entropy-decrement question (priority-2) — real, citable, but WRONG PAPER
The structural ineffectivity Wilson described is genuine and attributable — **Tao, "The logarithmically averaged Chowla
and Elliott conjectures,"** arXiv:1509.05422 (Forum Math Pi 4 (2016) e8), the entropy-decrement method: a good scale H
"cannot be specified in advance" ⇒ **log-averaging cannot be removed** (harmless for log-averages, fatal for natural
ones); same mechanism as Tao's Erdős-discrepancy resolution. **But it does NOT apply here.** Tao's Collatz Prop 1.17
does not use entropy decrement — "entropy" appears in 1909.03562 only in the *heuristic* Remark 1.15 (Shannon entropy of
Geom(2)=log4). **⚠️ Do NOT attribute our δ→0 to the entropy-decrement barrier;** the reason is simply that the true decay
is superpolynomial, established effectively-in-principle by the renewal argument.

## 5. Classical shelves — all VACUOUS at L~n=O(log q) (corroborating, demoted to context)
- **Korobov** (`Σ_N e(a bⁿ/m) ≤ √m(1+log m)`): nontrivial only for N ≳ √m = 3^{n/2}; our L~(4/3)n = O(log q) is
  *exponentially* below threshold ⇒ vacuous.
- **Differencing** (Banks–Shparlinski, arXiv:1606.07911, covers prime-power m): threshold only down to
  `N ≥ exp(log m/log₃loglog m) = m^{o(1)} ≫ log m` ⇒ still vacuous; saving degrades with m.
- **BGK / Konyagin lecture notes** (local `C:\Collatz\Bourgain-Konyagin\Konyagin_Lectures.pdf`, Hank read it): prime
  *field* full-subgroup sums needing |G|≥p^δ — two disqualifiers (prime field not prime power; full subgroup not a
  length-log-q initial segment). Off-object.
- **Erdős base-3 anchor** (Narkiewicz N(X)≤1.62X^{0.63}; verified n≤2·3⁴⁵): sharp control of 2^b mod 3ⁿ at L~n = the
  ternary-digits-of-2ⁿ problem, of which Erdős said "there is no method at our disposal to attack this" — the
  independent frontier signal that per-orbit-segment control at L~n is at/beyond the current frontier.

## 6. Bottom line
1. **Best-known bound = Tao 1909.03562 Prop 1.17**, `|E e(−ξ·Syrac/3ⁿ)| ≪_A n^{−A}`, uniform in n and ξ (3∤ξ). THE
   published home of our sum.
2. **Uniform power saving (fixed δ): FALSE / not expected.** δ→0 monotone-no-floor = correct fingerprint. Stop hunting δ_∞.
3. **Effective superpolynomial: KNOWN qualitatively (uniform), NOT written explicitly, NOT barred.** Constructive proof;
   C_A extractable from §7 (renewal + triangles + induction). No entropy decrement, no unlocated-scale obstruction.
4. **Framing correction:** entropy-decrement ineffectivity is real & citable (1509.05422) but belongs to
   Chowla/Erdős-discrepancy, NOT the Collatz estimate.
5. **Best route to an effective rate:** grind constants through §7 (Prop 7.1→7.3→7.8). For the *stronger* exp(−cm)
   fine-scale decay: expected but unproven (Remark 1.15) = a genuinely open improvement.

**Artifacts:** paper text `scratchpad/tao.txt` (Prop 1.17 ≈ line 984; eq. 1.29 ≈ 1005; §7 proof ≈ 3015–4750). Local
Konyagin PDF = off-object. Not at stake: LAMBDA/GRECURSION/CONTRACTION/SINGLEREC/BRIDGE2/CHANNEL_ID/R1–R30.
