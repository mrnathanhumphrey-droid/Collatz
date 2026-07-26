# RESULT — THE RE-CENTER: three functionals, Tao owns ℓ^∞, S_∞/channels are ours; both flags pulled from primary source (2026-07-26)

**Probes/reads:** `probe_calib.py` (T1); direct reads of `scratchpad/tao.txt` (Tao 1909.03562) and
`scratchpad/siegel.txt` (Maxwell Siegel, *(p,q)-adic Analysis and the Collatz Conjecture*, arXiv:2412.02902 =
the consolidated dissertation). Wilson's course-correction: the sup|π̂| chain was Tao's already-solved functional;
re-center onto the aggregate/channels. He asked both flags pulled with exact statements before banking.

## T1 — calibration (sup|π̂| decay vs Tao n^{−A})
On the certified LAMBDA sup values (k=3..15), fitting `sup(n) ~ C n^{−A}`:
- **Per-step `A_eff = ln(Srate)/ln(n/(n+1))` climbs 1.23 → 2.33** (trend +0.089/k); windowed power-law fit climbs
  1.38 → 2.15 monotone (all R²>0.996). **The climb IS the superpolynomial signature** (a fixed power gives constant A).
- **Pure-exp REJECTED:** Srate would be constant; it climbs 0.70 → 0.85. Stretched-exp best α≈0.43 (I get 0.43, not
  Wilson's 0.14 estimate — but α<1 corroborates near-power-law, not exp). This **contradicts exp(−cm) at these depths**
  (Tao Remark 1.15's heuristic prediction, which Tao explicitly does *not* prove — tao.txt line 979/1017).
- **Kills the tower-vs-factorial worry:** the useful regime is `A ≈ 2` (grind A=2,3,4 and STOP). Even a tower-in-A
  constant is tower(2.5) = nothing; the blow-up only bites at large A, which we never need.

## The three-functional sharpening (Wilson) — Parseval makes it a number
`‖π_k‖² = ⅓‖π_{k−1}‖² + 3^{−k} Σ_{3∤ξ}|π̂|²`, and Tao's `|π̂(ξ)| ≤ C_A k^{−A}` gives `‖π_k‖² ≲ k^{−2A}` — but the truth
is `‖π_k‖² = X_k/3^k ~ c·k·3^{−k}`. **Superpolynomial vs geometric: the sup does NOT give the aggregate.**
- **ℓ^∞ — sup|π̂|.** Tao Prop 1.17. = R66, CONTRACTION, GRECURSION, LAMBDA. **Closed, his turf.**
- **ℓ² — aggregate ‖π_k‖² = X_k/3^k.** Siegel's Parseval object, our **S_∞** and `A_r(0)`. Geometric, untouched by
  the sup bound.
- **ℓ⁴ / full transform — the channels.** `γ_r(k)=Σ_a|π̂(a)|² e(ak/3^r)` = Fourier transform of the power spectrum;
  `mean_k|γ_r(k)−1|² = Σ_{a≠0}|π̂|⁴ = U²(ρ_k) = 0.29754` level-invariant. **The channels need the spectral
  DISTRIBUTION — not its max, not its total. Nobody has an estimate for it.**
- **How vanishing the old handle was:** the sup mode's contribution to any channel is ~2^{−2m}, m≈4k/3, i.e.
  **~3^{−1.68k}** — against a total `Σ_a|π̂|² = X_k ~ k`. The sup→bridge→channels chain was controlling an O(k)
  quantity through a piece of size 3^{−1.68k}: not a suboptimal handle, a **vanishing** one. That is exactly why
  CONTRACTION/GRECURSION/LAMBDA kept resolving into "known/superpolynomial/δ→0 no floor" — they were Prop 1.17 all
  the way down.

## FLAG 1 (Tao Prop 1.14 — could TV reach the aggregate?) — NO, resolved by the source
tao.txt Prop 1.14 (line 946): `Osc_{m,n}(P(Syrac=Y)) ≪_A m^{−A}`, and (line 977) it literally equals
`dTV(Syrac, Syrac + Unif(3^mℤ/3^nℤ))` — so Wilson's read of its *type* (total-variation closeness to the geometric
model) is right. **But Remark 1.18 (line 995–997): "it is not difficult to use the triangle inequality to establish
|E e(…)| ≤ Osc_{n−1,n}(…) … Thus Proposition 1.17 and Proposition 1.14 are in fact equivalent."** So the TV/fine-scale
functional is the SAME STRENGTH as the sup bound — it does **not** reach the ℓ² aggregate. Door closed by Tao himself.

## FLAG 2 (Siegel p.92–93 — sharp constant or qualitative gap?) — NEITHER; the brief mis-read it
Two candidate locations, both refute the paraphrase "Siegel flags ‖π_k‖² decay as the open problem":
- **PDF 92–93 = diss 81–82 = "Connections to Tao (2019)"** (siegel.txt line 6459+). Siegel **cedes the decay to Tao**:
  "the central challenge Tao overcomes is obtaining explicit estimates for the decay of the characteristic function of
  the Syracuse Random Variables. In our terminology, Tao establishes decay estimates for the archimedean absolute value
  of the function **φ₃** … Tao's decay estimate is given in **Proposition 1.17**" (lines 6472–6490). His φ₃ (eq. 2.171)
  IS our π̂ / Tao's Syracuse char function. His stated open direction: "a natural next step … would be to study the
  characteristic function of χ_H for an **arbitrary semi-basic p-Hydra map**" (line 6491) — GENERALIZATION, not the
  sharp Collatz constant. And **Prop 2.18 / eq. 2.173–2.174** `φ_H(dt)=(1/p)Σ_j e^{−2πiβ_j t}φ_H(α_j t)` (line 6513–60)
  = **our single-level recursion (SINGLEREC), verbatim.**
- **Literal diss p.92–93 (PDF 103–104)** = Chapter 3 methods: **Dwork's Theorem** (rational power series via
  `Π_{p∈S} R_p > 1`) and Borel's observation — i.e. **the Mahler/house/Dwork RATIONALITY shelf**, exactly what the
  free-check (Π_{3∤a}π̂ rational?) feeds. Not an open-problem statement at all.

So Siegel flags **neither** the sharp constant nor a qualitative gap as open; he attributes the ℓ^∞ decay to Tao and
works a different axis (generalization + the (p,q)-adic Wiener–Tauberian program). Siegel's own words seal the object
identity: he renamed χ_H "the characteristic function of H" to "**numen**" precisely "in light of Tao's work — where,
in essence, one works with the characteristic function (in the probability-theoretic sense) of χ_H" (line 831–834), and
notes he "**has yet to grok Tao's probabilistic arguments**" (line 512) — so the two frameworks were never bridged, by
his own admission. ⚠️ **Correction to PHASE3 brief:** "Siegel flags ‖π_k‖² decay as the open problem (p.92–93)" is not
supported at either candidate page; Siegel cedes that decay to Tao.

## The re-centered map — one object, three functionals, three sources
| functional | object | who owns it | status |
|-----------|--------|-------------|--------|
| **ℓ^∞** sup\|π̂\| | R66 / CONTRACTION / GRECURSION / LAMBDA | **Tao Prop 1.17 ⟺ 1.14** | closed (superpolynomial, C_A extractable §7) |
| **ℓ²** ‖π_k‖²=X_k/3^k | S_∞, A_r(0), the domination/gap | Siegel's Parseval; **nobody sharp** | **OPEN — ours** |
| **ℓ⁴/full** γ_r(k) channels | CHANNEL_ID, v₃ HIERARCHY, U²(ρ) | **nobody** | **OPEN — ours, richest** |

**Two things survive the re-center intact:**
1. **The v₃ HIERARCHY is a different genre** — `γ_j(k)=X_j` for `j≤v₃(k)`, class means 5/3 & 2/3, from `ν_r mod 3=ν_1`
   + the tower. **Exact structural identities, not estimates.** Nothing in Tao's or Siegel's analytic machinery bears on
   them.
2. **The Galois observation is aggregate-side:** `Π_{3∤a}|π̂(a)| = |N(π̂(1))|` is the **geometric mean** of the
   spectrum (not the sup) — a spectral-distribution functional, the class the channels live in. On the right side of
   the line, and it plugs straight into the Dwork/Borel rationality shelf at diss p.92–93. **Keep it** (free-check
   pending: verify rational at k=2..6, report value + v₃).

**Net:** the sup was Tao's; the prize (S_∞ / the channels) is the aggregate + spectral-distribution, which neither Tao
(who proved only the ℓ^∞ decay) nor Siegel (who cedes the decay to Tao and generalizes elsewhere) computes. Re-center
the hunt on the channels. **Not at stake:** CHANNEL_ID, v₃ HIERARCHY, MEAN1, U², R1–R30. (Effective-C_A §7 grind =
Hank, in flight; free-check = pending.)
