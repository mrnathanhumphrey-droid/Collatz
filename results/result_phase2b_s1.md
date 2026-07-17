# Result — PHASE 2b Session One (agent requests A–E, C). Instrument work only; no proof committed here.

**Date:** 2026-07-16. Deliverables for Claude's D1 hand-derivation. **Headline: the toy's subdominant is EXACTLY λ₂ = −(1−λ)/(1+λ), living on a 2-state DIAGONAL ray; the earlier "r=0.831" was a mass-sequence artifact (the mass is pure Perron) — the true r(λ)=|λ₂|/λ₁ is q-flat, L-stable, and equals 0.60 at λ=1/2. L=4 (Request C) WALLS OUT for accurate gap resolution.**

## ⚠️ CORRECTION to Q6's "0.831 flat"
The toy's **total-mass** sequence `a_k = 1ᵀMᵏv₀` is **pure Perron**: every subdominant eigenvector is mean-zero (`1ᵀrᵢ = 0`, amplitude Aᵢ ≈ 0), so it is INVISIBLE to the total mass. ESPRIT/mass-detectors therefore fit noise — the "0.831 flat" (Q6) was that artifact. **The true subdominant is read from the RAW operator spectrum**, where λ₂ is a genuine (mean-zero) mode. This does not change Q6's qualitative conclusion (toy gapped, q-flat) — only the value (0.60, not 0.831) and the reason.

## Request A — constructor freeze (the object Claude derives against)
`M(q, g, λ)` on states `(a, b, γ)`, `a,b ∈ ⟨g⟩ mod q^L`, `γ ∈ ℤ/q^L`:
- **Phase group** `⟨g⟩ mod q^L`, order `ordn`. For `g=−1`: `{1, q^L−1}`, order **2 at every L** (never lifts).
- **Coordinates & weights (THE FROZEN FREE CHOICE):** v is **folded into the subgroup period** — one coordinate `δ = 1..ordn`, phase `g^{−δ}`, raw weight `λ^δ`, normalized `w = raw/Σraw`. For `g=−1`: `w = [1/(1+λ), λ/(1+λ)] = [u, s]` (flip u, stay s). **This fold is the "silent choice" named:** the phase coordinate is the *parity class*, not the raw v. Consequence: the Perron rate is `λ₁ = Σw² = u²+s² = (1+λ²)/(1+λ)²` (measured **5/9** at λ=1/2), which is DIFFERENT from the unfolded-ansatz `q₊ = (1−λ)/(1+λ) = 1/3`. **Derive against the folded λ₁.**
- **Gate / carry (frozen-phase handling — your explicit question):** transition `(a,b,γ) → (a·g^{−δ_a}, b·g^{−δ_b}, γ')` with weight `w_{δ_a}w_{δ_b}`; `T = a'−b' mod q^L`; **survive iff `(γ+T) ≡ 0 mod q`**, then `γ' = ((γ+T)//q) mod q^L`; else the branch dies (sub-stochastic). For frozen phases (⟨−1⟩, no tower) at levels ≥2, the digit condition still reads the SAME `(γ+T)≡0 mod q` — the phases just don't refine, so `a',b' ∈ {±1}` and `T ∈ {0, ±2a'}` only. **No special-casing: the gate is level-agnostic; the phase group's failure to lift is what keeps `T` in `{0,±2}`.**

## Request D — exact λ₂ at L=1, λ=1/2 (pins the polynomial)
| q | λ₁ | λ₂ (exact) | r=\|λ₂\|/λ₁ |
|---|---|---|---|
| 5 | 5/9 | **−1/3** | **3/5** |
| 7 | 5/9 | **−1/3** | **3/5** |

λ₂ = **−1/3 exactly** — a *rational* (degree 1 over ℚ), q-flat. Not 0.831. (Distinct \|eig\|: {5/9, 1/3, 0}.)

## Request E — λ₂ eigenvector localization (L=2, q=7, λ=1/2)
λ₂ = −1/3; its eigenvector is supported **entirely on the 2-state diagonal ray** `{(1,1,0), (−1,−1,0)}` (e=+1, γ=0), as the **antisymmetric** combination `(+1,−1)` (Perron is the symmetric `(+1,+1)`). The full diagonal-ray 2×2 block is
`[[s², u²], [u², s²]]` (stay-stay keeps state, flip-flip swaps to the mirror), eigenvalues `s²±u² = {λ₁, λ₂}`. **So the subdominant is the antisymmetric diagonal mode; λ₂ = s²−u² = −(1−λ)/(1+λ).** The skeleton search is confirmed to be this 2-dim diagonal subspace (no carry mixing needed at this order).

## Request B — λ-sweep r(λ)=\|λ₂\|/λ₁ (RAW spectrum), gate for the derived r(λ)
| λ | 0.30 | 0.40 | 0.50 | 0.60 | 0.70 |
|---|---|---|---|---|---|
| **q=7 (L=2)** | 0.8349 | 0.7241 | 0.6000 | 0.4706 | 0.3423 |
| **q=13 (L=2)** | 0.8349 | 0.7241 | 0.6000 | 0.4706 | 0.3423 |
| **q=7 (L=3)** | 0.8349 | 0.7241 | 0.6000 | 0.4706 | 0.3423 |

**q-FLAT (q=7≡q=13 to all digits) and L-STABLE (L=2≡L=3)** — the pre-registered q-independence CONFIRMED, no deviation. The curve is smooth and monotone-decreasing; the λ=1/2 value is exactly 3/5. (These are the values your derived r(λ) must hit.)

## Request C — the L=4 shot: WALLS OUT (reported, not extrapolated)
Reduced reachable operator built (n=54 / 1458 / **39366** at L=2/3/4; L=2,3 reproduce the standing gaps 2.9e-3, 1.0e-4 exactly — tool validated). **But the L=4 coalescing gap cannot be reliably resolved:** LM-ARPACK returns 1.67e-4, at/below its accuracy floor for a sub-1e-4 cluster in a 39366-dim operator (so **not a measurement** — it even reads non-monotone, a known ARPACK-cluster failure); **shift-invert at σ=1/3 (the accurate tool) times out** (LU factorization of the 39366-dim operator too heavy locally). **Per the guard: L=4 walls out — I report the wall and stop, no extrapolation substitute.**
- Standing reliable data: gap(L=1,2,3) = 0.889, 2.9e-3, 1.0e-4.
- 27^{−L} pattern (3.8e-6) and naive tower 2^{−2·3^{L−1}} (5.5e-17): **both remain unadjudicated at L=4** (a clean L=4 needs a direct solver — Lambda, or a smarter carry-skeleton reduction than the full reachable set).

## Not at stake
R1–R46, Phases 0/1/2a. This is Session-One instrument work: it freezes the constructor (A), hands D1 the exact λ₂=−1/3 (D) localized on the diagonal ray (E), the q-flat r(λ) sweep (B), and reports the L=4 wall (C). It also corrects Q6's 0.831→0.60 (mass-artifact).

_Reporting discipline: the 0.831→0.60 correction is disclosed as a correction (the mass sequence is pure Perron — subdominants mean-zero — so ESPRIT fit noise). λ₂=−1/3 and the diagonal 2×2 block are exact (read from the operator, not fit). The sweep q/L-flatness is verified across q=7,13 and L=2,3. C is reported as a WALL (ARPACK floor + shift-invert timeout), not filled with the 27^{−L} extrapolation._
