# BRIEF — THEOREM D1: the toy spectral gap `r(λ) = (1−λ²)/(1+λ²)`

**The first hand-derived spectral gap of the L3 program.** Committed 2026-07-16 by the derivation, then met the pre-published sweep table (five-for-five). Author: Nathan (derivation); Claude (probes/gates). This is the toy `M(q, −1, λ)` — the frozen constructor of `result_phase2b_s1.md` Request A. It is the *template* for the real theorem, not the real theorem.

Notation: `u = 1/(1+λ)` (flip weight), `s = λ/(1+λ)` (stay weight), `u + s = 1`.

---

## THEOREM D1. For `M(q, −1, λ)` as frozen — any odd prime `q`, any level `L`:

### (a) Exact invariance of the diagonal ray.
The 2-state ray `{(1,1,0), (−1,−1,0)}` is **exactly closed**:
- stay-stay (weight `s²`) holds the state — `T=0`, gate passes at `γ=0`, carry stays `0`;
- flip-flip (weight `u²`) swaps to the mirror — also `T=0`, gate passes, carry stays `0`;
- the two exit branches flip-stay / stay-flip (weight `us` each) produce `T = ±2` against `γ=0`, and **`±2 ≢ 0 mod q` for every odd `q`** — they die at the gate.

The block on the ray is exactly
```
[[ s²,  u² ],
 [ u²,  s² ]]
```
sub-stochastic with defect `2us` (the killed exits).

### (b) Spectrum (read off the 2×2 block).
- **Symmetric mode (Perron):** `λ₁ = s² + u² = (1+λ²)/(1+λ)²`  →  `5/9` at `λ=½` ✓
- **Antisymmetric mode (subdominant):** `λ₂ = s² − u² = (s−u)(s+u) = s − u = −(1−λ)/(1+λ)`

Because `s + u = 1`, the subdominant is **literally the first signed moment `m₋`** of the pair operator — degree one in the moments. That is why it is so clean. At `λ=½`: `λ₂ = −1/3` — matching Request D's exact algebraic answer at both q=5 and q=7.

### (c) The gap, closed form.
```
r(λ) = |λ₂| / λ₁ = (1−λ)(1+λ)/(1+λ²) = (1 − λ²)/(1 + λ²)
```

Against all five sweep points, **exact rational arithmetic** (curve published *before* the derivation touched it, per the gate design):

| λ | r(λ) exact | decimal | measured (q=7,13; L=2,3) |
|---|---|---|---|
| 0.3 | 91/109 | 0.834862 | 0.8349 ✓ |
| 0.4 | 21/29 | 0.724138 | 0.7241 ✓ |
| 0.5 | 3/5 | 0.600000 | 0.6000 ✓ |
| 0.6 | 8/17 | 0.470588 | 0.4706 ✓ |
| 0.7 | 51/149 | 0.342282 | 0.3423 ✓ |

**Five for five, both primes, both levels.** Not fit to the curve — the curve was published first.

### (d) Structure notes.
- **`r < 1` for all `λ ∈ (0,1]`** — the toy is gapped at every weight. Its only boundary is the degenerate `λ → 0` (pure deterministic flipping: perfect sign memory, `r → 1`). That is the toy's own "too rhythmic to mix," and it sits at a **different location** than the real map's `λ=½` resonance — because **the toy has no resonance: it has no tower.**
- **`λ₂ < 0` ⇒ oscillatory approach** — the sign-alternating genre R26/R27 measured on the real object, now derived in miniature.
- (zero evidential weight, filed for enjoyment) `λ = tan(θ/2)` gives `r = cos θ`; at the Syracuse weight `λ=½` the pair `(r, 2λ/(1+λ²)) = (3/5, 4/5)` — **the toy sits on the 3-4-5 triangle.**

---

## Remaining gap in D1 (stated honestly): MAXIMALITY.
Proving **no `e = −1` (off-diagonal) mode beats `|λ₂|`.** The mechanism is in hand — it is the carry skeleton's surviving job. The dangerous candidate is the **flip-flip pseudo-cycle at weight `u² = 4/9 > 1/3`**; the **carry kills it**: from `(a, −a, q−2a)`, flip-flip passes, but lands carry at `γ=0` in the mirror, where both `e=−1` continuations need `γ ± 2a ≡ 0` and fail — the mass is forced out through `us`-exits within two steps. Every `e=−1` walk dies or exits on that schedule.

Turning "every `e=−1` walk exits in ≤2 steps" into "spectral radius of the `e=−1` block `< 1/3`, all `λ`" is a **finite reachable-graph argument** — Nathan's pen, next session. One cheap assist queued first:

- **Request F** (pre-registered): numeric spectral radius of the `e=−1` sub-block **alone**, `L=1,2`, `λ = 0.3 / 0.5 / 0.7`. Pre-registration: **strictly below `|λ₂|` at every point** (expect `us`-to-`s²` scale). One number per case — makes the graph proof a **confirmation, not an exploration.**

---

## Why this brief matters to the real theorem.
The technique that cracked the toy — **find the exactly-invariant ray, diagonalize its block, then fight for maximality** — is now the **template walking into the real q=3 operator.** The D3 lead (`result_phase2b_s2.md`) shows the real operator has such a ray: one member of its coalescing pair is already closed-formed as the folded diagonal moment `Σ w_r²`. If a real invariant ray yielding `Σ w_r²` exists (and the six-digit L=3 match says it does), the real operator inherits this proof architecture: **one exactly-solvable ray + one dynamical partner + a gap statement between them.**
