# Result 84 — The family-distinguishing mod-9 offset

**Date:** 2026-07-14. **Verdict: H_ARTIFACT — the offset tracks the modulus (3^r), it is not a fixed 3². The Probe-83 lead is a normalization artifact and is removed.**

Probe `probe_84_mod9_offset.py`; data `result_84_data.csv`; log `result_84_log.txt`.

## Task A — c_{ℓ,ε} definition: no 3-power factor, but a top-layer twist phase

`c_{ℓ,ε} = 2^ε·(1+3^r)^ℓ mod 3^{r+1}` (result_78_FINAL 78.1 / result_81b). Every family is a **unit** (`v₃(c)=0`, verified r=2,3,4): **no** `3^{k(ℓ,ε)}` factor and **no** coset-representative choice differing by a multiple of 9. So the offset is *not* a crude bookkeeping ghost of a 3-power in c. **But** the family-defining twist `(1+3^r)^ℓ = 1 + ℓ·3^r + O(3^{2r})` carries a global phase: since **4 ≡ 1 mod 3**, `e_q(ℓ·3^r·4^u) = e_3(ℓ·4^u) = e_3(ℓ) = ω₃^ℓ` is *constant in u*, so it multiplies the whole chirp — and thus ĝ — by a global cube-root phase living at the **top 3-adic layer, order 3^r**. That predicts an offset scaling as 3^r, which Task B confirms.

## Task B — the decisive test: v₃(offset) vs r

Predicting level r from the fixed r=9 profile, the offset `pred − cert mod 3^{r+1}` is constant in b (a global phase) for every family, and:

| r | modulus 3^{r+1} | offsets (ℓ=0,1,2 / ε=0) | offsets (ε=1) | v₃(nonzero) |
|---|---|---|---|---|
| 2 | 27 | [18, 9, 0] | [9, 18, 0] | [2] |
| 3 | 81 | [0, 54, 27] | [0, 27, 54] | [3] |
| 4 | 243 | [0, 162, 81] | [0, 81, 162] | [4] |
| 5 | 729 | [0, 486, 243] | [0, 243, 486] | [5] |
| 6 | 2187 | [0, 1458, 729] | [0, 729, 1458] | [6] |
| 7 | 6561 | [0, 4374, 2187] | [0, 2187, 4374] | [7] |

**v₃(offset) = r at every level r=2..7.** The nonzero offsets are exactly `3^r·{1,2}` — the top modulus layer — not a fixed `3²`. The `9` Probe 83 reported was simply `3^{r=2}`. (The ℓ=0 offset is nonzero only at r=2 and vanishes for r≥3, so even the r=2 value was partly an extrapolation-boundary effect — the coincidence is complete.)

**Why the Probe-83 character check missed it:** that check ran on `offset/9` at the *single* level r=2 and asked whether `(ℓ,ε)→offset/9∈Z/3` is a character. It is not — but the reason is not hidden family structure; it is that the object isn't a fixed mod-9 quantity at all. Its 3-adic layer moves with r (`3^r`), and the ℓ,ε-pattern rotates with r (r=2: {2,1,0}; r≥3: {0,2,1}). A one-level check cannot see that. **Exactly the failure mode the pre-reg §2 flagged: a normalization ghost the character check was too coarse to catch.**

## Task C — not earned

H_CONST did not fire (the offset is r-dependent, tracking 3^r), so there is no r-independent family invariant `f(ℓ,ε)=±(2−ℓ)` to derive. The `±(2−ℓ)` pattern was the top-layer coefficient of the twist phase read at one level; it is not a spine member.

## Task D — the ε pointer: closed

From the definition, not the pattern: in `c_{ℓ,ε}=2^ε·(1+3^r)^ℓ`, ε toggles the `2^ε` factor — a **2-adic doubling**. The sibling `3x±1` map is `σ(r)=−r mod 3^k` (`K₋=σK₊σ`), a **negation**, unrelated to ×2. So ε is **not** the sibling map; the `f(ℓ,1)=−f(ℓ,0)` antisymmetry is two involutions coinciding (weak evidence). **Pointer closed** — and moot, since the offset is a normalization artifact.

## Consequences (per §6)

- **The Probe-83 lead is removed.** R83's disposition and STATE entry called the offset "genuine residual structure, not a normalization artifact." That is **overturned**: it is a top-layer (`3^r`) global phase from the `(1+3^r)^ℓ` twist — a property of the c_{ℓ,ε} normalization. The R83 disposition is annotated with a correction banner; the R83 STATE entry is corrected in place.
- **R81b certification amendment — FLAGGED, not written (for the operator):** the phrase "one r-independent 3-adic analytic function" should read "**r-independent in shape; the global phase (constant term) tracks the top modulus layer 3^r as the family-defining twist `(1+3^r)^ℓ` dictates**." The shape coefficients c₁,c₂,… do transfer across levels exactly (Probe 83 Task B); only the global phase is level-coupled, and now explained. This is a one-line refinement of R81b, surfaced for the operator to apply.

- **Untouched:** THEOREM_C_745 (7/45), Thms 78.1–78.3, and the R81b *shape* certification all stand. Nothing here bears on them.

_Reporting discipline: a cheap, decisive null — the offset dissolves as a modulus-tracking artifact (§5 outcome #2). The lead is removed where it lives. The v₃-vs-modulus check was run first per §4; no fitting was performed; the answer to Task A and Task D came from the definition, not the pattern._
