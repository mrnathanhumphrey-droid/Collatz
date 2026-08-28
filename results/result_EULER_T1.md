# RESULT — EULER-T1: inverse-tree `λ=(3+√21)/6` ↔ Plancherel `S_∞`. NOT SUPPORTED — natural candidates refuted by the exact floor; opposite sides of the structure/value split; forward↔inverse already disposed. (2026-08-28)

**Probe:** `probes/probe_euler_t1.py`. The prize target of the Euler hunt (CORPUS_LEDGER §6, T1): does the
backward-tree growth constant `λ` link to the forward Plancherel limit `S_∞`?

## Express
`λ = (3+√21)/6 = 1.263762616`, degree-2 (`3λ²−3λ−1=0`, checked 0). Since `[ℚ(λ):ℚ]=2`, the sharp algebraic
form of any relation is **`S_∞ = a + b·λ`, `a,b ∈ ℚ`**. Structurally-natural sub-candidates:
`S_∞ = 2 log λ` (motivated by `log λ = 0.234093 ≈ lim T_i`) and `S_∞ = λ − 1/λ` (clean algebraic combo;
note `1/λ = 3λ−3` from the min-poly, so `λ − 1/λ = 3 − 2λ` — the same number).

Hard constraints (copied from ledger; re-grep before penning): exact floor `2·T_20 = 0.473177`
(`result_SOLSTICE`) with `S_∞` **strictly above**; `S_16 = 2·T_15 = 0.471352` (`result_P6I`); bracket
`[0.4714, 0.478]`; `7/15 = 0.466667` EXCLUDED; point estimate `≈ 0.475`.

## Disprove — Test A: the natural candidates die on the exact floor
- `S_∞ = 2 log λ = 0.468186948` → **BELOW floor** by `−0.004990`. REFUTED.
- `S_∞ = λ − 1/λ = 3 − 2λ = 0.472474768` → **BELOW floor** by `−0.000702`. REFUTED.

Both are killed not by a fit but by a **proven exact lower bound** (`2·T_20`, a specific rational). The
closest natural λ-object, `3−2λ = 0.472475`, sits tantalizingly `0.0007` under the floor and just above
`S_16` — but "just under the floor" is *excluded*, full stop.

## Disprove — Test B: numerics cannot confirm anyway, and λ systematically undershoots
Low-complexity `a+b·λ` near the top of the bracket cluster **below** the floor (`3−2λ = 0.472475`,
`−2/7+3/5·λ = 0.472543`). In the region where `S_∞` actually lives — `[0.473177, 0.478]` — there is **no
low-complexity λ-expression**; the only simple hit is `a=19/40` (`b=0`, i.e. λ absent). So the natural
λ-combinations *undershoot* `S_∞`, and where the value is, λ contributes nothing simple. With a bracket of
width `0.0066` and a target that has **no closed form**, a numerical near-match would carry no information
regardless — **any real T1 relation must be structural, not numerical** ([[facts_from_data]], mirror
guardrail).

## Disprove — structural: opposite sides of the one established split, and already-disposed
- **`λ`** is the degree-2 Perron root of the *finite* backward-tree offspring matrix — an
  associated-graded / finite-place object: the **symmetric, law-governed side** of the structure/value split.
- **`S_∞`** is the forward Plancherel *limit* — infinite Mahler depth, doubly-exponential denominators,
  value in the tail, leading rational `7/15` overshot: the **pro-limit / archimedean value side** no
  framework delivers.

A clean `S_∞ = a+bλ` would drop the walled tail-value into a degree-2 number field — collapsing the
infinite Mahler depth into an algebraic-degree-2 value, i.e. **defeating the corpus's one established
split**. No mechanism does this. Moreover `duality_S_vs_D_verdict.md` already tested the forward↔inverse
bridge directly and **disposed it negative** (`D·S, D+S, D/S`, diagonal all fail; the one near-relation
`D_n(1)/S_1 = 1/3` is "a trivial coincidence"). T1 asks the growth constant of the inverse object to relate
to the limit of the forward object — across a split never bridged and a duality already killed.

## Verdict
**T1 NOT SUPPORTED (disposed-negative at the level testable).** We cannot *prove nonexistence* — `S_∞` has
no closed form to test an arbitrary `f(λ)` against — but: (1) both structurally-natural candidates are
**refuted by the exact floor**; (2) the value region admits **no low-complexity λ-expression** and numerics
cannot discriminate there anyway; (3) a relation would have to **defeat the structure/value split**, for
which there is no mechanism, and the sibling forward↔inverse duality is **already disposed**. The prize is
not here. `λ` lives on the finite/symmetric side; `S_∞`'s value lives in the tail — consistent with every
other framework in the corpus: **the value is delivered by nothing on the structured side.**

**Not at stake:** `λ`'s own algebra, `S_∞`'s floor/bracket, the MAHLER spine. T1 tested a cross-side link
and found none.
