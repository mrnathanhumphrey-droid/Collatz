# The Euler move — hunt the corpus for INTERNAL RELATIONSHIPS (2026-08-28, Wilson's side-quest)

## The idea (in Claude's framing)

Euler's leverage on `ζ(s) = Σ 1/n^s` did **not** come from grinding the asymptotics of the sum. It came from a
**second representation of the same object** — `ζ(s) = Π_p 1/(1−p^{−s})` — and from *equating two faces of one
thing* (the sum-face and the product-face), which is where the Fundamental Theorem of Arithmetic falls out. The
Basel value `ζ(2) = π²/6` likewise came from a **relation** (the series ↔ the `sin x / x` product), not from
summing harder.

**Wilson's observation:** across the entire Collatz corpus we have hammered the **asymptotic wall** — is `S_∞`
rational, what is its arithmetic type, how fast does the tail decay — and we have essentially **never** stepped
back to look at the **web of internal relationships** among the many constants, series, kernels, and structures we
have banked. We have been computing one face (the value) and never asking whether two of our *own* objects are
two faces of one thing.

**The target, stated as a slogan:** *find our `π²/6`* — a surprising internal identity linking two (or more)
banked objects of the corpus that we did not put there by construction. A "sum = product," a closed form relating
two constants, a generating identity, an unexpected coincidence that turns out to be a theorem. Not a new value —
a **relation**.

## Order of events (Wilson's instruction, verbatim intent)

1. **[done]** Write the idea down (this file + memory).
2. **[now]** Prepare to compact.
3. **[after compaction]** *Yap about it* — discuss the shape of the hunt, what "a relationship we can express"
   means, what would count as our `π²/6` vs a triviality-by-construction.
4. **[then]** A **methodical, systematic hunt through the corpus** for internal relationships we can **express**.
   First we EXPRESS (write the candidate identity), then we TEST (gate it exactly, honest-negative discipline).
5. **[only after]** If we find and confirm any real relation, *then* we go on (act on it).

Wilson's own caveat, kept on the record: *"I'm probably reaching and that's fine."* This is exploratory. Most
candidate relations will be trivial-by-construction or coincidental; the discipline is to separate a genuine
internal identity from an artifact of how we defined things.

## The launchpad — the objects to relate (inventory to be built during the hunt)

⚠️ **Do NOT recall exact values here** ([[copy_dont_recall]]). The hunt's FIRST step is to grep the exact banked
constants out of the `result_*.md` artifacts into one inventory table. Below is the *object list* (by name + where
it lives), not their values.

- **The Plancherel ladder / value.** `S_∞` and its floor; the exact S-ladder `S_1, S_2, …`; `T_i`, `Λ_i`;
  `7/15`, `10/21`, `7/45`. (`result_P6*`, `result_SOLSTICE`, EPS JSON.)
- **The inverse-tree constant.** `λ = (3+√21)/6`, its minimal polynomial `3λ²−3λ−1=0`, `1/λ`, the branch ratio.
  (`result_BRANCH_BIAS`, `result_MOD2K_STATIONARY`.)
- **The kernel / spectral objects.** `K = {±2:2, 0:−1}`, `K̂ = 4cos2θ−1`, `|D̃|² = 1/(5−4cosθ)`, the subdominant
  rate `ρ_1`, the box-dim `D`. (`result_P6G`, `result_LATTICE`, `result_MAHLER`.)
- **The denominator arithmetic.** Mersenne `2^{2·3^{k−1}}−1`, `ord(2 mod q)`, LTE `v_3`. (`result_DENOM_OBSTRUCTION`,
  `result_CTILDE_EXTEND`, `result_MAHLER`.)
- **The family law.** `(q,p)`-Hydra `S_{k+1}/S_k → q(p−1)/(p+1)`, boundary `q(p−1)=p+1`, `c̃_q=(q−3)/q`, the
  mirror `S_∞(2,3)`. (`result_PHYDRA_FAMILY`, `result_MIRROR`.)
- **The stopping-time / archimedean side.** `ε_S ≈ log 4`, `E[L⁻] = μ/(1−q_esc)` with `μ=log(4/3)`, `q_esc`,
  `log_2 3`, `ord_3(2)`. (`result_TWO_WALLS`, `result_EPS_LATTICE`, `caravenna_doney_attempt`, `closed_form_findings`.)
- **The measure-type / Galois objects.** GARSIA dim→1, the level-1 difference-Galois group `G_1 = 𝔾ₐ` (Lean-proven),
  the pro-unipotent tower. (`result_GARSIA`, `result_QDIFF2/3`, `QDiff2.lean`.)

## What would count (pre-registered, so we don't fool ourselves)

- **Real (our π²/6):** an identity relating two objects that were defined by *different* constructions (e.g. the
  inverse-tree `λ` and the Plancherel `S_∞`; or the archimedean `ε_S`/`log 4` and the finite-place denominators;
  or a product-over-`⟨2⟩`-cascade that equals a sum we already have). Must be EXPRESSIBLE as a clean statement and
  then gated exactly.
- **Trivial-by-construction (not it):** a relation that just restates how we defined an object (e.g. `S_{i+1}=2T_i`,
  `T_i = T_{i−1}+Λ_i` — these are definitions/telescopes, already banked, not discoveries).
- **Coincidence (suspect):** a numerical near-equality at low precision with no structural reason — treat as a lead
  to prove or kill, never trusted from digits ([[facts_from_data]], the mirror-guardrail).

The Euler analogue we're imitating: two *representations* of one object, equated. The most Euler-shaped thing in
our corpus to look for is a **product form** (over primes? over the `⟨2⟩`-cascade? over levels?) for something we
currently only have as a **sum/limit** — because that sum↔product bridge is exactly the move that has never been
tried here.

**Status:** IDEA CAPTURED. Not acted on. Hunt begins after compaction + discussion.
