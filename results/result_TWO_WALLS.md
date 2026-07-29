# RESULT — TWO WALLS / Gate 2 (Option 2): the ladder wall (E[L⁻]) is archimedean, the spectral wall (S_∞) is finite-place; they are two places of one idele. Gate 2 confirms (2026-07-28)

**Structural verdict (the "two walls are two places of one idele" reasoning) is penned and banked in the memory filesystem** — Wilson's `/areas/collatz-height-biextension.md`, section "THE TWO WALLS ARE TWO PLACES OF ONE IDELE." **This repo file carries the Gate 2 confirming DATA only, and points there for the reasoning — do not duplicate the derivation here (the two copies would drift).** Probe: `probes/probe_gate2_escape.py`.

## The question (Option 2)
Is the **Gelfond–Schneider** obstruction on the ladder side (no closed form for the mean descending-ladder height `E[L⁻]`, blocked by `log₂3 ∉ ℚ` in the Wiener–Hopf factorization) the *archimedean face* of the **infinite unipotent depth** that MAHLER proved on the spectral side (S_∞ = infinite-order Mahler; cyclotomic denominators `2^{2·3^{k−1}}−1`, `⟨2⟩`-order in `(ℤ/3^k)*`)?

**Penned verdict (see the /areas file):** not (A) shallow "both hard from `log₂3`" nor (B) "same number," but the sharp middle — **`ord₃(2)` (finite-place face) and `log₂3` (archimedean face) are the two places' readings of one idele "2," tied by the product formula `Σ_v log|2|_v = 0`.** The spectral wall is the finite-place shadow (accumulated `ord₃(2)` data up the tower = infinite depth); the ladder wall is the archimedean shadow (Gelfond–Schneider on `log₂3` at ∞). Dual faces, related by the product formula, but each obstruction is local to its place — killing one doesn't kill the other. Pre-registered predictions: both gates come back showing the ladder side archimedean, the spectral side finite-place.

## GATE 2 — the confirming data (escape probability across the (q,p)-Hydra family)
The ascending-ladder escape probability `q_esc = P(σ⁺<∞)`, via the confirmed identity `E[L⁻] = μ/(1−q_esc)`. **Analytic reduction:** each Syracuse step is `X = log q − v·log p`, so `S_n = n·log q − V_n·log p` (`V_n=Σvᵢ`) and

> **`S_n > 0 ⟺ V_n/n < log_p(q)`**  — so `q_esc = P(∃n: V_n/n < log_p q)` factors through the archimedean ratio `α = log_p(q)` and the p-geometric jump-law. No `⟨2⟩`-order / finite-place content anywhere.

MC confirmation (3M walks/member):
```
   (q,p)    α=log_p q   drift     q_esc     note
   (3,2)     1.58496   -0.2877   0.7137    Collatz — EXACT match to on-record 0.7137 (E[L⁻]=1.00466, μ=0.287682)
   (2.5,2)   1.32193   -0.4700   0.5433    p=2 slice, smooth/monotone in α
   (3.5,2)   1.80735   -0.1335   0.8556
   (3.9,2)   1.96347   -0.0253   0.9672
   (2,2)     1.00000   -0.6931   0.0000    α=1 (q=p) boundary: walk can't strictly ascend → q_esc=0
   (4,3)     1.26186   -0.2616   0.7689
   (5,4)     1.16096   -0.2390   0.8040
   (9,4)     1.58496   +0.3488   1.0000    same α as (3,2), DIFFERENT q_esc — the tell
```

## The (9,4) tell — the clean kill of thesis (B)
`(3,2)` and `(9,4)` share `α = log₂3 = log₄9 = 1.58496` **exactly**, yet `q_esc = 0.714` vs `1.000` — because the p-geometric jump-law differs (and `(9,4)` is supercritical). So **`q_esc = F(α, p)`: an archimedean `α` modulated by the p-jump-law, with zero `⟨2⟩`-multiplicative-order content, not algebraic, and unrelated to the `q(p−1)/(p+1)` graded rate** (which was the spectral/finite-place object). The escape constant is governed by `log_p q` at ∞ — exactly "governed by `log₂3` at the archimedean place, not `ord₃(2)` at the finite place."

## Verdict
- **Gate 2 CONFIRMS thesis (A): the ladder wall is archimedean, dual to the finite-place spectral wall.** Pre-registered prediction upheld, converted to a checked fact (the mirror move).
- **Gate 1 (the heavier 41-role / Wiener–Hopf read) is NOT needed** — per the pre-registered plan, only run on a wobble toward (B); Gate 2 didn't wobble. (Also flagged: the naive "41 divides `2²⁰−1`" is trivially universal — every odd prime divides some `2ⁿ−1` — so Gate 1's sharp form needs the actual factorization/Ostrowski structure, not a divisibility check.)
- **Structure-vs-value split, now confirmed at FOUR sites** (three falsifiable tests that held): biextension (structure homed, height walled) · p-Hydra law (graded functional equation real, value untouched) · MIRROR (involution preserves the graded, breaks the value, `0.475≠0.459`) · two-walls (ladder archimedean, spectral finite-place, Gate-2-confirmed). One fact: the shared structure lives in the associated-graded / finite place, the value lives in the extension / archimedean place, tied by the product formula but evaluated by neither.
- **Reasoning home:** `/areas/collatz-height-biextension.md` (memory filesystem). **Not at stake:** S_∞≈0.475, MAHLER, MIRROR, PHYDRA, GARSIA, DENOM, R1–R30. Neither value delivered — E[L⁻] stays Gelfond–Schneider-blocked at ∞, S_∞ stays infinite-unipotent-depth-blocked at the finite place.
