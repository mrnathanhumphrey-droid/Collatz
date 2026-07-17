# Result — L-A (no-return) and L-B (gauge factorization, k=0 exact) CLAUDE-VERIFIED on the real q=3 operator.

**Date:** 2026-07-16. Independent verification (requested) of the two skeleton lemmas previously banked as Nathan-reported. Probe `probes/probe_phase2b_LALB.py`, log `logs/probe_phase2b_LALB_log.txt`. Operator `M = build_M_gen(3, L, 2, [λ^δ])`, λ=½. State classes: Δ={(a,a,0)}, C (carried-diagonal)={(a,a,γ≠0)}, O (off-diagonal)={(a,b,γ): a≠b}. Edge src→dest iff M[dest,src]≠0.

**Headline: BOTH CONFIRMED. L-A — carried-diagonal never returns to Δ (one-step M[Δ,C]=0 exactly AND transitive BFS reaches 0 Δ-states, at L=2 and L=3). L-B — the k=0 co-invariant LEFT eigenvector ℓ₀ factors EXACTLY through the (ρ,γ) gauge quotient (within-orbit deviation 9.35e-16), is χ₀ (constant) on Δ and zero on carried-diagonal, and this holds for k=0 ONLY (k=1,3 deviation 1.0). Gauge = diagonal shift (a,b,γ)→(sa,sb,γ): 54 orbits at L=2, matching Nathan's "≤54 gauge classes."**

## L-A — no-return (C ↛ Δ). ✅ CLAUDE-VERIFIED at L=2 and L=3
| L | dim | \|Δ\| | \|C\| | \|O\| | one-step w[Δ←C] | w[Δ←O] | BFS from C: Δ reached / total reached |
|---|---|---|---|---|---|---|---|
| 2 | 324 | 6 | 48 | 270 | **0.0** | 3.9 | **0** / 126 |
| 3 | 8748 | 18 | 468 | 8262 | **0.0** | 12.0 | **0** / 3078 |

- **One-step:** `M[Δ, C] = 0` exactly — no carried-diagonal state transitions directly into Δ.
- **All-step:** BFS forward from *all* carried-diagonal states (48 at L=2, 468 at L=3) reaches **zero** Δ states (while spreading to 126 / 3078 states total). **No path C → Δ at any length.**
- **Precision of the claim confirmed:** Δ *does* receive flux from the off-diagonal sector O (w[Δ←O] = 3.9 / 12.0 ≠ 0) — the leak returns to Δ through O, not through C. So the no-return is specifically about **carried-diagonal**, exactly as Nathan stated; Δ is not blanket-closed. This is the correct and needed form: from C the dynamics reaches only C∪O, never Δ, so the co-invariant ℓ_k (=0 on C) never sees ℓ|_Δ from a C-transition.

## L-B — gauge factorization, k=0 exact. ✅ CLAUDE-VERIFIED (on the LEFT / co-invariant eigenvector)
Gauge = simultaneous diagonal phase shift `g_s: (a,b,γ) → (sa, sb, γ)`, s∈⟨2⟩; invariant label `(ρ = a·b⁻¹ mod 3^L, γ)`. **Total gauge orbits at L=2 = 54 — matches Nathan's "object reduced to ≤54 gauge classes."** Reachable-from-Δ = 17.

| k | c_k | LEFT eigvec within-orbit dev | gauge-invariant? |
|---|---|---|---|
| **0** | 0.343915 | **9.35e-16** | **YES — factors exactly** |
| 1 | 0.0794+0.2749j | 1.00 | no |
| 3 | −0.206349 | 1.00 | no |

- **k=0 factors exactly:** the co-invariant functional ℓ₀ (left eigenvector for c₀) is constant on each (ρ,γ) orbit to machine zero ⇒ it descends to the reduced (ρ,γ)-chain. **k=0 only** — k=1,3 are not gauge-invariant (dev 1.0), so the factorization is special to the trivial character, i.e. "k=0 exact."
- **Co-invariance structure confirmed:** |ℓ₀| = 1.000 (constant = χ₀) on Δ; |ℓ₀| ≤ 1.3e-15 (zero) on carried-diagonal C. Exactly Nathan's description (character on Δ, zero on diagonal-with-carry, dressed off-diagonal).
- **Method note (caught):** the first pass tested the RIGHT eigenvector, which is NOT gauge-invariant (dev 0.806) — the co-invariance mechanism lives on the LEFT eigenvector. Corrected; the banked probe tests the left (co-invariant) eigenvector.

## Provenance upgrade
- **L-A:** Nathan-reported → **CLAUDE-VERIFIED** (one-step + transitive BFS, L=2 and L=3).
- **L-B:** Nathan-reported → **CLAUDE-VERIFIED** for its concrete core (k=0 co-invariant eigenvector gauge-factorizes exactly; 54-orbit count matches). The full skeleton (three lemmas, keystone) now rests on two Claude-verified lemmas + the circulant-family completeness (already Claude-verified).

## Flag — the "21 closed classes"
Nathan's next-session object is "the 21 closed classes of the reduced (ρ,γ)-chain." This probe measured **17 reachable-from-Δ orbits at L=2** (of 54 total). 17 ≠ 21: these are different counts (reachable-from-Δ at L=2 vs. Nathan's "closed" classes, likely a communicating-class notion and/or a different level). **Not a contradiction — different quantities; the 21-closed-class enumeration is the next-session hand target and was not reproduced here.** The 54-total-orbit count (the reduction object) IS reproduced.

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, the H-gates, the J-refutation. No `r_q` value changes. The invariant that proves no-return as a *theorem* (the conserved quantity on the (ρ,γ)-chain) is next session's hand work; this probe verifies the no-return holds numerically at L=2,3 and that k=0 factors — it does not supply the invariant.

_Reporting discipline: L-A verified by two independent measures (one-step weight = 0 AND transitive BFS = 0) at two levels; the precise scope (C not O) is stated, since Δ does receive from O. L-B's first pass used the wrong (right) eigenvector and was corrected to the left/co-invariant one — disclosed. The 54-orbit count matches Nathan's ≤54; the 21-closed-class figure is flagged as a distinct, un-reproduced next-session quantity, not asserted._
