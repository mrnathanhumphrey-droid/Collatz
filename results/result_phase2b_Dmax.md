# Result — LEMMA D1-MAX (user's proof) BANKED + gated. THEOREM D1 COMPLETE. One decorative deviation (index law), lemma intact.

**Date:** 2026-07-16. Nathan wrote the pen-and-paper proof of **LEMMA D1-MAX** (acyclicity ⇒ nilpotence of the e=−1 carry graph), the one owed piece of THEOREM D1. This file banks the proof and records two pre-registered gates on its *new* falsifiable predictions (beyond Request F's L∈{1,2}). Probe `probes/probe_phase2b_Fcor.py`, log `logs/probe_phase2b_Fcor_log.txt`.

**Headline: the LEMMA is CONFIRMED — the e=−1 block is nilpotent (ρ₋=0) for q≥5 at L=1,2,3; the surviving T-alphabet is EXACTLY {2, q^L−2}; and the q=3 corollary is exact (ρ₋=s², maximality crossover precisely at λ=1/√2). ONE deviation: the nilpotency INDEX is not 2L — measured {2,4,5} at L={1,2,3}, so the proof's BOUND (index ≤ 2L) holds but the "= 2L exactly" flourish is false. The lemma needs only nilpotence, so THEOREM D1 closes.**

## The proof (Nathan; banked verbatim in structure)
Setting: in the e=−1 sector (b=−a), sector-preserving moves are stay-stay (s²) and flip-flip (u²), with T ∈ {2, q^L−2} always; gate γ ≡ −T mod q; carry γ′=(γ+T)/q exact (both branches stay < q^L).
- **Step 1 — cycle equation.** n surviving steps ⇒ `q^n·γ_n = γ_0 + Σ_{i<n} t_i q^i`, `t_i ∈ {2, q^L−2}`. A cycle forces `γ_0(q^n−1) = Σ t_i q^i`. Pure-A (all t=2): `γ_0 = 2/(q−1) ∈ ℤ ⟺ q=3`. Pure-B (all t=q^L−2): `γ_0=(q^L−2)/(q−1) ∈ ℤ ⟺ (q−1)|1` — never. Mixed: killed by the digit automaton.
- **Step 2 — digit automaton.** Base-q, d₀ least significant. Move A (t=2, needs d₀=q−2): rolls to 0, +1 carry to d₁, shift down, injects 0 at top; surviving paths have no deeper cascade (d₁=q−1 ⇒ next read 0 ∉ {2,q−2} = death) ⇒ **A adds +1 to the next read digit**. Move B (t=q^L−2, needs d₀=2): subtract 2 (no borrow), shift, injects 1 at top. After L steps all original digits flush; the digit read at step i is `r_i = inj_{i−L} + [step i−1 was A]`, inj ∈ {0(A),1(B)} ⇒ `r_i ∈ {0,1,2}`. Gate demands `r_i ∈ {2, q−2}`.
- **Step 3 — contradiction (q≥5).** q−2 ≥ 3 unreachable ⇒ every late read must be r_i=2 ⇒ inj=1 (step was B) AND step i−1 was A — "all B" and "all A" simultaneously. Contradiction. **No infinite surviving path ⇒ the sector is NILPOTENT.** ∎
- **Corollary (q=3).** Alphabet {2, q−2}={2,1}; r=1 IS reachable (A-inject 0 + A-carry 1) ⇒ the pure-A fixed point γ=1 survives at weight s². Maximality at q=3 holds for `s² < |λ₂| ⟺ 2λ² < 1 ⟺ λ < 1/√2` (covers λ=½). Toy-internal; the grid {5,7,13} never sees it.

## Gate P-F4 — nilpotence, T-alphabet, index (q=5,7; L=1,2,3; λ=½)
| q | L | dim | ρ₋ | T-alphabet | pred {2,q^L−2} | nilpotency index | 2L |
|---|---|---|---|---|---|---|---|
| 5 | 1 | 10 | 0 | {2, 3} | {2,3} ✓ | 2 | 2 |
| 5 | 2 | 50 | 0 | {2, 23} | {2,23} ✓ | 4 | 4 |
| 5 | 3 | 250 | 0 | {2, 123} | {2,123} ✓ | **5** | 6 |
| 7 | 1 | 14 | 0 | {2, 5} | {2,5} ✓ | 2 | 2 |
| 7 | 2 | 98 | 0 | {2, 47} | {2,47} ✓ | 4 | 4 |
| 7 | 3 | 686 | 0 | {2, 341} | {2,341} ✓ | **5** | 6 |

- **✅ CONFIRMED:** ρ₋=0 (nilpotent) and T-alphabet = {2, q^L−2} EXACTLY, at L=3 for both primes. Steps 1–3 of the proof are gate-validated one level beyond Request F.
- **⚠️ DEVIATION (reported as a deviation):** the nilpotency **index is {2, 4, 5}** for L={1,2,3} — **not 2L = {2,4,6}**. Power-norm ladder (identical prefix, one extra nonzero power per level): `L=1: 4.4e-1,0`; `L=2: 4.4e-1,2.0e-1,8.8e-2,0`; `L=3: …,9.8e-3,0`. The proof's **bound `index ≤ 2L` HOLDS** (2≤2, 4≤4, 5≤6) and `Msub^{2L}=0` holds — **so the lemma (nilpotence ⇒ ρ₋=0 ⇒ maximality) is untouched.** Only the "= 2L exactly / index-2L growth" refinement is false. **Correction of provenance:** that "= 2L" came from *Request F's* two-point fit (L=1,2 → 2,4); it is a small-window extrapolation and is retracted (see `result_phase2b_F.md` corrected). *The exact index is immaterial to the lemma and is NOT re-fit here (per the standing no-small-window rule).*

## Gate P-F5 — q=3 corollary (the toy exception): ρ₋ = s², crossover at λ=1/√2
| λ | ρ₋ (L=2, =L=3) | s²=λ²/(1+λ)² | \|λ₂\|=(1−λ)/(1+λ) | maximality | predicted (2λ²<1) |
|---|---|---|---|---|---|
| 0.30 | 0.053254 | 0.053254 ✓ | 0.538462 | HOLDS | holds ✓ |
| 0.50 | 0.111111 | 0.111111 ✓ | 0.333333 | HOLDS | holds ✓ |
| 0.60 | 0.140625 | 0.140625 ✓ | 0.250000 | HOLDS | holds ✓ |
| 0.70 | 0.169550 | 0.169550 ✓ | 0.176471 | HOLDS | holds ✓ |
| **0.7071** | **0.171571** | 0.171571 ✓ | **0.171578** | HOLDS (by 7e-6) | holds ✓ |
| 0.72 | 0.175230 | 0.175230 ✓ | 0.162791 | **FAILS** | fails ✓ |
| 0.80 | 0.197531 | 0.197531 ✓ | 0.111111 | **FAILS** | fails ✓ |

**✅ CONFIRMED exactly (7/7, both L).** `ρ₋(M(3,−1,λ)) = s²` to machine precision — the pure-A fixed point γ=1 at weight s² is the surviving e=−1 mode, as the corollary claims. The maximality crossover lands **precisely at λ=1/√2 ≈ 0.70711**: holds below (at 0.7071 by 7e-6), fails above (0.72). **λ=½ (Syracuse) is safely inside.** The q=3 toy is genuinely different from q≥5 — it has a surviving e=−1 cycle — but the ray's λ₂ still dominates at the Syracuse weight.

## ★ THEOREM D1 — COMPLETE
For M(q,−1,λ), **q≥5, any L**: the spectrum is exactly {λ₁=(1+λ²)/(1+λ)², λ₂=−(1−λ)/(1+λ)} on the diagonal ray, plus a **nilpotent** off-diagonal block (ρ₋=0). Hence **`r(λ)=(1−λ²)/(1+λ²) < 1` for all λ∈(0,1]** — gap open at every weight, prime, level. Derived, gated five-for-five on the sweep, maximality **proven** (Nathan) and gate-confirmed (Steps 1–3 at L=3; q=3 corollary exact). The first fully-closed spectral gap of the program.

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1(a)–(d), Request F (except its "index=2L" flourish, now corrected to ≤2L). No `r_q` value changes. The D3 lead / real q=3 object is untouched (next session's target: the T1-for-the-real-object invariant ray).

_Reporting discipline: the lemma's core (nilpotence ⇒ ρ₋=0 ⇒ maximality) and its T-alphabet + q=3 corollary are gate-confirmed and reported as passes. The one deviation — index ≠ 2L (measured {2,4,5}) — is reported AS a deviation, traced to Request F's two-point extrapolation, and NOT re-fit (small-window rule). The proof's actual bound (≤2L) holds, so THEOREM D1 is recorded as complete. q=3 corollary confirmed to machine precision including the 1/√2 boundary._
