# Result — PHASE 2b Session Two: D1 COMMITTED (first hand-derived gap) + D3 lead (real-object invariant ray identified) + requests F/G pre-registered

**Date:** 2026-07-16. Follows `result_phase2b_s1.md` (Session One instruments A–E, C). This session Nathan **committed** the toy's closed-form gap (the derivation, D1), then a single sandbox gate advanced D3 (the real q=3 object). Discipline held: the r(λ) formula was derived from the block, not fitted; the D3 gate was pre-registered and **came back split** (one arm refuted, one arm hit), reported as such.

**Headline: `r(λ) = (1−λ²)/(1+λ²)` — DERIVED (not fit), five-for-five against the pre-published sweep. First hand-derived spectral gap of the L3 program. Maximality owed (Request F closes it by confirmation). And D3 sharpened: one member of the real q=3 coalescing pair is now CLOSED-FORMED as the folded diagonal moment `Σ w_r² ≈ 1/3 + (2/3)·2^{−D}`, matched to six digits at L=3; an eigenvalue braid toward the EP observed.**

Full theorem statement + proof of (a)–(d): `BRIEF_D1_TOY_GAP.md` (canonical). This file records the D3 gate, the retraction, the instrument lesson, and the standing requests.

---

## D1 — COMMITTED (summary; full brief is canonical)
`M(q,−1,λ)`, `u=1/(1+λ)`, `s=λ/(1+λ)`, `u+s=1`:
- **(a)** diagonal ray `{(1,1,0),(−1,−1,0)}` exactly closed; block `[[s²,u²],[u²,s²]]`, defect `2us` (exits `T=±2 ≢ 0 mod q` die at the gate, every odd q, every L).
- **(b)** `λ₁ = s²+u² = (1+λ²)/(1+λ)²` (Perron, 5/9 at λ=½); `λ₂ = s²−u² = s−u = −(1−λ)/(1+λ)` (**= first signed moment `m₋`**, −1/3 at λ=½).
- **(c)** `r(λ) = |λ₂|/λ₁ = (1−λ²)/(1+λ²)`. Exact rationals hit all five sweep points: 91/109, 21/29, 3/5, 8/17, 51/149.
- **(d)** `r<1` ∀λ∈(0,1]; boundary only at degenerate `λ→0` (no tower ⇒ no resonance — the toy's boundary is elsewhere than the real map's λ=½); `λ₂<0` oscillatory; `λ=tan(θ/2) ⇒ r=cosθ`, and λ=½ ⇒ (3/5,4/5) on the 3-4-5 triangle.
- **OWED — maximality:** no `e=−1` mode beats `|λ₂|`. Mechanism in hand (flip-flip pseudo-cycle at `u²=4/9>1/3` killed by carry within 2 steps). Finite reachable-graph argument, Nathan's pen. **Request F** confirms.

---

## D3 lead — the fold-refinement gate (ONE sandbox gate; pre-registered; came back SPLIT)
**Premise (from the toy):** the pair spectrum is **moments of the folded weights** (D1's `λ₂ = m₋`). The real q=3 fold refines with level: `D = 2·3^{L−1}` folded coordinates. **Menu check:** do the measured coalescing-pair values match a folded-moment family?

**✗ REFUTED — the COHERENT character-fold reading.** `|F₁|² = 0.806` at L=3 is **above** the Perron — nonsense as a mode. The naive "pair spectrum = coherent fold-moments" dies at the cost of one run, as designed.

**✓ HIT (twice, with a twist) — the INCOHERENT folded moment.** `Σ w_r²` matches **one member of the measured coalescing pair** at both levels:
| L | Σ w_r² (folded) | measured pair-member | match |
|---|---|---|---|
| 2 | 0.343915 | 0.3439 | 4 digits |
| 3 | 0.333336 | 0.333336 | **6 digits** |

Expansion: **`Σ w² ≈ 1/3 + (2/3)·2^{−D}`, `D = 2·3^{L−1}`** — one eigenvalue of the pair **rides the folded diagonal moment**, pacing at the **tower clock `2^{−D}`** (R14's doubly-exponential scale, finally in a legitimate role).

**★ THE TWIST (braid, flagged as observation):** it is the **lower** member at L=2 (`+1.35e-2` above 1/3) and the **upper** member at L=3 (`−1.0e-4` below 1/3). The **dynamical partner crosses the diagonal mode between levels** — an eigenvalue **braid en route to the EP**, the exceptional-point signature arriving unprompted. Braids near EPs are textbook; **this braid needs the L=4 point to be more than two dots** (Request G).

**⇒ D3 sharpens concretely:** half the coalescing pair is **identified in closed form** (`Σ w_r²` on the folded diagonal, tower-clocked). The theorem's remaining job reduces to the **dynamical partner's approach to 1/3**. Next session's target is a **T1-for-the-real-object**: find the invariant subspace of the actual q=3 operator whose block yields `Σ w_r²` exactly — the same way the toy's diagonal ray yielded `s²+u²`. If that ray exists (six-digit L=3 match says it does), the real operator **inherits the toy's proof architecture**: one exactly-solvable ray + one dynamical partner + a gap between them.

---

## Retraction + instrument lesson (this session)
- **RETRACTED:** T2's "carry skeleton **is** λ₂" mechanism (Session One's structural ansatz). The subdominant is the **antisymmetric diagonal mode** `s²−u²` on the 2-state ray (E), **not** a skeleton eigenvalue. The skeleton's real job is **maximality** (killing the `e=−1` block), not producing λ₂.
- **INSTRUMENT LESSON (banked, durable):** **mass sequences are subdominant-blind — never ESPRIT a total mass again.** The toy's total mass `1ᵀMᵏv₀` is **pure Perron** (subdominants mean-zero, `1ᵀrᵢ=0`); ESPRIT on it fits noise. This is what produced Q6's spurious "0.831" (corrected to 0.60 in Session One). Read subdominants from the **raw operator spectrum**, not from any total-mass sequence.

---

## Standing requests for the agent (consolidated, pre-registered)
- **Request F** — spectral radius of the `e=−1` sub-block **alone**; `L=1,2`; `λ = 0.3 / 0.5 / 0.7`. **Pre-registered: strictly below `|λ₂|` at every point** (expect `us`-to-`s²` scale). Closes D1's maximality by **confirmation** (the graph proof becomes a check, not an exploration).
- **Request G** — the **L=4 gap** via E's localization trick applied to the **real** pair: find where the coalescing eigenvectors live at `L=2,3`, build the **reduced operator on that subspace** at `L=4`. (The full 39366-state solve **walled** in Session One — Request C; the reduced one should not.) **G carries double freight:** the rate law's 4th point (D3) **and** the braid's 3rd dot (turns two dots into a curve).

---

## Session ledger (honest)
- **D1 COMMITTED** — `r(λ)=(1−λ²)/(1+λ²)`, five-for-five against a curve published before the derivation, **first hand-derived spectral gap of the program.** Maximality owed pending F.
- **D3 ADVANCED** — one pair-member closed-formed (`Σ w_r²`, tower-clocked `2^{−D}`), braid observed (partner crosses 1/3 between L=2 and L=3), **invariant-ray target named** (T1-for-the-real-object).
- **One mechanism RETRACTED** (T2's skeleton-as-λ₂). **One instrument lesson BANKED** (mass sequences subdominant-blind).
- **The toy did its job:** the technique that cracked it — find the exactly-invariant ray, diagonalize the block, fight for maximality — is now the **template walking into the real theorem.**

## Not at stake
R1–R46, Phases 0/1/2a, Session One (A–E, C). No `r_q` value changes; no theorem below L3 touched. This is Session-Two derivation + one D3 gate.

_Reporting discipline: r(λ) derived from the 2×2 block then met the pre-published table (not fitted). The D3 gate was pre-registered and came back split — the coherent arm REFUTED (|F₁|²>Perron, nonsense), the incoherent arm HIT (six-digit L=3 match); both reported. The braid is flagged as an observation needing L=4 (Request G), not asserted as a law. Maximality is stated as OWED, not claimed. The T2 retraction and the mass-ESPRIT lesson are disclosed as corrections._
