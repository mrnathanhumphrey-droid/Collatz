# T_V_RECURSION — Phase 1 attempt to derive M(g) → M(g') recursion

**Date:** 2026-05-12. Seventh spectral probe (Phase 3b/3c follow-up from CROSS_FREQ_DISPOSITION H_CROSS_CLOSES_ON_ENLARGED_SPAN). Wilson reporting to Nathan.

This document attempts the Phase 1 derivation of the level-(n) → level-(n+1) recursion induced by Tao on the moment family

  M_n^{ab}(g, c) := Σ_{ξ ∈ (Z/3^n)^×, ξ ≡ c mod 3} e^{-2πi ξ·ẽ_g/3^n} · μ̂_n^a(ξ) · μ̂_n^{b*}(ξ·2^{-g} mod 3^n)

with ẽ_g := (1 - 2^{-g})/3 (well-defined 3-adic integer for g even ≥ 0; ẽ_0 = 0).

---

## Verdict at the end of Phase 1

> **The cross-freq materials' enlarged span V_M = span{M_n^{ab}(g, c) : g ∈ {0, 2, 4, ...}} does NOT close under Tao's iteration.** The recursion at level n+1 generates moments with extra phase offsets θ_{v,g} := 2^v·ẽ_g/3 (mod 3^n) that are NOT in the V_M family for generic (v, g). To close, V_M would have to be enlarged to a phase-parameterized family V'_M (with dim growing in n), which is incompatible with the brief's "finite truncation at g_max" framing.
>
> **Disposition (Phase 1 alone):** **H_M_RECURSION_UNDERSPECIFIED**, with substantive structural content — the obstruction is identified explicitly (phase offsets θ_{v,g}), not hand-waved. Phase 2–6 are not executable on V_M^{(g_max)} as the brief envisions.

---

## 1. Setup

We want to expand M_{n+1}^{ab}(g, c) via Tao's recursion and collect terms as a linear combination of {M_n^{a'b'}(g', c')} (over some index set).

By Tao's recursion + R66 class flow (CROSS_FREQ_PHASE1_EXPANSION §0):

  μ̂_{n+1}^a(ξ) = Σ_{v ∈ V_a} 2^{-v} e^{-2πi ξ·2^{-v}/3^{n+1}} μ̂_n(ξ·2^{-v} mod 3^n)

with V_+ = {2, 4, 6, ...} (even ≥ 2), V_- = {1, 3, 5, ...} (odd ≥ 1).

  μ̂_{n+1}^{b*}(η) = Σ_{v' ∈ V_b} 2^{-v'} e^{+2πi η·2^{-v'}/3^{n+1}} μ̂_n^*(η·2^{-v'} mod 3^n).

In the M_{n+1}^{ab}(g, c) integrand, the b-leg's argument is η = ξ·2^{-g} mod 3^{n+1}.

Substituting both legs into M_{n+1}^{ab}(g, c):

  M_{n+1}^{ab}(g, c)
    = Σ_{v ∈ V_a, v' ∈ V_b} 2^{-v-v'} · S_n(v, v', g, c, a, b)

with the inner level-(n+1) sum

  S_n(v, v', g, c, a, b)
    := Σ_{ξ ≡ c mod 3, ξ ∈ (Z/3^{n+1})^×}
         e^{-2πi ξ·D_{v,v',g}/3^{n+1}}
         · μ̂_n(ξ·2^{-v} mod 3^n) · μ̂_n^*(ξ·2^{-g-v'} mod 3^n)

where the consolidated phase exponent is

  **D_{v,v',g} := ẽ_g + 2^{-v} - 2^{-g-v'}**

(ẽ_g from the M_{n+1} definition; 2^{-v} from A_v(ξ); -2^{-g-v'} from A_{v'}^*(ξ·2^{-g}).)

---

## 2. Lift-fiber sum at level n+1 → level n

Decompose ξ ∈ (Z/3^{n+1})^× as ξ = u + j·3^n with u ∈ (Z/3^n)^× (u ≡ c mod 3) and j ∈ {0, 1, 2}. The μ̂_n arguments depend only on u (not j) because (u + j·3^n)·2^{-v} ≡ u·2^{-v} mod 3^n for any j.

The phase factor splits:

  e^{-2πi ξ·D/3^{n+1}} = e^{-2πi u·D/3^{n+1}} · e^{-2πi j·D/3}.

The j-sum gives 3 if 3 | D, else 0.

**Survival condition:** Σ_j e^{-2πi j·D_{v,v',g}/3} ≠ 0 ⟺ 3 | D_{v,v',g}.

When 3 | D: write D = 3·D̃ where D̃ ∈ Z[1/2] (well-defined mod 3^n). The surviving inner sum:

  S_n = 3 · Σ_{u ∈ (Z/3^n)^×, u ≡ c mod 3}
           e^{-2πi u·D̃/3^n} · μ̂_n(u·2^{-v} mod 3^n) · μ̂_n^*(u·2^{-g-v'} mod 3^n).

---

## 3. Unit shuffle and emergence of shift index G

Substitute s := u·2^{-v} mod 3^n (bijection on (Z/3^n)^×). Then u = s·2^v mod 3^n; second arg: u·2^{-g-v'} = s·2^{v-g-v'} = s·2^{-G} where

  **G := v' + g - v** ∈ Z (can be negative, zero, or positive).

Class flow: u ≡ c mod 3 ⟺ s·2^v ≡ c mod 3 ⟺ s ≡ c·2^{-v} ≡ c·(-1)^v mod 3. Define

  **c̃ := (-1)^v · c mod 3 ∈ {1, 2}** (v even: c̃ = c; v odd: c̃ = 3-c).

(Reminder: 2 ≡ -1 mod 3.)

The phase exponent transforms: -u·D̃/3^n = -s·2^v·D̃/3^n. Compute 2^v·D̃:

  2^v·D̃ = 2^v·D/3 = 2^v·(ẽ_g + 2^{-v} - 2^{-g-v'})/3
         = (2^v·ẽ_g + 1 - 2^{v-g-v'})/3
         = (2^v·ẽ_g + 1 - 2^{-G})/3
         = ẽ_G + 2^v·ẽ_g/3                   [using ẽ_G = (1 - 2^{-G})/3]

(Here G can be any integer; ẽ_G is well-defined in Q[1/2] regardless of sign, but its 3-adic valuation depends on G — see §4.)

Inserting this into the s-sum:

  S_n = 3 · Σ_{s ≡ c̃ mod 3, s ∈ (Z/3^n)^×}
           e^{-2πi s·(ẽ_G + θ_{v,g})/3^n}
           · μ̂_n(s) · μ̂_n^*(s·2^{-G} mod 3^n)

with the **phase offset**

  **θ_{v,g} := 2^v · ẽ_g / 3 = 2^v·(1 - 2^{-g})/9** ∈ Q[1/2] (well-defined mod 3^n).

Splitting μ̂_n = μ̂_n^+ + μ̂_n^- to extract classes:

  S_n = 3 · Σ_{a',b' ∈ {+,-}} M̃_n^{a'b'}(G, c̃; ẽ_G + θ_{v,g})

where the **generalized phase-twisted moment**

  **M̃_n^{a'b'}(G, c; φ) := Σ_{s ≡ c mod 3, s ∈ (Z/3^n)^×} e^{-2πi s·φ/3^n} · μ̂_n^{a'}(s) · μ̂_n^{b'*}(s·2^{-G} mod 3^n).**

When φ = ẽ_G (the canonical choice), M̃ reduces to M:

  M̃_n^{a'b'}(G, c; ẽ_G) = M_n^{a'b'}(G, c).

For φ ≠ ẽ_G, M̃ is a NEW moment object outside V_M.

---

## 4. The core obstruction: θ_{v,g} is NOT canceled by ẽ_G

The combined phase ẽ_G + θ_{v,g} would reduce M̃ back into V_M only if θ_{v,g} ≡ 0 mod 3^n or θ_{v,g} = ẽ_{G''} - ẽ_G for some G''.

**Claim 1.** θ_{v,g} ≢ 0 mod 3^n in general.

  v_3(θ_{v,g}) = v_3(2^v·(1 - 2^{-g})/9)
              = v_3(1 - 2^{-g}) - 2
              = v_3(2^g - 1) - 2
              = (1 + v_3(g/2)) - 2     [for g even ≥ 2, by LTE_3]
              = v_3(g/2) - 1.

For g = 2 (the dominant cross-frequency mode): v_3(θ_{v,2}) = -1, which means θ_{v,2} has a 3 in the denominator (not a 3-adic integer).

But θ_{v,g} appears in the exponent over 3^n, so it must be interpreted mod 3^n. Concretely: write D̃ = D/3 as an element of Z[1/2] / 3^n·Z[1/2], and 2^v·D̃ ∈ Z[1/2] mod 3^n. The "phase factor" e^{-2πi s·θ/3^n} is well-defined as long as θ ∈ Z[1/2] mod 3^n (clearing the 2-denominator: 2 is a unit mod 3, so s·θ ∈ Z[1/2] mod 3^n is fine).

**For our concrete cases, however**, the phase offset DOES contribute a non-trivial twist that doesn't reduce to ẽ_{G''} for any G''.

**Claim 2 (worked example).** For g = 2, v = 2, v' = 3 (so a = +, b = -, G = 3, c̃ = c):

  D_{v,v',g} = ẽ_2 + 2^{-2} - 2^{-5} = 1/4 + 1/4 - 1/32 = (8 + 8 - 1)/32 = 15/32.

  D̃ = D/3 = 5/32.

  2^v·D̃ = 4·5/32 = 20/32 = 5/8.

  ẽ_G = ẽ_3 = (1 - 1/8)/3 = (7/8)/3 = 7/24.

  θ_{v,g} = 2^v·D̃ - ẽ_G = 5/8 - 7/24 = (15 - 7)/24 = 8/24 = 1/3.

But ẽ_G' = (1 - 2^{-G'})/3 = (2^{G'} - 1)/(3·2^{G'}). For this to equal 1/3: 2^{G'} - 1 = 2^{G'} ⟹ -1 = 0. No solution.

So θ_{v=2, g=2} = 1/3 (a non-3-adic-integer rational) is NOT expressible as a difference of two ẽ_G's. The phase 5/8 = ẽ_3 + 1/3 is in NO ẽ_G family. Hence S_n(v=2, v'=3, g=2, ...) produces a generalized M̃ with phase 5/8 that is NOT in V_M.

**Claim 3.** As (v, v', g) vary over the surviving (3 | D) tuples, θ_{v,g} takes finitely many distinct values mod 3^n at each fixed n. Since 2^v cycles mod 3^n with period 2·3^{n-1}, the v-set yields at most 2·3^{n-1} distinct 2^v values mod 3^n, and θ_{v,g} = 2^v·ẽ_g/3 mod 3^n similarly takes at most 2·3^{n-1} distinct values for each fixed g.

So the closure space, if we insist on closing, has dimension growing with n — at level n, the relevant phase-twisted moments form a basis of size O(N_g · 3^{n-1}) per (a, b, c). This is **not a fixed dimension; it grows with n**. The operator T_V then doesn't have a uniform finite-dim truncation at any g_max independent of n.

---

## 5. Survival condition (3 | D) as a parity / class diagnostic

Reading 3 | D ⟺ ẽ_g + 2^{-v} - 2^{-g-v'} ≡ 0 mod 3.

Mod 3 we have 2 ≡ -1, so 2^{-v} ≡ (-1)^v, 2^{-g-v'} ≡ (-1)^{g+v'}. For g even (the relevant case), (-1)^g = 1, so 2^{-g-v'} ≡ (-1)^{v'}.

Mod 3 reduction of ẽ_g: for g = 2, ẽ_2 = 1/4 ≡ 1 mod 3 (inv(4) mod 3 = 1). For g = 4: ẽ_4 = 5/16 ≡ 5·inv(16) mod 3 = 5·1 = 2 mod 3. For g = 6: ẽ_6 = 21/64; v_3(ẽ_6) = v_3(21/64) = 1 (since v_3(21) = 1, v_3(64) = 0). So ẽ_6 ≡ 0 mod 3.

For g = 2 (ẽ_g ≡ 1 mod 3):
  D ≡ 1 + (-1)^v - (-1)^{v'} mod 3.
  - v even, v' even: 1 + 1 - 1 = 1 ≢ 0. **Kills.**
  - v even, v' odd: 1 + 1 - (-1) = 3 ≡ 0. **Survives.**
  - v odd, v' even: 1 + (-1) - 1 = -1 ≢ 0. **Kills.**
  - v odd, v' odd: 1 + (-1) - (-1) = 1 ≢ 0. **Kills.**

For g = 4 (ẽ_g ≡ 2 mod 3):
  - v even, v' even: 2 + 1 - 1 = 2 ≢ 0. Kills.
  - v even, v' odd: 2 + 1 - (-1) = 4 ≡ 1 ≢ 0. Kills.
  - v odd, v' even: 2 + (-1) - 1 = 0. **Survives.**
  - v odd, v' odd: 2 + (-1) - (-1) = 2 ≢ 0. Kills.

For g = 6 (ẽ_g ≡ 0 mod 3, v_3 ≥ 1, refined to ẽ_6 ≡ 3 mod 9):
  - All four parity combos have D mod 3 = 0 + (parity diff). Need a refined mod-9 check.

This already shows: the recursion DOES NOT couple same-class P-moments to same-class M(g)-moments at level n+1 → n. Instead:

- For g = 2: the surviving (v even, v' odd) pairs contribute to M_{n+1}^{+−}(g=2, c) (a=+, b=-).
- For g = 4: the surviving (v odd, v' even) pairs contribute to M_{n+1}^{-+}(g=4, c).

**The g=2 recursion thus produces (P^{+-})-class moments at level n+1, NOT (P^{++}) or (P^{--}).** But R76 §11 / cross_freq §3 said P^{+−} = 0 structurally (algebraic identity at g = 0). The cross_freq probe asserted this propagates to all g ≥ 2 even — but the survival pattern here suggests M^{+-}(g=2, c) could be NON-zero.

This means either:
- The cross_freq probe's claim "(1, 4)-direction preservation across all g" relies on summing over both (a,b) = (+,-) and (-,+) which DO have non-trivial g ≥ 2 contributions that cancel to give the (1, 4) projection. OR
- My derivation has missed a structural detail.

I lean toward the former: the cross_freq probe sums X̄_n(g, c) := Σ_{a,b} M_n^{ab}(g, c), and the (1, 4) structure comes from this CLASS-SUMMED moment, NOT from individual class-resolved moments.

This is consistent with how X̄_n(g, c) appears in the off-diagonal correction (CROSS_FREQ_PHASE1_EXPANSION §7).

---

## 6. The cascade: g' and v' relationship

Given (v ∈ V_a, v' ∈ V_b, g, n+1 → n recursion), the new shift index is G = v' + g - v.

For g = 2:
- (a=+, b=-) surviving: v ∈ {2, 4, 6, ...} (even ≥ 2), v' ∈ {1, 3, 5, ...} (odd ≥ 1). G = v' + 2 - v.
  - v = 2, v' = 1: G = 1 (odd, but ẽ_1 has v_3 < 0; ẽ_1 = (1 - 1/2)/3 = 1/6)
  - v = 2, v' = 3: G = 3 (odd, ẽ_3 = 7/24)
  - v = 2, v' = 5: G = 5 (odd, ẽ_5 = (1 - 1/32)/3 = 31/96)
  - v = 4, v' = 1: G = -1 (negative shift! ẽ_{-1} = (1 - 2)/3 = -1/3)
  - v = 4, v' = 3: G = 1
  - v = 4, v' = 5: G = 3
  - v = 6, v' = 1: G = -3
  - v = 6, v' = 3: G = -1
  - etc.

**G ranges over all ODD integers (positive and negative)** in the g=2, (+,-) recursion. Negative G means s·2^{-G} = s·2^{|G|}, a shifted argument.

But the cross_freq probe defined V_M with g ∈ {0, 2, 4, ...} — only EVEN non-negative. So the (g=2, (+,-)) recursion generates moments with ODD shifts G (and possibly negative), entirely outside V_M.

**This is a direct structural mismatch:** V_M as defined doesn't contain the moments the recursion produces.

For g = 4 (a=-, b=+): v ∈ V_- (odd), v' ∈ V_+ (even). G = v' + 4 - v. v_min = 1, v'_min = 2: G ∈ {2+4-1=5, 2+4-3=3, 2+4-5=1, ..., 4+4-1=7, ...}. All odd.

For g = 6 (mod-9 refinement needed, omitted): expect similar structure.

**Generalization:** for g ∈ {0, 2, 4, ...}, the surviving recursion produces moments with shifts G of OPPOSITE parity to g. So even-g moments produce odd-G moments under iteration, and vice versa.

So V_M = span{M^{ab}(g) : g even} would need to be enlarged to include M^{ab}(G) for G ODD as well — and possibly negative.

---

## 7. The negative-G symmetry

By the Hermitian symmetry of bilinear moments (P_n^{ab}(c) = P_n^{ba}(c)^* — and more generally the Hermitian symmetry of M_n^{ab}(g, c)):

  M_n^{ab}(-G, c) = M_n^{ba}(G, c̃')^* (with some c̃' related to c via a relabeling).

So negative-G moments reduce to positive-G with (a, b) swapped and complex conjugated. They don't add new dimensions if we already include all (a, b) and complex conjugates.

But the recursion produces odd-G moments, which are NOT in V_M = span{even-G}. So V_M doesn't close.

---

## 8. Phase 1 verdict

Three structural findings that make V_M as the brief proposes inadequate:

(F1) **The Tao recursion on M_{n+1}^{ab}(g, c) generates phase-twisted moments M̃_n^{a'b'}(G, c̃; ẽ_G + θ_{v,g}) where θ_{v,g} = 2^v·ẽ_g/3 mod 3^n.** For generic (v, g), θ_{v,g} ≠ 0 and θ_{v,g} ≠ ẽ_{G''} - ẽ_G for any G''. The phase-twisted moments are outside V_M.

(F2) **For each surviving (v, v') in the recursion, the new shift index is G = v' + g - v, which can have either parity.** For g = 2 with surviving (a=+, b=-) parity pairs, G is always ODD. So the recursion generates odd-shift moments outside V_M = span{even-g}.

(F3) **The closure space (if V_M is enlarged to include all phases and shifts) has dimension growing with n** (specifically 4 · 2 · O(N_g) · 3^{n-1} or larger), not a fixed dimension at any g_max. So a "finite truncation T_V at g_max" with level-uniform spectrum does not exist as the brief envisions.

**This is a structural obstruction.** The cross_freq probe's V_M definition closes only for the structural-collapse 2D piece (P_+, P_−)-projection of the **class-summed** moment X̄_n(g, c) = Σ_{a,b} M_n^{ab}(g, c). The CLASS-RESOLVED recursion on individual M_n^{ab}(g, c) does NOT close on V_M.

---

## 9. What might still be tractable

Two routes remain open, both BELOW the brief's "finite truncation T_V → λ_max → 1/2" framing:

**Route α: Class-summed X̄_n(g, c) recursion.** The (1, 4) projection of Off_lin onto (P_+, P_−) uses only X̄_n(c; g) = M_n^{++}(g, c) + M_n^{+-}(g, c) + M_n^{-+}(g, c) + M_n^{--}(g, c). If summing over (a, b) cancels the class-dependent obstruction (it might, since the parity-pattern survival in §5 is class-dependent), then X̄_n(g, c) might satisfy a cleaner recursion on the {g even} family. Worth checking explicitly. This would give a 2-or-3-dim operator on X̄(g) for g ∈ {0, 2, 4} — potentially with eigenvalue 1/2.

**Route β: Truncate at g = 0 only (recovers T_diag) + treat all g ≥ 2 as a perturbation.** This is what T_lead_2x2.py essentially does numerically. The "perturbation Off_n" is summarized by a single rate (empirically 1/2). Rigorous over Q only if the perturbation reduces to a known closed family — but per F1–F3, it doesn't.

Neither route yields the "T_V spectrum at g_max = 2, 4, 6" deliverable the brief asks for. Both would require fresh derivation that goes substantially beyond what cross_freq provided.

---

## 10. Relation to R77.6's branch-cut interpretation

R77.6's finding (E(z)'s singularity at z = 2 is a branch cut, not a simple pole) is CONSISTENT with the operator-level obstruction here. A finite-rank operator over Q on V_M^{(g_max)} with a discrete eigenvalue at 1/2 would be a simple-pole signature; the absence of such an operator (per Phase 1's H_M_RECURSION_UNDERSPECIFIED) is consistent with R77.6's branch-cut picture: the rate-1/2 behavior arises from a continuous spectral structure of the asymptotic operator, not from a finite-dim discrete eigenvalue.

The brief acknowledged this in (A4): "the probe's job is to test whether finite-truncation operators have discrete eigenvalues approaching 1/2, NOT to insist on a discrete eigenvalue at 1/2." Phase 1's finding is: such finite-truncation operators on V_M^{(g_max)} do not exist as the brief envisions — the natural closure is either much larger (V'_M with phase parameters), or operates only on the class-summed projection (Route α), or is not finite-dimensional at all (consistent with R77.6).

---

## 11. Disposition

**H_M_RECURSION_UNDERSPECIFIED.**

The recursion M(g) → M(g') cannot be derived in the form the brief envisions (closure on V_M^{(g_max)} as a finite Q-matrix at fixed g_max). The obstruction is the phase offset θ_{v,g} produced by the unit-shuffle step in the recursion, plus the parity-flipped shift G = v' + g - v that takes even-g moments to odd-G moments. Each is a substantive structural finding (not a hand-wave).

Phase 2–6 are NOT executable on V_M^{(g_max)} as defined by cross_freq materials. The probe halts at Phase 1.

---

## Adversarial check (A1, A2 in scope; A3–A6 not applicable since Phase 2+ blocked)

**(A1) Phase 1 fidelity.** Derivation traces to:
- CROSS_FREQ_PHASE1_EXPANSION.md §0–§5 (R77 sketch §5 + Tao recursion verbatim).
- LTE_3 for v_3(2^g - 1) (standard 3-adic).
- No invocation outside cross_freq + R77 + standard 3-adic number theory.

The new finding (phase offset θ_{v,g}) was NOT identified in cross_freq materials; it emerges from carrying the same level-(n+1) → level-(n) substitution one step further, applied to M_{n+1}(g, c) instead of to P_{n+1}^{ab}(c). The fact that this offset breaks V_M closure is the structural content of Phase 1.

**(A2) Exact rationals.** All quantities (ẽ_g, θ_{v,g}, D_{v,v',g}, ẽ_G) are over Q[1/2] (3-adic integers when relevant). The Phase 1 conclusion is structural (presence/absence of a closure), not numerical.

---

## Files referenced

- `CROSS_FREQ_PHASE1_EXPANSION.md` — sole-frequency Phase 1 (cross-freq at level n+1 → n with v ≠ v')
- `CROSS_FREQ_PHASE1_SPAN.md` — V_M definition
- `CROSS_FREQ_HIGHER_PAIRS.md` — Phase 3a sketch with heuristic cascade comment
- `cross_freq_compute.py` — verification at n = 2, 3 of cross_freq Phase 1 (confirms M(g≥2) ∉ span{P})
- `T_N_OFF_LIN_SPEC.md` — prior H_OFF_LIN_UNDERSPECIFIED disposition
- `result_77_T_lead_spectrum.md` §3, §6 — project's own open-ledger
- `result_77_6_generating_function.md` — branch-cut interpretation (consistent with Phase 1 finding)
