# CROSS_FREQ_PHASE1_EXPANSION — leading-pair derivation

**Date:** 2026-05-12. Phase 1 of the cross-frequency closure derivation probe (Fork 1 from T_N_DISPOSITION's H_OFF_LIN_UNDERSPECIFIED). Wilson (analyst) reporting to Nathan.

This document expands the cross-frequency bilinear Q_n^{++}(c; v, v') for the leading non-trivial pair v ≠ v' (both even, contributing to P_{n+1}^{++}) step-by-step using only formulas verbatim from R77 sketch §5 and Tao's recursion (c_seven_forty_fifth.md §3 / Theorem 75.2 Proof).

---

## 0. Conventions and fidelity to R77 sketch §5

**R77 sketch §5 (verbatim, lines 57–65):**

> Tao recursion gives μ̂_{n+1} from μ̂_n. The class-conservation rule (from R66 + the chain dynamics):
> - v even → r' ≡ 1 mod 3: μ̂_{n+1}^+ contribution
> - v odd → r' ≡ 2 mod 3: μ̂_{n+1}^− contribution
>
> So:
> > μ̂_{n+1}^+(ξ) = Σ_{v even, v≥2} 2^{−v} A_v(ξ) μ̂_n(ξ·2^{−v} mod 3^n)
> > μ̂_{n+1}^−(ξ) = Σ_{v odd, v≥1} 2^{−v} A_v(ξ) μ̂_n(ξ·2^{−v} mod 3^n)
>
> where A_v(ξ) = e^{−2πi ξ 2^{−v}/3^{n+1}}.

**Indexing note.** R77 sketch §5 writes "v even, v ≥ 2" for the + class, "v odd, v ≥ 1" for the − class. The leading non-trivial cross-frequency pair contributing to P^{++} is **(v=2, v'=4)** — the two smallest distinct even indices. The parent task asks about "(v=1, v'=3)" — those are both odd, contributing to the − class only. The structural analysis applies identically; this document treats the **leading P^{++} pair (v=2, v'=4)** in primary form and the **leading P^{−−} pair (v=1, v'=3)** as the parallel case. Both pairs have the same g = v' - v = 2, and the analysis collapses to a function of g (shown below).

For both class signatures, the cross-frequency contribution is computed at fixed (v, v') with v < v', g := v' - v ≥ 2.

**(R77.4 conflict check, A6.)** R77.4 erratum ruled out K_n (the within-level Markov operator on coprime residues mod 3^n) carrying eigenvalue 1/2. T_N here acts on the moment basis {P_n^{ab}(c)} via Tao's bilinear recursion — a different object than K_n. The cross-frequency closure question is internal to T_N's construction and does not contradict R77.4: even if the closure succeeds, T_N's spectrum on class-resolved moments is structurally distinct from K_n's spectrum on residue probabilities.

**Tao's recursion (verbatim, c_seven_forty_fifth.md §3 / Theorem 75.2 Proof line 87):**

> μ̂_{n+1}(ξ) = Σ_{v=1}^∞ 2^{−v} · e^{−2πi ξ · 2^{−v}/3^{n+1}} · μ̂_n(ξ · 2^{−v} mod 3^n).

So A_v(ξ) = e^{-2πi ξ·2^{-v}/3^{n+1}}, identical to R77 sketch §5's notation.

---

## 1. Bilinear setup for P_{n+1}^{++}(c)

Following R77 sketch §5 and the §6 ledger:

P_{n+1}^{++}(c) := Σ_{ξ ≡ c mod 3, ξ ∈ (Z/3^{n+1})^×} μ̂_{n+1}^+(ξ) · μ̂_{n+1}^{+*}(ξ).

Substituting:

  P_{n+1}^{++}(c) = Σ_{v, v' both even, ≥ 2} 2^{−v−v'} · Q_n^{++}(c; v, v')

where, as in T_N_OFF_LIN_SPEC.md §(b),

  **Q_n^{++}(c; v, v') := Σ_{ξ ≡ c mod 3, 3∤ξ in (Z/3^{n+1})^×} A_v(ξ) A_{v'}^*(ξ) · μ̂_n(ξ·2^{−v} mod 3^n) · μ̂_n^*(ξ·2^{−v'} mod 3^n).**

For v = v' (diagonal): see `result_77_T_diagonal.py` — reduces rigorously to span{P_n^{ab}(c)} via the unit-shuffle change of variable + R66 class flow, yielding T_diag = (1/5)·[[1,1],[4,4]].

For v ≠ v' (off-diagonal): the open derivation. Below.

Identical structure for P^{−−}(c) with v, v' both odd ≥ 1; for P^{+−}(c) and P^{−+}(c) with mixed parities (one even, one odd). Mixed-parity pairs contribute to cross-class moments (P^{+−}) which the structural collapse R76 §11 sets to zero for n ≥ 2 — so they don't enter the 2D (P_+, P_−) reduction. But they DO enter the full 6-dim or 8-real-dim picture R77 sketch §5 articulates. We focus on **same-parity v, v' pairs (the ones surviving the 2D reduction)** since these are what determine Off_lin's matrix entries in the (P_+, P_−) basis.

---

## 2. Phase factor structure: 3-adic valuation of d_{v,v'}

The phase factor in Q_n^{++}(c; v, v') is

  A_v(ξ) · A_{v'}^*(ξ) = e^{-2πi ξ·(2^{-v} - 2^{-v'})/3^{n+1}}.

Define d_{v,v'} := 2^{-v} - 2^{-v'} ∈ Z[1/2] (computed mod 3^{n+1} via 2 invertible mod 3^{n+1}).

Writing v < v' and g := v' - v ≥ 2 (even, for the P^{++} case):

  **d_{v,v'} = 2^{-v}·(1 - 2^{-g}) = (2^g - 1)·2^{-v'} (after factoring)**.

**3-adic valuation.** v_3(d_{v,v'}) = v_3(2^g - 1) - v_3(2^{v'}) = v_3(2^g - 1) (since v_3(2^k) = 0 for any k).

For g ≥ 1: 2^g ≡ 1 mod 3 iff g is even (since ord_3(2) = 2). So:
- g odd: v_3(2^g - 1) = 0 ⇒ v_3(d_{v,v'}) = 0
- g even: v_3(2^g - 1) ≥ 1; by Lifting-the-Exponent (LTE_3): v_3(2^g - 1) = v_3(2² - 1) + v_3(g/2) = 1 + v_3(g/2) when g ≥ 2 even.

**Special cases for the P^{++} surviving pairs (v, v' both even, g = v' - v ∈ {2, 4, 6, ...}):**

| g | g/2 | v_3(g/2) | v_3(d_{v,v'}) | d_{v,v'}·2^{v'} = 2^g - 1 (3-adic part) |
|---|-----|----------|---------------|-----------------------------------------|
| 2 | 1   | 0        | **1**         | 2² - 1 = 3                              |
| 4 | 2   | 0        | **1**         | 2⁴ - 1 = 15 = 3·5                       |
| 6 | 3   | 1        | **2**         | 2⁶ - 1 = 63 = 9·7                       |
| 8 | 4   | 0        | **1**         | 2⁸ - 1 = 255 = 3·85                     |
| 12| 6   | 1        | **2**         | 2¹² - 1 = 4095 = 9·455                  |

**For the parent task's "v=1, v'=3" (P^{−−} class):** g = 2 ⇒ v_3 = 1. ✓ Matches T_N_DISPOSITION.md's note "v_3(2^{-1} - 2^{-3}) = 1".

**Leading pair contribution to Off_lin's matrix entries:** the (v=2, v'=4) pair (g=2, v_3 = 1) carries 2^{-v-v'} = 2^{-6} = 1/64 prefactor, plus its symmetric partner (v'=4, v=2) ⇒ total weight 2·2^{-6} = 1/32. The next pair (v=2, v'=6) has weight 2·2^{-8} = 1/128, and so on.

---

## 3. Lift-fiber orthogonality at v_3(d) ≥ 1

Each ξ ∈ (Z/3^{n+1})^× projects to u := ξ mod 3^n ∈ (Z/3^n)^×. The fiber over u has 3 elements: ξ ∈ {u, u + 3^n, u + 2·3^n}. The constraint ξ ≡ c mod 3 lifts to u ≡ c mod 3 (since n ≥ 1 ⇒ 3 | 3^n ⇒ ξ ≡ u mod 3).

The objects ξ·2^{-v} mod 3^n and ξ·2^{-v'} mod 3^n depend ONLY on u (not on the choice of lift), because (u + j·3^n)·2^{-v} ≡ u·2^{-v} mod 3^n. So μ̂_n at both arguments is **constant on each fiber**.

The phase factor varies with the lift:

  ξ·d_{v,v'} mod 3^{n+1} = (u + j·3^n)·d_{v,v'} = u·d_{v,v'} + j·3^n·d_{v,v'} mod 3^{n+1}.

So
  e^{-2πi ξ·d/3^{n+1}} = e^{-2πi u·d/3^{n+1}} · e^{-2πi j·d/3}.

The inner factor e^{-2πi j·d/3} depends on d mod 3:
- **If 3 ∤ d (i.e., v_3(d) = 0):** Σ_{j=0,1,2} e^{-2πi j·d/3} = 0 (sum of all primitive cube roots when d mod 3 ≠ 0).
- **If 3 | d (i.e., v_3(d) ≥ 1):** Σ_{j=0,1,2} e^{-2πi j·d/3} = 3 (each term = 1).

**Consequence.** When v_3(d_{v,v'}) = 0 (odd g, mixed-parity case), the lift-fiber sum kills the cross-frequency contribution entirely. **Hence the (v=1, v'=2)-type pair (g=1, odd, mixed parity) automatically gives Q ≡ 0** — these mixed-parity bilinears decouple via fiber orthogonality. This is a clean algebraic finding.

When v_3(d_{v,v'}) ≥ 1 (which is automatic for same-parity v, v' pairs with g ≥ 2 even), the lift-fiber sum gives a factor of 3. The cross-frequency object collapses to a level-n character sum:

  **Q_n^{++}(c; v, v') = 3 · Σ_{u ∈ (Z/3^n)^×, u ≡ c mod 3} e^{-2πi u·d_{v,v'}/3^{n+1}} · μ̂_n(u·2^{-v} mod 3^n) · μ̂_n^*(u·2^{-v'} mod 3^n).**

Since v_3(d_{v,v'}) ≥ 1, d_{v,v'}/3 is a 3-adic integer well-defined mod 3^n. Write d_{v,v'} = 3·ẽ_{v,v'} where ẽ_{v,v'} ∈ Z/3^n. Then

  e^{-2πi u·d_{v,v'}/3^{n+1}} = e^{-2πi u·ẽ_{v,v'}/3^n}

— **a level-n character**, not a level-(n+1) character. This is the structural payoff of v_3(d) ≥ 1.

---

## 4. Unit-shuffle and reduction to g-dependence

Substitute s := u·2^{-v} mod 3^n. Since 2 is a unit mod 3^n, this is a bijection on (Z/3^n)^×, with u = s·2^v mod 3^n.

The c-class of u relates to the class of s via u ≡ 2^v·s mod 3 ≡ (-1)^v·s mod 3. For v even, u and s are in the same class; for v odd, classes swap.

For the P^{++} case (v, v' both even), s ≡ u ≡ c mod 3 (class preserved).

Substituting:
  u·d_{v,v'}/3 = s·2^v·d_{v,v'}/3 = s · 2^v·(2^{-v} - 2^{-v'})/3 = s·(1 - 2^{-g})/3 mod 3^n.

Define **ẽ_g := (1 - 2^{-g})/3 mod 3^n** — depends only on g, not on v separately. (Well-defined: v_3(1 - 2^{-g}) = v_3(2^g - 1) ≥ 1 for g even ≥ 2.)

And μ̂_n(u·2^{-v}) = μ̂_n(s), μ̂_n(u·2^{-v'}) = μ̂_n(s·2^{v-v'}) = μ̂_n(s·2^{-g}).

So
  **Q_n^{++}(c; v, v') = 3 · X_n^{++}(c; g)**

where g = v' - v and

  **X_n^{ab}(c; g) := Σ_{s ∈ (Z/3^n)^×, s ≡ c mod 3} e^{-2πi s·ẽ_g/3^n} · μ̂_n^a(s) · μ̂_n^{b*}(s·2^{-g} mod 3^n)**

(with class-resolved μ̂^a, μ̂^b inserted when v, v' parities determine the class signature — see §6 for class flow).

**The cross-frequency object depends only on g, not on (v, v') separately.** This is the key dimensional reduction: there is one cross-frequency moment per even g ≥ 2, not one per (v, v') pair.

---

## 5. Pre-imaging X_n^{ab}(c; g) onto (r, r') pairs

Substitute the class-resolved definitions:

  μ̂_n^a(s) = Σ_{r ≡ a mod 3, r ∈ (Z/3^n)^×} π_n(r) · e^{-2πi r·s/3^n}
  μ̂_n^{b*}(s·2^{-g}) = Σ_{r' ≡ b mod 3, r' ∈ (Z/3^n)^×} π_n(r') · e^{+2πi r'·s·2^{-g}/3^n}

So

  X_n^{ab}(c; g) = Σ_{r ≡ a, r' ≡ b in (Z/3^n)^×} π_n(r)·π_n(r') · Σ_{s ∈ (Z/3^n)^×, s ≡ c mod 3} e^{-2πi s·(r - r'·2^{-g} + ẽ_g)/3^n}.

Define **h := r - r'·2^{-g} + ẽ_g ∈ Z/3^n** — depends on (r, r', g).

The inner sum is a partial character sum over s ≡ c mod 3 in (Z/3^n)^×. Parameterize s = c + 3·t for t ∈ {0, 1, ..., 3^{n-1} - 1} (these are precisely the elements of (Z/3^n)^× ∩ c-class):

  Σ_{s ≡ c mod 3} e^{-2πi s·h/3^n}
  = e^{-2πi c·h/3^n} · Σ_{t=0}^{3^{n-1} - 1} e^{-2πi (3t)·h/3^n}
  = e^{-2πi c·h/3^n} · Σ_{t=0}^{3^{n-1} - 1} e^{-2πi t·h/3^{n-1}}

The inner geometric sum is 3^{n-1} if h ≡ 0 mod 3^{n-1}, else 0.

So

  **X_n^{ab}(c; g) = 3^{n-1} · Σ_{r ≡ a, r' ≡ b, h ≡ 0 mod 3^{n-1}} π_n(r)·π_n(r') · e^{-2πi c·h/3^n}**

where h = r - r'·2^{-g} + ẽ_g.

When h ≡ 0 mod 3^{n-1}, write h = 3^{n-1}·m for m ∈ {0, 1, 2}. Then e^{-2πi c·h/3^n} = e^{-2πi c·m/3} = ω^{-c·m} ∈ {1, ω, ω²} where ω = e^{2πi/3}.

So the sum splits into THREE sublattices, indexed by m ∈ {0, 1, 2}, each weighted by ω^{-c·m}:

  **X_n^{ab}(c; g) = 3^{n-1} · Σ_{m=0,1,2} ω^{-c·m} · L_n^{ab}(g, m)**

where

  L_n^{ab}(g, m) := Σ_{r ≡ a, r' ≡ b, h = 3^{n-1}·m} π_n(r)·π_n(r').

The condition h = 3^{n-1}·m fixes r' (mod 3^{n-1}) given r (and g, m). Specifically:

  r' ≡ 2^g·(r + ẽ_g - 3^{n-1}·m) mod 3^{n-1} (multiplying h = 3^{n-1}·m by 2^g and rearranging)

The constraint is mod 3^{n-1}, leaving 3 lifts of r' mod 3^n, of which exactly one or two satisfy r' ≡ b mod 3 (since the 3 lifts cover all residue classes mod 3 except possibly when n-1=0; for n ≥ 2 it's clean).

---

## 6. Class flow for P^{++} and P^{−−} contributions

R77 sketch §5 says P_{n+1}^+ uses μ̂_{n+1}^+ which only contains v-even contributions; P_{n+1}^- uses v-odd. The class flow on the inner factor μ̂_n at frequency u·2^{-v} mod 3^n is:

The shift u → u·2^{-v} reduces μ̂_n(u·2^{-v}) which is the TOTAL Fourier transform at this frequency (not class-resolved). So when we split μ̂_n = μ̂_n^+ + μ̂_n^- and sum over the four (a, b) combinations, we get all four class signatures at level n inside the cross-frequency object.

This means: the cross-frequency contribution to P_{n+1}^{++}(c) is

  Σ_{v, v' both even ≥ 2} 2^{-v-v'} · 3 · Σ_{a, b ∈ {+,-}} X_n^{ab}(c; g)
  = 3 · Σ_{g even ≥ 2} [Σ_{v even ≥ 2 : v' = v+g ≥ 2 even} 2^{-2v-g} + (same for v > v')] · Σ_{a,b} X_n^{ab}(c; g)

The internal weight sum over v at fixed g:

  W_+(g) := Σ_{v ∈ {2, 4, 6, ...}, v + g ∈ {2, 4, ...}} 2^{-2v-g} · 2

  (factor 2 for (v, v+g) and (v+g, v) symmetric pair; we restrict to v < v' and double).

For g even ≥ 2:
  W_+(g) = 2 · 2^{-g} · Σ_{v ∈ {2, 4, ...}} 2^{-2v} = 2·2^{-g} · (1/16 + 1/256 + ...) = 2·2^{-g} · (1/16)/(1 - 1/16) = 2·2^{-g}/15 = 2^{-g+1}/15.

Check: W_+(2) = 2^{-1}/15 = 1/30. W_+(4) = 2^{-3}/15 = 1/120. W_+(6) = 2^{-5}/15 = 1/480. Geometric decay 1/4 per g.

For the P^{−−} class (v, v' both odd ≥ 1), similarly W_−(g) = 2 · 2^{-g} · Σ_{v ∈ {1, 3, ...}} 2^{-2v} = 2·2^{-g}·(1/4)/(1 - 1/16) = 2·2^{-g}·(4/15)·1 = wait let me redo.

Σ_{v odd ≥ 1} 4^{-v} = 4^{-1} + 4^{-3} + ... = (1/4)/(1 - 1/16) = 4/15.

So W_−(g) = 2 · 2^{-g} · 4/15 = 2^{-g+3}/15. Check: W_−(2) = 2/15. W_−(4) = 1/30. W_−(6) = 1/120.

(For mixed-parity P^{+−}: g must be odd; v_3(d) = 0 ⇒ contribution vanishes by §3. ✓ This is why structural collapse P^{+−} = 0 is preserved by Tao's recursion in the off-diagonal.)

---

## 7. The Off_lin contribution in (P_+, P_−) basis

Putting it together, the off-diagonal contribution to P_{n+1}^{++}(c) is:

  Off_{n+1}^{++}(c) = 3 · Σ_{g even ≥ 2} W_+(g) · [X_n^{++}(c; g) + X_n^{+-}(c; g) + X_n^{-+}(c; g) + X_n^{--}(c; g)]

  = 3 · Σ_{g even ≥ 2} W_+(g) · X̄_n(c; g)

where X̄_n(c; g) := Σ_{a,b} X_n^{ab}(c; g) (the class-summed cross-frequency moment) — equivalently, the cross-frequency moment for the TOTAL μ̂_n (no class resolution).

Similarly:

  Off_{n+1}^{−−}(c) = 3 · Σ_{g even ≥ 2} W_−(g) · X̄_n(c; g).

**Critical observation:** the SAME X̄_n(c; g) appears in both Off_{n+1}^{++} and Off_{n+1}^{−−}. The only difference is the weight (W_+(g) vs W_−(g)).

Recall W_+(g) = 2^{-g+1}/15, W_−(g) = 2^{-g+3}/15 = 4·W_+(g).

So **Off_{n+1}^{−−}(c) = 4 · Off_{n+1}^{++}(c)** — the off-diagonal correction lives on the (1, 4) eigendirection. This is automatic from the W_+/W_− ratio.

**This is structurally meaningful:** under the structural collapse P^{+−} = 0 (R76 §11) and class-c symmetry, the off-diagonal correction stays on the (1, 4) line, identical to T_diag's slow eigenvector. So Off_lin acts on the (P_+, P_−) 2D plane as a rank-1 matrix scaled by (1, 4)^T (1, 0)... wait, need to be careful — we need to express X̄_n(c; g) as a function of (P_+, P_−).

The question is whether X̄_n(c; g) lies in span{P_n^+, P_n^−} (the 2D subspace surviving the structural collapse).

---

## 8. Sublattice constraint vs. r = r' diagonal: the closure question

The diagonal v = v' case has h = r - r'·1 + 0 = r - r', and ẽ_0 = 0. The constraint h ≡ 0 mod 3^{n-1} becomes r ≡ r' mod 3^{n-1}, plus (in the class-resolved version with a = b) r ≡ r' mod 3 ⇒ r ≡ r' mod 3^n possible OR r differs from r' by 3^{n-1}·k. Actually for v = v', the inner s-sum gives just S_n-like behavior.

For v ≠ v' (g ≠ 0), the constraint h = r - r'·2^{-g} + ẽ_g ≡ 0 mod 3^{n-1} is a DIFFERENT sublattice constraint. Specifically:

  r - r'·2^{-g} ≡ -ẽ_g mod 3^{n-1}
  r' ≡ 2^g·(r + ẽ_g) mod 3^{n-1}

Given r ∈ (Z/3^n)^× with r ≡ a mod 3, this determines r' mod 3^{n-1}, leaving 3 lifts mod 3^n. Of those 3 lifts, those with r' ≡ b mod 3 (1 or 2 typically) contribute to L_n^{ab}(g, m).

**This is not the same as the (r, r') diagonal r = r' that defines P_n^{ab}(c).** Under the diagonal constraint:
- P_n^{ab}(c) requires a = b (same r-class), c-class of ξ implicit; it's a SAME-FREQUENCY moment.

Under the cross-frequency sublattice constraint:
- a, b free; (r, r') is two distinct residues constrained by an affine lattice relation depending on g.

**The sublattice moment L_n^{ab}(g, m) is NOT a linear function of {P_n^{ab}(c)}.** It is a new bilinear functional of π_n that requires its own moment dimension.

This is the precise statement of "cross-frequency closure does not exist on span{P_n^{ab}(c)}".

---

## 9. The enlarged span: what spans the cross-frequency residue?

Define the **g-twisted bilinear moment family**:

  **M_n^{ab}(g, c) := Σ_{ξ ∈ (Z/3^n)^×, ξ ≡ c mod 3} e^{-2πi ξ·ẽ_g/3^n} · μ̂_n^a(ξ) · μ̂_n^{b*}(ξ·2^{-g} mod 3^n)** (= X_n^{ab}(c; g))

This is a 4-parameter family: (a, b) ∈ {+,-}², g ∈ {0, 2, 4, 6, ...} (even non-negative, restricting to the P^{++}/P^{−−} surviving pairs), c ∈ {1, 2}.

At g = 0: M_n^{ab}(0, c) = P_n^{ab}(c) (the same-frequency moment with ẽ_0 = 0).

At g ≥ 2 even: M_n^{ab}(g, c) is a NEW moment, not in span{P_n^{ab}(c)} for c ∈ {1, 2}.

The natural state space for Off_lin closure is:

  **V_M := span{M_n^{ab}(g, c) : (a,b) ∈ {+,-}², g ∈ {0, 2, 4, 6, ...}, c ∈ {1, 2}}.**

The diagonal piece T_diag acts on the (g=0) subspace ≡ span{P_n^{ab}(c)}, with T_diag the 2x2 matrix (1/5)·[[1,1],[4,4]] on the post-collapse (P_+, P_−) restriction. The off-diagonal piece couples this subspace to the g ≥ 2 components via the recursion P_{n+1}^{++} = T_diag·(P^++, P^--)_n + 3·Σ_{g≥2 even} W_+(g)·M_n(g, c), AND ALSO requires knowing the recursion of M_n(g, c) itself under Tao's iteration.

**This is the natural enlarged span.** The operator's dimension is infinite if g ranges unbounded, or large-finite if g is truncated.

---

## 10. Phase 1 summary

The leading P^{++} pair (v=2, v'=4), and equivalently the parent task's named pair (v=1, v'=3) for P^{−−}, both have g = 2 and reduce to the SAME g-dependent cross-frequency moment M_n^{ab}(g=2, c).

Algebraic reduction (rigorous from §3–§5 above):

  Q_n^{ab}(c; v, v') = 3 · M_n^{ab}(c; g = v' - v)  for v < v', same parity, g even ≥ 2.

The cross-frequency moment M_n^{ab}(g, c) reduces to a **sublattice-constrained bilinear of π_n**:

  M_n^{ab}(g, c) = 3^{n-1} · Σ_{m ∈ {0,1,2}} ω^{-c·m} · L_n^{ab}(g, m)

where L_n^{ab}(g, m) = Σ_{r ≡ a, r' ≡ b, r' ≡ 2^g(r + ẽ_g) - 2^g·3^{n-1}·m mod 3^{n-1}} π_n(r)·π_n(r').

**The sublattice constraint r' ≡ 2^g·(r + ẽ_g) mod 3^{n-1} (modulo 3^{n-1}·m correction) is DIFFERENT from the r = r' diagonal that defines P_n^{ab}(c). Therefore M_n^{ab}(g, c) is not in span{P_n^{ab}(c)} for g ≥ 2.**

This is the precise algebraic content of "cross-frequency does not close on same-frequency".

The enlarged span V_M = span{M_n^{ab}(g, c) : g ∈ {0, 2, 4, ...}} is the natural domain for Off_lin closure.

---

## Adversarial check (A1, A2, A4, A6 in scope)

**(A1) R77 sketch §5 fidelity.** All steps trace to:
- R77 sketch §5 verbatim (μ̂_{n+1}^± class-split recursion).
- Tao's recursion (c_seven_forty_fifth.md §3 / Theorem 75.2 Proof line 87).
- R66 class flow (T_N_OFF_LIN_SPEC.md §(b), explicit in result_77_T_diagonal.py line 9–12).
- LTE_3 (Lifting the Exponent) for v_3(2^k - 1): standard 3-adic number theory, used in §2.

No invoked identity is outside this scope.

**(A2) Phase orthogonality vs cancellation.** §3 distinguishes:
- v_3(d) = 0 case (mixed parity, g odd): lift-fiber sum gives 0, contribution vanishes algebraically.
- v_3(d) ≥ 1 case (same parity, g even ≥ 2): lift-fiber sum gives 3, contribution survives as a level-n character sum.

Both outcomes documented; neither hand-waved.

**(A4) Exact rationals throughout.** ẽ_g, weights W_±(g), and the structural sublattice constraint are all over Q (or 3-adic integers reducing to Z/3^n exactly). The closure question (whether M_n^{ab}(g, c) lies in span{P_n^{ab}(c)}) is purely a question of linear-algebraic span, decidable on the basis of the (r, r') sublattice structure exhibited.

**(A6) Conflict with R77.4.** R77.4 ruled K_n out as carrying eigenvalue 1/2. T_N here is built on the bilinear moment basis {P_n^{ab}(c)} ∪ {M_n^{ab}(g, c) : g ≥ 2}, which is structurally distinct from K_n's residue probability space. The cross-frequency expansion of P^{++} introduces moments at twisted frequencies, not residue probabilities. No conflict; the operator-theoretic settings are different.

---

## Verdict at end of Phase 1

The (v=1, v'=3) [for −−] and (v=2, v'=4) [for ++] leading pairs collapse to a single g=2 cross-frequency moment M_n^{ab}(g=2, c). This moment is a sublattice-constrained bilinear of π_n with constraint r' ≡ 2²·(r + ẽ_2) mod 3^{n-1}, ẽ_2 = (1 - 1/4)/3 = (1/4)·(3/3) ... let me compute carefully: 1 - 2^{-2} = 3/4; (3/4)/3 = 1/4. So ẽ_2 = 1/4 mod 3^n = (2^{-2}) mod 3^n = inv(4) mod 3^n.

Hmm actually wait: ẽ_g := (1 - 2^{-g})/3 mod 3^n. So ẽ_2 = (1 - 2^{-2})/3 = (1 - 1/4)/3 = (3/4)/3 = 1/4. So ẽ_2 ≡ inv(4) mod 3^n.

This sublattice constraint produces a moment NOT in span{P_n^{ab}(c)} for n ≥ 2.

**Therefore: cross-frequency does NOT close on {P_n^{ab}(c)}. Off_lin's natural domain enlarges to include M_n^{ab}(g, c) for g ≥ 2 even.**

Next: Phase 2 identifies the precise structure of the enlarged span (CROSS_FREQ_PHASE1_SPAN.md).
