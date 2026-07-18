# SL2_EMBEDDING_CANDIDATES — three candidate SL_2(ℝ) embeddings of Syracuse dynamics

**Date:** 2026-05-12. Phase 2 of the SL_2(ℝ)-embedding structural-compatibility probe.

This file enumerates candidate SL_2(ℝ) random-walk constructions that could plausibly place
Syracuse dynamics inside the Furstenberg-measure framework, with T_lead = (1/45)·[[7,9],[28,36]]
emerging as a derived rank-1 projection. For each candidate, the four (or five) framework gates
(G1-G5) from SL2_FRAMEWORK_HYPOTHESES are checked.

T_lead's load-bearing properties:
- T_lead = (1/45)·[[7,9],[28,36]]
- **det(T_lead) = (7·36 − 9·28)/45² = (252 − 252)/2025 = 0**
- rank(T_lead) = 1
- Spectrum {43/45, 0}, eigenvectors (1, 4) at 43/45 and (9, −7) at 0
- Outer-product form: T_lead = (1/45) · (1, 4)^T · (7, 9) (column · row, both integer)

---

## Candidate A: Direct lift of T_lead's outer-product factorization

### A.1 Construction

T_lead = (1/45) · u v^T where u = (1, 4)^T, v = (7, 9). Two natural lifts:

  **A.1a: Right-action lift.** Lift T_lead to act on the row-vector space; complete u to a 2×2
  matrix with column u and some second column u_⊥ to be chosen, scaled to unit determinant.

  **A.1b: Left-action lift.** Similarly with v as a row of a 2×2 matrix.

  **A.1c: Full completion.** Find a 2×2 matrix M such that:
    - M restricted to the (1, 4) direction acts by 43/45 (matches T_lead's leading eigenvalue),
    - M has determinant 1 (lives in SL_2(ℝ)).

  For (c), parametrize M = [[a, b], [c, d]] with ad - bc = 1, and require M · (1, 4)^T =
  (43/45) · (1, 4)^T. That gives a + 4b = 43/45 and c + 4d = 4·43/45 = 172/45.

  These are TWO linear constraints; the additional det = 1 constraint plus the (9, −7)-null
  constraint of T_lead would over-determine. Drop the null constraint (T_lead's null direction
  is not a structural feature of any natural lift — it's an artifact of T_lead's rank-1 status).

### A.2 Concrete example

A 1-parameter family of det-1 lifts satisfying M(1,4)^T = (43/45)(1,4)^T:

  M_t = [[a, b], [c, d]] with
    a + 4b = 43/45,    c + 4d = 172/45,    ad - bc = 1

Solving: let b = t (free), then a = 43/45 - 4t; let d = s (free), then c = 172/45 - 4s.
det = ad - bc = (43/45 - 4t)·s - t·(172/45 - 4s)
              = (43/45) s - 4ts - (172/45) t + 4st
              = (43/45) s - (172/45) t
              = (43 s - 172 t)/45

Setting det = 1: 43 s - 172 t = 45, i.e., s = (45 + 172 t)/43.

So a 1-parameter family of SL_2(ℝ) matrices completing the (1, 4)-eigenvalue 43/45 exists. Pick
e.g. t = 0:  s = 45/43, b = 0, a = 43/45, c = 172/45, d = 45/43.

  **M_0 = [[43/45, 0], [172/45, 45/43]]**

Check: det = (43/45)(45/43) - 0·(172/45) = 1. ✓
       M_0 (1, 4)^T = (43/45, 172/45 + 180/43)^T.

Wait — recompute: 172/45 + 4·(45/43) = 172/45 + 180/43.
Common denom 45·43 = 1935: 172·43/1935 + 180·45/1935 = (7396 + 8100)/1935 = 15496/1935.
Compare to 4·43/45 = 172/45 = 172·43/1935 = 7396/1935. **Not equal.** So M_0 (1, 4) ≠ (43/45)(1, 4)
— my parametrization was wrong. Re-derive.

  M (1, 4)^T = (a + 4b, c + 4d)^T should equal (43/45)(1, 4)^T = (43/45, 172/45)^T.

So a + 4b = 43/45 (correct) AND c + 4d = 172/45 (correct). Above I used the SAME parametrization
but the row second-coordinate check is c + 4d = 172/45 — let me re-test M_0:
  c + 4d = 172/45 + 4·(45/43) = 172/45 + 180/43

These have to equal 172/45. So 4·(45/43) must equal 0, i.e., d = 0. Then c = 172/45.
But then det = ad - bc = 0 - b·c = -b·c = -0·(172/45) = 0. **Forced det = 0.** The natural choice
d = 0 makes the lift rank-1 again — it's T_lead.

Let me redo: with b = t, d = s, the two eigenvalue constraints are independent of det. From
a + 4b = 43/45: a = 43/45 - 4t. From c + 4d = 172/45: c = 172/45 - 4s.
det = (43/45 - 4t)·s - t·(172/45 - 4s) = (43/45)s - 4ts - (172/45)t + 4ts = (43s - 172t)/45.
Setting det = 1: 43s = 45 + 172t, s = (45 + 172t)/43.

For t = 1: s = 217/43, a = 43/45 - 4 = -137/45, b = 1, c = 172/45 - 868/43 = (172·43 - 868·45)/(45·43)
        = (7396 - 39060)/1935 = -31664/1935, d = 217/43.

  **M_1 = [[-137/45, 1], [-31664/1935, 217/43]]**

Verify det:
  ad = (-137/45)(217/43) = -29729/1935
  bc = 1 · (-31664/1935) = -31664/1935
  ad - bc = -29729/1935 - (-31664/1935) = 1935/1935 = 1. ✓

Verify M_1 (1, 4)^T:
  first coord: -137/45 + 4 = -137/45 + 180/45 = 43/45. ✓
  second coord: -31664/1935 + 4·217/43 = -31664/1935 + 868/43 = -31664/1935 + 868·45/1935
              = -31664/1935 + 39060/1935 = 7396/1935 = 172/45. ✓

So M_1 ∈ SL_2(ℚ) is a det-1 lift of T_lead's leading eigenstructure. Entries are rational
(hence algebraic). **Candidate A passes G1 (det = 1) and G4 (algebraic entries in SL_2(ℚ)).**

### A.3 Gate checks for Candidate A

| Gate | Status | Comment |
|---|---|---|
| G1: det = 1 | ✓ | M_1 above; 1-param family in SL_2(ℚ) |
| G2: non-elementary | **FAILS — single matrix** | A single matrix M_1 generates a cyclic group (powers of M_1), which is contained in a 1-dim algebraic subgroup of SL_2(ℝ). For non-elementary, we need ≥ 2 generators with no common invariant set in P^1. A random walk driven by ONE matrix is degenerate (μ = δ_{M_1}, μ^{*n} = δ_{M_1^n}, ν is supported on the M_1-orbits of P^1 which is either a single point or two points or a single 1-parameter curve). **Single-matrix non-elementary is impossible.** |
| G3: moment | ✓ | Finite-support measure (just δ_{M_1}) has all moments. |
| G4: algebraic | ✓ | Entries in ℚ. |
| G5: Zariski-dense | FAILS | A single matrix generates a 1-dim subgroup. |
| **Verdict** | **G2 fails** | A direct lift via ONE matrix is structurally inadequate. Need a probability measure with at least two matrices to even ask the non-elementary question. |

### A.4 Multi-matrix variant

Could a probability measure μ supported on TWO or more matrices recover T_lead as the "averaged
random-walk operator"? Two natural sources of multiple matrices in Syracuse dynamics:

  - **Parity-branched lift.** Syracuse has two transition branches (v even → class +,
    v odd → class −). Each branch could give a matrix; μ is a 2-atom measure on
    {(M_+ with weight p), (M_- with weight 1-p)}, and T_lead is the **expectation matrix**
    E[M] = p M_+ + (1-p) M_-.

  - **Geometric-step branched lift.** Tao's recursion sums over v ∈ {1, 2, 3, ...} with weights
    2^{-v}; each v contributes a matrix M_v. μ is a geometric-distribution measure
    P(M = M_v) = 2^{-v}, and T_lead is again E[M] = Σ_v 2^{-v} M_v.

Both readings make μ a non-degenerate probability measure on a multi-atom support.
**HOWEVER**: in the Furstenberg framework, the "operator" of interest is NOT E[M] (the
expectation), but the random-walk product g_n ⋯ g_1 of i.i.d. samples. T_lead = E[M] gives the
**first-moment information** about μ, not the Furstenberg measure of μ.

This is a critical distinction:
- E[M] is a linear average; can be rank-1 even when individual M_i have det = 1.
- The Furstenberg measure ν on P^1 is the stationary measure of the SL_2 action; it depends on
  ALL of μ, not just E[M].

**So Candidate A in its "T_lead is the embedding" reading is structurally wrong.** T_lead is at
best the *expectation* of an SL_2(ℝ)-valued random matrix; the Furstenberg framework's input
is the *measure* μ, and the conclusion is about ν on P^1.

### A.5 Reinterpretation: T_lead as expectation of an SL_2 random variable

If T_lead = E[M] for some SL_2(ℝ)-valued random matrix M with distribution μ, the Furstenberg
framework would apply to μ, give ν on P^1 with polynomial Fourier decay, and the relation
between μ_n (Syracuse stationary on ℤ_3) and ν (Furstenberg on P^1) would need to be derived
separately.

Constructing μ such that E[M] = T_lead: choose μ supported on {M_+, M_-} with weights {p, 1-p}
satisfying p M_+ + (1-p) M_- = T_lead = (1/45)·[[7, 9], [28, 36]].

Pick p = 1/2 and parameterize: M_+ - M_- = (something) and (M_+ + M_-)/2 = T_lead. So
M_+ = T_lead + Δ, M_- = T_lead - Δ for any 2×2 matrix Δ.

Require M_+, M_- ∈ SL_2(ℝ):
  det(T_lead + Δ) = 1 AND det(T_lead - Δ) = 1.

  det(T_lead + Δ) = det(T_lead) + (cross terms) + det(Δ) = 0 + tr(adj(T_lead) Δ) + det(Δ).
  adj(T_lead) = [[36/45, -9/45], [-28/45, 7/45]]. So tr(adj(T_lead) Δ) = (36 Δ_11 − 9 Δ_21 − 28 Δ_12 + 7 Δ_22)/45.

  For BOTH det(T_lead ± Δ) = 1: the linear term in Δ must vanish (else one sign gets +, other -),
  and det(Δ) = 1.

  Linear term vanishing: (36 Δ_11 + 7 Δ_22 - 9 Δ_21 - 28 Δ_12) / 45 = 0, i.e.,
  36 Δ_11 + 7 Δ_22 = 9 Δ_21 + 28 Δ_12.

  And det(Δ) = 1.

These two constraints leave a 2-parameter family of Δ. Example: Δ_11 = 0, Δ_22 = 0, then
9 Δ_21 + 28 Δ_12 = 0, so Δ_12 = -(9/28) Δ_21. det(Δ) = -Δ_12 Δ_21 = (9/28) Δ_21² = 1, so
Δ_21 = √(28/9) = 2√7/3 (irrational!). With Δ_21 = 2√7/3, Δ_12 = -(9/28)(2√7/3) = -3√7/14.

So **M_+ = (1/45)·[[7, 9], [28, 36]] + [[0, -3√7/14], [2√7/3, 0]]** and
**M_- = (1/45)·[[7, 9], [28, 36]] - [[0, -3√7/14], [2√7/3, 0]]**.

Entries now contain √7 — algebraic of degree 2 over ℚ. Both M_+, M_- ∈ SL_2(K) for K = ℚ(√7).

### A.6 Gate checks for Candidate A (reinterpreted, equal-weight 2-atom)

| Gate | Status | Comment |
|---|---|---|
| G1: det = 1 for both atoms | ✓ | By construction. |
| G2: non-elementary | **Likely ✓ — needs check** | Two distinct matrices in SL_2(ℝ) generically generate a non-elementary group. Specific test below. |
| G3: moment | ✓ | Two atoms. |
| G4: algebraic (in SL_2(K)) | ✓ | K = ℚ(√7), entries in K. |
| G5: Zariski-dense | Likely ✓ — needs check | Two generic matrices generate Zariski-dense. |

**G2 specific check.** Strong irreducibility: M_+ and M_- have no common 1-dim invariant
subspace in ℝ² iff their commutator [M_+, M_-] = M_+ M_- M_+^{-1} M_-^{-1} is not the identity
and not a homothety. **In general, for irrational Δ entries**, this holds generically.
Proximality: M_+ has eigenvalues; tr(M_+) = (7 + 36)/45 + 0 = 43/45 (since Δ has zero diagonal),
det(M_+) = 1. Discriminant: (43/45)² - 4·1 = 1849/2025 - 4 = (1849 - 8100)/2025 = -6251/2025 < 0.
**M_+ has complex eigenvalues — it's an elliptic element, not hyperbolic.** Same for M_-.

If both M_+ and M_- are elliptic (rotations / conjugates of rotations), the generated group is
typically **inside SO(2, ℝ) or a conjugate thereof — i.e., compact**. Compact ⟹ not
non-elementary. **G2 likely FAILS for this specific construction.**

The trace constraint is forced by T_lead: tr(T_lead) = 43/45 < 2, so any det-1 lift with
zero-diagonal Δ has trace 43/45 and hence elliptic. To avoid elliptic, Δ must have nonzero
diagonal, redistributing trace.

Try Δ_11 = 1, Δ_22 = -1 (preserves the 2-param family's residual constraint up to recheck).
Then 36(1) + 7(-1) = 9 Δ_21 + 28 Δ_12, i.e., 29 = 9 Δ_21 + 28 Δ_12. Δ_12 Δ_21 - 1·(-1) = det(Δ),
det(Δ) = 1 ⟹ Δ_12 Δ_21 = 0. So one of {Δ_12, Δ_21} is zero. If Δ_12 = 0, Δ_21 = 29/9; if
Δ_21 = 0, Δ_12 = 29/28.

Take Δ_12 = 0, Δ_21 = 29/9. Then **Δ = [[1, 0], [29/9, -1]]**, all rational!

**M_+ = T_lead + Δ = (1/45)·[[7+45, 9], [28+145, 36-45]] = (1/45)·[[52, 9], [173, -9]]**
**M_- = T_lead - Δ = (1/45)·[[7-45, 9], [28-145, 36+45]] = (1/45)·[[-38, 9], [-117, 81]]**

Verify: tr(M_+) = (52 - 9)/45 = 43/45. det(M_+) = (52·(-9) - 9·173)/45² = (-468 - 1557)/2025 = -2025/2025 = **-1**.

Wait — det = -1, not 1! Re-derive. The linear-term-vanishing constraint gave the SUM of dets
equal to 2; not each det equals 1. Let me redo more carefully.

  det(T_lead + Δ) + det(T_lead - Δ) = 2 det(T_lead) + 2 det(Δ) = 0 + 2 det(Δ).

For both dets to equal 1: sum = 2, so det(Δ) = 1. ✓ already imposed.

  det(T_lead + Δ) - det(T_lead - Δ) = 2 · tr(adj(T_lead) Δ).

For both = 1: difference = 0, so tr(adj(T_lead) Δ) = 0. ✓ already imposed.

So my computation has det(M_+) = 1 by construction; let me redo the arithmetic.

  M_+ = T_lead + Δ. T_lead = (1/45)·[[7, 9], [28, 36]]. Δ = [[1, 0], [29/9, -1]].

  M_+ = [[7/45 + 1, 9/45 + 0], [28/45 + 29/9, 36/45 - 1]]
      = [[52/45, 9/45], [28/45 + 145/45, 36/45 - 45/45]]
      = [[52/45, 9/45], [173/45, -9/45]]
      = (1/45)·[[52, 9], [173, -9]]

  det(M_+) = (52·(-9) - 9·173)/45² = (-468 - 1557)/2025 = -2025/2025 = -1.

**det = -1, not 1.** The construction gave a det-(-1) matrix, not det-1. So my "both dets equal 1"
deduction is wrong — let me recheck.

Re-derive carefully:
  det(T_lead + Δ) = (T_lead + Δ)_{11} (T_lead + Δ)_{22} - (T_lead + Δ)_{12} (T_lead + Δ)_{21}
                  = T_11 T_22 + T_11 Δ_22 + Δ_11 T_22 + Δ_11 Δ_22 - T_12 T_21 - T_12 Δ_21 - Δ_12 T_21 - Δ_12 Δ_21
                  = det(T) + det(Δ) + (T_11 Δ_22 + Δ_11 T_22 - T_12 Δ_21 - Δ_12 T_21)

  Let L(Δ) := T_11 Δ_22 + Δ_11 T_22 - T_12 Δ_21 - Δ_12 T_21. Then det(T+Δ) = det(T) + det(Δ) + L(Δ).

For T = T_lead: T_11 = 7/45, T_22 = 36/45, T_12 = 9/45, T_21 = 28/45.
  L(Δ) = (7/45) Δ_22 + Δ_11 (36/45) - (9/45) Δ_21 - Δ_12 (28/45)
       = (7 Δ_22 + 36 Δ_11 - 9 Δ_21 - 28 Δ_12) / 45

For Δ = [[1, 0], [29/9, -1]]:
  L(Δ) = (7·(-1) + 36·1 - 9·(29/9) - 28·0) / 45 = (-7 + 36 - 29 - 0) / 45 = 0/45 = 0. ✓

  det(Δ) = 1·(-1) - 0·(29/9) = -1.

  det(T+Δ) = det(T) + det(Δ) + L(Δ) = 0 + (-1) + 0 = -1.

**det(Δ) = -1, not +1, so det(M_+) = -1.** My earlier deduction "det(Δ) = 1" was wrong; the
right constraint was det(M_+) = 1 ⟹ det(T) + det(Δ) + L(Δ) = 1 ⟹ 0 + det(Δ) + 0 = 1 ⟹
det(Δ) = +1.

Redo with Δ_11 = 1, Δ_22 = -1 forces Δ_12 Δ_21 - (-1) = det(Δ) = +1, so Δ_12 Δ_21 = 0.
If Δ_12 = 0, then det(Δ) = (1)(-1) - 0 = -1. ✗
If Δ_21 = 0, then det(Δ) = (1)(-1) - 0·Δ_12 = -1. ✗

**Diagonal (1, -1) forces det(Δ) = -1, incompatible with det(Δ) = +1.** Need diagonal (a, d)
with ad ≥ 1 to allow det(Δ) = +1. Try Δ_11 = 2, Δ_22 = 1: L = 7·1 + 36·2 - 9 Δ_21 - 28 Δ_12 = 79 - 9 Δ_21 - 28 Δ_12. Set L = 0: 9 Δ_21 + 28 Δ_12 = 79. det(Δ) = 2·1 - Δ_12 Δ_21 = +1 ⟹ Δ_12 Δ_21 = 1.

Two equations: 9 Δ_21 + 28 Δ_12 = 79, Δ_12 Δ_21 = 1. Let x = Δ_12, then Δ_21 = 1/x, so 9/x + 28 x = 79, i.e., 28 x² - 79 x + 9 = 0. Discriminant 79² - 4·28·9 = 6241 - 1008 = 5233. √5233 ≈ 72.34 (irrational).

So this Δ has irrational entries. To get rational Δ with det = +1 + nonzero diagonal preserving
L(Δ) = 0 is a Diophantine equation — sometimes solvable, sometimes not.

Try diagonal (3, 1): L = 7·1 + 36·3 - 9Δ_21 - 28Δ_12 = 115 - 9Δ_21 - 28Δ_12 = 0, so 9Δ_21 + 28Δ_12 = 115. det = 3 - Δ_12 Δ_21 = +1, so Δ_12 Δ_21 = 2. So 9Δ_21 + 28Δ_12 = 115 and Δ_12 Δ_21 = 2. With Δ_12 = x, Δ_21 = 2/x: 28x + 18/x = 115, 28x² - 115x + 18 = 0. Discriminant 115² - 4·28·18 = 13225 - 2016 = 11209. √11209 = 105.87... (not a perfect square — irrational).

Diagonal (4, 1): L = 7 + 144 - 9Δ_21 - 28Δ_12 = 151 - 9Δ_21 - 28Δ_12 = 0. det = 4 - Δ_12 Δ_21 = +1, so Δ_12 Δ_21 = 3. 9Δ_21 + 28Δ_12 = 151. With Δ_12 = x, Δ_21 = 3/x: 28x + 27/x = 151, 28x² - 151x + 27 = 0. Discriminant 151² - 4·28·27 = 22801 - 3024 = 19777. √19777 ≈ 140.6 — irrational.

Diagonal (5, 1): L = 7 + 180 - 9Δ_21 - 28Δ_12 = 187. 9Δ_21 + 28Δ_12 = 187, Δ_12 Δ_21 = 4. 28x² - 187x + 36 = 0. Discrim 187² - 4·28·36 = 34969 - 4032 = 30937. √30937 ≈ 175.9 — irrational.

The Diophantine equation 28 x² - K x + (a d - 1) · ... = 0 with K = (7 d + 36 a) tied to L = 0
generically gives irrational roots. Hence rational hyperbolic-preserving 2-atom lifts are
NON-GENERIC; algebraic-of-higher-degree lifts are required, expanding the algebraic number field
K = ℚ(√D) for D depending on the discriminant.

For our purposes: **a 2-atom lift of T_lead = E[M] with both atoms in SL_2(ℝ) exists, generally
with entries in a quadratic extension K = ℚ(√D) of ℚ.** Picking specific D values:

  D = 7 (from §A.5): K = ℚ(√7), entries in K.

### A.7 Phase 2 conclusion for Candidate A

A two-atom probability measure μ on SL_2(K) with E[M] = T_lead exists for K = ℚ(√7) (and many
other quadratic extensions). The Furstenberg framework's hypotheses G1, G3, G4 are passable; G2
(non-elementary) requires checking the trace of each atom.

**The construction in §A.5 with Δ = [[0, -3√7/14], [2√7/3, 0]] gives traces both equal to 43/45
which is less than 2 in absolute value, so BOTH atoms are ELLIPTIC.** Two elliptic generators in
SL_2(ℝ) almost always generate a non-elementary group BUT in this specific case the trace is
forced to be the same on both atoms (since Δ is traceless), and there's no asymmetry to break
proximality.

Whether the group generated by {M_+, M_-} is non-elementary depends on a more careful Tits-
alternative-style analysis. Heuristically: two elliptic elements of SL_2(ℝ) with the same trace,
related by an involution Δ → -Δ, generate either a finite group (if the elliptic angles are
commensurable with π) or a free non-abelian group (generically). In the latter case the group
is non-elementary.

**Conclusion for Candidate A.** A 2-atom μ on SL_2(K) for K = ℚ(√7) with E[M] = T_lead exists.
Whether {M_+, M_-} satisfies G2 (non-elementary) requires explicit trace-discriminant analysis,
and is GENERICALLY YES for a generic choice of Δ but specifically depends on whether the
elliptic angles are commensurable with π. **The deep issue is: even if the framework applies,
the resulting Furstenberg measure ν is on P^1, not on ℤ_3 or 𝕋. The link between ν and
μ_n is NOT given by the construction — T_lead is the EXPECTATION of M, not the projection of ν.**

The Furstenberg framework gives a measure ν such that ν = ∫ g_* ν dμ(g), which is a fixed
point of the AVERAGED projective action of g. The Syracuse measure μ_n is a fixed point of a
DIFFERENT operator — the Tao bilinear recursion. The two fixed-point equations are not the
same, so ν is not a priori μ_n. **The transfer (Gate T1) is the load-bearing question.**

---

## Candidate B: Tao-recursion 2-step matrix (the "natural" SL_2 candidate)

### B.1 Construction

Syracuse iteration is (3x + 1)/2^k. On the residue mod 3^n, this is a Markov chain on
(ℤ/3^n)^×. Tao's bilinear recursion gives Fourier-side dynamics:

    μ̂_{n+1}(ξ) = Σ_v 2^{-v} A_v(ξ) μ̂_n(ξ · 2^{-v} mod 3^n)

This is a **bilinear** recursion (the A_v factor depends on ξ; the μ̂_n is at a transformed
argument). It's not a linear SL_2 random walk; it's an integral operator with kernel that
depends on ξ.

Could the recursion be cast as an SL_2(ℝ) action on a 2D space? The natural 2D space is
(P_+, P_−) where P_a^{ab}(c) is the class-summed squared mass. The R76 §11 + R77 derivation
gives:

    (P_+, P_−)_{n+1} = T_lead · (P_+, P_−)_n + (higher-order corrections / cross-frequency)

T_lead is the within-level **expectation** of a 2x2 matrix-valued random variable indexed by v.
Specifically:

    M_v = (depends on v, and Tao recursion) and T_lead = Σ_v 2^{-v} E[M_v | v]

This is the same E[M] = T_lead structure as Candidate A, but with the v-distribution Geom(½) as
the underlying probability measure on a COUNTABLE-support set {M_v : v ∈ ℤ_{≥1}}.

### B.2 Are the M_v in SL_2?

Tao's recursion gives the matrix entries at v-step are W_+(g) = 2^{-g+1}/15 for g = v' - v, etc.
The transition operator at fixed v is essentially T_diag/2 plus an off-diagonal piece. Trace
and determinant of M_v:

  - For v = 1, 2, 3, ...: M_v's trace and det depend on the parity of v.
  - **No reason a priori for det(M_v) = 1.** In fact, T_diag = (1/5)·[[1, 1], [4, 4]] has det 0,
    so the "diagonal v-block" of M_v is rank-1. The full M_v at fixed v inherits rank-1 from
    T_diag's structure.

  - Indeed, T_lead = T_diag - (2/45)·[[1, 0], [4, 0]] (from T_LEAD_CORRECTED_SPECTRUM) is the
    SUM over v ≥ 1 of weighted M_v. Each individual M_v contribution within the Tao recursion
    is rank-1, since the only nonzero "image direction" is the (1, 4) class-mass ratio.

**The M_v at each v are individually rank-1, det = 0, not in SL_2(ℝ).** Hence the natural
v-branched lift puts the random matrices OUTSIDE SL_2 — not in it. G1 fails atom-by-atom.

### B.3 Could SL_2-membership be restored by composing two consecutive Tao steps?

Tao recursion has two natural "step" granularities:
  - one Syracuse iteration (× (3x+1)/2^v at one v)
  - two consecutive (or k-many) Syracuse iterations

The 2-step composition matrix M_{v, v'} (one step at v, another at v') is a product of two
rank-1 matrices. A product of two rank-1 matrices is rank-1 generally (unless they're aligned
in a specific way to give rank ≤ 2). Specifically: if M_v = (1, 4)^T u_v^T (column · row, with
u_v the per-step row vector), then M_{v'} M_v = (1, 4)^T u_{v'}^T · (1, 4) · u_v^T · (rank-1).

So 2-step compositions are still rank-1; det = 0; not in SL_2.

**Any finite product of M_v matrices is rank-1.** This is structural: the Syracuse-Tao bilinear
recursion has only ONE non-trivial direction of action (the (1, 4) class-ratio direction).
There's no second non-trivial direction for SL_2 to act on.

### B.4 Gate checks for Candidate B

| Gate | Status |
|---|---|
| G1: det = 1 | **FAILS** — individual M_v are rank-1 (det = 0); products still rank-1. |
| G2: non-elementary | N/A — G1 fails first. |
| G4: algebraic | ✓ (rational entries) but irrelevant since G1 fails. |
| **Verdict** | **G1 fails atom-by-atom.** The Tao recursion's 2D structure is fundamentally rank-1; doesn't lift to SL_2(ℝ). |

---

## Candidate C: Projective P^1 action (decouple from SL_2 membership)

### C.1 Construction

Hochman-Solomyak's machinery is fundamentally about the projective action on P^1(ℝ). Even if
T_lead isn't in SL_2(ℝ), perhaps its **projective action** [T_lead]_* on P^1 is well-defined
and dynamically interesting.

The projective action of a 2×2 matrix M ∈ GL_2(ℝ) on P^1(ℝ) is:
  [M](x : y) := (M·(x, y)^T)_{normalized}.

This is well-defined IF M is invertible (det ≠ 0). For T_lead with det = 0, **the projective
action is NOT defined as a map P^1 → P^1.** Instead, T_lead maps all of P^1 to the single point
[1 : 4] (the image of the rank-1 operator), except for the null direction [9 : -7] which maps
to 0 (undefined in P^1).

This is a **degenerate projective map**, not in PSL_2(ℝ). It cannot serve as a generator of a
projective random walk in the standard sense.

### C.2 Gate checks for Candidate C

| Gate | Status |
|---|---|
| G1 (projective version): well-defined map P^1 → P^1 | **FAILS** — T_lead collapses P^1 to a single point [1:4] (and breaks at [9:-7]). |
| **Verdict** | **G1 fails projectively too.** Rank-1 matrices don't admit non-degenerate projective actions. |

---

## Candidate D: Re-interpret the question — random walk in HIGHER dimension

### D.1 Construction

If the 2D operator T_lead is rank-1, perhaps the natural ambient dimension is higher. The
cross-frequency derivation lives on V_M = span{M_n^{ab}(g, c) : g ∈ {0, 2, 4, ...}}, an
infinite-dimensional Hilbert space (or its finite truncations V_M^{(g_max)} = finite-dimensional
truncation at index g ≤ g_max).

**Is V_M^{(g_max)} a natural ambient space for an SL_d(ℝ) random walk (with d = dim V_M^{(g_max)})?**
Possibilities:

  - **d = 2 (the (1, 4)-projection of V_M).** This is what T_lead lives on. As above, it's
    rank-1 in this projection — too low-dimensional for non-elementary action.

  - **d = 6, 7, ...** (V_M^{(g_max)} at g_max = 2, 4, ...). The T_V operator on V_M^{(g_max)}
    is the candidate higher-rank operator. From T_V_DISPOSITION: T_V doesn't close on
    V_M^{(g_max)} for any finite g_max — the cascade generates **odd-G** moments and **phase
    offsets θ_{v,g}** outside V_M. So V_M as currently constructed isn't a closed subspace
    under iteration.

  - **d = ∞**: the full closure under iteration is infinite-dimensional, with phase offsets θ in
    ℚ/ℤ ⊗ ℤ[1/2] and shift indices G in ℤ. The full closure is essentially an arithmetic
    group action on a larger space — possibly a torus 𝕋^∞ or a profinite ℤ_3 × ℤ_2 ⊗ (...)
    structure.

### D.2 He-de Saxcé alternative: torus random walk

The natural higher-dim setting is **𝕋^d for d > 1** (rather than P^1). He-de Saxcé's framework
applies. Specifically:

  - Tao recursion induces an action on ℤ/3^n × ℤ/2^v (a 2-d profinite quotient capturing both
    the 3-adic Fourier modes and the 2-adic valuation step). This is structurally a Markov
    chain on a 2-d torus.

  - SL_d(ℤ)-action on 𝕋^d? The Syracuse update (3x+1)/2^v involves multiplication by 3 and
    division by 2^v, which in the 3-adic + 2-adic joint setup is a 2×2 matrix
    [[3, 0], [0, 1/2^v]] (block-diagonal). **det = 3/2^v ≠ 1 in general** — so this isn't
    SL_2(ℤ) either.

  - However, multiplication-by-3 is an SL_1-action on the 3-adic factor, and division-by-2^v
    is an SL_1-action on the 2-adic factor. The product is an action on 𝕋^2 by diagonal
    matrices with non-unit determinant — not SL_2, not in He-de Saxcé.

### D.3 Gate checks for Candidate D

| Gate | Status |
|---|---|
| G1 (V_M version): T_V is rank-1 in (P_+, P_-) projection; higher-rank on V_M^{(g_max)} BUT not closed under iteration | **FAILS — V_M doesn't close** (T_V_DISPOSITION's H_M_RECURSION_UNDERSPECIFIED) |
| G1 (torus action version): Syracuse-update is diag(3, 1/2^v) — not in SL_2(ℤ) (det ≠ 1) | **FAILS — det isn't 1** |
| **Verdict** | **Both higher-dim candidates fail their det = 1 / closure prerequisite.** |

---

## Phase 2 verdict on candidates

| Candidate | G1 (det) | G2 (non-elem) | G4 (algebraic) | Verdict |
|---|---|---|---|---|
| A: direct lift via E[M] = T_lead, K = ℚ(√7) | ✓ (constructed) | likely ✓ (generic) | ✓ | Construction-dependent. T_lead is at best the FIRST MOMENT of a measure on SL_2; the framework gives ν on P^1 which is **not μ_n**. The transfer (T1) is the structural question. |
| B: Tao recursion v-branched M_v | **FAILS** (det = 0) | N/A | ✓ | Tao's recursion gives rank-1 atoms, not SL_2 atoms. |
| C: T_lead as projective action | **FAILS** (T_lead degenerate projectively) | N/A | ✓ | Rank-1 has no non-degenerate projective action. |
| D: higher-dim V_M / 𝕋^2 | **FAILS** (closure / det = 1) | N/A | ✓ | V_M doesn't close (per T_V_DISPOSITION); torus matrix has det ≠ 1. |

**Only Candidate A has any path through G1.** Its existence relies on an artificial 2-atom
construction with E[M] = T_lead in K = ℚ(√7); the construction is NOT natural to the Tao
recursion (which has rank-1 atoms, not SL_2 atoms).

**The natural "T_lead = E[M] for some non-natural μ" reframing** moves T_lead from a Syracuse-
intrinsic object to a constructed object whose connection to μ_n is NOT given.

This is the load-bearing finding of Phase 2: **the natural SL_2(ℝ) extensions of Syracuse
dynamics (Candidates B, C, D) all fail G1; Candidate A passes G1 only by giving up the natural
connection to Syracuse.**
