# GENERALIZATION_R78_SPECIFIC_FEATURES — Phase 4 A3 audit

**Purpose:** identify each R78-specific feature that makes a step of the Path 2 chain work, so we can judge which features transfer and which don't.

## Feature inventory

### Feature 1: Prime-power modulus on principal-unit coset

**Where it appears:** Step 1 (closed-form magnitude theorem), Step 2 (bijection).

**Why R78 needs it:** Cochrane's Theorem 2 requires polynomial-phase character sums on `Z/p^n`. The principal-unit coset `a ≡ 1 mod p` admits Postnikov decomposition.

**Generality:** This is a **broad class of problems** — covers most prime-power-moduli character sum literature.

**Transfer:** **GENERIC** — most candidates inhabit this setting natively.

---

### Feature 2: Cochrane Theorem 2 (polynomial-phase prime-power character sum)

**Where it appears:** Step 1.

**Why R78 needs it:** The closed-form magnitude `|G_p(a)| = √q` requires a polynomial-phase decomposition admitting Cochrane's saddle-point analysis.

**Generality:** Cochrane's theorem is a **published classical result** (2002, refined 2003 with Pinner). Available to anyone in the field.

**Transfer:** **GENERIC** — Candidates 1, 4 use it directly. Candidates 2, 3 need Postnikov-style substitution first (mostly works).

---

### Feature 3: Cubic phase degree

**Where it appears:** Step 2 (multi-saddle bijection), Step 3 (saddle exactness at r = J = 3), Step 4 (linear-in-second-digit structure after substitution).

**Why R78 needs it:** A cubic phase has the form `P(s) = α s³ + β s² + γ s + δ`. Its derivative `P'(s) = 3αs² + 2βs + γ` is QUADRATIC, with multiple roots — this gives the `p`-fold bijection in Step 2. Quadratic phases have a SINGLE saddle (no bijection); quartic-or-higher phases have too-many saddles (bijection ill-defined or over-counted).

The "linear-in-c_2" structure of Step 4 comes from: at saddle `s*`, the second base-p digit `c_2` of the parameter enters the cubic term `p³ · c_2 · s*` linearly. This is **algebraically specific to cubic phases** — at quartic phases, the second-digit dependence is QUADRATIC not LINEAR, breaking the clean Plancherel collapse.

**Generality:** Cubic phases are a specific structural class within the polynomial-phase literature. Quadratic phases (Heilbronn-on-coset, Candidate 3) FAIL the chain; quartic phases would also fail.

**Transfer:** **R78-SPECIFIC TO CUBIC PHASES.** This is the load-bearing structural feature.

---

### Feature 4: Saddle exactness at r = 3 (J_p = r alignment)

**Where it appears:** Step 3.

**Why R78 needs it:** The closed form `G_p(a) = √q · e_q(P_a(s*))` is EXACT only when the p-adic log truncation level `J_p` matches the polynomial-degree-times-r relationship. At r = 3 (q = p^4), `J_p = 3` matches the cubic-degree exact-saddle condition.

At r ≥ 4, the truncation level no longer matches and the saddle gives only an approximate phase — R78's Hensel correction `D(a)` quantifies the deviation.

**Generality:** The J = r alignment is a tunable parameter. Different problems achieve it at different r values:
- R78 cubic at r = 3 ↔ J = 3.
- Cochrane–Pinner cubic at n = 6 ↔ ⌈n/2⌉ = 3 (Candidate 1).
- Heath-Brown after Postnikov at r = 5, 6 (Candidate 2).

**Transfer:** **STRUCTURALLY ALIGNS** for cubic phases at the right r value; doesn't transfer to non-cubic phases (Candidate 3).

---

### Feature 5: Single-parameter principal-unit-coset parameterization

**Where it appears:** Step 2 (bijection a ↔ C_a), Step 4 (Inner-Plancherel collapse on c_2).

**Why R78 needs it:** The parameter `a` is a single base-p-decomposable index `a = 1 + p · α_1 + p² · α_2 + ...`. The Inner-Plancherel step collapses on `α_2` because the post-saddle phase is LINEAR in α_2 alone. If the parameter were two-dimensional `(a, b)`, the post-saddle phase would mix `(a_2, b_2)` and 1D Plancherel wouldn't apply.

**Generality:** R78's bilinear sum has ONE outer parameter (`a` ranging over the principal-unit coset). Heath-Brown's bilinear has TWO parameters (χ, a) — leading to Step 4 FAILURE on Candidate 2.

**Transfer:** **R78-SPECIFIC TO ONE-PARAMETER BILINEAR SETUPS.** Two-parameter bilinears (Heath-Brown, Banks-Shparlinski, joint χ × a) break Step 4.

---

### Feature 6: 1/sin grid identity on principal-unit coset of size p

**Where it appears:** Step 5.

**Why R78 needs it:** The outer sum after Inner-Plancherel is `Σ_{α ∈ Z/p} 1/|sin(π(1+pα)/p²)|`. This is bounded by `p + 2 log p ≤ 2p` via the cosecant grid identity.

**Generality:** Classical Dirichlet-kernel manipulation — applies to ANY sum of csc on a regular arithmetic progression of length p.

**Transfer:** **GENERIC** — transfers freely whenever the outer sum reduces to a cosecant grid.

---

### Feature 7: Truncated p-adic log phase structure

**Where it appears:** Step 1 (phase derivation), Step 3 (saddle position formula), Step 4 (post-saddle expansion).

**Why R78 needs it:** The phase `P_a(s) = ps − C_a · L_p(1+ps)` uses the p-adic log of `1 + ps` truncated at J_p terms. This SPECIFIC choice of phase function comes from the Collatz iteration `c · 4^u = c · (1+3)^u = c · exp(u · L_3(1+3))` rewriting.

**Generality:** Any Postnikov-decomposed character sum on the principal-unit coset has phase = (coefficient) × (truncated p-adic log). This covers:
- Generic exponential sums of the form `Σ e_q(c · g^u)` for `g ≡ 1 mod p` (Candidate 4).
- After Postnikov substitution: many character sums on `(Z/p^n)×` (Candidate 2 — but with two parameters mixing).

**Transfer:** **GENERIC for one-parameter Postnikov sums.** Distorted for two-parameter Postnikov-with-multiplicative-character sums.

---

## Load-bearing feature ranking

| Feature | Load-bearing for chain? | Transfers? | Failure mode if missing |
|---------|-------------------------|------------|-------------------------|
| 1. Prime-power modulus | Yes | Yes (generic) | Chain doesn't apply (out of scope) |
| 2. Cochrane Theorem 2 | Yes | Yes (classical) | Step 1 fails (no closed form) |
| **3. Cubic phase degree** | **Yes** | **R78-specific** | **Step 2 + Step 4 fail (no multi-saddle / wrong digit structure)** |
| 4. Saddle exactness at r = J | Yes | Aligns for cubic | Step 3 gives approximation only |
| **5. One-parameter principal-unit setup** | **Yes** | **R78-specific** | **Step 4 fails (2D Plancherel needed)** |
| 6. 1/sin grid identity | Yes | Yes (classical) | Step 5 needs different bound |
| 7. Truncated p-adic log phase | Yes | Mostly generic | Step 1 needs different decomposition |

## The two R78-specific features

**Features 3 (cubic phase) and 5 (one-parameter setup) are the TWO LOAD-BEARING R78-SPECIFIC FEATURES.**

The chain transfers iff a candidate has:
- **(3) Cubic-degree polynomial phase** after Postnikov substitution to principal-unit-coset coordinates. Quadratic fails (Heilbronn — Candidate 3); quartic-or-higher would also fail.
- **(5) Single-parameter principal-unit-coset bilinear structure**. Two-parameter (Heath-Brown — Candidate 2) fails.

Cochrane–Pinner cubic sums (Candidate 1) satisfy BOTH and the chain transfers. R78-variant Postnikov sums (Candidate 4) satisfy BOTH trivially. Outside this narrow specification, the chain fails.

---

## Implication for methods paper viability (A4 advance)

The chain's **distinct novelty over published Cochrane–Pinner machinery** is:
- The Inner-Plancherel step (Step 4) — exploiting the linear-in-second-digit structure to collapse the inner sum to a length-p Dirichlet kernel.
- This is **one technical lemma**, not a standalone methodology.

Steps 1, 2, 3, 5 are all classical / Cochrane–Pinner published machinery. Step 4 is the new wrinkle.

A methods paper would essentially be: "Cochrane–Pinner cubic exponential sums on principal-unit cosets admit a bilinear bound via an Inner-Plancherel collapse on the second p-adic digit of the parameter, when the phase is cubic and the parametrization is one-dimensional."

**Honest scoping:** this is a TECHNICAL LEMMA, not a methods paper. The scope is "cubic Postnikov phases on principal-unit cosets" — a specific class containing R78 and not much else. Specialists familiar with Cochrane–Pinner would view it as a "competent application + one minor technical wrinkle", not a new methodology.

**Conclusion:** Features 3 and 5 being load-bearing AND R78-specific means the chain transfers only to a narrow class of cubic-on-coset problems. The methods-paper hypothesis (H_FULL or H_PARTIAL_TEMPLATE with broad scope) is not supported.
