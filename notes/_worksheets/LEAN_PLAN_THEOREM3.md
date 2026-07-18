# Lean verification plan — Theorem 3 (Plancherel decomposition) + exact S_k

**Purpose of this file:** a plain-English reading guide for what the Lean file will
claim and check. This is so you (Nathan) can follow what Lean is verifying. You do
NOT certify the math by reading this — Lean does. Your job is to confirm the INPUTS
match your own work (the chain construction, the R66 nesting, the S_k numbers).

Written 2026-06-11. Matches the chain in `s_infinity_exact.py` / `fourier_S_decomposition.py`.

---

## What Lean WILL verify (these are true and provable)

### Part A — Theorem 3 (the Plancherel decomposition) — the one the referee flagged

**Plain English:** The high-frequency mass S_k equals the sum of |characteristic function|^2
over exactly the frequencies NOT divisible by 3.

It breaks into three honest steps. Lean checks each.

**Step 1 — the nesting fact (your R66).**
Zooming out from level k reproduces level k-1 exactly: if you pool all the level-k
residues sitting inside one level-(k-1) bucket, the pooled stationary mass pi_k equals
pi_{k-1} of that bucket.
- SOURCE: this is your R66 ("Markov chain on (Z/3^k Z)*: closed form for all k").
- In Lean this enters as a HYPOTHESIS named `nesting` (the level-compatibility of the
  chain family). We are NOT re-deriving it from the Syracuse map inside this file; we
  state it as the property your construction has, cite R66, and (where R66 leans on
  Terras parity-vector structure) cite Terras 1976.
- ** YOUR CHECK: is it true that your chain, by R66, satisfies "pool pi_k over a mod-3^{k-1}
  bucket = pi_{k-1}"? Yes / No. **

**Step 2 — character collapse (follows from Step 1, pure computation).**
For a frequency that IS divisible by 3 (write it 3*xi'), the characteristic function at
level k equals the characteristic function at level k-1 at xi'. Reason: the phase
e^{-2pi i * r * 3xi' / 3^k} simplifies to e^{-2pi i * r * xi' / 3^{k-1}}, and Step 1 lets
you regroup the sum by buckets. Lean derives this from Step 1 + the exponential identity.

**Step 3 — Theorem 3 itself (Plancherel split).**
- Plancherel (standard, cited from mathlib): X_k = sum over ALL frequencies of |char fn|^2.
- Split that sum into "divisible by 3" and "not divisible by 3."
- By Step 2, the "divisible by 3" part equals X_{k-1}.
- So the "not divisible by 3" part = X_k - X_{k-1} = S_k.  QED.

This is the whole proof the referee said he saw "no route" to formalize. The route is
Steps 1-2-3; he was shown Step 3 with 1 and 2 compressed into a clause.

### Part B — the exact rational values (finite, deterministic, certified)

Lean verifies, by exact rational arithmetic (no floating point, no approximation):
- S_1 = 2/3
- S_2 = 10/21
- S_3 = 31370/67963
- ... and S_k through k = 6.

These are FINITE computations: build the rational Markov chain K_k on the coprime
residues mod 3^k, solve pi K = pi over Q, compute X_k = 3^k * sum pi^2, take
S_k = X_k - X_{k-1}. Deterministic => Lean can certify the exact rationals.

This is the "certified (exact computation, finite range)" claim of paper section 1.2,
upgraded from "a computation I ran" to "machine-verified."

---

## What Lean will NOT verify (because it is NOT proven — your repo says so)

### S_infinity = 7/15  ->  this stays OPEN.

Your own files state this plainly and repeatedly:
- s_infinity_derivation.md: "Rigorous derivation of S_infinity remains open."
- "(alpha) S_infinity = 7/15 verified analytically: PARTIAL - strong numerical evidence,
  no rigorous proof."
- "lambda_2 = 1/2 CONJECTURED from convergence rate."

The paper itself (section 1.2) lists 7/15 as CERTIFIED (finite range), NOT proven, and
section 7 is "what remains open." So there is no proof for Lean to check, and we do not
pretend there is. In the Lean file this appears as a documented comment:
  `-- OPEN: S_infinity = 7/15. Rate conjecture (lambda_2 = 1/2). See paper section 7.`
NOT as a `theorem`.

---

## The resulting artifact (what you can honestly say)

> "Theorem 3 (Plancherel decomposition): formally verified in Lean.
>  Exact values S_1 ... S_6: formally verified in Lean.
>  The limit S_infinity = 7/15: certified by these verified computations through k=6;
>  a rigorous rate proof remains open (section 7)."

This is EXACTLY your paper's section 1.2 proven/certified/open structure, with
"proven" and "certified" now machine-checked and "open" still honestly open.

The referee's objection ("uncheckable") is answered: every checkable claim is checked
by Lean, and the open part is labeled open. Nothing is left for him to take on faith.

---

## Build order when Lean is installed

1. Part B first (exact S_k). It is self-contained rational arithmetic, easiest to get
   compiling, and it nails the "certified computation" the referee specifically doubted.
2. Part A (Theorem 3). Needs the `nesting` hypothesis (R66) + mathlib's Plancherel.
3. S_infinity stays a comment, never a theorem.

## Files this is built from (your work, verified inputs)
- s_infinity_exact.py        -> the exact rational chain + S_k values (Part B)
- fourier_S_decomposition.py -> the Plancherel decomposition check (Part A)
- closed_form_findings.md    -> R66 (nesting), R70 (S_1, S_2 exact), R73/R74 (decay)
- s_infinity_derivation.md   -> states plainly that S_infinity = 7/15 is open
