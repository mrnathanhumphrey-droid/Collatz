# RESULT — EULER-T2: `μ = H₁ − H₂` (drift = Shannon − collision entropy). TRUE for Collatz, DISPOSED as a discovery (trivial-by-construction ⊕ known criticality). (2026-08-28)

**Probe:** `probes/probe_euler_t2.py`, exact/numeric to 1e-12. First target of the Euler internal-relations
hunt (warm-up to shake down the express→disprove method). Candidate T2 from `notes/CORPUS_LEDGER.md §6`,
derived by pairing the stopping-time and entropy ledgers.

## Express
Halving law `P(v)=2^{−v}`, v≥1 (Geom ½). Shannon `H₁=−Σ P log P`, collision `H₂=−log Σ P²`,
drift `μ = E[v]·log2 − log3` (map `x↦(3x+1)/2^v`). Claim: **`μ = H₁ − H₂`**.

## Verify (exact)
`H₁ = log 4 = 1.386294361120`, `H₂ = log 3 = 1.098612288668`, `μ = log(4/3) = 0.287682072452`.
`μ − (H₁−H₂) = 0.00e+00`. Identity holds.

## Disprove-or-classify — the family-law test (same test that killed the cyclotomic-7)
Generalize to `(q,p)` with valuation law `r=1/p`. Closed forms:
`H₁ = −log(1−r) − (r/(1−r))log r`, `H₂ = log((1+r)/(1−r))`, so with `r=1/p`:
```
   H₁ − H₂ = (p/(p−1))·log p − log(p+1)        [PURE denominator/valuation-side function; no q]
   μ       = (p/(p−1))·log p − log q            [contains the multiplier q]
   ⟹  μ = H₁ − H₂   ⟺   q = p+1
```
Verified p=2,3,4,5,7: equality holds exactly at `q=p+1` and fails otherwise. So it **is** a family law —
but along the curve **`q = p+1`**, which is **NOT** the criticality curve `q=(p+1)/(p−1)`. The two coincide
**only at p=2** (both give q=3). Collatz sits where "q=p+1" and "critical" happen to be the same point.

## Decomposition — it factors into known pieces
- **(a) `μ = H₁ − log q` is DEFINITIONAL:** `H₁ = E[v]·log2 =` the mean log-denominator (since `−log P(v)=v·log2`),
  and `μ = E[v]log2 − log q` by definition. Verified `H₁ = E[v]·log2` and `H₁ − log3 = μ`.
- **(b) `H₂ = log q` at q=3 is already banked** — the collision-entropy criticality identity
  (`H₂(Geom½)=−log Σ4^{−v}=log3=log q`, CORPUS_LEDGER §1).

Hence `μ = H₁ − H₂ = (H₁−log q)_{definitional} + (log q − H₂)_{=0 by criticality}`.

## Verdict
**TRUE identity, DISPOSED as a discovery.** It is a clean *restatement* — "the Syracuse drift is the
Shannon-minus-collision entropy gap of the halving law" — worth stating once, but:
1. it lives entirely in the **leading/first-moment/entropy layer**; it carries **no tail information** and does
   not reach `S_∞≈0.475` (it is about `μ=log(4/3)`, a leading-order object);
2. it **factors into one trivial-by-construction identity + one already-banked criticality fact** — it moves no
   new information between the two faces, which is the bar for a real Euler identity;
3. its "family law" is `q=p+1`, a different curve from criticality — Collatz is just where they meet.

**Method validated:** express → generalize → find the pinning condition (`q=p+1`) → decompose into known
pieces → classify. A cross-thread numerical coincidence that looked Euler-shaped dissolved under
generalization. Exactly the discipline for T1.

**Mildly-interesting residue (not a discovery, logged):** `H₁ − H₂ = (p/(p−1))log p − log(p+1)` is a pure
function of the halving/denominator law alone (no `q`); its equality with the `q`-bearing drift is the
arithmetic pin `q=p+1`. That a finite-place/entropy quantity numerically *predicts* `log(p+1)` where the
archimedean drift substitutes `log q` is the kind of place↔place bookkeeping the two-walls thread tracks —
but here it is exact bookkeeping, not a new bridge.

**Not at stake:** S_∞, MAHLER spine, the structure/value split. T2 is a leading-layer restatement only.
