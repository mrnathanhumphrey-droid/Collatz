# Q-sweep Test 3 — multi-resolution orthogonal decomposition for qx+1

**Date:** 2026-05-04. **Outcome:** **DECOMP-UNIVERSAL**.

Tests whether R77.5's lift-residual decomposition
`R_k^(q) := pi_{k+1}^(q) − T_q(pi_k^(q))` is orthogonal to `T_q(V_k^(q))` in L²(Z/q^{k+1}) and whether the identity

    ‖R_k^(q)‖² · q^k = S_{k+1}^(q) / q       (where S_{k+1}^(q) := X_{k+1}^(q) − X_k^(q),  X_j^(q) := q^j · Σ pi_j^(q)²)

proved over Q at q=3 (R77.5 follow-up) extends to q ∈ {3, 5, 7, 11, 13}.

---

## 1. Theoretical input

The argument from `result_77_5_d_R_identity_check.md` is q-blind:
the lift map T_q : V_k^(q) → V_{k+1}^(q), defined by
`T_q(v)(r') := v(r' mod q^k) / q` for r' coprime in Z/q^{k+1}, satisfies

    ‖R_k^(q)‖² = Σ pi_{k+1}^(q)² − (1/q) Σ pi_k^(q)²

by **marginal consistency** of the projective Markov system Σ_{r' lifts of r} pi_{k+1}^(q)(r') = pi_k^(q)(r). Multiplying by q^k:

    ‖R_k^(q)‖² · q^k = q^k · Σ pi_{k+1}² − q^{k−1} · Σ pi_k²
                     = X_{k+1}^(q)/q − X_k^(q)/q  +  X_k^(q)/q − X_k^(q)/q  ... no, more directly:
                     = (X_{k+1}^(q) − X_k^(q)) / q = S_{k+1}^(q) / q.

So Test 3 is a regression check (orthogonality and identity must both hold by structure), not a discovery test.

## 2. Stage 0 — q=3 regression

Canonical values from `result_77_5_d_R_norms.csv`:

| k | canonical ‖R_k‖² | this run ‖R_k‖² | match |
|---|---|---|---|
| 1 | 10/189 | 10/189 | PASS |
| 2 | 31370/1835001 | 31370/1835001 | PASS |
| 3 | 5303542579979870/925406323431537423 | 5303542579979870/925406323431537423 | PASS |

**Stage 0 verdict:** PASS — q-generalized code reproduces canonical q=3 values exactly over Q.

## 3. Stage 1 — orthogonality at q ∈ {3, 5, 7, 11, 13}

For each (q, k), tested ⟨R_k^(q), T_q(v)⟩ = 0 over Q for three test vectors v ∈ V_k^(q):
- v_1 = pi_k^(q) itself
- v_2 = δ at the first coprime state in Z/q^k
- v_3 = balanced ±1 indicator on the first two coprime states

| q | k | N_k | N_{k+1} | ⟨R_k, T(pi_k)⟩ | ⟨R_k, T(δ)⟩ | ⟨R_k, T(±)⟩ | pass |
|---|---|---|---|---|---|---|---|
| 3 | 1 | 2 | 6 | 0 | 0 | 0 | PASS |
| 3 | 2 | 6 | 18 | 0 | 0 | 0 | PASS |
| 3 | 3 | 18 | 54 | 0 | 0 | 0 | PASS |
| 5 | 1 | 4 | 20 | 0 | 0 | 0 | PASS |
| 5 | 2 | 20 | 100 | 0 | 0 | 0 | PASS |
| 7 | 1 | 6 | 42 | 0 | 0 | 0 | PASS |
| 7 | 2 | 42 | 294 | 0 | 0 | 0 | PASS |
| 11 | 1 | 10 | 110 | 0 | 0 | 0 | PASS |
| 13 | 1 | 12 | 156 | 0 | 0 | 0 | PASS |

**Stage 1 verdict:** all PASS — orthogonality holds as exact rational equality at every tested (q,k).

## 4. Stage 2 — identity ‖R_k^(q)‖² · q^k = S_{k+1}^(q) / q

All comparisons over Q via fractions.Fraction.

| q | k | ‖R_k^(q)‖² (decimal) | ‖R_k‖²·q^k | S_{k+1}/q | identity? |
|---|---|---|---|---|---|
| 3 | 1 | 5.291005e-02 | 0.1587301587 | 0.1587301587 | PASS |
| 3 | 2 | 1.709536e-02 | 0.1538582268 | 0.1538582268 | PASS |
| 3 | 3 | 5.731042e-03 | 0.1547381361 | 0.1547381361 | PASS |
| 5 | 1 | 5.469010e-02 | 0.2734504764 | 0.2734504764 | PASS |
| 5 | 2 | 1.813455e-02 | 0.4533636626 | 0.4533636626 | PASS |
| 7 | 1 | 8.685784e-02 | 0.6080048599 | 0.6080048599 | PASS |
| 7 | 2 | 2.898455e-02 | 1.4202430479 | 1.4202430479 | PASS |
| 11 | 1 | 8.097870e-02 | 0.8907657283 | 0.8907657283 | PASS |
| 13 | 1 | 8.553257e-02 | 1.1119234044 | 1.1119234044 | PASS |

**Stage 2 verdict:** all PASS.

### Test 2 cache cross-check

X_j^(q) values from this run match Test 2's cached values (where overlap exists):

| q | k | X_k matches cache | X_{k+1} matches cache |
|---|---|---|---|
| 3 | 1 | True | True |
| 3 | 2 | True | True |
| 3 | 3 | True | True |
| 5 | 1 | True | True |
| 5 | 2 | True | True |
| 7 | 1 | True | True |
| 7 | 2 | True | True |
| 11 | 1 | True | True |
| 13 | 1 | n/a | n/a |

## 5. Stage 3 — outcome classification

**Outcome:** **DECOMP-UNIVERSAL**

Multi-resolution decomposition is q-universal.

## 6. Strategic implication

R77.5's multi-resolution / wavelet-style geometric framework — V_{k+1}^(q) = T_q(V_k^(q)) ⊕ W_k^(q) with R_k^(q) ∈ W_k^(q) — extends cleanly to the qx+1 family. The identity ‖R_k^(q)‖² · q^k = S_{k+1}^(q) / q is structural, inherited entirely from marginal consistency of the projective Markov system modulo q^k.

Consequences:
- The convergence S_∞^(q) (Test 2's open question) is equivalent to summability of ‖R_k^(q)‖² · q^k in the multi-resolution decomposition.
- Every analytical tool R77.5 provides for q=3 (orthogonal complement chain, wavelet-style basis, transfer-operator route on Ẑ_q^×) is available unmodified for any odd prime q.
- The q-blindness of marginal consistency means c_q rationality (if any) sits in the projective limit of pi_k^(q), not in q-specific arithmetic.

## 7. Output files

- `result_q_sweep_test_3.py` — q-generalized R_k computation script
- `result_q_sweep_test_3_norms.csv` — exact-rational table per (q, k)
- `result_q_sweep_test_3_decomposition.md` — this writeup
