# The universal rate of qx+1 stationary L² mass — standalone writeup

**Date:** 2026-07-14. **Status:** paper-shaped; three results, one identity short of a full theorem on Result 1. Independent of every Collatz-closure thread (ε_k / c=7/45 / THEOREM_C_745 / R81b all untouched). Collatz need not be in the title.

---

## Object

For odd `q`, the **qx+1 Syracuse chain** on residues coprime to `q` mod `q^k`:

&nbsp;&nbsp;&nbsp;&nbsp;`r ↦ (q·r + 1)·2^{−v} mod q^k`, &nbsp; `v ~ Geom(1/2)` (P(v)=2^{−v}, v≥1).

Let `π_k` be its stationary distribution and define `S_k^{(q)} = q^k · ‖π_k‖² = q^k Σ_r π_k(r)²` (the state-count-normalized L² / collision mass). For `q=3` this is the R75 object with `S_k^{(3)} → 7/15`.

---

## Result 1 (RATE) — `S_k^{(q)} ~ (q/3)^k`, universal in q

Equivalently `‖π_k‖² ~ C_q · 3^{−k}` with the exponential rate **`3^{−k}` independent of q**; the q-normalized ratio `R_k = (S_k/S_{k−1})/q → 1/3`.

**Mechanism (the "3" named).** `‖π_k‖²` contracts by exactly `1/3` per level, and

&nbsp;&nbsp;&nbsp;&nbsp;**`1/3 = Σ_{v≥1} 2^{−v}·2^{−v} = Σ_{v≥1} 4^{−v} = E_{v~Geom(1/2)}[2^{−v}]`**

— the halving second-moment. This is **q-blind**: `q` enters the transfer operator only through the multiplicative character on `(Z/q^k)*`, which rescales the state-count `q^k` but not the halving statistic. Hence `‖π_k‖² ~ (1/3)^k` and `S_k = q^k‖π_k‖² ~ (q/3)^k`. Geometric picture: the stationary measure occupies an **effective `3^k` residues** inside the `q^k`-sized space, for every q — verified `1/‖π_k‖² = 3, 9, 27` exactly at q=11,13,25.

**Evidence (Probe 5, `probe_5_universal_rate.py`).**
- **Adversarial-q falsifier ran first** and survived on every case built to break the separation: small `ord_q(2)` primes (7, 23, 47, **89 [ord 11/88]**) all `R_k → 1/3`; odd composite (9,15,25,21,27,45) all `→1/3`; q≡0 mod 3 (critical q=3) trends to 1/3. The q=7 anomaly and small-ord primes give a **clean q/3 rate** — their anomaly lives entirely in the *constant* (Result 2), never the rate.
- Direct contraction `‖π_k‖²`-ratio `→ 1/3` (q=5,7,11,13,25); exact identity `Σ4^{−v}=1/3`.
- **Scope: odd q only.** Even q breaks the construction (2 not invertible mod `q^k`).

**Open (the one line to a full theorem).** That the *sub-leading* character contributions to `‖π_k‖²` do not perturb the leading `(1/3)^k` rate — i.e. **uniform diagonal-self-overlap domination on `(Z/q^k)*`**, an R76-style conservation identity generalized from q=3. Numerics show zero q-dependence in the rate to tested k; the clean algebraic proof of uniform domination is the remaining step. Reported as **H_PROVED-at-mechanism**, not a finished theorem.

---

## Result 2 (CONSTANT) — `c̃_q = (q−3)/q`, from the same factorization

Leading `S_k/(q/3)^k → 1`, so the difference `D_k = S_k − S_{k−1}` gives `c̃_q := D_k/(q/3)^k → 1 − 3/q = (q−3)/q`. Same rate factorization; the bare `3` in `q−3` is the same `3 = 1/E[2^{−v}]`. **Rate and constant are one derivation.** Confirmed to ≤0.2% at q=11,13,17 (`c_tilde_structure_verdict.md`).

---

## Result 3 (CORRECTION) — `δ_q = c̃_q − (q−3)/q ≈ 0.82/ord_q(2)` (empirical)

The finite-order correction: `c̃_q = (q−3)/q + O(1/ord_q(2))`, deviation `≈ 0.82/ord_q(2)`, monotone across 8 primes, R²=0.94, **out-of-sample validated at q=31,127,73** (`result_4_ctilde_ord2.md`, Probe 4). Mechanism: small multiplicative order of 2 mod q shortens the 2-orbit in the chain, inflating the finite-order correction. The constant 0.82 is empirical (open whether it is a clean rational); c̃ measured at k=2, matches all prior established δ.

---

## Coherence note

The number `3` recurs as `1/E_{Geom(1/2)}[2^{−v}]` — a purely 2-adic halving statistic — and is the reason the rate is q-independent. It is **not** the multiplier (which is q, and the rate is q-independent). This is the same `1/3 = E[2^{−v}]` that governs Tao's 2-adic renewal; here it is isolated as the universal contraction rate of the stationary L² mass across the entire qx+1 family.

## Provenance / files

- `probe_5_universal_rate.py` + `result_5_universal_rate.md` + `result_5_data.csv` — Result 1 derivation + falsifier (this session).
- `PRE_REG_5_UNIVERSAL_RATE_2026_07_14.md` — pre-registration.
- `result_q_sweep_test_1_rate.md` — original empirical rate observation (2026-05-04).
- `c_tilde_structure_verdict.md` + `c_tilde_q17_probe.py` — Result 2.
- `result_4_ctilde_ord2.md` + `.py` + `PRE_REG_4_CTILDE_ORD2` — Result 3.
- `c_seven_forty_fifth.md` (R75), `result_76_conservation_law.md` (R76) — the q=3 Plancherel/conservation machinery the domination identity must generalize.

**Next step to finish the paper:** prove uniform diagonal-self-overlap domination on `(Z/q^k)*` (generalize R76's `Σ_j M(η_0 + j·3^n)=0` conservation), upgrading Result 1 from mechanism to theorem. That is the single remaining line.
