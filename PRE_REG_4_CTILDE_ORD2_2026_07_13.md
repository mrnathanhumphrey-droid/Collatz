# PROBE 4: c̃_q Deviation vs ord_q(2) (PRE_REG)

**Agent:** compute (Agent 1/2/3)
**Repo:** `C:\Collatz`
**Date locked:** 2026-07-13.
**Status:** pre-registration — fire only after reading §1.
**Track D, fully independent of Probes 81 / 1B / 3a — cheapest probe in the plan, can run concurrently.**
**Do not write to any file outside the deliverables in §6.**

---

## 0. The question

`c̃_q := lim_k S_k^{(q)} / (q/3)^k` (the renormalized qx+1 Plancherel-mass limit) fits **`c̃_q = (q−3)/q`** to 1% at q = 11, 13, 17. But two primes deviate:

- q = 5: δ ≈ 0.09 (Aitken suggests genuine, not finite-k).
- q = 7: δ ≈ 0.21 — **worse than q=5**, which kills any "finite-q correction" story (a correction would shrink with q).

The **primitive-root** explanation is already ruled out: q=17 has ord(2) = 8 ≠ q−1 yet fits cleanly (`c_tilde_structure_verdict.md`). This probe tests the next covariate: **ord_q(2)**, the multiplicative order of 2 mod q. q=7 is the **unique** small prime with ord = 3 (2, 4, 1), and it's the outlier — a sample of one you can't replicate within the current set. New primes at controlled orders break that degeneracy.

**Distinct from the c_∞ / Hecke arc.** `c̃_q` (this probe, the S_k/(q/3)^k limit) is NOT `c_∞(q)` (the Tao character moment in `probe_p73_p89_p97_2026_06_01.py` / the black-hole pre-reg). Do not conflate them or reuse those outputs as if they were c̃_q. If ord_q(2) turns out to explain BOTH, that's a bonus for a later probe — out of scope here.

---

## 1. Context (read first)

- `c_tilde_structure_verdict.md` — c̃_q = (q−3)/q confirmed q=11,13,17 within 1%; primitive-root ruled out via q=17.
- `c_tilde_q17_probe.py` — the existing c̃_q computation at q=17. **Extend this to the new primes.**
- `result_q_sweep_test_2_c_q.md` — the q-sweep; literal c_q = S_∞/q falsified, renormalized c̃_q exists universally; q/3 ratio universal.
- `result_qspectrum.md` — K_k^{(q)} spectrum trivial across q (q=7 anomaly is NOT spectral) — rules out a within-level-eigenvalue origin for δ, pointing the mechanism to the ratio/limit structure this probe measures.
- STATE.md Open pieces #10, #13 — the q-sweep follow-up notes ("cheap at k=2; minutes per q").

---

## 2. Object

For each prime q, the deviation from the closed form:

    δ_q := c̃_q − (q − 3)/q,      c̃_q := lim_k S_k^{(q)} / (q/3)^k

against the covariate `ord_q(2)` (multiplicative order of 2 mod q).

**Prime set (chosen so ord_q(2) is controlled and clean):**

| q | ord_q(2) | note |
|---|---|---|
| 7 | 3 | existing outlier (2³ = 8 ≡ 1) |
| 5 | 4 | existing mild deviation |
| 31 | 5 | Mersenne 2⁵−1 |
| 11 | 10 | existing clean fit (baseline) |
| 13 | 12 | existing clean fit (baseline) |
| 17 | 8 | existing clean fit (baseline) |
| 127 | 7 | Mersenne 2⁷−1 |
| 73 | 9 | 73 \| 2⁹−1 = 511 = 7·73 |

New compute: q ∈ {31, 127, 73} (ord 5, 7, 9). Recompute {5, 7, 11, 13, 17} for a consistent δ under one pipeline.

---

## 3. Hypotheses (pre-registered, mutually exclusive)

- **H_ORD** — |δ_q| is a decreasing function of ord_q(2) (small order → large deviation). Sharp prediction: q=31 (ord 5) has |δ| between q=7 (ord 3, large) and the clean fits (ord ≥ 8, |δ| ≲ 0.01); q=73 (ord 9) and q=127 (ord 7) fit (q−3)/q to ≈ 1%. Confirmed iff δ is monotone in ord and the new primes land as ordered.
- **H_ORD_OVER_Q** — the clean covariate is `ord_q(2)/q` (or `1/ord_q(2)`), not ord alone; regress and report which functional form linearizes δ.
- **H_NULL** — δ_q shows no clean dependence on ord_q(2); mechanism is elsewhere. (q=31, 127, 73 would then scatter regardless of order.)
- **H_FINITE_K** — apparent δ at the new primes is finite-k transient (Aitken/Richardson shrinks it toward 0). Re-test the limit before attributing δ to ord.

**Most-likely outcome:** H_ORD or H_ORD_OVER_Q — small multiplicative order of 2 collapses the 2-orbit structure inside the (Z/q)* Plancherel sum, and q=7's ord=3 is the extreme case. State the prior, don't let it relax §4.

---

## 4. Method

1. **Extend `c_tilde_q17_probe.py`** to q ∈ {31, 73, 127} and recompute {5, 7, 11, 13, 17} under the same pipeline. Compute S_k^{(q)} to the largest k the state space allows.
2. **Extrapolate the limit properly** — c̃_q = lim S_k/(q/3)^k via Aitken/Richardson, NOT a single k. Report c̃_q with an uncertainty and the k-range used. (STATE flags q=5's δ as Aitken-confirmed genuine; hold the new primes to the same standard.)
3. Compute δ_q = c̃_q − (q−3)/q and ord_q(2) (verify: 2^{ord} ≡ 1 mod q, minimal).
4. **Regress** δ_q on ord_q(2), on ord_q(2)/q, and on 1/ord_q(2). Report which linearizes; report R² and the residual for q=7 (the anchor outlier) under each.
5. Report the full (q, ord_q(2), c̃_q, (q−3)/q, δ_q) table sorted by ord.

**Compute note / bug guards:**
- State space scales with q^k (or the qx+1 analog used in `c_tilde_q17_probe.py`). Cheap at k=2–3; q=127 at high k is the heaviest cell — **size k to the extrapolation need, and report the k reached per q.** Do not silently use a smaller k for q=127 and compare its under-converged c̃_q against well-converged small-q values — that manufactures a fake δ. Match extrapolation quality across primes or flag the mismatch.
- ord_q(2) must be the **minimal** period (e.g. verify 2⁹ ≡ 1 mod 73 AND no smaller divisor of 9 works).
- Do NOT pull c_∞(q) values from the Hecke/black-hole arc as a substitute for c̃_q (§0).
- q=3 is a separate regime (c̃_3 = 7/15; the (q−3)/q form vanishes there) — exclude from the δ regression.

---

## 5. What NOT to do

- Do not attempt to derive (q−3)/q here (that's a separate theorem candidate) — this probe measures the *deviation's covariate*.
- Do not reuse or overwrite the c_∞ / Hecke arc outputs.
- Do not touch `lagarias_framework_synthesis.docx` or any external-facing file.

---

## 6. Deliverables

- `result_4_ctilde_ord2.py` — the probe (extends `c_tilde_q17_probe.py`).
- `result_4_ctilde_ord2.csv` — per q: ord_q(2), k-range, c̃_q (+ uncertainty), (q−3)/q, δ_q.
- `result_4_ctilde_ord2.md` — disposition: which hypothesis fired, the δ-vs-ord table sorted by ord, the regression form that linearizes (with R² and q=7 residual), and a one-line routing note (does ord_q(2) explain the q=5/q=7 anomalies, upgrading (q−3)/q toward a mechanism).
- `result_4_log.txt` — run log incl. per-q k reached and extrapolation method.

Append the disposition to `STATE.md` under a new dated entry. Do not rewrite existing STATE.md content.

---

## 7. Reporting discipline

The whole point is q=7 being worse than q=5 — if the new primes DON'T order by ord_q(2), that's H_NULL and it's a real result that redirects the mechanism hunt; report it straight. Do not declare H_ORD on 3 new points if q=127 or q=73 lands off-trend — call the miss. Match extrapolation quality across primes before comparing δ's, or the covariate story is confounded by k.
