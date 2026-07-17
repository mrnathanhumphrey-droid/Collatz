# PRE_REG: Probe 5 — Derivation of the Universal Rate S_k^{(q)} ~ (q/3)^k

**Agent:** the qx+1 / c̃_q paper agent (arithmetic side).
**Repo:** `C:\Collatz`
**Status:** pre-registration. Read §1 before firing. **§4 Task A (falsifier) is mandatory and comes first — before any proof.**
**Priority:** HIGH for the standalone paper; this is the one missing leg. Independent of every Collatz-closure thread.
**Filesystem:** read-only except §7 deliverables. Append to STATE.md only, one entry.

---

## 1. Context — the standalone paper, and the one gap

There is a complete three-result paper here about exponential sums / stationary measures on `(Z/q^k)*` for the qx+1 map `x → (qx+1)/2^v`. **Collatz need not be in the title.** The three results:

| # | result | status | source |
|---|---|---|---|
| **Rate** | `S_k^{(q)} ~ (q/3)^k`, universal in q | **empirical only — no derivation** | `result_q_sweep_test_1_rate.md`, STATE:197/273 |
| Constant | `c̃_q := lim S_k^{(q)}/(q/3)^k = (q−3)/q` | empirical, ≤0.2% at q=11,13,17 | `c_tilde_structure_verdict.md` |
| Correction | `δ_q = c̃_q − (q−3)/q ≈ 0.82/ord_q(2)` | empirical, R²=0.94, OOS at q=31,127,73 | `result_4_ctilde_ord2.md` (Probe 4) |

**Definitions (do not re-derive; read the sources).** `S_k^{(q)} = q^k · Σ_ξ π_k(ξ)²`, where `π_k` is the stationary distribution of Tao's qx+1 Syracuse chain on `(Z/q^k)*`. So `S_k = q^k · ‖π_k‖²`, and the rate claim `S_k ~ (q/3)^k` is **exactly** the statement

> **‖π_k‖² ~ C_q · 3^{−k}, with the rate `3^{−k}` INDEPENDENT of q.**

Read first: `result_q_sweep_test_1_rate.md` (the rate fact + exact S_k rationals through k=5, q∈{3,5,7,11,13}), `c_seven_forty_fifth.md` / `result_76_conservation_law.md` (the R75/R76 Plancherel + conservation machinery for q=3), `result_4_ctilde_ord2.md` (the correction pillar; **confirm it exists and says δ≈0.82/ord — it does, verified 2026-07-14**), and `c_tilde_structure_verdict.md`.

**This is the LEADING rate, not the subdominant ε_k swamp.** `‖π_k‖² ~ 3^{−k}` is the top-order L² decay — a clean transfer-operator contraction. It is a *different object* from the c=7/45 subdominant resonance (ε_k = S_k^{(3)} − 7/15, the 0.984 / period-9 / branch-cut saga that consumed two months). Do not import that difficulty here; the leading rate is a separate, lower tier. **This is why it should be derivable in ~a page.**

---

## 2. The claim under test (a DERIVATION target, not a measurement)

> **Theorem (target).** For the qx+1 Syracuse chain on `(Z/q^k)*` with `v ~ Geom(1/2)` halving, the stationary L² mass satisfies `‖π_k‖² = C_q · 3^{−k}·(1 + o(1))`, where the exponential rate `3^{−k}` is independent of q. Equivalently `S_k^{(q)}/S_{k−1}^{(q)} → q/3` for every q.

**The mechanism to make rigorous (separation of variables).** `S_k = q^k · ‖π_k‖²`. The `q^k` is the state-count / Plancherel normalization on `(Z/q^k)*`. The `3^{−k}` per-level contraction of `‖π_k‖²` is set by the **Geom(1/2) halving against the 3-adic modulus structure of the map** — and q enters the *operator* only through the multiplicative character, which affects the leading **constant** `C_q`, not the **rate**. The paper's own data already exhibits this separation: **q=7's anomaly lives entirely in the constant** (`c̃_7 ≈ 0.78` vs `(q−3)/q = 4/7`) while **its rate is still clean q/3**. The derivation must formalize an *observed* separation, not conjure one.

**Where the "3" must come from.** The target of the proof is to identify precisely why the contraction is `1/3` and not `1/q`, `1/2`, or `1/4`. State this explicitly; a proof that yields the right rate without saying where the 3 enters is incomplete. Candidate sources to test in the derivation: the `3n+1`→`3`-in-`qx+1` is a red herring (the rate is q-independent, so it is NOT the multiplier q); the `1/3` most likely traces to the **fixed 3-adic / mod-3 fiber structure of the halving map** shared across all q, or to a `(1 − 2·(1/2))`-type mass-conservation identity generalizing R76. Pin it.

---

## 3. Outcomes (pre-registered, mutually exclusive)

- **H_PROVED** — the separation goes through: `‖π_k‖²` decay factorizes into a q-normalization × a q-independent `3^{−k}` contraction, with the "3" identified. The rate becomes a theorem. Paper is complete.
- **H_BREAKS** — the derivation localizes the **exact step** where q re-enters the rate (i.e., the contraction is only *approximately* q-independent, or is q-independent only under a condition on ord_q(2) / q mod something). This is a real result: it says the universal rate is emergent/asymptotic, not exact, and names the coupling term. Report the step, do not paper over it.
- **H_FALSIFIED** — Task A finds a q where the empirical rate ≠ q/3. Then the "universal rate" is false as stated and the paper's first pillar is wrong. **This is the most valuable negative** and is why Task A runs first.

**Pre-registered most-likely:** the *fact* survives Task A (the rate is empirically rock-solid across 6 q including the q=7 anomaly — unlike this arc's dead structural priors, which were guesses; this one is a measured invariant). For the *derivation*, weight **H_PROVED at leading order but flag realistic H_BREAKS at one lemma** — the honest expectation is a proof that goes through modulo identifying the "3", with a real chance the q-independence is asymptotic rather than exact and needs a stated error term. Do not let the paper's momentum inflate a partial derivation into a clean theorem.

---

## 4. Tasks — in order, no reordering

### Task A — the falsifier FIRST (mandatory, before any proof)

**Prove-it-or-break-it demands you attack it before you defend it.** Sweep the empirical rate `S_k^{(q)}/S_{k−1}^{(q)} → q/3` on **adversarial q chosen where the separation is most likely to fail**, not on the comfortable q∈{5,7,11,13} that already confirm it:

1. **Small / pathological `ord_q(2)`** — the q=7 family. q with `ord_q(2)` unusually small relative to q (the case where the character/halving coupling is most resonant). If any breaks q/3, that's the coupling entering the rate.
2. **Composite and even q** — q=9, 15, 25, 21, and even q if the chain is well-defined there. The `q^k` state-count / `(Z/q^k)*` framing is cleanest for prime q; composite/even q is where "q^k states" and "3-adic modulus" assumptions can crack. Report where the construction itself breaks (that is a scope boundary, a finding).
3. **q ≡ 0 mod 3** (q=3, 9, 15, …) and q sharing structure with 2 — the degenerate slices. q=3 is the critical case (q/3=1); test whether the framing degenerates predictably.

**Method:** exact-rational S_k where feasible (reuse the `c_tilde_q17_probe.py` / `result_4_ctilde_ord2.py` exact-stationary machinery — do NOT build a second path), float sparse power-iteration with an exact cross-check at k=2 otherwise. Compute the ratio `S_k/S_{k−1}` and its approach to q/3 to as many k as the chain size allows. **Report the ratio, not a binary.** A q converging to `q/3 ± o(1)` confirms; a q converging to anything else falsifies.

**Deliverable of Task A is a table of `(q, ord_q(2), S_k/S_{k−1} → ?, verdict)`. Print it. If any q breaks q/3, STOP and report H_FALSIFIED — the derivation target is wrong and that is the finding.**

### Task B — the derivation (only after Task A confirms the fact)

Set up `S_k = q^k‖π_k‖²` via the R75/R76 Plancherel + conservation framework, generalized from q=3 to general q. Then attempt the separation:

1. Write the one-level transfer operator for `π_k → π_{k+1}` (or for the defect from uniform) on `(Z/q^k)* → (Z/q^{k+1})*`.
2. Show its action factorizes: a q-dependent normalization (state-count / character mass) × a **q-independent contraction** on the L² norm.
3. **Identify the "3".** Derive the `1/3` explicitly — from the halving-against-3-adic-fiber structure, or a mass-conservation identity generalizing `Σ_j M(η_0 + j·3^n) = 0` (R76). Name it.
4. If step 2 or 3 does not close, **localize the exact term where q re-enters the rate** and report it as H_BREAKS. A one-line "and here q leaks into the exponent" is worth more than a hand-wave.

**Decision rule for H_PROVED:** the rate `3^{−k}` must be derived q-independent **with the source of the 3 named**, not fitted or asserted. A derivation valid only for prime q, or only under a condition on ord_q(2), is **H_BREAKS with that condition stated**, not H_PROVED.

### Task C — only if H_PROVED: does the same argument yield the constant?

If the rate falls out, the leading constant of the same factorization *is* `c̃_q`. Check whether the argument delivers `c̃_q = (q−3)/q` (pillar 2) as the leading coefficient — that would unify pillars 1+2 into one derivation and is the paper's spine. The `(q−3) = q − 3` again puts a bare **3** against the multiplier; if the derivation shows why, note it. Do not fit; derive or report open. (The `0.82/ord` correction, pillar 3, stays empirical — out of scope for derivation here.)

---

## 5. Guards

- **Falsifier before proof** (§4 Task A). Attack adversarial q first. A proof built only on the confirming q is the fit-inside-the-window failure mode that killed ⌊r/2⌋+2 and H_EMPTY — in derivation clothes.
- **Leading rate ≠ subdominant ε_k.** `‖π_k‖² ~ 3^{−k}` is top-order and clean; do NOT drag in the c=7/45 / 0.984 / period-9 machinery. If you find yourself needing it, you've mis-scoped the object.
- **Name the 3.** A rate-right proof that doesn't say where `1/3` comes from is H_BREAKS, not H_PROVED.
- **Do not conflate `c̃_q` with `c_∞(q)`** (the Tao character moment / Hecke arc in the p73/p89/p97 black-hole probes). Different object; do not reuse those outputs. (Per PRE_REG_4 §warning.)
- **Reuse the exact-stationary machinery**, don't build a parallel one (`check_existing_convention`).
- **Exact arithmetic for the falsifier** where feasible; float only with a k=2 exact cross-check.
- **Do not touch** THEOREM_C_745, Thms 78.1–78.3, R81b, or the ε_k / subdominant thread. None are at stake; the rate is independent of all of them. State so.

---

## 6. What a null / break looks like, and why each is fine

1. **H_FALSIFIED** (Task A finds a q ≠ q/3) — the biggest result: the universal rate is false, pillar 1 dies, the paper reshapes around "rate q/3 for prime q with condition X." Cheap and decisive; that's why it runs first.
2. **H_BREAKS** (derivation localizes where q enters the rate) — a real refinement: the universal rate is asymptotic/conditional, with a named coupling term. Often *more* interesting than a clean proof, and it tells the paper exactly what to claim.
3. **H_PROVED** — the paper's first leg becomes a theorem and (via Task C) possibly unifies with the second. The target outcome; do not force it.

A partial derivation reported honestly as H_BREAKS is worth more than an H_PROVED that skipped the adversarial q or hand-waved the 3.

---

## 7. Deliverables

- `probe_5_universal_rate.py` — the Task A adversarial-q falsifier (exact/float S_k ratios) + any numerical checks supporting Task B.
- `result_5_universal_rate.md` — disposition. Must contain: the Task A adversarial-q table with per-q `S_k/S_{k−1} → ?` and verdict; the derivation writeup (Task B) OR the exact localized break-point; the source of the "3" if H_PROVED; the Task C constant check if reached; explicit statement that the ε_k thread / C_745 / 78.x are untouched.
- `result_5_data.csv` — per (q, k): S_k (rational or float), ratio, ord_q(2), verdict.
- `result_5_log.txt`

One dated STATE.md entry. If an intermediate verdict is printed and revised, remove the stale one.

---

## 8. Reporting discipline

This is a **derivation** probe: the target is a proof, and the falsifier is a q where the rate isn't q/3. Report the outcome that fired. A proof that works only on the confirming q is not a proof — it is the window-fit failure mode. If the derivation stalls, **name the step**; "open at the following identity" is a result. If the rate is only asymptotically q-independent, say so with the error term; do not round an approximate separation up to an exact theorem. The paper is stronger with an honest H_BREAKS than with an oversold H_PROVED.

---

## User additions
