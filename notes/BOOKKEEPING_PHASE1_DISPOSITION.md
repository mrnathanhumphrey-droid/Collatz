# Bookkeeping Phase 1 — Disposition

**Date:** 2026-05-11
**Pre-registration:** `BOOKKEEPING_PHASE1_PRE_REGISTRATION.md` (committed first)
**Phase 1a/b/c:** `TAO_PROOF_CONSTANT_MAP.md`, `TAO_BOOKKEEPING_TRACTABILITY.md`, `TAO_CA_LOOSENESS_PROJECTION.md`

---

## Disposition (one of four, no hedging): **INFEASIBLE**

The Tao §7.2–§7.4 bookkeeping route to an effective C_A for Nisoli closure is closed. Phase 2 is **not justified.**

---

## One-paragraph rationale

The proof structure is line-by-line trackable (Phase 1a produced a clean 42-entry constant map; OCR was usable; no fundamental ambiguities surfaced), and project expertise is sufficient to extract every named-unspecified constant — there are no BLOCKED-by-expertise entries (Phase 1b: 22 TRIVIAL / 13 MODERATE / 5 HARD / **0 BLOCKED**, with the 5 HARDs all in standard 2D local-CLT / large-deviations technique). Phase 1b therefore lands positive on two of the three pre-registered override checks. The third — Phase 1c looseness — lands unambiguously NEGATIVE. The iterated-cubic recursion in Case 3 of §7.4 (the `p_{i+1} ≤ 40A(1+p_i)³ + O(A)` iterated R = A²/ε times) forces P(A, ε) — and hence C_A = O(P^A) — to grow at least like exp(exp(A²)) under faithful bookkeeping, and at the optimistic floor A^{O(A)} under any reading. **At K = 10, no value of A produces a C_A small enough to satisfy Nisoli η = ε_K · M_3 < 1 with M_3 ≈ 800–1000** (even under the optimistic A^{20A} floor, requirement fails by 40+ orders of magnitude). The obstruction is structural to Tao's renewal-process method, not to bookkeeping; redoing §7 with Tao's method cannot produce a polynomial-in-A bound. Per pre-reg §3.4, INFEASIBLE applies whenever "projected C_A too loose at any K" — which is exactly what Phase 1c established. Pre-reg was NULL-favored; this disposition is consistent with the prior.

---

## Three pre-registered override checks — outcomes

| # | Check | Outcome | Where |
|---|---|---|---|
| 1 | Constant map clean (Phase 1a complete, no fundamental ambiguities) | ✅ POSITIVE | Phase 1a: 42 entries, no proof-structure obstruction |
| 2 | Tractability classification mostly TRIVIAL/MODERATE, HARDs have identifiable extraction paths, no BLOCKED-by-expertise | ✅ POSITIVE | Phase 1b: 22 TRIVIAL / 13 MOD / 5 HARD / 0 BLOCKED |
| 3 | Projected C_A satisfies Nisoli η < 1 at K ≤ 10 even under looseness pessimism | ❌ NEGATIVE | Phase 1c: fails at every K, every A, even under optimistic projection |

**Override requires THREE positive. Only TWO positive. Disposition stays NULL-aligned → INFEASIBLE.**

---

## Why INFEASIBLE rather than TRACTABLE_BUT_LOOSE

Pre-reg §3.4 distinguishes:

> **TRACTABLE_BUT_LOOSE:** 1a/1b clean, but projected C_A wouldn't satisfy η < 1 at verified K range. Phase 2 useful with adjusted expectations.
>
> **INFEASIBLE:** proof structure resists line-by-line tracking, OR projected C_A too loose at any K. Bookkeeping route closes.

Phase 1c established that the looseness is **not just bad at the verified K range; it is bad at any K ≤ 10 and indeed at any finite K** because the growth of C_A in A exceeds the K^A gain regardless of K. There is no Phase 2 with "adjusted expectations" that would salvage a usable number — the bound stays too loose for Nisoli closure no matter how carefully extracted, because the proof's renewal-process iteration is the source.

INFEASIBLE is the correct category. Phase 2 is **closed**, not "conditionally justified" or "shelved with adjustments."

---

## Highest-difficulty constants encountered

The HARDs from Phase 1b, in order of load-bearing-ness:

1. **C-14 — Lemma 2.2 (2D local-CLT for renewal process)** absolute constants. Cited by Tao §7.3 (Lemma 7.7 proof, end paragraph: "From Lemma 7.6 and Lemma 2.2 one has P(v_{[1,k−1]} = (j', s')) ≪ k^{−1} G_{k−1}(c((j', s') − (k−1)(4, 16)))"). The single largest reservoir of named-unspecified constants. Phase 2 would need to reprove Lemma 2.2 with effective constants for the specific (1, Pascal' / Hold) step distribution — a 10-15 page side-task using Bhattacharya & Rao 1976-style multivariate Berry-Esseen.
2. **C-28 — large-deviations rate** `c` in `P(j_{[1,k+P]} ≥ 0.9m) ≤ P · exp(−cm)` (Tao §7.4 Case 3, paragraph after Lemma 7.10 and before Lemma 7.7-application). Inherits structure from C-14 and is set by the gap `0.8 − (1/4)(log 9/log 2) ≈ 0.007`.
3. **C-36 — exponential-tail constant** `c` in Lemma 7.10's `exp(−c A²(1+p))` — same machinery as C-28, applied to a different gap.
4. **C-12 — Gaussian width parameter** `c` in argument of `G_{1+s}(c(j' − s/4))` of Lemma 7.7 (eq 7.48). Also inherits from C-14.
5. **C-40 — terminal Vinogradov accumulation** in `P(F) ≪ 10^{-A-2}` (Tao §7.4, paragraph after Lemma 7.10's union-bound, eq 7.66). The final assembly point where all 18 named-unspecified constants compound; cumulative ≪-drift is the largest single source of looseness.

The **dominant looseness source** is actually NOT in this list — it is the **iterated-cubic recursion in Case 3 of Prop 7.8's proof** (the `p′ ≤ p + 10·4·A·(1+p)³ + 10A/3 + 1` chain, last paragraph of §7.4 around eq 7.66 area). This is **explicit numeric** and so does not count as a HARD constant; the looseness it generates is "built into the proof method" rather than "absorbed into a Vinogradov." That is why Phase 1c is structurally negative even when Phase 1b is structurally positive.

---

## A1 OCR-discrepancy log

Phase 1a §6 lists 7 internal cross-checks; **no A1 discrepancies surfaced** in the load-bearing estimates (Lemma 7.4's strip and separation distance both at (1/10) log(1/ε), Lemma 7.7 first-passage bound, Geom(4) and EHold = (4, 16), Case 2/3 thresholds s ~ m/log² m, R = A²/ε from §7.4).

**However, live WebFetch cross-check against arXiv 1909.03562 was blocked by harness permissions** during Phase 1a. The disposition does not hinge on any of the OCR-mangled lines, but a full Phase-2 reading should re-verify directly against the published PDF before any bookkeeping work commits. This is logged as a Phase 1 limitation, not a finding.

---

## Phase 2 disposition — CLOSED

Phase 2 (actual bookkeeping) is **not justified** under the locked decision rule. Resources should not be spent on it.

**Alternative routes** (NOT decided here — surfaced for project's broader path-planning):

- **R77.3 bypass route** (already on file at `result_77_3_nisoli_bypass.md`): bound ‖T − T_K‖ directly without going through Prop 1.17. Most promising overall direction; does not depend on this Phase 2.
- **Different Fourier-analytic bound on the Syrac MGF** with polynomial-in-A constant: would require novel technique outside Tao's method. Out of scope for any bookkeeping pass.
- **Empirical-only closure** at fixed K via Plancherel + spectrum computation, never quantifying C_A: already in place for c = 7/45 through k = 6 (the |c − S_k/3| ≤ 0.0133 · (1/2)^k envelope). Not a rigorous proof.

The recommendation from this Phase 1 is therefore: **redirect effort to the R77.3 bypass route**, which does not need an effective C_A.

---

## Closure paragraph

R77.2's Stage 2 outcome (δ) — that Tao Prop 1.17 lacks an effective C_A — is upgraded here to: **the missing C_A cannot be recovered by a bookkeeping pass through §7.2–§7.4 to a value useful for Nisoli closure, even with perfect bookkeeping execution.** The bookkeeping route is technically possible (no expertise blocker) but produces a number too loose at any K, because Tao's renewal-process technique encodes a Case-3 iterated-cubic that forces super-exponential A-dependence regardless of constant-extraction quality. The Nisoli framework needs a polynomial-in-A bound on the Syrac Fourier coefficient; Tao's §7 method does not deliver one, and bookkeeping cannot create one. **Bookkeeping route closes; bypass routes remain open.**

---

End of Phase 1 disposition.
