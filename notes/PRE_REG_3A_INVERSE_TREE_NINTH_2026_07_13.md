# PROBE 3a: The 1/9 Inverse-Tree Decay (PRE_REG)

**Agent:** compute (Agent 1/2/3)
**Repo:** `C:\Collatz`
**Date locked:** 2026-07-13.
**Status:** pre-registration — fire only after reading §1.
**Track C (sibling), independent of Probes 81 / 1B — can run concurrently.**
**Do not write to any file outside the deliverables in §6.**

---

## 0. What this nails, and why it's on the spine

`duality_S_vs_D_verdict.md` reports that the 3x+1 single-basin inverse-tree Plancherel mass has consecutive-depth ratio `D_{n+1}(k)/D_n(k) → 1/9`, but marks it **"n=6 too small to confirm cleanly."** This probe extends the tree and settles whether the limit is **exactly 1/9**.

Two reasons it's worth a clean confirmation, not just a footnote:

1. **1/9 = 3⁻² is another spine entry.** D_n is an L² (Plancherel) mass, so a 1/9 decay is the *squared* version of a 1/3 mass-branching — the same 3² that appears as squared-class-mass in 7/45 = 7/(9·5) and as the mod-9 darkness onset. If it's exactly 1/9, the sibling thread joins the mod-9 spine.
2. **It anchors 3b** (the standalone result): *inverse-tree Plancherel mass is a basin/cycle-count detector*. 3x+1 (single conjectured basin) decays at 1/9; 3x−1 (three cycles {1,2}, {5,7,10,14}, {17..34}) stays elevated (`duality_S_vs_D_verdict.md` Test 2: ratio diverges to 10³–10⁴× by depth 6). Pinning the single-basin rate is step one of that statement.

**Scope of 3a: the single-basin 3x+1 decay rate only.** The cycle-count-detector generalization (3b) and the write-up (3c) are downstream and out of scope here.

---

## 1. Context (read first)

- `duality_S_vs_D_verdict.md` — the verdict quoting D_{n+1}/D_n → 1/9 (Agent 2, 3x+1) and the basin-fingerprint divergence (Agent 3, 3x−1). Test 1 table has the current n≤6 ratios.
- `duality_followup_verdict.md`, `duality_followup_check.py`, `duality_followup_data.csv` — matched-N control; ~95% of the raw sibling gap is sample-size, residual factor 0.2–4.
- `result_inverse_tree_residue.md` — Agent 2's 3x+1 single-basin D_n(k) tables, exact rationals, no value truncation.
- `agent3_inverse_tree_3xm1_Dn.py` / `agent3_inverse_tree_3xm1_Dn_v2.py` — Agent 3's 3x−1 three-basin builder (N_MAX = 10⁸ cap). For 3b later, not needed for 3a's single-basin rate.
- `sibling_3x_minus_1_symmetry_verdict.md` — forward K₋ = σK₊σ (S_k identical); the paradox 3a/3b resolve.

---

## 2. Object

For the **3x+1 single-basin inverse tree** (root 1, backwards map, exact rationals, no value truncation), the depth-n residue-level-k Plancherel mass `D_n(k)`, and the consecutive-depth ratio:

    ρ_n(k) := D_{n+1}(k) / D_n(k)

Target: `ρ_n(k) → 1/9` as n → ∞, for each fixed k ≥ 2 (k=1 is the degenerate fixed-point D_n(1) = 2/9, exclude from the limit test — it's already constant).

---

## 3. Hypotheses (pre-registered, mutually exclusive)

- **H_EXACT_NINTH** — `ρ_n(k) → 1/9` and the limit is **exactly** 1/9 (rational, from the 3-adic branching of the inverse map squared into the L² mass). Confirmed iff the exact-rational ρ_n(k) converges to 1/9 monotonically or with decaying oscillation, and the extrapolated limit matches 1/9 to the achieved precision.
- **H_OTHER_RATIO** — ρ_n(k) converges to a clean constant ≠ 1/9 (report it as an exact rational if the D_n are exact).
- **H_K_DEPENDENT** — the limit depends on k (would refute a single universal decay rate; interesting in its own right).
- **H_NONCONVERGENT** — no clean limit at reachable depth (report the trajectory and the obstruction).

**Most-likely outcome (state, then test honestly):** H_EXACT_NINTH, on the 3²-from-squared-branching argument. Do not let that relax the convergence bar in §4.

---

## 4. Method

1. **Extend the 3x+1 single-basin inverse tree** beyond n=6 using Agent 2's exact-rational construction (`result_inverse_tree_residue.md` method). Push n as deep as exact arithmetic allows within a sane compute budget; report the deepest n reached and why it stopped (vertex-count blowup, memory, time).
2. Compute `D_n(k)` as exact rationals for k = 2..k_max, n = 0..n_max. Report the tables.
3. Compute `ρ_n(k)` for each k; tabulate the sequence in n. Apply Aitken/Richardson extrapolation to estimate lim_n ρ_n(k) and its uncertainty.
4. **Exact-limit test.** Since D_n are exact rationals, test whether `ρ_n(k) − 1/9` → 0 with a clean rate. If the D_n satisfy an exact recurrence in n, derive the limit algebraically (that upgrades H_EXACT_NINTH from numerical to proved) — attempt this before falling back to extrapolation.
5. **k-coherence.** Report lim ρ_n(k) across k. A single k-independent 1/9 is the clean result; a k-dependent limit is H_K_DEPENDENT.

**Decision rule:** H_EXACT_NINTH fires iff the extrapolated lim_n ρ_n(k) = 1/9 within the achieved precision for **every** tested k ≥ 2, AND (strong form) an exact recurrence yields 1/9 algebraically. Numerical-only convergence to 1/9 fires the weak form; state which.

**Bug guards:**
- **k=1 excluded** from the limit test (degenerate fixed point D_n(1) = 2/9).
- Single-basin only — do NOT mix in Agent 3's three-basin 3x−1 masses here (that conflation is the basin fingerprint, which is 3b's subject, and would corrupt the single-basin rate).
- Use exact rationals as far as depth allows; only switch to float with an explicit precision note if forced, and never fit the limit from a single n.
- Matched-N is NOT relevant to 3a (that control was for the *sibling comparison*, not the single-basin rate) — do not truncate the 3x+1 tree to any external N.

---

## 5. What NOT to do

- Do not build 3b (cycle-count detector) or write 3c here — 3a is the single-basin rate only.
- Do not re-run the forward symmetry (proved) or the matched-N sibling control (done).
- Do not touch `lagarias_framework_synthesis.docx` or any external-facing file.

---

## 6. Deliverables

- `result_3a_inverse_tree_ninth.py` — the probe (extends Agent 2's builder).
- `result_3a_dn_tables.csv` — D_n(k) exact rationals (num/den) and ρ_n(k), per (n, k).
- `result_3a_inverse_tree_ninth.md` — disposition: which hypothesis fired, the ρ_n(k) → limit table with extrapolation, the exact-recurrence derivation if found (weak vs strong form stated), k-coherence, and a one-line routing note to 3b (does 1/9 confirmed motivate the basin-detector generalization).
- `result_3a_log.txt` — run log incl. deepest n and stop reason.

Append the disposition to `STATE.md` under a new dated entry. Do not rewrite existing STATE.md content.

---

## 7. Reporting discipline

If the limit is 1/9 numerically but no exact recurrence is found, say so — that's the weak form, not a theorem. If it's an ugly rational near 0.111 that is NOT 1/9, report the ugly rational; do not round it to 1/9 because the spine wants a 3². A k-dependent limit refutes the universal rate and is a finding, not a failure.
