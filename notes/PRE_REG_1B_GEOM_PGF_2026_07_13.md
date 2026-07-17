# PROBE 1B: Two-Root-of-Unity Audit of the Halving PGF (PRE_REG)

**Agent:** compute (Agent 1/2/3)
**Repo:** `C:\Collatz`
**Date locked:** 2026-07-13.
**Status:** pre-registration — fire only after reading §1.
**Fires BEFORE Probe 81** (cheaper; audits the headline denominator; sharpens Phase 2's target either way).
**Do not write to any file outside the deliverables in §6.**

---

## 0. What this audits and why it's first

The §CONJ block of `PRE_REG_81_FHAT_PHASE_2026_07_13.md` proposes a fully mechanistic origin for the constant:

    7/45  =  N(2 − ω) / ( 3² · (1 + 4) )

built entirely on the Geom(½) halving-weight generating function `G(z) = Σ_{v≥1} 2^{−v} z^v = z/(2−z)` evaluated at two roots of unity:

- `G(−1) = −1/3` → class mass (0, 1/3, 2/3) → squared 1:4 → **the 9**.
- `|G(ω)|² = 1/N(2−ω) = 1/7` → **the 7**.

**But v is NOT exactly Geom(½).** R68 (`lagarias_sinai_validation.md`) measured 0.5%–25% deviation at specific j; the unconditional ensemble mean is `v_2 = 2.102`, not the Geom(½) value 2.0 (`result_density_one_v2_bounds.md`). If the true v-measure moves `G(−1)` off −1/3 or `|G(ω)|²` off 1/7, the mechanistic story is idealized. This probe measures exactly that, at both constants at once — same computation, two evaluation points.

**This probe does NOT re-derive 7/45.** The rigorous value comes from exact-rational Plancherel (R75/R76), which computes the true stationary measure and does not visibly substitute Geom(½). What is at stake here is (a) whether the §CONJ *mechanism* reflects the true measure or only the idealization, and (b) a flag for a separate proof-audit of whether any step of `THEOREM_C_745.md` leans on Geom(½). Keep (b) as a motivated hypothesis, not an assumed conclusion.

---

## 1. Context (read first)

- `PRE_REG_81_FHAT_PHASE_2026_07_13.md` §CONJ — the conjecture being audited.
- `lagarias_sinai_validation.py` / `lagarias_sinai_validation.md` — R68 v ~ Geom(½) deviation measurements (0.5–25% at specific j). Existing tooling — extend, don't rebuild.
- `result_density_one_v2_bounds.md` — ensemble mean v_2 = 2.102; the TEST-B tautology caveat (do not repeat that error).
- `data/v_seq_N8388608.parquet` — 2,796,202 odd, coprime-to-3 starts in [3, 8388607], per-trajectory v-sequences. Primary data.
- R68 refined finding (STATE.md Open piece #6): `v_t` given `m_t mod 2^k` is arithmetic-deterministic; the Geom-like marginal arises because the trajectory measure is non-uniform mod 2^k for k ≥ 3. **This means "the v-distribution" is ambiguous — measure the marginal actually seen by the chain, and report which marginal was used.**

---

## 2. Object

Empirical PMF `P̂(v)` from the measured v-values. The empirical halving PGF:

    Ĝ(z) := Σ_{v≥1} P̂(v) · z^v

evaluated at two roots of unity:

    Ĝ(−1) = E[(−1)^v]      (target: −1/3)
    Ĝ(ω)  = E[ω^v],  |Ĝ(ω)|²   (target: 1/7),   ω = e^{2πi/3}

Also report the Geom(½) baseline `G(−1) = −1/3`, `|G(ω)|² = 1/7` alongside, and the deviation `Δ_{−1} = Ĝ(−1) − (−1/3)`, `Δ_ω = |Ĝ(ω)|² − 1/7`.

---

## 3. Hypotheses (pre-registered, mutually exclusive)

- **H_ROBUST** — both constants survive the true measure within tolerance (§4 decision rule) across all v-brackets. ⇒ 7/45 has a fully mechanistic, measure-robust origin; §CONJ upgrades from conjecture toward a **robustness theorem**, and Phase 2 gets a sharp target (derive R75/R76's 7 as N(2−ω)).
- **H_DRIFT** — one or both constants move outside tolerance. ⇒ the §CONJ generating-function mechanism is Geom(½)-idealized; the "7 = N(2−ω)" story is approximate, not the exact origin of the rigorous 7/45. Motivates the separate proof-audit of `THEOREM_C_745.md` for a hidden Geom(½) substitution (item (b), §0).
- **H_BRACKET_SPLIT** — constants hold in one v-regime and drift in another (e.g. hold at large v, drift in the low-v descent funnel, or vice versa). This is the outcome the retired Bohr signal warns about — do NOT aggregate it away.

---

## 4. Method

1. **Build P̂(v)** from `data/v_seq_N8388608.parquet`. Report N, support range, and the mean (expect ≈ 2.102 as a smoke check against `result_density_one_v2_bounds.md`; if it disagrees, STOP and reconcile).
2. **Evaluate** Ĝ(−1) and Ĝ(ω); report Ĝ(−1), Ĝ(ω) (complex), |Ĝ(ω)|², and the two deviations Δ.
3. **Bracket-stratify** — MANDATORY. Split by v-magnitude/trajectory-position brackets exactly as `result_bohr_probe_strat.md` did (the aggregate-vs-descent-funnel distinction that deflated the Bohr signal). Report Ĝ(−1), |Ĝ(ω)|² **per bracket**, not only pooled. A pooled pass masking a bracket split is H_BRACKET_SPLIT, not H_ROBUST.
4. **Marginal disclosure** — state which v-marginal was used (unconditional per-step? conditional on m mod 2^k? per-j?). Per R68's refined finding the marginal is measure-dependent; if cheap, report the two most natural marginals and whether the constants are marginal-robust.
5. **Decision rule / tolerance.** The Geom(½) idealization is exact; the question is how far the true measure moves the constants relative to the ~5% mean deviation already known.
   - **H_ROBUST fires** iff |Δ_{−1}| ≤ 0.01 AND |Δ_ω| ≤ 0.01 in **every** bracket.
   - Any single bracket outside → H_DRIFT (if pooled also fails) or H_BRACKET_SPLIT (if pooled passes but a bracket fails).
   - Report the raw deviations regardless; do not soften a 0.03 drift into "approximately holds."

**Bug guards:**
- Do NOT repeat the `result_density_one_v2_bounds.md` TEST-B tautology (mean_v > log₂3 restates "all trajectories reached 1"). This probe evaluates a PGF at roots of unity — a genuinely different object — but keep the caveat in view when interpreting.
- Do NOT aggregate across brackets before reporting (step 3 is load-bearing).
- ω is exact `e^{2πi/3}`; use exact `(-1 + i√3)/2`, not a rounded decimal, so |Ĝ(ω)|² isn't polluted at the 3rd digit.

---

## 5. What NOT to do

- Do not re-derive 7/45 (that's R75/R76, exact and separate).
- Do not audit `THEOREM_C_745.md`'s proof here — H_DRIFT only *motivates* that; it's a separate task.
- Do not touch `lagarias_framework_synthesis.docx` or any external-facing file.
- Do not fire Probe 81's phase machinery here — 1B is v-distribution only.

---

## 6. Deliverables

- `result_1b_halving_pgf.py` — the probe (extends `lagarias_sinai_validation.py`).
- `result_1b_pgf_data.csv` — per bracket (and pooled): N, mean v, Ĝ(−1), Re/Im Ĝ(ω), |Ĝ(ω)|², Δ_{−1}, Δ_ω.
- `result_1b_halving_pgf.md` — disposition: which hypothesis fired, the per-bracket table, the marginal used, and a one-paragraph routing statement — if H_ROBUST, hand Phase 2 the "derive R75/R76's 7 as N(2−ω)" target; if H_DRIFT/H_BRACKET_SPLIT, flag the `THEOREM_C_745.md` proof-audit and state how far the mechanism is idealized.
- `result_1b_log.txt` — run log incl. data source and exact ω used.

Append the disposition to `STATE.md` under a new dated entry. Do not rewrite existing STATE.md content.

---

## 7. Reporting discipline

Both outcomes are results. H_ROBUST → a robustness theorem and a mechanistic 7/45. H_DRIFT → you found a Geom(½) idealization in the spine's headline story before a referee did. Report the deviations as measured; a near-miss on the 100·Δ is a drift, not a "hold." Do not let the elegance of `7/45 = N(2−ω)/(3²·5)` relax the tolerance.
