# DWM_DISPOSITION — Davies-Wiseman-Milburn quantum trajectory vs Syracuse

**Date:** 2026-05-15
**Mode:** E. Self-adversarial. Honest per A/B/C menu. Avoiding the Belavkin probe's over-claim.
**Companion:** `DWM_VERBATIM.md`, `DWM_SYRACUSE_IDENTIFICATION.md`, `DWM_MOMENT_PREDICTIONS.md`

---

## 0. One-line verdict

**Outcome B (near fit) — DWM Kraus form `M_v = ⟨v|U|0⟩` is literally correct for Syracuse's T_j structure AND DWM natively admits countably-infinite outcome alphabets (no extension required). All 4 structural moment-row predictions match QUALITATIVELY. But two material gaps remain: row (b) "zero" requires a bath-averaging-at-innermost-insertion argument cleaner than a naive tower property, and rows (d), (f) numerical values 0.108 / 0.609 are NOT predicted by DWM as such — they require explicit Kraus-channel numerical computation in R77's (1,4) basis. The fit is structurally clean but the quantitative discriminator hasn't been computed.**

---

## 1. The two load-bearing questions answered

### Q1: Is `M_v = ⟨v|U|0⟩` literally correct for Syracuse's T_j?

**YES.** Per `DWM_SYRACUSE_IDENTIFICATION.md` §3:

- `M_v^{(j, b_{[1,j-1]})} f(ξ) = 2^{-v/2} · A_v^{(j)}(ξ, b_{[1,j-1]}) · f(ξ · 2^{-v} mod 3^n)` is a valid Kraus operator on ℋ_n.
- POVM resolution: `∑_{v≥1} (M_v^{(j)})† M_v^{(j)} = ∑_{v≥1} 2^{-v} · I = I` ✓
- T_j as DWM-channel-averaged-Kraus: `T_j(ρ) = ∑_v M_v^{(j)} ρ (M_v^{(j)})†` — matches Tao's recursion via the Geom(1/2) outer expectation.
- Stinespring dilation `U_j ∈ U(ℋ_n ⊗ ℓ²(ℕ_{≥1}))` exists (Stinespring theorem; CP-ness of T_j is automatic from positive-kernel structure).
- Explicit `⟨v|U_j|0⟩ f = √(2^{-v}) · A_v^{(j)} · σ_{-v} f = M_v^{(j)} f` ✓ verified.

### Q2: Does countably-infinite outcome extension exist in published DWM literature?

**YES — it's already there, no extension required.**

- Wiseman 1996 §2 eq. (7): `∑_r F_r(T) = 1` constitutes a POVM on the result-space — **cardinality unrestricted by axiom.**
- Plenio-Knight 1998 §IV.A eq. (51): `ρ(t) = ∑_{n=0}^∞ ρ_A^{(n)}(t)` — countably-infinite path-space.
- Davies-Lewis 1970 / Davies 1976 "instrument": σ-additive CP-map-valued measure on a measurable σ-algebra of results — accommodates finite, countable, and continuous cardinalities axiomatically.
- The "gap" the brief asked about — countably-infinite extension of Davies/Wiseman-Milburn — **does not exist as a gap**. The framework is already defined for arbitrary measurable outcome spaces. (BvHJ 2009 §5 is the 2-outcome binomial restriction; that's a SPECIFIC paper within the Belavkin family, not a Davies/Wiseman-Milburn limitation.)

**Net of Q1+Q2:** DWM passes both structural tests cleanly. The brand label "Davies-Wiseman-Milburn quantum trajectory" is correct for what the Belavkin probe was actually doing.

---

## 2. Consistent P1-P7 score across all four probes

Re-scored on consistent DWM-labels (P1 = abelian past, P2 = NC system, P3 = level-graded operators, P4 = adapted-history-coupled, P5 = row d non-zero, P6 = row f non-zero, P7 = Fubini constant):

| Probe | P1 | P2 | P3 | P4 | P5 | P6 | P7 | Total | Notes |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|---|---|
| HP 1984 (Hudson-Parthasarathy) | ✗ | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ | 1/7 | Non-commutative filtration, time-local Itô |
| AP 2006 (Attal-Pautrat toy QSC) | ✗ | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ | 1/7 | Same Fock structure as HP |
| AFL 1982 (Accardi-Frigerio-Lewis) | ✓ | ✓ | ✗ | ✗ | ✗ | ✗ | ✓ | 3/7 | *-hom transports kill rows d, f |
| Belavkin-as-DWM (this probe) | ✓ | ✓ | ✓ | ✓ | ⚠ | ⚠ | ✓ | 4–5/7 | Audit-corrected from "7/7 claim"; P5/P6 qualitative only |

**Audit re-score corrections vs prior dispositions:**
- AFL was scored "4/7" in BELAVKIN_DISPOSITION; correct is **3/7** (P3, P4 are AFL failures, not "partial").
- Belavkin was scored "7/7" in BELAVKIN_DISPOSITION; correct is **4–5/7** (P5, P6 are qualitative-only; P7 is generic across ergodic-CP, over-credited).

**Net structural advance:** DWM ≥ AFL > HP/AP (by Δ = 1-2 properties on P3+P4). The DWM advance is real but narrower than the Belavkin probe claimed.

---

## 3. Outcome A vs B vs C — explicit comparison

**Outcome A (clean fit).** Would require:
- ✓ Kraus form literally correct (verified §1 Q1)
- ✓ Countably-infinite outcomes natively admitted (verified §1 Q2)
- ✓ All structural rows match (qualitative ✓)
- ✗ All NUMERICAL row values predicted (0.108, 0.609, 6.347×10⁻²) — Fubini value is R77-internal not DWM-specific; rows (d), (f) require explicit Kraus computation not done this probe
- ⚠ Row (b) clean derivation (improved over Belavkin probe but still relies on bath-averaging operator argument, not pure σ-algebra tower)

**Outcome B (near fit, structural gaps).** Matches current state:
- Structural identification is clean (Q1 ✓, Q2 ✓)
- P1-P4 cleanly accommodated
- P5, P6 qualitative-only (Mode-E gap)
- P7 over-credited (R77-internal, not DWM-specific)
- Row (b) derivation cleaner than Belavkin probe's but still requires operator-algebra step

**Outcome C (doesn't fit).** Would require a structural mismatch DWM can't accommodate. **No such mismatch found** — Geom(1/2) outcomes, adaptive Kraus, running-sum coarsening all admissible.

**Verdict: B.**

The structural identification is clean; the framework is correctly named (DWM, not Belavkin-1992-QSDE); the countably-infinite extension question is non-existent (DWM admits it natively). But the QUANTITATIVE verification — does DWM's Kraus-channel computation give 0.108 for row (d)? — is the genuine framework-discriminator, and that's a Mode-E gap deferred to numerical computation.

---

## 4. Mode-E gaps remaining

| Gap | Description | Effort |
|---|---|---|
| DWM-V-G1 | Davies 1976 monograph Ch. 2 verbatim "instrument" definition (physical book, no open scan found) | 1 hour institutional library |
| DWM-V-G2 | Wiseman-Milburn 2010 book Ch. 3 / Ch. 5 verbatim Kraus form (Cambridge CUP, no open chapters) | 1 hour book pull |
| DWM-MP-G1 | Numerical Kraus-channel computation of row (d) = 0.108 in R77 (1,4) basis | 4-8 hours numpy |
| DWM-MP-G2 | Same for row (f) = 0.609 | 4-8 hours numpy |
| DWM-MP-G3 | Formal operator-algebra factorization for row (b) bath-averaging-at-innermost-insertion | 2-4 hours derivation |

None block the verdict (B). All are downstream verification.

---

## 5. Implications for c = 7/45

**No change.** Per `THEOREM_C_745.md`, the leading c = 7/45 coefficient is rigorous unconditional via R75 + R76 + R77 + R64.B + HR74. DWM identification is a structural-framework finding, not a leading-coefficient finding.

What DWM helps with:
- The framework brand label is now correctly identified (Davies-Wiseman-Milburn, not Belavkin).
- The countably-infinite outcome question is closed (DWM admits natively).
- Future numerical Kraus-channel computation (Mode-E gap MP-G1, MP-G2) could provide quantitative match on rows (d), (f) — converting Outcome B to Outcome A.

---

## 6. Mode-E discipline notes

- Verbatim pulled this session: Wiseman 1996 (arXiv quant-ph/0302080), Plenio-Knight 1998 (arXiv quant-ph/9702007), Belavkin 1992 (Project Euclid CMP open access).
- pypdf extracted cleanly from all three (16/26/134 pages); the Belavkin probe's "pypdf blocked" claim continues to be wrong.
- Verbatim NOT pulled: Davies 1976 monograph (physical book), Wiseman-Milburn 2010 (Cambridge book), Davies-Lewis 1970 (Project Euclid auth wall).
- The Plenio-Knight 1998 citation chain (lines 879-881, 3927-3930 of extracted) provides the Davies-foundational route at a secondary verbatim level.

---

## 7. Adversarial honesty: what I did NOT do

- I did NOT compute the numerical row (d) value 0.108 in R77 (1,4) basis. The prediction "non-zero" is qualitative only.
- I did NOT pull Davies 1976 verbatim (paywalled physical book).
- I did NOT pull Wiseman-Milburn 2010 verbatim (paywalled book; frontmatter WebFetch was permission-denied).
- I did NOT formally factor row (b) through operator-algebra at the innermost-insertion level; the §2 derivation in DWM_MOMENT_PREDICTIONS is a sketch, not a theorem.
- I did NOT inflate AFL's P-score (audit-corrected to 3/7); previous Belavkin disposition's "AFL: 4/7" was over-credited.
- I did NOT inflate the DWM score to 7/7; honest count is 4-5/7 with P5, P6 qualitative-only and P7 generic.

The verdict (B) is honest: structurally clean fit, but quantitative discriminators (row d, row f magnitudes) are open Mode-E gaps. The Belavkin probe's "Outcome A" was over-claim; this Outcome B is the corrected verdict.

---

## 8. Verdict in one paragraph

**Outcome B. Davies-Wiseman-Milburn quantum trajectory is the correct named framework for Syracuse's structural shape: discrete-time Kraus operators `M_v^{(j, b_{[1,j-1]})} = 2^{-v/2} · A_v^{(j)}(·, b_{[1,j-1]}) · σ_{-v}` on ℋ_n = L²((ℤ/3^n)*, π_n), with countably-infinite outcome alphabet v ∈ ℕ_{≥1} ~ Geom(1/2), abelian observation filtration 𝔅 generated by running sums b_{[1,k]}, and history-dependent (adaptive) phase coupling INSIDE the Kraus operator. The literal `M_v = ⟨v|U|0⟩` form is verified by explicit Stinespring dilation. The countably-infinite extension question that the brief asked about is a non-issue — DWM admits arbitrary measurable outcome cardinalities by axiom (Wiseman 1996 eq. 7 POVM resolution; Davies-Lewis 1970 σ-additive instrument). P1-P4 are clean structural matches (4/7); P5, P6 are qualitative-only ("non-zero generic" — true under DWM and under any non-degenerate Kraus framework); P7 (Fubini constancy 6.347×10⁻²) is R77-internal and generic across ergodic-CP, not DWM-specific. Row (b) ≈ 0 derivation is cleaner than Belavkin probe's tower-argument (uses bath-averaging at the innermost X̃ insertion, requires full-𝔅 centering doing the work), but still relies on operator-algebra factorization not σ-algebra-tower alone. The framework is structurally home but the quantitative discriminator (do explicit Kraus-channel computations give 0.108 / 0.609 in R77 (1,4) basis?) is the open Mode-E gap that would convert this to Outcome A. P1-P7 re-scored consistently: HP 1/7, AP 1/7, AFL 3/7 (was 4/7, audit-corrected), DWM 4-5/7 (was claimed 7/7, audit-corrected). Net: DWM is the right home, narrower margin than previously claimed.**

---

## 9. Files

- This file: `C:/Collatz/DWM_DISPOSITION.md`
- Verbatim: `C:/Collatz/DWM_VERBATIM.md`
- Identification: `C:/Collatz/DWM_SYRACUSE_IDENTIFICATION.md`
- Moment predictions: `C:/Collatz/DWM_MOMENT_PREDICTIONS.md`
- Source extractions: `C:/tmp/dwm_belavkin_1992.txt`, `C:/tmp/dwm_plenio_knight_1998.txt`, `C:/tmp/dwm_wiseman_1996_qtmt.txt`
- Source PDFs (cached, %PDF magic): tool-results paths in DWM_VERBATIM §0
- Prior arc: `BELAVKIN_ADVERSARIAL_AUDIT.md`, `BELAVKIN_DISPOSITION.md`, `BELAVKIN_SYRACUSE_IDENTIFICATION.md`, `AFL_DISPOSITION.md`, `QSC_DISPOSITION.md`, `AMALG_FREENESS_SETUP.md`, `C1_TAO_RECURSION_FORM.md`, `THEOREM_C_745.md`
