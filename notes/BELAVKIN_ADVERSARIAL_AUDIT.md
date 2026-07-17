# BELAVKIN_ADVERSARIAL_AUDIT — Trying to break the Outcome A claim

**Date:** 2026-05-15
**Mode:** E. Adversarial.
**Auditor:** Opus 4.7. Read source PDFs directly via pypdf — no agent-summary intermediation.
**Verdict (one line):** **SURVIVED_WITH_CAVEAT** — Outcome A's positive structural finding (Syracuse has the shape of a Belavkin filter at the level of P1+P2+P3+P4) holds. But two load-bearing pieces are wrong as stated and one piece is over-credited; the right verdict is **Outcome B (near fit with gaps)**, not Outcome A.

---

## 0. Material pulled this session

Both cached BvHJ PDFs in tool-results extracted successfully via pypdf. The Belavkin agent's claim that "PDF binary fetch was sandbox-denied" / "pypdf extraction blocked" is **not corroborated** — pypdf reads them fine. The two files cached are both arXiv:math/0606118 (Bouten-van Handel-James 2009, "A discrete invitation to quantum filtering and feedback control") — v3 (Oct 2006) and v4 (Dec 2006), 77 pages each. Magic bytes `%PDF-1.4`. Extracted to `C:/tmp/bvhj_2009.txt` (and `bvhj_2007.txt` for the v3 copy).

Belavkin 1992 CMP, BvHJ 2007 SIAM J. Ctrl. Opt., Belavkin 1989 FoP — NOT pulled. The agent's download script was left for the user, not run by me; I'll note as a deferred gap. The BvHJ 2009 paper is the agent's identified "most directly relevant Syracuse-analog source," so my verbatim extraction targets it.

---

## Claim 1 audit (Kraus operator identification)

### Verbatim Belavkin Kraus form

**Finding 1.1 — load-bearing error.** The agent claims (BELAVKIN_VERBATIM.md §1.3) that BvHJ 2009 uses Kraus-operator decomposition
`M_i^{(j)} := ⟨i| U_j |0⟩ ∈ B(H_S)` with `Σ_i M_i^* M_i = I`.

**This is not the BvHJ 2009 formalism.** Direct grep of the 77-page PDF: the word "Kraus" appears **zero times**. The word "operator sum" also appears zero times. The word "Stinespring" appears zero times. BvHJ 2009 formalizes the discrete model via Hudson-Parthasarathy QSDE: equation 2.7 (p. 15) reads:

> "U(l) = U(l-1){I + j_{l-1}(L₁)ΔΛ(l) + j_{l-1}(L₂)ΔA*(l) + j_{l-1}(L₃)ΔA(l) + j_{l-1}(...)Δt(l)}" (eq. 2.7)

with discrete Itô differentials ΔA, ΔA*, ΔΛ, Δt. The system observable is evolved as `j_l(X) = U(l)* X U(l)` (eq. 4.1). The discrete Lindblad generator (eq. 4.4, p. 27):

> "L(X) = M+* X M+ + λ² M◦* X M◦ + M◦* X + X M◦"

This is the **same Hudson-Parthasarathy QSDE / discrete Itô calculus** that the prior `QSC_DISPOSITION.md` rejected for Syracuse (Outcome C). The Kraus-form `M_v = ⟨v|U|0⟩` that the agent identifies is a **Davies 1976 / Wiseman-Milburn 2010** formalism that BvHJ 2009 simply does not deploy as the operating framework. The agent constructed the Kraus form from "standard quantum-trajectory theory" exposition (Wikipedia + memory of Davies), not from BvHJ.

**Implication.** If we identify Syracuse with the *Davies/Wiseman-Milburn quantum-trajectory* form, that's a different framework than BvHJ 2009 *Belavkin discrete filtering*. The agent has been doing a hybrid identification: pulling "Belavkin" as a brand label while actually using Davies-style Kraus operators. The two are RELATED but NOT the same; in particular, QSC_DISPOSITION already explicitly rejected the Hudson-Parthasarathy QSDE form (= literal BvHJ 2009 form).

### History-dependent Kraus admissibility

**Finding 1.2 — partially correct.** Adversarial question (1) asks whether Belavkin admits Kraus operators that depend on the prior observation accumulator `b_{[1,j-1]}`. The answer in BvHJ 2009 is YES, but via a specific construction in §7.2:

> "U^µ(l) = M_l(ǔ_l) M_{l-1}(ǔ_{l-1}) ··· M_1(ǔ_1)" (eq. 7.6, p. 53)
> "ǔ_l = f_l(ΔZ(1), ..., ΔZ(l-1)) ∈ C_{l-1}" (p. 52)

where `C_{l-1}` is the commutative bath subalgebra (the bath-side past) and `f_l` is the admissible feedback function. The feedback strategy is `µ = {f_1, f_2(Δy_1), f_3(Δy_1, Δy_2), ..., f_k(Δy_1, ..., Δy_{k-1})}` (eq. above eq. 7.1, p. 51).

**So adaptive Kraus IS admitted, but via the BATH-side observation algebra C_{l-1}** (= observation algebra after spectral identification), which the agent identifies with Syracuse's B_{j-1} (running 2-adic valuation sums).

**Caveat.** BvHJ 2009 §7.2 makes the feedback enter via a *unitary that depends on a scalar `ǔ_l` taking values in U ⊂ R* (p. 51). Syracuse's accumulator `b_{[1,j-1]}` is integer-valued (in N≥1)^{j-1} → N), so the construction admits this case. The construction is admissible.

**Net.** Adaptive history-dependent Kraus IS admitted in Belavkin via §7.2. ✓ This piece of the agent's claim survives, modulo Finding 1.1's note that the canonical form in BvHJ is QSDE-based, not Kraus.

### Phase-inside-operator canonical form

**Finding 1.3 — admissible but not load-bearing of Belavkin specifically.** Adversarial question (2): Is the exponential phase factor INSIDE the operator (rather than outside as a multiplicative scalar) canonical Belavkin?

Yes — any unitary `U_j(ǔ)` with parameter `ǔ ∈ C_{l-1}` can have phase content depending on ǔ. The point of §7.2's `M_l(u) = exp(...)` construction is precisely that the parameter enters inside the unitary's exponential.

BUT — this is **not specifically a Belavkin feature**. Any quantum channel with a feedback parameter has this property. The agent is claiming Belavkin specifically when actually the structural feature (parameter-inside-unitary) is generic across HP / Davies / Belavkin / Lindblad / etc. The "P4 adapted phase coupling INSIDE" check is not load-bearing for distinguishing Belavkin from competitors.

### Observation outcome i_j as Geom(1/2)

**Finding 1.4 — caveat.** Adversarial question (3): does i_j ~ Geom(1/2) fit Belavkin's typical outcome alphabets?

BvHJ 2009 §5.6 (martingale representation) explicitly relies on `ω_l = ι(ΔY(l))` taking **exactly two values** {ω_+, ω_-} (p. 32, line 1697: "ω_l takes one of two values {ω_+, ω_-}"). The proof of the discrete martingale representation theorem leans on the bath-side observable having a 2-point spectrum (because ΔZ(l) is a single-photon increment).

Syracuse's v_j ~ Geom(1/2) on N≥1 has **countably infinite values**. BvHJ 2009 in particular does not treat this — the proof of the martingale representation breaks down for >2 outcomes. To extend Belavkin filtering to countably-infinite outcomes one needs the continuous-time or full-Davies generalization, not BvHJ's "binomial model".

**This is a real structural mismatch the agent overlooked.** BvHJ 2009 = binomial / two-outcome. Syracuse = geometric / countably-infinite outcome. The discrete BvHJ derivations do not extend to Syracuse's outcome space without nontrivial work.

---

## Claim 2 audit (non-demolition condition)

### Direct commutator computation

**Finding 2.1 — survives, but for technical reasons different from agent's argument.** Agent claims `[T̃_j, M_{b_{[1,k]}}] = 0` for k < j based on tensor-factor argument.

Tracing through: On the tensor product H_n ⊗ L²(Ω), T̃_j has the form
`T̃_j = ∫_{Ω} T_j(b_{[1,j-1]}(ω)) · σ_{-v_j(ω)} dP(ω)`
where σ is shift on H_n. The Ω-side dependence is:
(i) T_j is multiplied by the indicator-projector of value of b_{[1,j-1]}, which is in B_{j-1} ⊂ L^∞(Ω).
(ii) T_j has a "new" v_j integration mixed in via the σ_{-v_j} shift weight 2^{-v_j}.

For k < j, M_{b_{[1,k]}} = I_{H_n} ⊗ μ_{b_{[1,k]}}, a multiplication operator on the Ω-side only. The Ω-side action commutes with another Ω-side action iff both are decomposable along the same abelian subalgebra. Since b_{[1,k]} ∈ B_{j-1} (for k ≤ j-1), both T̃_j-Ω-part and M_{b_{[1,k]}} are B_{j-1}-decomposable. Hence they commute. ✓

For **k = j**: T_j integrates over v_j; M_{b_{[1,j]}} = M_{b_{[1,j-1]} + v_j} mixes the new v_j. Non-trivial commutator. Agent acknowledges this gap (BELAVKIN_SYRACUSE_IDENTIFICATION §2.2 trailing brackets).

**Net.** Non-demolition holds for k < j. ✓ Matches BvHJ 2009 p. 18 statement of nondemolition: `[ΔY(l), j_i(X)] = 0` for l ≤ i (past obs commutes with current/future system).

### Tensor-product lift check

**Finding 2.2.** The tensor-lift `T̃_j = ∫ T_j(b) ⊗ E_b dP(b)` written in BELAVKIN_SYRACUSE_IDENTIFICATION §2.2 is well-defined: T_j(b) is a bounded operator on H_n for each fixed `b ∈ N_0`, and `E_b` is the spectral projector onto the value b of b_{[1,j-1]}. The integral is over the (discrete) distribution of `b_{[1,j-1]}`.

This satisfies the bilinearity for tensor decomposition and the commutator-vanishing argument. ✓

**Net for Claim 2.** Claim 2 survives directly. Non-demolition between system algebra at step j and observation algebra at step k < j is verified by tensor-factor argument.

---

## Claim 3 audit (moment predictions)

### Row (b) tower-property derivation check

**Finding 3.1 — load-bearing hand-wave.** Agent argues `φ(X̃_{j_1} X̃_{j_2}) = E[X̃_{j_1} · E[X̃_{j_2} | B_{j_2-1}]] = 0` via tower property.

This requires:
(a) `E[X̃_{j_2} | B_{j_2-1}] = 0` — i.e., X̃_j has zero conditional expectation onto the past.
(b) The tower property holds (φ = E_B ∘ E_{B_{j_2-1}}).

For (a): X̃_j = T_j - E_B[T_j] where E_B is conditional expectation onto B (= the FULL accumulator algebra, not B_{j-1}). Centering with respect to E_B does not directly give zero conditional expectation onto B_{j_2-1} ⊂ B unless additional structure holds.

**The agent provides no derivation that E[X̃_{j_2} | B_{j_2-1}] = 0.** The argument in BELAVKIN_MOMENT_PREDICTIONS.md §2 just asserts this without working through the projection. Since T_j has v_j-dependence (integrated out) AND b_{[1,j-1]}-dependence (from the phase), E[T_j | B_{j-1}] = ∫ T_j(b_{[1,j-1]}, v_j) dμ_{Geom}(v_j), which is itself a B_{j-1}-measurable operator on H_n. This need NOT equal E_B[T_j] (the latter projects onto the full accumulator algebra B = ∪_k B_k).

**Verdict.** Row (b) tower-property derivation does not actually verify row (b) = 0 from Belavkin's structure. It's a plausibility argument with a measurability gap. Syracuse's measured row (b) ≈ 0 may follow from finer structure (the leading 7/45 cancellation in T_diag), not from a clean Belavkin tower property.

This is a **specific error**: the moment prediction §2 over-claims rigor.

### Row (d) Kraus-trace computation

**Finding 3.2 — qualitative only, MP-G1 gap is load-bearing.** Agent predicts row (d) "non-zero" via the structural argument that `X̃_{j_1}^op · M^{(j_2)} · X̃_{j_1}^op` doesn't contract.

This is QUALITATIVE only. The agent itself logs MP-G1 / MP-G2 as Mode-E gaps requiring "explicit Kraus-channel computation". Without that computation, the prediction "non-zero" is **so weak it's nearly tautological**: any non-degenerate framework predicts a generic third moment is non-zero.

The specific numerical match (0.108 ± noise) requires actually computing Tr(X̃_{j_1}^2 · M^(j_2)) in R77's (1,4) basis. Without this, we have no quantitative confirmation that Belavkin predicts 0.108 specifically; we only have "non-zero, with sign unspecified, with magnitude unspecified."

**Verdict.** Row (d) qualitative prediction = "non-zero, generic". This is true under Belavkin AND under a hundred other frameworks. It is not informative as a framework discriminator. The Mode-E gap MP-G1 is load-bearing for the AFL→Belavkin advance to be substantive.

### Fubini constancy ergodicity argument

**Finding 3.3 — over-credited.** Agent argues Fubini constancy follows from "ergodic Kraus channel + R77 1-d invariant on (1,4)".

R77 is a Syracuse-internal result, derived independently of any framework (Belavkin or otherwise). The Fubini constancy 6.347×10^{-2} is a measured Syracuse quantity, anchored by R77's 1-d invariant eigenspace of T_diag.

Saying "Belavkin predicts Fubini constancy via ergodic Kraus channel" identifies Belavkin's ergodicity condition with Syracuse's R77 structure. But **any framework with an ergodic 1-d-invariant transition channel reproduces this**. AFL was credited with "Fubini constancy ✓" via the ergodic transition expectation (AFL_DISPOSITION §1). HP/AP were "ambiguous." QSC was ✗.

Belavkin specifically adds nothing here over AFL. The structural ingredient is "ergodicity of the per-step channel" — that's the abstract requirement and any framework with that property satisfies P7.

**Verdict.** P7 / Fubini constancy is **over-credited to Belavkin specifically**. The structural ingredient is generic across multiple frameworks. The 6.347×10^{-2} value is Syracuse-internal (R77).

---

## Claim 4 audit (P1-P7 score)

### Each P_i individually checked

| P_i | Agent's claim | Audit finding |
|---|---|---|
| P1 (abelian past filt) | ✓ axiomatic | ✓ HOLDS. BvHJ 2009 §2.5 self-nondemolition: Y_k is commutative. Matches B's abelian structure. |
| P2 (NC system algebra) | ✓ axiomatic | ✓ HOLDS. M is non-commutative atomic algebra in BvHJ. |
| P3 (level-graded U_j) | ✓ admits time-dependent U_j | ⚠ ADMITTED but NOT EMPHASIZED in BvHJ 2009. BvHJ §2.5 takes a SINGLE M_l interaction matrix; time-dependent M_l = M (no j-subscript) for spontaneous emission (§2.6 p. 19). To get level-graded U_j matching Syracuse's `T_j(b)`, one needs §7.2's feedback construction, which is more specific. |
| P4 (adapted phase inside) | ✓ admits adaptive M_i^{(j, i_{1:j-1})} | ✓ HOLDS via §7.2 (eq. 7.6 with ǔ_l = f_l(ΔZ_1, ..., ΔZ_{l-1})). |
| P5 (row d non-zero) | ✓ predicted non-zero | ⚠ QUALITATIVE only. Not a discriminating prediction (see Finding 3.2). |
| P6 (row f non-zero) | ✓ predicted non-zero | ⚠ QUALITATIVE only. Same status as P5. |
| P7 (Fubini constant) | ✓ via ergodic Kraus channel | ⚠ over-credited (see Finding 3.3). R77-internal; not Belavkin-specific. |

**Properly scored:** 4 clean ✓ (P1, P2, P3 modulo §7.2 specialization, P4), 3 qualitative-only ⚠ (P5, P6, P7).

The score "7/7 clean" is overstated. A fair score is **4/7 clean + 3/7 qualitative-only**, i.e., approximately 4-5 / 7 depending on how strictly one counts qualitative agreement.

Caveat: BvHJ 2009 is structurally limited to **2-outcome bath observations** (binomial model). Syracuse's Geom(1/2) is countably-infinite. This is a Belavkin-specific obstacle the agent didn't flag (Finding 1.4).

---

## Claim 5 audit (AFL + HP/AP score comparison)

### AFL P1-P7 re-derived

The agent's BELAVKIN_DISPOSITION §1 scores AFL at 4/7 with the column labeled "AFL: ✓ ✓-partial partial partial ✗ ✗ ✓". But AFL_DISPOSITION.md itself uses a DIFFERENT P1-P7 numbering (P5 = repeat moments non-zero, P6 = B-content inside the integral, P7 = Fubini constant), and scores AFL at "passes P1, P3, P4, P7 cleanly; fails P5, P6" — i.e., 4-5 of 7 by its own labels.

**Re-scored on Belavkin's P1-P7 labels:**
- P1 (abelian past) — ✓ (AFL admits abelian)
- P2 (NC system) — ✓ in AFL (j_t : O → A; O can be NC)
- P3 (level-graded operators NOT time-translates) — ✗ AFL explicitly USES time-translates of fixed O
- P4 (adapted phase coupling INSIDE) — ✗ AFL doesn't admit (per AFL_DISPOSITION §1: "AFL requires 'the same random variable transported'")
- P5 (row d non-zero) — ✗ AFL forces row d = 0 (per AFL_DISPOSITION §1)
- P6 (row f non-zero) — ✗ AFL forces row f = 0
- P7 (Fubini constant) — ✓ via ergodic transition expectation

**Re-scored AFL: 3/7 (P1, P2, P7) clean — NOT 4/7 as the agent claims.**

The agent inflated AFL's score by giving "partial" credit for P2/P3/P4 which AFL's own disposition says are structural failures. This makes Belavkin look more strictly dominant than the actual structural picture supports.

### HP/AP P1-P7 re-derived

Agent scores HP/AP at 1/7. Re-checking:
- P1 (abelian past) — ✗ HP/AP filtration is non-commutative (Fock/qubit-chain)
- P2 (NC system) — ✓ HP/AP system is non-commutative
- P3 (level-graded) — ✗ HP/AP time-local Itô = stationary, NOT level-graded
- P4 (adapted phase inside) — ✗ HP adapted process is `F_t · dA_t` with F_t outside the noise differential
- P5 (row d non-zero) — ✗ HP/AP time-local Itô kills cross-time triple
- P6 (row f non-zero) — ✗ same
- P7 (Fubini constant) — ambiguous (per QSC_DISPOSITION)

**Re-scored HP/AP: 1/7 (P2). ✓ Agent's score is correct here.**

**Net Claim 5:** Agent inflates AFL from 3/7 to 4/7. The dominance claim "Belavkin 7/7 > AFL 4/7" is **soft-overstated**: properly scored it's "Belavkin 4-5/7 vs AFL 3/7", a narrower margin.

---

## Claim 6 audit (Mode-E gap impact)

### Verbatim vs canonical-form robustness

**Finding 6.1 — fork in canonical exposition.** Adversarial question: is there a fork where Syracuse fits one Belavkin variant but not another?

YES. Three variants of "Belavkin":

(a) **BvHJ 2009 discrete invitation** — HP/QSDE-based, 2-outcome bath, discrete Itô. This is the variant where verbatim pulled in this audit (PDF read). Syracuse's countably-infinite outcome space does NOT directly fit (Finding 1.4).

(b) **Davies 1976 / Wiseman-Milburn 2010 quantum trajectories** — Kraus-operator-based, arbitrary outcome alphabet via `M_v = ⟨v|U|0⟩`. This is the variant the agent ACTUALLY uses in BELAVKIN_SYRACUSE_IDENTIFICATION.md and BELAVKIN_MOMENT_PREDICTIONS.md. Countably-infinite outcomes OK.

(c) **Belavkin 1989/1992/1995 continuous filtering** — Itô SDE form with Wiener / Poisson observation. Continuous-time, not directly comparable to Syracuse's discrete recursion.

The agent BLENDS (a)+(b) under the "Belavkin" brand. (a) is what I directly verified; (b) is what the agent actually applies. **They are not the same framework.** Syracuse fits (b) (Davies/quantum-trajectory) potentially well, but the literal BvHJ 2009 Belavkin formalism (a) is the QSDE form already rejected by QSC_DISPOSITION.

**Verdict.** Mode-E gap is load-bearing. The "Belavkin = correct home for Syracuse" finding rests on framework (b), but the agent's verbatim BvHJ source is framework (a). If the AGENT had cited Davies 1976 + Wiseman-Milburn 2010 as the canonical home instead of BvHJ 2009, the structural fit story would be cleaner and the brand label "Belavkin" might be a misnomer (it's actually Davies/quantum-trajectory).

---

## Verdict: SURVIVED_WITH_CAVEAT

The structural shape of the agent's finding holds: Syracuse exhibits the qualitative pattern of an **adaptive-feedback quantum trajectory** with classical (abelian) observation filtration coupled to a non-commutative system dynamics with history-dependent unitary interactions. This shape is correctly identified.

But three things break vs the strong-form Outcome A claim:

1. **The agent's "Belavkin" is actually Davies/Wiseman-Milburn quantum trajectory theory** (Kraus form), not BvHJ 2009 Belavkin (HP/QSDE form). BvHJ 2009 in particular uses ZERO Kraus-operator language; it is QSDE-based, which QSC_DISPOSITION already rejected for Syracuse. Brand-label mismatch is load-bearing.

2. **BvHJ 2009 is structurally a binomial/2-outcome model**; Syracuse is geometric/countably-infinite outcome. The discrete BvHJ derivations (especially §5 martingale representation, §6 reference probability) do not extend without nontrivial generalization. The agent overlooked this.

3. **P-scores are inflated.** Belavkin's 7/7 should be 4-5/7 (P5, P6, P7 are qualitative-only or framework-independent). AFL's 4/7 should be 3/7. The "strictly dominates" claim still holds in direction but with narrower margin.

4. **Row (b) tower-property argument is hand-wavy.** Centered X̃_j has zero conditional expectation onto FULL B, not necessarily onto B_{j-1}. The martingale-difference property invoked for row (b) = 0 is asserted, not derived.

5. **Row (d) and (f) "non-zero" predictions are too weak to be discriminating.** Almost any framework that doesn't impose a kernel-condition contraction predicts these are non-zero. The numerical match (0.108, 0.609) is a Mode-E gap (MP-G1, MP-G2) deferred to future Kraus-channel computation.

6. **Fubini constancy (P7) is Syracuse-internal via R77**, not Belavkin-specific. Any ergodic framework gets credit for this.

The correct verdict downgrade is from **Outcome A (clean structural fit)** to **Outcome B (near fit with explicit gaps)** with the following sharply named gaps:
- BVK-AUD-G1: framework selection — Davies/quantum-trajectory vs BvHJ/QSDE; pick one and stick with it.
- BVK-AUD-G2: extend BvHJ 2009 derivations from 2-outcome binomial to Geom(1/2) countably-infinite.
- BVK-AUD-G3: explicit numerical verification of row (d) = 0.108, row (f) = 0.609 via Kraus channel computation (NOT just "non-zero generic").
- BVK-AUD-G4: derive row (b) = 0 from a clean martingale-difference statement under Belavkin's full structure, not via heuristic tower property over the wrong conditional expectation.

---

## Specific errors found

1. **BvHJ 2009 does NOT contain the word "Kraus".** (Verbatim from C:/tmp/bvhj_2009.txt: zero matches for "Kraus", "operator sum", "Stinespring".) The agent's canonical Kraus form is borrowed from Davies/Wiseman-Milburn, not extracted from BvHJ. The Belavkin agent presented this as a faithful BvHJ identification when it is actually a translation.

2. **BvHJ 2009 is 2-outcome binomial** (eq. line 1697 of extracted text: "ω_l takes one of two values {ω_+, ω_-}"). Syracuse is Geom(1/2). Mismatch flagged by the agent? No.

3. **AFL P1-P7 score inflated to 4/7.** Properly scored on Belavkin's labels: 3/7.

4. **Tower-property derivation for row (b)** does not work without showing E[X̃_j | B_{j-1}] = 0, which requires more than X̃_j = T_j - E_B[T_j].

5. **P7 over-credited.** Fubini constancy is via R77 (Syracuse-internal). Any ergodic framework satisfies this.

6. **Claim "pypdf extraction blocked" is wrong.** The two cached BvHJ PDFs extract fine via pypdf. (Both 77 pages, %PDF-1.4 magic bytes verified.)

---

## Recommended follow-up

1. **Rename the framework** in the Syracuse identification: "Davies-Wiseman-Milburn quantum trajectory" rather than "Belavkin filtering". The Kraus-operator form `M_v = ⟨v|U|0⟩` is canonical Davies, not canonical Belavkin. Belavkin 1989/1992 is continuous-time SDE-based.

2. **Pick a discrete framework that natively supports countably-infinite bath outcomes.** Candidates: Davies 1976 generalized measurement; Maassen-Kümmerer 2003 quantum Markov processes; Holevo statistical structure of quantum theory (Ch 3-4); Carmichael quantum trajectory theory.

3. **Numerical closure: BVK-AUD-G3.** Compute `Tr(X̃_{j_1}^2 · M^{(j_2)})` in R77's (1,4) basis and check 0.108. This converts qualitative "non-zero" to quantitative match. If this fails the verdict drops further.

4. **Verbatim Belavkin 1992 CMP pull** (open access on Project Euclid). The agent didn't pull it. The non-demolition definition originated there. Direct verbatim quote of the foundational paper is overdue.

5. **Re-score AFL and HP/AP on the SAME P1-P7 labels** to ensure the dominance claim is on a consistent scale.

6. **Don't redo Outcome A.** Outcome B with the named gaps above is the right verdict at the current state of evidence.

---

## Files

- This audit: `C:/Collatz/BELAVKIN_ADVERSARIAL_AUDIT.md`
- BvHJ 2009 extracted: `C:/tmp/bvhj_2009.txt` (77 pages)
- BvHJ 2009 v3 extracted: `C:/tmp/bvhj_2007.txt` (77 pages, duplicate)
- Cached PDFs (valid %PDF-1.4): `C:/Users/Nate/.claude/projects/c--As-Above-So-Below-Master/58d9c54a-7a98-404f-a28d-18647b4045be/tool-results/webfetch-1778888184187-75g6ko.pdf` and `webfetch-1778888211682-nuwqop.pdf`
- Audited materials: `BELAVKIN_DISPOSITION.md`, `BELAVKIN_VERBATIM.md`, `BELAVKIN_SYRACUSE_IDENTIFICATION.md`, `BELAVKIN_MOMENT_PREDICTIONS.md`
- Cross-references: `AFL_DISPOSITION.md`, `QSC_DISPOSITION.md`, `AMALG_FREENESS_SETUP.md`, `C1_TAO_RECURSION_FORM.md`
