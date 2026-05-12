# T_N_OFF_LIN_SPEC — Phase 1 specification extraction for Off_lin

**Date:** 2026-05-12. Phase 1 of the T_N = T_diag + Off_lin construction probe, Fork 1 from R76 §11 INCONCLUSIVE.
Wilson (analyst) reporting to Nathan.

---

## Headline

> **R77 sketch §5 articulates Off_lin as a PROCEDURE in a different state space (6-dim cross-frequency bilinear), not as a 2x2 matrix on (P_+, P_−). The reduction to a 2x2 over Q on the (P_+, P_−) basis requires a closure step that R77 sketch §5 / §6 / §7 do not provide and that the project-internal documents (T_lead_spectrum.md §3, T_lead_2x2.py, result_77_T_diagonal.py) explicitly flag as "more careful analysis needed" / "Open".**

The diagonal piece T_diag = (1/5)·[[1, 1], [4, 4]] is rigorously derivable from R77 sketch §5's procedure because the v = v' diagonal collapses each μ̂_n(ξ·2^{−v}) μ̂_n^*(ξ·2^{−v}) into |μ̂_n^a(ξ·2^{−v})|², which after summing over ξ projects back onto the same {P_n^{ab}(c)} basis (cover map ξ → ξ·2^{−v} is a coprime-unit shuffle).

The off-diagonal v ≠ v' piece does NOT collapse onto {P_n^{ab}(c)} because μ̂_n(ξ·2^{−v}) μ̂_n^*(ξ·2^{−v'}) is a bilinear at TWO DIFFERENT frequencies, with phase factor A_v(ξ) A_{v'}^*(ξ) = e^{−2πi ξ (2^{−v} − 2^{−v'})/3^{n+1}}. After summing over ξ, this produces NEW objects M_n^{ab}(c; v, v') that are not linear functions of {P_n^{ab}(c)} under any closure R77 sketch identifies.

**Without that closure, Off_lin is not a well-defined 2x2 matrix over Q.** The verbal heuristic in T_lead_spectrum.md §3 (λ_2 = P(v=1) = 1/2) is acknowledged as unspecified ("more careful analysis needed", §3 line 67). The numerical fit in T_lead_2x2.py is data-driven, not derivation-driven.

This matches the pre-registered failure mode **H_OFF_LIN_UNDERSPECIFIED**.

---

## (a) Verbatim text of R77 sketch §5

R77 sketch §5 is titled **"Tao recursion induces T_lead"**. The full content:

> Tao recursion gives μ̂_{n+1} from μ̂_n. The class-conservation rule (from R66 + the chain dynamics):
> - v even → r' ≡ 1 mod 3: μ̂_{n+1}^+ contribution
> - v odd → r' ≡ 2 mod 3: μ̂_{n+1}^− contribution
>
> So:
> > μ̂_{n+1}^+(ξ) = Σ_{v even, v≥2} 2^{−v} A_v(ξ) μ̂_n(ξ·2^{−v} mod 3^n)
> > μ̂_{n+1}^−(ξ) = Σ_{v odd, v≥1} 2^{−v} A_v(ξ) μ̂_n(ξ·2^{−v} mod 3^n)
>
> where A_v(ξ) = e^{−2πi ξ 2^{−v}/3^{n+1}}.
>
> The level-n+1 class-resolved moments P_{n+1}^{ab}(c) are bilinear in μ̂_{n+1}^a, μ̂_{n+1}^b*, hence quadratic in μ̂_n. Substituting Tao's recursion expresses them as quadratic forms in {P_n^{ab}(c)}.
>
> This defines the **transfer operator T** acting on the 6-dim vector P_n = (P_n^{++}(1), P_n^{++}(2), P_n^{−−}(1), P_n^{−−}(2), Re P_n^{+−}(1), Im P_n^{+−}(1), Re P_n^{+−}(2), Im P_n^{+−}(2)) (8 reals, with constraints from total mass).
>
> S_n = linear combination of P_n entries; R_n similarly.

That is the complete content of R77 sketch §5.

---

## (b) Diagnosis — is this concrete enough to build T_N as a 2x2 over Q?

**No.** Three concrete obstructions:

### Obstruction 1 — §5's claim "quadratic forms in {P_n^{ab}(c)}" is the OPEN CLOSURE, not a derivation

The §5 sentence

> "Substituting Tao's recursion expresses them as quadratic forms in {P_n^{ab}(c)}."

is the **conclusion** R77 sketch §5 asserts but does NOT derive. Expanding Tao's recursion in P_{n+1}^{ab}(c) gives:

  P_{n+1}^{++}(c) = Σ_{ξ ≡ c, 3∤ξ} |μ̂_{n+1}^+(ξ)|²
                  = Σ_{ξ ≡ c, 3∤ξ} (Σ_{v even} 2^{−v} A_v(ξ) μ̂_n(ξ·2^{−v}))
                                  · (Σ_{v' even} 2^{−v'} A_{v'}^*(ξ) μ̂_n^*(ξ·2^{−v'}))
                  = Σ_{v, v' both even} 2^{−v−v'} · Σ_{ξ ≡ c, 3∤ξ} A_v(ξ) A_{v'}^*(ξ) · μ̂_n(ξ·2^{−v}) μ̂_n^*(ξ·2^{−v'})

Define the inner cross-frequency object:

  Q_n^{++}(c; v, v') := Σ_{ξ ≡ c, 3∤ξ} A_v(ξ) A_{v'}^*(ξ) · μ̂_n(ξ·2^{−v}) μ̂_n^*(ξ·2^{−v'})

Then

  P_{n+1}^{++}(c) = Σ_{v, v' both even} 2^{−v−v'} Q_n^{++}(c; v, v')

The case v = v' gives, by A_v · A_v^* = 1 and Tao's class-shuffle:

  Q_n^{++}(c; v, v) = Σ_{ξ ≡ c, 3∤ξ} |μ̂_n(ξ·2^{−v})|²

The substitution u := ξ·2^{−v} is a unit shuffle on (Z/3^{n+1})^×; combined with R66's class flow (v even ⇒ class preserved; v odd ⇒ class flipped), this collapses Q_n^{++}(c; v, v) onto a linear combination of {P_n^{++}(c'), P_n^{−−}(c'), P_n^{+−}(c')} entries — this is exactly the derivation in `result_77_T_diagonal.py`. The diagonal v = v' piece IS a closed linear function of {P_n^{ab}(c)}.

The case v ≠ v' is qualitatively different:

  Q_n^{++}(c; v, v') = Σ_{ξ ≡ c, 3∤ξ} e^{−2πi ξ (2^{−v} − 2^{−v'})/3^{n+1}} · μ̂_n(ξ·2^{−v}) μ̂_n^*(ξ·2^{−v'})

This is a CROSS-FREQUENCY bilinear with a NON-TRIVIAL phase factor e^{−2πi ξ d_{v,v'}/3^{n+1}} where d_{v,v'} := 2^{−v} − 2^{−v'} ∈ Z/3^{n+1}. The two μ̂_n arguments are DIFFERENT lattice points (ξ·2^{−v} ≠ ξ·2^{−v'}), so this does NOT reduce to a sum of |μ̂_n(u)|² objects under any change of variable. It produces a new bilinear object:

  R_n^{++}(c; d) := Σ_{u, u': 3∤u, 3∤u'} (substitution-dependent) μ̂_n^+(u) μ̂_n^{+*}(u')  ⋅  (phase depending on d)

These cross-frequency μ̂_n^a(u) μ̂_n^{b*}(u') bilinears do NOT lie in span{P_n^{ab}(c)}. They are 3^n × 3^n × (constraint) dimensional in principle; only the diagonal u = u' projects back onto {P_n^{ab}(c)}.

**The §5 sentence "expresses them as quadratic forms in {P_n^{ab}(c)}" is therefore false as written for v ≠ v'.** It would require an additional closure identity — something like a Plancherel-type collapse showing Σ_{u ≠ u'} (cross-freq) ≡ linear combination of {P_n^{ab}(c)} — and R77 sketch does not provide it.

### Obstruction 2 — §5's "6-dim vector" / "8 reals with constraints" is not the 2x2 (P_+, P_−) basis

R77 sketch §5 explicitly says T acts on a 6-dim vector (8 reals with constraints), enumerating the basis:

> P_n = (P_n^{++}(1), P_n^{++}(2), P_n^{−−}(1), P_n^{−−}(2), Re P_n^{+−}(1), Im P_n^{+−}(1), Re P_n^{+−}(2), Im P_n^{+−}(2))

This is a 6-dim or 8-real-dim object. The 2-dim (P_+, P_−) reduction R76 §11 obtains is by IMPOSING the structural collapse identities (P^{+−}(c) = 0; P^{++}(1) = P^{++}(2); P^{−−}(1) = P^{−−}(2)), which collapse the 6 dims to 2.

That collapse is established for the **state vector** at level n ≥ 2 (rigorous from R76 §11), but for **T to act on the 2-dim subspace as a 2x2 matrix**, we need T to preserve this subspace. That is:

  Claim: If P_n^{+−}(c) = 0 and P_n^{++}(1) = P_n^{++}(2) and P_n^{−−}(1) = P_n^{−−}(2), then T(P_n) has the same structure.

This is true for T_diag (the diagonal piece) — `result_77_T_diagonal.py` derivation explicitly uses the class-c symmetry assumption. It is **not derived** anywhere in R77 sketch / R76 §11 / T_lead_spectrum.md for Off_lin (the off-diagonal piece).

In other words: even if Obstruction 1 were resolved (cross-frequency bilinears closed onto {P_n^{ab}(c)}), the 2x2 reduction requires showing that the resulting 6×6 operator's restriction to the (P^{++}(c=1) = P^{++}(c=2), P^{−−}(c=1) = P^{−−}(c=2), P^{+−} = 0) subspace is invariant. R77 sketch §5 asserts this implicitly ("S_n = linear combination of P_n entries; R_n similarly") but does not establish it.

### Obstruction 3 — the project's own honest documentation flags this

`result_77_T_lead_spectrum.md` §3 ("Conjectured exact rate λ_2 = 1/2") attempts to extract λ_2 = 1/2 from the v ≠ v' off-diagonal sum and explicitly stops, with two passages worth quoting verbatim:

> Next contributing terms: v = 1, v' = 3 (or v = 3, v' = 1), with weight 2^{−1−3} = 1/16. The phase 2^{−1} − 2^{−3} mod 3^{n+1} has **3-adic valuation = 1** (computed: 2^{−1} − 2^{−3} = (4 − 1)/8 = 3/8, hence v_3 = 1 in 3-adic). This 1-step 3-adic gap means the phase character at level n+1 reduces to a level-n character with non-trivial sum.
>
> Working through this: the **leading off-diagonal eigenvalue of the v ≠ v' contributions is** ~ **2·(1/4)·(weight) = 1/2** when summed over the leading bilinear couplings. The factor 2 comes from the (v, v') ↔ (v', v) symmetry; the (1/4) from 2^{−1−1}.
>
> **More precisely:** the leading off-diagonal term is the (v=1, v'=1) coincidence on cross-frequency, which contributes weight P(v=1)² = 1/4 with sign +1 (not the trivial diagonal). Combined with the (1, 4) eigenvector projection: **λ_2 = 4·(1/4) = 1**? — no, more careful analysis needed.
>
> Actually the cleanest derivation: **λ_2 = P(v=1) = 1/2** because at each level k → k+1, the v=1 contribution is the "fresh" perturbation that, once integrated through Plancherel, gives a contraction by 1/2.

The text directly admits "more careful analysis needed" and offers two contradictory candidate derivations (4·(1/4)=1 vs P(v=1)=1/2) without selecting. This is the project's own statement that Off_lin is unspecified.

The ledger §6 in the same document places this on the OPEN side:

> ### Open (analytical work for fully rigorous closure)
> - ✗ Off-diagonal exact bilinear-sum analysis to confirm λ_2 = 1/2 from Tao's recursion

The "Off-diagonal exact bilinear-sum analysis" is exactly the closure step needed to realize Off_lin as a 2x2 matrix. The project's own ledger says it's open.

`T_lead_2x2.py` numerically fits a 2x2 T using observed (P_+, P_−) deviations at k = 2..5. That is fit-from-data, not derived-from-recursion. The eigenvalues output by the fit are conjectural numerical values, not exact Q values; and §11's H_T_N_DIFFERENT_EIGENVALUE / H_T_N_CONFIRMS_RATE_HALF distinction requires exact Q to decide.

`result_77_2_T_N_construction.py` builds a 1-dim "T_N" as kappa_N := delta_+(N+1) / delta_+(N) (the level-N ratio) — that is also fit-from-data, and acknowledged as such in the script's docstring header.

---

## (c) If explicit: write the matrix entries

Not applicable — Phase 1 finds the spec is not concrete enough to build Off_lin as a 2x2 matrix over Q. See (d).

---

## (d) Diagnosis: the specific gap

Off_lin's 2x2 matrix entries over Q are NOT specified by R77 sketch §5. The gap is structural:

**Gap 1 — Cross-frequency closure.** §5 says off-diagonal bilinears collapse onto {P_n^{ab}(c)}, but the v ≠ v' terms involve cross-frequency μ̂_n(u) μ̂_n^*(u') with u ≠ u'. These do NOT collapse onto same-frequency {P_n^{ab}(c)} under any identity R77 sketch provides. A closure would require either:
- A Plancherel-type formula expressing Σ_{u ≠ u'} μ̂_n(u) μ̂_n^*(u') (with character weight) as a linear combination of {P_n^{ab}(c)} entries.
- A separately-articulated bilinear-pair recursion stating that the cross-frequency object is itself representable in the same-frequency basis.

Neither is provided in R77 sketch §5–§9.

**Gap 2 — 2D subspace invariance for Off_lin.** Even if Gap 1 were closed, Off_lin (as a 6x6 operator on the full P_n vector) must preserve the (P^{+−} = 0, class-c symmetric) 2D subspace for the 2x2 reduction to be a valid restriction. This is asserted but not derived.

**Gap 3 — Selection of which leading term gives the 1/2.** T_lead_spectrum.md §3 offers two contradictory candidate paths (λ_2 = 4·(1/4) = 1 vs λ_2 = P(v=1) = 1/2). Both are heuristic. Neither is derived to the matrix-entry level.

---

## (e) Can the gap be bridged in one session?

**Gap 1 alone is a substantive derivation** — it requires deriving a Plancherel-type closure for cross-frequency bilinears, which is mathematically of the same scale as the original derivation of T_diag (which took its own dedicated session in R77 §1). Sketch-level: this would mean computing, for each pair (v, v', class, c), the explicit contribution of Q_n^{ab}(c; v, v') in terms of {P_n^{ab}(c)} and any leftover cross-frequency objects, then checking that the cross-frequency leftovers are zero, decouple, or are themselves expressible in {P_n^{ab}(c)}.

This is exactly the "off-diagonal exact bilinear-sum analysis" T_lead_spectrum.md §6 lists as open. R77 sketch §10 estimates this at "1-2 hours of focused implementation", but R77 sketch's §10 estimate refers to running the IMPLEMENTATION (i.e., assuming the closure exists, build T_N matrices at level N = 1, 2, 3 numerically and check spectrum). The DERIVATION of the closure itself is treated as already-done in §5's claim — which is precisely the asserted-but-not-derived step.

Bridging Gap 1 in this probe would mean DERIVING the cross-frequency closure — i.e., doing the open analytical work R76 §11 + R77 sketch §6 + T_lead_spectrum.md §6 all flag as open. That is reconstruction, not verification. The parent task is explicit:

> Phase 1 is gating just like in the prior R76 §11 probe. Don't reconstruct what isn't specified; only execute what is.

So the bridge is wide, not narrow.

**Alternate bridges considered and rejected:**

1. **Numerical fit from (P_+, P_−) data at k = 2..5.** This is what T_lead_2x2.py already does (numerical eigenvalues converging to ≈ 1/2). It is not a Q-derivation — it produces approximate floats, not Fraction entries. Decision rules (H_T_N_CONFIRMS_RATE_HALF vs H_T_N_DIFFERENT_EIGENVALUE) require Q. This bridge would force a degraded probe that the pre-registration explicitly excludes.

2. **Build T as a 6x6 over Q from level-n+1/level-n data via Markov stationary distributions.** Possible mechanically (build_markov_rational gives π_n; char_classes_exact gives all 8 entries of {P_n^{ab}(c)} over Q[ω]; finite difference between levels gives T entries that fit the recursion locally). But:
  - This is "fit empirical operator from a single transition", not "derive Off_lin from §5's procedure".
  - The fitted operator depends on the specific transition (k=2→3, k=3→4, etc.) — different transitions yield different "Off_lin"s, none of which is "the" Off_lin §5 describes.
  - The result is conjectural in the same way as the kappa_N flavor of result_77_2_T_N_construction.py — and that exact pattern was previously dispositioned as R77.3-falsified (see CANDIDATE_A_DISPOSITION.md + R_K_DISPOSITION.md).
  
  This is reconstruction. Rejected.

3. **Use one of the two heuristic derivations in T_lead_spectrum.md §3 (4·(1/4) = 1 OR P(v=1) = 1/2) as a chosen ansatz and proceed.** Pre-registration explicitly excludes choosing between competing heuristics:
  > Don't reconstruct what isn't specified; only execute what is.

  The two candidate derivations are mutually inconsistent in §3 itself ("4·(1/4) = 1? — no, more careful analysis needed"). Selecting one is invention, not verification.

---

## (f) Disposition: H_OFF_LIN_UNDERSPECIFIED

R77 sketch §5's articulation of Off_lin is at the level of:
1. A procedure (substitute Tao's recursion into the bilinear moment definition).
2. A claim that the result is "quadratic forms in {P_n^{ab}(c)}".
3. A statement that the operator acts on a 6-dim or 8-real-dim space (not 2-dim).

The procedure produces, at the v ≠ v' level, **cross-frequency bilinears that do NOT lie in span{P_n^{ab}(c)} without an additional closure step**. The claim "quadratic forms in {P_n^{ab}(c)}" is the open closure, not a derivation. T_lead_spectrum.md §3 attempts the closure heuristically and stops at two contradictory candidates with "more careful analysis needed". T_lead_2x2.py provides a numerical fit only.

**No 2x2 matrix over Q for Off_lin is specifiable from R77 sketch §5 alone.** Constructing one would require either:
- (Inventive) DERIVING the cross-frequency-to-same-frequency closure that §5 asserts but does not provide. This is the "off-diagonal exact bilinear-sum analysis" the project's own ledger (T_lead_spectrum.md §6) classifies as open work.
- (Degrading) Fitting Off_lin numerically from data, which produces approximations not Q values, and cannot decide H_T_N_CONFIRMS_RATE_HALF over Q.

Both options are excluded by the parent task's pre-registration gate.

**Disposition: H_OFF_LIN_UNDERSPECIFIED.** Phase 2 cannot proceed.

This is the exact failure mode the pre-registration anticipated:

> **H_OFF_LIN_UNDERSPECIFIED:** R77 sketch §5 describes Off_lin in terms that don't pin down a specific 2x2 matrix over Q. The probe blocks at Phase 1 because reconstructing Off_lin would mean inventing rather than verifying. Same failure mode as the R_K and R76 §11 probes.

---

## Adversarial check (Phase 1 A1 only — others gated on Phase 2)

**(A1) R77 sketch §5 fidelity.** The complete content of §5 is quoted verbatim in section (a). The judgment in (b)–(d) is based on:
- §5 itself (verbatim).
- T_lead_spectrum.md §3 (the project's own attempt to derive Off_lin's leading eigenvalue, with its own "more careful analysis needed" admission).
- T_lead_spectrum.md §6 ledger (project's classification of "off-diagonal exact bilinear-sum analysis" as Open).
- result_77_T_diagonal.py (the rigorous derivation of T_diag, which uses the diagonal v = v' collapse but does not extend to v ≠ v').

No reconstruction of R77 sketch §5 beyond restating what it says.

---

## Files referenced

- `result_77_sketch.md` §5 — primary source for Off_lin's articulation (procedural, not matrix-level)
- `result_77_T_lead_spectrum.md` §2, §3, §6 — companion document; explicitly flags off-diagonal as open
- `result_77_T_diagonal.py` — rigorous T_diag derivation (rests on v = v' collapse)
- `T_lead_2x2.py` — numerical 2x2 fit (data-driven, not derivation-driven)
- `result_77_2_T_N_construction.py` — 1-dim T_N flavor (kappa_N from data ratios)
- `R76_S11_DISPOSITION.md` — prior INCONCLUSIVE disposition that scheduled this Fork 1 probe
- `R_K_DISPOSITION.md` — prior H_R_K_INTRACTABLE; same reconstruction-vs-verification failure mode
