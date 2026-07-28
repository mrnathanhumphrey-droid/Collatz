# RESULT — SCOPING: S_∞ as a biextension height — the STRUCTURE has an abstract home, the HEIGHT (the value) is walled; and a live PORT lead (2026-07-28)

**Type:** literature scoping + structural assessment (no new compute; builds on the pen "S_∞ = height pairing" derivation).
Question: is S_∞ *literally* the archimedean height of a 𝔾_m-biextension, and does the named formalism close the open axiom
(the local–global functional equation / the value)? Verdict: **the pairing STRUCTURE (axioms 1–3) has a genuine abstract home;
the HEIGHT that attaches the real number 0.4749 (axiom 4) is walled by the same infinite-unipotent-depth obstruction — no
existing framework attaches a number to our q-difference / infinite-order object.** But a specific PORT is worth trying.

## The three papers, pulled and read
- **Nakayama, "Asymptotics of local height pairing," [arXiv:2512.22551](https://arxiv.org/abs/2512.22551) (Dec 2025).**
  Real, current, exactly on biextension height pairings — but **exclusively geometric**. The 𝔾_m-biextension (Bloch–Seibold)
  requires *two families of generically-homologically-trivial algebraic cycles on a smooth proper family X→S over a smooth
  curve*, codim m+n=d+1. Hain's biextension = a MHS with `Gr⁰=ℤ, Gr⁻¹=H, Gr⁻²=ℤ(1)`, a ℂ×-biextension over intermediate
  Jacobians; archimedean height = a real number via the real Hodge structure / Green currents. **No difference / q-difference /
  dynamical / abstract input anywhere.**
- **"Hardouin's biextension" is NOT in Nakayama** — it's an external citation `[Har05, Prop 4.4.1]`, used *geometrically* (to
  relate Hain's and Bloch's constructions, generalizing the Poincaré line bundle). `Har05` = **Charlotte Hardouin's 2005 PhD
  thesis, "Structure Galoisienne des extensions itérées de modules DIFFÉRENTIELS"** — iterated extensions of **differential**
  modules. **Differential, not q-difference** — the exact axis-mismatch P10 flagged. The earlier read ("named as the
  difference-equation generalization") was an over-reading; its origin is differential-Galois, its use in Nakayama is geometric.
- **Eskandari, "On blended extensions in filtered Tannakian categories," [arXiv:2307.15487](https://arxiv.org/abs/2307.15487).**
  The abstract home. Operates in **filtered ABELIAN categories** (v4: "an abelian rather than tannakian category with a weight
  filtration" — even weaker than Tannakian, so **no finite-rank obstruction at the structural level**). Explicitly treats the
  **unipotent** case (maximal unipotent radicals). BUT its output is a **homological classification of extensions — NOT a
  real-number height.** No regulator/period to ℝ/ℂ; no geometry needed for the structure, geometry needed only for applications.

## The split (this is the finding)
- **Axioms 1–3 (bilinear, symmetric, lifting-independent = the Annihilation Lemma) — HOMED.** The abstract biextension of a
  3-step weight filtration (Grothendieck SGA7 *extensions panachées*; Eskandari's filtered-abelian generality) genuinely houses
  our pairing. Our `M=D(I+N)` gives exactly the weight-filtered object. The hand-derivation was right about the *structure*, and
  it needs neither geometry nor finite-rank Tannakian input.
- **Axiom 4 / the HEIGHT (the value 0.4749) — WALLED.** Every framework that attaches an actual real number requires a
  realization our object lacks: geometric (cycles over a curve — Nakayama/Bloch, we don't have), differential-Hodge (Hardouin —
  wrong axis), or Tannakian-finite-rank (and S_∞ has **no finite-rank difference module at any order**, proven — `result_MAHLER.md`).
  Eskandari's abstract framework classifies the extension but **stops before the number**. The number-attaching functor is
  precisely the piece that does not transfer to the q-difference / infinite-order setting. Same wall as P8/P9/P10, in
  biextension language.

## THE PORT LEAD (worth trying — the first import that RETRODICTED a measurement)
Nakayama's paper is about **asymptotics of local height under UNIPOTENT degeneration** — structurally our `T_i → S_∞` tower,
with `N` the unipotent monodromy. The finite-rank obstruction only bites in the LIMIT, which is exactly what "asymptotics"
studies:
- **Each finite level is finite-rank** (`build_markov_q(3,i)` is a finite operator) ⟹ the biextension/height applies cleanly at
  level i; `T_i` = the local height at level i; `S_∞ = 2·lim T_i` = the asymptotic.
- **Geometric home exists** (P10 missed it): q-difference modules live on the Tate curve `ℂ*/q^ℤ` (Sauloy–Ramis); Mahler `z↦z³`
  is the geometry of the 3-power isogeny of 𝔾_m; the archimedean height there is the **Néron/theta local height** (the theta
  functions Holmes–de Jong use for the unipotent case).
- **Sharp retrodiction:** Nakayama Thm 0.1 = `h(t) − ⟨W,Z⟩·log|t|` bounded. Translate `log|t|` → level i: `T_i ~ ⟨W,Z⟩·i +
  bounded`. `T_i` **converges** ⟹ the non-archimedean Néron leading coeff `⟨W,Z⟩ = 0`, value = the bounded archimedean
  remainder. **This matches P9's independent finding that the 3-adic part is TAME** (`v₃(S_k)` bounded, value archimedean). The
  imported framework retrodicts our measured place-decomposition — the first route all session to do so, not just share a name.
- **Concrete gate + payoff:** is `T_i = ⟨C_ρ, K⟩` literally Nakayama's Green-current archimedean height (our `{2,−1,2}` kernel =
  his Green current)? If yes ⟹ a **closed form for `T_i`** (Néron/theta / Bernoulli) that we currently lack — a major gain
  *independent of rationality*. If no ⟹ the 𝔾_m-isogeny avatar is the wrong geometry.
- **Confidence flags:** SOLID = finite levels finite-rank; Nakayama = asymptotics-under-unipotent-degeneration; place-split
  matches measured 3-adic tameness. UNVERIFIED = that our renewal maps to the 𝔾_m-isogeny geometry and `T_i` is a theta height
  (Mahler `z↦z³` is the p-power *isogeny*, not a translation — the exact object may be solenoid/Drinfeld-flavored, not the clean
  Tate curve). CEILING = likely **re-expresses** rather than decides rationality (a limit of Néron heights being rational is its
  own hard question); the realistic win is the closed form for `T_i`, not a rationality proof.

## Net
- **Bank the split, not "literally the height pairing."** S_∞'s pairing IS a biextension in the abstract SGA7 / filtered-abelian
  sense (Eskandari homes axioms 1–3, no geometry/finite-rank needed). It is **NOT** established to be *the height* — the
  real-number realization (axiom 4) is walled, exactly by the infinite-unipotent-depth obstruction, in every framework with a
  name for it. Do NOT cite "Hardouin's biextension closes axiom 4": Hardouin is *differential*, Nakayama is *geometric*, and
  neither ingests our q-difference infinite-order object.
- **The live lead is the PORT** (cast the level-tower as a unipotent degeneration; ask if `T_i` is a Néron/theta local height on
  the 𝔾_m-isogeny geometry), because Nakayama's asymptotic framework retrodicts our measured 3-adic tameness. Next gate: check
  whether `⟨C_ρ, K⟩` is Nakayama's Green-current height at finite level (on our data), then a targeted read on the Mahler-`z↦z³`
  geometric realization.
- **Not at stake:** `result_MAHLER.md` (infinite-order, proven), GARSIA, DENOM, SOLSTICE, R1–R30. 7/15 unaffected (floor 0.473177).

Sources: [Nakayama arXiv:2512.22551](https://arxiv.org/abs/2512.22551) · [Eskandari arXiv:2307.15487](https://arxiv.org/abs/2307.15487) · [Hain–Eskandari, Jumps in the Archimedean Height, arXiv:1701.05527](https://arxiv.org/pdf/1701.05527) · Hardouin, "Structure Galoisienne des extensions itérées de modules différentiels," PhD thesis 2005 [Har05].
