# Cycle equation, α_det closed forms, and qx+1 generalization (2026-05-29)

Probe arc from translating a desktop-produced "cycle equation in the project machinery," then
pushing into α_det closed forms and the qx+1 generalization. All numbers verified in code.

## 1. The cycle equation — corrected convention (verified)

For a Syracuse cycle of L odd elements (m_{i+1} = (3m_i+1)/2^{k_i}, m_L = m_0):

  (2^K − 3^L) m_0 = Σ_{i=0}^{L-1} 3^{L-1-i} 2^{S_i},   S_i = k_0+···+k_{i-1},  S_0 = 0,  K = Σ k_i.

The exponent is **S_i = cumulative halvings strictly BEFORE step i** (S_0=0, max S_{L-1}=K−k_{L-1}),
NOT the through-step sum s_i = k_0+···+k_i (which wrongly injects a 2^K term).
- Verified: general telescoping identity `2^K m_L = 3^L m_0 + Σ 3^{L-1-i} 2^{S_i}` holds EXACTLY for
  12 structured+random trajectories. Trivial cycle {1,4,2,1} (L=1,k_0=2): (2²−3)·1 = 1 = 2^{S_0}=2^0. ✓
  (Desktop's s_i convention gave 4 ≠ 1 — the off-by-one.)
- File: `_cycle_eq_check.py`.

## 2. Literature check — the translation is CANONICAL (no new cycle leverage)

- Cycle equation: standard (Lagarias 1985 survey, Eliahou 1993). The corrected S_i convention IS the
  canonical one.
- The `(3/4)` rewrite (Route 1): canonical identity; `(3/4)=2²/3`, `log(4/3)` is the Lagarias-school
  per-step factor. An identity ⇒ no constraint by itself.
- Prefix-tree at mod 2^k (Route 2): KNOWN object ("residues mod 2^k evolve by a binary tree").
- Binomial j-distribution counting (Route 3): = Eliahou's composition count C(K−1,L−1) of admissible
  halving-vectors. Standard.
- **Binding obstruction is TRANSCENDENCE**: continued-fraction bound on log₂3 controlling |2^K−3^L|.
  Eliahou 1993: any nontrivial cycle has ≥ 17,087,915 elements. m-cycle bound (Simons–de Weger 2005 →
  Hercher): no m-cycle for m ≤ 91. The prefix/F̂_p machinery sits on the SAME Gelfond–Schneider wall
  (irrational log₂3) that blocks the project's closed forms — it can re-derive the counting, not sharpen
  the transcendence. So the translation adds nothing to the cycle bound.
- Project's own prior `result_cycle_obstruction.md` (2026-05-05) reached the same place (framework silent
  on cycles). NOTE: that doc's bound figures are garbled — Eliahou is 1.7×10⁷ (not 1.5×10⁸); Simons–de
  Weger is an m-cycle (circuit-count) bound, not "1.7×10¹⁰".

## 3. α_det closed forms (the genuinely productive thread)

α_det = log(6)/log(4/3) = 6.2282625189596… Equivalent closed forms (all exact to 60 digits):

  α_det = log_{4/3}(6) = (1 + log₂3)/(2 − log₂3) = log_{4/3}2 + log_{4/3}3 = (3·log_{4/3}3 + 1)/2.

- PSLQ: no independent integer relation (genuine transcendental ratio). CF = [6,4,2,1,1,1,2,71,3,1,168,…]
  (big partial quotients ⇒ no quadratic structure).
- **Cross-link:** the denominator `2 − log₂3 = 0.41504` is EXACTLY the depth-walk drift from the
  char-fn/cycle saddle calc (see `SYRAC_CHARFN_PEAK_DERIVATION.md`). Same constant, two roles.

## 4. qx+1 generalization

α_det^(q) = log(2q)/log(4/q) = (1 + log₂q)/(2 − log₂q). The "4" = 2^{E[v]}, E[v]=2 (Geom 1/2).

| q | drift 2−log₂q | α_det^(q) | regime |
|---|---|---|---|
| 3 | +0.415 | +6.228 | converge |
| 5 | −0.322 | −10.32 | diverge |
| 7 | −0.807 | −4.716 | diverge |
| 11 | −1.459 | −3.056 | diverge |

- **Pole at q = 4 = 2^{E[v]}**: α_det's singularity IS the converge/diverge boundary (q/4 < 1 ⇔ q=3).
- **q=3 is the unique odd q with α_det > 0.** For q≥5 it's negative (expansion signature).
- File: `probe_3adic_cycle_2026_05_29.py` and inline searches.

## 5. Project Cramér θ(q): clean form + critical point

θ(q) = root of `q^{−θ} = 2^{1−θ} − 1`. Substituting y = 2^{−θ}, x = log₂q gives the CLEAN form:

  **y^x = 2y − 1.**

- θ's nontrivial root (y<1) exists **iff x>2 (q>4)** — born exactly at the q=4 pole where α_det diverges.
  θ(3)=0 (only trivial root), θ(5)=0.349, θ(7)=0.626, θ(11)=0.804.
- θ and α_det share the q=4 critical point but are NOT related by a clean identity (PSLQ null;
  θ·(x−2), θ/(x−2), −α·θ all vary across q). Two independent order parameters for one transition.

## 6. Universal 2-adic / q-adic asymmetry (the "different mods" payoff)

The qx+1 Syracuse map n → (qn+1)/2^v:
- **Mod 2^k: single-valued** (v determined by residue) for ALL odd q — the deterministic prefix tree.
- **Mod q^L: multivalued** (frac = 1.000) for ALL odd q — intrinsically a transfer-operator/averaged
  object. 2 is the *division* prime; q is the *multiplication* prime. The split is structural, universal.
- **Finer:** q-adic mixing degree = ord₂(q^L). FULL (= φ(q^L)) iff 2 is a primitive root mod q:
  q=3, q=9 full; **q=7 is half (ord₂ mod 49 = 21 = φ/2)**. This retroactively explains why this session's
  3-adic transfer-operator/dlog machinery was so clean — q=3 is in the full-mixing class (2 primitive
  root mod 3ⁿ). For q=7 the qx+1 transfer operator carries a 2-fold coset structure.
- File: `probe_3adic_cycle_2026_05_29.py` + inline.

## Net

The desktop translation is faithful and (after the off-by-one fix) correct, but canonical — no new
cycle bound. The productive yield is the **α_det closed-form family + its qx+1 generalization** (pole at
q=4=2^{E[v]}, q=3 unique converging odd q), the **clean θ form y^x=2y−1**, and the **universal
2-adic/q-adic asymmetry with ord₂(q) fine structure**. Open/next: does the qx+1 char-fn Plancherel mass
reflect the ord₂(q) mixing degree (q=3 full vs q=7 half-coset)? — the F̂_p-saturation probe.
