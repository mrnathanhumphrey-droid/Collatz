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

## 7. qx+1 char-fn Plancherel mass does NOT see ord₂(q) — symmetry-protected null (RESOLVED)

The teed-up open question. Built the qx+1 Syracuse offset char-fn over Z/q^n
(X_n^(q) = Σ_{j} q^{j-1} 2^{-S_j} mod q^n — the 3→q generalization of the validated q=3 offset
builder), partitioned the units (Z/q^n)* by cosets of the cyclic subgroup H = ⟨2⟩, and measured the
Plancherel mass Σ|μ̂(ξ)|² per coset. Prediction: q=7's half-mixing (index-2 ⟨2⟩) would SUPPRESS the
dark (non-⟨2⟩) coset. **Falsified — exact null:**

| q | mixing | idx | mass in ⟨2⟩ | mass in dark coset | maxH/maxOut |
|---|---|---|---|---|---|
| 3 | full | 1 | 1.000000 | 0 | (no dark coset) |
| 5 | full | 1 | 1.000000 | 0 | (no dark coset) |
| 7 | half | 2 | **0.500000** | **0.500000** | **1.0000** |

For q=7 the split is EXACTLY 50/50 with identical per-coset maxima.

**Mechanism (verified to 1e-16).** The offset distribution is real ⇒ μ̂(−ξ) = conj(μ̂(ξ)) ⇒
|μ̂(−ξ)| = |μ̂(ξ)|. And −1 ∉ ⟨2⟩ mod 7 (6 ∉ {1,2,4}), so negation ξ→−ξ is a magnitude-preserving
BIJECTION between ⟨2⟩ and the complementary coset (verified: complement == −1·⟨2⟩ exactly). This forces
the 50/50 split independent of any Syracuse-specific structure. For q=3,5, −1 ∈ ⟨2⟩ (2≡−1 mod 3;
4=2² mod 5), so conjugate pairs stay inside the single full coset.

**Consequences.**
- The dark coset at q=7 is a PURE CONJUGATE REPLICA — zero independent magnitude information.
  Plancherel mass (a magnitude / 2nd-moment observable) is symmetry-protected and structurally cannot
  see ord₂(q). ~~**ord₂(q) is a PHASE phenomenon, not an amplitude one.**~~
  **[OVERGENERALIZED — REFINED IN §8: this holds only when −1 ∉ ⟨2⟩ (q=7,23). When −1 ∈ ⟨2⟩ (q=17),
  half-mixing IS amplitude-visible in Plancherel mass. ord₂(q) is NOT purely phase.]**
- Refines §6 / B6: the powers-of-2 transfer operator at q=7 reaches only ⟨2⟩, but since the complement
  is its conjugate mirror, every |μ̂| value also appears in ⟨2⟩ ⇒ the 1-D reduction stays
  **magnitude-complete** at q=7. So the q=3 transfer-op cleanliness (A2's exact 2.78e-17 FFT match) was
  NOT actually contingent on full mixing — it survives half mixing.
- The only amplitude-level trace of ord₂(q): the global argmax bounces between cosets for q=7 (its
  equal-magnitude conjugate twin lives in the other coset), vs trivially always in ⟨2⟩ for q=3,5.
- File: `probe_qx1_coset_plancherel_2026_05_29.py` (+ inline conjugation-symmetry check).

## 8. Refinement of §7: half-mixing IS amplitude-visible — q=7 was symmetry-protected (the −1∈⟨2⟩ switch)

§7's "ord₂(q) is a PHASE phenomenon, Plancherel mass can't see it" was **overgeneralized from q=7.**
The true control is whether **−1 ∈ ⟨2⟩** (verified `probe_qx1_neg1_coset_2026_05_29.py`):

| q | mixing | −1∈⟨2⟩ | mass ⟨2⟩ / dark | maxH/maxDark | reading |
|---|---|---|---|---|---|
| 11 | full (idx 1) | yes | 1.000 / — | — | no dark coset |
| 7  | half | NO (ord₂=3 odd) | 0.50000 / 0.50000 | 1.000 | symmetry-locked |
| 23 | half | NO (ord₂=11 odd) | 0.50000 / 0.50000 | 1.000 | symmetry-locked |
| 17 | half | YES (2⁴=−1 mod 17) | **0.5451 / 0.4549** | 1.16–1.20 | GENUINE asymmetry, mass-visible |

(q=17 split is stable across n=2,3,4 at ≈0.54506.)

**Mechanism.** Conjugation μ̂(−ξ)=conj(μ̂(ξ)) acts on the ⟨2⟩-cosets via multiplication by −1:
- **−1 ∉ ⟨2⟩** (q=7,23; ord₂ odd): −1 sends ⟨2⟩ to the *other* coset ⇒ dark coset = conjugate MIRROR of
  ⟨2⟩ ⇒ equal mass FORCED (50/50), and the dark coset carries **zero independent information** (pure
  replica). Half-mixing is invisible to mass AND there is no hidden phase signal — the effective frequency
  space is just ⟨2⟩ with conjugate-symmetric extension.
- **−1 ∈ ⟨2⟩** (q=17; ord₂ even, 2^{ord₂/2}=−1): −1 fixes EACH coset ⇒ no symmetry relates ⟨2⟩ to the
  dark coset ⇒ masses unconstrained ⇒ the powers-of-2 resonance genuinely concentrates in ⟨2⟩ and the
  dark coset is **suppressed (54.5 vs 45.5)** — half-mixing shows up directly in Plancherel mass.

**Corrected conclusion.** Half-mixing (ord₂(q)<φ) generically DOES register in amplitude/Plancherel mass
(q=17). The phantom "phase-only" reading came from q=7, where −1∉⟨2⟩ makes the dark coset a conjugate
replica that is mass-locked by reality of P. So there is **no hidden phase-only half-mixing signal to
chase**: for −1∉⟨2⟩ the dark coset is informationally empty; for −1∈⟨2⟩ it is already amplitude-visible.
File: `probe_qx1_neg1_coset_2026_05_29.py`.

## 9. Multi-coset generalization (q=31, index 6) + the general level-count LAW

q=31: ord₂ mod 31 = 5 ⇒ ⟨2⟩={1,2,4,8,16}, index 6 ⇒ 6 cosets; −1 ∉ ⟨2⟩ (ord₂ odd). Quotient
Q=(Z/31^n)*/⟨2⟩ ≅ Z/6, conjugation (×−1) = its order-2 element ⇒ the 6 cosets fall into 3 conjugate-equal
pairs. Result (`probe_qx1_multicoset_2026_05_29.py`), stable across n=2,3,4:

| pair | per-coset mass | pair total | reading |
|---|---|---|---|
| ⟨2⟩-pair (resonant, contains powers of 2) | 0.228461 | 0.456922 | heaviest |
| middle | 0.151084 | 0.302168 | |
| low | 0.120455 | 0.240910 | most suppressed |

Within each pair the two cosets are equal to 1e-17 (reality of P). THREE DISTINCT levels — a graded
hierarchy, ⟨2⟩-pair on top. (Per-coset masses are converged constants in n — intrinsic limiting coset
masses; no closed form yet, cf. q=17's 0.5451.)

**GENERAL LAW (derived + confirmed for every prime q=3..47).** The number of distinct Plancherel
coset-mass levels = the number of ORBITS of the conjugation involution (×−1) acting on the cyclic
quotient Q = (Z/q^n)*/⟨2⟩ (|Q| = index):

  #levels = index       if −1 ∈ ⟨2⟩   (conjugation trivial on Q ⇒ index fixed points)
          = index / 2   if −1 ∉ ⟨2⟩   (conjugation = order-2 element ⇒ index/2 swapped pairs)

Confirmed 14/14 primes: q=7→1, q=17→2, q=23→1, q=31→3, q=41→2, q=43→3, q=47→1, … The Syracuse dynamics
ACHIEVES the maximum allowed by symmetry — every conjugation-orbit class comes out distinct, no accidental
coincidences. **This unifies §7 (q=7: 1 level), §8 (q=17: 2 levels), §9 (q=31: 3 levels):** the ⟨2⟩
powers-of-2 resonance always tops a graded ladder whose length is fixed by how −1 sits relative to ⟨2⟩.
File: `probe_qx1_multicoset_2026_05_29.py` (+ inline q=3..47 law sweep).

## Net

The desktop translation is faithful and (after the off-by-one fix) correct, but canonical — no new
cycle bound. The productive yield is the **α_det closed-form family + its qx+1 generalization** (pole at
q=4=2^{E[v]}, q=3 unique converging odd q), the **clean θ form y^x=2y−1**, and the **universal
2-adic/q-adic asymmetry with ord₂(q) fine structure**. The teed-up open question — does the qx+1 char-fn
Plancherel mass reflect ord₂(q) mixing? — is **RESOLVED (§7 + §8)**: it depends on whether −1 ∈ ⟨2⟩.
For −1 ∉ ⟨2⟩ (q=7,23) the dark coset is the conjugate mirror, mass-locked to 50/50 and informationally
empty (§7). For −1 ∈ ⟨2⟩ (q=17) the dark coset is genuinely suppressed (54.5/45.5) — half-mixing IS
amplitude-visible (§8, correcting §7's "phase-only" overgeneralization). No hidden phase-only signal exists.
The full multi-coset picture (§9): the number of distinct Plancherel coset-mass levels = #orbits of (×−1)
on the quotient (Z/q^n)*/⟨2⟩ (= index if −1∈⟨2⟩, else index/2), with the ⟨2⟩-resonance always heaviest —
a single law confirmed for all primes q=3..47 (q=31 gives a 3-level graded ladder).
