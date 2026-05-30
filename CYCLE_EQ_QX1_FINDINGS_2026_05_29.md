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

## 10. Closed forms for the converged coset-mass constants: leading-order rational + non-elementary exact

Pursued closed forms for the limiting per-coset masses (q=17: 0.5450561…; q=31: 0.2284611/0.1510842/0.1204547).

**(a) Not elementary.** Geometric-truncation-independent (A_MAX 100=300); converged to ~7 digits (q=17
oscillates n=5/6 ⇒ Richardson limit of the asymmetry A_χ/A_triv → 0.0901126). Continued fractions don't
terminate, mpmath `identify`=None, and PSLQ rejects low-degree algebraic: `pslq[1,x,x²]`=`pslq[1,x,x²,x³]`
=None (height ≤1e5–1e6). So: not rational, not low-degree algebraic at available precision.

**(b) Leading-order rational mechanism (the structure that IS clean).** The masses are set by
A_χ = Σ_ξ χ(ξ)|μ̂(ξ)|² for the `index` Dirichlet characters χ trivial on ⟨2⟩ (E_c = (1/index)Σ_χ χ̄(c)·
A_χ/A_triv). The self-similar recursion X = 2^{−a}(1+qX') with χ(2)=1 forces the χ-value of the difference's
unit part to χ(1−2^Δ), Δ=b−a (iid Geom gaps), giving the leading-order closed form
  **V_χ = Σ_{Δ≠0} (2^{−|Δ|}/3) χ(1−2^Δ)**  (rational; periodic in Δ mod ord₂(q) ⇒ finite geometric sum).
- q=17 (χ=Legendre mod 17): **V = 76/765** = 0.0993464.
- q=31 (Z/6 characters): only the order-3 characters j=2,4 survive at leading order (V_1=V_3=V_5=0 — the
  conjugation character j=3 vanishes, consistent with −1∉⟨2⟩); leading levels 0.2527/0.1882/0.0591.
V_χ reproduces the QUALITATIVE structure (⟨2⟩-coset heaviest, correct #levels and ordering) for every q.

**(c) The gap = non-elementary corrections.** V misses the true value (q=17: 0.0993 vs 0.0901, ~9%;
q=31 low level 0.059 vs 0.120) because the deepening collisions (a=b, and a≠b with ord₂(q)∣Δ) carry shift
terms t=(2^{−Δ}−1)/q that the leading sum drops. The exact A_χ obeys a transfer recursion whose phase
e(−ζ(1−2^Δ)/q^n) does NOT close into a finite character system — it is the SAME boundary-layer object as
the char-fn peak (Arc A). Conclusion: the exact coset constants are **non-elementary**, with V_χ their
explicit leading-order rational shadow — directly parallel to Arc A's c_∞ (clean leading structure,
non-elementary exact amplitude).
Files: `probe_coset_closedform_2026_05_29.py`, `probe_coset_mechanism_2026_05_29.py`,
`probe_coset_precision_2026_05_29.py`, `probe_coset_q31_leading_2026_05_29.py`.

## 11. Collision-operator decomposition: λ = 1/3 universal, asym = (2 g(χ)/(q−3))·c_∞

Refining §10. The asymmetry decomposes into intrinsic objects via the collision transfer:

  A_χ / A_triv  →  g(χ) · c_∞ · (1 − λ) / (qλ − 1)

with **λ** the limiting q-adic collision rate (lim Σ_y P(y)² ratio across levels), **c_∞** the deep-collision χ̄-moment of the leading q-digit of (X−Y), and g(χ) the Gauss sum (=√17 for Legendre mod 17, since 17≡1 mod 4).

**(a) λ = 1/3 EXACTLY (universal).** From q=17 FFT at n=3,4,5,6:
  λ(3→4)=0.3333526461, λ(4→5)=0.3333341293, λ(5→6)=0.3333334049.
Each step ~17× closer to 1/3 (O(q⁻ⁿ) corrections vanish in the limit). Mechanism:
λ = P(a=b for two iid Geom(2)) = Σ_a 2^{−2a} = (1/4)/(1−1/4) = **1/3** — the dominant deepening event
at every q-adic level, with a≠b corrections (ord₂(q)|Δ) contributing only finite-n boundary terms.
This is universal: λ_∞=1/3 for every qx+1.

**(b) Clean rational factor 2/(q−3).** With λ=1/3, (1−λ)/(qλ−1) = (2/3)/((q−3)/3) = **2/(q−3)** —
explicit closed form. Combined with the Gauss sum:

  **asym(q,χ) = (2 g(χ) / (q − 3)) · c_∞.**

For q=17: asym = (2·√17 / 14) · c_∞ = (√17/7) · c_∞. Consistency: g(χ)·c·(1−λ)/(qλ−1) reproduces
the FFT asymmetry to 1e-16 at every n we checked (n=4,5,6).

**(c) c_∞ ≈ 0.152988994 — the irreducible non-elementary content.** Shanks-extrapolated c sequence
(c(n−1) at n=3,4,5,6: 0.15324792, 0.15300533, 0.15298871, 0.15298900). PSLQ on c against
[1, c, c², c³] = None at maxcoeff 1e6; identify = None. The c_∞ piece is non-elementary — same
boundary-layer character as the char-fn peak (Arc A's c_∞) — but **it is the ENTIRE non-elementary
content of the asymmetry**. Everything else around it is closed-form clean.

**Result.** §10 said "non-elementary." §11 says: **the non-elementary content is concentrated in a
single intrinsic constant c_∞ per character; structurally the asymmetry is asym = (2 g(χ)/(q−3))·c_∞,
fully clean except for c_∞.** Universal rate λ=1/3 and the clean rational factor 2/(q−3) hold for every
qx+1; c_∞ depends on q and χ. File: `probe_coset_collision_op_2026_05_29.py`.

## 12. Depth-resolved c(m): a rational ladder toward c_∞, with c(0) = 19/127 EXACT

Refining §11. From a single level-n FFT, c(m) at every depth m=0..n−1 is extractable via the
autocorrelation P_D = IFFT(|μ̂|²). Each c(m) is a finite weighted character sum ⇒ **rational at every
finite m** (computable in exact fraction arithmetic). Empirical extraction at q=17, n=6, plus exact
fraction computation:

| m | c(m) | denominator |
|---|---|---|
| 0 | **19/127** (exact) | 127 (= 2⁷−1, Mersenne) |
| 1 | 0.153178230055 | ~10³⁴ (rational, ladder explodes) |
| 2 | 0.153247920779 | rational |
| 3 | 0.153005331692 | rational |
| 4 | 0.152988709047 | rational |
| 5 | 0.152988999414 | rational |
| ∞ | ≈ 0.152988994 (Shanks) | limit of rationals with diverging denominators |

**c(0) = 19/127 closed form** (verified exact via probe_c_exact_rationals_2026_05_29.py): from §11
notation, c(0) = N(0)/T(0) with N(0) = V_χ = 76/765 (the leading-order sum from §10's V_χ mechanism!)
and T(0) = 508/765 (the depth-0 termination probability), giving c(0) = 76/508 = 19/127. So **V_χ from §10
is literally the numerator of the depth-0 rung**; the 9% gap between V_χ-prediction and the true asym is
exactly the gap c(0) ≠ c_∞ — the depth-correction that ripples through deeper rungs.

**c(1) is rational** with ~10³⁴-digit numerator and denominator (verified to 1e-13 vs FFT); subsequent
c(m) explode further. So the entire **rational ladder c(0), c(1), c(2), ..., c_∞** lives in Q with each
rung exactly computable, but denominators grow super-exponentially ⇒ the limit c_∞ is generically
non-elementary (matches PSLQ failures up to degree 4 / height 1e7 on the Shanks-extrapolated value).

**Final picture.** asym = (2 g(χ)/(q−3))·c_∞ with:
- **2/(q−3)** rational (from λ=1/3, §11)
- **g(χ)** Gauss sum (√q for q≡1 mod 4 + Legendre)
- **c_∞** = limit of an explicit rational ladder with anchor c(0) = 19/127 (= V_χ/T_χ from §10's V_χ)
  and depth-corrections that compound super-exponentially — the irreducible boundary-layer residual,
  same flavor as Arc A's c_∞.

For q=17 Legendre: asym = (√17/7)·c_∞, c(0) = 19/127, c_∞ ≈ 0.152988994.
File: `probe_c_inf_depth_extrap_2026_05_29.py`, `probe_c_exact_rationals_2026_05_29.py`.

**Negative result: no V/T-style miracle at c(1) (`probe_c1_structure_2026_05_30.py`).** Hunted for a
depth-1 analog of c(0)=19/127. Decomposing c(1) = (2¹³⁶X + Y)/(2¹³⁶X′ + Y′), the leading-order ratio
X/X′ matches FFT c(1) to 3e-13 with structure 1466497186566969389 / 9573796394196615981 — factorization
7·17·67·(15-digit) / 3⁶·(17-digit). The 17 (=q) and 3⁶ (deepening rate to power 6) are recognizable, but
there is **no analog of c(0)'s factor-of-4 reduction** that took 76/508 → 19/127. The 17-term weighted
sum at depth 1 (collapsed to 9 via the χ(−1)=+1 symmetry N(−s)=N(s)) does not telescope into a clean
small rational. ⇒ **c(0)=19/127 cleanliness is depth-0-specific** — at depth 0 the formula reduces to
N(0)/T(0) directly with no k-sum to integrate, and a coincidental gcd of 4 produces the small denominator.
At depth ≥ 1 the rational ladder is genuinely giant.

## 13. Wirsching ✗ + K=2 truncation washes out π; 17 closed-form mod-q deepening rates as byproduct

Tested two routes toward c_∞ via the deep-collision shift distribution π.

**(a) Wirsching's φ doesn't match π** (`probe_wirsching_check_2026_05_30.py`). Built φ_17 numerically as
the fixed point of W_17 f(x) = (17/16) ∫_{17x−16}^{17x} f (the qx+1 analog of Wirsching's predecessor
density operator). φ_17 is near-uniform on [0,1] (max 1.0625), giving near-Haar 17-digit distribution and
Legendre χ̄-moment **−0.0323**, vs c_∞ ≈ 0.1530 and c(0) = 0.1496. Self-convolution φ_17 * φ_17^rev also
fails (moment 0.0022). Wirsching's predecessor-density framework is structurally too smooth — it doesn't
crack π.

**(b) K=2 truncated kernel under Haar prior washes out π's structure** (`probe_aerial_dye_cameras_
2026_05_30.py` + `probe_truncation_diagnosis_2026_05_30.py`). Built the 289×289 deepening kernel
K(s, s′) on Z/q² with Haar marginal on higher digits. Leading eigenvalue came out to **1/q = 0.0588**, not
the true 1/3. Diagnostic confirms: row sums of K depend **only on s mod q** (std=0 within each mod-q
class). Under Haar prior the chain mixes to Haar in one step ⇒ leading eigenvalue = avg row sum = 1/q.
The true 1/3 survival rate requires the non-Haar higher-digit marginal — which is itself a property of π,
i.e. self-referential. Resolution requires (i) lazy q-adic simulation (state as digits-on-demand) or
(ii) self-consistent power iteration with on-demand higher-digit sampling.

**(c) Byproduct: 17 closed-form mod-q deepening rates.** The row sums by mod-q class are:

| r | rate | r | rate | r | rate |
|---|---|---|---|---|---|
| 0 | **0.193637** | 6 | 0.042361 | 12 | 0.056131 |
| 1 | 0.051782 | 7 | 0.039468 | 13 | 0.026117 |
| 2 | 0.074965 | 8 | 0.044920 | 14 | 0.067438 |
| 3 | 0.067438 | 9 | 0.044920 | 15 | 0.074965 |
| 4 | 0.026117 | 10 | 0.039468 | 16 | 0.051782 |
| 5 | 0.056131 | 11 | 0.042361 |   |   |

Each rate is **rate(r) = Σ_{δ ≡ −r mod q} P_D(δ)** for the explicit Δ = 2⁻ᵃX* − 2⁻ᵇY* distribution — a
finite weighted character sum, closed-form rational like c(0). The χ(−1)=+1 symmetry is directly visible
(rate(r) = rate(17−r) for all r∈[1,16]; rate(0) standalone). Mean = 1/q exactly. These 17 rates are
**conserved scalars of the chain** — water-test diagnostics that survive any finite truncation. They
don't determine π, but they constrain it. (For comparison: π's mod-q marginal is NOT equal to these rates
normalized — that would be the Haar π, which gives χ̄-moment 0, not 0.153.)

**Cumulative status of the deep-π chase:** §10 leading rational V_χ → §11 universal λ=1/3 + clean rational
2/(q−3) → §12 c(m) rational ladder, c(0)=19/127, no V/T miracle at depth ≥1 → §13 Wirsching ✗ + K=2-with-
Haar washes out π. The empirical/algebraic route is genuinely exhausted; cracking c_∞ requires either
(a) lazy q-adic Markov simulation (computational, gives arbitrary precision but no closed form), or
(b) genuinely new theory (Doob h-transform with explicit collision probability ψ, characterization of π
as identifiable measure). Both are real research projects.

Files: `probe_wirsching_check_2026_05_30.py`, `probe_aerial_dye_cameras_2026_05_30.py`,
`probe_truncation_diagnosis_2026_05_30.py`.

## 14. Doob h-transform: σ chain, multiplicative walk, π is uniform-on-⟨2⟩-coset (q=17)

The Doob route delivered a real structural identification of π for q=17.

**Setup.** D_m = X_m − Y_m for two iid Syracuse chains. The recursion
  D_{m+1} = (2⁻ᵃ − 2⁻ᵇ)(1 + qY_m) + q·2⁻ᵃ·D_m
gives, for collision at level m (D_m = q^m σ_m), the next-level collision condition:
v_q(2⁻ᵃ − 2⁻ᵇ) ≥ m+1. For m ≥ 1 this dominantly requires a = b (prob 1/3 = λ; subdominant
ord₂(q) | (a−b) corrections are O(q⁻ᵐ) and vanish in the limit). **In the dominant a=b regime,
σ_{m+1} = 2⁻ᵃ·σ_m — a multiplicative random walk on (Z/q)* by Geom(2)-distributed powers of 2.**

**Consequence: σ_∞ mod q is uniform on the ⟨2⟩-coset of σ_1.** The cumulative exponent
S_m = a_2 + … + a_m is a sum of iid Geom(2) draws; S_m mod ord₂(q) → uniform on Z/ord₂(q)
by ergodicity. So σ_m mod q = σ_1 · 2⁻ˢᵐ mod q → uniform on the orbit of σ_1 under ×⟨2⟩ =
σ_1's ⟨2⟩-coset.

**For q=17 specifically: −1 ∈ ⟨2⟩ (since 2⁴ = 16 = −1 mod 17).** Both ⟨2⟩-cosets are
negation-closed, and Legendre χ is constant on each (+1 on ⟨2⟩=QR, −1 on the other).
Therefore for ANY measure ρ on (Z/17)*:
  **E_ρ[Legendre] = 2·P_ρ(⟨2⟩) − 1.**
Combined with σ_m being uniform within its coset:
  **c_∞ = 2·P(σ_∞ ∈ ⟨2⟩) − 1.**

**Numerical verification** (`probe_doob_sigma_chain_2026_05_30.py`, 5M Monte Carlo pairs):

| depth m | n samples | χ̄-moment | P(unit-part ∈ ⟨2⟩) | 2P−1 | match |
|---|---|---|---|---|---|
| 1 | 1,120,424 | +0.154165 | 0.577082 | +0.154165 | ✓ exact |
| 2 | 372,736 | +0.154329 | 0.577165 | +0.154329 | ✓ exact |
| 3 | 124,755 | +0.152371 | 0.576185 | +0.152371 | ✓ exact |
| 4 | 41,350 | +0.149601 | 0.574800 | +0.149601 | ✓ exact |
| 5 | 13,935 | +0.144887 | 0.572443 | +0.144887 | ✓ exact |
| 6 | 4,590 | +0.166885 | 0.583442 | +0.166885 | ✓ exact |

The identity χ̄(coset) = 2·P − 1 holds to machine precision at every depth — confirming the
constancy-on-coset structure of Legendre for q=17 + that the σ chain is well-defined.

**Identification of π.** The deep-collision shift distribution π on Z_q (q=17) is:
  π = (uniform Haar on a ⟨2⟩-coset of Z/q)·(binary distribution on which coset).
That is, **π's only non-elementary content is a single binary probability P_∞ = P(σ_∞ ∈ ⟨2⟩) ≈ 0.5765.**
The "complicated measure on Z_q" reduces to a coin flip whose bias is non-elementary. The
empirical/algebraic π ladder rungs (P_0 = 73/127 at depth 0; P_m → P_∞ at depth m) all sit inside
this binary structure.

**Limits of the result.** (a) Specific to q with −1 ∈ ⟨2⟩ AND a quadratic character (q=17). For
q=7,23 (−1 ∉ ⟨2⟩) or q=31 (index 6, multi-character structure), the reduction is different:
π's content distributes over multiple cosets/characters. (b) P_∞ itself remains non-elementary;
the Doob framework rephrases the question (binary classification of σ_∞) but doesn't crack it.
Computing P_∞ exactly still requires either lazy q-adic simulation OR explicit construction of
the conditioning kernel.

**Where this leaves the arc.** §10–§12 algebra → §13 truncation washes out → §14 Doob reveals
π is uniform-on-coset with a single binary bit of non-elementary content. The empirical/structural
chase is now genuinely complete: we know π's *shape* (uniform on ⟨2⟩-coset), we know the
*mechanism* (multiplicative walk in the dominant collision regime), and we've identified the *one*
irreducible constant (the binary probability ≈ 0.5765). Cracking that one constant in closed form
remains the open theoretical problem.

File: `probe_doob_sigma_chain_2026_05_30.py`.

## 15. CORRECTION to §14 + spectral decomposition of σ chain

The structural claim in §14 — "in dominant a=b regime σ_{m+1} = 2⁻ᵃσ_m preserves coset ⇒ P_m = P_1
for m ≥ 1 modulo O(2⁻¹³⁶) corrections" — is **incomplete**. Direct test (`probe_metaphor_
experiments_2026_05_30.py`, experiment D): c(2) at FFT n=3 boundary vs n=6 interior gives
**identical** value 0.153247920779. So c(m) variation across m is **real structure, not
finite-n artifact**. Empirically:

| m | c(m) | c(m) − c(1) | comment |
|---|---|---|---|
| 0 | 19/127 = 0.149606 | −0.00357 | exact rational, boundary jump |
| 1 | 0.153178 | 0 | exact rational, 10³⁴ denom |
| 2 | 0.153248 | +7e-5 | ≠ c(1), real variation |
| 3 | 0.153005 | −2e-4 | |
| 4 | 0.152989 | −2e-4 | |
| 5 | 0.152989 | −2e-4 | converged toward c_∞ |

**What went wrong in §14:** the chain interior (X_m, Y_m) at depth m is fresh iid Syracuse from
new Geom sequences (a_{m+1}, a_{m+2}, …) — *independent* of σ_1's history. So σ_m's coset isn't
inherited from σ_1; it's set by the fresh interior pair's depth-conditional structure. P_m is a
*stationary* limit of independent depth-conditional shifts, not a propagation of a depth-1 boundary.
The closed-form Doob shortcut collapses; the actual picture is closer to "iid Syracuse difference
coset distribution at depth m, conditioned on collision-then-decollision pattern."

**Spectral decomposition (experiment A) IS clean and real.** The within-coset multiplicative
transfer operator T f(σ) = E_{a~Geom(1/2)}[f(2⁻ᵃσ)] on ⟨2⟩ ≅ Z/8 diagonalizes via characters
χ_k(σ) = ω^(k·log₂σ), ω = e^(2πi/8). Eigenvalues are closed-form:

  λ_k = E[χ_k(2⁻ᵃ)] = ω⁻ᵏ / (2 − ω⁻ᵏ).

Explicit values: λ_0 = 1 (Haar invariant); **λ_4 = −1/3 EXACTLY** (sign character on Z/8);
λ_{1,7} = ±i·conjugate pair with |λ| = 1/√(5 − 2√2) ≈ 0.679; λ_{2,6} with |λ| = 1/√5; λ_{3,5}
with |λ| = 1/√(5 + 2√2). **Coincidence note:** |λ_4| = 1/3 matches the empirical cross-depth
collision rate. These are formally different objects (within-coset mixing vs cross-depth survival)
but the same number — possibly a structural relationship worth investigating.

**Other experiments confirmed:** (B) temperature sweep p = 0.1, 0.3, 0.5, 0.7, 0.9 gives smooth
c(0)(p) curve from −0.13 to +0.80 through 19/127 at p=1/2. (C) chain is **irreversible**:
K(1,9)/K(9,1) = 64, detailed balance fails ⇒ P_∞ requires full transfer-op spectral structure.

**Where this leaves the π chase.** The §14 Doob framework gave a useful framing (the χ̄ = 2P−1
identity is a tautology of Legendre's coset-constancy) and identified λ_4 = −1/3 as the within-coset
spectral signature. But the simple "P_m constant for m ≥ 1" prediction failed empirically. The
actual mechanism for cross-coset transitions at depth m ≥ 1 is not the rare subdominant ord₂(q^m)|Δ
events I posited — it's the **fresh independence of (X_m, Y_m) from σ_1 at every depth**, which I
missed. P_∞ is genuinely a *limiting* stationary, not a depth-1 rational.

File: `probe_metaphor_experiments_2026_05_30.py`.

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
a single law confirmed for all primes q=3..47 (q=31 gives a 3-level graded ladder). §10 found a clean
leading-order rational shadow V_χ but the exact values are non-elementary. §11 sharpens that: the
asymmetry decomposes as asym = (2 g(χ)/(q−3))·c_∞ with **λ=1/3 universal**, isolating the
non-elementary content in a single intrinsic c_∞ per character.
