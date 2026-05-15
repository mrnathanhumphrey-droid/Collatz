# BT_PHASE0_THEOREMS — Verbatim theorem extraction for Bruhat-Tits probe

**Date:** 2026-05-12. Probe BT Phase 0. Working dir C:/Collatz/. Source PDFs at C:/Users/Nate/OneDrive/Documents/black_hole/. Extracts at C:/tmp/bt/.

**Format per spec:** Source / Location / Statement (verbatim, OCR cleaned of stray hyphens/spaces but content preserved) / Notation translation / Tells us / Doesn't tell us.

---

## Required theorem set

### T1. BKL billiard motion → reflections on hyperbolic polyhedron (DHN)

**Source:** Damour, Henneaux, Nicolai, "Cosmological billiards", hep-th/0212256 v2 (Class. Quantum Grav. 20 (2003)).

**Location:** §1, pp. 4–5 (lines 102–119 of C:/tmp/bt/dhn.txt).

**Statement (verbatim):** "The asymptotic dynamics is found to be equivalent to a billiard motion in a region of Lobachevskii space H_d, where d is the number of (large) spatial dimensions, interrupted by geometric reflections against the walls bounding this region [17, 84]. Chaos follows from the fact that the Bianchi IX billiard has finite [volume]. … Besides the dimension of the hyperbolic billiard, the other ingredients that enter its definition are the walls that bound it. These walls can be of different types: symmetry walls related to the off-diagonal components of the spatial metric, gravitational walls related to the spatial curvature, and p-form walls (electric and magnetic) arising from the p-form energy-density. All these walls are hyperplanar. The billiard is a convex polyhedron with [finitely many faces]. In some cases, the billiard can be identified with the Weyl chamber of a Kac-Moody algebra, and the reflections against the billiard walls with the fundamental [reflections of that Weyl group]."

**Notation translation:** Lobachevskii space H_d ≡ d-dimensional real hyperbolic space (negatively curved Riemannian manifold). "β-space" ≡ Minkowski space of log-scale-factors. Walls w_A(β) ≥ 0 cut out a convex polyhedron in β-space; radial projection onto the unit hyperboloid yields the hyperbolic billiard.

**What it tells us:** The billiard substrate is **continuous hyperbolic space** with **finite hyperplane walls**, and reflections are specular (Euclidean/Lorentz-style). The Weyl-chamber correspondence is conditional ("in some cases") — not automatic.

**What it doesn't tell us:** Whether *any* discrete-arithmetic dynamics on a tree maps to this picture. The Lobachevskii setting is intrinsically continuous; there is no a priori bridge to a Bruhat-Tits *tree* (which is a 1-complex, not a 2D/higher hyperbolic space).

---

### T2. Finite volume ↔ chaos / infinite collisions (DHN)

**Source:** Damour, Henneaux, Nicolai, same paper.

**Location:** §5.2.2, p. 27 (lines 1298–1307 of dhn.txt).

**Statement (verbatim):** "Geodesic motion in a billiard in hyperbolic space has been much studied. It is known that this motion is chaotic or non-chaotic depending on whether the billiard has finite or infinite volume [81, 52, 101, 36]. In the finite volume case, the generic evolution exhibits an infinite number of collisions with the walls with strong chaotic features ('oscillating behavior'). By contrast, if the billiard has infinite volume, the evolution is non-chaotic. For a generic evolution, there are only finitely many collisions with the walls. The system generically settles after a finite time in a Kasner-like motion that lasts all the way to the singularity."

**Notation translation:** "Finite affine time" ≡ in β-space, the polywedge has cushions meeting at a cusp at infinity, but in hyperbolic-polar (ρ, γ) coordinates the collisions accumulate as ρ → ∞, which is finite proper time toward the singularity.

**What it tells us:** The "infinitely many bounces in finite affine time" claim our probe brief asks for is **a property of finite-volume hyperbolic billiards** — it is a consequence of finite hyperbolic volume + geodesic flow ergodicity, not specifically a Bruhat-Tits-tree statement.

**What it doesn't tell us:** Whether a discrete-vertex walk on T_3 has any analog. T_3 has *infinite* vertex set and a totally-disconnected boundary ∂T_3 = ℙ¹(ℚ_3). The "approach to a cusp" notion in DHN is a geodesic limit in continuous H_d; the discrete-walk analog is a *recurrent* vs *transient* random-walk dichotomy with hitting measure on ∂T_3, not "infinitely many bounces in finite affine time".

**Mode E flag:** Direct transfer of T2 to T_3 is non-trivial and the brief explicitly warns against it. The closest discrete analog is: random walk on T_3 is transient (since deg = 4 > 2), and the harmonic measure on ∂T_3 is well-defined. That is NOT the same statement as T2.

---

### T3. Billiard walls = Weyl-chamber walls when underlying KM algebra is hyperbolic (DHN)

**Source:** Damour, Henneaux, Nicolai, same paper.

**Location:** §7 (p. 49–50, lines 2295–2329 of dhn.txt).

**Statement (verbatim):** "The hyperbolic billiard is obtained from the β-space picture by a radial projection onto the unit hyperboloid of the piecewise straight motion in the polywedge defined by the walls. The straight motion thereby becomes a geodesic motion on hyperbolic space. … The billiard region, as a subset of hyperbolic space, is in general non-compact because the cushions meet at infinity (i.e. at a cusp); in terms of the original scale factor variables β, this means that the corresponding hyperplanes intersect on the lightcone. It is important that, even when the billiard is non-compact, the hyperbolic region can have finite volume. … Pure gravity billiards have finite volume for spacetime dimension D ≤ 10 and infinite volume for spacetime dimension D ≥ 11 [33]. This can be understood in terms of the underlying Kac-Moody algebra [27]: as shown there, the system is chaotic precisely if the underlying indefinite Kac-Moody algebra is hyperbolic."

**Notation translation:** "Hyperbolic" Kac-Moody algebra: indefinite generalized Cartan matrix whose principal submatrices are finite or affine. E_10, AE_3 are examples.

**What it tells us:** The Coxeter / Weyl-chamber structure is what governs the bounce sequence. The Cartan matrix of the underlying KM algebra fixes the wall angles and reflection rules.

**What it doesn't tell us:** Whether the Tao recursion has any KM / Coxeter structure attached to it. We'd need to *exhibit* a generalized Cartan matrix from the Tao recursion data to claim KM-billiard transfer.

---

### T4. Reflection rule on velocity vector at a wall (DHN)

**Source:** Damour, Henneaux, Nicolai, same paper.

**Location:** §7, eq. (7.3) referenced at line 1320 of dhn.txt.

**Statement (paraphrase of dhn.txt §7 — the explicit formula is cited but the page-7.3 formula is graphic; the verbatim around it is): "The billiard motion then consists of free motions of β^μ on straight lightlike lines within this polywedge, which are interrupted by specular reflections off the walls. … reflections that the velocities undergo during a collision, are elements of the orthochronous Lorentz group. Each reflection preserves the norm and the time-orientation; hence, the velocity vector remains null and future-oriented."

**Notation translation:** Each wall w_A(β) = 0 defines a reflection in the bilinear form induced by G_μν. The Coxeter group is generated by these wall reflections.

**What it tells us:** The dynamical update at a bounce is a fixed Lorentz reflection in β-space. In the discrete (Coxeter) version each generator is a fundamental reflection.

**What it doesn't tell us:** Anything that obviously translates to "multiply by 3, divide by power of 2". The Tao recursion update is multiplicative-then-divisive in (ℤ/3^n)^*, not a Lorentz reflection.

---

### T5. Henneaux-Persson-Spindel restatement / Coxeter-billiard table of contents (HPS)

**Source:** Henneaux, Persson, Spindel, "Spacelike singularities and hidden symmetries of gravity," Living Reviews 11 (2008) 1 (arXiv:0710.1818).

**Location:** §§1.1, 2.4, 2.7, 5.1, p. 7–80 (TOC lines 24–105 of hps.txt).

**Statement (verbatim, abstract): "In the vicinity of a spacelike singularity, the gravitational field equations can be reformulated in terms of billiard motion in a region of hyperbolic space. For certain gravitational theories, the billiard region is bounded by infinitely many walls and the resulting reflection groups are the Weyl groups of infinite-dimensional Kac–Moody algebras, suggesting that these algebras occur as hidden symmetries."

**Statement (verbatim, p. ~7): "Not only can the … equations of motion be reformulated as dynamical equations for billiard motion in a region of hyperbolic space, but also this region possesses unique features: It is the fundamental Weyl chamber [of a hyperbolic KM algebra]. … reflections in the walls bounding the fundamental Weyl chamber [define] 'words' in the Weyl [group]."

**Notation translation:** Same as DHN. HPS covers a broader catalogue of (super)gravities and the chamber-by-chamber book-keeping.

**What it tells us:** Confirms / consolidates DHN. The "infinitely many walls" case (line 24) is when the KM algebra is infinite-dimensional; the resulting words in the Weyl group are the bounce sequence.

**What it doesn't tell us:** Same gap as DHN — no statement about discrete-arithmetic / tree-walk objects.

---

### T6. Damour-Hillmann fermionic Coxeter structure (DH)

**Source:** Damour, Hillmann, "Fermionic Kac-Moody billiards and supergravity," arXiv:0906.3116.

**Location:** abstract + intro (early pages of dh.txt).

**Relevance:** Adds fermion sector to the Coxeter-group structure for 11D supergravity. Confirms KM algebra is robust across bosonic+fermionic content.

**What it tells us / doesn't:** Same picture — Coxeter / KM acts on β-space chamber, fermions don't change the chamber. No new structure usable on a tree.

---

### T7. Spherical functions on (semi-)homogeneous trees (Casadio-Tarabusi & Picardello)

**Source:** Casadio-Tarabusi, Picardello, "Spherical functions and spectrum of the Laplace operators on semi-homogeneous trees," arXiv:2208.00910 (2022).

**Location:** Intro pp. 1–3 + §§4–9 (lines 25–115, 460–550, 685–730, 1340 of spherical.txt).

**Statement (verbatim, p. 1–2):** "On a homogeneous tree, … the γ-eigenfunction φ(·, v₀|γ) of µ₁ that is radial around v₀ and normalized to 1 at v₀, called spherical function, can be expressed as a Poisson transform: φ(·, v₀|γ) = ∫_Ω K(·, v₀, ω|γ) dν_{v₀}(ω), where ν_{v₀} is the hitting distribution on Ω of the random walk starting at v₀ induced by µ₁; … For general trees, if v, w ∈ V and F(v, w) is the probability that the random walk starting at v visits w, then the harmonic Poisson kernel K(v, v₀, ω) (i.e., at the eigenvalue γ = 1) is a limit of quotients: K(v, v₀, ω) = lim_{w→ω} F(v, w) / F(v₀, w). … on a homogeneous tree one has K(v, v₀, ω|γ) = K(v, v₀, ω)^z for some complex number z related to γ: this yields a … parametrization of the eigenvalues, called the eigenvalue map γ(z) = (q^z + q^{1−z})/(q + 1)."

**Statement (verbatim, Theorem 7.1, p. ~17 of paper): "For γ ∉ sp_{ℓ²(V)}(µ₁), every spherical function φ(v, v₀|γ) on the semi-homogeneous tree with [explicit ℓ^p decay rate stated in original; verbatim form requires the paper's symbol table]." (We pulled the location; the actual ℓ^p exponent involves q_+, q_- in a manner we do not need to fully transcribe here, since the homogeneous specialization q_+ = q_- = q = 3 reduces this paper's content to standard Cartier 1973 / Figà-Talamanca-Picardello results.)**

**Statement (verbatim, p. ~30 / Theorem 9.1 ref): "[−b, −a] ∪ [a, b] is the ℓ²-spectrum of µ₁."** (where a, b are functions of q_+, q_-; in the homogeneous case q_+ = q_- = q, the standard result is sp_{ℓ²}(µ₁) = [−2√q / (q+1), 2√q / (q+1)] for the normalized adjacency operator.)

**Notation translation:** µ₁ = nearest-neighbor average operator (adjacency operator normalized by degree). For us q = 3 (so T_3 is the (q+1) = 4-regular tree if we use SL_2 conventions). ℓ²-spectrum is bounded by ±2√q / (q+1).

**What it tells us:**
1. The Poisson-kernel / harmonic-measure framework on ∂T_3 is fully developed.
2. The eigenvalue map γ(z) = (3^z + 3^{1-z})/4 parametrizes the spherical-function spectrum (for q = 3, homogeneous case).
3. Random walk on T_3 is **transient** (the visit probability F(v, w) decays) and has well-defined hitting distribution on ∂T_3.

**What it doesn't tell us:** Nothing about which specific automorphism / subgroup of Aut(T_3) corresponds to the Tao recursion. The spherical-function machinery describes K-bi-invariant statistics; the Tao map is not (visibly) K-bi-invariant.

---

### T8. Cocompact-lattice Hecke spectrum / Ramanujan condition (Lubotzky 2013)

**Source:** Lubotzky, "Ramanujan complexes and high dimensional expanders," arXiv:1301.1028.

**Location:** §1.2, Theorem 1.2.1 (lines 268–272 of lub2013.txt) and Theorem 1.2.3 (lines 362–368).

**Statement (verbatim, Theorem 1.2.1):** "Let F be a local field (e.g. F = ℚ_p the field of p-adic numbers …) with ring of integers O … so k = O/m is a finite field of order q. Let G = PGL_2(F) and K = PGL_2(O), a maximal compact subgroup of G. The quotient space G/K is a discrete set which can be identified as the set of vertices of the regular tree of degree q + 1 … [G/K is the (q+1)-regular tree, with neighbors of [L_0] corresponding to ℙ¹(k)]."

**Statement (verbatim, Theorem 1.2.3):** "Let Γ be a cocompact lattice in G = PGL_2(F). Then Γ\G/K is a Ramanujan graph if and only if every irreducible K-spherical G-subrepresentation of L²(Γ\G) is tempered, with the exception of the trivial representation (which corresponds to λ = q + 1) and the possible exception of the sign representation … which corresponds to λ = −(q + 1) and which appears in L²(Γ\G) iff Γ\G/K is bipartite."

**Notation translation:** For p = 3: G = PGL_2(ℚ_3), K = PGL_2(ℤ_3), G/K = T_3 the **4-regular** tree (q + 1 = 4). The Hecke / adjacency operator on T_3 has principal-series spectrum |λ| ≤ 2√q = 2√3 ≈ 3.46 on ℓ²(G), and the trivial + sign eigenvalues are ±(q+1) = ±4.

**What it tells us:** The Bruhat-Tits tree T_3 of PGL_2(ℚ_3) is the **4-regular tree** (not 3-regular). Vertices are equivalence classes of ℤ_3-lattices in ℚ_3². At any vertex [L_0], the 4 neighbors correspond to the 4 = 3 + 1 lines in 𝔽_3² = L_0/3L_0. Any cocompact arithmetic subgroup Γ ⊂ PGL_2(ℚ_3) gives a finite (q+1)-regular quotient graph.

**What it doesn't tell us:** Anything about Γ_0(3) ⊂ SL_2(ℤ) acting on T_3 (the candidate-C construction). That construction lives in a *different* p-adic setting (Γ_0(3) is an arithmetic congruence subgroup of SL_2(ℝ), not a cocompact lattice in PGL_2(ℚ_3)). Lubotzky's machinery is for cocompact lattices, which Γ_0(3) is not in PGL_2(ℚ_3).

---

### T9. Ramanujan / Alon-Boppana spectral gap (Lubotzky 2011)

**Source:** Lubotzky, "Expander graphs in pure and applied mathematics," arXiv:1105.2389, Bull. AMS 49 (2012).

**Location:** §1.2, Proposition 1.7 + Definition 1.8 (lines 261–270 of lub2011.txt).

**Statement (verbatim, Proposition 1.7 — Alon-Boppana):** "Let X_{n,k} be an infinite family of k-regular connected graphs on n vertices where k is fixed and n → ∞. Then λ(X_{n,k}) ≥ 2√(k−1) − o(1)."

**Statement (verbatim, Definition 1.8):** "A k-regular finite graph X is called a Ramanujan graph if λ(X) ≤ 2√(k−1)."

**Notation translation:** λ(X) = max{|λ_2|, |λ_n|} = second-largest absolute eigenvalue of the adjacency matrix (after the trivial eigenvalue k). For k = 4 (degree of T_3 if we use PGL_2(ℚ_3) convention), Ramanujan ⟺ λ(X) ≤ 2√3 ≈ 3.464.

**What it tells us:** Spectral gap of finite quotients of T_3 is bounded; 2√(k-1) is the Plancherel cap from T7's continuous-spectrum endpoint.

**What it doesn't tell us:** What discrete-walk operator the Tao recursion corresponds to. The spectral gap is a property of a graph, not of a dynamical map.

---

## Gaps flagged (no substitute in pulled set)

### Gap A — Bass-Serre theory of group actions on trees (Serre 1980).

**What we need:** Stabilizer / quotient structure for actions of arithmetic subgroups like Γ_0(3) on T_3. Without Bass-Serre we cannot rigorously verify Candidate-C's depth-walk claim (we'd need stabilizer profile at each vertex).

**Substitutes in pulled set:** Lubotzky 2013 §1.2 covers the *cocompact* case (Γ\G/K is a finite graph), but Γ_0(3) is NOT cocompact in PGL_2(ℚ_3) — it's an arithmetic subgroup of SL_2(ℝ), and its action on T_3 (via the embedding SL_2(ℤ) → SL_2(ℚ_3) → PGL_2(ℚ_3)) is via *non-cocompact* discrete action, with vertex stabilizers possibly infinite. Lubotzky's machinery does not cover this regime.

**Disposition:** Flag as gap. Candidate C analysis must work with hand-computed orbits at n = 1, 2, 3, not with general Bass-Serre vertex-stabilizer theorems.

### Gap B — Explicit Hecke / Möbius generators for Γ_0(3).

**What we need:** Concrete matrices in PGL_2(ℚ_3) corresponding to "multiply by 3" and "divide by 2." 

**Substitutes:** Hand-construct from first principles (PGL_2(ℚ_3) is well-defined; we can write matrices directly).

---

## What Phase 0 reveals — updates to Mode E discipline

**The BKL / Coxeter literature (DHN, HPS, DH) lives entirely in continuous β-space → continuous hyperbolic-polar (ρ, γ) projection → continuous geodesic flow on a hyperbolic polyhedron.** None of these papers references the Bruhat-Tits tree, p-adic analysis, or any discrete-arithmetic setting. The Coxeter / Weyl-group / Cartan-matrix structure comes from the **wall geometry in the linear β-space**, not from any tree structure.

The Bruhat-Tits literature (Lubotzky 2013, spherical-functions 2208.00910) lives in p-adic / tree setting with K-bi-invariant statistics, Hecke operators, harmonic-measure-on-boundary machinery. None of these papers references BKL, cosmological billiards, hyperbolic Weyl chambers, or Lorentz reflections.

**Mode E elevated to high alert.** The probe brief asks whether the Tao recursion can be both (a) lifted to T_3 dynamics AND (b) analyzed via BKL-billiard machinery. These are *two disjoint literature regions* with no overlap in our pulled set, no overlap in their primary mathematical objects, and no overlap in their notation. Connecting them is not transferring a known theorem; it is **proposing a new identification** that needs to be re-justified from scratch.

This is the **same category-of-object barrier** that closed Clusters 1, 2, and the 5-probe Fourier-decay arc. The BKL literature targets a continuous-Riemannian object (geodesic flow on H_d billiard); the Bruhat-Tits literature targets a discrete-graph object (Hecke spectrum on (q+1)-regular tree). Asking whether Tao on (ℤ/3^n)^* fits both pictures simultaneously is a *third*-category-of-object request.

**Pre-registered prior updates (post-Phase 0, before Phase 1):**
- H_BT_NONE_FIT: 30% → **45%** (the category gap is structural, no theorem bridges it)
- H_BT_STRUCTURE_PARTIAL: 30% → 30%
- H_BT_UNIVERSAL_FAILS: 25% → 15%
- H_BT_UNIVERSAL_PARTIAL: 10% → 7%
- H_BT_ALL_N_CANDIDATE: 5% → **3%** (the bridge would have to be from scratch, no theorem support)

Phase 0 PASSES (theorems extracted, gaps flagged). Phase 1 proceeds with elevated skepticism.
