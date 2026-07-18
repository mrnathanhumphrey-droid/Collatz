# FG Candidate H — Varopoulos-Saloff-Coste-Coulhon heat-kernel on locally compact groups

**Status:** UNVERIFIABLE-PHASE-0 — book paywalled, NOT in corpus (per operational note).

**Secondary citation reconstruction:** the only mention in the FG corpus is in Saloff-Coste 2004 bibliography (line 4480: "Astérisque 175, SMF" — reference to Saloff-Coste's earlier paper, NOT the Varopoulos-Saloff-Coste-Coulhon book proper).

The Varopoulos-Saloff-Coste-Coulhon framework's primary content (per general knowledge of the literature, not extractable from corpus):

- Heat-kernel p_t(x, y) on a locally compact (unimodular) group G with sub-Laplacian Δ from left-invariant Hörmander vector fields.
- Decay rate p_t(e, e) ~ t^{-D/2} where D is the volume-growth dimension of G (D = 0 for compact, D = n for nilpotent of step ≤ 2 with n-dim, etc.).
- Polynomial / exponential heat-kernel bounds on volume-growth + algebraic structure.

**For Syracuse on (Z/3^n)*:**
- The chain heat-kernel is the n-step transition probability K_n^t(x, y). For each fixed n, the state space is finite, so all heat-kernel decay theorems are trivially asymptotic-zero (finite-state ergodic).
- In the inverse-limit Z_3^*, the "heat kernel" would be the n → ∞ limit transition kernel. The volume growth of Z_3^* is **polynomial of dimension 0** (it's compact). VSC's machinery for polynomial-growth groups gives no quantitative content for compact groups.

Even *if* we had the VSC book statements verbatim, the heat-kernel framework targets **non-compact unimodular groups with polynomial or sub-exponential volume growth**, with the heat-kernel asymptotic being the load-bearing output. (Z/3^n)* and its inverse limit Z_3^* are **compact**, so VSC's results are vacuous (the heat kernel on a compact group equilibrates to Haar measure 1/|G|, no rate analysis needed in the VSC sense).

---

## Disposition H: **BLOCKER / UNVERIFIABLE-PHASE-0**, with structural prediction:

If statements were available, the most likely disposition would be **NO_FIT** on the compactness / volume-growth mismatch:

- h_H.group: locally compact + polynomial-or-exponential volume growth. Syracuse profinite is **compact** → vacuous case.
- h_H.heat_kernel: target object is p_t (continuous-time heat kernel from sub-Laplacian). Syracuse's K_n is discrete-time Markov, no canonical sub-Laplacian without additional structure.

**Phase 3 prediction: STRUCTURALLY_BLOCKED** on the compact-vs-noncompact category mismatch. Heat-kernel machinery is for "large groups, small heat-kernel diffusion"; Syracuse profinite is "small group, fast equilibration".

---

## Note: this is the "wrong" candidate per the brief's honest prior; VSC heat-kernel was the candidate most-flagged as risky for non-applicability, and the structural analysis confirms it.
