# Lagarias-class manifestations in q=3 bridge equation framework (Result 34)

**Purpose.** Canonical enumeration of every observable identified as
"Lagarias-class" or potentially-so across Results 1–33. For each item:
statement, reduction, open piece, bypasses, connections.

**Underlying classification.** "Lagarias-class" = the open piece reduces to
the **forward trajectory measure** on the Syracuse map (Lagarias 1985,
Tao 2022, equivalent to trajectory-measure invariance for the Collatz
3n+1 problem). Distinct from:
- Renewal-theoretic / fluctuation-theory open pieces (Path A, ladder mean E[L⁻])
- Finite-N empirical questions (asymptote tightness, double reversals)
- Cross-q questions (qx+1 universality of multiplier)
- Transcendence-class obstructions (Gelfond-Schneider on Path C)

## Headline count

**1 underlying Lagarias-class question** in q=3 framework:
forward trajectory measure on Syracuse map.

**4 manifestations** of that single question (A, D, F, K all reduce to it).

**1 likely-Lagarias-class adjacent** (B): reduces to a moderate-deviation
rate function under the same trajectory measure; not independently confirmed
to be a separate question.

**1 separate hard problem** (E): transcendence-class, not trajectory-measure.

**Other open pieces** (C, G, H, I, J): not Lagarias-class — finite-N,
renewal-theoretic, cross-q, or uncharacterized.

## A. Per-j W_j → ⟨σ_S | j⟩

**Statement.** P(j) closed form via Result 17 chain machinery; conditional
Wald W_j on entry classes m_j = (4^j−1)/3 reduces to ⟨σ_S | j⟩ via
conservation identity:

W_j = ⟨σ_S | j⟩ − (log N − 1)/log(4/3) − 1 + log(m_j)/log(4/3)

**Open piece.** ⟨σ_S | j⟩ closed form.

**Bypasses falsified (5 attempts):**
- Cramer-Esscher (gap +0.118)
- Geom-tilted Geom(1/2) (gap +0.197)
- Last-step inverse-density mod-3 mixing (gap +0.068)
- Uniform-mod-3 inverse-tree growth (gap −10)
- Single-eigenvalue spectral (Result 32, hump-shape pattern + λ_j non-invariant)

**Status.** Lagarias-class confirmed.

## B. σ-quantile-conditional w_q(q) closed form

**Statement.** Per-band Esscher tilt parameter w_q satisfies w_q ↔ E_band(q)
via the exact Esscher inversion (Result 25). 10-band measurement gave
w_q non-symmetric across q.

**Open piece.** Closed form for w_q(q) function.

**Bypasses falsified:**
- w_q ≈ c·z_q linear (Result 26, asymmetric pattern at 15 quantile points)
- Local CGF curvature (Result 28, doesn't flatten asymmetry)

**Status.** Likely-Lagarias-class via moderate-deviation rate function on
the same forward trajectory measure as A. Not formally confirmed as a
separately-Lagarias question. May reduce to A under deeper analysis
(rate function dependence on trajectory measure).

## C. K_full → K_h asymptote

**Statement.** Aggregate K_full(N) = quartile-weighted average of
K_eff_band(q, N). Per-band E_band closes via Esscher per-step + Result 32
algebraic Cov[T,V|band] correction.

**Open piece.** Precise asymptote (K_h exactly vs ~10.46) — requires N > 2⁴²
to discriminate empirically.

**Status.** **NOT Lagarias-class.** Per Result 33, boundary correction
closes within ±0.7 with no monotone trend; bulk closure is structural.
Asymptote tightness is finite-N empirical question, not trajectory-measure.

## D. ε(σ) per-orbit asymptote

**Statement.** ε(σ) = ε_S · log(6)/log(2) − ⟨α_det⟩, with
ε_S = Σ_j P(j)·[W_j − log(m_j)/log(4/3) + 1].

**Open piece.** Same as A (per-j W_j → ⟨σ_S | j⟩).

**Status.** **NOT a separate Lagarias-class piece** — different observable,
same underlying open question. Reduces to A.

## E. E[L⁻] closed form (i.i.d. ladder mean)

**Statement.** Strict descending ladder mean for Syracuse log-walk under iid
P(v=k) = 2^(−k). Empirical: E[L⁻] = 1.00456 ± 0.00061 nats (10⁷ orbits).

**Open piece.** Closed-form derivation.

**Bypasses falsified:**
- Rational-in-u Path C (Gelfond-Schneider obstruction via log₂ 3 transcendence)
- Esscher-duality identity (1480σ falsification)
- All algebraic candidates (e.g. (7/2)·log(4/3), at >7σ each)

**Status.** **NOT Lagarias-class.** Separate hard problem — transcendence-class
obstruction in standard fluctuation theory, not a trajectory-measure question.

## F. Inverse-tree depth distribution at m_j with mod-3 propagation

**Statement.** ⟨σ_S | j⟩ = expected depth of pre-image of m_j in the inverse
tree, weighted by mod-3 forward-propagation. Equivalent reformulation of A.

**Status.** **NOT a separate Lagarias-class piece** — equivalent to A.

## G. Closed form for b ≈ 2.275 and X-shape parameters (skew, kurt)

**Statement.** Per-band fit K_eff_band(q) = K_h + b·X_q where X_q is a
shape function of σ-distribution. Empirical b ≈ 2.275; X has specific
skew/kurt structure.

**Open piece.** Closed form for b and X-shape from Geom(1/2) odd-step
structure under Esscher tilt.

**Status.** Renewal-theoretic for Geom(1/2) under Esscher tilts; may close
within fluctuation-theory framework without trajectory-measure information.
**Probably NOT Lagarias-class.** Not yet investigated systematically.

## H. q=0.875 double-reversal mechanism

**Statement.** K_eff_band(q=0.875) double-transitions across N: settles at
plateau 14.63 (Result 23), with ξ_X / K_q875 decoupling documented.

**Open piece.** Mechanism for double-transition.

**Status.** **NOT Lagarias-class.** Probably finite-N artifact in band-edge
dynamics, not trajectory-measure question.

## I. qx+1 Cramér multiplier C ≈ 5/2 (q-independent)

**Statement.** conv_rate(j; q) ≈ A(q)·exp(−θ(q)·log(q)·j) with empirical
q-independent multiplier ≈ 5/2 in C across q ∈ {5, 7, 9}.

**Open piece.** Closed-form derivation; structural reason for q-independence.

**Status.** **Separate problem from q=3 Lagarias-class** — cross-q phenomenon
in qx+1 family, not q=3 trajectory measure.

## J. Per-octave β oscillation peak at j ≈ 21–22

**Statement.** E[v] correction explains baseline; per-octave variation
(~0.3 nats) does not.

**Open piece.** Mechanism for residual per-octave ~0.3 nat oscillation.

**Status.** Open; **not yet connected to Lagarias-class.** Could be finite-N
sampling structure or reflect deeper trajectory-measure modulation.

## K. Inverse-tree non-uniform residue stationary distribution

**Statement.** Stationary distribution of residues mod-3 in inverse-tree.

**Reduction.** Closed form via leading eigenvector (Result 23). Inverse-tree
measure is structurally distinct from forward trajectory measure
(anti-correlation r ≈ −0.20 at attractor structure).

**Status.** **CLOSED for inverse-tree.** The forward trajectory measure
(which IS the Lagarias-class object underlying A) remains open. So K's
"open piece" is precisely the same as A's.

## Summary table

| Item | Observable | Status | Reduces to |
|:----:|------------|:------:|:----------:|
| A | Per-j W_j → ⟨σ_S\|j⟩ | **Lagarias-class** | (root) |
| B | w_q(q) closed form | Likely-Lagarias adjacent | A (probable) |
| C | K_full → K_h asymptote | NOT Lagarias | finite-N |
| D | ε(σ) asymptote | NOT separate | A |
| E | E[L⁻] | NOT Lagarias | transcendence |
| F | Inverse-tree depth at m_j | NOT separate | A |
| G | b, X-shape constants | Probably NOT Lagarias | renewal theory |
| H | q=0.875 double-reversal | NOT Lagarias | finite-N |
| I | qx+1 multiplier 5/2 | NOT q=3 Lagarias | cross-q separate |
| J | Per-octave β peak | Unclassified | open |
| K | Inverse-tree residue dist | CLOSED (inverse) / A (forward) | A |

## Synthesis: 1 underlying Lagarias-class question

Items A, D, F, K-forward all reduce to the **forward trajectory measure**
on the Syracuse map. This is the single open question:

**What is the structural form of the forward trajectory measure on the
Syracuse map T(m) = (3m+1)/2^v_2(3m+1) restricted to convergent orbits?**

Equivalently: what is the distribution of the visit measure over residue
classes (or other natural partitions) along orbits, beyond the per-step
Esscher-tilt characterization (which captures only the marginal v_t
distribution)?

B (w_q closed form) likely reduces to a moderate-deviation rate function
under the same trajectory measure, but this connection has not been
formally established. If confirmed, B joins A/D/F/K as a manifestation of
the single underlying question. If independent, B is a second Lagarias-class
question.

Standard fluctuation theory + Esscher tilt + algebraic σ-identity (Result 32)
captures everything except trajectory-measure dependence at the residue-class
or sub-stratum level.

## Outside the q=3 framework

E (E[L⁻]) and I (qx+1 multiplier) are open but distinct from q=3
trajectory measure. They sit in different problem classes:
- E: transcendence in path-integral identities
- I: universality across qx+1 family (different problem class entirely)

## Implications for v3.6

- Bridge equation: 4 of 4 constants closed except for trajectory-measure
  manifestations of constant 3 (per-j W_j ⟺ ⟨σ_S|j⟩).
- Constant 4: bulk closes structurally; boundary closes within ±0.7
  (Result 33). Not Lagarias-class.
- Constants 1, 2 closed; constant 4 effectively closed; constant 3 reduces
  to trajectory-measure question (A).
- B remains as the one secondary potentially-Lagarias-class piece pending
  confirmation that w_q reduces to A.

## Implications for qx+1 cross-q comparison

q=3 baseline established: 1 (or 2) underlying Lagarias-class question with
4–5 manifestations. Same enumeration for q ∈ {5, 7, 9}:
- Item I (multiplier 5/2) is q-independent, suggests cross-q universality
- Per-q analogs of A (and possibly B) need separate measurement
- Renewal-theoretic items (E, G) may have q-independent character or not —
  open

## Files

- All Results 1–33 in `closed_form_findings.md`
- `route_b_v_orbit_test.md` (Result 32)
- `delta_k_band_decomposition.md` (Result 33)
- This catalog: `lagarias_class_catalog.md`
