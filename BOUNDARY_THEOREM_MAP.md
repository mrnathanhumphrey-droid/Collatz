# BOUNDARY THEOREM — theorem → gate bridge

Maps the master write-up (`D:/Resolve Research/Collatz Documents/Papers/boundary_theorem_master.md`,
Thm 1–6) to this repo's verification record. Wilson's pen writes the theorems; these are the probes/gates that
judge them. Status key: **PROVEN** (transcribe from source) · **OWED** (pen derivation remaining) ·
**NUMERICAL** (pre-registered, commit-stamped). Live chronological log: `STATE.md`.

| Thm | Statement (short) | Status | Gate / probe | Result doc |
|---|---|---|---|---|
| **1** | Kinematic spectrum = {c_k}, c_k = Σw²ω^{kδ}; closed-form ℓ_k; rest nilpotent | PROVEN | Real-T1, gate 18/18 (commit `1fe5a64`) | result_phase2b_H / LALB |
| **2** | spec M = {c_k} ∪ {0} ∪ spec M_tower; partner = ρ(M_tower) | PROVEN | PARTNER-CHAR, G0 @ 1e-14 both L | result_phase2c0_* |
| **3** | Crossing: c₀(λ)=(1−λ)/(1+λ) moves, ρ(S)=1/3 weight-free ⟹ collide at λ=½ | PROVEN | G1 arc + **G-D0.3** gate (probe_trackD0_gate, λ=0.4/0.5/0.6) | result_trackD0_gate |
| | └ numerical companion: partner pinned 1/3 within 6e-4 while c₀ tracks | NUMERICAL | D-1 crossing rider (probe_trackD1_fork) | result_trackD_fork |
| **4** | Ladder **labeling** total {±1,±2,±4}+DC (626/626 L3); triangular selection 3^{L−1−j}\|κ; **NOT a reduction** (purity 0.07–0.32) | PROVEN (totality+rule); scope from measured purity | **J1** (probe_judge_J1) + 2c0-G3 | **result_judge.md §J1** |
| **5** | σ(θ)=(1/3)((1+e^{iθ})/2)²; bulk=½-flux, edge=Lebesgue restriction; dominant mode = d within couplings ≤3.6e-5 | PROVEN (sketch rigor); **Lemma D OWED** (closed-form d + march) | D2-a ½-flux (exact), D2-b/c edge, **D2-e** ladder matrices, **J2** couplings | result_trackD2a / result_trackD2e / result_judge.md §J2 |
| **6** | Coalescence at (3,½): condensation 0.681→0.900→0.987, phase→2π/3^{L−1}, gap (1/3)sin²(θ/2) @2% by L4; braiding; overlap→1/defect 17→189 | NUMERICAL spine complete; **rides on Lemma D** | G4 braid, D-1 EP recon, W4 witnesses, **J2** arm + **rider** (12-digit L4 doublet) | result_phase2d_G4 / result_trackD1 / result_W4 / result_judge.md §J2 |
| **7** | **Blindness:** spec(M₋)=spec(M₊) exactly — sign-blind; the distributional→pointwise barrier, located & proven | PROVEN (exact perm similarity); **lemma corrected** | **B1** — M₋=P M₊ P (pair-swap) EXACT 0.0; partner 0.34682666/0.33323630 reproduced | **result_blindness_B1.md** |

**Thm 7 mechanism note (walk-back #25):** the intertwiner is the **pair-swap** P:(a,b,γ)↦(b,a,γ), **not** the
committed negation Σ:(a,b,γ)↦(−a,−b,−γ). Σ fails (M₋≠ΣM₊Σ, 0.25–0.26) on carry-floor ⌊·⌋ vs modular-negation
non-commutation — the same breakage as the J-involution (walk-back #14); the durable constraint "invariance must
act on the carry as an INTEGER map, not modular" selects P (an integer relabel that sends T→−T = the +↔− flip).
The claim is unchanged and cleaner (exact permutation similarity).

## OWED (Wilson's pen — the two remaining derivations)
- **Lemma D** — closed-form d(gf, L) (exact lattice sum) with σ as its limit and the correction law reproducing
  the phase march 0.705 → 0.940 → 0.9928. Feeds Thm 5 limit + Thm 6 condensation.
- **Lemma B** — diagonal-dominance bound for dominant modes (empirically 2.5e-7 → 3.6e-5) from the selection
  structure. Feeds Thm 5's "dominant mode = d" claim.

## Judge (Probe J) — the verdict, unsealed 2026-07-18 (commit `d7692ab`)
- **J1 completeness (L=2,3, dense):** totality ✅ (all 626 modes → {±1,±2,±4}+DC, no 4th family); jury modes
  captured 1e-3→1e-2; **strong pre-reg FAILS** for subdominant (purity 0.07–0.32, resid→mode-scale) ⟹ ladder =
  complete labeling + dominant-mode theory, **not** a spectral reduction. (Thm 4's scope clause is this result.)
- **J2 arm (L=4, SpMV, no dense eig):** PASS ×3 — (i) |d(1,4)|/σ(θ₁)=1.00059 ∈(1,1.0201), dressing
  +2.01%→+0.059%; (ii) phase 0.9928 ∈(0.9434,1); (iii) couplings 3.6e-5 « L3. Cross-check: reproduces banked
  block-6 L4 doublet to 1.6e-4.
- **Rider (12-digit L=4 doublet split):** LANDED, **12 digits** (converged it 259, Δ=7.5e-14). Members
  0.320422712770+0.075242317692j & 0.320222549235+0.075251807019j; **split 2.00388e-4, ratio-to-L3 0.0758 ≈
  ×0.076** (block-splitting ledger confirmed). Partner 0.333499901322 = g4 ρ₄ to 12 digits. §6 entry complete.

## THREAD 3 — the constant (7/15, 7/45) and the loss ledger (Theorem F + the off-diagonal derivation)
The paper's constant, welded to the frozen instrument and its derivation reduced to a finitely-computable ledger.
Shell sequence S_k (= corpus per-scale L² cost, backward convention S_k = 3^k(a_k − a_{k−1}/3) = A_k − A_{k−1}).

| item | statement | status | gate / probe | result doc |
|---|---|---|---|---|
| **Theorem F (flat level)** | S_k → S∞ = 7/15 = 3·(7/45); three faces (spectral ρ=1/3, renewal driftless, Fourier non-decaying shells) | **welded** (S₁=2/3, S₂=10/21 exact; asymptote L-truncated at finite L) | **R1** re-weld + boundary contrast; **R2** A∞=(3/2)S∞=7/10 | result_thread3_R1 / result_thread3_R2 |
| A∞ = 7/10 (amplitude) | c₀-mode overlap of the independent pair | value holds; **ℓ₀-route DEAD** (2531/4095≠7/10) → it's the Jordan coupling (3/2)·β, not isolated ℓ₀ | **R2-B** (walk-back on crown route) | result_thread3_R2 §R2-B |
| **Theorem S (secular)** | S∞ = 3·g·φ_tow·ψ_kin; diverging g/Δ × vanishing Δ cancel | **mechanism gated** (Δ-cancellation L=3 ratio 0.974, L=4 0.9916); **L→∞ closed form OWED** | **R3** (L=4 product discarded: near-EP underconvergence) | result_thread3_R3 |
| edge-density law | 7/15 = L→∞ 1/θ edge residue (Dirichlet integral) | **NOT numerically confirmable** (finite-L band discrete, EP wall) — must close symbolically (Ĝ) | **R4** (scaling laws refuted at finite L) | result_thread3_R4 |
| dark-state selection | some band modes unread by the agreement functional | **readout zero** ⟨1\|r⟩=0 (DC-free), NOT a symmetry (no involution commutes; P dies #29) | **D1** (walk-back #29) | result_darkstate_D1 |
| deviation law | S_k = 7/15 + d_k | candidate 1/(5·21^{k−1}) **DEAD** at k=3 (wrong sign); true d_k signs +,+,−,−,−,− (overshoot) | **R5** (exact S₁..S₆ from Basic.lean) | result_deviation_R5 |
| **off-diagonal ledger** | S_k = S_{k−1} + OffDiag_k (diagonal replicates via 3×⅓=1); Σ_{k≥2} OffDiag = −1/5 | **R6-A GATE PASS**: OffDiag₂ = −4/21 derived from v≠v′ Ramanujan sums; diagonal replicates exactly k=2,3 | **R6** | result_offdiag_R6 |
| **channel engine (⟨4⟩-orbit law)** | OffDiag_k = (2/3)Σ_m 4^{−m}C_k(m), C_k(m) = twisted ⟨4ᵐ⟩-orbit character sum of μ_{k−1}, period 3^{k−1} | **R7 FULL GATE PASS k=2,3,4,5** (engine = frozen); sign-flip DERIVED (C_k(1) crosses zero); odd gaps ≡0; Mersenne denoms 4^P−1 | **R7** | result_engine_R7 |
| **uniform kill + strata** | correlation lives in μ non-uniformity (uniform μ ⟹ C_k≡0); OffDiag_k = (2/3)Σ_j W_j C̄_k(j), W_j k-independent | **R8 ALL PASS**: uniform kill C_k≡0 exact; band count m-independent (3^{k−1}, 2·3^{k−1}); ledger↔deviation weld exact; W_j closed form; overshoot = C̄_k(0) crosses zero k=3→4 | **R8** | result_strata_R8 |
| **collision identity (γ on ℤ₃)** | S_K = 2Σ_m 4^{−m}γ_{K−1}(τ_m), γ_n(τ)=collision density, γ_n(0)=X_n; whole campaign = one stationary γ_∞ | **R9 FULL PASS**: identity exact K=2..6 (from μ tables); γ_n(0)=X_n weld (τ=0 = qx+1 corpus); DC self-similarity C_k(DC)=3S_{k−1}; off-DC bounded; P_n=S_n, ⟨γ_n⟩=3S_{n+1}/2 | **R9** | result_gamma_R9 |
| **spectral form (character ledger)** | S = partial sum of frozen layers; OffDiag_{r+1}=2Λ_r, Λ_r=Σ_{ord χ=3^r}\|ν̂\|²/(4χ(4)−1); S_∞=7/15⟺Σ_{r≥1}Λ_r=−1/10 | **R10 FULL PASS**: 2Λ_r=OffDiag_{r+1} exact r=1..5 (char side only); layer mass=S_r weld; Mersenne trace 3^r/(4^{3^r}−1); r=2 excess 1.935×, r≥3 non-equidist. only | **R10** | result_charledger_R10 |
| **chirp kernel + two-campaign merger** | 2Λ_r=ε_{r+1}−ε_r (7/15 value ⟺ corpus ε-rate); Λ_r=⟨μ̂,K_rμ̂⟩, K_r=U*diag(w)U chirp-conjugated | **R11 A/C/D PASS**: U μ̂=ν̂ + shell block-diag (layer=shell mass = theorem of β); quadratic form=Λ_r; weld 2Λ_r=d_{r+1}−d_r exact, extended to new exact Λ₆; B: U sparse flat-on-support (not dense-flat); E: value=angle-half cancellation | **R11** | result_chirpkernel_R11 |
| **support law + closed loop + lobes** | U(k,ξ)=0 unless k≡ξ mod 3; angular moments = C-table entries (loop closed); Λ_r=L_r−M_r lobe difference | **R12 A/B PASS + F certified**: support law {k≡ξ mod3} exact; Σ\|θ̂\|²e(km)=C_{r+1}(m)/3 exact (no 5th coordinate); ε_k=d_k byte-equal (exact k≤8); Λ₁..Λ₇ exact; \|Re w\|-moment L_r+M_r stabilizes; classes balanced c₁=c₂ | **R12** | result_lobes_R12 |
| **lobe constant + ψ-existence + transport** | Re w=15/(2D)−½ (mean 0, ‖w‖²=1/15); chirp β = renewal (t′=β(2⁻ᵛ4ᵗ)); ψ = transport fixed point | **R13 D PASS, C→resolved by R14**: renewal-in-orbit gate exact n=1..5 (β IS the dynamics); ψ NOT uniform (osc ~2.5% deficit, depletion near trivial char); ψ-existence decider was open at r≤7, **RESOLVED in R14** | **R13** | result_psi_R13 |
| **deviation-field retarget (ψ-resolution)** | γ_n=1+Σ_{r≤n}A_r (A_r=C_{r+1}/3); bounded ⟹ no non-uniform ψ; Λ_r=S_r⟨δ_r,Re w⟩+Λ_r^unif | **R14 R13-C RESOLVED**: no non-uniform ψ (limiting-shape dead), object=deviation field δ_r; retarget identity exact; ⚠️Flag A (Λ^unif≠0 at r=2, 52%), Flag B (sign=two-band not near-0 lobe), ± reframe gated (odd lattice, no self-conj primitives) | **R14** | result_deviation_field_R14 |
| **endpoint split + bulk correlation** | Λ_r=S_r·b_r+Λ_r^unif; Λ^unif=R10-C trace (measure-free); bulk b_r = clean object | **R15 A PASS, B forced**: endpoint identity exact r=2..7 (anchors byte-exact); no self-conj primitive angle (odd) ⟹ Λ^unif=sampling residual; −1/10 is ~99% uniform baseline BUT this LOCALIZES difficulty into 1% not reduces (tail 2.59× answer, opp sign, 72% cancel); ⚠️#34 (5/3), #35 (trapezoid) killed | **R15** | result_endpoint_R15 |
| **transport recursion (the operator)** | θ_r = one step from μ_{r−1}+Geom; b_r = output of PASS-gated chain; crux = tower contraction | **R16 A PASS**: transport step exact r=2..6 (b_r not unexplained; b₁=0, no initial condition); dim δ_r=3^{r−1}−1. ⚠️its "q-sign oscillation q=3-critical" = R16-C **design defect** (off critical surface, superseded R23-A); ⚠️its "CRUX = SAME as R5 qx+1 step" **RETIRED #42** (R23) — R5 needs uniform decay, this is 94%-on-2-values localized | **R16** | result_transport_R16 |
| **slow mode (transport symbol + QSD)** | \|D\|²=1/(5−4cos(πξ/3^r)) linearized symbol; δ = quasi-stationary field | **R17 A/D**: closed form = linearized slow-mode symbol ((i)(ii)(iii) hold; exact transport non-diagonal, flagged); self-map+invariance PASS; ⟨\|D\|²⟩_δ≈1/3 flat (δ contracts at mean, QSD balance vs source); angular width holds (broad, not localized) | **R17** | result_slowmode_R17 |
| **regimes + roughness + branch + max-coeff** | settle ρ↔√(1/3) conflation; δ roughness; branch factorization; the Prop-1.17 quantity | **R18 A/C/D + B/E**: exact ratio ≈**1/2** not √(1/3), not geometric (osc @r=6); ρ≈0.988=signed-**envelope** decay = different object (my 0.984-as-√⅓ **corrected**); δ **equipartitioned across orders, NOT high-frequency** (adjective corrected, broad/non-decaying holds); **branch T=U₊D₊+U₋D₋ PASS** (DC 1/3,2/3; wts 1/15,4/15; #39 not incurred); **max/typical GROWS 1.34→4.25 ⟹ triangle/max bound DEAD** (Prop-1.17 gap concrete); R85 rung-1 exact-dead/float-cheap. **Crux refined: not norm contraction but summable fixed-m C-table decay (R12-B), max route confirmed dead** | **R18** | result_regimes_R18 |
| **typical or exceptional (the decider)** | is the additive spike the same obstruction as the fixed-m coefficients? | **R19 DECIDER (B): fixed low-harmonic A_r(m) m∈{1..4} are TYPICAL/depleted, O(1) osc, NOT growing** ⟹ equipartition/regularity route **SURVIVES R18-D**; within-stratum spikes real (C, 40×@r=7) on isolated members (refines R18-B to stratum-average). ⚠️R19-A "⟨2⟩/2-adic resonance" + R19-D "A-side orthogonal spike" **both corrected by R20** (vacuous/definitional); conclusion survives on B+R20 | **R19** | result_exceptional_R19 |
| **the thin window** | is the owed control uniform in m, or thin? + settle m=9; is the route one computation? | **R20 E PASS + A/B/C/D**: **Λ_r=Σ_m 4^{−m}A_r(m)=OffDiag/2 exact ⟹ deviation-field = R7 engine = collision-γ, ONE computation**; running sum **saturates by m≈3–4**, window O(1) ⟹ owed = summable bound on **thin O(1) window m≲r** (tail ≤(4/3)4^{−r}S_r); **m=9,27 OSCILLATE (R13-C read an oscillation) + weight-suppressed** ⟹ tension dissolved; additive slow mode migrates to **trivial char at (2/3)^r** (R19-A ⟨2⟩ retracted vacuous); **A_r(3^{r−1})=−S_r/2 VACATED** (R19-D artifact) | **R20** | result_window_R20 |
| **the ratio law (paper-abstract form)** | plainest reading of γ; independent gate; does f exist / is the weight on a typical part? | **R21 A GATE PASS + B/C/D**: **γ_r(τ_m)=3^r·ρ_r(4^{−m})**, ρ=law of X′/X by **group division** (independent of all machinery) reproduces γ all m + welds to **7/30**; theorem = **Σ 4^{−m}f(4^{−m})=7/30**, f=Haar ratio-density (#41 not incurred); **f ≈1 in bulk, SINGULAR at identity u=1 (=X_r→∞) but the 4^{−m} weight sits on the SMOOTH region (dev 0.11 vs bulk 0.44) ⟹ singularity weight-suppressed**; ⚠️R21-C **MISS** (argmax_7=2⁸ not 2⁷, peak near-degenerate — Prop 1.17 flat peak; additive thread closed by citation) | **R21** | result_ratio_R21 |
| **is f stratum-only** | can the theorem collapse to a stratum profile Σ_j W_j F(j)? | **R22 GATE FAILS — honest negative, §3 stratum-reduction VOID**: within-stratum std/mean stable ≈0.24 in j=0 (not →0); **the 3 j=0 orbit points γ(τ₁,τ₂,τ₄)→3 distinct limits ≈0.717/0.476/0.868 (spread 0.39)** ⟹ f u-dependent within strata. Diagnostics survive: **C** Haar sum=1 exact (binning ok); **B** means F(j) clean+stable (F(0)=2/3, diffs→7/15); **E** geometric weight **95% in j=0, ≈75%/19% on m=1/m=2** ⟹ theorem set by two orbit values f(τ₁),f(τ₂). **Confirms R20 thin-window, retires stratum-mean reduction** | **R22** | result_strata2_R22 |
| **critical family + f extrapolation** | is 7/15=M₄/M₃ a real closed form (q-sweep)? do f(τ₁),f(τ₂) lock? | **R23 TRIPLE NEGATIVE**: **Conj-2 FALSIFIED** — S_∞(q)≠(3q²+1)/(2q(q²+1)) at q=5,7,11,13 (builder validated vs q=3 byte-gate + Wilson's exact S₁; q=7 misses 5×) ⟹ 7/15=M₄/M₃ is a q=3 coincidence; **f(τ₁),f(τ₂) NOT locked at r=10** (γ(τ₁)↑0.723 spread 0.087; γ(τ₂)↓0.4731 **< 10/21** ⟹ **Conj-1 falsified**); closed-form hopes dead (no Mersenne, no 10/21, **functional eq proved non-existent**). ⚠️**#42 retires "same as R5"**; R16-C=design defect. **Owed = genuine convergence estimate, no shortcut** | **R23** | result_critical_family_R23 |
| **subcritical scaling law** | can 7/15 be reached from OUTSIDE criticality (positive-term, no cancellation)? | **R24 POSITIVE** (strongest evidence yet): critical sum is CONDITIONALLY convergent (no finite prefix determines it — the R23-B scatter was model misspecification); step off to λ=½+ε ⟹ clean positive geometric, exact rate **ρ=3(1−λ)/(1+λ) CONFIRMED**; **ε·X_∞→7/40 & amplitude (1−ρ)X_∞→7/15 both monotone** (0.602→0.471 @ε=0.1→0.01). ⚠️not 6-digit-locked (build wall r=10, small-ε tail-dominated); reframes target as **continuity of C(λ) at λ=½** (unproved). Trap avoided: extrapolate X_∞ not f(τ_m) | **R24** | result_subcritical_R24 |
| **the spectral gap (the gatekeeper)** | does the gap survive as ε→0 (⟺ C continuous at ½ ⟺ route closes)? | **R25** — drop X_∞, use **C(λ)=lim S_r/ρ^r** (ρ=exact leading eigenvalue); theorem = **boundary value C(½)=7/15**. Deep build r→14. **p₁₄(ε)→7/15 monotone** (supports it) but plateau unconverged small ε. ⭐**GATEKEEPER |λ₂|/ρ (1st look): HEALTHY large ε (0.69,0.57 @ε=0.1,0.05, real) but NARROWS toward ε→0** (plateau won't form @r=14, complex-pair oscillation emerges @ε≈0.02) — **leans |λ₂|/ρ→1 = gap shutting @½ = the period-9 pair**. UNRESOLVED (r=14 too shallow). Next: **direct transfer-operator 2nd eigenvalue** (finite matrix, no extrap) | **R25** | result_gap_R25 |

**Thread-3 net:** the constant's derivation is reduced to **Σ_{k≥2} OffDiag_k = −1/5**, a finitely-computable
exact-rational ledger. The **channel engine (R7)** now derives every term through k=5 from first principles:
OffDiag_k = (2/3)Σ_m 4^{−m}C_k(m), with C_k the twisted ⟨4ᵐ⟩-orbit character sum of μ_{k−1} (period 3^{k−1}),
engine == frozen for k=2,3,4,5. The two-sign tail mechanism is derived, not observed: the gap-2 sign flip
(−1/6 → +2/147) is **C_k(1) crossing zero** (C₂(1)=−1, C₃(1)=+4/49); odd-gap channels vanish identically
(conjugate-kill, verified); a positive DC/self-orbit class of collapsing Mersenne weight 1/(4^{3^{k−1}}−1) races
the negative bulk to produce the −,−,+,+ overshoot. The **uniform kill + strata (R8)** sharpen this: the whole
correlation lives in μ's non-uniformity (uniform μ ⟹ C_k ≡ 0 exact), the affine-band pair-count is m-independent
(3^{k−1}, 2·3^{k−1}), and OffDiag_k = (2/3)Σ_j W_j C̄_k(j) with **W_j k-independent** (closed form x/(1−x)−x³/(1−x³),
x=4^{−3^j}) — so the limit law is entirely C̄_∞(j), and the −,−,+,+ overshoot is **C̄_k(0) (the dominant W₀=20/63
stratum) crossing zero between k=3 and k=4**. The **collision identity (R9)** completes the reduction: the whole
campaign is **one stationary function γ_∞ on ℤ₃** — S_K = 2Σ_m 4^{−m}γ_{K−1}(τ_m) exact K=2…6, the τ=0 line welds
to the qx+1 corpus (γ_n(0)=X_n), DC self-similarity is literal (C_k(DC)=3S_{k−1}), and off-DC columns are bounded
(divergence confined to τ=0). The **spectral form (R10)** recasts the same content multiplicatively: S is a partial
sum of **frozen character layers** Λ_r = Σ_{ord(χ)=3^r}|ν̂(χ)|²/(4χ(4)−1), with OffDiag_{r+1}=2Λ_r exact (r=1…5),
the layer mass welding to S_r, and the Mersenne (4ᴺ−1) denominators exposed as the character-group trace of the
weight 1/(4ω−1). **Still owed (pen):** Σ_{r≥1}Λ_r = −1/10 in closed form — the within-layer distribution of
|ν̂(χ)|² over the χ(4)-angles (equivalently γ_∞(τ), the stationary C̄_∞(j), and Theorem S's L→∞ limit — one object,
three coordinates: additive γ on ℤ₃, stratum C̄, multiplicative Λ).
Walk-backs this arc: #25 (P not Σ), #26/#27 (shell convention — killed, backward pinned), #29 (dark = readout not
symmetry); crown ℓ₀-route retracted; deviation candidate 1/(5·21^{k−1}) killed.

## Verification ethos (§8)
Derive-blind / machine-gate / adjudicate; magnitudes pre-registered before the runs that judge them; the
walk-back ledger (~32 logged retractions, public; Thread-3 adds #31 dense-flat kernel, #34 5/3 end-lobe, #35
endpoint-atom/trapezoid, #42 "same object as R5's qx+1 step" (RETIRED corpus-wide R23 — R5 needs uniform-in-freq
decay, here 4^{−m} localizes 94% onto 2 values = structure R5 lacks), plus two self-corrections in R18 —
ρ≈0.984-as-√(1/3) conflation and the "high-frequency" adjective (δ is broad/equipartitioned); killed conjectures
(R23): Conj-1 f(τ₂)=10/21 and Conj-2 7/15=M₄/M₃ (critical-family q-sweep), and the ρ↦ρ′ functional-equation route
(proved non-existent) — #32/#33/#36/#37/#39/#41 pre-registered but NOT incurred); instrument law near the EP (dense/direct or within-block
power iteration; NO ARPACK/shift-invert; exact-rational; no rate-fitting). Honest negatives ARE the win —
the 27^{−L} law (44× miss), the ±4-seat withdrawal, J1's strong form, the crown ℓ₀-route (2531/4095), the
deviation candidate 1/(5·21^{k−1}), and R4's finite-L edge-density scalings are all in the record.
