# Collatz residue-class structural analysis

**Status (2026-07-21):** Live state in [`STATE.md`](STATE.md). The repo has accreted research threads on top of the original prefix-decomposition / Tao-bridge / qx+1-Cramér work — see "Recent landmarks" below.

## Table of Contents

- [Latest (2026-07-21) — the 7/15 constant campaign (Thread 3, R1–R10)](#latest-2026-07-21)
- [Latest (2026-07-16) — the qx+1 L3 spectral-gap campaign](#latest-2026-07-16)
- [Latest (2026-07-14) — fourth paper-shaped result](#latest-2026-07-14)
- [At a glance — where the project stands](#at-a-glance)
- [Recent landmarks](#recent-landmarks-post-2026-05-02)
- [Index by mathematical field](#index-by-mathematical-field)
- [TL;DR — what's here](#tldr--whats-here)
- [How to navigate this repo](#how-to-navigate-this-repo)
- [Documents](#documents)
- [Experiment index](#experiment-index-organized-by-theme)
- [Data files](#data-files)
- [Reproduction smoke check](#reproduction-smoke-check)
- [Open follow-ups](#open-follow-ups)
- [Outreach packages on Desktop](#outreach-packages-on-desktop)

---

<a id="latest-2026-07-21"></a>
> **Latest (2026-07-21) — the 7/15 constant campaign (Thread 3): the loss constant reduced to one stationary object, verified exact through k=6.**
> A parallel thread derives the marginal (q=3) loss constant **7/15 = 3·(7/45)** — the flat per-scale level of the shell sequence `S_k` (`S₁=2/3`, `S₂=10/21`, `S∞=7/15`). Wilson derives (pen); agents gate (pre-registered, exact-rational, no fitting). Live log [`STATE.md`](STATE.md); theorem→gate bridge [`BOUNDARY_THEOREM_MAP.md`](BOUNDARY_THEOREM_MAP.md). The arc R1–R10 walked the constant down to a **single finitely-computable object, seen in three coordinates**:
> - **The off-diagonal ledger (R6).** `S_k = S_{k−1} + OffDiag_k`; the diagonal channel replicates exactly (3-adic lift ×3 × 2-adic self-collision ×⅓ = 1), and `Σ_{k≥2} OffDiag_k = −1/5` is the constant's target. `OffDiag₂ = −4/21` **derived** from the v≠v′ Ramanujan sums, not read off the S-table.
> - **The channel engine (R7) — FULL GATE PASS k=2..5.** `OffDiag_k = (2/3)Σ_m 4^{−m}C_k(m)`, `C_k(m)` the twisted ⟨4ᵐ⟩-orbit character sum of `μ_{k−1}` (period `3^{k−1}`). Engine == frozen at every level; the gap-2 sign-flip (−1/6 → +2/147) is **derived** as `C_k(1)` crossing zero; odd-gap channels ≡ 0; the Mersenne `(4^P−1)` denominators arrive.
> - **Uniform kill + strata (R8) — ALL PASS.** Uniform `μ ⟹ C_k ≡ 0` exact (the correlation lives entirely in `μ`'s non-uniformity); the affine-band pair-count is m-independent (`3^{k−1}`, `2·3^{k−1}`); `OffDiag_k = (2/3)Σ_j W_j C̄_k(j)` with **k-independent** weights `W_j = x/(1−x)−x³/(1−x³)`, `x=4^{−3^j}` — so the whole limit law is `C̄_∞(j)`, and the −,−,+,+ overshoot is `C̄_k(0)` crossing zero between k=3 and k=4.
> - **The collision identity (R9) — FULL PASS.** `S_K = 2·Σ_m 4^{−m}γ_{K−1}(τ_m)` exact **K=2..6** (from the μ tables), `τ_m = (4^{−m}−1)/3 ∈ ℤ₃`, `γ_n` a collision density with `γ_n(0)=X_n` welding the τ=0 line to the qx+1 corpus; DC self-similarity `C_k(DC)=3·S_{k−1}` literal; off-DC columns bounded (divergence confined to τ=0). **The whole campaign is one stationary function `γ_∞` on ℤ₃.**
> - **The spectral form (R10) — FULL PASS.** Multiplicatively, S is a partial sum of a **frozen character series**: `OffDiag_{r+1} = 2Λ_r` exact **r=1..5**, `Λ_r = Σ_{ord(χ)=3^r}|ν̂(χ)|²/(4χ(4)−1)` (character side only), so `S_∞=7/15 ⟺ Σ_{r≥1}Λ_r = −1/10`. The order-`3^r` character mass welds to the shell mass `S_r` (a genuine cross-thread weld); the `(4ᴺ−1)` denominators are the character-group trace of the weight `3^r/(4^{3^r}−1)`; the within-layer excess is a hand-matched 1.935× at r=2, pure non-equidistribution from r=3 on.
> - **The chirp kernel + two-campaign merger (R11) — A/C/D PASS.** `2Λ_r = ε_{r+1}−ε_r` (ε_k = S_k−7/15) **welds the 7/15 value campaign to the corpus's old ε-rate campaign** (ρ≈0.984, period-9.2) — sum and tail of one frozen sequence. The chirp unitary `U(k,ξ)=3^{−r}Σ_z e((kβ(z)−ξz)/3^r)`, `β(z)=dlog₄(1+3z)`, transports `μ̂ → ν̂` exactly and is **shell block-diagonal** (so layer mass = shell mass is a theorem of the bijection), and `Λ_r = ⟨μ̂, K_rμ̂⟩` with `K_r = U*diag(w)U` — an additive quadratic form with the kernel written down. (Honest partial, walk-back #31: U is a *sparse* unitary, flat-on-support but not dense-flat — the "maximally-spreading" reading corrected. Angular: Λ_r is carried entirely by near-cancellation between angle-halves.)
> - **Support law + closed loop + lobes (R12) — A/B PASS, F certified.** The chirp support law `U(k,ξ)=0 unless k≡ξ mod 3` is exact (U = two conjugate dense-flat blocks split by the frequency's ± class); every angular moment of the profile is a banked C-table entry (`Σ|θ̂(k)|²e(km/3^r)=C_{r+1}(m)/3`), so the four representations (γ, strata, layers, kernel) **close on each other with no fifth coordinate**. The corpus's exact ε-table is certified byte-identical to `d_k=S_k−7/15` (exact through k=8), extending the exact ledger to **Λ₇** (signs `−,−,+,+,+,−,+`). The `|Re w|`-moment `L_r+M_r` stabilizes (~0.072) — one weak-* moment converges.
> - **Lobe constant + ψ-existence + transport (R13) — D PASS.** The weight is closed-form `Re w = 15/(2D)−½` (mean 0, `‖w‖²=1/15`, the Mersenne `4²−1` = denominator of 7/15). **The chirp β is the renewal itself**: `t′=β(2⁻ᵛ4ᵗ)` and `X_{n+1}=1+3·2⁻ᵛX_n` both reproduce the frozen measure exactly (R13-D). ψ is *not uniform* (an oscillating ~2.5% deficit, depletion concentrated near the trivial character x≈0).
> - **Deviation-field retarget (R14) — R13-C resolved.** γ_n(τ_m) is a partial sum of a frozen series (`A_r=C_{r+1}(m)/3`); boundedness (R9-C) forces `A_r` to have no nonzero limit, so **no non-uniform ψ exists** — the limiting-shape reading is ruled out, and the object is the deviation field `δ_r`. The exact retarget is `Σ_r [S_r⟨δ_r,Re w⟩ + Λ_r^unif] = −1/10`. Two honest flags corrected the session: the uniform-baseline term `Λ_r^unif` isn't zero at r=2 (52% of Λ₂), and the sign of Λ_r is a two-band (x≈0 vs x≈½) competition, not the near-zero lobe alone.
> - **Endpoint split + bulk correlation (R15) — A gate PASS, B forced.** The corrected retarget `Λ_r = S_r·b_r + Λ_r^unif` is exact (r=2…7), with `Λ_r^unif` the measure-free R10-C trace. The ± reframe is gated: since `N=3^r` is odd, **neither x=0 nor x=½ is a primitive angle** (no self-conjugate fixed point), so `Λ_r^unif` is a primitive-sampling residual, not a literal endpoint atom. The −1/10 is ~99% the r=1 uniform baseline (`Λ₁=−2/21`) — but that **localizes the difficulty into 1%, it does not reduce it**: the deviation-field bulk sum `Σ_{r≥2} S_r·b_r ≈ −9.85e−4` is a small residue of a large sign-alternating series (tail 2.59× the answer, opposite sign, 72% cancellation) — the *hardest* convergence to establish. (Walk-backs #34/#35: the 5/3 two-band and endpoint-atom mechanisms both died here.)
> - **Transport recursion (R16) — A gate PASS.** One transport step from `μ_{r−1}+Geom(2)` reproduces the next layer exactly (r=2…6), so **b_r is the output of a fully PASS-gated operator chain** (β=renewal, U certified, D explicit), not an unexplained sequence — and `δ₁=0` is forced, so the deviation field has no initial condition. The crux is now named and shared: the transport is a **tower map between growing spaces**, so the theorem is a **uniform contraction estimate on the tower — the same obstruction as R5's open step in the qx+1 paper. Two campaigns, one open problem.** (q-sweep: the sign oscillation is a q=3 criticality feature; q≠3 grows monotonically with no fixed point.)
> - **Slow mode (R17) — convention + quasi-stationary.** The transport's slow-mode symbol is `|D(ξ)|²=1/(5−4cos(πξ/3^r))` (the linearized shift; conjugation-symmetric, circle-average exactly 1/3 = Σ p_v²) — with the honest caveat that the *exact* transport is non-diagonal (it scales frequencies), so this is the leading symbol, not the exact operator. The deviation field is a **broad quasi-stationary field**: `⟨|D|²⟩_{δ_r} ≈ 1/3` (flat), so δ contracts at the mean rate 1/3 and is replenished by the source (‖δ‖²→7/15 steady state), with a fixed angular width (not localized). The theorem's difficulty is thus the *signed* phase residual (the bulk b_r), not the magnitude.
> - **Regimes + roughness (R18) — conflation settled, two adjectives corrected, one gate PASS.** The exact-regime term-ratio is **≈ 1/2** (below √(1/3)=0.577) and *not* geometric — it declines then begins a sign oscillation at r=6; the float-era `ρ≈0.988` is the **envelope decay of the signed sequence's oscillation**, a different object (a prior conflation of the two, corrected here). The deviation field is **equipartitioned across order-strata (≈1/r each, mildly coarse-loaded), not high-frequency** — the *broad, non-decaying* picture holds but the "high-frequency" adjective does not. The **branch factorization `T = U₊D₊ + U₋D₋` is certified** (v-parity split of the renewal; DC weights 1/3 & 2/3, symbol weights 1/15 & 4/15; #39 not incurred). And the concrete **Prop-1.17 quantity** — `max|μ̂_r|/typical` — **grows 1.34 → 4.25**, so a triangle-inequality/max-coefficient bound provably cannot close: the level-bound `r^{−A}` is the wrong species for the shell demand.
> - **Typical or exceptional (R19) — the decider.** The fixed low-harmonic coefficients `A_r(m)`, `m∈{1,2,3,4}` — the ones `b_r=⟨δ_r,Re w⟩` actually couples to, since Re w is smooth — are **typical or depleted (O(1), oscillating, not growing)**. So R18-D's additive max-death does **not** kill the equipartition / regularity route. (Within strata there *are* isolated spikes, 40× at r=7, refining R18-B's equipartition to a stratum-average.) *[Two R19 framings were corrected in R20: "pinned to ⟨2⟩ / 2-adic resonance" is vacuous — 2 is a primitive root, so ⟨2⟩ is the whole unit group; and the "A-side spike at `m*=3^{r−1}`" was the definitional `−S_r/2`. The decider result stands.]*
> - **The thin window (R20) — the route collapses to R7.** The weight is exactly geometric: `w(k)=1/(4χ_k(4)−1)=Σ_{m≥1}4^{−m}e(−mk/3^r)`, so **`Λ_r = Σ_{m≥1} 4^{−m}A_r(m) = OffDiag_{r+1}/2` exactly — the deviation-field route, the R7 channel engine, and the collision-γ identity are one computation, not three.** The running sum **saturates by m≈3–4** and the window coefficients are O(1), so the owed control is a **summable-in-r bound on a thin O(1) window `m≲r`** (tail `≤(4/3)4^{−r}S_r`). The m=9 tension is dissolved — `A_r(9)` turns over at r=6 and `A_r(27)` sign-alternates (R13-C sampled an oscillation), and both are weight-suppressed out of `Λ_r`. The additive slow mode migrates to the **trivial character** at `(2/3)^r`; `A_r(3^{r−1})=−S_r/2` is definitional (vacated).
> - **The ratio law (R21) — the paper-abstract form.** The plainest reading of γ: from `(1+3a)(1+3τ)≡(1+3a′)`, the ratio `X′/X ≡ 4^{−m}`, so **`γ_r(τ_m)=3^r·ρ_r(4^{−m})`** and the theorem is **`Σ_{m≥1} 4^{−m} f(4^{−m}) = 7/30`**, where `f` is the Haar density of the ratio of two i.i.d. Syracuse values. An **independent gate** — computing ρ by group division `X′·X^{−1}`, touching none of the prior machinery — reproduces γ for all m (r=2…6) and welds to 7/30. The density is `≈1` in the bulk but has a **growing spike at the identity `u=1` (`=X_r→∞`)** — yet the `4^{−m}` weight sits on the *smooth* part (deviation 0.11 vs bulk 0.44), so **the singularity is weight-suppressed out of the theorem.** (R21-C: a pre-registered miss — argmax₇ is 2⁸ not 2⁷; under Prop 1.17's superpolynomial decay the additive peak is near-degenerate, so the additive-peak thread is closed by citation.)
> - **Is f stratum-only? (R22) — no; an honest negative that sharpens the target.** A pre-registered gate on whether the theorem could collapse to a stratum profile `Σ_j W_j F(j)=7/30`. **It cannot:** the within-stratum spread of `f` is stable at ≈0.24 (not →0), and the three j=0 orbit points `γ(τ₁),γ(τ₂),γ(τ₄)` head to **three distinct limits (≈0.717, 0.476, 0.868)** — f is genuinely u-dependent within a stratum, so the stratum-reduction is **void.** The diagnostics survive and sharpen the target: the stratum *means* are clean (`F(0)=2/3`, differences → 7/15), and the geometric weight is **95% in j=0** with **≈75%/19% on m=1/m=2** — so the constant is set by two specific orbit values, `f(τ₁)≈0.717` and `f(τ₂)≈0.476`. This **confirms R20's thin-window reading and retires the stratum-mean route.**
> - **Critical family + f-extrapolation (R23) — a clean triple negative.** Two conjectured closed forms and one proposed method, all killed in one probe. **(1)** `S_∞(q) = (3q²+1)/(2q(q²+1))` (which would make `7/15 = M₄/M₃`, a ratio of renewal-step moments) is **falsified** by running the critical renewal at `q=5,7,11,13` — S_∞ does not track the smooth curve (q=7 misses 5×), so `7/15=M₄/M₃` is a q=3 coincidence (builder validated against the q=3 byte-gate *and* Wilson's exact `S₁(q=5),S₁(q=7)`). **(2)** `f(τ₂)=10/21` is **falsified** — pushed to r=10 in float, `γ(τ₂)` fell *through* 10/21 to 0.4731 and is still dropping; and `f(τ₁)` is still climbing (Aitken spread 0.087) — **neither dominant value is locked at r=10.** **(3)** the self-map functional-equation idea is **proved non-existent** (`ρ↦ρ′` isn't well-defined). Net: the closed-form hopes are dead, and — with walk-back **#42** retiring "same object as R5's step" — the honest state is that the owed step is a **genuine convergence estimate with no shortcut**, on a problem whose 94%-on-two-values localization is real structure Tao's uniform-decay problem lacks.
> - **Subcritical scaling → the gatekeeper (R24–R26).** Stepping off criticality (`λ=½+ε`, q held at 3) turns the shell sequence positive-geometric with **exact rate `ρ=3(1−λ)/(1+λ)`** confirmed; the scaling law `ε·X_∞→7/40` and amplitude `→7/15` recast the theorem as **continuity of `C(λ)=lim S_r/ρ^r` at the boundary λ=½** (`C(½)=7/15`). That continuity holds **iff the spectral gap survives** as ε→0, and Wilson's **`|λ₂|/ρ = 2λ²`** is confirmed to 4 digits (0.6841 vs 0.6840 at ε=0.1) — so `|λ₂|/ρ(½)≈½<1`, gap survives, route closes. (R25-C's "gap narrowing→1" reading was **retracted** in R26 as a small-ε artifact.)
> - **The recurrence + the gap operator (R27–R29).** The "exact closed recurrence" jackpot is **dead** (4 exact solves disagree ⟹ no finite recurrence), but the clean solves pin **`|λ₂|≈0.5` real** ({3,4}=0.503, = R18-A), so the dominant subdominant mode is **not** a near-circle complex pair; period-9 demoted to a minor `|λ|≈0.19` mode; ≥3 modes. Wilson's **gap-correlation operator (R27-E) is CERTIFIED (R28)** — `X_r=Σ_d P(d)R_{r−1}(d)` closes exactly (r=2…7), Diagonal Flatness exact, and the **odd-gap conjugate-kill is DERIVED** (odd-`d` correlations vanish to machine zero, `ord₃(2)=2`), giving R7's odd-gap vanishing a mechanism. But the bound is **mechanism-level, not rigorous**: channel mixing is O(1) (R28-C), and the finite gap **matrix** `3·diag(κ)K` **fails its sanity gate** (leading eig 1.08≠1, spurious mode; κ₂ drifts) — so `|λ₂|≈½` survives only on the converging numerical lines (R18-A, `2λ²`, R27, and R29-D's exact Λ-ratios straddling ½), not yet a rigorous finite eigenproblem.
> - **The digit probability (R30) — the cleanest statement, and a corrected target.** Wilson **retires his own √-cancellation requirement** (`max|μ̂|≲3^{−r/2}` was the triangle route, aimed 3 orders too high): the theorem needs only `Σ_r|A_r(m)|<∞`, so `O(r^{−1−δ})` suffices and the decay is observed geometric. The object recasts **exactly** as a conditional-digit probability — with `A_r(m)=γ_{r−1}(τ_m)(3q_r(m)−1)`, `q_r=p_r/p_{r−1}`, the theorem **is** `q_r(m)→1/3` (the next 3-adic digit of the ratio equidistributes). **R30-A GATE PASS:** the identity is certified exactly (m=1..4, r=2..7), a **new weld** joining the collision-γ ledger and the primitive-layer character coefficient `A_r(m)=Σ_{ord χ=3^r}|ν̂(χ)|²χ(4^{−m})` for the first time at m≥1 (the m=0 case is the R10-B/R28 weld `A_r(0)=S_r`). But the digit process is **not** the clean 3-state Markov chain that would make the gap free: the per-channel rate is neither a clean ½ nor m-independent (**B**, channel mixing), the excess over uniform **alternates** d=1↔d=2 by parity of r (**C**, the |λ₂| complex pair in digit space), and `q_r` depends on **which** digits matched (mod-27 spread 0.06–0.12) though that dependence **decays** with depth (**D**, asymptotically depth-Markov). All three point back to the same `|λ₂|≈½` gap, now wearing the guise of the homogenization rate of the conditional-digit distribution — the same open problem, in the regime coupling/renewal arguments are built for.
> - **M-reality promoted, and the constant retro-de-risked (R80).** The older `c=7/45` corpus (R74–R77) rests on `S_{n+1}=−2·M_{n+1}(1+3^n)` (Thm 76.3), whose proof had one empirical link — `Im M=0`, "verified numerically + class-symmetry." That link is now **Lemma 76.0**, an unconditional two-liner: the frequency set `{3∤ξ}` is `−1`-closed and fixed-point-free, and `π` real gives `μ̂(−ξ)=μ̂(ξ)*`, so reindexing `ξ→−ξ` yields `M(η)=M(η)*`. Needs only that `π` is a real measure — **not** class-symmetry (a red herring). Gate: `max|Im M|≤2.3×10⁻¹⁷` over all η to k=7. And the value never depended on it: conservation + Hermitian symmetry already give `S=−2·Re M(1+3^n)` unconditionally.
> - **The subordination warp (W1/W3).** A candidate derivation of `|λ₂|=½`: our random-scale renewal as the resolvent of a fixed-scale Bernoulli operator, dressed by the affine `+1` shift `A_*`. The **W3 gate is decisive and lands safe** — the unshifted core has *no* spectral gap (`|λ₂|→1`), the affine phase *creates* a stable gap (`|λ₂|=0.970`), so `A_*` is a **genuine spectrum-moving deformation, not a conjugation**: `½` is not an artifact and the four `|λ₂|≈½` lines are not reversed, with the shift confirmed load-bearing. But `½` is a **second-moment** quantity and does not appear in this first-moment picture (the `f(½)=1/3` step conflates the Bernoulli-convolution factor with the permutation resolvent), so the value-derivation must be built on the second-moment operator — the R29 growing-spaces object.
> - **What's owed (pen):** *sharpened* by R18–R20 — the whole program is the **R7 channel engine** `Λ_r=Σ_{m≥1}4^{−m}A_r(m)`, and the owed control is a **summable-in-r bound on the O(1) window coefficients `A_r(m)`, `m≲r`** — logarithmically thin, not uniform in m; the tail `m>r` is `≤(4/3)4^{−r}S_r` and the high-m modes oscillate + weight-suppress. R19-B and R20-A have measured the dominant part (`m≤4`) as O(1) through r=7; **conditional on that O(1) persisting, `Σ_r 4^{−m}A_r(m)` converges and `Λ_r` is summable — the theorem — for the weight-carrying harmonics** (measured, not proved). R23 killed the closed-form shortcuts (no Mersenne `M₄/M₃`, no `10/21`, no self-map functional equation) and **retired the "same as R5's qx+1 step" framing (#42)** — R5 needs *uniform* Fourier decay, whereas here 94% of the constant sits on two orbit values, real structure R5 lacks. So the owed step is a **genuine convergence estimate with no shortcut** (bulk sum ≈ −9.85e−4), and the honest open question is those two values `f(τ₁),f(τ₂)` — which R23 showed are **not even converged at r=10**. ~35 honest walk-backs logged (Thread-3 adds #31/#34/#35/#42, two R18 self-corrections, two R19 framings corrected in R20, killed conjectures Conj-1/Conj-2/functional-eq in R23; #32/#33/#36/#37/#39/#41 pre-registered but not incurred); the deviation candidate `1/(5·21^{k−1})` and the crown ℓ₀-route are both in the kill record.

<a id="latest-2026-07-16"></a>
> **Latest (2026-07-16) — the qx+1 paper's Result 1 enters its final step: a structured campaign on the spectral gap L3.**
> The standalone **qx+1 universal-rate paper** (`‖π_k‖² ~ C_q·3^{−k}`) has one open step — **L3: the spectral gap `r_q < 1` for `d = ord_q(2) ≥ 3`, with `r_q = 1` at the critical `d = 2` (⟺ q=3)**. It is now under a **6-phase, falsifier-first campaign** ([`STATE.md`](STATE.md), [`L3_DEFINITIONS.md`](L3_DEFINITIONS.md), [`PHASE1_WORKSHEET.md`](PHASE1_WORKSHEET.md)).
> - **Object frozen (Phase 0).** One `r_q`, three welded coordinates: **(M)** `build_M` pair operator (`r_q = |λ₂|/λ₁`), **(A)** renewal `A(z) = Σ S₀(i) zⁱ`, **(c)** `c_k = 3^k‖π_k − lift(π_{k−1})‖²`. Dictionary `c_k = (3/q)^k S₀(k)` and mass identity `Σ_k M^k v₀ = ‖π_k‖²` are exact; **`r_q` gate-validated `0.62` (q=5) / `0.39` (q=7)** against exact references (≤6e-10). The `3/q` confound (3/7 = 0.4286 vs r₇ = 0.39) raised and killed.
> - **Substrate proved (Phase 1).** Five lemmas — **FORGET / ONE-STEP / INTERTWINE / REFINE / PYTHAGORAS** — proved and machine-verified to ≤8e-16 (FORGET exact at **q = 1093**, the Wieferich prime: the substrate is s-blind). Pythagoras `X_k = (3/q)X_{k−1} + c_k` (R7/R42) is now *derived*, not measured.
> - **Boundary asset.** `build_M` is **genuinely defective at q = 3** (`cond(R) = ∞`, top-eigenvector overlap 0.998) — a real Jordan block, the exceptional point the d = 2 clause must break against.
> - **Platform placed (R41–R44).** Our `r_q = 1` boundary = Siegel's `σ_H = 1` degenerate map; the object is Tao's Syracuse random variable / Siegel's self-similar measure, read as a pipeline (Chang 2-adic input → Siegel transform → Tao 1st-moment / Nathan 2nd-moment `r_q`). The order-reciprocity "rig" refuted (a power-of-2 / Catalan coincidence).
> - **Honest walk-backs banked:** the single-address refinement operator was retired as the L3 object (first-moment ⇒ gapless — Gate G0 caught the drift); the 2/3-echo / 7/45-drip recursion died (no injection term).
> - **★ Phase 2 — two proven results (the entrance exam's kinematic half is closed).** **THEOREM D1** (the toy): `r(λ) = (1−λ²)/(1+λ²)`, maximality via nilpotence of the e=−1 block ([`BRIEF_D1_TOY_GAP.md`](BRIEF_D1_TOY_GAP.md), [`result_phase2b_Dmax.md`](results/result_phase2b_Dmax.md)). **THEOREM Real-T1** (the real q=3 operator — the program's *second proven result, first on the real operator*): the exact eigenvalues are the **twisted autocorrelations of the halving weights**, `c_k = R_k(0) = Σ_δ w_δ² ω^δ`, with closed-form left eigenvectors on zero-carry — gate-verified 18/18 at L=3 ([`result_phase2b_T1.md`](results/result_phase2b_T1.md)). The planned "invariant hunt" was **superseded** — the R1 STOP (no fishing) is what surfaced the closed form.
> - **Phase 2c — the limit theorem `ρ(M_tower,L)→1/3` (in progress, 2026-07-17).** The **dynamical partner** — the one object left in the entrance exam — is now pinned as `ρ(M_tower)`, the Perron of the carry-tower (γ≠0) principal submatrix, exactly (Probe G0). Its defect from the mean-field is frozen: an exact sector kernel `R_{k_in}(s)·N_{k_in−k_out}/D` (k=0 recovers the gate-confirmed E-FORM) and a **depth-selection rule** — resolving the carry to `mod 3^j` excites only the characters with `3^{L−1−j}\|k`, by exact orthogonality on the principal-unit tower (Probe 2c0). A **Collatz–Wielandt corrector chain** brackets the partner and shrinks: rung-1 `β*=3/5` (bracket `0.486 → 9/49`, kills the (odd,v₃=0) survival dichotomy), rung-2 the trit `τ∈ℤ/3` with the residual in closed form `r = 4/9 − κ(pop)·W⁻(τ)` — Wilson's blind derivation vindicated on all four pre-registered predictions ([`result_phase2c3_gate.md`](results/result_phase2c3_gate.md)). **Open — the crux:** the contraction estimate. Amplitude decay is ruled out (`|c_k|/c₀→1`; raw defect mass doesn't discriminate q3/q7), so the mechanism is the selection rule × the cascade tax `3^{−j}`. Then Phase 4 (the bound). USER-written, on the frozen object.

<a id="latest-2026-07-14"></a>
> **Latest (2026-07-14) — a fourth paper-shaped result, and the Tauberian thread closes.**
> - **★ NEW STANDALONE PAPER — qx+1 universal rate.** `S_k^(q) ~ (q/3)^k`, universal in q, **DERIVED at mechanism** ([`QX1_UNIVERSAL_RATE_WRITEUP_2026_07_14.md`](QX1_UNIVERSAL_RATE_WRITEUP_2026_07_14.md)). The "3" is named: `1/3 = Σ4^{-v} = E_{v~Geom(1/2)}[2^{-v}]`, the halving second-moment — **q-blind** (q enters only via the character / state-count `q^k`, never the halving statistic). Three pillars: rate, constant `c̃_q=(q-3)/q`, correction `δ_q≈0.82/ord_q(2)`. Adversarial-q falsifier ran FIRST and survived (small `ord_q(2)`, odd composites, `q≡0 mod 3`). **Scope: odd q only.** **Independent of every Collatz-closure thread below.** *Open — the single remaining line:* uniform diagonal-self-overlap domination on `(Z/q^k)*` (R76 conservation generalized from q=3) to upgrade mechanism → theorem. Reported as proved-at-mechanism, **not** as theorem.
> - **c = 7/45 is now provenance-hole-free.** The load-bearing 1:4 squared-class-mass ratio — previously asserted without a result file — is reconstructed and **CONFIRMED elementary + exact** (`result_64B.md`): class ± = v-parity ⇒ `P(v even)=1/3` ⇒ `(1/3)²:(2/3)²=1:4`. A one-line Geom(½) parity identity. The R75+R76+R77+R64.B chain has no remaining gap.
> - **Tauberian/BGT corpus FULLY DISPOSED.** The last surviving PARTIAL (Bingham–Ostaszewski candidate E) closes **NO_FIT**: `L(k)=|ε_k|·2^k` does not reach a second plateau — for k≥10 it grows **geometrically** toward `2·0.984=1.968`, which is outside the regular-variation framework entirely. Three regimes (plateau ≤6 / period-9 transition 7–9 / escalation ≥10); **multi-regime obstruction confirmed STRUCTURAL.** Byproduct: independently re-confirms the subdominant rate ρ≈0.984.
> - **Four clean negatives** (structural priors 0-for-5 this arc — the falsifier-first protocol is why): inverse-tree 1/9 **dead in both measures — do not re-add**; `7=Φ_p(2)` refuted (2-point coincidence, breaks ~0.65× at p≥7); the mod-9 offset is a normalization artifact; `v₃(c_k)` closed form NULL.
> - **Owed:** the R85 rung-1 operator-DFT chirp identity (a positive) is **single-r (n=6) — owed an r=5 (n=8) extension before it is trusted.** Same small-window profile this arc has killed three times.

---

<a id="at-a-glance"></a>
## AT A GLANCE — where the project stands (2026-07-16)

> **qx+1 paper (Result 4), final step in progress:** its universal rate `‖π_k‖² ~ C_q·3^{−k}` is proved at mechanism; the last open step is the **spectral gap L3** (`r_q < 1` for `d ≥ 3`), now a 6-phase campaign with **Phases 0–1 closed** (object frozen + welded, `r_q` gate-validated 0.62/0.39, five substrate lemmas proved). See the 2026-07-16 "Latest" block above, [`STATE.md`](STATE.md), [`L3_DEFINITIONS.md`](L3_DEFINITIONS.md), [`PHASE1_WORKSHEET.md`](PHASE1_WORKSHEET.md).

**Four paper-shaped results in hand** (the fourth, qx+1, is standalone — see "Latest" above):

1. **Leading c = 7/45 RIGOROUS UNCONDITIONAL** ([`THEOREM_C_745.md`](THEOREM_C_745.md)). `S_k = 3^k · ‖d_k‖² → 7/15` proved via R75 Plancherel × R76 conservation × R77 T_diag × R64.B class-mass × HR74 algebraic identity. Paper-shaped. Independent of all operator-valued probability framework questions.

2. **Syracuse = Davies-Wiseman-Milburn quantum trajectory, NUMERICALLY VERIFIED to 6 sig digits** ([`FRAMEWORK_IDENTIFICATION.md`](notes/FRAMEWORK_IDENTIFICATION.md), [`DWM_MP_G1_RESULT.md`](notes/DWM_MP_G1_RESULT.md)). DWM cross-Kraus form `M̃_{v,v'}^{(j, b_prior)} · f(ξ) = phase·σ_{-(v+v')}·f(ξ)` reproduces Syracuse's measured moments **exactly to 6 sig digits across all 4 scalar reductions** for both 3-alternating (0.108) and 4-alternating (0.609) moments.
3. **Dark-subspace classification of Syracuse's adaptive Kraus family** ([`R3_DARK_SUBSPACE_STRUCTURAL.md`](notes/R3_DARK_SUBSPACE_STRUCTURAL.md), integrates Phases 1+2+4). Three sub-results: (a) the full adaptive Kraus family at level n is irreducible at finite n (dim(A') = 1 at n=2, 3, via SVD of commutator-stack); (b) the **TRUE D_W = 3-fiber-zero-mean subspace is EXACTLY dark under the j ≥ 2 sub-family** (machine-epsilon leakage at n=2, 3, 4 verified directly), with j = 1 the unique mixing event — structurally forced by `x_{j≥2} ≡ 0 mod 9` vs `x_1 ≡ ±1 mod 3`; (c) the per-step channel `L|_{D_W}` for j ≥ 2 has **closed-form below-commutant spectrum** `λ_below(n) = 0.5/|1 − 0.5·e^{iπ/3^{n−1}}|`, **j-independent for j ≥ 2**, **verified across 8-point grid (n, j) ∈ {(3,2), (3,3), (4,2), (4,3), (4,4), (5,2), (5,3), (5,4)} via matrix-free ARPACK** + full-eigvals cross-check at n=5 (rel err < 1.3e-5 throughout). λ_below(n) → 1 as n → ∞. **Route B partial closure (2026-05-16)** ([`ROUTEB_PERIOD9_IDENTIFICATION.md`](notes/ROUTEB_PERIOD9_IDENTIFICATION.md)): the (class, b_prior mod M) Markov chain at M=18 reproduces the same eigenvalue, with period 9.5 ≈ empirical PADE 9.2 — **the empirical period-9 CC pair IS Phase 4's L|_{D_W} below-commutant eigenvalue.** Magnitude residual 9% (0.898 vs 0.984) closes structurally at effective level n ≈ 3.91; full magnitude closure via direct ε_k fitting is data-limited at k=13 ([`ROUTEB_PRIME_EPS_FIT_RESULT.md`](notes/ROUTEB_PRIME_EPS_FIT_RESULT.md)).

**ε_14, ε_15, ε_16 computed (2026-05-17):** ε_14 = +3.588e-03, ε_15 = +4.161e-03, ε_16 = +4.685e-03. Hadamard radius |ε_k|^{-1/k}: 1.5655 (k=13) → 1.4951 (k=14) → 1.4412 (k=15) → 1.3982 (k=16), inward by 0.043-0.071/step (geometrically decelerating). Envelope ratios: 1.217 → 1.160 → 1.126 (decelerating ~4 pct/step). Empirical growth rate |ε_k|^{1/k}: 0.639 → 0.669 → 0.694 → 0.715 (approaching predicted asymptote 0.984 from below; linear extrapolation reaches asymptote around k≈27). Sign pattern still + + − − − − − − − + + + + + + +. Period-9 oscillation predicts next sign flip around k=19 (3 more ε_k away).

**Subdominant-rate question (open, structural boundary mapped — R3 does NOT close the 2.9% gap):**

The 2.9% gap between T_lead's exact spectrum {43/45, 0} and empirical Hadamard 0.984 + period-9.2 CC oscillation remains open. R3 dark-subspace classification is structurally clean but does not provide the closure mechanism. The c=7/45 subdominant rate has **no finite-truncation discrete-eigenvalue closure available**. Every natural finite-rank operator over Q has been probed and either:
- Has trivial spectrum (K_k = {1, 0, ..., 0} via 3-fiber row-equality + marginal consistency).
- Has continuous-on-circle spectrum (U_n, Phi_omega = T^ω ∘ U_n: continuous distribution at radius 0.319 or 0.587, no discrete CC pair).
- Is identically zero on trivial-twist class projections at any modulus (U_n → W_n exact, kills all class averages by cube-root cancellation).
- Doesn't recover T_lead's 43/45 (which is itself a class-resolved coherent-summation phenomenon, not a primitive eigenvalue).
- Has been ruled out by Tauberian BLOCKER (20-PDF corpus, Mode H circular) or Nisoli budget infeasibility (M_3''=24.4 blown 18× under realistic Tao C_A). **[2026-07-14: the Tauberian corpus is now FULLY disposed — the last PARTIAL (BGT candidate E) closed NO_FIT on the pre-registered escalation branch; no BGT theorem accommodates the multi-regime asymptotic. See "Latest" above.]**

**Best finite-truncation closure available:** T_lead's eigenvalue **43/45 = 1 − Σ_g W_+(g) = 1 − 2/45** over Q on (P_+, P_-) class-resolved space (`T_LEAD_CORRECTED_DISPOSITION.md`). First positive algebraic spectral result. Within-level rate-carrier, but Nisoli closure inequality `|K|·K^{-A}·M_3'' < 1` fails at realistic A under Tao C_A bookkeeping.

**Surviving productive directions:**
1. **Continue ε_n exact extension** via R77.7 v2 modular CRT (~3-10hr per coefficient). Empirical-discriminative.
2. **V'_M with phase parameters** (T_V Route B, 5-10 sessions). Substantial theoretical reconstruction.
3. **Paper-grade documentation** of the structural boundary combined with Results 1 + 2 above.

See [`SESSION_2026_05_15_STRUCTURAL_BOUNDARY.md`](notes/SESSION_2026_05_15_STRUCTURAL_BOUNDARY.md) for the full session writeup tying together morning's paper-shaped results + evening's structural exhaustion mapping.

**Independent results, paper-grade:**
- **Bilinear bound** `|S_partial(r)| ≤ 2√N` at r ≤ 3 (Path 2), polylog-free `2√p·√N` at r ≥ 4 via Hensel (2026-05-11). See `HENSEL_DISPOSITION.md`.
- **F̂_p Plancherel saturation** `|F̂_p^full(ξ)| = p^{(r+3)/2}` verified across 33 cells (2026-05-11). See `FHAT_THEOREM_VERIFICATION_RESULTS.md`.
- **Prefix-decomposition + Tao bridge** at k = 4..14 with `s_mean(r) ≈ α_det(r) + K_h · log(N/f(N))`, K_h = 3/log(4/3), slope = 1.000 ± 0.005 across 40 verification cells (2026-05-02). See `writeup.md`, `tao_bridge_findings.md`.
- **qx+1 Cramér convergence law** `q^(-θ) = 2^(1-θ) − 1` exact at q ∈ {5, 7, 9, 11} with q=5 match to 0.01% (2026-05-02). See `experiments/16_cramer_root.py`.

**Original status (2026-05-02):** Bridge result to Tao 2022 documented. Three findings consolidated: (a) prefix-decomposition theorem at modular resolutions k = 4..14; (b) `s_mean(r) ≈ α_det(r) + K_h · log(N/f(N))` with slope = 1.000 ± 0.005 at K_h = 3/log(4/3) across two independent observables (σ and first-passage), four modular resolutions, and two data scales; (c) qx+1 Cramér convergence law at q ∈ {5, 7, 9, 11} with q=5 match to 0.01%.

---

## Recent landmarks (post-2026-05-02)

Threads layered on top of the architectural overview below. **For the current state of any of these, read [`STATE.md`](STATE.md) first** — it is the live document and supersedes any drift in this README.

### Leading c=7/45 RIGOROUS UNCONDITIONAL + Syracuse = DWM quantum trajectory NUMERICALLY VERIFIED — 2026-05-15 (morning)

Two paper-shaped results landed in the same session:

- **Leading c = 7/45 RIGOROUS UNCONDITIONAL** ([`THEOREM_C_745.md`](THEOREM_C_745.md)). 8 sections, full hypotheses verbatim, 6-step proof chain via R75 Plancherel × R76 conservation × R77 T_diag × R64.B class-mass × HR74. Independent of all probability-framework questions. The Hasebe-Saigo 2014 monotone-independence overlay was interpretive only — D3 audit confirmed the derivation never depended on HS 2014 Thm 3.4.
- **Syracuse = Davies-Wiseman-Milburn quantum trajectory, numerically verified to 6 sig digits** ([`FRAMEWORK_IDENTIFICATION.md`](notes/FRAMEWORK_IDENTIFICATION.md), [`DWM_MP_G1_RESULT.md`](notes/DWM_MP_G1_RESULT.md)). 6-probe framework-identification arc (H1' → D2 Tier 1 → BMT/bigraph → HP/QSC → AFL → Belavkin/DWM) closed at DWM with adaptive Kraus `M_v^{(j, b_{[1,j-1]})} = 2^{-v/2} A_v^{(j)} σ_{-v}`. DWM-MP-G1+G2 numerical match across all 4 scalar reductions (sum_entries / tr_π / delta_1 / vac_π) for both 3-alt (0.108) and 4-alt (0.609) moments.
### Structural boundary mapped: finite-truncation discrete-eigenvalue paths exhausted — 2026-05-15 (post-compact)

Five new probes targeting the c=7/45 subdominant rate's missing discrete-eigenvalue carrier all confirm R77.6's continuous-spectrum reading + close the structural picture:

- **K_k structural lemma** ([`K_STRUCTURE_RESULT.md`](notes/K_STRUCTURE_RESULT.md)). K_k has spectrum {1, 0, 0, ..., 0} EXACTLY with Jordan chain length k. K_k maps W_{k-1} → 0 exactly. K_k converges to stationary in exactly k Markov steps via rank pattern N_{k-1} → N_{k-2} → ... → 1. R77.4 erratum's "|λ_2| ≈ 10⁻³ growing with k" was numerical noise.
- **U_n → W_n structural lemma + Phi_omega continuous spectrum** ([`INTERLEVEL_U_PROBE_RESULT.md`](notes/INTERLEVEL_U_PROBE_RESULT.md)). Fourier-side Tao transfer U_n: V_n^Fourier → V_{n+1}^Fourier maps V_n entirely into W_n exactly (3rd-root-of-unity phase cancellation). Twisted endomorphism Phi_omega = T^ω ∘ U_n on V_n has top |eigenvalue| converging to 0.319 (ω_3) or 0.587 (ω_3²) but arguments continuously distributed in arcs — no discrete CC pair at θ = 2π/9.2.
- **Bilinear T_M (V_n^M truncation + tensor V_n ⊗ V_n*)** ([`D1_T_M_NEGATIVE_RESULT.md`](notes/D1_T_M_NEGATIVE_RESULT.md)). Two attempted constructions of the bilinear pair-correlation operator both give max |eig| ≈ 0.345, NOT recovering T_lead's 43/45. T_lead's 43/45 is a class-resolved coherent-summation phenomenon at the (P_+, P_-) projection, NOT a primitive eigenvalue of any natural finite-truncation operator.
- **Option III: mod-9 / mod-27 class projection** ([`T_M_class_mod9_spectrum.py`](probes/T_M_class_mod9_spectrum.py)). Trivial-twist projections vanish identically at any modulus 3^k because U_n → W_n exact + each class contains integer-many 3-fibers. Only character-twisted projections give non-trivial structure (which Probe 2 already explored).

Combined with prior-session work (Tauberian 20-PDF BLOCKER 2026-05-13, T_V V_M non-closure 2026-05-12, Nisoli at 43/45 budget-blown 2026-05-12, R77.6 branch-cut + PADE_NUMERICAL z=2 refuted 2026-05-12), the verdict is:

> **No finite-rank operator over Q at finite truncation carries the c=7/45 subdominant rate as a discrete eigenvalue.** T_lead's 43/45 = 1 − Σ_g W_+(g) = 1 − 2/45 is the deepest finite-rank closure available, and even it doesn't close c=7/45 rigorously (Nisoli inequality fails under realistic Tao C_A).

Full session writeup: [`SESSION_2026_05_15_STRUCTURAL_BOUNDARY.md`](notes/SESSION_2026_05_15_STRUCTURAL_BOUNDARY.md).

### Bilinear bound on the R78 wall — DELIVERED 2026-05-11

The dual-side Plancherel bilinear character sum `|S_partial(r)| ≤ C · √N` (Route 2 of the c=7/45 closure landscape) has rigorous proof:

- **r ≤ 3:** strict `|S_partial(r)| ≤ 2√N`, family-level at p ∈ {3, 5, 7, 11}, no tradition ingredients (Path 2 + pushback reconstruction CONFIRMS). See `PATH2_DISPOSITION.md`, `PATH2_BILINEAR_FROM_CLOSED_FORM.md`.
- **r ≥ 4:** polylog-free `|S_partial(r)| ≤ 2√p · √N` via **Hensel-lifted closed form**. The truncated p-adic log gives a saddle equation `1+ps = C_a` linear in s, exactly solvable in Z_p; "Hensel correction collapse" is literal digit extraction `s*(r) = (C_a − 1)/p mod p^{r-1}`. Phase identity `P_a(s*) = M(C_a − 1)` where `M(y) = y − (1+y)·log(1+y) = Σ_{j≥2} (−1)^{j-1}/(j(j-1)) · y^j` falls out algebraically. Triple-verified on closed form (10 cells at 1e-15, max rel dev 6.4e-15 at p=7,r=6) + independent re-derivation CONFIRMS. See `HENSEL_DISPOSITION.md`.

Empirical `K_max/√N ≈ 2.0` at r = 8..20 (`r79b_S_partial_empirical.md`) matches the rigorous constant level. Six literature candidates (Milićević, Banks-Shparlinski, Petrow-Young, Garcia-Young, Pascadi, DFI 1995) + the Polymath8 chain ruled out by object-shape mismatch before Path 2 closed; see `BURGESS_LITERATURE_FINDINGS.md`.

### F̂_p family-level Plancherel saturation — VERIFIED 2026-05-11

The candidate theorem `|F̂_p^full(ξ)| = p^{(r+3)/2}` is verified across 33 cells (primes p ∈ {3..31}, r ∈ {1..6}); mpmath at 50 digits confirms exact algebraic equality to 1e-49 at (p=5, r=3). Standalone result independent of c=7/45 closure. See `FHAT_THEOREM_VERIFICATION_RESULTS.md`, `QX1_FAMILY_THEOREM_ATTEMPT.md`.

### Eighth + ninth probes — **obstruction-separation finding** — 2026-05-12

After the seven-probe trajectory and Tauberian framework arc, two more probes ran at the corrected rate parameters and **separated two obstructions that the seven-probe framing was conflating:**

- **Eighth probe (T_lead at corrected rate)** — `T_LEAD_CORRECTED_DISPOSITION.md`: cross-freq machinery delivers **T_lead = (1/45)·[[7, 9], [28, 36]]** over Q with **spectrum {43/45, 0}**. Eigenvalue **43/45 ≈ 0.9556** on (1, 4). Closed-form origin: 1 − Σ_g W_+(g) = 1 − 2/45. M_3'' = ||(I − T_lead)^{−1}|| = **24.426 exact**. First positive algebraic spectral result of the trajectory. R77.3's rate-1/2 falsification reframes: wrong target, not framework failure.
- **Ninth probe (Nisoli closure at λ=43/45)** — `NISOLI_CLOSURE_CORRECTED_DISPOSITION.md`: closure inequality `|K_bil|·K^{−A}·M_3'' < 1` **FIRES at (r=3, A=3, K=6) with product 0.679 under optimistic C_A = 1.** But under realistic Tao bookkeeping C_A ≥ A^{O(A)}: at A=3, C_3 ≥ 27 blows the K^3 budget by 18×. **No realistic Tao-C_A delivers closure** — H_A_EXTRACTION_HARD.

**The c=7/45 closure landscape goes from "three independent obstructions" (2026-05-11) to "single isolated obstruction" (2026-05-12):**

| Obstruction | Rate-specific? | Status |
|---|---|---|
| (1) Operator-theoretic (discrete eigenvalue at target rate) | rate-specific | **LIFTED at corrected rate** (T_lead 43/45 exact, M_3'' 24.43 exact) |
| (2) Tao Prop 1.17 effective C_A | rate-invariant | **INFEASIBLE — single remaining obstruction** |
| (3) Bilinear \|K\| | independent | DELIVERED 2026-05-11 |

**Single load-bearing unblocker: novel polynomial-in-A Fourier bound on |μ̂_n(ξ)| outside Tao's method.**

### Seven-probe spectral trajectory + Tauberian framework arc — 2026-05-12

Continuation session mapped the structural boundary of c=7/45 closure beyond the bilinear bound. Seven probes (T_3 → R_k → Candidate A → R76 §11 → T_N → cross-frequency closure → T_V iteration → Tauberian scoping) plus the R77.6 Padé re-read converged on a single conclusion:

**Nisoli framework is STRUCTURALLY INAPPLICABLE** — it requires a discrete eigenvalue of a resolvent at rate 1/2, but no Q-constructable finite-rank operator carries one. Rate-1/2 lives in continuous spectrum / branch-cut endpoint structure, NOT a discrete eigenvalue.

First positive structural advance: cross-frequency closure exists on enlarged span V_M = span{M_n^{ab}(g, c)} parameterized by g = v' − v. But V_M doesn't close under iteration n → n+1 (phase + parity obstructions). See `CROSS_FREQ_DISPOSITION.md`, `T_V_DISPOSITION.md`.

The Tauberian framework arc is the live direction: Flajolet-Sedgewick Ch. VI singularity analysis with Chevalier 2507.15394 Thm 1.16 (n^{M − 3/2} via meromorphic h with pole of order M) as cleanest single-theorem candidate. Single-theorem selection is gated on ε_7 exact-rational compute (R77.7 re-fire with new modular-arithmetic solver in flight). See `TAUBERIAN_SCOPING_DISPOSITION.md`, `SESSION_DISPOSITIONS_2026_05_12.md`.

Literature bundle: `C:/Users/Nate/Documents/burgess/literature/` — 73 PDFs across 7 math-field lots, master index at `literature/INDEX.md`.

### c = 7/45 closed-form thread (R75–R79.x)

Plancherel-side derivation of the trajectory measure's structural constant via Tao's Syracuse Markov chain on (Z/3^k)*.

- **R75 / `c_seven_forty_fifth.md`** — derivation of c = 7/45. Rate-1/2 envelope on |ε_n|·2^n appeared stable at k=2..6 (~0.04). **WALKED BACK 2026-05-05** — at k=7 the envelope jumps to 0.150 (4× the supposed plateau); ρ ≈ 0.984 single-pair model also falsified at k=11. Current best fit: order-3 linear recurrence on ε_2..ε_11 with ρ_slow ≈ 0.827 (real). The asymptotic rate question is open; structural form not yet identified.
- **R76 / `result_76_conservation_law.md`** — `Σ_j M_{n+1}(η_0 + j·3^n) = 0` conservation law; reduces rate question to scalar sequence R_n.
- **R77 / `result_77_T_lead_spectrum.md`** through R77.6 / `result_77_6_generating_function.md` — operator-shape attempt: R77.3 falsified the 3-mode geometric ansatz over Q; R77.4 envelope fits gave verdict (M); R77.4 erratum / `result_77_4_K_spectrum_erratum.md` showed K_k itself has no eigenvalue near 1/2 (rate operator is inter-level, not within-level); R77.6 generating-function probe found branch-cut signature at z=2 (type indeterminate at N=5).
- **R77.7 (RE-FIRE IN FLIGHT, `R77_7_V2_*`)** — original k=7 ε-extension killed at ~8.5 hr (`result_77_7_status.md`). New solver design (modular sparse linear solve mod p + CRT recovery + rational reconstruction) targets <2hr; gates the Tauberian framework arc's [3/3] Padé refresh.
- **R78 / `result_78.md`, `result_78_extended.md`, R79 / `result_79.md`, R79b / `r79b_S_partial_empirical.md`** — Path-A obstruction map (Cochrane / van der Corput / direct band-l¹ / band-spectral): all internal subroutes closed for analytical closure of Kalafatelis eq 190 within Tao's framework. **Resolved 2026-05-11 via direct construction (Path 2 + Hensel):** see "Bilinear bound" landmark at top of this section — strict `|S_partial| ≤ 2√N` at r ≤ 3, polylog-free `2√p · √N` at r ≥ 4 via Hensel-lifted closed form.
- **Seven-probe spectral trajectory (2026-05-12)** — `SESSION_DISPOSITIONS_2026_05_12.md`, `READING_A_SCOPING_*`, `CANDIDATE_A_*`, `R76_S11_*`, `T_N_*`, `CROSS_FREQ_*`, `T_V_*`, `TAUBERIAN_SCOPING_*`. After R77.4 erratum's "K_k itself has no eigenvalue near 1/2," the M_3 obstruction was attempted via Reading A scoping → Candidate A construction → R76 §11 verification → T_N construction → cross-frequency closure → T_V spectrum → Tauberian scoping. Conclusion: **Nisoli framework is structurally inapplicable** — no Q-constructable finite-rank operator carries a discrete eigenvalue at rate 1/2; rate-1/2 lives in continuous spectrum / branch-cut endpoint structure. Three positive findings within: (1, 4) eigendirection structurally forced by R64.B's squared class-mass ratio; cross-frequency closure on enlarged span V_M = span{M_n^{ab}(g, c)} (g=0 IS T_diag); R76 §11's empirical P^{+−} = 0 upgraded to rigorous algebraic identity via lift-fiber orthogonality.
- **Tauberian framework arc (2026-05-12, OPEN)** — Flajolet-Sedgewick Ch. VI singularity analysis with Chevalier 2507.15394 Thm 1.16 as cleanest single-theorem candidate. Today's delta diagnostic (`DELTA_DIAGNOSTIC_*`) + Padé extension (`PADE_EXTENSION_*`) both confirm n=2..6 is pre-asymptotic in strong sense: leading "1/30·(1/2)^n" is fast-transient, true asymptotic rate gated on ε_7 from R77.7 v2 (in flight) plus structural reading from PADE_NUMERICAL_* (in flight, uses numerical ε_7..ε_13).

The rate-1/2 rigorous-proof gate is no longer "the single open piece" — bilinear bound delivered the explicit C·√N piece. Remaining open: the analytic structure of E(z) at z=2 (or wherever the true leading singularity is), via Tauberian framework arc.

### Joint 2-3-adic Bohr empirical positive — RETIRED 2026-05-05

Original claim (`result_bohr_probe.md`, 2026-05-04): structured non-CRT-independence between Z/2^a and Z/3^b residues of N=10⁷ Syracuse iterates, z=16.5 at k=20 (a=5, b=4).

**Deflated 2026-05-04/05** by bracket-stratification probe (`result_bohr_probe_strat.md`): per-bracket χ²/df at v ∈ (10⁶, 10⁹] is 0.95 (z = −0.99) and at v > 10⁹ is 0.94 (z = −1.18) — statistically CRT-independent within ±2σ. The original aggregate signal was driven by the v ≤ 100 descent funnel (low-v trajectories transiting toward 1), not by joint structure at the scales relevant to D_emp.

**Status:** retired as a load-bearing structural object. Do not cite as a closure path or as evidence of R58/R60 gap structure.

### qx+1 sweep + sibling 3x±1 study

- **q-sweep / `result_q_sweep_test_2_c_q.md`** — literal hypothesis c_q = S_∞^{(q)}/q falsified for q ≥ 5 (S_k diverges geometrically as (q/3)^k); but renormalized c̃_q := lim S_k^{(q)} / (q/3)^k exists universally.
- **c̃_q structure / `c_tilde_structure_verdict.md` + `c_tilde_q17_probe.py`** — c̃_q = (q − 3)/q confirmed at q=11, 13, 17 within 1%. q=17 (where 2 is NOT a primitive root, like q=7) cleanly fits the formula, ruling out the non-prim-root explanation for q=7's anomalous +0.21 deviation.
- **Sibling 3x±1 forward / `sibling_3x_minus_1_symmetry_verdict.md`** — K_- = σK_+σ proved with σ(r) = −r mod 3^k. Implies S_n^{3x−1} = S_n^{3x+1} as exact rationals; all R76/R77 derived quantities transfer; c=7/45 is automatic for the 3x−1 system by symmetry.
- **Sibling 3x±1 inverse-tree / `duality_S_vs_D_verdict.md` + `duality_followup_verdict.md`** — D_n^((x±1)/3) tables (Agent 2 single-basin from 1, Agent 3 three-basin from {1,5,17} cycles). Raw 10³–10⁴× difference is ~95% sample-size artifact after matched-N control; residual structural difference factor 0.2-4. No clean forward-backward duality D = f(S).

The c̃_q = (q − 3)/q observation is a publishable theorem candidate independent of c = 7/45 closure status.

---

## Index by mathematical field

### Analytic number theory
- Plancherel decomposition on Z/3^k coprime classes (R75; `c_seven_forty_fifth.md`, `result_76_conservation_law.md`)
- **Bilinear bound on R78 wall — DELIVERED 2026-05-11.** Strict `|S_partial(r)| ≤ 2√N` at r ≤ 3 (family-level p ∈ {3, 5, 7, 11}); polylog-free `2√p · √N` at r ≥ 4 via Hensel-lifted closed form (`PATH2_DISPOSITION.md`, `HENSEL_DISPOSITION.md`, `PATH2_BILINEAR_FROM_CLOSED_FORM.md`)
- **F̂_p family-level Plancherel saturation theorem — VERIFIED 2026-05-11.** `|F̂_p^full(ξ)| = p^{(r+3)/2}` across 33 cells, primes 3-31, r ∈ {1..6}; mpmath 50-digit confirms 1e-49 at p=5, r=3 (`FHAT_THEOREM_VERIFICATION_RESULTS.md`, `QX1_FAMILY_THEOREM_ATTEMPT.md`)
- Sign-invariance theorem K_- = σK_+σ (`sibling_3x_minus_1_symmetry_verdict.md`)
- (q−3)/q closed-form candidate for c̃_q at q ∈ {11, 13, 17} (`c_tilde_structure_verdict.md`, `c_tilde_q17_probe.py`)
- q-spectrum probe across q ∈ {3, 5, 7, 11, 13} (`result_qspectrum.md`)
- **Cross-frequency bilinear closure on V_M — POSITIVE STRUCTURAL FINDING 2026-05-12.** Cross-frequency bilinears reduce to enlarged span V_M = span{M_n^{ab}(g, c)} parameterized by g = v' − v; g=0 IS span{P_n^{ab}(c)}; mixed-parity vanishing upgrades R76 §11's empirical P^{+−} = 0 to rigorous algebraic identity via lift-fiber orthogonality. R77 sketch §5's "quadratic forms in {P_n^{ab}(c)}" assertion is FALSE as stated (`CROSS_FREQ_DISPOSITION.md`, `cross_freq_compute.py`)

### Arithmetic dynamics
- Syracuse Markov chain construction at k=1..7 (`K_full.npz`, `K_derived_v2.npz`)
- α_det deterministic prefix algorithm (`experiments/01_alpha_decomposition.py`)
- Inter-level renormalization R̃ operator probe at k=4..7 (`result_R_operator_spectrum.md`, `result_renormalization_spectrum.md`)
- Bridge to Tao 2022 leading term (`tao_bridge_findings.md`)
- **T_diag = (1/5)·[[1, 1], [4, 4]] over Q** — spectrum {0, 1} via char poly λ² − λ = 0; (1, 4) is the λ=1 eigenvector (the **conserved quantity is the squared class-mass direction**, structurally forced by R64.B's class fractions (1/3)² : (2/3)² = 1:4) (`R76_S11_VERIFICATION.md`, `R76_S11_DISPOSITION.md`, `analytical_abc_derivation.md`)
- **W_k filtration via φ_n bilinear pair-form moment — H_CANDIDATE_A_FALSIFIES_F2 2026-05-12.** c_{n, k} := ⟨φ_n, lift_n(R_k)⟩ over Q computed at n=1..6: 15 of 21 are exactly 0/1, only diagonal k=n−1 nonzero; structurally forced by K_n's coset support — φ_n ∈ W_{n−1} by construction. Rate-1/2 does not live in W_k filtration (`CANDIDATE_A_DISPOSITION.md`, `candidate_a_compute.py`)
- **T_V on V_M does not close under iteration — 2026-05-12.** F1 phase obstruction: θ_{v,g} = 2^v · ẽ_g / 3 generically not expressible as ẽ_{G''} − ẽ_G. F2 parity: incoming g ∈ {2, 4} produces only odd outgoing G; V_M = span{g ∈ {0, 2, 4, ...}} doesn't contain odd-G moments (`T_V_DISPOSITION.md`, `T_V_RECURSION.md`, `t_v_compute.py`)
- **Seven-probe spectral trajectory mapping the structural boundary of operator-spectral framework** — T_3 (falsified) → R_k (intractable) → Candidate A (F2) → R76 §11 (inconclusive) → T_N (off_lin underspecified) → cross-frequency (positive on V_M) → T_V (recursion underspecified) → Tauberian scoping (ambiguous, framework right). Combined verdict: no Q-constructable finite-rank operator carries discrete eigenvalue at rate 1/2 (`SESSION_DISPOSITIONS_2026_05_12.md`)

### Probability / Markov chain theory
- Stationary distribution Plancherel mass S_k (`S_k_recursion.csv`, `s_infinity_exact.py`)
- Convergence rate ε_k = S_k − 7/15 trajectory at k=2..11 (`result_epsilon_6.md` … `result_epsilon_11.md`); extended to k=13 (slow oscillating mode at ρ ≈ 0.984, period ≈ 9.2 in k-space)
- Direct K_k top-10 eigenvalue spectrum at k=5,6,7 (`result_eigenvalue_spectrum.md`)
- Order-3 linear recurrence characterization on ε_2..ε_11, ρ_slow ≈ 0.83 (`result_renormalization_recurrence_fits.csv`, `result_renormalization_spectrum.md`); WALKED BACK 2026-05-06 as window-unstable
- Cycle obstruction null result (`result_cycle_obstruction.md`)
- **R77.7 v2 — modular CRT solver for exact-rational π_7 — COMPLETED 2026-05-12.** ε_7 = -1.175236830374320×10⁻³ as exact Fraction (~4485 bits), computed in 39.4 min with 151 primes + 1 witness verification (0/1458 reconstruction failures, sum(π)=1). **~13× speedup vs original R77.7 killed at 8.5hr.** |ε_7|·2^7 = 0.1504 confirms the 4.7× envelope jump at exact precision — R76 §11's "(1/30)·(1/2)^n + O((1/4)^n)" conjecture algebraically refuted (`R77_7_V2_*`, `result_77_7_v2.py`, `experiments_output/result_77_7_eps_exact_through_k7_v2.json`)

### Hierarchical Bayesian statistics (Collatz σ)
- Hierarchical NB GLM with Stan for σ vs log n (`fit.py`, `experiments/nb2_glm.stan`)
- Bonacorsi-Bordoni 2026 NB GLM replication (`experiments/04_head_to_head_nb_glm.py`, `experiments/06_bb_replication.py`, `experiments/06b_bb_pathfinder.py`)

### Ergodic theory / measure theory
- Connection to Lagarias 1985 2-adic equidistribution (`literature_check.md`)
- Tao 2022 bridge at 40 verification cells (`tao_bridge_findings.md`)
- TA.1 N-stability of σ structural offset across N ∈ 2²⁵..2³² (`experiments/36_TA1_sigma_offset_N_sweep.py`)
- Unconditional ensemble mean v_2 = 2.102 vs Geom(½) prediction 2.0 across 2.8M trajectories (`result_density_one_v2_bounds.md`)

### Analytic combinatorics / Tauberian theory (NEW 2026-05-12)
- **R77.6 Padé analysis of E(z) = Σ ε_n z^n** — diagonal [n/n] poles drift monotonically toward z=2 from above the real axis with convergence ratio ~0.67 (slower than exponential → consistent with branch-cut, NOT simple pole) (`result_77_6_generating_function.md`)
- **Tauberian scoping probe** — Flajolet-Sedgewick Ch. VI singularity analysis is the right level of abstraction; Chevalier 2507.15394 Thm 1.16 (n^{M − 3/2} via meromorphic h with pole of order M) cleanest single-theorem candidate; Newman-Zagier excluded (Dirichlet, requires pole). Chevalier Thm 1.14 pure √-singularity FALSIFIED at leading order (n^{-3/2} prediction has growing product across n=2..6) (`TAUBERIAN_SCOPING_DISPOSITION.md`, `tauberian_verify.py`)
- **Delta diagnostic** — H_DELTA_IRREGULAR. All five pre-registered ansatze (geometric, power-law, log, two-term Prony, oscillating·(1/2)^n) fail held-out n=6 with residuals 1.8× to 4.0×; |δ_n/ε_n| ∈ [0.5, 3.0] indicates pre-asymptotic regime (`DELTA_DIAGNOSTIC_DISPOSITION.md`, `delta_diagnostic.py`)
- **Padé extension probe** — H_AMBIGUOUS within n=2..6 window; rules out H_COMPLEX_SECONDARY (no off-axis primary poles) and H_PURE_SIMPLE_POLE (ratios r_n = |ε_n|/|ε_{n-1}| not monotone-approaching 0.5; +0.035, −0.018, −0.030, −0.068 with accelerating downward deviations). Consistent with branch-cut at z=2 with negative subleading coefficient (`PADE_EXTENSION_DISPOSITION.md`)
- **Padé numerical extension — LANDED 2026-05-12: H_TWO_SINGULARITIES_VISIBLE.** Extending Padé budget to numerical ε_7..ε_13 REFUTES z=2 as the leading singularity. R77.6's monotone-to-z=2 pattern BREAKS at [3/3] (first approximant with numerical ε_7): closest pole drifts 0.076 → 0.070 → **0.079** (reverses), then [4/4] jumps inward to z=1.681, [5/5] gives complex-conjugate pair 0.829 ± 0.330i. Hadamard radius from |ε_n|^{1/n} at n=10..13: 2.06 → 1.81 → 1.66 → **1.57** (monotone inward). Sign pattern + + − − − − − − − + + + + has single zero-crossing at n=9→10, period ≈ 9-10 matching STATE.md's slow oscillating mode. Verdict: leading singularity is at |z| ≈ 1.5-1.7 (complex pair plausible), z=2 is sub-leading. R77.6's z=2 reading was a fast-transient fingerprint visible only in n=2..6. Tauberian framework arc redirects to multi-singularity (Flajolet-Sedgewick §VI.4 / §VI.5) (`PADE_NUMERICAL_DISPOSITION.md`, `pade_numerical.py`)
- **Literature bundle** — 18 Tauberian / analytic combinatorics PDFs at `C:/Users/Nate/Documents/burgess/literature/tauberian/`, indexed by Hank's curated reading order (Borwein survey → Haggstrom basics → Flajolet-Sedgewick Ch. VI → Chevalier 2507.15394 + 2504.16233 modern guide)

### qx+1 generalization
- Cramér convergence law q^(−θ) = 2^(1−θ) − 1 exact (`experiments/16_cramer_root.py`, `experiments/17_cramer_dual_verification.py`)
- Per-prime decomposition q ∈ {3, 5, 7, 9, 11} (`experiments/10_q_decomposition.py`, `experiments/12_q_convrate_analytical.py`, `experiments/13_cross_q_unification.py`)
- Cycle classification for q ∈ {5, 7, 11, 13} (`experiments/22_q5_cycle_detection.py`, `experiments/29_qx1_cycle_classification.py`, `experiments/36_q5_fourth_cycle_search.py`)
- Sibling 3x±1 forward-symmetry and inverse-tree asymmetry (`sibling_3x_minus_1_symmetry_verdict.md`, `duality_S_vs_D_verdict.md`, `duality_followup_verdict.md`)

### Numerical methods
- Exact-rational stationary computation k ≤ 5 (Gaussian elimination over Q; `lifting_operator_spectral.py`)
- Float64 power iteration with scipy.eigs cross-check at k=6, 7 (`result_epsilon_6.py`, `result_epsilon_7.py`, `result_epsilon_7_verify.py`)
- Matrix-free power iteration for k=8..11 (sparse + Krylov; `result_epsilon_8.py` … `result_epsilon_11.py`)
- Exact-rational K_n verification at n=1..6 via `candidate_a_compute.py` (Markov chain stationary vector over Q using `fractions.Fraction`)
- **R77.7 v2 modular CRT solver (2026-05-12, in flight)** — sparse linear solve mod p for ~20-30 primes + Chinese Remainder Theorem recovery + rational reconstruction; targets <2hr at k=7 vs original 8.5hr-killed O(N³) Fraction Gauss elimination
- **Hensel digit-extraction verification** — `hensel_approach_a_verify_fast.py` (numpy vectorized) verified the closed-form Hensel saddle s*(r) = (C_a − 1)/p mod p^{r-1} at 10 cells across p ∈ {3, 5, 7, 11} and r ∈ {4, 5, 6}, max rel dev 6.4e-15 at (p=7, r=6)
- **Cross-frequency closure verification** — `cross_freq_compute.py` confirmed at n=2, 3: augmented rank 6 (n=2) and 7 (n=3) vs P-only rank 1, demonstrating V_M has 5-6 new dimensions beyond span{P_n^{ab}(c)}
- **T_V iteration obstruction verification** — `t_v_compute.py` demonstrates both F1 phase obstruction and F2 parity obstruction explicitly with worked examples
- FFT-based Plancherel mass computation (`result_epsilon_10.py` cross-validation chain)

---

## TL;DR — what's here

**Current focus (2026-07-17).** The standalone **qx+1 universal-rate paper** — `‖π_k‖² ~ C_q·3^{−k}`, derived at mechanism — is in its final step: the spectral-gap **L3** (`r_q < 1` for `d = ord_q(2) ≥ 3`, `= 1` at the critical `d = 2 ⟺ q = 3`), under a 6-phase falsifier-first campaign. Phases 0–1 closed; Phase 2b delivered **two proven results** (THEOREM D1 toy gap + THEOREM Real-T1, the kinematic half of the q=3 boundary). **Phase 2c — the limit theorem `ρ(M_tower,L)→1/3`** is now underway (Wilson derives, agents gate): the dynamical partner is pinned as `ρ(M_tower)` exactly (Probe G0), the defect object is frozen with a gate-confirmed sector kernel + an exact `3^{L−1−j}\|k` depth-selection rule (Probe 2c0), and a **Collatz–Wielandt corrector chain** is shrinking the bracket — rung-1 `β*=3/5` (2.6×, kills the bad-cell dichotomy) and rung-2 (the trit `τ`: residual `= 4/9 − κ·W⁻(τ)`, Wilson's blind derivation gate-vindicated on all four predictions). The open crux is the contraction estimate; amplitude-decay is ruled out, leaving selection-rule × cascade-tax `3^{−j}`. First hand-derived gap: the toy `r(λ) = (1−λ²)/(1+λ²)` ([`BRIEF_D1_TOY_GAP.md`](BRIEF_D1_TOY_GAP.md)). Live log: [`STATE.md`](STATE.md).

**Repo layout:** artifacts are foldered by type — `probes/` (scripts), `results/` (`result_*.md`), `notes/` (worksheets/briefs), `logs/`, `outputs/` (numeric data); the repo root holds only this README, `STATE.md`, and the active L3 campaign docs. See [How to navigate this repo](#how-to-navigate-this-repo).

---

The rest of this section describes the project's founding result — the **3x+1 prefix decomposition**:

For odd integers n, the total stopping time σ(n) under the Collatz map has structure indexed by the residue r = n mod 2^k via a *deterministic prefix*. Tracking the symbolic state (a, c) such that orbit value = a·m + c (where m is the integer tail of n) under the Collatz iteration until a becomes odd terminates at a_final ∈ {3^j : 1 ≤ j ≤ k}. The 2^(k−1) odd residue classes therefore collapse onto exactly k distinct conditional distributions of σ.

α_det(r) := prefix_steps(r) + K_h · log(a_final(r) / 2^k), where K_h = 3/log(4/3) ≈ 10.4282, predicts:

1. **σ-intercepts per class** (the original definition; verified at k = 4..12).
2. **First-passage-time means per class** for the orbit reaching ≤ f(N) for *any* threshold f, with offset matching Tao 2022's (5.15) leading term `K_h · log(N/f(N))` to ≤ 1 step (1%-trimmed mean) across 40 verification cells.

The bridge to Tao 2022 is structurally `s_mean(r; f) ≈ α_det(r) + K_h · log(N/f(N)) + ε`, where ε is small, observable-dependent, and stable across modular resolution and N.

---

## How to navigate this repo

### Repository layout (foldered by artifact type; audited & cleaned 2026-07-17)

> **New:** [`RESEARCH_ARC.md`](RESEARCH_ARC.md) is the top-level map of the *whole* program — every result (proven / open / falsified / retracted), the May→July chronology, the field-of-math taxonomy, and the replication path. Start there for orientation.

| Folder | Contents |
|---|---|
| **`papers/`** | The 4 finished papers (Paper 1–4, author's own PDFs) — see [`papers/README.md`](papers/README.md). |
| **`probes/`** | All analysis scripts. **Run from the repo root:** `python probes/<name>.py`. Live L3-campaign probes use relative paths + `probes/paths.py` (repo-root anchor); ~214 legacy probes hardcode paths (superseded). |
| **`results/`** | `result_*.md` write-ups (one per probe / gate). |
| **`notes/`** | The exploratory record — **131 capstones** (`*_DISPOSITION`, `*_verdict`, `*_RESULT`, canonical synthesis) at top; **`notes/_worksheets/`** holds the 321 subordinate scaffolds. See [`notes/README.md`](notes/README.md). |
| **`logs/`, `outputs/`** | Run logs / generated numeric artifacts. |
| **`archive/`** | Dead-end explorations + superseded root docs (incl. the May `PHASE1/2/4_DARK` ancillary) — kept for provenance, out of the live tree. See [`archive/README.md`](archive/README.md). |
| **`references/`, `data/`, `figures/`, `__lean_check/`** | Literature (PDFs gitignored), input data, plots, and the Lean 4 verification project. |
| **Root** | `README.md`, `RESEARCH_ARC.md`, `STATE.md`, and the live campaign docs: `L3_DEFINITIONS.md`, `PHASE1_WORKSHEET.md`, `PHASE3_*.md`, `BRIEF_D1_TOY_GAP.md`, `THEOREM_C_745.md`, `QX1_UNIVERSAL_RATE_WRITEUP_2026_07_14.md`. |

The 125 MB of regeneratable `.npz` stationary caches and all third-party literature PDFs are **untracked** (see `.gitignore`); the tracked repo is ~62 MB. Author-draft folders (`_Author_Papers/`, `JNT Submission/`, …) are kept on disk but untracked.

### If you want the headline result
1. Read [`writeup.md`](notes/writeup.md) — Result 1, Result 3, and the subsection "α_det predicts mean first-passage time and matches Tao (5.15) at the per-class level" within Result 3.
2. Then [`tao_bridge_findings.md`](notes/tao_bridge_findings.md) for the TA.1/TA.2/TA.3 tightening of ε(N).

### If you want to verify a specific claim
- Each substantive claim in `writeup.md` traces to a `findings.md` entry and one or more numbered experiments. Use the experiment index below to locate the script. Run it; CSV outputs land in `experiments_output/`.

### If you want to check the prior-art positioning
- [`literature_check.md`](notes/literature_check.md) — audit identifying Terras 1976 Lemma 4 as the asymptotic predecessor of the prefix decomposition, plus connections to Sinai 2003, Tao 2022, Bonacorsi & Bordoni 2026.

### If you want to see the audit trail with sanity-check protocol applied per finding
- [`findings.md`](notes/findings.md) — chronological log, append-only.
- [`agent2_findings.md`](notes/agent2_findings.md) — trajectory-measure deep dive (v=4/v=10 spike mechanism, q=5 trajectory v, MGF preservation across q).
- [`compute_threads_findings.md`](notes/compute_threads_findings.md) — σ-record extension to OEIS A006877 b-file, prefix-tail mechanism analysis, q=5 cycle search.

### If you want the qx+1 generalization (companion to the 3x+1 work)
- The qx+1 Cramér convergence-rate result lives separately. Code is mixed in among the experiments below (10–22 range and 28/29 range). For prior consolidations see `findings.md` entries dated 2026-05-01 onward and the auto-memory file `project_collatz_qx1.md` (in the user's external memory).

---

## Documents

**7/15 constant campaign — Thread 3 (2026-07-21):**

| File | Contents |
|---|---|
| [`BOUNDARY_THEOREM_MAP.md`](BOUNDARY_THEOREM_MAP.md) | Theorem→gate bridge: Thm 1–7 (boundary theorem) + the Thread-3 constant arc (Theorem F, off-diagonal ledger, channel engine, strata, collision identity, spectral form) mapped to verifying probes/commits. |
| [`results/result_offdiag_R6.md`](results/result_offdiag_R6.md) | **R6** — the interference ledger; `OffDiag₂ = −4/21` derived from v≠v′ Ramanujan sums; `Σ_{k≥2} OffDiag = −1/5`. |
| [`results/result_engine_R7.md`](results/result_engine_R7.md) | **R7** — the ⟨4⟩-orbit channel engine; `OffDiag_k = (2/3)Σ_m 4^{−m}C_k(m)` == frozen, k=2..5; sign-flip derived. |
| [`results/result_strata_R8.md`](results/result_strata_R8.md) | **R8** — uniform kill (`C_k≡0`), m-independent band count, ledger↔deviation weld, k-independent stratum weights `W_j`. |
| [`results/result_gamma_R9.md`](results/result_gamma_R9.md) | **R9** — the collision identity `S_K = 2Σ_m 4^{−m}γ_{K−1}(τ_m)` exact K=2..6; `γ_n(0)=X_n` corpus weld; one stationary `γ_∞` on ℤ₃. |
| [`results/result_charledger_R10.md`](results/result_charledger_R10.md) | **R10** — the spectral form; `OffDiag_{r+1}=2Λ_r` (frozen character layers) exact r=1..5; `S_∞=7/15 ⟺ Σ_{r≥1}Λ_r=−1/10`. |
| [`results/result_chirpkernel_R11.md`](results/result_chirpkernel_R11.md) | **R11** — the chirp kernel; `2Λ_r=ε_{r+1}−ε_r` welds the value & ε-rate campaigns; U shell block-diagonal; `Λ_r=⟨μ̂,K_rμ̂⟩`; exact ledger extended to Λ₆. |
| [`results/result_lobes_R12.md`](results/result_lobes_R12.md) | **R12** — support law `U(k,ξ)=0 unless k≡ξ mod 3`; closed loop (angular moments = C-table entries); ε-table certified = d_k (exact k≤8); Λ₇ exact; `\|Re w\|`-moment stabilizes. |
| [`results/result_psi_R13.md`](results/result_psi_R13.md) | **R13** — closed-form weight `Re w=15/(2D)−½`; chirp β *is* the renewal (`t′=β(2⁻ᵛ4ᵗ)`, gate); ψ not uniform (depletion near x≈0). |
| [`results/result_deviation_field_R14.md`](results/result_deviation_field_R14.md) | **R14** — R13-C resolved: no non-uniform ψ (γ_n = partial sum, bounded); object = deviation field δ_r; retarget `Σ[S_r⟨δ_r,Re w⟩+Λ_r^unif]=−1/10`; 2 flags. |
| [`results/result_endpoint_R15.md`](results/result_endpoint_R15.md) | **R15** — endpoint identity `Λ_r=S_r·b_r+Λ_r^unif` exact (anchors byte-exact); no self-conjugate primitive angle (odd lattice); the −1/10 baseline localizes (not reduces) the difficulty; bulk sum ≈−9.85e−4. |
| [`results/result_transport_R16.md`](results/result_transport_R16.md) | **R16** — transport step exact (b_r = PASS-gated operator output, δ₁=0 no initial condition). *(Its "same as R5" crux retired #42, and its q-oscillation reading was a design defect — off the critical surface — both corrected in R23.)* |
| [`results/result_slowmode_R17.md`](results/result_slowmode_R17.md) | **R17** — slow-mode symbol `\|D\|²=1/(5−4cos(πξ/3^r))` (linearized; exact transport non-diagonal, flagged); δ quasi-stationary (⟨\|D\|²⟩≈1/3, contracts at mean rate, width holds); self-map + invariance PASS. |
| [`results/result_regimes_R18.md`](results/result_regimes_R18.md) | **R18** — exact ratio ≈1/2 not √(1/3) (ρ≈0.988=signed-envelope decay, conflation settled); δ **equipartitioned across orders, not high-frequency** (adjective corrected); **branch `T=U₊D₊+U₋D₋` PASS** (#39 not incurred); **max/typical grows 1.34→4.25 ⟹ triangle/max bound dead** (Prop-1.17 gap concrete). |
| [`results/result_exceptional_R19.md`](results/result_exceptional_R19.md) | **R19** (the decider) — **fixed low-harmonic `A_r(m)` are typical/depleted, not exceptional ⟹ equipartition route survives R18-D**; within-stratum spikes real (40× @ r=7). *(R19-A ⟨2⟩ and R19-D A-side spike framings corrected in R20; the decider stands.)* |
| [`results/result_window_R20.md`](results/result_window_R20.md) | **R20** — **`Λ_r=Σ_m 4^{−m}A_r(m)=OffDiag/2` exact ⟹ deviation-field = R7 engine = collision-γ, one computation**; running sum saturates by m≈3–4 ⟹ owed = summable bound on **thin O(1) window m≲r**; m=9/27 oscillate (R13-C read an oscillation) + weight-suppressed; additive slow mode → trivial char at (2/3)^r; `A_r(3^{r−1})=−S_r/2` vacated. |
| [`results/result_ratio_R21.md`](results/result_ratio_R21.md) | **R21** — **`γ_r(τ_m)=3^r·ρ_r(4^{−m})`, theorem = `Σ4^{−m}f(4^{−m})=7/30`** (Haar ratio-density); independent **group-division gate PASS** (#41 not incurred); f singular at identity (`=X_r`) but **weight sits on the smooth part ⟹ singularity weight-suppressed**; R21-C pre-registered miss (argmax₇=2⁸, near-degenerate peak / Prop 1.17). |
| [`results/result_strata2_R22.md`](results/result_strata2_R22.md) | **R22** — **stratum-only GATE FAILS (honest negative)**: 3 j=0 orbit points → 3 distinct limits (0.717/0.476/0.868), within-stratum std/mean ≈0.24 stable ⟹ **§3 stratum-reduction void**; diagnostics survive (F(0)=2/3, diffs→7/15, Haar sum=1 exact; weight 95% j=0, 75%/19% on m=1/m=2) ⟹ **confirms R20 thin-window, retires stratum-mean route**. |
| [`results/result_critical_family_R23.md`](results/result_critical_family_R23.md) | **R23** — **triple negative**: Conj-2 `S_∞(q)=(3q²+1)/(2q(q²+1))` **falsified** (q=5,7,11,13 miss; 7/15=M₄/M₃ is q=3 coincidence; builder validated vs byte-gate + Wilson's S₁); Conj-1 `f(τ₂)=10/21` **falsified** (γ(τ₂) fell through it at r=10); f(τ₁),f(τ₂) **not locked**; functional-eq **proved non-existent**; **#42 retires "same as R5"** ⟹ owed = genuine convergence estimate, no shortcut. |

**qx+1 paper — Result 1 / L3 spectral-gap campaign (2026-07-16):**

| File | Contents |
|---|---|
| [`STATE.md`](STATE.md) | **Live research log** — supersedes any drift here. Per-result entries (R1–R46) + the L3 campaign gates (G0/G0b/G0c/G0c′/G1) with pre-registrations, verdicts, and walk-backs. |
| [`L3_DEFINITIONS.md`](L3_DEFINITIONS.md) | The frozen L3 object (Phase 0): one `r_q`, three welded coordinates (build_M / renewal `A(z)` / `c_k`); dictionary + mass identity; the L3 statement and boundary clause. |
| [`PHASE1_WORKSHEET.md`](PHASE1_WORKSHEET.md) | The five substrate lemmas (Phase 1): FORGET / ONE-STEP / INTERTWINE / REFINE / PYTHAGORAS, proved and machine-verified (incl. q=1093). |
| [`result_phase2a_recon.md`](results/result_phase2a_recon.md) | Phase 2a boundary recon (Q1–Q6): names the theorem = **L→∞ coalescence**; both 2's (phase ⟨2⟩ + weight ½) load-bearing = marginality of the actual 3x+1 map. |
| [`result_phase2b_s1.md`](results/result_phase2b_s1.md) | Phase 2b Session 1 (instruments A–E, C): the toy `M(q,−1,λ)` frozen + hand-solvable; diagonal-ray localization; L=4 walls out. |
| [`BRIEF_D1_TOY_GAP.md`](BRIEF_D1_TOY_GAP.md) | **THEOREM D1** (complete) — the first hand-derived spectral gap of the program: `r(λ) = (1−λ²)/(1+λ²)`, derived then met the pre-published sweep five-for-five; maximality proven via nilpotence. |
| [`result_phase2b_s2.md`](results/result_phase2b_s2.md) | Phase 2b Session 2: D1 committed + the **D3 lead** (real q=3 invariant ray `Σw_r² ≈ 1/3 + (2/3)·2^{−D}`, six-digit match; eigenvalue braid toward the EP); requests F/G. |
| [`result_phase2b_F.md`](results/result_phase2b_F.md) | Request F: the e=−1 sub-block is **nilpotent** (ρ₋=0) — D1's maximality closes with maximal margin. |
| [`result_phase2b_Dmax.md`](results/result_phase2b_Dmax.md) | **LEMMA D1-MAX** (acyclicity/nilpotence proof) + gates ⇒ **THEOREM D1 COMPLETE**; q=3 corollary exact (crossover at λ=1/√2). |
| [`result_phase2b_H.md`](results/result_phase2b_H.md) | Real q=3 Δ-channel operator: **H_EXACT + H_CIRC** confirmed (circulant family complete to L=3, 18/18); the toy is the **parity quotient**. |
| [`result_phase2b_LALB.md`](results/result_phase2b_LALB.md) | Skeleton lemmas verified: **L-A** (no-return C↛Δ) and **L-B** (k=0 co-invariant eigenvector gauge-factorizes exactly). |
| [`result_phase2b_R.md`](results/result_phase2b_R.md) | Probe R: R1 hit the **STOP** (modular (e,γ) closure not clean — gauge broken); R2/L=4 deferred (local wall). Superseded by Real-T1. |
| [`result_phase2b_J.md`](results/result_phase2b_J.md) | Swap-involution J refuted (walk-back #14); the extracted constraint (invariance must act on the carry as an **integer map**); board consolidation. |
| [`result_phase2b_T1.md`](results/result_phase2b_T1.md) | **★ THEOREM Real-T1** (proven; program's 2nd result, 1st on the real operator): the exact eigenvalues are the **twisted autocorrelations of the halving weights** `c_k = Σ_δ w_δ² ω^δ`; closed-form eigenvectors gate-verified 18/18 at L=3. |
| [`QX1_UNIVERSAL_RATE_WRITEUP_2026_07_14.md`](QX1_UNIVERSAL_RATE_WRITEUP_2026_07_14.md) | Standalone qx+1 universal-rate paper draft: `S_k^(q) ~ (q/3)^k`, derived at mechanism; three pillars (rate, constant, correction). |
| [`THEOREM_C_745.md`](THEOREM_C_745.md) | The `c = 7/45` rigorous unconditional result (`S_k = 3^k‖d_k‖² → 7/15`) — the q=3 marginal fixed point. |

**Original 3x+1 prefix-decomposition / Tao-bridge work:**

| File | Contents |
|---|---|
| [`writeup.md`](notes/writeup.md) | Canonical result document. Result 1 (slope universality + non-monotone β oscillation), Result 2 (tail shape), Result 3 (prefix decomposition + Tao bridge subsection), Related Work (B&B comparison + Pathfinder caveat), Limitations. |
| [`findings.md`](notes/findings.md) | Append-only chronological audit trail. Every empirical finding gets sanity-check protocol entries (sampling bias / definition / finite-N / parity / numerical precision). ~600 lines. |
| [`agent2_findings.md`](notes/agent2_findings.md) | Trajectory-measure characterization: q=3 trajectory v moments, MGF preservation across q, m mod 32/2048/131072 pushforward (mechanism for v=4/v=10 spikes), q=5 unconditional trajectory v. |
| [`compute_threads_findings.md`](notes/compute_threads_findings.md) | σ-record class-fraction analysis (T1.1/T1.2/T1.5/T1.6/TB.2). Gaussian-tail Gumbel mechanism for prefix-class σ-record fractions, replacing earlier exponential-θ guess. |
| [`tao_bridge_findings.md`](notes/tao_bridge_findings.md) | TA.1 N-stability of σ structural offset (constant ≈ −2.45 across N = 2²⁵..2³²), TA.2 trim-quantile sweep (q* = 1.18% drives gap to 0 at √N), TA.3 parametric fit (gap ≈ −2.35 + 0.486·log(threshold)). |
| [`closed_form_findings.md`](notes/closed_form_findings.md) | Closed-form derivations for the bridge structural constants. ⟨α_det⟩ = log(6)/log(4/3) DERIVED exactly. ε(σ) and slope-on-log(threshold) ruled out as having clean closed forms; trace back to either Lagarias trajectory measure (open) or finite-N μ_β characterization (TA.1 follow-up). |
| [`literature_check.md`](notes/literature_check.md) | Prior-art audit. Terras 1976 Lemma 4, Sinai 2003, Lagarias 1985, Tao 2022, B&B 2026. |
| [`one_sheet_lin.py`](probes/one_sheet_lin.py) / [`one_sheet_yosef.py`](probes/one_sheet_yosef.py) | PDF generators for one-sheet summaries (Lin: 3x+1; Yosef: qx+1). |

---

## Experiment index (organized by theme)

Numbered by creation order; some numbers collide because the project ran two parallel agents. Filename disambiguates.

### Stage 1–4 pipeline (original Bayesian fit)
| Script | Purpose |
|---|---|
| [`generate.py`](probes/generate.py) | Numba memoized σ / syracuse / max_excursion / residues for n ∈ [1, N]. Outputs `data/main_N{N}.parquet`. |
| [`generate_q.py`](probes/generate_q.py) | qx+1 generalization: writes `data/q_main_q{q}_N{N}.parquet`. |
| [`analyze.py`](probes/analyze.py) | Stage 2 EDA: σ vs log n by mod-16 class, v-distribution, residual tails. |
| [`stage3_prep.py`](probes/stage3_prep.py) | Stage 3 input: odd-only filter, class index, uniform stratified subsample. |
| [`fit.py`](probes/fit.py) | Stage 3 hierarchical Stan fit (k=6, k=10). Outputs to `fits/{tag}/`. |
| [`diagnose.py`](probes/diagnose.py) | Stage 4 posterior summary, GPD on tails, posterior tail probabilities. Outputs to `stage4_results/{tag}/`. |

### Per-class structure (3x+1)
| # | Script | Purpose |
|---|---|---|
| 01 | [`experiments/01_alpha_decomposition.py`](experiments/01_alpha_decomposition.py) | Per-class OLS α(r) vs predicted α_det(r) at given k. |
| 02 | [`experiments/02_moment_universality.py`](experiments/02_moment_universality.py) | Higher per-class moments (variance, skew, kurtosis) vs prefix prediction. |
| 03 | [`experiments/03_n_scaling.py`](experiments/03_n_scaling.py) | μ_β scaling N ∈ {2²⁰..2²⁵}. |
| 05 | [`experiments/05_cfinal_ks_analysis.py`](experiments/05_cfinal_ks_analysis.py) | Within-a_final c_final substructure via KS tests. |
| 07 | [`experiments/07_anderson_darling.py`](experiments/07_anderson_darling.py) | Distributional clustering of per-class σ residuals via Anderson-Darling. |
| 08 | [`experiments/08_all_n_decomposition.py`](experiments/08_all_n_decomposition.py) | Decomposition extended to all n (odd ∪ even). |
| 09 | [`experiments/09_multi_stat_decomposition.py`](experiments/09_multi_stat_decomposition.py) | Decomposition for σ, syracuse, odd_steps, even_steps, log(max_excursion). |
| 24 | [`experiments/24_k_sweep_alpha_decomposition.py`](experiments/24_k_sweep_alpha_decomposition.py) | k-sweep at N=2²⁷ for k ∈ {4..12}; noise-floor ratio band. |

### B&B NB GLM replication
| # | Script | Purpose |
|---|---|---|
| 04 | [`experiments/04_head_to_head_nb_glm.py`](experiments/04_head_to_head_nb_glm.py) | Frequentist NB GLM head-to-head (M0..M4). |
| 06 | [`experiments/06_bb_replication.py`](experiments/06_bb_replication.py) | Bayesian NB GLM via cmdstanpy NUTS. |
| 06b | [`experiments/06b_bb_pathfinder.py`](experiments/06b_bb_pathfinder.py) | Pathfinder VI fallback (used when NUTS multi-chain locked at N_train=500K). |
| — | [`experiments/nb2_glm.stan`](experiments/nb2_glm.stan) | Shared Stan model. |

### Trajectory measure
| # | Script | Purpose |
|---|---|---|
| 15 / 15b | [`experiments/15_step_variance.py`](experiments/15_step_variance.py), [`experiments/15b_step_variance_unconditional.py`](experiments/15b_step_variance_unconditional.py) | Conditional vs unconditional v variance. |
| 25 | [`experiments/25_trajectory_measure.py`](experiments/25_trajectory_measure.py) | High-resolution v-distribution at N_start=10⁸, T=200; v=4 and v=10 spike characterization. |
| 27 | [`experiments/27_m_residue_pushforward.py`](experiments/27_m_residue_pushforward.py) | m mod 32 / 2048 / 131072 pushforward (mechanism for v=4/v=10/v=16 spikes). Agent 2. |
| 28 | [`experiments/28_per_octave_trajectory_E_v.py`](experiments/28_per_octave_trajectory_E_v.py) | Per-octave trajectory E[v] for the K(E[v]) closed-form prediction of β_local. |
| 28 | [`experiments/28_q5_trajectory_measure.py`](experiments/28_q5_trajectory_measure.py) | q=5 trajectory v-distribution on convergent orbits. Agent 2. |
| 29 | [`experiments/29_v_step_correlation.py`](experiments/29_v_step_correlation.py) | Lag-1 autocorrelation of v along Syracuse trajectories. |

### β oscillation and N-extension
| # | Script | Purpose |
|---|---|---|
| 26 | [`experiments/26_mu_beta_n_extension.py`](experiments/26_mu_beta_n_extension.py) | Streaming OLS at N up to 2³². Non-monotone β oscillation finding. |
| 27 | [`experiments/27_beta_oscillation_diagnostic.py`](experiments/27_beta_oscillation_diagnostic.py) | Per-octave β_local + top-K outlier exclusion (record-σ hypothesis test). |

### First-passage / Tao bridge
| # | Script | Purpose |
|---|---|---|
| 23 | [`experiments/23_sigma_fiber_cardinality.py`](experiments/23_sigma_fiber_cardinality.py) | σ-fiber cardinality (Avenue A diagnostic; cryptographic hardness ruled out). |
| 30 | [`experiments/30_first_passage_alpha_det.py`](experiments/30_first_passage_alpha_det.py) | First Spearman ρ = 1.0 finding for s_median vs α_det at k=8. |
| 31 | [`experiments/31_first_passage_replication.py`](experiments/31_first_passage_replication.py) | Replication at k=8/10/12 × 4 thresholds. |
| 32 | [`experiments/32_alpha_det_K_calibration.py`](experiments/32_alpha_det_K_calibration.py) | K-recalibration test for s_median (slope < 1 mechanism diagnosis). |
| 33 | [`experiments/33_alpha_det_K_calibration_mean.py`](experiments/33_alpha_det_K_calibration_mean.py) | s_mean version: slope = 1 at K_h with raw mean; trim-1% match to Tao. |
| 34 | [`experiments/34_alpha_det_K_calibration_mean_k_sweep.py`](experiments/34_alpha_det_K_calibration_mean_k_sweep.py) | k=8/10/12 mean replication. |
| 35 | [`experiments/35_alpha_det_full_bridge.py`](experiments/35_alpha_det_full_bridge.py) | Full bridge: 4 k × 5 observables (σ + 4 first-passage thresholds) × 2 N. 40-cell verification. |
| 36 | [`experiments/36_TA1_sigma_offset_N_sweep.py`](experiments/36_TA1_sigma_offset_N_sweep.py) | TA.1: σ offset N-stability at N = 2²⁵..2³². |
| 37 | [`experiments/37_TA2_trim_quantile_sweep.py`](experiments/37_TA2_trim_quantile_sweep.py) | TA.2: trim quantile sweep finding q* = 1.18%. |
| 37 | [`experiments/37_alpha_det_k16_verification.py`](experiments/37_alpha_det_k16_verification.py) | Higher-k verification at k=16. Agent 2. |
| 38 | [`experiments/38_TA3_parametric_fit.py`](experiments/38_TA3_parametric_fit.py) | TA.3: parametric fit `gap ≈ −2.35 + 0.486 · log(threshold)`. |
| 39 | [`experiments/39_overshoot_at_first_passage.py`](experiments/39_overshoot_at_first_passage.py) | First-passage overshoot ⟨log(f/v*)⟩ — empirical ≈ 0.298 nats, constant in f. Rules out the K_h·⟨overshoot⟩ explanation for the slope on log(f). |
| 40 | [`experiments/40_K_eff_decomposition.py`](experiments/40_K_eff_decomposition.py) | Parity decomposition of K_eff into odd-σ slope + halving compensation; user's hypothesis K_eff = 9.31 + 0.63 ruled out empirically. |

### qx+1 generalization (companion)
| # | Script | Purpose |
|---|---|---|
| 10 | [`experiments/10_q_decomposition.py`](experiments/10_q_decomposition.py) | qx+1 prefix decomposition at k=6 for q ∈ {3, 5, 7, 11}. |
| 10b | [`experiments/10b_q_partial_correlation.py`](experiments/10b_q_partial_correlation.py) | Partial correlation diagnostics (j-slope vs log(m)). |
| 12 | [`experiments/12_q_convrate_analytical.py`](experiments/12_q_convrate_analytical.py) | log(conv_rate) vs j slope fits. |
| 13 | [`experiments/13_cross_q_unification.py`](experiments/13_cross_q_unification.py) | C ≈ 5/2 universal-multiplier hypothesis test (rejected at q=11). |
| 14 | [`experiments/14_conv_rate_vs_N.py`](experiments/14_conv_rate_vs_N.py) | conv_rate(N) decay exponent. |
| 16 | [`experiments/16_cramer_root.py`](experiments/16_cramer_root.py) | Exact Cramér root: `q^(−θ) = 2^(1−θ) − 1`. |
| 17 | [`experiments/17_cramer_dual_verification.py`](experiments/17_cramer_dual_verification.py) | Dual j-slope and N-decay verification of θ(q). |
| 18 | [`experiments/18_q7_x_binning_diagnostic.py`](experiments/18_q7_x_binning_diagnostic.py) | Pooled-X vs per-class diagnostic at q=7. |
| 19 | [`experiments/19_bahadur_rao.py`](experiments/19_bahadur_rao.py) | Bahadur-Rao 1/√L sub-exponential prefactor test (rejected). |
| 20 | [`experiments/20_m_selection_test.py`](experiments/20_m_selection_test.py) | m-selection partial correlation diagnostic. |
| 21 | [`experiments/21_two_term_fit.py`](experiments/21_two_term_fit.py) | Two-term fit `f(X) = A·X^(−θ) + B` (B → 0; rejected). |
| 22 | [`experiments/22_q5_cycle_detection.py`](experiments/22_q5_cycle_detection.py) | Floyd cycle detection at q=5; non-trivial cycle landings 0.12%. |
| 29 | [`experiments/29_qx1_cycle_classification.py`](experiments/29_qx1_cycle_classification.py) | qx+1 cycle classification at q ∈ {5, 7, 11, 13}. Agent 2. |
| 36 | [`experiments/36_q5_fourth_cycle_search.py`](experiments/36_q5_fourth_cycle_search.py) | q=5 fourth-cycle search. Agent 2. |

### Auxiliary
| # | Script | Purpose |
|---|---|---|
| 30 | [`experiments/30_sigma_records_prefix_analysis.py`](experiments/30_sigma_records_prefix_analysis.py) | σ-record class-fraction Gaussian-Gumbel analysis. Agent 2. |

---

## Data files

`data/main_N{N}.parquet`:

| N (= 2^k or scientific) | rows | size |
|---|---|---|
| 2²⁰ = 1,048,576 | 1.0M | 7.5 MB |
| 4,194,304 | 4.2M | 27 MB |
| 8,388,608 | 8.4M | 52 MB |
| 10,000,000 | 10M | 62 MB |
| 16,777,216 | 16.8M | 102 MB |
| 33,554,432 = 2²⁵ | 33.6M | 201 MB |
| 134,217,728 = 2²⁷ | 134M | 786 MB |

Schema: `n, sigma, syracuse, odd_steps, even_steps, max_excursion, is_record, res_mod_16, res_mod_64, res_mod_256`.

`data/q_main_q{q}_N{N}.parquet`: qx+1 versions for q ∈ {3, 5, 7, 9, 11, 13} at varying N up to 10⁹.

Larger σ caches at N ∈ {2²⁸, 2³⁰, 2³²} are built in-memory by experiment 36 and not persisted (~1, 4, 17 GB int32 respectively).

---

## Reproduction smoke check

If these three pass, the work is intact:

1. **Data generation:** `python generate.py --N 1048576` should finish in ~1.5 s. σ at n=27 should be 111.
2. **Prefix algorithm:** Run by hand on residue r = 21 starting from state (a=64, c=21). Expected: 7 steps, terminating at (a_final=3, c_final=1).
3. **Bridge result:** `python experiments/35_alpha_det_full_bridge.py` should report slope at K_h ≈ 1.000 ± 0.005 across all 40 cells, with offset gaps matching the table in `tao_bridge_findings.md`.

For the headline N-extension finding: `python experiments/36_TA1_sigma_offset_N_sweep.py` reports gap ≈ −2.45 across N = 2²⁵ → 2³², stable to 0.01.

Compute requirements:
- All experiments at N ≤ 2²⁷ run in seconds-to-minutes on a 16-thread CPU.
- N=2³⁰ sigma cache: ~14 s, 4 GB RAM.
- N=2³² sigma cache: ~55 s, 17 GB RAM.

---

## Open follow-ups

Not load-bearing for the bridge claim, but adjacent and worth pursuing:

- **Closed-form derivation of the −2.35 structural constant** (TA.3 intercept). The N-stability data (TA.1) supports it being a structural invariant of the σ distribution; analytical form unresolved.
- **Closed-form derivation of the +0.486 ≈ 1/2 slope on log(threshold)** (TA.3). Hints at a √f scaling but the underlying mechanism is not isolated.
- **Trim-quantile interpretation at larger N.** TA.2's q* = 1.18% does not match a clean log^(−c) N exponent; might be an N=2²⁷ artifact. Re-test at N=2³².
- **q=5 / q=7 first-passage analog.** Does the bridge `s_mean ≈ α_det^(q) + K_q · Δlog` hold for qx+1 with K_q derived from θ(q)? Generalizes the bridge across the qx+1 family.
- **Bonacorsi-side HMC validation** at full N=10⁷ on the NB GLM replication (deferred to N. Bonacorsi at Columbia).

---

## Outreach packages on Desktop

- `collatz_for_lin_2026-05-01/` and `.zip` — Lin (Maryland) review package: writeup, findings, one-sheet PDF, data sample.
- `collatz_for_bonacorsi_2026-05-01_v4_pathfinder/` — B&B framework replication artifacts.
- `collatz_qx1_for_yosef_2026-05-01/` — qx+1 Cramér derivation package for the probability-theory audience.

These predate the Tao bridge result. Re-zipping with the updated `writeup.md` and `tao_bridge_findings.md` would refresh them.
