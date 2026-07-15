# The universal rate of qx+1 stationary L² mass — standalone writeup

**Date:** 2026-07-14, **substantially revised 2026-07-15 (R6–R10).** **Status:** paper-shaped; three results, **now all flowing from one mechanism**. One genuinely open input (Result 1's overlap bound). Independent of every Collatz-closure thread (ε_k / c=7/45 / THEOREM_C_745 / R81b all untouched). Collatz need not be in the title.

> **2026-07-15 revision summary.** Three things changed and one died:
> - **R7 — the object was misidentified.** The primitive is `M_k = q^k‖d_k‖²`, **not** `q^k‖π_k‖²`. See Object below. The q-sweep and the q=3 thread (R74/R75/R76) have been computing **one object by two routes** since 2026-05-03.
> - **R8 — the mechanism is named.** `π_k` **is** a q-adic self-similar measure; the "3" is `1/Σ_v p_v²`, its correlation dimension. Domination = an **overlap bound**, not a character estimate.
> - **R9/R10 — Result 3 was WRONG and is rewritten.** `δ_q ≈ 0.82/ord_q(2)` is refuted by **2.55×10¹³**. The law is **`δ_q = 2^{1−ord_q(2)}·(q−3)/q`**.
> - **R6 — the previously-advertised "one line to a theorem" (generalize R76's conservation) is DEAD.** Conservation ports for free *and is insufficient*; q=3 is structurally unique. See Result 1.
> - **R11 — the collision count.** Result 3's leading term is now an **exact identity**; the prefactor 2 is derived (and the earlier ordered-pair explanation was wrong).
> - **R12 — Tao 2019 read.** His Prop 1.17 does **not** apply to our object; but the contrast names our criticality. See the abstract line below.

---

## ★ The one-line statement (abstract material)

> **`H₂(Geom(1/2)) = −log Σ_{v≥1} 4^{−v} = log 3 = log q` — exactly, at q = 3.**
> The stationary measure's **correlation dimension** is `D₂ = H₂/log q = log 3 / log q`, so
> ### **`D₂ = 1 ⟺ q = 3` — the Collatz multiplier sits exactly at the critical point of its own family.**

The `3` is `1/Σ_v p_v² = 1/E[2^{−v}] = 2²−1`: a property of the **halving law alone**, carrying no information about the multiplier. The Collatz multiplier merely *happens* to equal it. This is the entire content of "3 is special," in standard vocabulary — and it makes the Bernoulli-convolution template (Erdős/Solomyak/Hochman: a one-parameter family with a phase transition located by L² methods) **the same kind of statement**, not an analogy.

**Positioning against Tao 2019** (arXiv:1909.03562): same object, adjacent question, **different Rényi order.** Tao's Remark 1.15 counts the same addresses by **Shannon** entropy — `H₁ = log 4` ⇒ `4^{n+o(n)}` tuples into `3ⁿ` slots, a *comfortable surplus*, which is why fine-scale mixing holds for him. Our L² mass is governed by **collision** entropy — `H₂ = log 3` ⇒ `3ⁿ` effective addresses into `3ⁿ` slots, **exactly saturated**. ⚠️ **Tao's Prop 1.17 (`|μ̂_n(ξ)| ≪_A n^{−A}`, pointwise, uniform) does NOT give our result** — summed over the ~3ⁿ admissible frequencies it yields `≪ 3ⁿ·n^{−2A}`, which *grows*. Cite him for the architecture; state the non-applicability explicitly, or a reader will assume otherwise.

**Five independent derivations of the same phase boundary:** R6 (`(q−1)/2 = 1 ⟺ q=3`) · R7 (`Σ(q/3)^j` diverges iff q≤3) · R7 (`X_k` grows linearly at q=3) · R8 (domination fails by a factor of k at q=3) · **R12 (`D₂ = 1 ⟺ q=3`)**.

---

## Object

For odd `q`, the **qx+1 Syracuse chain** on residues coprime to `q` mod `q^k`:

&nbsp;&nbsp;&nbsp;&nbsp;`r ↦ (q·r + 1)·2^{−v} mod q^k`, &nbsp; `v ~ Geom(1/2)` (P(v)=2^{−v}, v≥1).

Let `π_k` be its stationary distribution. **The primitive object is the level-incremental deviation** (R74's, generalized):

&nbsp;&nbsp;&nbsp;&nbsp;`‖d_k‖² := Σ_{r'} (π_k(r') − π_{k−1}(parent(r'))/q)²`, &nbsp; **`M_k := q^k‖d_k‖² = Σ_{gcd(ξ,q)=1}|μ̂_k(ξ)|²`**

For `q=3` this is exactly R75/R76's `S_k`, with `M_k^{(3)} → 7/15` (verified: its deviations reproduce R76 §10's ε_2, ε_3, ε_4 exactly).

**⚠️ Naming correction (R7).** Earlier drafts used `X_k := q^k‖π_k‖²` and called it `S_k^{(q)}`. **`X_k` is a different object** — it is the *cumulative sum* of `M_k`:

&nbsp;&nbsp;&nbsp;&nbsp;**`X_k − X_{k−1} = M_k`, exactly, at every q** — one line from the projection property, since each parent has exactly `q` coprime children whose masses sum to the parent's:
&nbsp;&nbsp;&nbsp;&nbsp;`‖d_k‖² = ‖π_k‖² − (1/q)‖π_{k−1}‖²` &nbsp;⇒&nbsp; `q^k‖d_k‖² = X_k − X_{k−1}`. ∎

This is **Pythagoras**: the lifted parent is the orthogonal projection of `π_k` onto the level-(k−1) σ-algebra and `d_k` is its complement, so `‖π_k‖² = (1/q)‖π_{k−1}‖² + ‖d_k‖²` — *(what level k−1 already knew) + (new information at level k)*. `X_k` is the accumulated energy of that decomposition.

**Use `M_k`, not `X_k`.** `M_k` carries the clean `(q/3)^k` rate at **every** odd q **including the critical q=3**; `X_k` carries it only for `q ≥ 5` (inherited by geometric summation) and **degenerates at q=3**, where the series diverges linearly and `X_k ≈ (7/15)·k`. The object that degenerates at the case of interest must not carry the headline.

---

## Result 1 (RATE) — `S_k^{(q)} ~ (q/3)^k`, universal in q

Equivalently `‖π_k‖² ~ C_q · 3^{−k}` with the exponential rate **`3^{−k}` independent of q**; the q-normalized ratio `R_k = (S_k/S_{k−1})/q → 1/3`.

**Mechanism (the "3" named).** `‖π_k‖²` contracts by exactly `1/3` per level, and

&nbsp;&nbsp;&nbsp;&nbsp;**`1/3 = Σ_{v≥1} 2^{−v}·2^{−v} = Σ_{v≥1} 4^{−v} = E_{v~Geom(1/2)}[2^{−v}]`**

— the halving second-moment. This is **q-blind**: `q` enters the transfer operator only through the multiplicative character on `(Z/q^k)*`, which rescales the state-count `q^k` but not the halving statistic. Hence `‖π_k‖² ~ (1/3)^k` and `S_k = q^k‖π_k‖² ~ (q/3)^k`. Geometric picture: the stationary measure occupies an **effective `3^k` residues** inside the `q^k`-sized space, for every q — verified `1/‖π_k‖² = 3, 9, 27` exactly at q=11,13,25.

**Evidence (Probe 5, `probe_5_universal_rate.py`).**
- **Adversarial-q falsifier ran first** and survived on every case built to break the separation: small `ord_q(2)` primes (7, 23, 47, **89 [ord 11/88]**) all `R_k → 1/3`; odd composite (9,15,25,21,27,45) all `→1/3`; q≡0 mod 3 (critical q=3) trends to 1/3. The q=7 anomaly and small-ord primes give a **clean q/3 rate** — their anomaly lives entirely in the *constant* (Result 2), never the rate.
- Direct contraction `‖π_k‖²`-ratio `→ 1/3` (q=5,7,11,13,25); exact identity `Σ4^{−v}=1/3`.
- **Scope: odd q only.** Even q breaks the construction (2 not invertible mod `q^k`).

**Mechanism, properly named (R8 — supersedes the hand-wave above).** `π_k` **is a q-adic self-similar measure.** Iterating from any `r_0`, `r_k = q^k·r_0·2^{−A_k} + Σ_{m=1}^k q^{m−1}2^{−S_m}` — **mod `q^k` the `r_0` term vanishes identically**, so `π_k` is exactly the law of `Σ_m q^{m−1}2^{−S_m}`. That is an IFS `T_v(x)=(qx+1)/2^v` with weights `p_v = 2^{−v}/Z`, and **q-adically every map contracts by exactly `1/q`** (`2^v` is a unit). Writing the L² mass by address `a=(v_1..v_k)`, `p_a = ∏_i p_{v_i}`:

&nbsp;&nbsp;&nbsp;&nbsp;**`‖π_k‖² = Σ_a p_a²` [DIAGONAL] `+ Σ_{a≠a', val(a)=val(a')} p_a p_{a'}` [OVERLAPS]**, &nbsp; `DIAGONAL = (Σ_v p_v²)^k → (1/3)^k`

**So the "3" is `1/Σ_v p_v²` — the address measure's own participation ratio — and `D₂ = log3/log q` is the measure's correlation dimension.** This makes the constant a *known species* of object (cf. Bernoulli convolutions) rather than a bespoke identity. Gated at machine zero: address enumeration vs power iteration agree to **0.000e+00** (q=3,k=2) and ≤1.1e−16 at (3,3), (5,2), (7,2).

**Open (restated correctly) — an OVERLAP bound, not a character estimate.** `‖π_k‖² = (1/3)^k` exactly iff distinct addresses land on distinct residues. Overlaps can only *increase* the sum (Cauchy–Schwarz per fiber), so **`C_q ≥ 1` is forced** — which is why every measured `δ_q > 0`; that is a theorem, not an observation. The open step is therefore:

&nbsp;&nbsp;&nbsp;&nbsp;**off-diagonal collision mass = `O(3^{−k})`** — a counting problem in `(Z/q^k)*`, not a Fourier estimate.

Measured status: **confirmed for q ≥ 5** (`offdiag/diag` converges: 0.210 / 0.361 / 0.00209 / 0.00068 at q=5/7/11/13). **At q=3 it fails by exactly a factor of k** (`ratio_k ≈ 0.4655·k`, and 0.4655 ≈ 7/15) — which is not a defect but *the critical behaviour itself* (`‖π_k‖² ~ (7/15)·k·3^{−k}`).

**Why this should be tractable** (unlike the archimedean Bernoulli-convolution analogue, open since Erdős): the metric is **non-archimedean** — q-adic balls are nested or disjoint, never partially overlapping; all maps contract by *exactly* `1/q` with no distortion; collisions are **exact algebraic coincidences** (`2^{−S}` matching mod `q^k`), not near-misses, so no transversality is needed; and the coding is **triangular** (digit m depends only on `2^{−S_m}`), so it can be analysed digit by digit.

**⛔ DEAD ROUTE — do not retry (R6).** Earlier drafts advertised "one line to a theorem: generalize R76's conservation identity to `(Z/q^k)*`." **That framing is wrong twice over.** (i) Conservation **ports verbatim and for free** — R76 Thm 76.1 is a complete-character-sum vanishing (`Σ_j e^{2πi rξj/q} = 0` unless `q|rξ`; `gcd(r,q)=1` ⇒ `q∤ξ`) that **never uses q=3**; confirmed to ≤2.9e−15 including composite q=9, 15. If it were the missing step the paper would already be done. (ii) It is **insufficient** — R76 Thm 76.3's leading-mode collapse `S_{n+1} = −2M_{n+1}(1+3^n)` needs Lemma 76.2's pairing to leave *one* unknown, and the pairing ports fine (M is palindromic on inverse pairs at every q≥5) but the **count** does not: conservation is 1 equation in `(q−1)/2` unknowns, and

&nbsp;&nbsp;&nbsp;&nbsp;**`(q−1)/2 = 1 ⟺ q = 3`.**

**Machinery triage (complete): R74 ✓ ports · R75 (Plancherel) ✓ ports · R76 Thm 76.1 ✓ ports · R76 Thm 76.3 ✗ (q=3 only).**

**★ Byproduct — a second, structural sense in which q=3 is critical.** It is the unique odd q for which conservation determines the leading mode, with a one-line proof. Combined with the geometric-series divergence at `(q/3)=1` (Object, above), the paper now has a **phase-boundary** justification for "q=3 is critical" — precisely the register the Bernoulli-convolution template (Erdős/Solomyak/Hochman) uses — rather than only the observation `M_k^{(3)} → 7/15`.

---

## Result 2 (CONSTANT) — `c̃_q = (q−3)/q`, from the same factorization

Leading `S_k/(q/3)^k → 1`, so the difference `D_k = S_k − S_{k−1}` gives `c̃_q := D_k/(q/3)^k → 1 − 3/q = (q−3)/q`. Same rate factorization; the bare `3` in `q−3` is the same `3 = 1/E[2^{−v}]`. **Rate and constant are one derivation.** Confirmed to ≤0.2% at q=11,13,17 (`c_tilde_structure_verdict.md`).

---

## Result 3 (CORRECTION) — `δ_q = 2^{1−ord_q(2)}·(q−3)/q`

> **⚠️ This result was REWRITTEN 2026-07-15 (R9/R10). The previous claim — `δ_q ≈ 0.82/ord_q(2)`, "R²=0.94, OOS-validated at q=31,127,73" — is REFUTED by a factor of 2.55×10¹³. It must not be cited.**

**The law.** `δ_q := c̃_q − (q−3)/q` and, with `d := ord_q(2)`, `x := 2^{−d}`, `M := ord_{q²}(2)`:

&nbsp;&nbsp;&nbsp;&nbsp;**`ratio_2 := offdiag_2/diag_2 = [(1+x)/(1−x)·(1−2^{−M})/(1+2^{−M}) − 1]·(1 + ε_q)`** &nbsp;→&nbsp; **`2^{1−ord_q(2)}·(1+ε_q)`**
&nbsp;&nbsp;&nbsp;&nbsp;**`δ_q = ratio_2 · (q−3)/q`** &nbsp; (exact; verified at q=41, 47)

with `ε_q ≈ 0.007` for `q ≥ 31` (erratic below). **The bracketed leading term is an EXACT IDENTITY — a theorem, verified in exact rational arithmetic (Fraction equality, not tolerance).** `ε_q` is open.

**DERIVED — the count (R11).** The structural fact: `value(v_1,v_2) = 2^{−v_2} + q·2^{−(v_1+v_2)} mod q²`, and **the second term carries a factor `q`, so it needs `A = v_1+v_2` only mod `d`** (since `2^d ≡ 1 mod q`). Hence

&nbsp;&nbsp;&nbsp;&nbsp;**`v_1 → v_1 + d` leaves the value EXACTLY unchanged — `v_1` is only ever determined mod `d`.**

(Verified by exact integer equality at q=11, 13, 17, 31, 41, 47.) So collisions are **not rare events that "cost" a period shift — they are structural**: every value-bucket contains a whole **geometric tower** in `v_1`. Summing the tower (`G_c = 2^{−c}/(1−x)`, `H_c = 4^{−c}(1+2^{−M})/((1−x²)(1−2^{−M}))`, `Σ_c 4^{−c} = (1−x²)/3`, `P2 = (1/3)(1+2^{−M})/(1−2^{−M})`) gives the bracketed identity above.

> **★ The prefactor 2 is the tower's CROSS-TERM: `(1+x)/(1−x) − 1 = 2x/(1−x)`.**

**⚠️ Correction to earlier drafts:** the previously-stated reasons — "the cheapest collision costs a full period shift" and "`offdiag` sums ordered pairs ⇒ counted twice ⇒ prefactor 2" — **were wrong.** They produced the right number for the wrong reason. The 2 comes from the geometric tower's cross-term, and the collisions are structural, not cheap. `ratio_1 = 0` (the k=1 bijection onto `⟨2⟩`, measured at machine zero on 12 primes) remains correct and is what forces the tower to first appear at k=2.

`1/ord_q(2)` has **no mechanism** — nothing in the structure produces a reciprocal.

**Evidence — predictions pre-committed before the run, on primes far outside the fitted range (ord 3–12):**

| q | ord_q(2) | measured `ratio_2` | `2^{1−ord}` predicts | miss |
|---|---|---|---|---|
| 41 | 20 | 1.9199e−06 | 1.9073e−06 | **0.7%** |
| 47 | 23 | 2.4011e−07 | 2.3842e−07 | **0.7%** |
| 59 | 58 | **0** (−2.5e−16, sub-eps) | 6.9e−18 | consistent |

Free 2-parameter fit `ratio_2 = a·c^{−ord}` on the *old* primes recovers **`c = 2.01704`** against a mechanism that named **2** before fitting; the large-q primes pin the prefactor at **2.0132, 2.0142**.

`δ_q = 2^{1−ord_q(2)}·(q−3)/q` checked against the old primes, which never informed it: **q=73 → 0.1% · q=41 → 0.7% · q=47 → 0.8% · q=127 → 0.8% · q=31 → 4% · q=11 → 8% · q=17 → 12% · q=7 → 47% · q=13 → 49% · q=5 → 84%** — sub-1% at large q, degrading at small q. **The small-q degradation is `ε_q` (family b), whose solvability condition is genuinely arithmetic — see "Still owed" below. It is NOT an `O(1/q)` decay** (R11): it plateaus at ~0.7% across q=31→47 rather than continuing to shrink.

**★ Result 3 carries Result 2's `(q−3)/q` factor — Results 2 and 3 are not independent.** With `X_k − X_{k−1} = M_k` linking Results 1↔2 (Object) and the `3 = 1/Σ_v p_v²` naming (Result 1), **all three results now flow from the single self-similar-overlap mechanism.**

**Still owed (flagged, not claimed) — narrowed by R11.** The leading term is now an exact identity and the prefactor 2 is derived. What remains is **family (b)**: collisions with `v'_2 = v_2 + jd` and a compensating `v_1`, which exist iff

&nbsp;&nbsp;&nbsp;&nbsp;`2^{−A'} ≡ 2^{−A} + j·t·2^{−v_2} (mod q)` &nbsp; is solvable — i.e. iff the RHS lands in `⟨2⟩`, where `t = (2^d−1)/q mod q`.

That is a **genuinely arithmetic** condition, which is why small q is erratic. Family (b) contributes `ε_q`: measured `resid/family(a)` = 0.058 (q=11), **0.382 (q=13)**, 0.097 (q=17), 0.0068 (q=31), 0.0066 (q=41), 0.0071 (q=47) — positive everywhere (consistent with the proved `δ_q > 0`), erratic at small q, then **flat at ~0.7% across q=31→47**.

**⚠️ Correction:** earlier drafts called this an `O(1/q)` correction. **That is unsupported** — it plateaus rather than decays. Restated as `(1+ε_q)`. Honest limit: three large-q points cannot distinguish a plateau from slow decay. Mechanically, family (b) is *also* `O(x)` (its pairs cost `x^j`), which would make `ε_q` a constant rather than `~1/q`.

**Why the old claim survived as long as it did — a methodological note worth keeping.** Recomputing the old fit gives `δ = −0.08903 + 0.81886/ord` with **linear R² = 0.94045**, reproducing the published slope 0.82 and R² 0.94 to three digits — **while being off by 121.5× at q=13, 53.4× at q=11, 24.3× at q=73.** δ is convex and spans 350×, so a linear fit is dominated by its two largest points while the rest cluster near the origin where `1/ord` is also small. **Linear R² cannot discriminate monotone candidates across that range.** Worse, the fitted law has intercept `−0.089`, so it predicts **δ < 0 for ord > 9.2** — contradicting the *proved* positivity `δ_q > 0`. And the original "OOS validation" used q=31, 127, 73 with ord **5, 7, 9** — **interpolation inside the fitted range 3–12.** It tested the fit; it never tested the functional form.

---

## Coherence note

The number `3` recurs as `1/E_{Geom(1/2)}[2^{−v}]` — a purely 2-adic halving statistic — and is the reason the rate is q-independent. It is **not** the multiplier (which is q, and the rate is q-independent). This is the same `1/3 = E[2^{−v}]` that governs Tao's 2-adic renewal; here it is isolated as the universal contraction rate of the stationary L² mass across the entire qx+1 family.

## Provenance / files

- `probe_5_universal_rate.py` + `result_5_universal_rate.md` + `result_5_data.csv` — Result 1 derivation + falsifier (this session).
- `PRE_REG_5_UNIVERSAL_RATE_2026_07_14.md` — pre-registration.
- `result_q_sweep_test_1_rate.md` — original empirical rate observation (2026-05-04).
- `c_tilde_structure_verdict.md` + `c_tilde_q17_probe.py` — Result 2.
- `result_4_ctilde_ord2.md` + `.py` + `PRE_REG_4_CTILDE_ORD2` — Result 3.
- `c_seven_forty_fifth.md` (R75), `result_76_conservation_law.md` (R76) — the q=3 Plancherel/conservation machinery the domination identity must generalize.

**Next step to finish the paper (revised 2026-07-15).** ~~Generalize R76's conservation~~ — **DEAD, see Result 1's ⛔ DEAD ROUTE.** The remaining input is a **collision-counting bound**, not a conservation collapse:

&nbsp;&nbsp;&nbsp;&nbsp;**bound the off-diagonal collision mass by `O(3^{−k})` on `(Z/q^k)*`**, upgrading Result 1 from mechanism to theorem.

Concretely tractable next moves, in order:
1. ~~**The full collision count at k=2.**~~ **DONE (R11).** Family (a) is counted: the leading term is an exact identity and the prefactor 2 is derived (the geometric tower's cross-term). **What remains is the family-(b) count** — the compensating-`v_1` solvability condition `2^{−A'} ≡ 2^{−A} + j·t·2^{−v_2} (mod q)`, solvable iff the RHS lands in `⟨2⟩`. Counting that would deliver `ε_q` and turn Result 3 into a theorem outright. Still the cheapest real win on the board (family (a) took 0.18 s).
2. **Read Tao 2019** (*Almost all orbits of the Collatz map attain almost bounded values*) before deriving anything — it proves Fourier decay for `Syrac(Z/3^n)` and may already contain most of the k-fold overlap estimate. Re-deriving it would be this arc's characteristic mistake.
3. **p-adic Bernoulli convolution literature** — a direct sibling never covered by the 2026-05-04 lit dive (which looked at archimedean BCs and Siegel). Non-archimedean self-similar measures are exactly this object.
4. **Does R77/`T_lead` port?** Mod q the chain gives `r_out ≡ 2^{−v}`, valued in `⟨2⟩` of size `ord_q(2)`. At q=3, `⟨2⟩={±1}` and `2^{−v} ≡ (−1)^v` — **R64.B/R66's v-parity class rule is the q=3 shadow of an `ord_q(2)`-class structure.** Porting R77 would complete the machinery triage.

**Read Siegel before publishing** ((p,q)-adic Analysis and Collatz, USC 2022 / Springer 2024–25) — closest sibling, not a duplicate.
