# Phase 3 — Re-entry Recursion Worksheet (scaffold for L3)

**Purpose.** L3 is the one open step: prove `d = ord_q(2) ≥ 3 ⟹ r_q < 1` (spectral gap), with `r_3 = 1` at the boundary. This document assembles the **level-`m → m+1` re-entry recursion** for the cross-cell correlation mass, with the two known collapse points — the phase-spread factor (fails at `d=2`) and the shift-coupling constant `s_R13` (fails at `s≥2`) — placed **explicitly and separately**, so that any closed form you write for the recursion can be checked against both fixed points at a glance. It is a scaffold, not a proof. Everything is tagged **[DERIVED]** (exact, on paper + gate-validated), **[NUM]** (measured), or **[OPEN]** (you supply this).

Written 2026-07-16, after R13/R14/R20/R35/R36. Companion to `PHASE3_SPECTRAL_GAP_BRIEF.md` (the L1–L4 skeleton) and `PHASE3_L3_REDTEAM.md` (the attack surface).

---

## 0. Notation (fixed once)

- `q` odd prime, `d = ord_q(2)`, `H = ⟨2⟩ ⊆ F_q^*`, `|H| = d`.
- `s := v_q(2^d − 1) ≥ 1` (the LTE constant). `s = 1` for every prime `< 1093`; `s ≥ 2` first at `q = 1093, 3511` [NUM, R34/R35].
- `s_R13 := (2^d − 1)/q mod q ∈ F_q`. Note `s_R13 = 0 ⟺ q² | 2^d − 1 ⟺ s ≥ 2`. **This is the single scalar that carries the `s≥2` boundary.** [DERIVED, R13]
- Address `= (v_1, …, v_k)`, `v_i ≥ 1`, weight `p_{v_i} = 2^{−v_i}`, `Σ_v p_v = 1`, `Σ_v p_v² = 1/3`.
- Suffix sums `S_m := v_{k−m+1} + … + v_k` (so `S_1 = v_k`, `S_k = v_1+…+v_k`). [DERIVED, R14/R20]
- `f(k) := ‖π_k‖²`. Value of an address mod `q^k` is `value = Σ_{m=1}^{k} q^{m−1} 2^{−S_m} (mod q^k)` [DERIVED, R14].
- `x_j := 2^{−d q^{j−1}}` — the tower ratios, doubly-exponentially separated (`x_1 = 2^{−d}`, `x_2 = 2^{−dq}`, …) [DERIVED, R14].

---

## 1. The object: the increment IS the fresh off-diagonal mass

Split a correlation pair by its **leading** coordinate `(v_1, v'_1)`:

- `v_1 = v'_1` (leading factor `Σ_v p_v² = 1/3`): the pair reduces to a depth-`(k−1)` pair with the *identical* collision requirement. This is the **diagonal branch** and the entire source of `λ_1 = 1/3`.
- `v_1 ≠ v'_1`: a fresh **off-diagonal** contribution `g(k)`.

So the exact one-step recursion is [DERIVED, R8]:

```
    f(k) = (1/3)·f(k−1) + g(k),        g(k) = off-diagonal (v_1 ≠ v'_1) mass surviving to depth k.
```

Normalize by the diagonal rate. Let `F(k) := f(k)/(1/3)^k = 1 + within(k) + cross(k)`. Then

```
    F(k) = F(k−1) + G(k),        G(k) := 3^k · g(k).
```

`within(k)` is **k-flat** (`≈ 2·2^{−d}`, every tower past the first is `≤ 2^{−dq}` ≈ nonexistent) [DERIVED, R14]. Therefore the increment

```
    c_k := cross(k) − cross(k−1) = G(k)   (up to the doubly-exponential within-tail, R14/R33 Attack 5).
```

> **`c_k = G(k) = 3^k g(k)` is the object of the recursion. `r_q = lim_k c_{k+1}/c_k` is the subdominant eigenvalue** [NUM, matches operator `M`, R30/R32]. Result 1 ⟺ `r_q < 1`.

---

## 2. The cascade gates (what "survives to depth `k`" means)

A pair with leading shift `v_1 ≠ v'_1` collides mod `q^k` iff `Σ_{m=1}^{k} q^{m−1} T_m ≡ 0 (mod q^k)`, `T_m := 2^{−S_m} − 2^{−S'_m}`. Dividing by `q` successively gives **`k` gates**, one per q-adic digit [DERIVED, R20]:

| gate | exact condition | pass-rate [NUM, R20] |
|---|---|---|
| **1** | `T_1 ≡ 0 (mod q)` ⟺ `S'_1 = S_1 + j_1·d` (shift is a multiple of `d`) | `≈ 1/d` |
| **2** | `2^{−S'_2} ≡ 2^{−S_2} + j_1·s_R13·2^{−S_1} (mod q)` | `≈ 1/q` |
| `m ≥ 3` | `W_{m−1} + T_m ≡ 0 (mod q)`, `W_{m−1} := (U_{m−2}+T_{m−1})/q` | `≈ 1/q`, **stabilized** (`L2≈L3≈L4`) |

Three facts to hold onto:

- **Gate 2 is where `s_R13` enters, and it is the ONLY place a single power of `s_R13` appears linearly.** [DERIVED, R13/R20] Set `s_R13 = 0` and gate 2 degenerates to `2^{−S'_2} ≡ 2^{−S_2}`, i.e. `S'_2 ≡ S_2` — the coupling to the level-1 shift `j_1` vanishes. This is fixed point 2 (§5).
- **The level-`≥3` gate is a definition, not yet a formula.** `W_2 = (U_1 + T_2)/q` needs the 2nd-order expansion of `2^{−jd} mod q³`. **[OPEN — R20 §"honest scope"]** You supply this; the analog of `s_R13` at the next q-power (call it `s^{(2)} := v_{q}((2^d−1)/q^{s})`-type second digit) should appear here.
- **Pass-rates are `≈1/q` and q-uniform, NOT `d`-dependent.** [NUM, R20/R36] R36 measured the *conditional* re-collision ratio `ρ = r_2/r_1 ≈ 0.75` **uniformly across q=3,5,7,11,13** — it does **not** discriminate `q=3`. ⇒ **The counting/rate side of the recursion is blind to the boundary. Both fixed points must live on the weight/phase side.** This is the load-bearing empirical constraint (§6).

---

## 3. The re-entry recursion, factored

Write the fresh off-diagonal mass as a **transfer of the previous level's mass through one prepended coordinate pair**. Group the surviving pairs by their **phase class** — the residue of the leading phase `2^{−S_1} mod q` (equivalently the `H`-coset structure of the shift). Then, schematically [DERIVED structure / OPEN weights]:

```
    G(k)  =  3 · Σ_{phase classes φ}  κ_φ · ω_φ(d, s_R13) · G(k−1|φ)
             └┬┘  └─────┬─────┘   └──┬──┘   └────┬────┘
       diagonal    rate factor    phase/weight    previous-level mass
       normalization  (§2 gates)     factor        in class φ
```

Two multiplicative factors per level, and **the two fixed points live in different ones**:

- **`κ_φ` = rate factor.** How many `(v_1, v'_1)` shifts feed phase class `φ` while passing the gate. Governed by the gate pass-rates (`1/d` at level 1, `1/q` after). **q-uniform, boundary-blind** [NUM, R36]. This is *not* where the gap opens or closes.
- **`ω_φ(d, s_R13)` = phase/weight factor.** The `2^{−v}`-weighted amplitude that phase class `φ` contributes — i.e. how much *mass* (not count) re-enters. This is where `d` (via the phase spread of `2^{−v} mod q` over `H`) and `s_R13` (via gate 2's coupling) both live.

> **The recursion's asymptotic ratio is `r_q = 3 · ρ(κ·ω)`**, the spectral radius of the one-level transfer operator `(κ·ω)` — which is exactly the subdominant eigenvalue of `M|Krylov(v₀)` [NUM]. The claim `r_q < 1` for `d ≥ 3` is a statement that **`3·κ·ω` contracts off the diagonal** — and by the bullet above, all the `d`-dependence is inside `ω`.

---

## 4. Fixed point 1 — `d = 2` (q=3): the phase-spread factor saturates

At `d = 2`, `H = ⟨2⟩ = {1, −1}` and `2^{−v} ≡ (−1)^v (mod 3)` — the phase takes **only two values**. Consequences:

- The phase class `φ` collapses to **parity of `v`**. A shift `v_1 ≠ v'_1` lands in the collision-supporting class whenever `v_1, v'_1` share parity — a **fixed positive fraction at every depth**, not a decaying one.
- So `ω_φ` does **not** decay with `k`: the off-diagonal mode stays aligned with the diagonal, `λ_2 = λ_1 = 1/3`, and `3·κ·ω = 1`. **`r_3 = 1`, `cross(k)` linear.** [NUM/DERIVED, R8/R32]
- **The collapse is NOT in `κ`.** R36 is the direct check: `κ`'s conditional ratio `ρ ≈ 0.75` is the same at `q=3` as at `q=5,7,11,13`. If the boundary lived in the rate factor, `q=3` would be an outlier in `ρ`; it is in the pack. **⇒ fixed point 1 is a property of `ω` alone.** [NUM, R36]

**Constraint FP1:** a correct `ω(d, ·)` must satisfy `3·κ·ω → 1` as `d → 2`, and this must be visibly a **phase-count** collapse (`H` disperses `2^{−v} mod q` over only `|H|=2` values ⇒ the "no-collision" fraction `Pr[d ∤ m]` vanishes), *not* a rate collapse.

---

## 5. Fixed point 2 — `s_R13 = 0` (s ≥ 2): the shift-coupling collapses to an index shift

At `s ≥ 2` (`q = 1093, 3511, …`), `s_R13 = 0`, and gate 2 (§2) degenerates:

```
    2^{−S'_2} ≡ 2^{−S_2} + j_1·s_R13·2^{−S_1}  ⟶  2^{−S'_2} ≡ 2^{−S_2}   (mod q).
```

The coupling of level 2 to the level-1 shift `j_1` vanishes. What this does to the mass [NUM, R35]:

- Pre-registered guess (mine) was **inflation** — depth-2 mass jumping to depth-1 scale. **WRONG.** Measured `cross(2)|_{q=1093} = exactly 0` — the folded cells (`ord_{q²}(2) = d`, so the address collapses harder) are **injective at depth 2**, colliding nothing.
- Correct reading: the collapse **folds addresses into the same cell rather than colliding distinct ones**, so the overlap **onset is delayed by one level**: `cross(2)|_{s=2} = cross(1)|_{s=1} = 0`, and the whole structure index-shifts down, `cross(k)|_{s=2} ≈ cross(k−1)|_{s=1}`.
- **Consequence: the gap survives, index-shifted. The `s≥2` side condition is likely REMOVABLE** — L3's statement simplifies (no side condition) rather than gaining one. Caveat: established at depth 2 only; `cross(3)|_{s=2}` (≈ `d³q²` cells) is infeasible, so the shift-down is *suggested*, not proven [NUM, R35].

**Constraint FP2:** a correct `ω(d, s_R13)` must have `s_R13` entering such that `s_R13 → 0` produces an **index shift `k → k−1`** in `G(k)` (onset delay), **not** a factor blowup and **not** `r_q → 1`. Concretely: `G(k; s_R13=0) = G(k−1; s_R13≠0)` to leading order. The `s_R13` dependence is a *phase offset*, not an *amplitude*.

---

## 6. The two constraints, side by side (the check for your closed form)

Any closed form `r_q = 3·κ·ω(d, s_R13)` you derive must pass **both** simultaneously:

| | trigger | must produce | must NOT be | evidence |
|---|---|---|---|---|
| **FP1** | `d = 2` | `3·κ·ω → 1` (divergence), via **phase-count** saturation in `ω` | a change in `κ` (rate) | R36: `ρ` q-uniform |
| **FP2** | `s_R13 = 0` | **index shift** `k→k−1` in `G(k)` (onset delay); `r_q` unchanged | inflation, or `r_q → 1` | R35: `cross(2)=0` |

And the three sanity anchors from the arc:

- **A1.** `r_q` has **no elementary closed form** — `3/q` is refuted (`r_5≈0.62 > 3/5`, `r_7≈0.38 < 3/7`) [NUM, R28]. So `ω` is a **character sum** (R13 §"the route"), not a rational function of `d, q`. Expect square-root cancellation, not a formula. Do not fit six points.
- **A2.** The mechanism is **`d=2` specifically**, not "`H` proper." At `q=5`, `2` is primitive (`H` full) yet `r_5 < 1` [NUM]. `ω` must gap for *all* `d ≥ 3` including small `d` (`d=3` at `q=7`), where subgroup-sum bounds are vacuous. So `ω`'s contraction is **qualitative non-degeneracy for `d ≥ 3`**, not an asymptotic bound.
- **A3.** `P_coll(d) = (2/3)/(2^d−1)` does **NOT** match the measured level-1 rate `r_1` [NUM, R36] — do not build `ω` on it. Use the measured/derived gate rates (`1/d` then `1/q`), which are the validated numbers.

---

## 7. `ω` as a character sum (the concrete target for L3)

R13 already reduced the family-(b) mass to a character sum, and that is exactly `ω`. Expanding the `⟨2⟩`-membership indicator over the `(q−1)/d` characters trivial on `H` [DERIVED, R13]:

```
    1[y ∈ ⟨2⟩] = (d/(q−1)) · Σ_{χ : χ|_H = 1} χ(y).
```

Applied to gate 2's condition `y = 2^{−S_2} + j_1·s_R13·2^{−S_1}`, the family-(b) weight is

```
    ω  =  (d/(q−1)) · Σ_{χ|_H = 1}  Σ_{v, v', j_1}  p_v p_{v'} · χ(2^{−S_2} + j_1·s_R13·2^{−S_1}) · [weights].
```

- **Principal character** (`χ = 1`) gives the "expected" density `d/(q−1)` — the smooth part.
- **Non-principal characters** give the fluctuation (the `14×` swing between `q=17` and `q=41` at identical `d/(q−1)`, R13). This is the erraticness of `ε_q`, and it is **structural** (Jacobi/Kloosterman species), not noise.

> **L3 in this register:** show the character sum `ω` has **spectral radius `< 1/3` for `d ≥ 3`** (gap) and `= 1/3` at `d = 2` (boundary). The `d=2` boundary is anchored in literature (Konyagin small-subgroup non-cancellation, the `|G|=2` example `S(1,{±1}) = 2cos(2π/q) = |H| + O(q^{−2})`; Siegel pole at `q∈{1,3}`). The **`d≥3` positive side has no off-the-shelf route** (bounds vacuous for small `d`) — it is the novel content, and it must come from the **structure of this specific character sum**, exploiting that `s_R13 ≠ 0` for `d ≥ 3, s = 1` makes gate 2 genuinely couple (unlike FP2).

---

## 8. Derived vs. what you supply (honest ledger)

**Derived / gate-validated (lean on these freely):**
- The one-step diagonal recursion `f(k) = (1/3)f(k−1) + g(k)`; `c_k = G(k) = 3^k g(k)` [R8].
- The value formula, suffix-sum grading, tower ratios `x_j` [R14].
- Gate 1 (`S'_1 = S_1 + j_1 d`) and gate 2 (`2^{−S'_2} ≡ 2^{−S_2} + j_1 s_R13 2^{−S_1}`) as **exact iff conditions** [R13/R20, ~2M pairs, zero failures].
- The character-sum reduction of `ω` [R13].
- Level rates `1/d` (L1), `1/q` (L≥2, stabilized) [R20]; `ρ ≈ 0.75` q-uniform [R36].
- Both fixed points as measured: FP1 `r_3=1` phase collapse [R32]; FP2 `cross(2)|_{1093}=0` index shift [R35].

**[OPEN] — the pen-and-paper work (yours):**
1. **The level-`≥3` explicit gate** — expand `W_2 = (U_1+T_2)/q` via `2^{−jd} mod q³`. Expect the next-digit analog of `s_R13` to surface. Without this the recursion is exact only for the `k=1→2` step and *structurally extrapolated* (rates stabilize) beyond.
2. **The functional form of `ω(d, s_R13)`** — the character sum of §7, bounded for `d ≥ 3`. This is L3 proper. Strategy menu (from the brief): **S1** direct off-diagonal norm bound `< 1/3` via the `d∤m` fraction; **S2** renormalization `g(k) ≤ ρ_q·(1/3)·g(k−1)`, `ρ_q < 3`; **S3** finite-dimensional descent on the amplitude-carrying (Krylov) subspace. The recursion here is the raw material for all three.
3. **Confirm FP2 is truly removable** — the index-shift `cross(k)|_{s=2} = cross(k−1)|_{s=1}` is depth-2 evidence only. A symbolic argument (not compute) that the shifted cascade reproduces the same `r_q` would retire the side condition for good.

---

## 9. One-line summary

`c_k = 3^k g(k)` obeys a re-entry recursion `G(k) = 3·Σ_φ κ_φ·ω_φ(d,s_R13)·G(k−1|φ)`; the **rate factor `κ` is q-uniform** (R36) so **both boundaries live in the phase/weight factor `ω`** — `d=2` saturates its phase-count (FP1, `r_3=1`), `s_R13=0` index-shifts it (FP2, gap survives). `ω` is a character sum (R13); L3 = "this character sum has spectral radius `< 1/3` for `d ≥ 3`," boundary-anchored in Konyagin/Siegel, positive side novel.
