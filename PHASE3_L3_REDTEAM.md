# L3 Red-Team Brief — the spectral gap for `d ≥ 3`

**Purpose.** A brief written to be *attacked*. L3 is the one unproven step in Result 1 (everything else is gate-validated or citable — see `PHASE3_SPECTRAL_GAP_BRIEF.md`). This document states the claim as sharply as possible, isolates every assumption, and ranks the **attack surface** so a red team can go straight at the weak joints. Written 2026-07-16. `[NUM]` = numerically established; `[ASSUME]` = load-bearing assumption not yet proven; `[GAP]` = the actual hole.

---

## 1. The claim, stated three ways (attack the statement first)

Let `q` be an odd prime, `d = ord_q(2)`, `H = ⟨2⟩ ⊆ F_q^*`, `π` the q-adic self-similar measure of the IFS `T_v(x)=(qx+1)2^{-v}` (weights `p_v=2^{-v}`), `π_k` its image mod `q^k`.

- **(Analytic)** The transfer operator `L` of `π`, restricted to the Krylov subspace of the initial correlation vector, has a spectral gap for `d ≥ 3`: its second eigenvalue `λ_2` satisfies `|λ_2| < λ_1 = 1/3`. Define `r_q := |λ_2|/λ_1`; then `r_q < 1`.
- **(Probabilistic)** With `X, X'` iid `∼ π`, and `f(k) := P(X ≡ X' mod q^k) = ‖π_k‖²`, the normalized excess `cross(k) := f(k)/(1/3)^k − 1 − within(k)` is **bounded in k** for `d ≥ 3` (and grows linearly for `d=2`).
- **(Boundary)** `d = 2 ⟺ q = 3 ⟺ 2 ≡ −1 (mod q)`; there the gap **closes** (`r_3 = 1`, linear growth). The theorem is precisely the implication `d ≥ 3 ⟹ r_q < 1`, **uniform in the truncation depth** used to define `L`.

**First attack:** is the Krylov restriction legitimate, or does it quietly discard mass that matters? (See §6, Attack 3.)

---

## 2. Non-negotiable facts — any proof MUST reproduce all of these

A proposed proof that violates any of these is wrong on arrival:

1. `λ_1 = 1/3 = Σ_v p_v²` exactly, q-independent. `[NUM/exact, R8]`
2. Gap closes **exactly** at `d = 2` (q=3), open for **every** `d ≥ 3` tested (q=5,7,11,13,19,29). `[NUM]` So the argument must (a) use `d ≥ 3`, and (b) **break at d = 2** — a proof that doesn't visibly fail at d=2 is proving the wrong thing.
3. `r_q` has **no elementary closed form** (`3/q` refuted: `r_5≈0.62>3/5`, `r_7≈0.38<3/7`). `[NUM, R28/R30]` So the proof yields a **bound/existence**, never a formula.
4. The mechanism is **NOT** "H is a proper subgroup." At q=5, `2` is a primitive root (`H = F_5^*` full) yet `r_5 < 1`. `[NUM]` Any argument keyed to `H ⊊ F_q^*` is dead (kills q=5).
5. The near-`λ_1` "tower" eigenvalues carry **exactly zero amplitude** (orthogonal to Krylov(v₀)). `[NUM, R32]` They are real eigenvalues but dynamically invisible; the gap is `λ_1 → r_q·λ_1`, not `λ_1 → tower`.
6. Subgroup exponential-sum bounds are **vacuous for small d** (Konyagin Thm 1.7/3.3 need `|H|>√q`/`q^δ`; Thm 1.8 gives *no* cancellation for `|H| ≪ log q`). `[gather]` So no off-the-shelf `|S(a,H)| < |H|` will carry it.

---

## 3. Where the contraction comes from (the mechanism, as far as we understand it)

**The difference process + Lifting-the-Exponent.** Write `X = 2^{-v}(1 + qX_1)` (self-similarity, `X_1` iid copy). For two copies with leading coordinates `v, v'`:
```
   D := X − X' = (2^{-v} − 2^{-v'}) + q·(2^{-v}X_1 − 2^{-v'}X_1').
```
- **Diagonal branch `v = v'`** (prob `Σp_v² = 1/3`): `D = 2^{-v}·q·(X_1 − X_1')`, so `v_q(D) = 1 + v_q(X_1−X_1')` — valuation advances by 1 and **recurses**. This is the entire source of `λ_1 = 1/3`: `f(k) = (1/3)f(k−1) + [off-diagonal]`.
- **Off-diagonal branch `v ≠ v'`** (prob `2/3`): let `m = v'−v`. The leading term `2^{-v}−2^{-v'} = 2^{-v}(1 − 2^{-m})` has
  ```
     v_q(1 − 2^{-m}) = v_q(2^m − 1) = 0                if d ∤ m   (NO collision mod q),
                     = s + v_q(m/d)                   if d | m,  s := v_q(2^d − 1) ≥ 1  (LTE).
  ```
  So **a differing pair collides to depth k only if its shift `m` is a multiple of `d` that is q-adically deep**: `q^{k−s} | (m/d)`. Deep shifts are tower-rare (geometric coordinate weights), and this rarity is the contraction.

**Why `d` is the knob (the heart — and the softest part).** At `d = 2`, `H = {1, −1}` and `2^{-v} ≡ (−1)^v (mod q=3)`: the phase takes only **two** values, so `v ≠ v'` still collides mod 3 whenever `v, v'` share parity — collisions are *too easy*, the off-diagonal doesn't decay relative to the diagonal, `r_3 = 1`. At `d ≥ 3` the phase `2^{-v} mod q` spreads over `≥ 3` values, the `d ∤ m` case (no collision) captures a positive fraction, and the surviving collisions require deep shifts → strict contraction. **This is a heuristic, not a proof** — see Attack 1.

**Dual description (R6).** q=3 is *also* the unique prime where the Fourier conservation law `Σ_{j} M(η₀ + j q^k) = 0` **determines** the leading off-diagonal mode (`(q−1)/2 = 1`: one equation, one unknown), forcing `λ_2 = λ_1`. For `q ≥ 5`, `(q−1)/2 ≥ 2` leaves the mode **underdetermined** — free to be smaller. R6 verified the modes are genuinely *not* forced equal for q≥5. **But "underdetermined ⟹ strictly smaller" is exactly the missing implication.**

---

## 4. Candidate proof strategies (pick one to harden, or break)

- **(S1) Direct off-diagonal norm bound.** Bound the operator norm of `L` on the off-diagonal (Krylov-minus-Perron) subspace by `< 1/3`, using the LTE tower structure of §3: sum the off-diagonal collision mass over first-difference position and shift `m`, show it is a strict `θ_q < 1` fraction of the diagonal per renormalization step. *Cleanest if it works; the d-dependence must enter through the `d∤m` fraction.*
- **(S2) Renormalization / two-scale.** From `f(k) = (1/3)f(k−1) + g(k)`, show the off-diagonal `g(k)` is itself contracted: `g(k) ≤ ρ_q · (1/3)·g(k−1)` with `ρ_q < 3`. Equivalent to a Lyapunov/Doeblin–Fortet inequality for `L`. Ties directly to Ruelle quasi-compactness (L1/L2 gather).
- **(S3) Finite-dimensional descent.** The tower modes have zero amplitude (fact 5); prove that the *amplitude-carrying* subspace is finite-dimensional (or admits a uniform-in-L basis), reducing L3 to a finite eigenvalue inequality per q — which then needs a q-uniform argument for the `< 1/3` bound. *Riskiest, but closest to what the numerics actually compute.*

---

## 5. Ruled-out routes (do not spend time here)

- Subgroup exponential sums for strictness (vacuous, fact 6).
- "H proper ⟹ additively non-closed ⟹ contraction" (false at q=5, fact 4).
- Any finite linear recurrence for `c_k` (none exists — the tower, R23/R30).
- A closed form for `r_q` (fact 3).
- Importing the Bernoulli-convolution Pisot/Salem phase boundary (vacuous in the ultrametric non-overlap-in-the-usual-sense setting — gather).

---

## 6. THE ATTACK SURFACE — ranked. Hit these.

**Attack 1 (the crux — `d=2` vs `d≥3` threshold).** The mechanism in §3 is a *heuristic* for why `d≥3` contracts. **Is it actually true that `d ≥ 3 ⟹ strict contraction`, or could some `d ≥ 3` fail?** We have only 6 primes `[NUM]`. Red-team: (a) find a candidate `q ≥ 5` where the "deep-shift rarity" argument might *not* beat the diagonal — e.g. large `s = v_q(2^d−1) > 1` (Wieferich-like: `q=1093, 3511` have `q² | 2^{q−1}−1`; do any give large `s` and a near-1 `r_q`?). (b) Is the `d∤m`-fraction argument actually monotone in `d`, or is there a `d` where it's degenerate? **This is the most likely place the theorem is *false as stated* or needs a side condition.**

**Attack 2 (uniformity in L / the limit operator).** Everything numeric is at finite truncation depth `L`; the operator `M_L` is only faithful to `k ≤ L`, and `λ_2(M_L)/λ_1` was **L-dependent and value-biased** (`q=7`: 0.475 at L=2 vs 0.39 true). `[ASSUME]` The claim is about the `L→∞` limit `L`. Red-team: does the gap *persist* in the limit, or could `λ_2(M_L) → λ_1` as `L→∞` for some `q≥5` (gap closing in the limit, like q=3 does)? What guarantees quasi-compactness of the limit operator (`[GAP]`, = Hole 5)? A skeptic should demand the essential spectral radius of `L` is proven `< 1/3`.

**Attack 3 (the Krylov restriction).** We *defined away* the tower by restricting to `Krylov(v₀)` (fact 5). `[ASSUME]` Is the tower amplitude *exactly* zero in the limit, or only `~10⁻¹³` at finite L (i.e. we're reading a truncation artifact)? If the tower has tiny-but-nonzero amplitude that *grows* with L, it could re-enter and there's no gap. Red-team: prove the orthogonality `⟨tower, Krylov(v₀)⟩ = 0` structurally, or exhibit its L-scaling.

**Attack 4 (reducibility).** `M` is **reducible** (R33: many transient SCCs feed one dominant SCC). `[NUM]` The Perron/gap argument was applied to "the dominant class." Red-team: do the transient SCCs contribute long-lived transients that masquerade as slow decay (fattening `r_q` toward 1)? Is the dominant-SCC restriction spectrally clean, or does the reducible structure leak?

**Attack 5 (the `within` sleight-of-hand).** We use "`cross` increments = `total` increments" because `within(k)` is flat. `[NUM]` It is flat only to `~10⁻³⁰` (R33) — a doubly-exp tail, not exact. Red-team: could the `within` tail, summed over *all* k, contribute an O(1) shift that changes whether `cross` is bounded? (Almost certainly not — geometric-of-doubly-exp — but state the bound.)

**Attack 6 (does `s > 1` break the LTE tower?).** §3 assumed `s = v_q(2^d−1) = 1` implicitly in "collision needs `q^{k−s}|m/d`". For primes with `s ≥ 2` the tower is shifted. Red-team: recompute the boundary and the contraction for `s ≥ 2` primes; does the `d=2` characterization survive, or is it really an `s`-and-`d` joint condition?

---

## 7. Falsifiable checks a red team can run (cheap, decisive)

- **Wieferich/`s≥2` primes:** compute `r_q` (best-effort, operator L=2 modal) for `q = 1093, 3511` (if feasible) or any prime with `s = v_q(2^d−1) ≥ 2`, and for a spread of `d`. If any `q ≥ 5` shows `r_q → 1`, **Attack 1 wins and L3 is false as stated.**
- **L-limit probe:** for `q = 5, 7`, push the operator to L=3 (q=5 L=3 ≈ 1.25M — feasible) and confirm `λ_2/λ_1` (amplitude-selected) *converges* rather than drifting toward 1. Tests Attack 2/3.
- **Tower amplitude vs L:** track the tower mode's `|A|` at L=1,2,3 (q=3,5). If it is not monotonically → 0, Attack 3 has teeth.

---

## 8. Honest one-paragraph summary for the red team

We are confident `d ≥ 3 ⟹ r_q < 1` **empirically** (6 primes, direct + operator, cross-validated) and we know the boundary is `d = 2` (q=3), doubly-characterized (phase-too-coarse; conservation-determines-mode). We have a **mechanism sketch** (difference process + LTE tower: collisions need q-adically deep shifts, rare for `d ≥ 3`) and three proof strategies (S1 norm bound / S2 renormalization / S3 finite descent). The **theorem could still be false or need a side condition** at `s ≥ 2` / Wieferich primes (Attack 1+6), the **limit/uniformity is assumed not proven** (Attack 2+3 = Hole 5), and the **`d=2 vs d≥3` threshold is a heuristic, not a proof** (Attack 1). Break Attack 1 first: if the deep-shift-rarity argument doesn't cleanly beat the diagonal for *every* `d ≥ 3`, the whole edifice needs a side condition — and that is the single most valuable thing to learn.
