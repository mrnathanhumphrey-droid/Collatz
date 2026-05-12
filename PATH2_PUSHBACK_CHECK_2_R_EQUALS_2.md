# PATH2 Pushback — Check 2: r=2 derivation mechanism

**Adversarial frame:** Is the r=2 bound `|S_partial| ≤ (1 + log(p)/p) · √N` a **special case** of the r=3 Inner-Plancherel argument, or a **separate construction** with tradition ingredients that would mean the "right side of the structural divide" framing applies only at r=3 and above?

## Disposition

> **r=2_SEPARATE_CONSTRUCTION** with ingredients: triangle inequality on the phase + harmonic-decay bound on `|1̂(p·a)|`. NO tradition ingredients used (no Cauchy-Schwarz, no smooth amplification, no square-free factorization). Right side of the structural divide preserved.

## Step-by-step trace of the r=2 argument

At r=2: q=p^3, period=p², N=p, |supp|=p (singleton classes when partitioned by s*).

### Step A — Setup (shared with r=3)

Cochrane factorization (T78.4_p, rigorous, p-blind) gives `F̂_p(p·a) = p · e_q(c) · G_p(a)` and `|G_p(a)| = √q` on support.

### Step B — Phase formula at r=2

From the truncated p-adic log expansion at r=2 (J_p = 2 for p ≥ 5, J_p = 3 for p = 3):

`P_a(s*) ≡ −p²·s*²/2 mod p³` (for p ≥ 5, J_p = 2)
`P_a(s*) ≡ −p²·s*²/2 + (3-specific correction) mod p³` (for p = 3, J_p = 3)

In both cases:
`e_q(P_a(s*)) = e_p(−s*²/2) · ε_p`

where `ε_p` is an a-independent root-of-unity Gaussian factor. **Phase 2 empirical verification confirmed `r2_uniform_factor = True` at all 4 tested cells** (p ∈ {3,5,7,11}, r=2; PATH2_FAMILY_EXTENSION_VERIFICATION.csv):

| p | r=2 ε_p | ε_p magnitude |
|---|---|---|
| 3 | 0.866 + 0.500i = e^{iπ/6} | 1.000 |
| 5 | -1.000 + 0i = e^{iπ} | 1.000 |
| 7 | 0 + 1.000i = e^{iπ/2} | 1.000 |
| 11 | 0 - 1.000i = e^{-iπ/2} | 1.000 |

The ε_p factor is a-independent (uniform across all a in support), so:

`G_p(a) = √q · ε_p · e_p(−s*(C_a)²/2)`

### Step C — The bilinear sum at r=2

`T_p := Σ_{a ∈ supp} 1̂(p·a) · e_q(P_a(s*(C_a)))`
     `= ε_p · Σ_{a ∈ supp} 1̂(p·a) · e_p(−s*(C_a)²/2)`

|T_p| ≤ |ε_p| · |Σ_a 1̂(p·a) · e_p(−s*²/2)|

At r=2, the bijection a ↔ s* is essentially `s* = (C_a − 1)/p mod p` with `c_2 = 0` (no c_2 dimension at r=2 — |supp| = p exhausted by the s* index). So the inner Plancherel structure of r=3 (collapse on c_2 dimension) is **trivial / absent** at r=2.

### Step D — Direct triangle + |1̂| decay bound

At r=2, no phase cancellation across s* is exploited. The bound proceeds by:

`|T_p| ≤ Σ_{a ∈ supp} |1̂(p·a)|`

This is just triangle on the phase. Now bound `Σ |1̂(p·a)|`:

For a ∈ {1, 1+p, 1+2p, ..., 1+(p-1)p}, the kernel `1̂(p·a) = Σ_{u=0}^{p-1} e_q(p·a·u) = Σ_u e_{p²}(a·u)`, which is a length-p Dirichlet kernel mod p². Magnitudes:

`|1̂(p·a)| = sin(π·a·p / p²) / sin(π·a / p²) = sin(πa/p) / sin(πa/p²)`

For a = 1+pα with α ∈ {0,...,p-1}:
- α=0: |1̂| = sin(π/p)/sin(π/p²) ≈ p
- α ≥ 1: |1̂| ≈ sin(π/p)/sin(πα/p) (the π/p² correction is sub-leading)

The α=0 term contributes ~p; the α ≥ 1 sum is the cosecant grid sum ≈ 2 log(p) (Check 1 §F asymptotic).

Hence `Σ |1̂(p·a)| ≤ p + 2·log(p) + O(1)`, and:

`|T_p|(r=2) ≤ p + 2·log(p) ≤ 2N for p ≥ 3` (since 2·log(p) ≤ p for p ≥ 3).

Converting to |S_partial|: at r=2, `|S_partial| = (p/√q) · |T_p| = (p/(p·√p)) · |T_p| = |T_p|/√p = |T_p|/√N`. So:

`|S_partial|(r=2) ≤ (p + 2·log(p))/√N · (1/p) · √p ... `

Wait, let me redo more carefully. With q = p^3, N = p, √q = p^{3/2}:

`|K_p|(r=2) = (p/√q) · |T_p| = (p / p^{3/2}) · |T_p| = |T_p|/√p = |T_p|/√N`

So `|K_p|(r=2) ≤ (p + 2 log p)/√N = √p · (1 + 2 log p / p) ≤ √N · (1 + 2 log p / p)`.

For p ≥ 3, 2 log p / p ≤ 0.73, so `|K_p|(r=2) ≤ 1.73·√N < 2√N`.

(Note: PATH2_DISPOSITION.md line 15 reports "(1 + log p/p)·√N" with constant ≤ 2 — slight off-by-factor-of-2 in the log term, same as the Check 1 doc arithmetic discrepancy. Doesn't affect the bound shape or final constant ≤ 2.)

## Ingredients audit at r=2

1. **Cochrane factorization (T78.4_p):** rigorous, p-blind. ✓ NOT tradition.
2. **Saddle exactness at r=2 (T78.6_p phase):** verified empirically (uniform Gaussian factor ε_p, r2_uniform_factor = True in Phase 2 CSV). ✓ NOT tradition.
3. **Triangle inequality on phase:** elementary. ✓ NOT tradition.
4. **Harmonic-decay bound on |1̂(p·a)|:** classical Dirichlet-kernel asymptotic / cosecant sum. ✓ NOT tradition.

**No Cauchy-Schwarz, no smooth amplification, no square-free factorization.** Right side of structural divide preserved at r=2.

## Comparison with r=3 mechanism

| | r=2 | r=3 |
|---|---|---|
| Phase structure | e_p(−s*²/2)·ε_p (quadratic in s*, no c_2) | e_{p²}(−s*²/2)·e_p(s*³/6 − c_2·s*) |
| Inner-Plancherel collapse | not needed (1 element per s* class) | required (p elements per s* class, c_2 dimension) |
| Bound mechanism | triangle + harmonic decay on |1̂| | triangle + Inner-Plancherel + cosecant grid identity |
| Final shape | √N · (1 + 2 log p/p) | √N · (1 + 2 log p/p) (Check 1) |
| Final constant | ≤ 2 | ≤ 2 |

Both mechanisms land at the same bound `≤ 2√N`, but via DIFFERENT paths:
- r=2: direct |1̂| decay (harmonic/cosecant bound on |1̂(p·a)|).
- r=3: Inner-Plancherel reduces to a cosecant grid sum on |D_p(a_0(s*), p²)|.

These coincide in shape because at r=2 the Inner-Plancherel structure is degenerate (Inner(s*) reduces to a single term, namely 1̂(p·a(s*,0))), so r=2 IS a degenerate special case of r=3 — but the bound mechanism is simpler.

## Adversarial frame revisited

**Q:** Is r=2 a special case of r=3 (then no separate construction issue), or is it separate-with-tradition-ingredients (then "structural divide" framing applies only at r=3)?

**A:** r=2 is a **degenerate special case of r=3 with simpler mechanics**: the Inner-Plancherel reduction at r=2 has trivial inner sum (one c_2 value per s*-class since the c_2 dimension is absent), so the bound reduces to `|T_p| ≤ Σ |1̂(p·a)|`, bounded by triangle + cosecant-grid.

This is **not tradition-ingredient construction** — no Cauchy-Schwarz, no smooth amp, no square-free factoring. So the "right side of the structural divide" framing applies at r=2 AND r=3, with r=2 being the simpler case.

## Verdict

> **r=2_SEPARATE_CONSTRUCTION** in mechanism (triangle + |1̂| harmonic decay, distinct from r=3's Inner-Plancherel), but **with no tradition ingredients used**. Right-side-of-divide framing preserved at r=2.

**Pre-registration partial walk-back trigger?**

Per the locked rules: "If Check 2 = SEPARATE_CONSTRUCTION_WITH_INGREDIENTS → partial walk-back on r=2 framing". 

The condition is "with tradition ingredients". The mechanism is separate but does NOT use tradition ingredients (no Cauchy, no smooth amp, no square-free). So the trigger condition is NOT met. **No walk-back on r=2 framing.**

## Implication for Tao email framing

The r=2 result is the simpler case — same bound `≤ 2√N`, derived via direct triangle + harmonic decay (no Inner-Plancherel needed). The "right side of structural divide" claim is preserved at r=2 (no tradition ingredients).

For Tao: it's accurate to say "the bound holds at both r=2 and r=3 family-level, with r=2 reducing to elementary Dirichlet-kernel decay and r=3 requiring an Inner-Plancherel + cosecant-grid argument; neither uses tradition Type I/II ingredients."
