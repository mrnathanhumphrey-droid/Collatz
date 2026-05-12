# PATH2 Pushback — Check 4: Translation-invariant ingredients audit

**Adversarial frame:** The "tradition" (amplification + Type I/II + L-function moment) relies on:
- (a) smooth amplitudes
- (b) square-free moduli
- (c) Cauchy-Schwarz-affordable halving

The Path 2 claim is "right side of the structural divide" — explicitly NOT using these. Does Phase 3's argument silently use any of (a), (b), or (c)?

## Disposition

> **NO_TRADITION_INGREDIENTS** at r=2 and r=3. The argument uses Cochrane factorization, saddle exactness, Plancherel orthogonality (exact identity, NOT Cauchy halving), triangle inequality, and the cosecant grid identity. None of these are tradition ingredients.

**Caveat at r ≥ 4:** the Hensel-triangle fallback used at r ≥ 4 uses an extra triangle inequality on the deviation D(a), which produces the loose log N factor (Check 3 finding). This is NOT a tradition ingredient (it's an elementary triangle, not Cauchy), but the looseness of the resulting bound is one reason the "tight √N" picture at r ≥ 4 remains open.

## Detailed audit of each ingredient

### (a) Smooth amplitudes — NOT USED

The Path 2 argument operates on the **sharp Dirichlet kernel** `1̂(p·a) = Σ_{u=0}^{N-1} e_q(p·a·u)` throughout. At no step is `1̂` replaced by a smoothed analog (e.g., convolution with a Gaussian, Beurling-Selberg majorant, or partial-summation-induced smooth weight).

The "1/sin grid identity" / cosecant grid sum operates on `|1̂(p·a)| = sin(π·a/p)/sin(π·a/p²)` — the exact closed-form magnitude of the sharp kernel. The asymptotic `Σ csc(πα/p) ~ (2p/π) log p` is a discrete-lattice cosecant sum, not a smooth integral.

The Inner-Plancherel collapse exploits the **discrete orthogonality** `Σ_{c_2 ∈ Z/p} e_p(c_2·k) = p·𝟙(k≡0)` — also sharp, not smoothed.

**Verdict (a): No smooth amplitudes used.** ✓

### (b) Square-free moduli — NOT USED

The modulus throughout is `q = p^{r+1}` — a prime power, NOT square-free. The argument is specifically for prime-power q.

At no step is the bound improved by factoring q into coprime square-free components. The Chinese Remainder Theorem decomposition into coprime moduli (e.g., when q = q_1 · q_2 with gcd(q_1, q_2) = 1) is NOT invoked. There is no "square-free-essential" structure used.

The Cochrane factorization itself operates at prime-power level — it's a structure that **doesn't exist** at square-free moduli (the principal-unit cyclic subgroup of order p^r in (Z/p^{r+1})* requires the prime-power setting). So the entire argument is intrinsically prime-power-native.

**Verdict (b): No square-free factoring used.** ✓

### (c) Cauchy-Schwarz halving — NOT USED IN THE ACCEPTED CHAIN

The PATH2_BILINEAR.md doc explored Cauchy-Schwarz approaches (Attempts B, D) but **abandoned them** because they gave bounds AT or WORSE THAN trivial. From line 414-417: "**This Cauchy-Plancherel approach can't beat trivial.** The phase cancellation has to come from a different mechanism."

The accepted chain (Attempt G+, line 491+) uses **triangle inequality** on the outer sum:
`|T_p| ≤ Σ_{s*} |Inner(s*)|`

and bounds each `|Inner(s*)|` **pointwise via an exact identity** (Inner-Plancherel collapse on c_2 via Z/p orthogonality), THEN sums the magnitudes with the cosecant grid identity.

The Z/p orthogonality used in the Inner-Plancherel collapse is an **EXACT identity**:
`Σ_{c_2 ∈ Z/p} e_p(c_2·k) = p · 𝟙(k ≡ 0 mod p)`

NOT Cauchy-Schwarz. There's no `|sum|² ≤ N · Σ|.|²` style halving. The collapse is a discrete delta function, exploited because the inner phase is LINEAR in c_2 (which makes the c_2-sum a trivial geometric/orthogonality sum).

**Verdict (c): No Cauchy-Schwarz halving used in the accepted chain.** ✓

## Pre-registered "tradition" ingredient — is the linear-in-c_2 collapse a hidden tradition use?

**Counter-argument to consider:** the "linear-in-c_2 collapse via orthogonality" looks structurally similar to the moment-method-with-orthogonality used in some tradition arguments (e.g., Burgess uses orthogonality of multiplicative characters). Is it tradition under a different name?

**No, and here's why:**

1. The orthogonality used is **additive character orthogonality on Z/p** (`Σ e_p(c·k) = p·δ_{k≡0}`), not multiplicative character orthogonality. Multiplicative character orthogonality on (Z/q)* is what Burgess/Heath-Brown leverage. The additive version is much more elementary (it's just discrete Fourier inversion).

2. The collapse is "free" (no inequality involved) BECAUSE the inner phase is linear in c_2 — a structural fact about the saddle phase. If the phase were quadratic-in-c_2 (e.g., at r ≥ 4 with Hensel correction), the collapse would NOT be exact, and one would need bilinear/Cauchy methods (then yes, tradition).

3. The argument's "save" comes from **the structure of the saddle phase**, not from amplification or moment method. The saddle phase is structurally linear in c_2 at r=3 (proved by direct computation of P_a(s*) mod p^4). This is a feature of the Cochrane Prop 4 + saddle setup, not a Cauchy-Schwarz / Fourier-amplification trick.

**Conclusion:** the linear-in-c_2 collapse is NOT a hidden tradition ingredient. It is an exact identity following from the structure of P_a(s*).

## Where the line is drawn

The Path 2 argument is "right side of structural divide" because it uses:

✓ Cochrane truncated p-adic log (structural identity, p-blind)
✓ Plancherel saturation / Gauss-sum magnitude (Parseval identity, exact)
✓ Saddle-exactness at r=3 (polynomial identity, verified by Phase 2)
✓ Linear-in-c_2 phase structure (algebraic identity from saddle computation)
✓ Z/p additive character orthogonality (trivial discrete identity)
✓ Triangle inequality on outer sum (elementary)
✓ Cosecant grid sum asymptotic (classical Hardy-Littlewood-level, NOT moment-method)

It does NOT use:

✗ Smooth-amplitude majorants (Beurling-Selberg, partial summation)
✗ Square-free factorization / CRT decomposition
✗ Multiplicative character orthogonality (Burgess / Heath-Brown / cubic-character moments)
✗ Cauchy-Schwarz halving on the bilinear (Attempts B, D abandoned)
✗ Type I / Type II sum decomposition
✗ L-function moment estimates
✗ Vinogradov mean-value theorem

## The r ≥ 4 fallback uses triangle, not Cauchy

At r ≥ 4 the argument falls back to:
`|S_partial| ≤ |S_partial(lead)| + Σ_a |1̂(p·a)| · max|D(a)|`

This is triangle inequality on `S_partial = S_lead + Σ 1̂ · D`. NOT Cauchy-Schwarz. The looseness comes from the Hensel-triangle wasting structural cancellation in `Σ 1̂ · D` — but the looseness is via triangle, not Cauchy halving.

So the r ≥ 4 bound, while loose (Check 3), still doesn't use tradition ingredients.

## Verdict

> **NO_TRADITION_INGREDIENTS** at r ≤ 3 (the saddle-exact regime). The argument operates entirely via Cochrane / saddle / orthogonality / triangle / cosecant — none of which are smooth amplification, square-free essential, or Cauchy halving.

**At r ≥ 4:** still no tradition ingredients, but the Hensel-triangle fallback gives a loose bound (Check 3). Tightening this at r ≥ 4 to remove the log factor would require either explicit Hensel-lifted closed form (open) or a structural inequality on `Σ 1̂·D` that respects D's class-mean-zero structure — neither uses tradition.

## Pre-registration trigger check

Per pre-reg: "If Checks 1, 4, OR 6 fails → full walk-back." 

Check 4 passes: NO_TRADITION_INGREDIENTS. **No walk-back trigger from Check 4.**

## Implication for Tao email framing

The "right side of structural divide" framing is preserved. The Path 2 argument can be honestly characterized as:

> "Phase 3 derives `|S_partial| ≤ 2√N` at r ≤ 3 family-level via:
> (i) Cochrane Prop 4 + truncated p-adic log → explicit polynomial form of P_a(s*) mod p^4 at r=3;
> (ii) Inner-Plancherel collapse via Z/p additive character orthogonality (exact identity, exploits linear-in-c_2 phase structure);
> (iii) cosecant grid sum identity on the resulting Dirichlet-kernel magnitudes.
>
> The chain does NOT invoke smooth amplitudes, square-free factorization, Cauchy-Schwarz halving, or L-function moment estimates. The bound is therefore on the 'right side' of the Pascadi-diagnosed structural divide — not relying on the tradition Type I/II ingredients."

This framing is accurate and defensible.
