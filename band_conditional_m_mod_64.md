# Band-conditional m ≡ 21/53 mod 64 split (Result 54)

**Status.** Outcome (c) with clean structural picture. The simple closed form
E_V = 5 + 2·p_21(B) fails by up to 0.13 in extreme bands. Three layers of
band-conditional structure: (1) p_21(B) split varies monotonically;
(2) within m ≡ 21 mod 64, the V distribution itself shifts with band
(deeper-mod structure); (3) gap-conditional p_21 drops dramatically for
long gaps (confirms the gap-mixing mechanism).

The open piece is at minimum 2D: (band, mod-2^k for k ≥ 7) — not just 1D
p_21(B). But the structural mechanism is identified: gap length determines
m's mod-2^k bit pattern at next visit, which determines V.

## Step 1-3: Direct p_21(B) measurement vs predicted E_V

Walking 500K orbits at N=2³⁶, 2.40M r=21 visits.

| band | n | p_21 | p_53 | pred E_V | emp E_V | gap |
|------|---:|---:|---:|---:|---:|---:|
| 0–25 | 560,231 | 0.5366 | 0.4634 | 6.0731 | **6.1990** | **−0.126** |
| 25–50 | 574,525 | 0.4803 | 0.5197 | 5.9605 | 5.9632 | −0.003 |
| 50–75 | 603,876 | 0.4428 | 0.5572 | 5.8855 | 5.8410 | +0.045 |
| 75–95 | 515,934 | 0.4150 | 0.5850 | 5.8301 | 5.7439 | +0.086 |
| 95–100 | 144,429 | 0.3937 | 0.6063 | 5.7874 | **5.6707** | **+0.117** |

The simple closed form predicts E_V from p_21 alone, assuming
**E[V | m ≡ 21 mod 64] = 7** (pure shifted Geom(1/2) starting at V=6).

Gap up to 0.126 — significantly larger than bootstrap noise. **The simple
prediction fails.**

## Step 4: Why — within m ≡ 21 mod 64, V distribution is band-dependent

| band | P(V=6) | P(V=7) | P(V=8) | P(V=9) | P(V≥10) | ⟨V⟩ − 6 |
|------|---:|---:|---:|---:|---:|---:|
| pred uniform | 0.5000 | 0.2500 | 0.1250 | 0.0625 | 0.0625 | 1.000 |
| 0–25 | 0.4494 | 0.2351 | 0.1424 | 0.0667 | **0.1065** | **1.235** |
| 25–50 | 0.4995 | 0.2498 | 0.1201 | 0.0579 | 0.0727 | 1.006 |
| 50–75 | 0.5235 | 0.2567 | 0.1095 | 0.0547 | 0.0557 | 0.900 |
| 75–95 | 0.5555 | 0.2543 | 0.0993 | 0.0487 | 0.0422 | 0.792 |
| 95–100 | **0.5827** | 0.2544 | 0.0898 | 0.0407 | 0.0323 | **0.704** |

**Within m ≡ 21 mod 64, the geometric tail differs by band**:
- Low-σ (0–25): heavier upper tail, P(V≥10) = 0.107 vs uniform 0.0625
- High-σ (95–100): concentrated at V=6, P(V=6) = 0.583 vs uniform 0.5

This is the **j-distribution at r=21 mod 64 visits being band-dependent**.
Equivalently: m mod 2^k for k > 6 distribution depends on band — not just
m mod 64.

The deeper-mod visit measure carries band-conditional structure beyond the
mod-64 split.

## Step 5: p_21(B) curve and uniform-cylinder deviation

Uniform-cylinder prediction: p_21 = 0.5 in every band.
Empirical p_21 ranges 0.394 to 0.537.

| band | p_21 | deviation |
|------|---:|---:|
| 0–25 | 0.537 | **+0.037** |
| 25–50 | 0.480 | −0.020 |
| 50–75 | 0.443 | −0.057 |
| 75–95 | 0.415 | −0.085 |
| 95–100 | 0.394 | **−0.106** |

Monotone decreasing in band index (matches direction predicted by
gap-mixing mechanism). The signature is roughly linear: p_21(B) ≈ 0.50 −
0.18·B_index where B_index ∈ {−1, −0.5, 0, +0.5, +1} would map to bands.

But this is only PART of the open piece — Layer 2 (within-class V) is the
other part.

## Step 6: Gap-conditional p_21 — dramatic long-gap effect

Pooled across all bands, p_21 stratified by gap quintile:

| gap range | n | p_21 |
|---|---:|---:|
| [1, 4] | 462,789 | 0.491 |
| [4, 8] | 428,010 | 0.523 |
| [8, 15] | 453,699 | 0.513 |
| [15, 28] | 401,469 | 0.482 |
| **[28, 207]** | 393,309 | **0.274** |

For long gaps (≥ 28 steps), p_21 drops to 0.274 — orbits arrive at r=53 mod 64
dramatically more often. **Confirms the gap-mixing mechanism**: long gap = more
deterministic-residue traversals = m's bits more "mixed" = arrival lands at
the more probable mod-2^k residue (which empirically is r=53 mod 64).

In band 50–75 specifically, p_21 drops from 0.50 (typical gap) to 0.22
(gap ≥ 32). Strong gap-conditional dependence within band.

## Verdict — outcome (c) with structural picture

The mechanism is partially correct:
- p_21(B) varies as predicted (monotone decreasing with band)
- Long gaps drive p_21 down (gap-mixing mechanism confirmed)
- BUT the simple closed form E_V = 5 + 2·p_21 fails by 0.13 in extreme bands

The reason: within m ≡ 21 mod 64, the V distribution is itself band-dependent
because the **j-distribution (m mod 2^k for k > 6) is non-uniform per band**.

## What's the open piece, exactly

Not 1D p_21(B). At minimum **2D**: the band-conditional joint distribution
of (m mod 64, m mod 2^k) at r=21 visits, for k > 6.

Equivalently, given Result 50/52's closed-form arithmetic V = 5 + v_2(2+3h)
where h = (m−21)/32:

> **The open piece is the band-conditional distribution P(h | r=21 visit, band)
> over h ∈ {0, 1, 2, ..., 2^(k−5) − 1} for arbitrary k.**

Under uniform h: V follows shifted Geom(1/2) starting at V=5.
Under band-conditional h-distribution: V deviates from uniform Geom in
band-specific ways.

The j-distribution at r=21 mod 64 (the V≥6 sub-cylinder) is non-uniform per
band; the geometric tail depths reflect this.

## For v3.6 / Chang correspondence

The Lagarias-class open piece is:

> **The band-conditional measure on the 2-adic cylinder {m ≡ 21 mod 32}**:
> P(m mod 2^k | r=21 visit, σ-band) as k → ∞, with non-trivial structure at
> all k ≥ 6. Not 1D-reducible to p_21(B); minimum 2D over band × deeper mod.

The mechanism is identified (gap-conditional bit-mixing on the deterministic
31-residue spectrum). The 1D p_21(B) curve is the marginal projection at
mod 64; deeper-mod structure persists at all higher k.

Given Result 47's mixed (b)+(d) verdict (renewal at gap level + G→V coupling),
this Result confirms the G→V coupling operates by gap-conditional bit-mixing.
Long gaps mix the bits more thoroughly, producing different mod-2^k visit
measures than short gaps.

## Reduction for the unified framework

The unified framework requires capturing:
1. **Renewal at visit level** (Result 49, confirmed)
2. **Closed-form arithmetic V = f(m mod 2^k)** (Result 50/52, confirmed)
3. **Band-conditional visit measure on cylinder** (this Result, open at deeper mods)

The deterministic 31-residue Markov chain's mixing rate determines (3). If
that mixing has closed form (a property of the specific Collatz residue
dynamics on Z/2^k), the entire trajectory measure closes.

**The remaining mathematical question is the asymptotic mixing of m mod 2^k
on the deterministic spectrum, conditional on σ-band.**

## Files

- `experiments/78_band_conditional_split.py`
- `experiments_output/78_band_conditional_split_log.txt`
- `experiments_output/78_band_conditional_split.csv`

Compute: 1.2s (500K orbits at N=2³⁶, 2.40M r=21 visits).
