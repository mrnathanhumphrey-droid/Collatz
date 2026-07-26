# RESULT — PROBE GRECURSION: the closed-form high-v₂ recursion is EXACT, but |G|≠O(1) and the decay is NOT counting (2026-07-26)

**Probe:** `probe_grecursion.py`. Wilson's high-v₂ recursion in closed form on the lattice ξ=3^i2^j:
`G(i,j) := 2^j π̂_k(3^i2^j)`, `G(i,j) = Σ_{b<j} e(2^b/3^{k−i}) G(i+1,b)` (+ wrap a>j), `|π̂(3^i2^j)| = 2^{−j}|G(i,j)|`.
Run order (Wilson's): #3 pin m(k) & test the rate first; #1 gate the recursion; #2 profile |G|.

## #1 — the closed form is EXACT, wrap negligible (VALIDATED)
Full recursion vs π̂ on the lattice at (i=0,j=m): **rel ~2e-16 all k=3..12.** Truncated (genuine-halving b≥0 only):
rel `~2^{−k}` and shrinking (7.7e-2 @k3 → 1.9e-6 @k12). So the unit-modulus reindexing holds and the wrap terms
(a>j) are down by ~2^{−k} against the main term, exactly as derived. The framework is correct.

## #3 — the decay is NOT a counting fact; |G(0,m)| GROWS (refutes the item-3 hope)
| k | ξ*=2^m | m | dm | S=\|π̂(ξ*)\| | Srate | 2^{−dm} | \|G(0,m)\|=2^m·S | Grate |
|---|--------|---|----|-----------|-------|---------|----------------|-------|
| 3 | 8 | 3 |  | 0.25224 |  |  | 2.018 |  |
| 5 | 32 | 5 | +1 | 0.12927 | 0.730 | 0.500 | 4.137 | 1.461 |
| 7 | 256 | 8 | +2 | 0.07587 | 0.789 | 0.250 | 19.42 | 3.158 |
| 9 | 1024 | 10 | +1 | 0.04803 | 0.789 | 0.500 | 49.18 | 1.577 |
| 10 | 3¹⁰−2¹² | 12 | +2 | 0.03828 | 0.797 | 0.250 | 156.8 | 3.188 |
| 12 | 515057 | 14 | +1 | 0.02646 | 0.828 | 0.500 | 433.5 | 1.657 |

`S = 2^{−m}|G|` is an **identity**, so `Srate = 2^{−dm}·Grate`. The observed 0.80/step is therefore **half counting,
half growth of G**: dm=+1 gives `2^{−1}·1.6 ≈ 0.80`, not the `2^{−1}=0.50` that a bounded G would give. **Wilson's
item-3 collapse ("decay = where the argmax sits") is refuted** — the counting factor `2^{−m}` is real but exactly half
the story.

- **m(k) = 3,4,5,6,8,9,10,12,13,14** for k=3..12. m−k = 0,0,0,0,1,1,1,2,2,2 — rises by 1 every ~3 levels ⟹
  **m ~ 4k/3** asymptotically (dm/dk → 4/3), i.e. the argmax climbs *faster* than the grid.
- **|G(0,m)| is NOT O(1)** — it grows ~geometrically, Grate → ~1.6–2.0/step (mean ~2.0 over the 1,1,2 dm-pattern).
  Solving `2^{−m}|G| = 0.80^k` with m~4k/3 gives **|G(0,m)| ~ 2.0^k** — a genuinely growing object, not a constant.

## #2 — |G(0,j)| profile: monotone growth, no plateau
`|G(0,j)| = 2^j|π̂(2^j)|` climbs monotonically across j (k=12: 0.00→18837 over j=1..26). For j>m the 2^j factor
simply dominates the ~typical π̂, so |G| grows trivially past the argmax; the meaningful object is |G(0,m)| at the
argmax, and that grows too. There is no bounded regime.

## Verdict — the content stayed in the ⟨2⟩-orbit sum; it did not collapse into counting
The closed form is exact and the decay factorizes cleanly as `sup|π̂| = 2^{−m(k)} · |G(0,m)|`, but **both factors are
nontrivial**: the counting factor `2^{−m}` (m~4k/3) and the **growth of |G(0,m)| ~ 2^k**. The true rate 0.80/step is
the ratio of these two competing geometric growths. So:
- The Heilbronn/BGK ⟨2⟩-orbit-sum shelf is still the right target — but for the **growth rate of |G(0,m)|**
  (`|Σ_b e(2^b/3^{k−i}) G(i+1,b)|` compounded down the tree), NOT for a boundedness statement. Wilson's hope that the
  orbit sum is O(1) is refuted: it grows ~2^k, and the ⟨2⟩-orbit phases supply only *partial* cancellation (enough to
  hold Grate below 2, not enough for O(1)).
- The sharpened open object: `|G(0,m)| ~ C·λ^k` with **λ ≈ 2.0** measured — pin λ and the constant C, and then
  `sup|π̂| ~ C·2^{−m(k)}·λ^k` with m(k) the argmax-location law. The decay rate 0.80 = `λ·2^{−dm/dk} = 2.0·2^{−4/3}`.
- **NOT a counting fact.** The cancellation content is real and lives in G's growth exponent λ. Pen territory: bound λ
  (or the per-level orbit-sum growth) via the subgroup-cancellation machinery, now that the regime is unit-coefficient
  orbit sums mod 3^{k−i}.

**Not at stake:** SINGLEREC, BRIDGE2, CONTRACTION (this refines its "arg at feed points" seam into the explicit G-tree),
MAXMODE2, MEAN1, HIERARCHY, CHANNEL_ID, R1–R30. Cheap (2.3s).
