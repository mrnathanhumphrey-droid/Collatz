# RESULT — C-TILDE-EXTEND: c̃_q=(q−3)/q holds ⟺ ord(2 mod q) large; q=5,7,31 structural, dev ~ 2^{−ord} (2026-07-27)

**Probe:** `probes/probe_ctilde_extend.py` (extends `c_tilde_q17_probe.py`). `c̃_q := lim_k S_k^{(q)}/(q/3)^k`, the
renormalized qx+1 Plancherel-mass constant. Conjecture `c̃_q = (q−3)/q` held at q=11,13,17 but q=5,7 deviated — **open:
finite-k transient or structural?** Now resolved, and the differentiator identified.

**Gate:** float sparse power-iteration vs the exact-rational build at q=17,k=2 — match to **2.2×10⁻¹⁶**.

## PART B — q=5, q=7 deviations are STRUCTURAL, not finite-k (resolves the open question)
Pushed to k=5 (q=5) and k=4 (q=7):
```
 q=5  (q-3)/q=0.400:  c̃(k)= 0.533, 0.492, 0.490, 0.4877, 0.4880   -> ~0.488, dev +0.088 PERSISTENT
 q=7  (q-3)/q=0.571:  c̃(k)= 0.857, 0.782, 0.783, 0.7811          -> ~0.781, dev +0.210 PERSISTENT
```
c̃_5 → 0.488 and c̃_7 → 0.781 are **stable limits ≠ (q−3)/q**. The deviations do **not** shrink with k ⟹ **structural**.

## PART A — large primes, and a surprise: q=31 breaks "large q → (q−3)/q"
```
 q    ord(2 mod q)   (q-3)/q     c̃_q(k=2)    dev
 19       18         0.84211     0.842112    +0.00001   exact
 23       11         0.86957     0.870542    +0.00098   0.1%
 29       28         0.89655     0.896552    +0.00000   exact
 31        5         0.90323     0.961938    +0.05871   *** big deviation ***
 37       36         0.91892     0.918919    +0.00000   exact
```
**q=31 (a large prime) deviates by +0.059** — as much as the small primes. It is Mersenne (`31=2⁵−1`), `ord(2 mod 31)=5`.

## PART C — THE DIFFERENTIATOR: ord(2 mod q), monotone, dev ~ O(2^{−ord})
Deviation from `(q−3)/q` sorted by `ord(2 mod q)` — **perfectly monotone decreasing**:
```
 ord( 3) q= 7:  +0.21029       ord(11) q=23:  +0.00098
 ord( 4) q= 5:  +0.09221       ord(12) q=13:  +0.00056
 ord( 5) q=31:  +0.05871       ord(18) q=19:  +0.00001
 ord( 8) q=17:  +0.00722       ord(28) q=29:  +0.00000
 ord(10) q=11:  +0.00154       ord(36) q=37:  +0.00000
```
- **`c̃_q = (q−3)/q` holds (to <0.1%) iff `ord(2 mod q)` is large** (≥10 ⟹ dev ≤ 0.15%; primitive-root q=19,29,37 ⟹ dev ~1e-5).
  It **fails structurally for small-order primes**: q=7 (ord 3), q=5 (ord 4), q=31 (ord 5) — the Mersenne-like / small-order cases.
- **`dev·(2^{ord}−1)` is O(1)** (1.4–2.3 across all q), so **`dev ~ C·2^{−ord(2 mod q)}`** — the deviation is governed by
  `2^{ord}−1`, the **same `2^M−1` structure** as the current program's denominator theorem and the `4ᵏ−1=2^{2k}−1` Poisson
  unification (`M = ord₂` is exactly the channel index there). A genuine bridge between the qx+1 c̃_q arc and the 7/15 machinery.

## Revises the prior reading
The writeup previously concluded "q=17 non-prim-root still hits (q−3)/q, so non-prim-root status is not the differentiator."
Correct on the negative, but the **positive** differentiator is now pinned: it is the **size of `ord(2 mod q)`**, not
primitive-root status. q=17 (ord 8) is moderate → small dev 0.7%; q=31 (ord 5) is small → big dev 5.9%. `2` being a primitive
root is just the extreme `ord=q−1` (largest possible) ⟹ smallest dev.

## Net
- q=5,7 deviations **structural** (c̃_5→0.488, c̃_7→0.781), not finite-k. `(q−3)/q` is **not universal**.
- **New law:** `c̃_q = (q−3)/q` asymptotically as `ord(2 mod q) → ∞`; deviation `~ 2^{−ord(2 mod q)}`, largest for small-order
  (Mersenne-like) primes. q=31 is a large-prime structural exception (ord 5).
- The `2^{ord}−1` scaling ties this arc to the current program's `2^M−1` / `4ᵏ−1` denominator structure.
- **Open (Wilson's pen):** the exact prefactor `C(q)` in `dev = C(q)·2^{−ord}` (the O(1) residual 1.4–2.3), and an analytic
  derivation of the `q/3` ratio + the `(q−3)/q` limit from Tao's Plancherel framework (still the standing open item).
- **Not at stake:** the current 7/15 / SOLSTICE / GARSIA program (q=3 is the boundary, ord₂ handled by its own machinery).
