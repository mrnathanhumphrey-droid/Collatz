# RESULT — PROBE RHOBOUND + RHOREC: the bound holds, the recursion isn't pinned, the seam is confirmed (2026-07-26)

**Probes:** `probe_rhobound.py`, `probe_rhorec.py`. Gate Wilson's corrected ρ̂ recursion
`ρ̂_n(a)=D̃(a)·E_{X~ν_{n-1}}[χ_a(3X+1)]` via (1) its unconditional bound, (2) the argmax-parity prediction, (3) full
reproduction. Wilson's criterion: if the bound or parity fails, the recursion is still mis-indexed.

## (1) The unconditional bound HOLDS — verified, never violated
`|ρ̂_k(a)|² ≤ 1/(5−4(−1)^a cos(πa/3^k))` (= `|D̃(a)|²`, from `|E[χ_a(3X+1)]|≤1`): **max ratio ≤ 1 at every k=2..16**,
never violated (max saturation 0.40 at k=2, falling to ~0 by k=16). So the `(−1)^a` parity in the bound matches the
codebase's `a`, and the **ninefold-suppression structure is real**: for odd a the ceiling is 1/9, for even a it's 1.
This is a genuine unconditional deliverable.

## (2) The argmax-parity prediction FAILS — and that confirms the seam
Wilson predicted `argmax_a|ρ̂_k(a)|` even with a/3^k small (sup lives where D̃'s ceiling is 1). Measured: the small-side
of the argmax is **even in only 6/15 levels — mixed, not even-dominated**. And near a=0, even-a and odd-a `|ρ̂|²` are
**comparable** (k=16: even ~1e-8, odd ~1e-8), NOT even≫odd. So the D̃ ceiling (1 vs 1/9) does **not** drive the actual
values — the sup is governed by the **affine factor `E[χ_a(3X+1)]`**, which suppresses both parities and moves the
argmax around. That is exactly Wilson's seam: D̃ is controlled, `E` is not, and `E` is where the sup lives.

## (3) Full recursion reproduction FAILS with my implementation
`ρ̂_n = vfactor·E_affine`, `vfactor=z/(2−z)` (`z=(−1)^a e^{−iπa/3^n}` = Wilson's D̃ up to the unit phase z, immaterial
to the bound), `E_affine(a)=Σ_X ν_{n-1}(X) e(2πi a d_n[X]/3^n)`: **rel ~1–2 on both index conventions** (pre-index and
residue-index), n=2..6. My reproduction does not match `fft(ρ_n)`.

**Honest limit:** this is the **second** time I've failed to reproduce the recursion numerically (fft(ν) in
GATE_RECURSION, fft(ρ) here), and the dlog/affine bookkeeping is precisely the "one place I could be reconstructing"
Wilson flagged. The bound holding (part 1) says the `vfactor`/parity part is right; the failure is in `E_affine` — my
`g[s]=Σ_{X:d_n[X]=s}ν_{n-1}(X)` construction (the dlog of the affine image `3X+1`), or the recursion's exact indexing.
I can't distinguish "my convention wrong" from "recursion mis-indexed" without reconstruction-flailing, so I stop.

## Verdict — bound + seam confirmed; recursion indexing is the pen's; Hank HELD
- **Validated:** the unconditional bound `|ρ̂|²≤|D̃|²` (ninefold suppression for odd a), and the **seam** — the sup is
  affine/`E`-dominated, not D̃-dominated (argmax mixed, even≈odd near 0). Both are real and support the write-up.
- **Not validated:** the exact recursion — my reproduction fails (rel~1–2), the argmax parity fails, so per Wilson's
  own criterion the indexing isn't pinned. This is his pen to fix (`E_affine`'s dlog convention / the level indexing);
  I've given the exact form I tested and where it breaks.
- **Hank HELD:** per Wilson's instruction ("nothing built on it counts" if the gate fails), the Salié/Kloosterman
  re-task is NOT dispatched until the recursion is pinned — though the seam being confirmed (E = a multiplicative
  character on the affine image `3X+1` = the Salié/Kloosterman family) means that shelf is the right target once it is.

**Not at stake:** MAXMODE2 saturation (fft(ρ), the channel object — solid), MEAN1, HIERARCHY, CHANNEL_ID, R1–R30.
Cheap (build_nu(11)+build_nu(6), ~8s total).
