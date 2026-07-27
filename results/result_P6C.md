# RESULT — P6C: NO support separation — Wilson's level-1 negative extends to all levels (2026-07-26)

**Probe:** `probe_p6c.py`. Wilson's candidate mechanism for the `3∤m → cross-parity` law (P6B) was a **mod-3 support
separation**: if `supp(ν_e) ⊂ {s≡α mod 3}` and `supp(ν_o) ⊂ {s≡β mod 3}`, α≠β, then same-parity correlations force
`3|m` and cross-parity forces `m≡β−α≢0` — exactly the P6B table. He checked LEVEL 1 by hand and it fails (even-a
images = {1,4,7} mod 9 = base-4 dlogs {0,1,2} = all classes). Question: does a separation appear at a finer
modulus / higher level? Graded the certified, gated `ρ_e, ρ_o` (P6B) by base-4 dlog `s` — Wilson's exact coordinate.

## The check — mod 3 and mod 9, levels 2..5
Max share of the parity-specific L1 mass carried by any single residue class (1.0 = disjoint support = separation):

| level | s mod 3 (ρ_e / ρ_o) | s mod 9 (ρ_e / ρ_o) |
|-------|---------------------|---------------------|
| j=2   | 0.57 / 0.62         | 0.36 / 0.27         |
| j=3   | 0.43 / 0.53         | 0.21 / 0.26         |
| j=4   | 0.46 / 0.39         | 0.24 / 0.20         |
| j=5   | 0.35 / 0.41         | 0.13 / 0.18         |

(gates `ρ_e+ρ_o==ρ_full`: 1e-2 at j=2 tightening to 4e-4 at j=5 — SINGLEREC truncation floor; the qualitative
read is robust across j=3,4,5.)

## Verdict — NEGATIVE, and it TIGHTENS: the support fills in, opposite of separation
No level, no modulus shows concentration. Max class-share sits at **0.13–0.62** and **drifts toward uniform as j
grows** (j=5, s mod 9 → ~1/9 = equidistribution). Wilson's level-1 negative is not a low-level accident — ν_e and ν_o
spread across all base-4-dlog classes at every level, and increasingly so. **The `3∤m → cross` mechanism is NOT a
mod-3 (or mod-9) support separation in the base-4 coordinate.** Per Wilson's own branch ("if a separation shows, the
collapse is a corollary; if not, it's a separate build"), this routes to **P6D — the collapse tested directly**.
Not at stake: P6/P6B findings (d₁=cross-parity, law 3∤m), P1LVL, BRIDGE2, CHANNEL_ID, dichotomy, R1–R30. Cheap.
