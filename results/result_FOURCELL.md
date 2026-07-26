# Probe FOURCELL -- Wilson's direct four-cell decomposition + the rate question -- **ANSWER: outcome 2 (conspiracy), with two clean anchors. NO single cell decays cleanly at the excess rate 0.89 -- the individual cells carry oscillatory sub-structure (two-step rates scatter 0.14-3.2). What IS clean: (a) V1 (the c=1 cell, rotation-1 coherence) is POSITIVE at all 13 levels r=4..16 -- the one sign-definite term, and V1+W2 (rotations 1,2) is the steady-positive 74% majority; (b) the parity ALTERNATION is confined to the ROTATION-0 cells U0,U3 (U0 = clean period-2 -+-+..., exactly Wilson's 3^r mod 4 = 1/3 boundary reading). THE MECHANISM the rates expose: the excess decays SLOWER (0.89) than its own components (~0.78-0.80 envelope) because the rotation-0 group U0+U3 is NEGATIVE early (cancelling V1+W2) and crosses to POSITIVE around r=12 -- the excess's 0.89 is the tail of an EASING CANCELLATION between rotation groups, not a true asymptotic rate. Wilson's instinct confirmed and generalized: not just E[u], but EVERY isolable component decays faster than the excess, so none is the asymptotic carrier -- the carrier is the resolution of the rotation-0-vs-rest cancellation. This is the harder outcome, but the pen has a clean anchor: V1 = p1 E[v|c=1] > 0 at every level.**

**Date:** 2026-07-25. Probe `probes/probe_fourcell.py`, log `logs/fourcell_run.log`. Direct decomposition q_r(1)-1/3 = U0+U3+V1+W2 (each cell its OWN rotation c mod 3: c=0,3->u, c=1->v, c=2->w), U_c = p_c E[own-coh|c]. Table saved scratchpad/fourcell.json.

## The four terms (sum = q-1/3, verified) and their sign patterns
| r | U0 | U3 | V1 | W2 | q-1/3 |
|---|---|---|---|---|---|
| 4 | -4.77e-3 | +1.90e-3 | +4.98e-3 | -2.10e-4 | +1.90e-3 |
| 8 | -1.38e-3 | -1.36e-5 | +2.43e-3 | -2.71e-5 | +1.01e-3 |
| 12 | -4.13e-4 | +4.45e-4 | +5.19e-4 | +8.79e-5 | +6.40e-4 |
| 16 | -6.40e-5 | +1.73e-4 | +2.04e-4 | +1.05e-4 | +4.18e-4 |

Sign across r=4..16: **U0: `-+-+-+-+-+-+-`** (clean period-2), U3: `+-+---+-+-+-+`, **V1: `+++++++++++++`** (all +), W2: `-+----+++++++` (+ from r=10), excess: `+++++++++++++`. Wilson's r=16 split confirmed: U0+U3 = +1.09e-4 (26%), V1+W2 = +3.09e-4 (74%).

## The rate question -- ANSWERED: none clean at 0.89, the excess rate is a cancellation transient
Two-step per-step rates (r=16): U0 0.54, U3 2.16, V1 0.59, W2 0.92, **excess 0.891** (excess is stable ~0.90 across r=12..16: 0.907, 0.916, 0.907, 0.900, 0.891 -- slight DOWNdrift). Individual cells scatter 0.14-3.2 = **oscillatory sub-structure (real, not float noise), no clean per-cell geometric rate.** Envelope over r=4..16: V1+W2 decays ~0.80, U0+U3 ~0.77 (magnitude) -- **both components decay FASTER than the excess (0.89).**
**How the sum decays slower than its parts:** early V1+W2 (+4.77e-3) and U0+U3 (-2.87e-3) partially CANCEL to the excess (+1.90e-3); U0+U3 crosses zero neg->pos around **r=12** (r=11: -2.16e-4 -> r=12: +3.25e-5); late they REINFORCE (+3.09e-4 + 1.09e-4). The excess's 0.89 is the tail of this easing cancellation; asymptotically (once fully reinforcing) it should steepen toward the ~0.78-0.80 component rate -- so **the 0.89 d1-rate quoted throughout may be transitional, not asymptotic** (tentative -- ~1 decade of data; flag, do not pin).

## Structure that IS clean (the pen's anchors)
1. **V1 = p1 E[v|c=1] > 0 at every level** -- the single sign-definite term. The rotation-1 coherence of the c=1 cell. Attacking "E[v|c=1] > 0" is attacking the robust anchor.
2. **The parity/boundary alternation is the ROTATION-0 cells (U0, U3)** -- confirmed Wilson's 3^r mod 4 boundary reading (U0 clean -+-+). Not dynamical; a lattice-alignment effect at the cell boundaries j 3^r/4, and it decays fast (~0.54).
3. **The sign is a rotation-group conspiracy:** V1+W2 (rot 1,2) steady positive; U0+U3 (rot 0) neg->pos crossover ~r=12. The excess is carried by the crossover, not by any single cell.

## Wilson's error-class note (his own flag), confirmed
Both prior predictions failed by averaging over the coordinate holding the information (c independent of coherence). The four-cell form removes that averaging, and the answer is a coordinate-resolved conspiracy -- E[u|c=1] was 20x the other cells because the (c,coherence) coupling IS the structure, and the direct per-cell terms make it explicit.

## Status
FOURCELL: direct decomposition verified (U0+U3+V1+W2 = q-1/3). **Outcome 2 (conspiracy):** no single cell decays cleanly at 0.89; individual cells carry oscillatory sub-structure; the excess's slow-vs-components rate = an EASING CANCELLATION (rotation-0 group U0+U3 neg->pos at r~12) -> the 0.89 rate is likely TRANSITIONAL (flag). Clean anchors: **V1 = p1 E[v|c=1] > 0 all 13 levels** (sign-definite); parity alternation = rotation-0 boundary (Wilson's 3^r mod 4 confirmed). Pen target: E[v|c=1] > 0 (the anchor) + the rotation-0 crossover (the conspiracy). Not at stake: CHANNEL_ID identity, CARRYLEMMA identity, R1-R30, R80-R82, all Thread-3. commit pending.
