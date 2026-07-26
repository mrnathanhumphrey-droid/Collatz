# Probe CARRYLEMMA -- Wilson's carry identity + attenuation mechanism -- **THE LEMMA IS EXACT (verified as Fractions at r=4,5: q-1/3 == E_mu<dpi_x, T_c dpi_4x> EXACTLY; float rel <= 2e-13 to r=16; u+v+w=0 at machine zero; the survival conditioning IS fully absorbed into mu = nu(x)nu(4x) -- the bookkeeping stands). BOTH PRE-REGISTERED PREDICTIONS FAIL, and the failure is the delta-2 answer: (P1) E[u] is NOT ~4(q-1/3) -- it ALTERNATES SIGN with r (-0.043, +0.023, -0.025, ..., -0.0012 at r=16), is 3-25x LARGER than the excess, and at r=16 has the WRONG sign (ratio -0.74). The independence assumption fails completely: the (c,u) coupling is not a correction, it IS the mechanism. (P2) cells are SKEWED (P1=0.262 > P2=0.243 at r=16), the v-term is LIVE and carries ~58% of the excess (T2=+2.43e-4 vs T1=+1.75e-4). THE STRIKING FACT: the raw coherence u alternates sign with r, yet the carry-rotated excess is steadily POSITIVE at every r=4..16 -- the rotation is not attenuating a fixed positive coherence, it is PHASE-ALIGNING an alternating one. Wilson's coordinate-mismatch caution is confirmed as the operative difficulty: the (c,u) coupling is strong and structured -- the x2/x3 tension in its fourth hat, now measured.**

**Date:** 2026-07-25. Probe `probes/probe_carrylemma.py`, log `logs/carrylemma_run.log`. Exact Fractions r=4,5; float r=4..16 (fresh build_nu(16)).

## Gates -- the lemma survives everything
- **Identity EXACT:** q_r(1) - 1/3 = E_mu<dpi_x, T_{c(x)} dpi_{4x}> -- equal AS FRACTIONS at r=4,5; float rel 5e-15..2e-13 at r=4..16. mu(x) prop nu_low(x) nu_low(4x mod 3^r); pi_x(d) = nu_hi(x+d 3^r)/nu_low(x); c = floor(4x/3^r); lift event D' = c+D (mod 3) [4 = 1 mod 3].
- **u+v+w = 0:** machine zero (1.4e-17) at every x, every r (dpi mean-zero on Z/3).
- **Two-term decomposition:** T1 + T2 == q-1/3 exact (Fractions) and to 2e-13 (float), where T1 = E[(1_{c in {0,3}} - 1_{c=2}) u], T2 = E[(1_{c=1} - 1_{c=2}) v].
- **Survival conditioning fully absorbed** into the mu-weighting -- no Doob/QSD needed. Deltas (1),(3) closed as claimed.

## P1 FAILED -- and the failure is the answer to delta (2)
| r | q-1/3 | E[u] | E[u]/4(q-1/3) |
|---|---|---|---|
| 4 | +1.90e-3 | -4.33e-2 | -5.69 |
| 5 | +1.76e-3 | +2.34e-2 | +3.32 |
| 8 | +1.01e-3 | -1.54e-2 | -3.81 |
| 12 | +6.40e-4 | -4.26e-3 | -1.66 |
| 16 | +4.18e-4 | -1.24e-3 | -0.74 |

**E[u] ALTERNATES SIGN with r** (period-2, decaying ~0.5-0.6/level -- same alternating family as AC-E's fiber antisymmetry; flagged as an observation, not an identification) while the excess stays steadily positive. The unrotated coherence is 3-25x the excess and often the wrong sign. **c is NOT independent of u.** E[u|c] at r=16: c=0: -2.6e-4, c=1: -4.9e-3 (dominant negative), c=2: -2.7e-4, c=3: +7.0e-4 -- the c=1 cell carries a large negative coherence that the rotation maps AWAY from the diagonal.

## P2 FAILED -- the v-term is live and carries the majority
Cells at r=16: P0..P3 = 0.2474, **0.2619**, **0.2431**, 0.2476 (P1 > P2 skew). T2/(q-1/3) = **+0.58** -- the v-term (rotation-1 coherence weighted by the c=1-vs-c=2 cell asymmetry) carries ~58% of the excess at r=16; T1 carries ~42% (and T1 was NEGATIVE for r=4..11, flipping positive only at r>=12 -- historically the v-term carried the whole sign).

## What replaces the attenuation picture
The mechanism is NOT "coherence attenuated by 4, sign preserved." It is: **an alternating raw coherence u_r, a skewed carry-cell measure, and a strong structured (c,u,v) joint law that rotates whatever sign u has into a steady positive excess.** Thirteen consecutive levels of positive excess out of an alternating u is the actual phenomenon -- the carry rotation does real, sign-producing work. Wilson's caution ("if u is strongly correlated with c, that correlation is precisely the arithmetic-vs-multiplicative coupling, and it won't be soft") is CONFIRMED as the operative difficulty. Delta (2) is now measured: the joint law is strongly coupled, cell-skewed, with the v-channel majority.

## Status
CARRYLEMMA: **lemma VERIFIED EXACT** (Fractions r=4,5; 2e-13 to r=16; conditioning absorbed; deltas 1,3 closed). **Attenuation-by-4 mechanism REFUTED** (P1: E[u] alternates, wrong sign at r=16, |u| >> excess; P2: cells skewed, v-term = 58% of excess). **Delta-2 answered by measurement: strong structured (c,u) coupling = the x2/x3 tension, fourth hat.** The target is now the SIGN of T1+T2 with both terms explicit mu-weighted cell statistics -- exact language, measured through r=16, mechanism open. Not at stake: R1-R30, R80-R82, CHANNEL_ID, all Thread-3. commit pending.
