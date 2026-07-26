# Probe CHEB -- gating Wilson's pen (Chebyshev form + digit coupling) -- both VERIFIED, two sharp additions: (1) the k=1 covariance positivity is carried by a MINORITY of base positions (only ~16% of t have Cov_t>0; the average is + because their magnitudes win) => any proof must be GLOBAL/weighted, pointwise-Chebyshev is dead on arrival; (2) "more digits coupled => smaller excess" holds for lag-1 dominance (7x) but is NOT monotone beyond it (k=4, k=8 spike = the MOON oscillatory channels). Sandwich obstruction (2^-T = modular inverse, non-monotone) banked as named. Carries shelf dispatched to Hank.

**Date:** 2026-07-25. Probe `probes/probe_cheb.py`, log `logs/cheb_run.log`.

## (1) Chebyshev/covariance form -- VERIFIED (machine precision)
With a_j = rho_r(t+jm'), b_j = rho_r(t+k+jm'), m' = 3^{r-1}:
   sum_t Cov_j(a,b) = (1/3)p_r(k) - (1/9)p_{r-1}(k) = A_r(k)/3^{r+1}
verified rel <= 4e-12 at r=12,14,16, k=0,1,2. k=0: ZERO negative variances (Cauchy-Schwarz; the m=0 theorem in one
line; strict unless fiber-constant). Wilson's delimitation confirmed: same inequality, degenerate (variance) at
lag 0, genuine covariance at lag 1, no proof to inherit -- and now we know why. Depletion-vs-divergence characters
confirmed: gamma(0) diverges (enhanced), gamma(1) -> 0.738 < 1 (depleted; conjecture = depletion eases monotonically).

**NEW (constraint on any proof): at k=1 only ~16% of base positions t have Cov_t > 0** (0.158/0.160/0.161 at
r=12/14/16, stable). The average is positive because the minority's magnitudes outweigh the ~84% negative majority.
=> positivity is a GLOBAL, weighted phenomenon; any pointwise or per-position argument is refuted by the data.

## (2) x4 as nearest-neighbour digit operation -- VERIFIED
digit_k(4X) = [digit_k(X) + digit_{k-1}(X) + carry] mod 3: 0 failures / 2000 random X < 3^30 (4X = X + 3X).
4^2 = (1,2,1)_3 exactly; 4^k = (1+3)^k binomial coupling of k+1 digits (carries beyond k=2). So
q_r(1) - 1/3 = conditional correlation of two TOP 3-adic digits, coupled by the CARRY -- the translation stands.

## Channel ladder vs coupling width -- lag-1 dominance HOLDS, monotone law does NOT
|A_16(k)|/S_16, 3-nmid k: k=1: 1.94e-3 >> k=2: 1.57e-4; but k=4: 2.61e-4 > k=2, k=8: 2.15e-4 -- spikes at the MOON
oscillatory channels. 3|k rows fold to lower-level content (k=9: 4.96e-4, different animal). Lag 1 IS distinguished
by 7x (structural reason = shortest digit coupling, confirmed); "excess decreasing in coupling width" is NOT a law.

## The sandwich obstruction -- banked as named (Wilson)
X = sum_k 3^k 2^{-T_k} and 2^{-T} is a MODULAR INVERSE -- not monotone in T under any order on Z/3^r. The
product-geometric (v,v') measure is FKG-perfect with nothing monotone to act on. Structural obstruction, not a
time-out. Sandwich OPEN.

## Status
CHEB: covariance identity verified (4e-12); m=0-free-by-Cauchy-Schwarz confirmed (0 negative variances); ~16%
pointwise-positive fact NEW (kills pointwise routes; proof must be weighted/global); digit recursion verified
(0/2000); lag-1 dominance 7x confirmed, monotone-width law refuted; sandwich obstruction named; carries shelf
(Holte / Diaconis-Fulman) dispatched to Hank -- first shelf whose conclusions have the right TYPE (signed digit
distributions, not magnitude bounds). Not at stake: R1-R30, R80-R82, CHANNEL_ID identity (d1 = A_r(1)/S_r), all
Thread-3. Hank pending.
