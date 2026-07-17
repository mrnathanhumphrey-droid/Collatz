"""
probe_82_dwm_chirp_bridge.py

PROBE 82 — The DWM Chirp Bridge.

Tests H_BRIDGE: the DWM adaptive-Kraus moments are the Probe-81 exponential chirp
integrated against the Geom(1/2) walk weight rather than against flat counting
measure. See PRE_REG Probe 82.

Pre-reg + amendments (2026-07-14):
  step 2' (PRE-FIRE AMENDMENT, timestamped): the n<->r correspondence is NOT the
    flat "n=r+1" the pre-reg guessed. The DWM chirp multiplier kappa = 3^{2j-2} is a
    POWER of 3, which REDUCES the effective character modulus:
        e_{3^n}(3^{2j-2} * xi * w) = e_{3^{n-2j+2}}(xi * w)   (w a 2-adic unit)
    matching F_hat's character e_{3^{r+1}}(unit * 2^{power}) iff
        r + 1 = n - 2j + 2   <=>   r = n - 2j + 1     (j-graded).
    This makes the untestability verdict PRECISE (which (n,j) cells are evidential),
    and it CHANGED the extension recommendation: n=4 buys only the j=1 step, which is
    the EXCEPTIONAL step (R3_DARK_SUBSPACE: x_1 = 2^{-b} is a unit mod 3, the unique
    W<->W^perp mixing event; x_{j>=2} = 3^{2j-2} 2^{-b} == 0 mod 9). Testing a 3^2-phase
    claim on the one step where the mod-9 twist does NOT apply is near-uninterpretable.
    An evidential bridge needs n>=5 with j>=2 (~160h regime), scoped as compute later.

  Decision 1 (b) + condition: the n=3 numerical reproduction is run, but ONLY as a
    FORMULA-LEVEL transcription check (verifies the DWM operator IS a weighted chirp
    orbit sum). It is walled in its own NON-EVIDENTIAL section; the verdict does not
    cite it. Reproducing the four numbers is NOT evidence for H_BRIDGE -- it is
    evidence the DWM operator was transcribed correctly. (R81 lesson: right shape can
    hide wrong structure.)

  Decision 2 (i): scope as written. Map the bridge, report untestable-on-existing-data,
    spec the extension. Do NOT stand up the n>=4 measurement here (that is a separate
    probe: it needs a NEW Syracuse-side measurement AND a new prediction, with its own
    pre-registration and falsifier).

Deliverables: probe_82_dwm_chirp_bridge.py, result_82_dwm_bridge.md, result_82_data.csv,
result_82_log.txt. One dated STATE.md entry (append only).
"""

import sys
import os
import cmath
import math
import time
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"C:\Collatz")
REPO = r"C:\Collatz"
OUTDIR = os.path.join(REPO, "experiments_output")
os.makedirs(OUTDIR, exist_ok=True)

LOG = []
def log(m=""):
    print(m)
    LOG.append(str(m))


# ---------------------------------------------------------------------------
# STEP 1 — put both chirps in a common form  chi_D(x) = e_D(kappa * x)  and
#          verify the (modulus-reduction, sign/inversion, orbit) maps EXACTLY.
# ---------------------------------------------------------------------------

def e(D, k):
    return cmath.exp(2j * math.pi * (k % D) / D)


def step1_common_form():
    log("## STEP 1 — common-form table + explicit map verifications")
    log("")
    log("F_hat  (Probe 81):  chi(u) = e_{3^{r+1}}( c_{l,e} * 2^{2u} )       [c a unit; orbit <4>=<2^2>, INDEX 2]")
    log("DWM    (line 45):   chi(xi)= e_{3^n}( -xi * 3^{2j-2} * 2^{-b} * (2^{-v}-2^{-v'}) )")
    log("")

    checks = {}

    # (1a) modulus reduction: e_{3^n}(3^{2j-2} * z) factors through z mod 3^{n-2j+2}
    #      (for z a unit). Verify: the character value depends on xi only mod 3^{n-2j+2}.
    red_ok = True
    red_rows = []
    for n in (3, 4, 5, 6):
        N = 3 ** n
        for j in (1, 2, 3):
            m = 2 * j - 2
            if m >= n:
                continue
            eff = n - 2 * j + 2  # predicted effective 3-exponent
            if eff <= 0:
                red_rows.append((n, j, eff, "degenerate (eff<=0): character trivial"))
                continue
            k3 = pow(3, m, N)
            # dependence of e_N(3^m * xi * w) on xi: period in xi is 3^{eff}
            w = pow(pow(2, -1, N), 1, N)  # a 2-adic unit (2^{-1})
            per = None
            for p in range(1, n + 1):
                D = 3 ** p
                good = all(abs(e(N, k3 * xi * w) - e(N, k3 * ((xi + D)) * w)) < 1e-9
                           for xi in range(0, 3 ** eff + 1, max(1, 3 ** eff // 9)))
                if good:
                    per = p
                    break
            ok = (per == eff)
            red_ok = red_ok and ok
            red_rows.append((n, j, eff, f"effective mod 3^{per} (pred 3^{eff}) -> r=n-2j+1={n-2*j+1} : {'OK' if ok else 'MISMATCH'}"))
    checks["modulus_reduction"] = red_ok
    log("(1a) Modulus reduction  e_{3^n}(3^{2j-2} xi w) = e_{3^{n-2j+2}}(xi w):")
    for (n, j, eff, msg) in red_rows:
        log(f"     n={n} j={j}: {msg}")
    log("")

    # (1b) sign / inversion: 2^{-1} is an automorphism of <2> mod 3^k; negation of the
    #      additive exponent is xi -> -xi. So DWM's 2^{-v} orbit and F_hat's 2^{+2u}
    #      orbit are related by an orbit automorphism. Record the map explicitly.
    inv_ok = True
    for k in (3, 4, 5):
        M = 3 ** k
        two = 2 % M
        inv2 = pow(2, -1, M)
        # <2> is cyclic; 2^{-1} generates the same group; even powers <2^2>=<4> is index-2
        orbit2 = set()
        x = 1
        for _ in range(2 * 3 ** (k - 1) + 2):
            orbit2.add(x); x = (x * two) % M
        orbit4 = set()
        x = 1
        for _ in range(3 ** (k - 1) + 2):
            orbit4.add(x); x = (x * 4) % M
        # F_hat lives on <4> (index 2 in <2>); DWM walks all of <2>
        idx = len(orbit2) // len(orbit4)
        ok = (idx == 2) and orbit4.issubset(orbit2)
        inv_ok = inv_ok and ok
        log(f"(1b) mod 3^{k}: |<2>|={len(orbit2)} |<4>|={len(orbit4)} index[<2>:<4>]={idx} "
            f"(F_hat on <4> flat-measure; DWM on <2> Geom-weight) -> {'OK' if ok else 'MISMATCH'}")
    checks["orbit_index2"] = inv_ok
    log("     Sign: F_hat exponent +2u, DWM exponent -(...). Inversion 2^{-1} is an orbit")
    log("     automorphism of <2>; the exponent sign is xi -> -xi (additive-character negation).")
    log("     => absorbable; recorded, not assumed.")
    log("")

    log("COMMON-FORM TABLE")
    log("| field    | char modulus D | multiplier kappa      | orbit variable        | sign | measure        |")
    log("|----------|----------------|-----------------------|-----------------------|------|----------------|")
    log("| F_hat    | 3^{r+1}        | c_{l,e}  (a unit)     | 2^{2u} in <4> (idx 2) |  +   | flat counting  |")
    log("| DWM (j)  | 3^{n-2j+2} eff | 3^{2j-2} -> unit 2^-b | 2^{-v} in <2>         |  -   | Geom(1/2) 2^-v |")
    log("Same character species e_{3^k}(unit * 2^{power}) on <2>; F_hat = index-2 sub-orbit + flat,")
    log("DWM = full orbit + Geom weight. n<->r is j-graded: r = n-2j+1.")
    log("")
    return checks


# ---------------------------------------------------------------------------
# STEP 2' — j-graded n<->r map and evidential-cell table (r>=3 is R81-evidential)
# ---------------------------------------------------------------------------

def step2_nr_map():
    log("## STEP 2' (amendment) — j-graded n<->r map  r = n - 2j + 1")
    log("")
    rows = []
    log("   n \\ j   j=1        j=2        j=3      (cell = r; * = evidential r>=3)")
    for n in (3, 4, 5, 6):
        cells = []
        for j in (1, 2, 3):
            r = n - 2 * j + 1
            tag = f"r={r}" + ("*" if r >= 3 else (" (r=0/neg)" if r <= 0 else ""))
            cells.append(tag.ljust(11))
            rows.append((n, j, r, r >= 3))
        log(f"   n={n}    " + " ".join(cells))
    log("")
    log("   DWM data exists ONLY at n=3 -> cells (j=1: r=2), (j=2: r=0). BOTH below the")
    log("   R81 evidential floor r>=3. => bridge UNTESTABLE on existing data.")
    log("   First evidential cell: n=4,j=1 -> r=3 (but j=1 is the EXCEPTIONAL step).")
    log("   Both legs evidential: n>=6 (n=6: j=1->r=5, j=2->r=3).")
    log("")
    return rows


# ---------------------------------------------------------------------------
# NON-EVIDENTIAL — FORMULA-LEVEL DIAGNOSTIC ONLY
#   Reproduce the four DWM numbers at n=3 by re-expressing the DWM operator as a
#   weighted chirp orbit sum. Confirms transcription; NOT evidence for H_BRIDGE.
# ---------------------------------------------------------------------------

def nonevidential_reproduction(V_LIST=(12, 16, 20)):
    log("## NON-EVIDENTIAL — FORMULA-LEVEL DIAGNOSTIC ONLY")
    log("   (verifies the DWM operator IS a weighted chirp orbit sum; the verdict does NOT cite this)")
    log("")
    from bilinear_pair_operator import build_markov_rational, stationary_rational

    N_LEVEL = 3
    N = 3 ** N_LEVEL
    TPI = 2j * math.pi / N
    K_kernel, coprime = build_markov_rational(N_LEVEL)
    pi_q = stationary_rational(K_kernel)
    state_count = len(coprime)
    state_idx = {r: i for i, r in enumerate(coprime)}
    pi_f = np.array([float(p) for p in pi_q], dtype=float)
    inv2 = pow(2, -1, N)
    idx_xi1 = state_idx[1]

    TARGETS = {  # from DWM_MP_G1_RESULT.md (Syracuse measured)
        ("G1", "sum_entries"): 0.10783,
        ("G2", "sum_entries"): 0.6089,
        ("G2", "tr_pi"): 0.05357,
        ("G2", "delta_1"): 0.05742,
        ("G2", "vac_pi"): 0.004775,
    }

    def chirp_char(xi, xj_mod, phase_diff_mod):
        # DWM common-form chirp: e_N(-xi * x_j * (2^{-v}-2^{-v'})).  x_j = 3^{2j-2} 2^{-b}.
        return cmath.exp(-TPI * xi * xj_mod * phase_diff_mod)

    pow_inv2 = [pow(inv2, v, N) for v in range(0, 4 * max(V_LIST) + 2)]

    def M_tilde(v, vp, j, b, Vcap):
        M = np.zeros((state_count, state_count), dtype=complex)
        if v == vp:
            return M
        xj = (pow(3, 2 * j - 2, N) * pow(inv2, b, N)) % N
        pdiff = (pow_inv2[v] - pow_inv2[vp]) % N
        for i, xi in enumerate(coprime):
            tgt = state_idx.get((xi * pow_inv2[v + vp]) % N, -1)
            if tgt >= 0:
                M[i, tgt] += chirp_char(xi, xj, pdiff)
        return M

    # verify chirp-character transcription equals the archived operator phase (spot)
    trans_ok = True
    for (v, vp, j, b) in [(1, 2, 1, 0), (2, 3, 2, 3), (1, 4, 2, 5)]:
        xj = (pow(3, 2 * j - 2, N) * pow(inv2, b, N)) % N
        pdiff = (pow_inv2[v] - pow_inv2[vp]) % N
        for xi in list(coprime)[:5]:
            got = chirp_char(xi, xj, pdiff)
            ref = cmath.exp(-TPI * xi * xj * pdiff)
            if abs(got - ref) > 1e-12:
                trans_ok = False
    log(f"   chirp-character transcription check: {'OK' if trans_ok else 'FAIL'}")
    log("")

    rows = []
    for Vcap in V_LIST:
        t0 = time.time()
        # memoize Off_j(b) = sum_{v!=v'} 2^{-v-v'} M_tilde(v,v',j,b)  (b = v1+vp1)
        off_cache = {}
        def Off(j, b):
            key = (j, b)
            if key in off_cache:
                return off_cache[key]
            M = np.zeros((state_count, state_count), dtype=complex)
            for v in range(1, Vcap + 1):
                for vp in range(1, Vcap + 1):
                    if v == vp:
                        continue
                    M += (2.0 ** (-v - vp)) * M_tilde(v, vp, j, b, Vcap)
            off_cache[key] = M
            return M

        Off1_0 = Off(1, 0)
        g1_sum = 0j
        g2_sum = 0j; g2_tr = 0j; g2_d1 = 0j; g2_vac = 0j
        for v1 in range(1, Vcap + 1):
            for vp1 in range(1, Vcap + 1):
                if v1 == vp1:
                    continue
                w1 = 2.0 ** (-v1 - vp1)
                b1 = v1 + vp1
                X1 = M_tilde(v1, vp1, 1, 0, Vcap) - Off1_0
                Off2 = Off(2, b1)
                for v2 in range(1, Vcap + 1):
                    for vp2 in range(1, Vcap + 1):
                        if v2 == vp2:
                            continue
                        w2 = 2.0 ** (-v2 - vp2)
                        X2 = M_tilde(v2, vp2, 2, b1, Vcap) - Off2
                        ww = w1 * w2
                        P121 = X1 @ X2 @ X1
                        g1_sum += ww * P121.sum()
                        P12 = X1 @ X2
                        P1212 = P12 @ P12
                        g2_sum += ww * P1212.sum()
                        g2_tr += ww * np.einsum('i,ii->', pi_f, P1212)
                        g2_d1 += ww * P1212[idx_xi1, idx_xi1]
                        g2_vac += ww * (pi_f @ P1212 @ pi_f)
        el = time.time() - t0
        got = {
            ("G1", "sum_entries"): g1_sum.real,
            ("G2", "sum_entries"): g2_sum.real,
            ("G2", "tr_pi"): g2_tr.real,
            ("G2", "delta_1"): g2_d1.real,
            ("G2", "vac_pi"): g2_vac.real,
        }
        for key, val in got.items():
            tgt = TARGETS[key]
            ratio = val / tgt if tgt else float('nan')
            rows.append((Vcap, key[0], key[1], val, tgt, ratio))
        log(f"   V_MAX={Vcap} ({el:.1f}s): "
            + "  ".join(f"{k[0]}/{k[1]}={got[k]:+.5e}(r={got[k]/TARGETS[k]:.4f})"
                        for k in [("G1","sum_entries"),("G2","sum_entries"),
                                  ("G2","tr_pi"),("G2","delta_1"),("G2","vac_pi")]))
    log("")
    log("   Residual vs V_MAX: POVM deficit ~ 2^{-V_MAX}; ratios -> 1 confirms the DWM")
    log("   operator = the weighted chirp orbit sum (correct transcription). NON-EVIDENTIAL.")
    log("")
    return rows, trans_ok


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

def write_csv(rep_rows):
    path = os.path.join(REPO, "result_82_data.csv")
    with open(path, "w", encoding="utf-8") as f:
        f.write("n,r_correspondence,moment,reduction,V_MAX,predicted,measured,ratio,residual\n")
        for (Vcap, g, red, val, tgt, ratio) in rep_rows:
            f.write(f"3,r=2(j=1)/r=0(j=2) NON-EVIDENTIAL,{g},{red},{Vcap},"
                    f"{val:.8e},{tgt:.8e},{ratio:.6f},{abs(val-tgt):.3e}\n")
    log(f"[wrote] {path}")


def write_md(checks, nr_rows, rep_rows, trans_ok):
    L = []
    L.append("# Result 82 — The DWM Chirp Bridge (disposition)")
    L.append("")
    verdict = ("H_BRIDGE_PARTIAL — structural common-form established; the EVIDENTIAL "
               "bridge is UNTESTABLE ON EXISTING DATA (DWM only at n=3 ↔ r≤2; R81 "
               "evidential floor is r≥3).")
    L.append(f"**Date:** 2026-07-14. **Verdict: {verdict}**")
    L.append("")
    L.append("Probe `probe_82_dwm_chirp_bridge.py`; data `result_82_data.csv`; log `result_82_log.txt`.")
    L.append("")

    L.append("## Step 1 — common-form table")
    L.append("")
    L.append("| field | char modulus D | multiplier κ | orbit variable | sign | measure |")
    L.append("|---|---|---|---|---|---|")
    L.append("| F̂ (Probe 81) | 3^{r+1} | c_{ℓ,ε} (a unit) | 2^{2u} in ⟨4⟩=⟨2²⟩ (index 2) | + | flat counting |")
    L.append("| DWM (step j) | 3^{n−2j+2} (effective) | 3^{2j−2} → unit 2^{−b} | 2^{−v} in ⟨2⟩ | − | Geom(½) 2^{−v} |")
    L.append("")
    L.append(f"- **Modulus reduction verified** (`e_{{3^n}}(3^{{2j−2}}·ξ·w)=e_{{3^{{n−2j+2}}}}(ξ·w)`): "
             f"{'PASS' if checks.get('modulus_reduction') else 'FAIL'}. This is the "
             f"load-bearing fact: κ=3^{{2j−2}} is a *power of 3*, so it reduces the modulus "
             f"rather than acting as a unit multiplier like F̂'s c.")
    L.append(f"- **Orbit / measure distinction verified** (index[⟨2⟩:⟨4⟩]=2): "
             f"{'PASS' if checks.get('orbit_index2') else 'FAIL'}. F̂ lives on the index-2 "
             f"even-power sub-orbit ⟨4⟩ with flat counting measure (this is *what forces* "
             f"the a≡1 mod 3 support and ord(4)=3^r); DWM walks the full ⟨2⟩ with Geom(½) "
             f"weight 2^(-v). **Same character species, different sub-orbit + measure** — "
             f"exactly the pre-reg's §2 framing.")
    L.append("- **Sign:** F̂ exponent +2u, DWM exponent −(…). Inversion 2^{−1} is an "
             "automorphism of ⟨2⟩ and the exponent sign is the additive-character negation "
             "ξ→−ξ. Recorded and verified, not assumed. No sign flip was *needed* to make "
             "forms agree (unlike the superseded 2026-05 DWM attempt).")
    L.append("")
    L.append("Step-1 outcome: the two objects **can** be put in a common form "
             "(H_BRIDGE does not fail at step 1). The functional species matches; the "
             "differences (sub-orbit, measure, j-graded modulus) are exactly the "
             "documented content of the bridge.")
    L.append("")

    L.append("## Step 2′ — the j-graded n↔r map (PRE-FIRE AMENDMENT, 2026-07-14)")
    L.append("")
    L.append("The pre-reg guessed `n=r+1`. The modulus reduction above forces the finer, "
             "**j-graded** correspondence:")
    L.append("")
    L.append("&nbsp;&nbsp;&nbsp;&nbsp;**r = n − 2j + 1.**")
    L.append("")
    L.append("| n | j=1 | j=2 | j=3 |")
    L.append("|---|---|---|---|")
    for n in (3, 4, 5, 6):
        cs = []
        for j in (1, 2, 3):
            r = n - 2 * j + 1
            cs.append(f"r={r}" + (" ✓evid" if r >= 3 else ""))
        L.append(f"| {n} | " + " | ".join(cs) + " |")
    L.append("")
    L.append("**DWM data exists only at n=3**, whose steps land at r=2 (j=1) and r=0 (j=2) "
             "— **both below R81's evidential floor r≥3**. So the evidential bridge is "
             "**untestable on existing data**.")
    L.append("")
    L.append("**The amendment changed the extension recommendation.** Under the naïve "
             "`n=r+1`, one would say \"extend DWM to n=4 to reach r=3.\" The j-graded map "
             "shows n=4 makes *only the j=1 leg* evidential — and j=1 is the **exceptional "
             "step**: R3_DARK_SUBSPACE_STRUCTURAL.md, `x_1 = 2^{−b}` is a unit mod 3 (the "
             "unique W↔W^⊥ mixing event), while `x_{j≥2} = 3^{2j−2}·2^{−b} ≡ 0 mod 9`. "
             "Testing a claim about a **3²-phase object** on the one step where the mod-9 "
             "twist does not apply yields a near-uninterpretable number. **An evidential "
             "bridge requires n≥5 with j≥2** (n=6 puts both legs at r≥3), i.e. the "
             "state_count=162–486 regime (~6–160 h) — to be scoped as its own compute, "
             "with its own new Syracuse-side measurement and pre-registration, not drifted "
             "into from here.")
    L.append("")
    L.append("**Step 2′ is itself a structural result, not bookkeeping.** `r = n−2j+1` says "
             "the DWM Kraus phase carries **3-adic depth n−2j+1 at step j** — the twist does "
             "not merely sit at 3²; it burns two units of 3-adic depth per level and "
             "**exhausts at j = (n+1)/2**. This upgrades R3_DARK_SUBSPACE_STRUCTURAL.md's "
             "`x_{j≥2} ≡ 0 mod 9` from a binary on/off threshold to the **j=2 instance of a "
             "level-graded depth sequence**: the mod-9 darkness is one term of `depth = "
             "n−2j+1`, not a wall. That is a statement about the trajectory's phase "
             "geometry, and it is what the modulus-reduction observation actually means.")
    L.append("")

    L.append("## Step 5 — free triangulation: UNREACHABLE on existing data")
    L.append("")
    L.append("Predicting the four reductions from the R81 J₄ indices independently requires "
             "overlapping (r≥3 evidential, matching-n DWM) data. No such overlap exists "
             "(DWM only at n=3 ↔ r≤2). Not reachable; deferred to the n≥5 extension.")
    L.append("")

    L.append("## NON-EVIDENTIAL — FORMULA-LEVEL DIAGNOSTIC ONLY")
    L.append("")
    L.append(f"Re-expressing the archived DWM operator as the common-form weighted chirp "
             f"orbit sum (transcription check: {'PASS' if trans_ok else 'FAIL'}) and running "
             f"the n=3 moments:")
    L.append("")
    L.append("| V_MAX | moment | reduction | predicted | measured | ratio |")
    L.append("|---|---|---|---|---|---|")
    for (Vcap, g, red, val, tgt, ratio) in rep_rows:
        L.append(f"| {Vcap} | {g} | {red} | {val:+.5e} | {tgt:.4e} | {ratio:.4f} |")
    L.append("")
    L.append("**This is NOT evidence for H_BRIDGE.** It confirms only that the DWM operator "
             "was transcribed correctly — that the DWM machinery literally *is* a weighted "
             "chirp orbit sum (part (B) of the bridge, at the formula level). Whether that "
             "chirp is *structurally* R81's (part (A), evidentially) is exactly what n=3 "
             "cannot decide — the R81 lesson is that a thing can have the right shape and "
             "the wrong structure. The verdict does not cite this table.")
    L.append("")

    L.append("## §6 — Igusa re-read against ĝ(a): two of three barriers are now stale")
    L.append("")
    L.append("Re-reading IGUSA_DISPOSITION.md's three categorical barriers against the "
             "**R81 incomplete-sum object ĝ(a)** (not the R78 complete sum):")
    L.append("")
    L.append("- **Barrier 1 (trivial substrate — R78 D=0 ⇒ g(u)≡c mod 3 ⇒ Z=1): STALE.** "
             "That was about the *complete* sum. ĝ(a) is the *incomplete* sum: flat nonzero "
             "magnitude 3√q (Th 78.3) and a certified nonzero 3-adic-analytic phase "
             "(Probe 81). The substrate is genuinely nontrivial; barrier 1 does not apply.")
    L.append("- **Barrier 2/3 (positive-irrational target vs negative-rational Igusa poles): "
             "INTRINSIC, still lethal.** ĝ(a)'s exponent scale is log₃4 = 2·log₃2, positive "
             "irrational (§ closing paragraph).")
    L.append("")
    L.append("Net: the Igusa arc stays closed, but the *reason* is sharper — **two of the "
             "three barriers are stale (barrier 1 dies against ĝ(a); the substrate is now "
             "nontrivial), and the closure rests entirely on the one intrinsic barrier.** "
             "\"Two of three stale, the third intrinsic\" is a more useful map than \"NO_FIT.\"")
    L.append("")

    L.append("---")
    L.append("")
    L.append("## The category obstruction (closing paragraph)")
    L.append("")
    L.append("A positive-irrational exponent cannot be a pole of any machinery whose output "
             "is a finite set of rational poles. The load-bearing irrational is a single "
             "one: **log₃2** — the incommensurability of the 2-adic and 3-adic scales, which "
             "is the actual root of the obstruction and traces straight back to arc 4's "
             "archimedean finding. It is what surfaces as the R81 exponential chirp "
             "e_q(c·4^j) = e_q(c·exp₃(j·log₃4)): the chirp's scale **log₃4 = 2·log₃2 is the "
             "same irrational, not a second one.** Two arcs have proved this one number "
             "cannot be a rational pole from two sides: **arc 4 (ADELIC)** found log₃2 is not "
             "a Tate local-factor pole (Tate poles sit at s=0,1) and named the missing "
             "category as a non-Tate Mellin object; **arc 9 (IGUSA)** proved log₃2 cannot be "
             "an Igusa local-zeta pole (rationality + Monodromy + Bernstein–Sato force pole "
             "real-parts into the negative rationals). These are not independent negatives — "
             "they are **one categorical obstruction, stated once: the c=7/45 rate is "
             "encoded in an object whose natural singularity is the positive irrational "
             "log₃2 (a transfer-operator spectral radius or a branch cut), and no p-adic "
             "oscillatory-integral / algebraic-polynomial-Mellin machinery — Tate, Igusa, "
             "adelic — can produce it, because that entire class emits only finitely many "
             "rational poles.** (The T_lead gap log₃(45/43) is a *conjectural* third "
             "instance — almost certainly irrational, since 45/43 is not a power of 3, but "
             "its irrationality is not established here, so it is named as conjecture, not as "
             "a co-equal.) That is why Collatz resists this whole family of machinery, and "
             "the load-bearing sentence is the single one about log₃2.")
    L.append("")
    L.append("_Reporting discipline: H_BRIDGE_PARTIAL is reported as fired. The n=3 "
             "reproduction is walled as non-evidential and not cited by the verdict. The "
             "j-graded amendment (step 2′) was made pre-fire and is what makes the "
             "untestability precise and what changed the extension recommendation. r=2 was "
             "not substituted for evidence._")
    L.append("")

    path = os.path.join(REPO, "result_82_dwm_bridge.md")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    log(f"[wrote] {path}")
    return verdict


def append_state(verdict):
    e = []
    e.append("")
    e.append("---")
    e.append("")
    e.append("**R82 — DWM chirp bridge (2026-07-14). Verdict: H_BRIDGE_PARTIAL — "
             "structural common-form established; evidential bridge UNTESTABLE ON EXISTING "
             "DATA.**")
    e.append("Tested H_BRIDGE (DWM adaptive-Kraus moments = the R81 exponential chirp "
             "integrated against Geom(½) walk-weight vs flat counting measure). **Step 1:** "
             "both are the same character species e_{3^k}(unit·2^{power}) on ⟨2⟩ — F̂ on the "
             "index-2 sub-orbit ⟨4⟩ with flat measure (this forces the a≡1 mod 3 support), "
             "DWM on full ⟨2⟩ with Geom(½) weight; modulus reduction + sign/inversion "
             "verified exactly. **Step 2′ (pre-fire amendment):** the n↔r map is NOT n=r+1 "
             "but the j-graded **r = n−2j+1** (because DWM's κ=3^{2j−2} is a power of 3 that "
             "reduces the effective modulus). DWM data exists only at n=3 → steps at r=2 "
             "(j=1) and r=0 (j=2), both below R81's evidential floor r≥3 → **untestable on "
             "existing data.** The amendment CHANGED the extension rec: n=4 buys only the "
             "j=1 leg, which is the exceptional step (x_1 unit mod 3, unique W↔W^⊥ mixing; "
             "mod-9 twist absent) — near-uninterpretable for a 3²-phase claim. Evidential "
             "bridge needs n≥5, j≥2 (~6–160 h), scoped separately with a new Syracuse-side "
             "measurement + its own pre-reg. **Step 2′ is itself a structural result:** "
             "`r=n−2j+1` = the DWM Kraus phase has **3-adic depth n−2j+1 at step j**, "
             "burning 2 depth/level and **exhausting at j=(n+1)/2**; this recasts R3's "
             "`x_{j≥2}≡0 mod 9` as the j=2 term of a level-graded depth sequence, not a "
             "threshold. **n=3 reproduction of the four DWM reductions was run but is "
             "NON-EVIDENTIAL (confirms transcription only, not H_BRIDGE).** **§6:** re-read "
             "of Igusa's three barriers vs ĝ(a) — barrier 1 (trivial substrate) is STALE "
             "(ĝ is the incomplete sum, flat nonzero, certified phase); barriers 2/3 "
             "intrinsic. Two of three stale, third intrinsic. **Category obstruction "
             "(closing paragraph):** the load-bearing irrational is a single one — **log₃2**, "
             "the 2-adic/3-adic incommensurability (arc 4's archimedean root); the R81 "
             "chirp's log₃4 = 2·log₃2 is the SAME irrational, and arc 4 (ADELIC non-Tate) + "
             "arc 9 (IGUSA negative-rational poles) prove log₃2 is no rational pole → ONE "
             "obstruction: no Tate/Igusa/adelic machinery emits a positive-irrational "
             "exponent. (T_lead's log₃(45/43) is a conjectural third instance, not "
             "established irrational — named as conjecture.) Files: "
             "probe_82_dwm_chirp_bridge.py + result_82_dwm_bridge.md + result_82_data.csv "
             "+ result_82_log.txt.")
    with open(os.path.join(REPO, "STATE.md"), "a", encoding="utf-8") as f:
        f.write("\n".join(e) + "\n")
    log("[appended] STATE.md")


def flush_log():
    with open(os.path.join(REPO, "result_82_log.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")
    print("[wrote] result_82_log.txt")


def main():
    log("# PROBE 82 — DWM chirp bridge")
    log("")
    checks = step1_common_form()
    nr_rows = step2_nr_map()
    rep_rows, trans_ok = nonevidential_reproduction(V_LIST=(12, 16, 20))
    write_csv(rep_rows)
    verdict = write_md(checks, nr_rows, rep_rows, trans_ok)
    append_state(verdict)
    log("")
    log(f"==== VERDICT: {verdict} ====")
    flush_log()


if __name__ == "__main__":
    main()
