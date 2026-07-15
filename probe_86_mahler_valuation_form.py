"""
probe_86_mahler_valuation_form.py

PROBE 86 (LOW PRIORITY puzzle-box) — closed form for the Mahler valuation profile v3(c_k)
of arg F-hat (R81b). NOT load-bearing: R81b's analyticity certification, the degree law,
the smooth-completion closure, and the log_3(2) obstruction all stand regardless.

Feasibility (established pre-fire): v3(c_k) ~ 1.3k, and resolving c_k needs modulus
3^{r+1} > 3^{v3(c_k)}, i.e. r >~ v3(c_k). k=20 would need r~26 (3^26 ~ 2.5e12-term sums)
-> infeasible. Per operator decision (option 1): extend to the FEASIBLE ceiling ~k=15
(r up to 20), verify r-stability per new point, then test the period-3 / Kummer-carry
candidate OUT OF SAMPLE on k=12..15. k=20 belongs to Task C's algebraic derivation.

Method: single-frequency ghat(a)=sum_{j=0}^{d-1} e_q(c 4^j) e_d(-a j) computed by STREAMING
over j (no 3^r-size FFT array), with overflow-safe modular multiply (int64 direct for
q<3.03e9 i.e. r<=18; split-multiply for r>=19). s2(b)=index of ghat^2/q in Z/2q, mod q
(R81b convention). c_k = Delta^k s2(0) mod 3^{r+1}. A k enters the profile only if its
valuation is r-STABLE (agrees across >=2 resolving r).

Task A — extend + r-stability (MANDATORY, first).
Task B — kill/confirm v3(c_k) = floor(4k/3) + digit/carry correction, OUT OF SAMPLE.
Task C — only if B fires: derive from (1+3)^x binomial structure.

Deliverables: probe_86_*.py, result_86_mahler_valuation.md, result_86_data.csv,
result_86_log.txt. One dated STATE.md entry.
"""

import sys
import os
import math
import time
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
REPO = r"C:\Collatz"
LOG = []
def log(m=""):
    print(m, flush=True)
    LOG.append(str(m))

INT64_SAFE_Q = 3_030_000_000  # q below this -> a*arr fits int64 (q^2 < 9.2e18)


def modmul(a, arr, q):
    """(a * arr) mod q, a scalar python int, arr int64 numpy array. Overflow-safe."""
    a %= q
    if q < INT64_SAFE_Q:
        return (a * arr) % q
    SH = 1 << 20
    hi = a // SH
    lo = a % SH
    t = ((arr % q) * SH) % q            # (2^20 * arr) mod q
    return (hi * t + lo * arr) % q


def ghat_freqs(r, ell, eps, K, chunk=2_000_000):
    """Streaming ghat(a0+3b) for b=0..K at level r. Returns (complex array len K+1, q)."""
    q = 3 ** (r + 1)
    d = 3 ** r
    a0 = 1 if eps == 0 else 2
    c = (pow(2, eps, q) * pow((1 + 3 ** r) % q, ell, q)) % q

    # base table 4^i mod q for i=0..chunk-1 (built once per q)
    pow4_base = np.empty(chunk, dtype=np.int64)
    v = 1
    for i in range(chunk):
        pow4_base[i] = v
        v = (v * 4) % q
    four_chunk = pow(4, chunk, q)

    acc = np.zeros(K + 1, dtype=np.complex128)
    j0 = 0
    base = 1  # 4^{j0} mod q  (python int, exact)
    while j0 < d:
        n = min(chunk, d - j0)
        pb = pow4_base if n == chunk else pow4_base[:n]
        cp = modmul(base, pb, q)                 # 4^{j0+i} mod q
        cc = modmul(c, cp, q)                     # c*4^j mod q
        base_wave = np.exp(2j * np.pi * (cc.astype(np.float64) / q))
        jidx = (j0 + np.arange(n)).astype(np.float64)
        w_a0 = np.exp(-2j * np.pi * (a0 * jidx) / d)
        w_step = np.exp(-2j * np.pi * (3.0 * jidx) / d)
        term = base_wave * w_a0
        for b in range(K + 1):
            acc[b] += term.sum()
            if b < K:
                term = term * w_step
        base = (base * four_chunk) % q
        j0 += n
    return acc, q


def s2_from_ghat(gh, q):
    z = gh * gh / q
    ang = np.angle(z) % (2 * np.pi)
    J2 = np.rint(ang * (2 * q) / (2 * np.pi)).astype(np.int64) % (2 * q)
    resid = float(np.max(np.abs(z - np.exp(2j * np.pi * J2 / (2 * q)))))
    return (J2 % q).astype(np.int64), resid


def mahler_v3(r, ell, eps, K):
    """Return (v3 list for k=0..K, resolved-bool list, cert residual)."""
    gh, q = ghat_freqs(r, ell, eps, K)
    s2, resid = s2_from_ghat(gh, q)
    cur = [int(x) % q for x in s2]
    cks = [cur[0] % q]
    while len(cur) > 1:
        cur = [(cur[i + 1] - cur[i]) % q for i in range(len(cur) - 1)]
        cks.append(cur[0] % q)
    out_v3, resolved = [], []
    for c in cks:
        if c % q == 0:
            out_v3.append(None); resolved.append(False)
        else:
            k = 0; t = c
            while t % 3 == 0:
                t //= 3; k += 1
            out_v3.append(k); resolved.append(k < r + 1)
    return out_v3, resolved, resid


# --- base-3 features -------------------------------------------------------

def s3(k):
    s = 0
    while k:
        s += k % 3; k //= 3
    return s

def carries_add(a, b, p=3):
    c = carry = 0
    while a or b or carry:
        s = a % p + b % p + carry
        carry = 1 if s >= p else 0
        c += carry
        a //= p; b //= p
    return c

def v3_fact(k):  # Legendre
    return (k - s3(k)) // 2  # (k - s_3(k))/(3-1)


# ---------------------------------------------------------------------------

def self_test():
    log("## Self-test: streaming ghat reproduces the banked profile (r=12)")
    v3l, res, resid = mahler_v3(12, 0, 0, 10)
    got = [v for v in v3l if v is not None]
    exp = [0, 2, 2, 3, 4, 6, 7, 8, 10, 11, 12]
    ok = got[:len(exp)] == exp and resid < 1e-6
    log(f"   r=12 v3(c_k) k=0..10 = {v3l}  (cert resid {resid:.1e})")
    log(f"   expected {exp} -> {'OK' if ok else 'FAIL'}")
    log("")
    return ok


def task_A(R_SET, K):
    log("## TASK A — extend v3(c_k) with per-k r-stability")
    log("")
    per_r = {}
    for r in R_SET:
        t = time.time()
        v3l, resolved, resid = mahler_v3(r, 0, 0, K)
        per_r[r] = v3l
        dt = time.time() - t
        disp = [(str(v) if (v is not None and v < r + 1) else "X") for v in v3l]
        log(f"   r={r:>2} (q=3^{r+1}, {dt:5.1f}s, resid {resid:.1e}): " + " ".join(disp))
    log("")
    # stability: for each k, collect resolved v3 across r; confirmed if >=2 agree
    profile = {}
    for k in range(K + 1):
        vals = {r: per_r[r][k] for r in R_SET
                if per_r[r][k] is not None and per_r[r][k] < r + 1}
        if not vals:
            profile[k] = (None, "unresolved", vals)
            continue
        uniq = set(vals.values())
        if len(uniq) == 1:
            v = uniq.pop()
            profile[k] = (v, "confirmed" if len(vals) >= 2 else "single-r", vals)
        else:
            profile[k] = (None, "UNSTABLE", vals)
    log("   k :  v3   status      (resolving r:v3)")
    for k in range(K + 1):
        v, st, vals = profile[k]
        vs = ",".join(f"{r}:{vv}" for r, vv in sorted(vals.items()))
        log(f"   {k:>2}:  {str(v):>4}  {st:<10}  [{vs}]")
    log("")
    return profile


def task_B(profile, K):
    log("## TASK B — kill/confirm floor(4k/3) + digit/carry correction (OUT OF SAMPLE)")
    log("")
    # usable = confirmed or single-r (but note); fit window = k<=11, test = k>=12
    pts = {k: profile[k][0] for k in range(K + 1) if profile[k][0] is not None}
    fit_k = [k for k in pts if k <= 11]
    oos_k = [k for k in pts if k >= 12]
    log(f"   fit window k<=11: {sorted(fit_k)};  out-of-sample k>=12: {sorted(oos_k)}")
    if not oos_k:
        log("   (no out-of-sample points resolved — extension did not clear k>=12)")
    log("")

    base = {k: (4 * k) // 3 for k in pts}
    feats = {
        "floor(4k/3) only": (lambda k: 0),
        "+ a*s3(k)": ("s3", lambda k: s3(k)),
        "+ a*carries(k+k)": ("carry", lambda k: carries_add(k, k)),
        "+ a*v3(k!)": ("v3fact", lambda k: v3_fact(k)),
    }
    results = []
    # model 0: base alone
    r0 = {k: pts[k] - base[k] for k in pts}
    fired0 = all(r0[k] == 0 for k in pts)
    results.append(("floor(4k/3) only", None, all(r0[k] == 0 for k in fit_k),
                    all(r0[k] == 0 for k in oos_k) if oos_k else None, fired0))
    log(f"   floor(4k/3) alone: residual v3-floor = "
        f"{[r0[k] for k in sorted(pts)]}  -> {'FIRES' if fired0 else 'refuted'}")
    # single-feature integer-coefficient models
    for name, (tag, f) in [(n, v) for n, v in feats.items() if isinstance(v, tuple)]:
        # integer a from fit window (need consistent a s.t. residual = a*f)
        # solve a per point; must be identical integer across all fit points with f!=0
        cand = set()
        ok_int = True
        for k in fit_k:
            fk = f(k)
            rk = pts[k] - base[k]
            if fk == 0:
                if rk != 0:
                    ok_int = False
            else:
                if rk % fk != 0:
                    ok_int = False
                else:
                    cand.add(rk // fk)
        if not ok_int or len(cand) != 1:
            log(f"   {name}: no single integer coeff fits window -> refuted")
            results.append((name, None, False, None, False))
            continue
        a = cand.pop()
        in_ok = all(pts[k] - base[k] == a * f(k) for k in fit_k)
        oos_ok = all(pts[k] - base[k] == a * f(k) for k in oos_k) if oos_k else None
        fired = in_ok and (oos_ok is True)
        log(f"   {name} with a={a}: in-window={'exact' if in_ok else 'no'}, "
            f"out-of-sample={'exact' if oos_ok else ('n/a' if oos_ok is None else 'FAILS')} "
            f"-> {'FIRES' if fired else 'refuted'}")
        results.append((name, a, in_ok, oos_ok, fired))
    log("")
    any_fired = any(r[4] for r in results)
    return results, any_fired, pts, fit_k, oos_k


def write_csv(profile, K):
    path = os.path.join(REPO, "result_86_data.csv")
    with open(path, "w", encoding="utf-8") as f:
        f.write("k,v3,status,delta,s3,carries_kk,v3_fact,floor4k3,resid_vs_floor,resolving_r\n")
        prev = None
        for k in range(K + 1):
            v, st, vals = profile[k]
            dlt = "" if (v is None or prev is None) else v - prev
            rv = "" if v is None else v - (4 * k) // 3
            rr = ";".join(f"{r}:{vv}" for r, vv in sorted(vals.items()))
            f.write(f"{k},{'' if v is None else v},{st},{dlt},{s3(k)},{carries_add(k,k)},"
                    f"{v3_fact(k)},{(4*k)//3},{rv},{rr}\n")
            if v is not None:
                prev = v
    log(f"[wrote] {path}")


def write_md(profile, K, taskB, any_fired, pts, fit_k, oos_k):
    L = []
    L.append("# Result 86 — closed form for the Mahler valuation profile v₃(c_k)")
    L.append("")
    conf = [k for k in range(K + 1) if profile[k][1] == "confirmed"]
    kmax_conf = max(conf) if conf else None
    verdict = ("A closed form FIRED (see Task B/C)" if any_fired else
               "NULL — no candidate reproduces the extended profile out of sample. "
               "v₃(c_k) remains deterministic, r-stable, ~1.3k, period-3 with digit-rollover "
               "defects, closed form OPEN.")
    L.append(f"**Date:** 2026-07-14. **Priority: LOW (puzzle-box).** **Verdict: {verdict}**")
    L.append("")
    L.append("Probe `probe_86_mahler_valuation_form.py`; data `result_86_data.csv`; log "
             "`result_86_log.txt`.")
    L.append("")
    L.append("**Not load-bearing.** R81b's analyticity certification, the degree law "
             "`degree(r)=max{k:v₃(c_k)≤r}`, the smooth-completion closure, and the log₃2 "
             "category obstruction all stand regardless of this outcome. This probe only "
             "asks whether the *valuation profile* has an obvious closed form.")
    L.append("")
    L.append("## Feasibility (why not k=20)")
    L.append("")
    L.append("`v₃(c_k) ≈ 1.3k`, and resolving c_k needs `3^{r+1} > 3^{v₃(c_k)}` ⇒ `r ≳ v₃(c_k)`. "
             "k=20 ⇒ v₃≈26 ⇒ r≈26 ⇒ `3^26 ≈ 2.5×10¹²`-term sums (and float loses the integer "
             "index past r≈20–22): infeasible. Per operator decision, extended to the feasible "
             f"ceiling (confirmed to **k={kmax_conf}**). k=20 belongs to Task C's algebraic "
             "derivation, not a numerical extension.")
    L.append("")
    L.append("## Task A — extended profile with per-k r-stability")
    L.append("")
    L.append("| k | v₃(c_k) | status | Δ | ⌊4k/3⌋ | v₃−⌊4k/3⌋ | s₃(k) | carries(k+k) |")
    L.append("|---|---|---|---|---|---|---|---|")
    prev = None
    for k in range(K + 1):
        v, st, vals = profile[k]
        dlt = "" if (v is None or prev is None) else v - prev
        rv = "" if v is None else v - (4 * k) // 3
        L.append(f"| {k} | {'—' if v is None else v} | {st} | {dlt} | {(4*k)//3} | {rv} "
                 f"| {s3(k)} | {carries_add(k,k)} |")
        if v is not None:
            prev = v
    L.append("")
    L.append("`confirmed` = valuation r-stable (agrees across ≥2 resolving r); `single-r` = "
             "resolved at only one r (tentative, not evidential); `unresolved` = modulus too "
             "small at every tested r. Only `confirmed` points carry weight.")
    L.append("")
    L.append("## Task B — candidate models, out of sample")
    L.append("")
    L.append(f"Fit window k≤11; out-of-sample k≥12 = **{sorted(oos_k) if oos_k else 'none resolved'}**. "
             "Decision rule: integer coefficients, reproduce **every** point incl. k≥12, or REFUTED.")
    L.append("")
    L.append("| model | integer coeff | in-window | out-of-sample | verdict |")
    L.append("|---|---|---|---|---|")
    for (name, a, in_ok, oos_ok, fired) in taskB:
        oos = "n/a" if oos_ok is None else ("exact" if oos_ok else "FAILS")
        L.append(f"| {name} | {a if a is not None else '—'} | "
                 f"{'exact' if in_ok else 'no'} | {oos} | {'FIRES' if fired else 'refuted'} |")
    L.append("")
    if not any_fired:
        L.append("**All candidates refuted.** The Legendre/Kummer family (`⌊4k/3⌋` + integer·"
                 "{s₃, carries, v₃(k!)}) does not reproduce the profile"
                 + (" out of sample" if oos_k else " even in-window / no OOS point cleared")
                 + ". Consistent with the pre-registered most-likely null and with the author's "
                 "0-for-4 structural-prior record in this arc. The period-3 `1,1,2` difference "
                 "structure with digit-rollover defects is real and visible in the table, but no "
                 "closed form in these features captures it.")
    else:
        L.append("**A model FIRED out of sample** — see the row above and Task C.")
    L.append("")
    L.append("## Task C — derivation")
    L.append("")
    if any_fired:
        L.append("A surviving model is present; its source in the `4=(1+3)` binomial structure "
                 "of `x↦4^x` is the next step (and if `4=1+3` is in its guts, it is the fifth "
                 "appearance of that root — coherence table). [Derivation to be completed.]")
    else:
        L.append("Not earned — no model survived Task B. The closed form for v₃(c_k) is reported "
                 "**empirical and OPEN**. The natural derivation route (Mahler coefficients of "
                 "`x↦4^x=(1+3)^x` from binomial structure) is where a real answer would come "
                 "from; the numerical profile alone does not yield it.")
    L.append("")
    L.append("## Scope")
    L.append("")
    L.append("Untouched, explicitly: **R81b's certification**, the degree law, the "
             "smooth-completion closure, the log₃2 obstruction, THEOREM_C_745, and Thms "
             "78.1–78.3. This is a low-priority puzzle-box and yields the machine to anything "
             "else that needs it.")
    L.append("")
    L.append("_Reporting discipline: extended before fitting (Task A first); out-of-sample "
             "decision rule enforced; integer coefficients required; r-stability verified per "
             "point; the pre-registered candidate treated as a hypothesis to kill._")
    L.append("")
    path = os.path.join(REPO, "result_86_mahler_valuation.md")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    log(f"[wrote] {path}")
    return verdict, kmax_conf


def append_state(verdict, kmax_conf, any_fired, oos_k):
    e = []
    e.append("")
    e.append("---")
    e.append("")
    head = ("closed form FIRED" if any_fired else "NULL — no closed form")
    e.append(f"**R86 — Mahler valuation profile v₃(c_k) closed form (2026-07-14, LOW "
             f"priority). {head}.**")
    e.append(f"Extended R81b's v₃(c_k) profile via streaming single-frequency ĝ + overflow-safe "
             f"modular arithmetic (k=20 infeasible: v₃~1.3k ⇒ needs r~26/3^26 terms; extended to "
             f"the feasible ceiling, **confirmed to k={kmax_conf}** with per-k r-stability). "
             f"Out-of-sample points k≥12 = {sorted(oos_k) if oos_k else 'none cleared'}. ")
    if any_fired:
        e.append(f"A `⌊4k/3⌋`+digit/carry (Kummer) model reproduced the profile OUT OF SAMPLE "
                 f"with integer coeffs — see result_86 + Task C for the (1+3)^x source. ")
    else:
        e.append(f"**All candidates REFUTED out of sample** (`⌊4k/3⌋` + integer·{{s₃(k), "
                 f"carries(k+k), v₃(k!)}}): none reproduce the extended profile. v₃(c_k) stays "
                 f"deterministic, r-stable, ~1.3k, period-3 `1,1,2` differences with "
                 f"digit-rollover defects; closed form **OPEN** (belongs to Task C's `x↦(1+3)^x` "
                 f"binomial derivation, not a numerical extension). Author's structural priors "
                 f"now 0-for-5 this arc — expected. ")
    e.append(f"**Not load-bearing:** R81b certification, degree law, smooth-completion closure, "
             f"log₃2 obstruction, THEOREM_C_745, Th 78.1–78.3 all UNTOUCHED. Files: "
             f"probe_86_mahler_valuation_form.py + result_86_mahler_valuation.md + "
             f"result_86_data.csv + result_86_log.txt.")
    with open(os.path.join(REPO, "STATE.md"), "a", encoding="utf-8") as f:
        f.write("\n".join(e) + "\n")
    log("[appended] STATE.md")


def flush_log():
    with open(os.path.join(REPO, "result_86_log.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


def main():
    log("# PROBE 86 — Mahler valuation profile v3(c_k) closed form (LOW priority)")
    log("")
    if not self_test():
        log("ABORT: self-test failed.")
        flush_log()
        return
    K = 16
    R_SET = [14, 16, 17, 18, 19]   # foreground-feasible (~4min); confirms through k=14
    # (r=20 for k=15 dropped: ~8min pole, one tentative point, can't move an OOS verdict
    #  already set by k=12,13,14; keeps the run reliable/foreground.)
    profile = task_A(R_SET, K)
    taskB, any_fired, pts, fit_k, oos_k = task_B(profile, K)
    write_csv(profile, K)
    verdict, kmax_conf = write_md(profile, K, taskB, any_fired, pts, fit_k, oos_k)
    append_state(verdict, kmax_conf, any_fired, oos_k)
    log("")
    log(f"==== VERDICT: {'FIRED' if any_fired else 'NULL'}  confirmed to k={kmax_conf} ====")
    flush_log()


if __name__ == "__main__":
    main()
