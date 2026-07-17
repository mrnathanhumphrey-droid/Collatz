"""
PROBE R1 -- the (e, gamma) class digraph of the real q=3 operator (handoff for the
invariant hunt). Direct methods only.

Operator: build_M_gen(3, L, 2, [lam^delta]), lam=1/2, L=2 and L=3.
Coordinates: GENERATOR CONVENTION 2^e = b*a^{-1} mod 3^L, e in Z/D (D=|<2>|); gamma in Z/3^L.
  e=0 <=> a=b (diagonal). Move type Delta_e = (delta_a - delta_b) mod D; e' = e + Delta_e.
  T = a*2^{-delta_a} - b*2^{-delta_b}; gate (gamma+T)==0 mod 3; gamma' = (gamma+T)//3.
  (T depends on the representative a, so gate-survival can vary within a class --
   the closure clean-union check below is a REAL test.)

(a) project state digraph -> (e,gamma) classes (54 nodes L=2, 486 L=3), edges tagged
    (Delta_e, gamma->gamma', weight).
(b) CLOSURE = forward closure of diag-carry classes {e=0, gamma!=0} vs SUPPORT.
    Expect 21/33 (L=2), 171/315 (L=3). Verify the state-level closure is a CLEAN UNION
    of whole classes; if not -> STOP, report loudly (gauge broken).
(c) handoff table: one row per class [e, gamma base-3, closure flag, out-edges (Delta_e, g->g')].
"""
import numpy as np

from probe_phase2a_q2b_q6 import build_M_gen, subgroup

Q = 3


def setup(L):
    qL = Q ** L
    sub = subgroup(2 % qL, qL)
    D = len(sub)
    dlog = {}
    x = 1 % qL
    for e in range(D):
        dlog[x] = e
        x = (x * 2) % qL
    states = [(a, b, g) for a in sub for b in sub for g in range(qL)]
    sidx = {s: i for i, s in enumerate(states)}
    cls = {}                                   # state -> (e, gamma)
    for (a, b, g) in states:
        e = dlog[(b * pow(a, -1, qL)) % qL]
        cls[(a, b, g)] = (e, g)
    raw = [0.5 ** d for d in range(1, D + 1)]
    w = np.array(raw) / sum(raw)
    geninv = pow(2, -1, qL)
    mult = [((pow(geninv, d, qL)), w[d - 1]) for d in range(1, D + 1)]  # (2^{-d}, w_d), d=1..D
    return qL, sub, D, states, sidx, cls, mult


def class_edges(L):
    """Return: nodes, state_adj (src->set dest), class_out[(e,g)] = set of (De, g', wt-bucket)."""
    qL, sub, D, states, sidx, cls, mult = setup(L)
    state_adj = {s: set() for s in states}
    class_out = {}                              # (e,g) -> dict{(De, g'): weight}
    for (a, b, g) in states:
        e = cls[(a, b, g)][0]
        for da in range(1, D + 1):
            ga, wa = mult[da - 1]
            ap = (a * ga) % qL
            for db in range(1, D + 1):
                gb, wb = mult[db - 1]
                bp = (b * gb) % qL
                T = (ap - bp) % qL
                if (g + T) % Q == 0:
                    gp = ((g + T) // Q) % qL
                    dest = (ap, bp, gp)
                    state_adj[(a, b, g)].add(dest)
                    De = (da - db) % D
                    key = (e, g)
                    d = class_out.setdefault(key, {})
                    d[(De, gp)] = d.get((De, gp), 0.0) + wa * wb
    return qL, D, states, cls, state_adj, class_out


def forward_closure_states(seed, state_adj):
    seen = set(seed); frontier = list(seed)
    while frontier:
        nxt = []
        for s in frontier:
            for d in state_adj[s]:
                if d not in seen:
                    seen.add(d); nxt.append(d)
        frontier = nxt
    return seen


def run(L):
    qL, D, states, cls, state_adj, class_out = class_edges(L)
    all_classes = sorted(set(cls.values()))
    nC = len(all_classes)
    # class membership
    members = {}
    for s, c in cls.items():
        members.setdefault(c, []).append(s)
    # seed = diag-carry STATES {e=0, gamma!=0}
    seed = [s for s in states if cls[s][0] == 0 and s[2] != 0]
    reached = forward_closure_states(seed, state_adj)
    reached_classes = set(cls[s] for s in reached)
    # CLEAN-UNION check: every reached class fully reached; no partial
    partial = []
    full = []
    for c in reached_classes:
        if all(s in reached for s in members[c]):
            full.append(c)
        else:
            partial.append(c)
    clean = (len(partial) == 0)
    closure = sorted(reached_classes)
    support = sorted(set(all_classes) - reached_classes)
    return dict(L=L, D=D, qL=qL, nClasses=nC, nClosure=len(closure), nSupport=len(support),
                clean=clean, nPartial=len(partial), closure=closure, support=support,
                full=set(full), partial=set(partial),
                class_out=class_out, all_classes=all_classes)


def base3(g, ndig):
    d = []
    for _ in range(ndig):
        d.append(g % 3); g //= 3
    return ''.join(str(x) for x in reversed(d))


def emit_table(res, path):
    L = res['L']; ndig = L
    full = res['full']; partial = res['partial']

    def flag(c):
        return 'F' if c in full else ('P' if c in partial else '0')   # F=full-closure P=partial(gauge-split) 0=support
    rows = sorted(res['all_classes'], key=lambda c: ({'F': 0, 'P': 1, '0': 2}[flag(c)], c[0], c[1]))
    lines = [f"# (e,gamma) class table  L={L}  D={res['D']}  classes={res['nClasses']}  "
             f"reached_classes={res['nClosure']} (full={len(full)} PARTIAL={len(partial)}) support={res['nSupport']}  "
             f"clean_union={res['clean']}  [convention 2^e=b*a^-1 mod 3^L; gamma UNTRANSFORMED (modular gauge)]",
             "# flag: F=full-closure  P=PARTIAL (gauge-split, closure NOT clean)  0=support",
             "# columns: e | gamma(base3) | flag | out_edges: (De,g->g') ..."]
    for c in rows:
        e, g = c
        outs = res['class_out'].get(c, {})
        edstr = ' '.join(f"({De},{g}->{gp})" for (De, gp) in sorted(outs.keys()))
        lines.append(f"{e}\t{base3(g,ndig)}\t{flag(c)}\t{edstr}")
    with open(path, "w", encoding="utf-8") as f:
        f.write('\n'.join(lines) + '\n')


def main():
    print("# PROBE R1 -- (e,gamma) class digraph; generator convention 2^e = b*a^{-1} mod 3^L\n")
    for L in [2, 3]:
        r = run(L)
        exp = {2: (21, 33), 3: (171, 315)}[L]
        ok = (r['nClosure'], r['nSupport']) == exp
        print(f"## L={L}: {r['nClasses']} classes (D={r['D']}, gamma in Z/{r['qL']}). "
              f"CLOSURE={r['nClosure']} / SUPPORT={r['nSupport']}  (expected {exp[0]}/{exp[1]}) "
              f"[{'MATCH' if ok else 'MISMATCH'}]")
        if not r['clean']:
            print(f"   *** STOP: closure is NOT a clean union of classes -- {r['nPartial']} partial "
                  f"classes (gauge broken). Skeleton pauses. ***")
        else:
            print(f"   clean union of whole classes: YES ({r['nClosure']} closure classes, "
                  f"all members reached; 0 partial)")
        path = f"outputs/class_table_L{L}.tsv"
        emit_table(r, path)
        print(f"   handoff table -> {path}\n")


if __name__ == "__main__":
    main()
