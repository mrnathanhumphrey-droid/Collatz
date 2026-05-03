"""Phase 1 inverse Collatz tree build — depth 30, parquet + checkpoint MD."""
import sys
import time
import math
import importlib.util
import random
from pathlib import Path

import numpy as np
import polars as pl

# Reuse existing prefix-decomposition utility
spec = importlib.util.spec_from_file_location(
    "alpha_decomp",
    r"C:\Collatz\experiments\01_alpha_decomposition.py",
)
alpha_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(alpha_mod)
deterministic_prefix = alpha_mod.deterministic_prefix
LOG_FACTOR_ODD = alpha_mod.LOG_FACTOR_ODD

OUT_DIR = Path(r"C:\Collatz\inverse_tree")
OUT_DIR.mkdir(parents=True, exist_ok=True)
DEPTH_MAX = 50

def is_eligible(n):
    if n % 3 != 1:
        return False
    m = (n - 1) // 3
    return m > 1 and m % 2 == 1

def forward_T(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1

def forward_sigma(n):
    s = 0
    while n != 1:
        n = forward_T(n)
        s += 1
    return s

def make_prefix_cache(k):
    M = 1 << k
    cache = {}
    def get(r):
        if r not in cache:
            steps, a_f, c_f = deterministic_prefix(r, M)
            alpha = steps + LOG_FACTOR_ODD * math.log(a_f / M) if a_f > 0 else 0.0
            cache[r] = (steps, a_f, c_f, alpha)
        return cache[r]
    return get, M

# BFS from root
t0 = time.perf_counter()
layers = [[1]]
seen = {1}
for d in range(DEPTH_MAX):
    nxt = []
    for n in layers[d]:
        c1 = 2 * n
        if c1 not in seen:
            seen.add(c1); nxt.append(c1)
        if is_eligible(n):
            c2 = (n - 1) // 3
            if c2 not in seen:
                seen.add(c2); nxt.append(c2)
    layers.append(nxt)
build_time = time.perf_counter() - t0

# OEIS A005186 cross-check d=0..20
OEIS = [1,1,1,1,1,2,2,4,4,6,6,8,10,14,18,24,29,36,44,58,72,91,113,143,179,227,287,366,460,578,732]
oeis_pass = True
oeis_mm = None
for d in range(min(len(OEIS), DEPTH_MAX + 1)):
    if len(layers[d]) != OEIS[d]:
        oeis_pass = False
        oeis_mm = (d, len(layers[d]), OEIS[d])
        break

# Build full table with prefix decomposition for k in {6,8,10}
caches = {k: make_prefix_cache(k) for k in (6, 8, 10)}
rows = []
for d in range(DEPTH_MAX + 1):
    for n in layers[d]:
        parent = None if n == 1 else forward_T(n)
        cl = 2 * n
        cr = ((n - 1) // 3) if is_eligible(n) else None
        row = {"n": n, "depth": d, "parent": parent, "child_left": cl, "child_right": cr}
        for k in (6, 8, 10):
            get_p, M = caches[k]
            r = n % M
            steps, a_f, c_f, alpha = get_p(r)
            row[f"r_{k}"] = r
            row[f"ell_{k}"] = steps
            row[f"a_final_{k}"] = a_f
            row[f"c_final_{k}"] = c_f
            row[f"alpha_det_{k}"] = alpha
        rows.append(row)

df = pl.DataFrame(rows)
parquet_path = OUT_DIR / "tree_d50.parquet"
df.write_parquet(parquet_path, compression="zstd")
parquet_size = parquet_path.stat().st_size

# Forward-sigma verification on 10^4 random sample
random.seed(42)
sample_size = min(10000, len(rows))
sample = random.sample(rows, sample_size)
sigma_pass = True
sigma_mm = None
for row in sample:
    s = forward_sigma(row["n"])
    if s != row["depth"]:
        sigma_pass = False
        sigma_mm = (row["n"], s, row["depth"])
        break

# Growth-rate OLS on log(layer_count) over d=10..30
ds = np.arange(10, DEPTH_MAX + 1)
ys = np.log(np.array([len(layers[d]) for d in ds], dtype=float))
slope, intercept = np.polyfit(ds, ys, 1)
yhat = intercept + slope * ds
ss_res = float(np.sum((ys - yhat) ** 2))
ss_tot = float(np.sum((ys - ys.mean()) ** 2))
r_sq = 1 - ss_res / ss_tot

total_nodes = sum(len(L) for L in layers)

# Schema lines
schema_lines = [f"  {c}: {t}" for c, t in zip(df.columns, df.dtypes)]

md_lines = [
    "# Phase 1 Checkpoint — Inverse Collatz Tree (d ≤ 50)",
    "",
    "## Wall time",
    f"Tree construction: {build_time:.3f}s",
    "",
    "## Node counts",
    f"Total nodes (d=0..{DEPTH_MAX}): **{total_nodes:,}**",
    "",
    "| depth | n_nodes |",
    "|------:|--------:|",
    *[f"| {d} | {len(layers[d]):,} |" for d in (0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)],
    "",
    f"## Empirical growth rate (window d=10..{DEPTH_MAX})",
    f"- slope = {slope:.4f}",
    f"- intercept = {intercept:.4f}",
    f"- R² = {r_sq:.6f}",
    f"- asymptote log(4/3) ≈ 0.2877 (finite-window slope expected lower)",
    "",
    "## OEIS A005186 cross-check (d=0..20)",
    f"Result: **{'PASS' if oeis_pass else 'FAIL'}**",
]
if not oeis_pass:
    d_mm, got, want = oeis_mm
    md_lines.append(f"First mismatch at d={d_mm}: got {got}, OEIS expects {want}")
md_lines += [
    "",
    "## Forward-σ verification",
    f"Sample size: {sample_size:,}",
    f"Result: **{'PASS' if sigma_pass else 'FAIL'}**",
]
if not sigma_pass:
    n_mm, got, want = sigma_mm
    md_lines.append(f"First mismatch: n={n_mm}, forward σ={got}, stored depth={want}")
md_lines += [
    "",
    "## Parquet output",
    "Path: `C:/Collatz/inverse_tree/tree_d50.parquet`",
    f"Size: {parquet_size/1024:.1f} KiB ({parquet_size:,} bytes)",
    "Schema:",
    "```",
    *schema_lines,
    "```",
    "",
    "## Surprises / anomalies",
    "(none flagged at this checkpoint)",
]

md_path = OUT_DIR / "phase1_checkpoint.md"
md_path.write_text("\n".join(md_lines), encoding="utf-8")

print(f"[ok] wall={build_time:.3f}s total={total_nodes} OEIS={oeis_pass} sigma={sigma_pass} parquet={parquet_size/1024:.1f}KiB")
print(f"[wrote] {md_path}")
