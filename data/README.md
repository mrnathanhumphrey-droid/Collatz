# data/

Parquet files are git-ignored (regeneratable, ~1 GB total).

To recreate the data used by experiments and the writeup:

```bash
python ../generate.py --N 10000000 --no-vseq      # 62 MB,  for 04_head_to_head_nb_glm.py
python ../generate.py --N 33554432 --no-vseq      # 200 MB, for k=6 alpha decomposition + diagnose.py
python ../generate.py --N 134217728 --no-vseq     # 786 MB, for k-sweep + N-scaling + c_final tests
```

Each takes a few seconds with the numba-memoized cache pass.

Smaller N values (2^20, 2^22, 2^23, 2^24) are also referenced in 03_n_scaling.py:

```bash
python ../generate.py --N 1048576 --no-vseq
python ../generate.py --N 4194304 --no-vseq
python ../generate.py --N 8388608 --no-vseq
python ../generate.py --N 16777216 --no-vseq
```
