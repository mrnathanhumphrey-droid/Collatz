#!/usr/bin/env bash
# Tomorrow morning batch: q=5, k=8, N=10^7 partial correlation analysis.
#
# At k=6, odd_steps_in_prefix and prefix_total_steps and log(a_final) are
# 100% collinear (even_steps_in_prefix is forced to k by construction).
# At k=8, the prefix algorithm allows non-trivial joint variation in
# (odd_steps, even_steps); even_steps is NOT forced constant for all
# residue classes once we move beyond minimal prefix. This run tests
# whether log(a_final), odd_steps, and even_steps survive joint partial
# correlation — three possible outcomes:
#   1. odd + even both significant: two-axis structural finding
#   2. only odd significant: clean one-axis (odd-step count) result
#   3. log(a_final) significant but odd alone is not: multiplicative
#      product effect not captured by count alone
#
# Also reports c_final partial correlation, which trended +0.32 at k=6
# (potential second axis: multiplicative growth hurts, additive offset helps).

set -e
cd /c/Collatz

# Generate q=5, N=10^7 data — already done if file exists
if [ ! -f data/q_main_q5_N10000000.parquet ]; then
  echo "[gen] q=5 N=10^7 data missing; generating ..."
  python generate_q.py --q 5 --N 10000000 --max_steps 100000 --max_value 1e18
fi

# Run partial correlation at k=8 (256 odd residue classes mod 512)
echo "[run] partial correlation, q=5, k=8, N=10^7"
python experiments/10b_q_partial_correlation.py --q 5 --N 10000000 --k 8 \
  2>&1 | tee experiments_output/11_q5_k8_partial.log

echo
echo "[done] log saved to experiments_output/11_q5_k8_partial.log"
echo
echo "Look for:"
echo "  1. Does even_steps_in_prefix have non-zero variance? (collinearity collapse)"
echo "  2. Multivariate regression: which of odd_steps / even_steps / log_a_final / log_c_final"
echo "     survive partial correlation at p<0.05?"
echo "  3. R^2 of full model vs the k=6 single-predictor R^2 = 0.87"
