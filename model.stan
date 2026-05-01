// Hierarchical Bayesian model for Collatz total stopping time sigma(n),
// grouped by residue class mod 2^k (k configurable, default 10 -> 512 odd
// classes). Bulk likelihood is Normal; the GPD tail layer is fit
// post-hoc in Python on the residuals from this fit (see diagnose.py).
//
// Non-centered parameterization for class-level alpha, beta. Vectorized
// + reduce_sum likelihood for parallel evaluation across threads.

functions {
  real partial_sum_lpdf(array[] real y_slice,
                        int start, int end,
                        array[] int class_idx,
                        vector log_n,
                        vector alpha,
                        vector beta,
                        real phi) {
    real lp = 0.0;
    int n = end - start + 1;
    vector[n] mu;
    for (i in 1:n) {
      int j = start + i - 1;
      mu[i] = alpha[class_idx[j]] + beta[class_idx[j]] * log_n[j];
    }
    return normal_lpdf(y_slice | mu, phi);
  }
}

data {
  int<lower=1> N;                      // observations
  int<lower=1> K;                      // residue classes
  array[N] int<lower=1, upper=K> class_idx;
  vector[N] log_n;                     // log(n)
  array[N] real y;                     // sigma(n) as real for reduce_sum
  int<lower=1> grainsize;              // reduce_sum grainsize
}

parameters {
  // Hyperparameters
  real mu_alpha;
  real mu_beta;
  real<lower=0> tau_alpha;
  real<lower=0> tau_beta;

  // Class-level offsets (non-centered)
  vector[K] alpha_z;
  vector[K] beta_z;

  // Bulk residual scale (global; class-specific can come later)
  real<lower=0> phi;
}

transformed parameters {
  vector[K] alpha = mu_alpha + tau_alpha * alpha_z;
  vector[K] beta  = mu_beta  + tau_beta  * beta_z;
}

model {
  // Hyperpriors. Heuristic for ODD-n filter: 3 / (ln 4 - ln 3) = 10.43
  mu_alpha    ~ normal(0.0, 5.0);
  mu_beta     ~ normal(10.4, 1.5);
  tau_alpha   ~ normal(0.0, 3.0);
  tau_beta    ~ normal(0.0, 0.3);
  phi         ~ normal(0.0, 100.0);  // empirical residual std ~64

  // Class effects
  alpha_z     ~ std_normal();
  beta_z      ~ std_normal();

  // Likelihood, parallelized via reduce_sum
  target += reduce_sum(partial_sum_lpdf, y, grainsize,
                       class_idx, log_n, alpha, beta, phi);
}

generated quantities {
  // Per-draw class-level slope/intercept are saved automatically (alpha, beta).
  // Posterior predictive tail probabilities are computed in Python from
  // these draws (see diagnose.py) so the Stan output stays small.
}
