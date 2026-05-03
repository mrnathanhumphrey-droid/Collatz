// Bayesian NB2-GLM with optional hierarchical random effects.
// Used by experiments/06_bb_replication.py for the apples-to-apples
// comparison against Bonacorsi-Bordoni's framework.
//
// Specification:
//   y_i | mu_i, phi  ~  NegBinomial2(mu_i, phi)
//   log(mu_i)        =  X_i^T beta + (use_re ? u_{re_idx_i} : 0)
//   u_j              ~  Normal(0, sigma_u)        // hierarchical
//   beta             ~  Normal(0, 5)
//   phi              ~  Gamma(2, 0.1)              // weakly informative
//   sigma_u          ~  HalfNormal(0, 2)
//
// Posterior predictive log score on a held-out set is computed in
// generated quantities. NB2 parameterization: mean = mu, variance =
// mu + mu^2 / phi.

functions {
  real partial_sum_lpmf(array[] int y_slice,
                        int start, int end,
                        matrix X,
                        array[] int re_idx,
                        vector beta,
                        vector u,
                        real phi,
                        int use_re) {
    int n = end - start + 1;
    vector[n] log_mu;
    for (i in 1:n) {
      int j = start + i - 1;
      log_mu[i] = X[j] * beta;
      if (use_re == 1) log_mu[i] += u[re_idx[j]];
    }
    return neg_binomial_2_log_lpmf(y_slice | log_mu, phi);
  }
}

data {
  int<lower=1> N_train;
  int<lower=1> N_test;
  int<lower=1> P;                              // number of fixed-effect columns (incl. intercept)
  matrix[N_train, P] X_train;
  matrix[N_test, P] X_test;
  array[N_train] int<lower=0> y_train;
  array[N_test] int<lower=0> y_test;
  int<lower=0, upper=1> use_re;                // 1 = include hierarchical random effect
  int<lower=1> N_re_levels;                    // number of RE levels (set to 1 if use_re=0)
  array[N_train] int<lower=1> re_idx_train;
  array[N_test] int<lower=1> re_idx_test;
  int<lower=1> grainsize;
}

parameters {
  vector[P] beta;
  real<lower=0.1> phi;
  vector[N_re_levels] u_raw;
  real<lower=0> sigma_u;
}

transformed parameters {
  vector[N_re_levels] u = sigma_u * u_raw;     // non-centered hierarchical
}

model {
  beta ~ normal(0, 5);
  phi ~ lognormal(log(5), 1);  // mass concentrated 1-25, no zero-boundary pathology
  u_raw ~ std_normal();
  sigma_u ~ normal(0, 2);

  target += reduce_sum(partial_sum_lpmf, y_train, grainsize,
                       X_train, re_idx_train, beta, u, phi, use_re);
}

generated quantities {
  // Posterior predictive log score on held-out test set
  vector[N_test] log_lik_test;
  for (i in 1:N_test) {
    real log_mu = X_test[i] * beta;
    if (use_re == 1) log_mu += u[re_idx_test[i]];
    log_lik_test[i] = neg_binomial_2_log_lpmf(y_test[i] | log_mu, phi);
  }
}
