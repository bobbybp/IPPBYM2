// ============================================================
// REVIEW CODE: flexible area-level BYM2 for honest model
// selection between M3 (spatial only) and M4 (spatial+covariates),
// and for exact leave-one-area-out cross-validation.
//
//   - X          : N_subs x P covariate matrix (P may be 0)
//   - w[j]       : likelihood weight (1 = used, 0 = held out for CV)
//   - spatial    : 1 = include BYM2 field, 0 = no spatial effect
//   - log_offset : 0 for the paper's count model, log(pop) for a rate model
//
// generated quantities returns log_lik[j] = poisson_log_lpmf(y[j]|.)
// for EVERY area, so a held-out area's log predictive density can be
// scored even though it was excluded from the likelihood.
// ============================================================
data {
  int<lower=1> N_subs;
  array[N_subs] int<lower=0> y;

  int<lower=0> P;
  matrix[N_subs, P] X;

  vector<lower=0, upper=1>[N_subs] w;
  int<lower=0, upper=1> spatial;
  vector[N_subs] log_offset;

  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=N_subs> node1;
  array[N_edges] int<lower=1, upper=N_subs> node2;
  real<lower=0> scaling_k;

  real alpha_prior_mean;
  real<lower=0> alpha_prior_sd;
}
parameters {
  real alpha;
  vector[P] beta;
  real<lower=0> sigma_total;
  real<lower=0, upper=1> rho;
  vector[N_subs] theta;
  vector[N_subs] epsilon;
}
transformed parameters {
  vector[N_subs] phi;
  vector[N_subs] log_lambda;

  if (spatial == 1)
    phi = (sqrt(1 - rho) * epsilon + sqrt(rho / scaling_k) * theta) * sigma_total;
  else
    phi = rep_vector(0.0, N_subs);

  log_lambda = alpha + log_offset + phi;
  if (P > 0) log_lambda += X * beta;
}
model {
  alpha ~ normal(alpha_prior_mean, alpha_prior_sd);
  if (P > 0) beta ~ normal(0, 1);

  if (spatial == 1) {
    sigma_total ~ normal(0, 1);
    rho ~ beta(0.5, 0.5);
    epsilon ~ std_normal();
    target += -0.5 * dot_self(theta[node1] - theta[node2]);   // pure ICAR (no double prior)
    sum(theta) ~ normal(0, 0.001 * N_subs);
  } else {
    // keep parameters proper but inert
    sigma_total ~ normal(0, 1);
    rho ~ beta(0.5, 0.5);
    epsilon ~ std_normal();
    theta ~ std_normal();
  }

  // weighted Poisson likelihood (held-out areas get weight 0)
  for (j in 1:N_subs)
    target += w[j] * poisson_log_lpmf(y[j] | log_lambda[j]);
}
generated quantities {
  vector[N_subs] log_lik;
  for (j in 1:N_subs)
    log_lik[j] = poisson_log_lpmf(y[j] | log_lambda[j]);
}
