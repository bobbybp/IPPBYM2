// ============================================================
// RECOMMENDED MODEL — IPP-BYM2, three-tier, per-capita, power-law kernel
//   log lambda_j = log(pop_j) + alpha + b_NDVI NDVI + b_health health + phi_j
//   catchment kernel: k_hj = (1 + gamma_h * d_hj)^(-b)     [gravity / power-law]
//   gamma_h ~ lognormal(mu_g, sig_g)  (hierarchical; NO decay-covariate regression)
//   b        shared tail exponent
//   Three tiers: exact (tier1), region-integrated-at-centroid (tier2), missing (tier3)
// ============================================================
data {
  int<lower=1> J; int<lower=1> H;
  int<lower=0> N1; array[N1] int<lower=1,upper=H> h1; array[N1] int<lower=1,upper=J> j1; vector<lower=0>[N1] d1;
  int<lower=0> N2; array[N2] int<lower=1,upper=H> h2; array[N2] int<lower=1,upper=J> j2; array[N2] int<lower=0> n2;
  array[H] int<lower=0> nmis;
  matrix[H,J] dist_mat;
  vector[J] NDVI; vector[J] health; vector[J] log_pop;
  int<lower=0> N_edges; array[N_edges] int<lower=1,upper=J> node1; array[N_edges] int<lower=1,upper=J> node2;
  real<lower=0> scaling_k;
  vector<lower=0,upper=1>[H] w_hosp;
}
parameters {
  real alpha; real beta_NDVI; real beta_health;
  real mu_g; real<lower=0> sig_g;
  vector<lower=0>[H] gamma_h;
  real<lower=0> b;
  real<lower=0> sigma_total; real<lower=0,upper=1> rho;
  vector[J] theta; vector[J] epsilon;
}
transformed parameters {
  vector[J] phi = (sqrt(1-rho)*epsilon + sqrt(rho/scaling_k)*theta)*sigma_total;
  vector[J] log_lambda = log_pop + alpha + beta_NDVI*NDVI + beta_health*health + phi;
  matrix[H,J] logk;
  for (h in 1:H) for (j in 1:J) logk[h,j] = -b * log1p(gamma_h[h]*dist_mat[h,j]);
}
model {
  alpha ~ normal(-6,3); beta_NDVI ~ normal(0,1); beta_health ~ normal(0,1);   // intercept ~ log per-capita rate
  mu_g ~ normal(0,1); sig_g ~ normal(0,1);
  gamma_h ~ lognormal(mu_g, sig_g);
  b ~ lognormal(log(1.5), 0.5);
  sigma_total ~ normal(0,1); rho ~ beta(0.5,0.5);
  epsilon ~ std_normal(); sum(epsilon) ~ normal(0,0.001*J);
  target += -0.5*dot_self(theta[node1]-theta[node2]); sum(theta) ~ normal(0,0.001*J);

  for (h in 1:H) for (j in 1:J) target += -w_hosp[h]*exp(log_lambda[j] + logk[h,j]);          // compensator
  for (i in 1:N1) target += w_hosp[h1[i]]*(log_lambda[j1[i]] - b*log1p(gamma_h[h1[i]]*d1[i]));
  for (m in 1:N2) target += w_hosp[h2[m]]*n2[m]*(log_lambda[j2[m]] + logk[h2[m],j2[m]]);
  for (h in 1:H) if (nmis[h]>0) {                                                     // tier 3 missing
    vector[J] lp; for (j in 1:J) lp[j] = log_lambda[j] + logk[h,j];
    target += w_hosp[h]*nmis[h]*log_sum_exp(lp);
  }
}
generated quantities {
  vector[J] lambda = exp(log_lambda);                 // expected cases (with pop)
  vector[J] log_rr = alpha + beta_NDVI*NDVI + beta_health*health + phi;  // log per-capita relative risk
  vector[H] hosp_expected;
  vector[H] catch_radius;                              // km at which kernel = 0.5
  matrix[H,J] origin_prob;
  for (h in 1:H) {
    vector[J] lp; for (j in 1:J) lp[j] = log_lambda[j] + logk[h,j];
    hosp_expected[h] = sum(exp(lp));
    origin_prob[h,] = to_row_vector(softmax(lp));
    catch_radius[h] = (pow(2.0, 1.0/b) - 1.0) / gamma_h[h];
  }
}
