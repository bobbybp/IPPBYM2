# =============================================================================
# 03_origin_recovery.R
# Held-out cross-validation of case-origin recovery (Table 3). Each hospital
# with recorded origins is withheld in turn; the model is refitted without its
# cases (so its catchment scale comes from the prior); and the inferred origin
# distribution is compared with the truth by total variation distance (TVD) and
# multinomial log predictive score. Benchmarked against uniform and
# population-proportional allocation.
# =============================================================================
source("R/00_setup.R")
suppressMessages({ library(rstan) })
rstan_options(auto_write = TRUE); options(mc.cores = 4)

sm <- stan_model("stan/ipp_bym2_loho.stan")   # IPP-BYM2 with a hospital weight vector
area_km2 <- as.numeric(units::set_units(sf::st_area(sf::st_transform(sf::st_make_valid(map), 32750)), km^2))
population <- area_data$pop_dens * area_km2

stan_data <- function(w) list(
  J = J, H = H, N1 = length(t1$h), h1 = t1$h, j1 = t1$j, d1 = t1$d,
  N2 = length(n2), h2 = h2, j2 = j2, n2 = n2, nmis = n_missing,
  dist_mat = dist_mat, NDVI = NDVI_s, health = health_s, log_pop = log_pop,
  N_edges = N_edges, node1 = node1, node2 = node2, scaling_k = scaling_k, w_hosp = as.numeric(w))

tvd <- function(truth, p) 0.5 * sum(abs(truth / sum(truth) - p))
lsc <- function(truth, p) sum(truth * log(pmax(p, 1e-10)))
p_pop <- population / sum(population); p_unif <- rep(1 / J, J)

res <- data.frame()
for (hout in 1:7) {                                 # 7 hospitals with recorded origins
  w <- rep(1, H); w[hout] <- 0
  f <- sampling(sm, data = stan_data(w), chains = 4, cores = 4, iter = 2000, warmup = 1000,
                seed = 500 + hout, refresh = 0, control = list(adapt_delta = 0.9, max_treedepth = 12))
  pr <- colMeans(as.matrix(f, pars = "origin_prob")[, grep(sprintf("^origin_prob\\[%d,", hout), colnames(as.matrix(f, pars="origin_prob")))])
  pr <- pr / sum(pr); truth <- n_obs[hout, ]
  res <- rbind(res, data.frame(hospital = hospital$name[hout], n = sum(truth),
    tvd_model = tvd(truth, pr), tvd_pop = tvd(truth, p_pop), tvd_unif = tvd(truth, p_unif),
    ls_model = lsc(truth, pr),  ls_pop = lsc(truth, p_pop),  ls_unif = lsc(truth, p_unif)))
}
w <- res$n
cat("\n### Held-out origin recovery (Table 3) ###\n")
cat(sprintf("Case-weighted mean TVD  model %.3f | pop-proportional %.3f | uniform %.3f\n",
            sum(res$tvd_model * w) / sum(w), sum(res$tvd_pop * w) / sum(w), sum(res$tvd_unif * w) / sum(w)))
cat(sprintf("Pooled log-score        model %.0f | pop-proportional %.0f | uniform %.0f\n",
            sum(res$ls_model), sum(res$ls_pop), sum(res$ls_unif)))
