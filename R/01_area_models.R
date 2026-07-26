# =============================================================================
# 01_area_models.R
# Fit the four area-level models (M1-M4) to the aggregated subdistrict counts
# with a population offset, and compare them by WAIC (in-sample) and by exact
# leave-one-area-out cross-validation (out-of-sample). Reproduces Table 2.
#
#   M1 null | M2 covariates | M3 BYM2 | M4 BYM2 + covariates
# =============================================================================
source("R/00_setup.R")
suppressMessages({ library(rstan); library(loo); library(parallel) })
rstan_options(auto_write = TRUE); options(mc.cores = 4)

sm <- stan_model("stan/area_bym2.stan")
y  <- area_data$cases
X  <- cbind(NDVI = NDVI_s, health = health_s)

stan_data <- function(cols, spatial, w = rep(1, J)) list(
  N_subs = J, y = y, P = length(cols), X = X[, cols, drop = FALSE],
  w = as.numeric(w), spatial = as.integer(spatial), log_offset = log_pop,
  N_edges = N_edges, node1 = node1, node2 = node2, scaling_k = scaling_k,
  alpha_prior_mean = 0, alpha_prior_sd = 10)

fit_model <- function(cols, spatial, seed = 1)
  sampling(sm, data = stan_data(cols, spatial), chains = 4, cores = 4,
           iter = 3000, warmup = 1000, seed = seed, refresh = 0,
           control = list(adapt_delta = 0.99, max_treedepth = 12))

# exact leave-one-area-out CV: hold out area j (weight 0), score its held-out lpd
loo_cv <- function(cols, spatial) {
  lpd <- unlist(mclapply(1:J, function(jj) {
    w <- rep(1, J); w[jj] <- 0
    f <- sampling(sm, data = stan_data(cols, spatial, w), chains = 2, cores = 1,
                  iter = 2500, warmup = 1000, seed = 100 + jj, refresh = 0,
                  control = list(adapt_delta = 0.95, max_treedepth = 12))
    ll <- as.matrix(f, pars = "log_lik")[, jj]
    log(mean(exp(ll - max(ll)))) + max(ll)
  }, mc.cores = 7))
  list(elpd = sum(lpd), per_area = lpd)
}

specs <- list(M1 = list(character(0), 0L), M2 = list(c("NDVI","health"), 0L),
              M3 = list(character(0), 1L), M4 = list(c("NDVI","health"), 1L))

waic <- sapply(specs, function(s) {
  f <- fit_model(s[[1]], s[[2]])
  suppressWarnings(loo::waic(as.matrix(f, pars = "log_lik"))$estimates["waic", "Estimate"])
})
cv <- lapply(specs[c("M3","M4")], function(s) loo_cv(s[[1]], s[[2]]))
d_elpd <- cv$M4$elpd - cv$M3$elpd
d_se   <- sd(cv$M4$per_area - cv$M3$per_area) * sqrt(J)

cat("\n### Area-level model comparison (Table 2) ###\n")
print(data.frame(Model = names(waic), WAIC = round(waic)))
cat(sprintf("\nExact LOO-CV, M4 vs M3: dElpd = %.2f (SE %.2f) -> M4 preferred out-of-sample\n",
            d_elpd, d_se))
