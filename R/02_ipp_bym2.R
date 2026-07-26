# =============================================================================
# 02_ipp_bym2.R
# Fit the full IPP-BYM2 model: population-adjusted BYM2 intensity, gravity
# catchment kernel, and three location-resolution groups (exact / subdistrict /
# unlocated) in a single likelihood. Produces the parameter estimates, the
# per-person relative-risk surface and exceedance probabilities, hospital
# catchment summaries, calibration, and the inferred origins of unlocated cases.
# =============================================================================
source("R/00_setup.R")
suppressMessages({ library(rstan); library(dplyr) })
rstan_options(auto_write = TRUE); options(mc.cores = 4)
dir.create("output", showWarnings = FALSE)

stan_data <- list(
  J = J, H = H,
  N1 = length(t1$h), h1 = t1$h, j1 = t1$j, d1 = t1$d,
  N2 = length(n2), h2 = h2, j2 = j2, n2 = n2, nmis = n_missing,
  dist_mat = dist_mat, NDVI = NDVI_s, health = health_s, log_pop = log_pop,
  N_edges = N_edges, node1 = node1, node2 = node2, scaling_k = scaling_k)

fit <- sampling(stan_model("stan/ipp_bym2.stan"), data = stan_data,
                chains = 4, cores = 4, iter = 3000, warmup = 1000, seed = 1234,
                refresh = 100, control = list(adapt_delta = 0.97, max_treedepth = 13))
saveRDS(fit, "output/fit_ipp_bym2.rds")
print(check_hmc_diagnostics(fit))

## ---- parameter estimates (Table 4) ----
pars <- c("beta_NDVI", "beta_health", "b", "sigma_total", "rho")
print(round(summary(fit, pars = pars, probs = c(0.05, 0.95))$summary[, c("mean","5%","95%","Rhat")], 3))

## ---- per-person relative risk + exceedance probability (Table 5 / Fig) ----
RR <- exp(as.matrix(fit, pars = "log_rr")); RR <- RR / rowMeans(RR)   # mean-1 normalised
risk <- data.frame(Subdistrict = area_data$Subdistrict,
                   RR = round(colMeans(RR), 2),
                   RR_lo = round(apply(RR, 2, quantile, 0.05), 2),
                   RR_hi = round(apply(RR, 2, quantile, 0.95), 2),
                   EP = round(colMeans(RR > 1), 3))
risk <- risk[order(-risk$RR), ]
cat("\n### Per-person relative risk and exceedance probability ###\n"); print(risk, row.names = FALSE)
write.csv(risk, "output/relative_risk.csv", row.names = FALSE)

## ---- hospital catchment (mean distance of inferred origins) ----
op <- as.matrix(fit, pars = "origin_prob")
catch <- sapply(1:H, function(h) {
  p <- colMeans(op[, grep(sprintf("^origin_prob\\[%d,", h), colnames(op))]); p <- p / sum(p)
  sum(p * dist_mat[h, ])
})
cat("\n### Hospital mean catchment distance (km) ###\n")
print(data.frame(hospital = hospital$name, mean_km = round(catch, 2)))

## ---- calibration and inferred origins of the unlocated cases ----
cat("\n### Hospital calibration (observed vs predicted totals) ###\n")
print(data.frame(hospital = hospital$name, observed = hospital$total_cases,
                 predicted = round(colMeans(as.matrix(fit, pars = "hosp_expected")), 1)))
p8 <- colMeans(op[, grep("^origin_prob\\[8,", colnames(op))])
cat("\n### Inferred origins of the unlocated cases (expected count per subdistrict) ###\n")
print(data.frame(Subdistrict = area_data$Subdistrict,
                 expected = round(n_missing[8] * p8, 1))[order(-p8), ], row.names = FALSE)
