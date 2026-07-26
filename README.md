# IPP-BYM2: mapping dengue risk from incomplete hospital surveillance data

Code accompanying the paper *"Mapping dengue risk from incomplete hospital
surveillance data in Makassar, Indonesia."*

The method is a hierarchical Bayesian inhomogeneous Poisson process (IPP) with a
population-adjusted Besag–York–Mollié (BYM2) intensity and a distance-based
gravity catchment kernel. It maps dengue risk per person while modelling hospital
catchment behaviour and inferring the residential origins of cases that are
recorded at coarse or missing spatial resolution.

## Repository layout

```
R/
  00_setup.R            Read data and build all model inputs
  01_area_models.R      Area-level models M1-M4 (WAIC + exact LOO-CV)  -> Table 2
  02_ipp_bym2.R         Full IPP-BYM2 model: parameters, per-person risk map,
                        exceedance probabilities, catchment, calibration,
                        inferred origins                                -> Tables 4-5, Figs
  03_origin_recovery.R  Held-out cross-validation of origin recovery    -> Table 3
stan/
  area_bym2.stan        Area-level BYM2 model (covariates/spatial toggles) for M1-M4
  ipp_bym2.stan         Full IPP-BYM2 model
  ipp_bym2_loho.stan    IPP-BYM2 with a hospital-hold-out weight (for origin recovery)
data/                   Data are NOT included (restricted); see data/README.md
```

## Requirements

R (>= 4.2) with `rstan`, `loo`, `sf`, `spdep`, `readxl`, `dplyr`, and `units`.
Models are written in Stan and compiled through `rstan`.

## Reproducing the analysis

1. Obtain the data (restricted; see `data/README.md`) and place the files in
   `data/`.
2. From the repository root, run in order:

   ```r
   Rscript R/01_area_models.R      # area-level model comparison
   Rscript R/02_ipp_bym2.R         # full model + risk map + origins
   Rscript R/03_origin_recovery.R  # held-out origin-recovery cross-validation
   ```

   Each script sources `R/00_setup.R`, which builds the model inputs from the
   data in `data/`.

## Model summary

For subdistrict *j* with population *P_j* and covariates **X**_j, the per-person
log relative risk uses a BYM2 random effect,

```
log lambda_j = log P_j + alpha + X_j' beta + phi_j .
```

A case resident in *j* is observed at hospital *h* with a gravity kernel
`k(d) = (1 + gamma_h d)^(-b)`, giving expected cases `mu_hj = lambda_j k(d_hj)`.
The IPP log-likelihood sums log-intensities over observed cases minus the
compensator `sum_h sum_j mu_hj`. Cases are handled at three resolutions:
geolocated (exact distance), subdistrict-only, and unlocated (origin
marginalised over subdistricts). See the paper for full details.

## Data availability

Individual case data are not redistributed owing to ethical-approval and
data-sharing restrictions. See the paper's *Data availability* statement.
