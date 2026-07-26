# =============================================================================
# 00_setup.R
# Read the data and build all inputs for the area-level and IPP-BYM2 models.
# All case data are read from the individual patient file; population density
# from the BPS workbook; geometry (adjacency, distances) from the shapefile.
#
# Run scripts from the repository root, e.g.  Rscript R/00_setup.R
# =============================================================================
suppressMessages({
  library(readxl); library(sf); library(spdep); library(dplyr)
})

DATA_DIR <- "data"                                   # edit if your data live elsewhere
IND <- file.path(DATA_DIR, "Individual dengue data.xlsx")
WB  <- file.path(DATA_DIR, "Dengue Spatial Data 5d.xlsx")            # population density only
SHP <- file.path(DATA_DIR, "shapefile", "Makassar_no_sangkarrang.shp")

H <- 8L; J <- 14L                                    # hospitals, subdistricts

## ---- subdistrict covariates (individual file) + population density (BPS) ----
area_data <- read_excel(IND, "total case by subdistrict") |>
  transmute(sub_id = as.integer(subdistrict_id),
            Subdistrict = trimws(Subdistrict),
            cases = as.integer(`Number of cases`),
            NDVI = as.numeric(NDVI),
            health = as.numeric(health_facility_density)) |>
  left_join(read_excel(WB, 1) |>
              transmute(sub_id = as.integer(sub_id), pop_dens = as.numeric(Pop_dens)),
            by = "sub_id") |>
  arrange(sub_id)

## ---- hospital table (individual file) ----
hospital <- read_excel(IND, "total case by hospital") |>
  transmute(hosp_id = as.integer(hosp_id), name = trimws(hosp_name),
            total_cases = as.integer(total_cases),
            lat = as.numeric(lat), lon = as.numeric(lon)) |>
  arrange(hosp_id)

## ---- case-level data -> three location-resolution groups ----
cases <- read_excel(IND, "Individual data combined") |>
  transmute(hosp = as.integer(Hosp_id), sub = as.integer(Subdistrict_id),
            dist = as.numeric(Distance))

# geolocated cases (hospitals with coordinate records): exact travel distance
coord_cases <- cases |> filter(hosp <= 5, is.finite(dist), dist > 0)
# subdistrict-only cases (H6, H7): aggregate to (hospital, subdistrict) counts
region <- as.data.frame(table(hosp = cases$hosp[cases$hosp %in% c(6, 7)],
                              sub  = cases$sub [cases$hosp %in% c(6, 7)]))
region <- region[region$Freq > 0, ]
# unlocated cases (H8): total only
n_missing <- integer(H); n_missing[8] <- hospital$total_cases[8]
# full known-origin table (used by area models and origin-recovery scoring)
n_obs <- matrix(0L, H, J)
n_obs[] <- as.integer(table(factor(cases$hosp, 1:H), factor(cases$sub, 1:J)))

## ---- geometry: adjacency, BYM2 scaling, hospital->centroid distances ----
shp <- st_cast(st_make_valid(st_read(SHP, quiet = TRUE)), "MULTIPOLYGON", warn = FALSE)
W   <- listw2mat(nb2listw(poly2nb(shp, queen = TRUE), style = "B", zero.policy = TRUE))
edges <- which(W == 1, arr.ind = TRUE); edges <- edges[edges[, 1] < edges[, 2], , drop = FALSE]
node1 <- edges[, 1]; node2 <- edges[, 2]; N_edges <- nrow(edges)
Qp <- solve(diag(rowSums(W)) - W + matrix(1 / J, J, J))
scaling_k <- exp(mean(log(diag(Qp))))                # BYM2 scaling factor (Riebler et al. 2016)

map    <- shp |> left_join(area_data, by = c("NAMOBJ" = "Subdistrict")) |> arrange(sub_id)
map_utm <- st_transform(st_make_valid(map), 32750)
hosp_utm <- st_transform(st_as_sf(hospital, coords = c("lon", "lat"), crs = 4326), 32750)
dist_mat <- as.matrix(units::drop_units(
  st_distance(hosp_utm, suppressWarnings(st_centroid(map_utm))))) / 1000   # H x J, km

## ---- standardised covariates and population offset ----
NDVI_s   <- as.numeric(scale(area_data$NDVI))
health_s <- as.numeric(scale(area_data$health))
area_km2 <- as.numeric(units::set_units(st_area(map_utm), km^2))
population <- area_data$pop_dens * area_km2
log_pop  <- log(population) - mean(log(population))  # centred population offset

# convenient tier objects
t1 <- list(h = coord_cases$hosp, j = coord_cases$sub, d = coord_cases$dist)
h2 <- as.integer(as.character(region$hosp))
j2 <- as.integer(as.character(region$sub))
n2 <- as.integer(region$Freq)

message(sprintf("[setup] cases: %d coord + %d region + %d unlocated = %d | edges %d",
                length(t1$h), sum(n2), n_missing[8],
                length(t1$h) + sum(n2) + n_missing[8], N_edges))
