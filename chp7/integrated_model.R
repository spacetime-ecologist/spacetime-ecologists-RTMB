library(sf)
library(viridisLite)
library(rnaturalearth)
library(fmesher)
library(terra)
library(RTMB)
# install RTMBconvenience for logspace_sub:
# remotes::install_github("calbertsen/RTMBconvenience/RTMBconvenience")
library(RTMBconvenience)

#--------------------------------------------------------------------------------
# read in data, manipulate for analysis
#--------------------------------------------------------------------------------

data <- read.csv("data/red_snapper_data.csv")
sf_data <- st_as_sf(data, coords = c("Lon", "Lat"), crs = st_crs("EPSG:4326"))
sf_data$Year <- factor(sf_data$Year)
extent <- st_read("data/red_snapper_extent.shp", crs = st_crs("EPSG:4326"))
coast_sf <- ne_coastline(scale = 10, returnclass = "sf")

bathy <- rast("data/bathy.grd")
sf_bathy <- st_as_sf(as.data.frame(bathy / 1000, xy = TRUE),
    coords = c("x", "y"), crs = st_crs("EPSG:4326")
)

# Add bathymetry at sample points
nearest <- RANN::nn2(st_coordinates(sf_bathy),
    st_coordinates(sf_data),
    k = 1
)
sf_data$bathy <- st_drop_geometry(sf_bathy)[nearest$nn.idx]

# create mesh
mesh <- fm_mesh_2d(st_coordinates(sf_data),
    plot.delay = NULL, cutoff = 0.5
)
# Create matrices in INLA
spde <- fm_fem(mesh, order = 2)
# create projection matrix from vertices to samples
A_is <- fm_evaluator(mesh, loc = st_coordinates(sf_data))$proj$A
# Create extrapolation grid
cellsize <- 0.1
sf_fullgrid <- st_make_grid(extent, cellsize = cellsize)
sf_grid <- st_make_valid(st_intersection(sf_fullgrid, extent))
# Add bathymetry
nearest <- RANN::nn2(st_coordinates(sf_bathy),
    st_coordinates(st_centroid(sf_grid)),
    k = 1
)
sf_grid <- st_sf(sf_grid,
    "bathy" = st_drop_geometry(sf_bathy)[nearest$nn.idx]
)
# create projection matrix from vertices to grid
A_gs <- fm_evaluator(mesh,
    loc = st_coordinates(st_centroid(sf_grid))
)$proj$A
a_g <- st_area(sf_grid)

# plot(st_geometry(sf_grid))
# plot(sf_data[1], add = TRUE)

# Define covariates
Q_formula <- ~ 0 + Data_type
X_formula <- ~ 0 + poly(bathy, 2, raw = TRUE) + Year

# Make Q-matrix
Q_ij <- model.matrix(Q_formula, sf_data)
Q_ij <- Q_ij[, c("Data_typeCount", "Data_typeEncounter")]

# Make X-matrix for data
X_ik <- model.matrix(X_formula, sf_data)

# Make X-matrix for grid
sf_grid$Year <- factor(2014, levels = sort(unique(sf_data$Year)))
frame0 <- model.frame(formula = X_formula, data = sf_data)
terms0 <- terms(frame0)
xlevels <- .getXlevels(terms0, frame0)
terms1 <- delete.response(terms0)
frame1 <- model.frame(terms1, sf_grid, xlev = xlevels)
X_gk <- model.matrix(terms1, frame1)

# scale as per kaskr's fix for non-converging model, see:
# https://github.com/James-Thorson/Spatio-temporal-models-for-ecologists/issues/2

X_ik[, 1:2] <- X_ik[, 1:2] * 100

#--------------------------------------------------------------------------------
# Estimate in RTMB
#--------------------------------------------------------------------------------

data <- list(
    "c_i" = sf_data$Response,
    "e_i" = as.numeric(factor(sf_data$Data_type,
        levels = c("Encounter", "Count", "Biomass_KG")
    )) - 1,
    "X_ik" = X_ik,
    "X_gk" = X_gk,
    "Q_ij" = Q_ij,
    "A_is" = A_is,
    "A_gs" = A_gs,
    "M0" = spde$c0,
    "M1" = spde$g1,
    "M2" = spde$g2
)

par <- list(
    "ln_tau" = 0,
    "ln_kappa" = 0,
    "ln_phi" = 0,
    "finv_power" = 0,
    "gamma_k" = numeric(ncol(X_ik)),
    "eta_j" = numeric(ncol(Q_ij)),
    "omega_s" = numeric(nrow(spde$c0))
)

f <- function(par) {
    getAll(data, par, warn = FALSE)
    phi <- exp(ln_phi)
    power <- 1 + plogis(finv_power)
    logmu_g <- A_gs %*% omega_s + X_gk %*% gamma_k
    logmu_i <- A_is %*% omega_s + Q_ij %*% eta_j + X_ik %*% gamma_k
    Q <- (exp(4 * ln_kappa) *
        M0 + 2 * exp(2 * ln_kappa) * M1 + M2) * exp(2 * ln_tau)
    jnll <- 0
    # Pr(random coefficients)
    jnll <- jnll - dgmrf(omega_s, 0, Q, TRUE)
    # Pr(observations|fixed and random effects)
    for (i in 1:length(c_i)) {
        # Bernoulli
        if (e_i[i] == 0) {
            if (c_i[i] > 0) {
                jnll <- jnll - logspace_sub(log(1.0), -1 * exp(logmu_i[i]))
            } else {
                jnll <- jnll - -1 * exp(logmu_i[i])
            }
        }
        # Poisson
        if (e_i[i] == 1) {
            jnll <- jnll - dpois(
                c_i[i],
                exp(logmu_i[i]), TRUE
            )
        }
        # Tweedie
        if (e_i[i] == 2) {
            jnll <- jnll - dtweedie(
                c_i[i], exp(logmu_i[i]),
                phi, power, TRUE
            )
        }
    }
    REPORT(logmu_g)
    REPORT(jnll)
    jnll
}

obj <- MakeADFun(f, par, random = "omega_s")
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt # --> book solution 25102.1 when scaling
sdr <- sdreport(obj, getJointPrecision = TRUE)
sdr
