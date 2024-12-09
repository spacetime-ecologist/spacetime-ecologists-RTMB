library(sf)
library(fmesher)
library(rnaturalearth)
library(RTMB)
# note, also requires stars package

#--------------------------------------------------------------------------------
# read in data and manipulate for analysis
#--------------------------------------------------------------------------------

# ozone data
ozone <- st_read("2019_ozone.csv",
    options = c("X_POSSIBLE_NAMES=X", "Y_POSSIBLE_NAMES=Y"),
    crs = st_crs("+proj=longlat +datum=WGS84")
)
ozone$Daily.Max.8.hour.Ozone.Concentration <-
    as.numeric(ozone$Daily.Max.8.hour.Ozone.Concentration)
ozone <- subset(ozone, Date == "07/01/2019")

# density data
pop_dens <- st_read("population_density.csv",
    options = c("X_POSSIBLE_NAMES=X", "Y_POSSIBLE_NAMES=Y"),
    crs = st_crs("+proj=longlat +datum=WGS84")
)
pop_dens$Dens2000 <- as.numeric(pop_dens$Dens2000)
pop_dens$Dens2020 <- as.numeric(pop_dens$Dens2020)
pop_dens_stars <- stars::st_rasterize(pop_dens)
# plot(pop_dens_stars["Dens2000"])

# Define states
state_set <- c(
    "Virginia", "North Carolina",
    "South Carolina", "Georgia", "Florida", "Maryland",
    "Delaware", "District of Columbia", "New Jersey",
    "Pennsylvania", "New York"
)
states_sf <- ne_states(c("United States of America", "Canada"), return = "sf")
states_sf <- states_sf[pmatch(state_set, states_sf$name), ]
domain_sf <- st_union(states_sf)

# Get pop-dens for nearest sample
index <- st_nearest_feature(ozone, pop_dens)
ozone$pop_dens <- pop_dens[index, "Dens2020"]
# plot(ozone$pop_dens)

#--------------------------------------------------------------------------------
# unequal distance 2d autoregressive
#--------------------------------------------------------------------------------

# reate extrapolation grid
cellsize <- 0.25
grid <- st_make_grid(states_sf, cellsize = cellsize)
grid <- st_intersection(grid, domain_sf)

# create mesh
mesh <- fm_mesh_2d(st_coordinates(st_centroid(grid)), refine = TRUE, cutoff = 0.2)

# create matrices in INLA
spde <- fm_fem(mesh, order = 2)

# create projection matrix from vertices to samples
A_is <- fm_evaluator(mesh, loc = st_coordinates(ozone))$proj$A

# create projection matrix from vertices to grid
A_gs <- fm_evaluator(mesh, loc = st_coordinates(st_centroid(grid)))$proj$A

# formula
formula <- ~ 0 + log(Dens2000)
options(na.action = "na.pass")
DF_i <- data.frame(pop_dens[st_nearest_feature(ozone, pop_dens), ])
X_ij <- model.matrix(formula, data = DF_i)
DF_g <- data.frame(pop_dens[st_nearest_feature(grid, pop_dens), ])
X_gj <- model.matrix(formula, data = DF_g)

# Expansion rates
a_g <- as.numeric(st_area(grid))
d_g <- DF_g$Dens2020

#--------------------------------------------------------------------------------
# fit in RTMB
#--------------------------------------------------------------------------------

data <- list(
    "y_i" = ozone$Daily.Max.8.hour.Ozone.Concentration, "A_is" = A_is,
    "A_gs" = A_gs, "M0" = spde$c0, "M1" = spde$g1,
    "M2" = spde$g2, "X_ij" = X_ij, "X_gj" = X_gj, "a_g" = a_g, "d_g" = d_g
)
par <- list(
    "beta0" = 0, "ln_tau" = 0, "ln_kappa" = 0, "log_sig" = 0,
    "b_j" = rep(0, ncol(X_ij)), "omega_s" = rnorm(mesh$n)
)

f <- function(par) {
    getAll(data, par, warn = FALSE)
    y_i <- OBS(y_i)
    range <- sqrt(8) / exp(ln_kappa)
    sigE <- 1 / sqrt(4 * pi * exp(2 * ln_tau) * exp(2 * ln_kappa))
    Q <- exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2
    jnll <- 0
    # Pr(random coefficients):
    jnll <- jnll - dgmrf(omega_s, 0.0, Q, TRUE, scale = 1 / exp(ln_tau))
    omega_i <- A_is %*% omega_s
    phi_i <- numeric(nrow(X_ij))
    yhat_i <- numeric(length(y_i))
    phi_i <- X_ij %*% b_j
    # Pr(data|fixed + random effects):
    for (i in 1:length(y_i)) {
        yhat_i[i] <- exp(beta0 + omega_i[i] + phi_i[i])
        # shape = 1/cv^2, scale = mean*cv^2
        jnll <- jnll - dgamma(
            y_i[i],
            shape = 1 / exp(2 * log_sig),
            scale = yhat_i[i] * exp(2 * log_sig),
            log = TRUE
        )
    }
    jnll
}

obj <- MakeADFun(f, par, random = "omega_s")

opt <- nlminb(obj$par, obj$fn, obj$gr)
opt

sdr <- sdreport(obj,
    bias.correct = TRUE,
    getJointPrecision = TRUE,
    bias.correct.control = list(sd = TRUE)
)
sdr
