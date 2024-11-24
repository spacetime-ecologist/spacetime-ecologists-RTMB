library(sf)
library(rnaturalearth)
library(fmesher)
library(RTMB)

#--------------------------------------------------------------------------------
# read in data and manipulate for analysis
#--------------------------------------------------------------------------------

out <- st_read("data/samples_3520.csv",
    options = c("X_POSSIBLE_NAMES=X", "Y_POSSIBLE_NAMES=Y")
)
st_crs(out) <- "+proj=longlat +datum=WGS84"
out$year <- as.numeric(out$Year)
out$species_total <- as.numeric(out$SpeciesTotal)
out <- subset(out, State %in% c("BritCol", "Alaska", "Yukon"))
out <- subset(out, Year %in% 2019)

sf_states <- ne_states(c("United States of America", "Canada"), return = "sf")
sf_states <- sf_states[pmatch(c("Brit", "Alas", "Yukon"), sf_states$name_en), ]
borders <- st_polygon(list(matrix(c(-180, 0, -180, 80, -50, 80, -50, 0, -180, 0),
    byrow = TRUE, ncol = 2
)))
sf_states <- st_intersection(
    st_sfc(borders, crs = "+proj=longlat +datum=WGS84"),
    sf_states
)
sf_states <- st_union(sf_states)

#--------------------------------------------------------------------------------
# set up spde approximation
#--------------------------------------------------------------------------------

xy_i <- st_coordinates(out)
mesh <- fm_mesh_2d(xy_i, refine = TRUE, cutoff = 0.5)

# create matrices in fmesher / INLA
spde <- fm_fem(mesh, order = 2)

# create projection matrix from vertices to sample locations
A_is <- fm_evaluator(mesh, loc = xy_i)$proj$A

# create extrapolation grid
cellsize <- 1
grid <- st_make_grid(sf_states, cellsize = cellsize)

# create projection matrix from vertices to grid
A_gs <- fm_evaluator(mesh, loc = st_coordinates(st_centroid(grid)))$proj$A

#--------------------------------------------------------------------------------
# fit spatial GLMM using SPDE approximation via RTMB
#--------------------------------------------------------------------------------

data <- list(
    "c_i" = out$species_total, "A_is" = A_is, "A_gs" = A_gs,
    "M0" = spde$c0, "M1" = spde$g1, "M2" = spde$g2
)
par <- list(
    "beta0" = 0, "ln_tau" = 0, "ln_kappa" = 0,
    "omega_s" = numeric(nrow(spde$c0))
)

f <- function(par) {
    getAll(data, par, warn = FALSE)
    range <- sqrt(8) / exp(ln_kappa)
    sigE <- 1 / sqrt(4 * pi * exp(2 * ln_tau) * exp(2 * ln_kappa))
    Q <- exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2
    jnll <- 0
    jnll <- jnll - dgmrf(omega_s, 0.0, Q, TRUE, scale = 1 / exp(ln_tau))
    omega_i <- A_is %*% omega_s
    jnll <- jnll - sum(dpois(c_i, exp(beta0 + omega_i), TRUE))
    logN_g <- beta0 + A_gs %*% omega_s
    REPORT(logN_g)
    REPORT(Q)
    REPORT(range)
    REPORT(sigE)
    jnll
}

obj <- MakeADFun(f, par, random = "omega_s")
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt # book solution 196.619
sdr <- sdreport(obj, bias.correct = TRUE)
sdr
