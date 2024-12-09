# dgmrf() scale bug fixed by fishfollower and kaskr
# see https://github.com/kaskr/RTMB/issues/38

library(fmesher)
library(sf)
library(rasterVis)
library(RTMB)

#--------------------------------------------------------------------------------
# read in data and format for analysis
#--------------------------------------------------------------------------------

pollock <- readRDS("pollock.rds")
pollock <- st_as_sf(pollock,
    coords = c("Long", "Lat"),
    crs = "+proj=longlat +datum=WGS84"
)
pollock <- st_transform(pollock,
    crs = st_crs("+proj=utm +zone=2 +datum=WGS84 +units=km")
)
coldpool <- readRDS("coldpool.rds")
survey_domain <- readRDS("survey_domain.rds")
survey_domain <- st_sfc(survey_domain, crs = "+proj=longlat +datum=WGS84")
survey_domain <- st_transform(survey_domain, crs = st_crs(pollock))

# Country shapefiles for plotting
sf_maps <- rnaturalearth::ne_countries(
    return = "sf",
    scale = 10, country = "united states of america"
)
sf_maps <- st_transform(sf_maps, crs = st_crs(pollock))

# Make triangulated mesh
mesh <- fm_mesh_2d(st_coordinates(pollock), cutoff = 100, refine = TRUE)
# Create matrices in INLA
spde <- fm_fem(mesh, order = 2)
# create projection matrix from vertices to samples
A_is <- fm_evaluator(mesh, loc = st_coordinates(pollock))$proj$A
# Create extrapolation grid
cellsize <- 25
grid <- st_make_grid(survey_domain, cellsize = cellsize)
grid <- st_intersection(grid, survey_domain)
# create projection matrix from vertices to grid
A_gs <- fm_evaluator(mesh, loc = st_coordinates(st_centroid(grid)))$proj$A
year_set <- min(pollock$Year):max(pollock$Year)

#--------------------------------------------------------------------------------
# fit with RTMB
#--------------------------------------------------------------------------------

data <- list(
    "n_t" = length(year_set),
    "a_g" = as.numeric(st_area(grid)),
    "z_g" = st_coordinates(st_centroid(grid))[, 2],
    "b_i" = pollock$Wt,
    "a_i" = pollock$AreaSwept_ha / 100, # convert hectares to km^2
    "t_i" = pollock$Year - min(pollock$Year) + 1,
    "A_is" = A_is,
    "A_gs" = A_gs,
    "M0" = spde$c0,
    "M1" = spde$g1,
    "M2" = spde$g2
)

par <- list(
    "beta_t" = rep(0, data$n_t),
    "ln_tauO" = log(1),
    "ln_tauE" = log(1),
    "ln_kappa" = 1,
    "ln_phi" = log(1),
    "logit_rhoE" = 0,
    "finv_power" = 0,
    "omega_s" = rep(0, mesh$n),
    "epsilon_st" = matrix(0.1, nrow = mesh$n, ncol = data$n_t)
)

f <- function(par) {
    getAll(data, par, warn = FALSE)
    jnll <- 0
    rhoE <- 1 - plogis(logit_rhoE)
    n_i <- nrow(A_is)
    n_g <- nrow(A_gs)
    Q <- exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2
    # Pr(random effects):
    # spatial effect
    jnll <- jnll - dgmrf(omega_s,
        mu = 0, Q = Q,
        log = TRUE, scale = 1 / exp(ln_tauO)
    )
    # temporally evolving AR-1 st field two ways (ctl = 0 or 1)
    # first, using built in separable functions:
    if (ctl == 0) {
        f1 <- function(x) {
            dgmrf(x,
                mu = 0, Q = Q,
                log = TRUE, scale = 1 / exp(ln_tauE)
            )
        }
        f2 <- function(x) {
            dautoreg(x,
                mu = 0, phi = rhoE,
                log = TRUE, scale = 1 / exp(ln_tauE) / sqrt(1 - rhoE^2)
            )
        }
        jnll <- jnll - dseparable(f1, f2)(epsilon_st)
    }
    # second, space-time "by-hand:"
    if (ctl == 1) {
        jnll <- jnll - dgmrf(epsilon_st[, 1],
            mu = 0, Q = Q,
            log = TRUE, scale = 1 / exp(ln_tauE) / sqrt(1 - rhoE^2)
        )
        for (t in 2:n_t) {
            jnll <- jnll - dgmrf(epsilon_st[, t],
                mu = rhoE * epsilon_st[, t - 1],
                Q = Q, log = TRUE, scale = 1 / exp(ln_tauE)
            )
        }
    }
    # latent density field
    omega_g <- A_gs %*% omega_s
    epsilon_gt <- A_gs %*% epsilon_st
    ln_d_gt <- matrix(0, nrow = n_g, ncol = n_t)
    for (t in 1:n_t) {
        for (g in 1:n_g) {
            ln_d_gt[g, t] <- beta_t[t] + omega_g[g] + epsilon_gt[g, t]
        }
    }
    # Pr(data | fixed and random effects):
    omega_i <- A_is %*% omega_s
    epsilon_it <- A_is %*% epsilon_st
    bhat_i <- numeric(n_i)
    for (i in 1:n_i) {
        bhat_i[i] <- a_i[i] *
            exp(beta_t[t_i[i]] + omega_i[i] + epsilon_it[i, t_i[i]])
        jnll <- jnll - dtweedie(
            b_i[i],
            bhat_i[i],
            exp(ln_phi),
            1 + plogis(finv_power),
            TRUE
        )
    }
    # derived quantities
    b_t <- zmean_t <- numeric(n_t)
    for (t in 1:n_t) {
        for (g in 1:n_g) {
            b_t[t] <- b_t[t] + a_g[g] * exp(ln_d_gt[g, t])
        }
        for (g in 1:n_g) {
            zmean_t[t] <- zmean_t[t] + z_g[g] * a_g[g] * exp(ln_d_gt[g, t]) / b_t[t]
        }
    }
    REPORT(ln_d_gt)
    REPORT(bhat_i)
    REPORT(b_t)
    ADREPORT(b_t)
    REPORT(zmean_t)
    ADREPORT(zmean_t)
    jnll
}

data$ctl <- 0 # 0 -> using st random effects using dseparable(), 1 -> 'by-hand'
obj <- MakeADFun(f, par, random = c("omega_s", "epsilon_st"))

opt <- nlminb(obj$par, obj$fn, obj$gr,
    control = list(trace = 1, eval.max = 1e4, iter.max = 1e4)
)
opt

sdr <- sdreport(obj)
sdr
