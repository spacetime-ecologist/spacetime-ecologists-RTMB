library(sf)
library(fmesher)
library(splines2)
library(rnaturalearth)
library(RTMB)

#--------------------------------------------------------------------------------
# read in data and format for analysis
#--------------------------------------------------------------------------------

ozone <- st_read("2019_ozone.csv",
    options = c("X_POSSIBLE_NAMES=X", "Y_POSSIBLE_NAMES=Y"),
    crs = st_crs("+proj=longlat +datum=WGS84")
)
ozone$Daily.Max.8.hour.Ozone.Concentration <-
    as.numeric(ozone$Daily.Max.8.hour.Ozone.Concentration)
ozone$Date <- as.POSIXlt(ozone$Date, tryFormats = "%m/%d/%Y")
ozone$Julian <- as.numeric(julian(ozone$Date,
    origin = as.POSIXct("01-01-2019", tryFormats = "%m-%d-%Y")
)) / 365
ozone <- ozone[which(ozone$Daily.Max.8.hour.Ozone.Concentration > 0), ]

# Define states
state_set <- c(
    "Virginia", "North Carolina", "South Carolina",
    "Georgia", "Florida", "Maryland", "Delaware", "District of Columbia",
    "New Jersey", "Pennsylvania", "New York"
)
states_sf <- ne_states(c("United States of America", "Canada"),
    return = "sf"
)
states_sf <- states_sf[pmatch(state_set, states_sf$name), ]
domain_sf <- st_union(states_sf)

# Create extrapolation grid
cellsize <- 0.25
grid <- st_make_grid(states_sf, cellsize = cellsize)
grid <- st_intersection(grid, domain_sf)

# create mesh
mesh <- fm_mesh_2d(st_coordinates(st_centroid(grid)),
    refine = TRUE, cutoff = 0.2
)
# Create matrices in INLA
spde <- fm_fem(mesh, order = 2)

season_df <- 6
season_plots <- 12

A_is <- fm_evaluator(mesh, loc = st_coordinates(ozone))$proj$A
A_gs <- fm_evaluator(mesh, loc = st_coordinates(st_centroid(grid)))$proj$A
Afull_gs <- A_gs
for (r in 2:season_plots) Afull_gs <- rbind(Afull_gs, A_gs)

# formula
formula <- ~ 0 + mSpline(Julian,
    knots =
        seq(0, 1, length = season_df + 1)[2:season_df],
    degree = 3, intercept = TRUE, periodic = TRUE,
    Boundary.knots = c(0, 1)
)
options(na.action = "na.pass")
X_ij <- model.matrix(formula, data = data.frame("Julian" = c(ozone$Julian)))
X_gj <- model.matrix(formula,
    data =
        data.frame("Julian" = rep(seq(0, 1, length = season_plots + 1)[1:season_plots],
            each = length(grid)
        ))
)

#--------------------------------------------------------------------------------
# fit with RTMB
#--------------------------------------------------------------------------------

data <- list(
    "y_i" = ozone$Daily.Max.8.hour.Ozone.Concentration,
    "A_is" = A_is,
    "A_gs" = Afull_gs,
    "M0" = spde$c0,
    "M1" = spde$g1,
    "M2" = spde$g2,
    "X_ij" = X_ij,
    "X_gj" = X_gj
)
par <- list(
    "ln_tau" = 0,
    "ln_kappa" = 0,
    "log_sig" = 0,
    "beta_j" = rep(0, ncol(X_ij)),
    "xi_sj" = matrix(0, ncol = ncol(X_ij), nrow = mesh$n)
)

f <- function(par) {
    getAll(data, par, warn = FALSE)
    jnll <- 0
    # derived variables
    range <- sqrt(8) / exp(ln_kappa)
    sigE <- 1 / sqrt(4 * pi * exp(2 * ln_tau) * exp(2 * ln_kappa))
    # set up Q
    Q <- exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2
    # Pr(random coefficients), independent spatial field for each time slice j:
    for (j in 1:ncol(xi_sj)) {
        jnll <- jnll - dgmrf(xi_sj[, j], 0, Q, TRUE, scale = 1 / exp(ln_tau))
    }
    xi_ij <- A_is %*% xi_sj # project to datum locations
    phi_i <- X_ij %*% beta_j
    yhat_i <- numeric(length(y_i))
    # Pr(data | fixed + random effects):
    for (i in 1:length(y_i)) {
        yhat_i[i] <- exp(sum(xi_ij[i, ] * X_ij[i, ]) + phi_i[i])
        # shape = 1/CV^2,  scale = mean*CV^2
        jnll <- jnll - dgamma(
            x = y_i[i],
            shape = 1 / exp(2 * log_sig),
            scale = yhat_i[i] * exp(2 * log_sig),
            log = TRUE
        )
    }
    # extrapolation
    xi_gj <- A_gs %*% xi_sj
    phi_g <- X_gj %*% beta_j
    yhat_g <- rep(0, nrow(A_gs))
    for (g in 1:nrow(A_gs)) {
        yhat_g[g] <- exp(sum(xi_gj[g, ] * X_gj[g, ]) + phi_g[g])
    }
    REPORT(yhat_g)
    REPORT(yhat_i)
    REPORT(range)
    REPORT(sigE)
    jnll
}

obj <- MakeADFun(f, par, random = "xi_sj")
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt
sdr <- sdreport(obj)
sdr
