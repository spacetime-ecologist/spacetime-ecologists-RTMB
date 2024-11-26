library(sf)
library(fmesher)
library(rnaturalearth)
library(RTMB)

#--------------------------------------------------------------------------------
# read in data and format for analysis
#--------------------------------------------------------------------------------

ice <- read.csv("data/Ice.csv")

# project data
sf_ice <- st_as_sf(ice, coords = c("Longitude", "Latitude"))
st_crs(sf_ice) <- "+proj=longlat +datum=WGS84"
sf_ice <- st_transform(sf_ice,
    crs = st_crs("+proj=laea +lat_0=90 +lon_0=-30 +units=km")
)
sf_pole <- st_point(c(0, 90))
sf_pole <- st_sfc(sf_pole,
    crs = "+proj=longlat +datum=WGS84"
)
sf_pole <- st_transform(sf_pole,
    crs = st_crs("+proj=laea +lat_0=90 +lon_0=-30 +units=km")
)
sf_pole <- st_buffer(sf_pole, dist = 3000)
sf_ice <- st_intersection(sf_ice, sf_pole)

# country shapefiles for plotting
sf_maps <- ne_countries(
    return = "sf", scale = 10,
    continent = c("north america", "europe", "asia")
)
sf_maps <- st_transform(sf_maps,
    crs = st_crs("+proj=laea +lat_0=90 +lon_0=-30 +units=km")
)
sf_maps <- st_union(sf_maps)

# shapefile for water
sf_water <- st_difference(st_as_sfc(st_bbox(sf_maps)), sf_maps)

# create mesh
xy_i <- st_coordinates(sf_ice)
mesh <- fm_mesh_2d(xy_i, plot.delay = NULL, cutoff = 200)
# Create matrices in INLA
spde <- fm_fem(mesh, order = 2)
# create projection matrix from vertices to samples
A_is <- fm_evaluator(mesh, loc = xy_i)$proj$A
# Create extrapolation grid
cellsize <- 50
sf_grid <- st_make_grid(sf_pole, cellsize = cellsize)
# Restrict to water
grid_i <- st_intersects(sf_water, sf_grid)
sf_grid <- sf_grid[unique(unlist(grid_i))]
# Restrict to 3000 km from North Pole
grid_i <- st_intersects(sf_pole, sf_grid)
sf_grid <- sf_grid[unique(unlist(grid_i))]
# create projection matrix from vertices to grid
A_gs <- fm_evaluator(mesh, loc = st_coordinates(st_centroid(sf_grid)))$proj$A
# Area for each grid cell
a_g <- rep(cellsize^2, length(sf_grid))
# sum(a_g) < pi*3000^2 due to missing land

# set up data and parameters
rmatrix <- function(nrow, ncol, sd = 1) matrix(sd * rnorm(nrow * ncol), nrow = nrow)
n_factors <- 2
n_years <- max(sf_ice$Year) - min(sf_ice$Year) + 1

#--------------------------------------------------------------------------------
# fit with RTMB
#--------------------------------------------------------------------------------

data <- list(
    "y_i" = sf_ice$Ice,
    "t_i" = sf_ice$Year - min(ice$Year) + 1,
    "A_is" = as.matrix(A_is),
    "A_gs" = as.matrix(A_gs),
    "a_g" = a_g,
    "M0" = spde$c0,
    "M1" = spde$g1,
    "M2" = spde$g2
)

par <- list(
    "beta0" = 0,
    "ln_tau" = 0,
    "ln_kappa" = log(1),
    "ln_sigma" = 0,
    "L_ft" = rmatrix(nrow = n_factors, ncol = n_years),
    "omega_s" = rnorm(nrow(spde$c0)),
    "epsilon_sf" = rmatrix(nrow = nrow(spde$c0), ncol = n_factors)
)

f <- function(par) {
    getAll(data, par, warn = FALSE)
    jnll <- 0
    Q <- exp(2 * ln_tau) * (exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2)
    jnll <- jnll - dgmrf(omega_s, 0, Q, TRUE)
    for (t in 1:ncol(epsilon_sf)) {
        jnll <- jnll - dgmrf(epsilon_sf[, t], 0, Q, log = TRUE)
    }
    epsilon_it <- A_is %*% epsilon_sf %*% L_ft
    omega_i <- A_is %*% omega_s
    ypred <- numeric(length(y_i))
    for (i in 1:length(y_i)) {
        ypred[i] <- beta0 + omega_i[i] + epsilon_it[i, t_i[i]]
    }
    jnll <- jnll - sum(dnorm(y_i, ypred, exp(ln_sigma), TRUE))
    jnll
}

map <- list("L_ft" = matrix(1:prod(dim(par$L_ft)), nrow = n_factors))
map$L_ft[lower.tri(map$L_ft)] <- NA
map$L_ft <- factor(map$L_ft)
par$L_ft[lower.tri(par$L_ft)] <- 0

TapeConfig(atomic = "disable") # nifty trick for use if RTMB hangs

obj <- MakeADFun(f, par,
    map = map,
    random = c("omega_s", "epsilon_sf")
)

opt <- nlminb(obj$par, obj$fn, obj$gr,
    control = list(trace = 1, eval.max = 1e4, iter.max = 1e4)
)
opt # -9766.518 is book solution

TapeConfig(atomic = "enable") # set it back
