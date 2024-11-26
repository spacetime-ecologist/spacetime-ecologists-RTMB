library(mvtnorm)
library(stars)
library(RTMB)

#--------------------------------------------------------------------------------
# simulate data
#--------------------------------------------------------------------------------

set.seed(101)
SD_km <- 100
peak_km <- 2 # km
best_elevation <- 1
SD_logelev <- 0.5
n_indiv <- 500

get_elevation <- function(loc) { # simulation elevation
    e0 <- mvtnorm::dmvnorm(c(0, 0), mean = c(0, 0), sigma = SD_km^2 * diag(2))
    peak_km * mvtnorm::dmvnorm(loc, mean = c(0, 0), sigma = SD_km^2 * diag(2)) / e0
}

get_density <- function(loc) { # simulate density
    elev <- get_elevation(loc)
    exp(-1 * (log(elev / best_elevation))^2 / SD_logelev^2)
}

# function for rejection sampling for locations of individuals
max_density <- optim(
    fn = get_density, par = c(0, 100),
    control = list(fnscale = -1)
)$value

simulate_location <- function(...) {
    loc <- NULL
    while (is.null(loc)) {
        samp <- runif(n = 2, min = -200, max = 200)
        D <- get_density(samp)
        rand <- runif(n = 1, min = 0, max = max_density)
        if (rand < D) loc <- samp
    }
    loc
}
loc_i <- t(sapply(1:n_indiv, FUN = simulate_location))

#--------------------------------------------------------------------------------
# format data for glm
#--------------------------------------------------------------------------------

samples <- data.frame("x" = loc_i[, 1], "y" = loc_i[, 2])
samples <- st_as_sf(samples, coords = c("x", "y"))

# Get count in each grid cell
grid_size <- 10
grid <- st_make_grid(st_bbox(c(xmin = -200, xmax = 200, ymin = -200, ymax = 200)),
    cellsize = grid_size
)
grid_i <- st_intersects(samples, grid)
n_i <- tapply(rep(1, nrow(samples)),
    INDEX = factor(unlist(grid_i), levels = 1:length(grid)),
    FUN = sum
)

# Convert to a data frame
data <- data.frame(st_coordinates(st_centroid(grid)),
    "N" = ifelse(is.na(n_i), 0, n_i)
)

data$elev <- get_elevation(data[, c("X", "Y")])

#--------------------------------------------------------------------------------
# fit with RTMB
#--------------------------------------------------------------------------------

formula <- ~ log(elev) + I(log(elev)^2)
X_ij <- model.matrix(formula, data = data)
data <- list("y_i" = na.omit(data)$N, "X_ij" = X_ij)
par <- list("b_j" = rep(0, ncol(X_ij)))

f <- function(par) {
    getAll(data, par, warn = FALSE)
    y_i <- OBS(y_i)
    n_i <- nrow(X_ij)
    log_mu <- X_ij %*% b_j # note: could be a loop
    jnll <- -sum(dpois(y_i, exp(log_mu), TRUE))
    REPORT(log_mu)
    jnll
}

obj <- MakeADFun(f, par)
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt # book solution 867.6721
sdr <- sdreport(obj)
sdr

obj$simulate() # --> can be used to create simulation residuals etc.
