library(sf)
library(rnaturalearth)
library(Matrix)
library(RTMB)

#--------------------------------------------------------------------------------
# set up data for analysis
#--------------------------------------------------------------------------------

sf_states <- ne_states(c("United States of America", "Canada"), return = "sf")
sf_states <- sf_states[pmatch(c("Brit", "Alas", "Yukon"), sf_states$name_en), ]
borders <- st_polygon(list(matrix(c(
    -180, 0, -180, 80,
    -50, 80, -50, 0, -180, 0
), byrow = TRUE, ncol = 2)))
sf_states <- st_intersection(st_sfc(borders,
    crs = "+proj=longlat +datum=WGS84"
), sf_states)
sf_states <- st_union(sf_states)

samples <- st_read("samples_3520.csv",
    options = c("X_POSSIBLE_NAMES=X", "Y_POSSIBLE_NAMES=Y")
)
st_crs(samples) <- "+proj=longlat +datum=WGS84"
samples$year <- as.numeric(samples$Year)
samples$species_total <- as.numeric(samples$SpeciesTotal)
samples <- st_intersection(samples, sf_states)
samples <- subset(samples, Year %in% 2019)

# bin into grids
cellsize <- 1
sf_fullgrid <- st_make_grid(sf_states, cellsize = cellsize)
sf_grid <- st_make_valid(st_intersection(sf_fullgrid, sf_states))

grid_i <- as.integer(st_intersects(samples, sf_fullgrid))
plotgrid <- data.frame(st_coordinates(st_centroid(sf_grid)))
n_x <- 1 + round(diff(range(plotgrid$X)) / cellsize)
n_y <- 1 + round(diff(range(plotgrid$Y)) / cellsize)
xy_i <- as.matrix(expand.grid(1:n_x, 1:n_y))[grid_i, ]

#--------------------------------------------------------------------------------
# code up the .hpp files as functions
#--------------------------------------------------------------------------------

conditional_dist_ll <- function(epsilon_xy, rho, sigma2, n_x, n_y) {
    # next line needed for just in time compiler nonsense
    # and needed to be within this function specifically
    "[<-" <- ADoverload("[<-")
    Q_yy <- matrix(0, n_y, n_y)
    diag(Q_yy) <- 1 + rho^2
    for (y in 2:n_y) {
        Q_yy[y - 1, y] <- -rho
        Q_yy[y, y - 1] <- -rho
    }
    # calculate probability
    ans <- 0
    tmp_y <- numeric(n_y)
    Q0_yy <- matrix(Q_yy * (1 - rho^2) / sigma2, n_y, n_y)
    Q1_yy <- Q_yy / sigma2
    for (x in 1:n_x) {
        for (y in 1:n_y) {
            if (x == 1) tmp_y[y] <- epsilon_xy[1, y]
            if (x >= 2) tmp_y[y] <- epsilon_xy[x, y] - rho * epsilon_xy[x - 1, y]
        }
        if (x == 1) ans <- ans + dmvnorm2(tmp_y, Q0_yy)
        if (x >= 2) ans <- ans + dmvnorm2(tmp_y, Q1_yy)
    }
    ans
}

joint_dist_ll <- function(epsilon_xy, rho, sigma2, n_x, n_y) {
    Q_xx <- matrix(0, n_x, n_x)
    diag(Q_xx) <- 1 + rho^2
    Q_xx[cbind(1:(n_x - 1), 2:n_x)] <- -rho # upper diagonal
    Q_xx[cbind(2:n_x, 1:(n_x - 1))] <- -rho # lower diagonal
    Q_yy <- matrix(0, n_y, n_y)
    diag(Q_yy) <- 1 + rho^2
    Q_yy[cbind(1:(n_y - 1), 2:n_y)] <- -rho # upper diagonal
    Q_yy[cbind(2:n_y, 1:(n_y - 1))] <- -rho # lower diagonal
    # calculate probability
    n_z <- n_x * n_y
    Q_zz <- kronecker(Q_yy, Q_xx)
    Q_zz <- Q_zz / sigma2
    epsilon_z <- as.vector(epsilon_xy)
    dmvnorm2(epsilon_z, Q_zz)
}

dmvnorm2 <- function(x, Q, give_log = TRUE) {
    n_x <- length(x)
    logres <- 0
    logdet_Q <- determinant(Q)$modulus[1]
    Tmp_x <- Q %*% x
    logres <- logres + 0.5 * logdet_Q
    logres <- logres - 0.5 * sum(x %*% Tmp_x)
    if (give_log) {
        return(logres)
    } else {
        return(exp(logres))
    }
}

#--------------------------------------------------------------------------------
# fit via RTMB
#--------------------------------------------------------------------------------

par <- list(
    "beta0" = 0,
    "log_sig2" = 0,
    "logit_rho" = 0,
    "epsilon_xy" = matrix(0, n_x, n_y)
)
data <- list(
    "c_i" = as.numeric(samples$SpeciesTotal),
    "xy_i" = xy_i
)

f <- function(par) {
    getAll(data, par, warn = FALSE)
    jnll <- 0
    sig2 <- exp(log_sig2)
    rho <- 1 - plogis(logit_rho)
    n_x <- nrow(epsilon_xy)
    n_y <- ncol(epsilon_xy)
    # calculate probability of random coefficients:
    # conditional in dimension y, joint in dimension x
    if (ctl == 1) {
        jnll <- jnll - conditional_dist_ll(epsilon_xy, rho, sig2, n_x, n_y)
    }
    # using kronecker product of precision in x and y
    if (ctl == 2) {
        jnll <- jnll - joint_dist_ll(epsilon_xy, rho, sig2, n_x, n_y)
    }
    # using built-in functions
    if (ctl == 3) {
        f1 <- function(x) {
            dautoreg(x,
                mu = 0, phi = rho,
                log = TRUE
            )
        }
        f2 <- function(x) {
            dautoreg(x,
                mu = 0, phi = rho,
                log = TRUE
            )
        }
        jnll <- jnll - dseparable(f1, f2)(epsilon_xy,
            scale = sqrt(sig2) / sqrt(1 - rho^2) / sqrt(1 - rho^2))
    }
    # predict densities
    D_xy <- matrix(0, nrow = nrow(epsilon_xy), ncol = ncol(epsilon_xy))
    for (x in 1:nrow(epsilon_xy)) {
        for (y in 1:ncol(epsilon_xy)) {
            D_xy[x, y] <- exp(beta0 + epsilon_xy[x, y])
        }
    }
    # Pr(data|fixed and random coefficients):
    for (i in 1:length(c_i)) {
        jnll <- jnll - dpois(c_i[i], D_xy[xy_i[i, 1], xy_i[i, 2]], TRUE)
    }
    REPORT(D_xy)
    ADREPORT(D_xy)
    jnll
}

# fit model conditional in y, joint in x:
data$ctl <- 1
obj <- MakeADFun(f, par, random = "epsilon_xy")
opt1 <- nlminb(obj$par, obj$fn, obj$gr)
opt1
sdr <- sdreport(obj)
sdr

# fit model using Kronecker product of precision in both dimensions
TapeConfig(atomic = "disable") # trick for use if RTMB hangs doing matrix algebra
data$ctl <- 2
obj <- MakeADFun(f, par, random = "epsilon_xy")
opt2 <- nlminb(obj$par, obj$fn, obj$gr)
opt2
sdr <- sdreport(obj)
sdr
TapeConfig(atomic = "enable")

# fit model using built-in functions
data$ctl <- 3
obj <- MakeADFun(f, par, random = "epsilon_xy")
opt3 <- nlminb(obj$par, obj$fn, obj$gr)
opt3
sdr <- sdreport(obj)
sdr
