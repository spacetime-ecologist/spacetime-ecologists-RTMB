library(RTMB)
library(sf)
library(numDeriv)
library(mvtnorm)
library(pracma)
library(plotrix)
library(stars)
library(sem)
source("shared_functions/build_ram.R") # --> needed for build_ram()

#--------------------------------------------------------------------------------
# create preference function, simulate the tracks w/taxis,
# plot preference gradient that changes through time
#--------------------------------------------------------------------------------

n_t <- 1000
beta_t <- seq(0.5, -0.5, length = n_t)
gamma_t <- rep(0, n_t)

get_preference <- function(beta, gamma, s) {
    beta * sqrt(s[1]^2 + s[2]^2) + gamma * (s[2])
}

# Simulation function
simulate_track <- function(n_t, beta_t, gamma_t, get_preference) {
    s_t <- array(NA, dim = c(n_t, 2), dimnames = list(NULL, c("x", "y")))
    # Euler approximation
    for (t in seq_len(n_t)) {
        if (t == 1) s_t[t, ] <- c(0, 0)
        if (t >= 2) {
            gradient <- grad(
                f = get_preference, x0 = s_t[t - 1, , drop = FALSE],
                beta = beta_t[t - 1], gamma = gamma_t[t - 1]
            )
            s_t[t, ] <- s_t[t - 1, ] + gradient + rmvnorm(n = 1, sigma = diag(2))
        }
    }
    s_t
}

#--------------------------------------------------------------------------------
# simulate a path

set.seed(4)
x_iz <- simulate_track(
    n_t = n_t, beta_t = beta_t,
    gamma_t = gamma_t, get_preference = get_preference
)
y_iz <- x_iz + mvtnorm::rmvnorm(n = nrow(x_iz), sigma = diag(2))
deltaT_i <- c(NA, rep(1, nrow(y_iz) - 1))
t_i <- c(0, cumsum(deltaT_i[-1]))

#--------------------------------------------------------------------------------
# set up covariance function

make_covar <- function(s2_z, RAM, RAMstart_z, nrows, ncols) {
    I_rr <- diag(nrows)
    if (nrow(RAM) > 0) { # structural equation covariance
        Rho <- Gamma <- matrix(0, nrows, nrows)
        for (zI in 1:nrow(RAM)) {
            if (RAM[zI, 4] > 0) { # RAM(zI,3) in C++ becomes RAM[zI,4] in R (1-indexed)
                if (RAM[zI, 1] == 1) {
                    Rho[RAM[zI, 2], RAM[zI, 3]] <- s2_z[RAM[zI, 4]]
                }
                if (RAM[zI, 1] == 2) {
                    Gamma[RAM[zI, 2], RAM[zI, 3]] <- s2_z[RAM[zI, 4]]
                }
            } else {
                if (RAM[zI, 1] == 1) {
                    Rho[RAM[zI, 2], RAM[zI, 3]] <- RAMstart_z[zI]
                }
                if (RAM[zI, 1] == 2) {
                    Gamma[RAM[zI, 2], RAM[zI, 3]] <- RAMstart_z[zI]
                }
            }
        }
        L_rc <- I_rr - Rho
        L_rc <- solve(L_rc)
        L_rc <- L_rc %*% Gamma
        Cov_rr <- L_rc %*% t(L_rc)
    } else { # assemble loadings matrix covariance
        L_rc <- matrix(0, nrows, ncols)
        ctr <- 1
        if (ncols > 0) {
            for (i in 1:nrows) {
                for (j in 1:ncols) {
                    x_iz[i, ]
                    if (i >= j) {
                        L_rc[i, j] <- s2_z[ctr]
                        ctr <- ctr + 1
                    } else {
                        L_rc[i, j] <- 0.0
                    }
                }
            }
        }
        # diagonal and equal covariance matrix
        Cov_rr <- I_rr * exp(s2_z[ctr])
        # add factor-model covariance matrix
        Cov_rr <- Cov_rr + L_rc %*% t(L_rc)
    }
    Cov_rr
}

#--------------------------------------------------------------------------------
# set up simplest model

set.seed(4)
x_iz <- simulate_track(
    n_t = n_t, beta_t = beta_t, gamma_t = gamma_t,
    get_preference = get_preference
)
y_iz <- x_iz + mvtnorm::rmvnorm(n = nrow(x_iz), sigma = diag(2))
deltaT_i <- c(NA, rep(1, nrow(y_iz) - 1))
t_i <- c(0, cumsum(deltaT_i[-1]))

# Define drift as a function of time
formula <- ~0
X_ij <- as.matrix(model.matrix(formula, data = data.frame("t_i" = t_i)))

# Build inputs
data <- list(
    "y_iz" = y_iz,
    "deltaT_i" = deltaT_i,
    "error2_i" = rep(1, nrow(y_iz)),
    "X_ij" = X_ij,
    "n_factors" = 0,
    "RAM" = matrix(nrow = 0, ncol = 4)
)

par <- list(
    "sigma2_z" = log(1),
    "x_iz" = data$y_iz,
    "beta_jz" = matrix(0, nrow = ncol(data$X_ij), ncol = 2)
)

# drop some data
which_include <- seq(1, nrow(data$y_iz), length = 50)
data$y_iz[-which_include, ] <- NA

f <- function(par) {
    getAll(data, par, warn = FALSE)
    n_i <- nrow(y_iz)
    n_z <- ncol(y_iz)
    I_zz <- diag(n_z)
    gamma_iz <- X_ij %*% beta_jz
    V_zz <- make_covar(sigma2_z, RAM, sigma2_z, n_z, n_factors)
    Gsum_iz <- matrix(0, n_i, n_z)
    Gsum_iz[1, ] <- gamma_iz[1, ] + x_iz[1, ]
    jnll <- 0
    # Pr(random coefficients)
    for (i in 2:n_i) {
        if (!is.na(deltaT_i[i])) {
            Vi_zz <- deltaT_i[i] * V_zz
            jnll <- jnll - RTMB::dmvnorm(
                x_iz[i, ], x_iz[i - 1, ] + gamma_iz[i - 1, ], Vi_zz, TRUE
            )
        }
        Gsum_iz[i, ] <- Gsum_iz[i - 1, ] + gamma_iz[i, ] # accumulate covar sums
    }
    # Pr(data|fixed,random par)
    for (i in 1:n_i) {
        if (!any(is.na(y_iz[i, ]))) {
            S_zz <- I_zz * error2_i[i]
            jnll <- jnll - RTMB::dmvnorm(y_iz[i, ], x_iz[i, ], S_zz, TRUE)
        }
    }
    REPORT(gamma_iz)
    ADREPORT(Gsum_iz)
    REPORT(V_zz)
    jnll
}

obj <- MakeADFun(f, par, random = "x_iz")
obj$fn()
obj$gr()
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt # 313.7647 is book answer
sdr <- sdreport(obj)
sdr

#------------------------------------------------------------------------------
# drift as a function of time

formula <- ~ splines::bs(t_i, 5)
data$X_ij <- as.matrix(model.matrix(formula, data = data.frame("t_i" = t_i)))
par$beta_jz <- matrix(0, nrow = ncol(data$X_ij), ncol = 2)

# Refit model
obj <- MakeADFun(f, par, random = "x_iz")
obj$fn()
obj$gr()
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt # 285.5892 is book answer
sdr <- sdreport(obj)
sdr

#------------------------------------------------------------------------------
# Northern fur seal demo

DF <- read.csv("data/FSdata_2016.csv")
DF <- st_as_sf(DF, coords = c("longitude", "latitude"), crs = "+proj=longlat")
DF <- st_transform(DF, crs = "+proj=utm +datum=WGS84 +units=km +zone=2")

DF <- subset(DF, tripno == 1 & dbid == 818)
DF$duration <- c(NA, diff(as.POSIXlt(DF$gmt)))
DF$duration[-1] <- ifelse(DF$tripno[-1] == DF$tripno[-nrow(DF)] &
    DF$dbid[-1] == DF$dbid[-nrow(DF)], DF$duration[-1], NA)
DF$t_i <- c(0, cumsum(DF$duration[-1]))

# Build drift covariances
formula <- ~ splines::bs(t_i, 3)
X_ij <- as.matrix(model.matrix(formula, data = data.frame("t_i" = t_i)))

# Build inputs
data <- list(
    "y_iz" = st_coordinates(DF),
    "deltaT_i" = DF$duration,
    "error2_i" = (DF$error_semi_major / 1000)^2,
    "X_ij" = X_ij, "n_factors" = 0,
    "RAM" = matrix(nrow = 0, ncol = 4)
)
par <- list(
    "sigma2_z" = log(1),
    "x_iz" = data$y_iz,
    "beta_jz" = matrix(0, nrow = ncol(data$X_ij), ncol = 2)
)

# drop some data
which_include <- seq(1, nrow(data$y_iz), length = 20)
data$y_iz[-which_include, ] <- NA

# Refit model
obj <- MakeADFun(f, par, random = "x_iz")
obj$fn()
obj$gr()
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt # 133.0738 is book answer
sdr <- sdreport(obj)
sdr

#------------------------------------------------------------------------------
# Simulate multiple tracks

# Parameters
gamma_t <- beta_t <- seq(0.5, -0.5, length = n_t)
beta_t <- ifelse(beta_t > 0, 0, beta_t)
gamma_t <- ifelse(gamma_t < 0, 0, gamma_t)

# Simulate four tracks
n_tracks <- 4
set.seed(101)
x_iz <- NULL
for (track in 1:n_tracks) {
    s <- simulate_track(
        n_t = n_t, beta_t = beta_t,
        gamma_t = gamma_t, get_preference = get_preference
    )
    colnames(s) <- paste0(c("x", "y"), track)
    x_iz <- cbind(x_iz, s)
}
y_iz <- x_iz + array(rnorm(prod(dim(x_iz))), dim = dim(x_iz))

# Build inputs
data <- list(
    "y_iz" = y_iz,
    "deltaT_i" = c(NA, rep(1, nrow(y_iz) - 1)),
    "error2_i" = rep(1, nrow(y_iz)),
    "X_ij" = array(dim = c(nrow(y_iz), 0)),
    "n_factors" = 0,
    "RAM" = matrix(nrow = 0, ncol = 4)
)
par <- list(
    "sigma2_z" = log(1),
    "x_iz" = data$y_iz,
    "beta_jz" = array(0, dim = c(ncol(data$X_ij), ncol(data$y_iz)))
)

# drop some data
which_include <- seq(1, nrow(data$y_iz), length = 20 * n_tracks)
data$y_iz[-which_include, ] <- NA

# Refit model
obj <- MakeADFun(f, par, random = "x_iz")
obj$fn()
obj$gr()
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt # 1860.966 is book answer
sdr <- sdreport(obj)
sdr

#--------------------------------------------------------------------------------
# Switch to factor model

data$n_factors <- 2
par$sigma2_z <- rep(0.1, 2 * ncol(y_iz))
obj <- MakeADFun(f, par, random = "x_iz")
obj$fn()
obj$gr()
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt # 1790.638 is book answer
sdr <- sdreport(obj)
sdr

#--------------------------------------------------------------------------------
# Structural equation model

# Specify SEM
text <- "
  y1 -> y2, y
  y1 -> y3, y
  y1 -> y4, y
"

SEM_model <- sem::specifyModel(
    text = text,
    exog.variances = TRUE,
    endog.variances = TRUE,
    covs = colnames(y_iz)
)
RAM <- build_ram(SEM_model, colnames(y_iz)) # uses Jim's helper

# Build with RAM
data$n_factors <- 0
data$RAM <- as.matrix(RAM[, 1:4])
data$RAM[data$RAM[, 1] == 2, 4] <- 2
par$sigma2_z <- rep(0.1, max(data$RAM[, 4]))
obj <- MakeADFun(f, par, random = "x_iz")
obj$fn()
obj$gr()
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt # 1846.172 is book answer
sdr <- sdreport(obj)
obj$report()$V_zz
