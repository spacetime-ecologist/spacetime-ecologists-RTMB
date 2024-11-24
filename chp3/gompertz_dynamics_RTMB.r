library(RTMB)
library(diagram)
# ! Note also needs Matrix library
# ! NOTE, I can only get SAR Gompertz from book to match if I solve on logit_rho

#--------------------------------------------------------------------------------
# read in data, manipulate for analysis
#--------------------------------------------------------------------------------

index_t <- read.csv(file = "data/Biomass_index.csv")
index_t <- array(index_t[, 2], dimnames = list(index_t[, 1]))
brange <- range(index_t) * c(0.8, 1.2)
prange <- c(0.5, 2)
y <- log(index_t[-1])
x <- log(index_t[-length(index_t)])
xpred <- seq(log(brange[1]), log(brange[2]), length = 1000)

#--------------------------------------------------------------------------------
# fit CAR Gompertz using RTMB
#--------------------------------------------------------------------------------

data <- list(
    "log_b_t" = log(index_t),
    "log_bnew_z" = xpred,
    "simulate_t" = rep(0, length(index_t))
)

par <- list(
    "log_d0" = 0,
    "log_sdp" = 1, "log_sdo" = 1,
    "alpha" = 0,
    "rho" = 0,
    "log_d_t" = rep(0, length(index_t))
)

f <- function(par) {
    getAll(data, par, warn = FALSE)
    log_b_t <- OBS(log_b_t)
    jnll <- 0
    jnll <- jnll - dnorm(log_d_t[1], log_d0, exp(log_sdp), TRUE) # initialize
    for (t in 2:length(log_b_t)) { # Pr(random coefficients)
        jnll <- jnll - dnorm(
            log_d_t[t],
            alpha + rho * log_d_t[t - 1], exp(log_sdp),
            TRUE
        )
    }
    # Pr(data|fixed + random)
    jnll <- jnll - sum(dnorm(log_b_t, log_d_t, exp(log_sdo), TRUE))
    jnll
}

obj <- MakeADFun(f, par, random = "log_d_t")
obj$fn()
obj$gr()
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt # book solution -13.77209
sdr <- sdreport(obj)
sdr
obj$simulate()

#--------------------------------------------------------------------------------
# fit SAR Gompertz using RTMB
#--------------------------------------------------------------------------------

data <- list(
    "log_b_t" = log(index_t),
    "log_bnew_z" = xpred
)

par <- list(
    "log_delta" = 0,
    "log_sdp" = 1,
    "log_sdo" = 1,
    "alpha" = 0,
    "logit_rho" = 0,
    "eps_t" = rep(0, length(index_t))
)

Q_ar1 <- function(rho, n_t) {
    Q <- matrix(0, n_t, n_t)
    for (t in 1:n_t) {
        Q[t, t] <- 1 + rho^2
        if (t >= 2) Q[t - 1, t] <- -rho
        if (t >= 2) Q[t, t - 1] <- -rho
    }
    Q <- Q / (1 - rho^2) # scale
    Q <- as(Q, "sparseMatrix")
    Q
}

f <- function(par) {
    getAll(data, par)
    log_b_t <- OBS(log_b_t)
    n_t <- length(log_b_t)
    rho <- 1 - plogis(logit_rho)
    jnll <- 0
    # Pr(random coefficients)
    Q <- Q_ar1(rho, n_t)
    jnll <- jnll - dgmrf(eps_t, 0.0, Q, TRUE)
    log_d_t <- numeric(n_t)
    for (t in 1:n_t) {
        log_d_t[t] <- log_delta * rho^t + alpha / (1 - rho) + exp(log_sdp) * eps_t[t]
    }
    # Pr(data|fixed + random effects):
    jnll <- jnll - sum(dnorm(log_b_t, log_d_t, exp(log_sdo), TRUE))
    jnll
}

obj <- MakeADFun(f, par, random = "eps_t")
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt # book solution -13.91269 when I solve book version on logit_rho
sdr <- sdreport(obj)
sdr
obj$simulate()
