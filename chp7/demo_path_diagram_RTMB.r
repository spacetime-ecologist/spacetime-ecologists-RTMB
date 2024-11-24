library(mvtnorm)
library(stars)
library(diagram)
library(sem)
library(RTMB)
source("shared_functions/build_ram.r")
source("shared_functions/make_covar.r")

#--------------------------------------------------------------------------------
# simulate data and define SEM
#--------------------------------------------------------------------------------

names <- c("Tree cover", "Temperature", "log(Density)")
box.type <- rep("square", length(names))

# construct M
M2 <- M1 <- array(0,
    dim = c(
        length(names),
        length(names)
    ),
    dimnames = list(names, names)
)
M1["Temperature", "Tree cover"] <- -1
M1["log(Density)", "Temperature"] <- 1

set.seed(101)
n_obs <- 100
C <- rnorm(n_obs, mean = 0, sd = 1)
T <- M1["Temperature", "Tree cover"] * C + rnorm(n_obs, mean = 0, sd = 0.1)
logD <- 0 + M1["log(Density)", "Temperature"] * T
N <- rpois(n_obs, exp(logD))

# define model and convert to RAM
text <- "
  C -> T, b1
  T -> logD, b2
  logD <-> logD, NA, 0.01
"

SEM_model <- sem::specifyModel(
    text = text, exog.variances = TRUE,
    endog.variances = TRUE, covs = c("C", "T", "logD")
)

RAM <- build_ram(SEM_model, c("C", "T", "logD"))

#--------------------------------------------------------------------------------
# estimate using RTMB
#--------------------------------------------------------------------------------

data <- list(
    y_iz = cbind(C, T, N),
    RAM = as.matrix(RAM[, 1:4]),
    familycode_z = c(0, 0, 1),
    RAMstart = ifelse(is.na(RAM[, 5]), 0, as.numeric(RAM[, 5]))
)

par <- list(
    x_iz = data$y_iz,
    beta_j = rep(1, max(data$RAM[, 4]))
)

f <- function(par) {
    getAll(data, par, warn = FALSE)
    n_i <- nrow(y_iz)
    n_z <- ncol(y_iz)
    # define process variance:
    V_zz <- make_covar(beta_j, RAM, RAMstart, n_z, as.integer(0))
    jnll <- 0
    # Pr(random coefficients)
    for (i in 1:n_i) {
        jnll <- jnll - RTMB::dmvnorm(x_iz[i, ], 0, V_zz, TRUE)
    }

    # Pr(data|fixed,random par)
    for (i in 1:n_i) {
        for (z in 1:n_z) {
            if (familycode_z[z] == 1) {
                jnll <- jnll - dpois(y_iz[i, z], exp(x_iz[i, z]), TRUE)
            }
        }
    }
    REPORT(V_zz)
    jnll
}

map <- list("x_iz" = factor(cbind(NA, NA, 1:nrow(data$y_iz))))
obj <- MakeADFun(f, par, random = "x_iz", map = map)
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt # book solves to 176.9859
sdr <- sdreport(obj)
sdr
