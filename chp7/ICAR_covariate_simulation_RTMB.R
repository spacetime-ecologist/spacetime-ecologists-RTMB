library(sf)
library(Matrix)
library(RTMB)
# note: also requires 'igraph' package
source("shared_functions/rmvnorm_prec.R")
set.seed(1) # needed to match book solution, which otherwise is not reproducible

#--------------------------------------------------------------------------------
# read in data and manipulate for analysis
#--------------------------------------------------------------------------------

sf_area <- rnaturalearth::ne_countries(
    scale = 10,
    country = "Finland", return = "sf"
)

# extract main island
sf_area <- st_cast(st_geometry(sf_area), "POLYGON")[[1]]
sf_area <- st_sfc(sf_area, crs = "+proj=longlat +datum=WGS84")

# make intro grid
cellsize <- 0.5
sf_fullgrid <- st_make_grid(sf_area, cellsize = cellsize)
sf_grid <- st_make_valid(st_intersection(sf_fullgrid, sf_area))
sf_grid <- sf_grid[st_area(sf_grid) > (0.2 * max(st_area(sf_grid)))]

# get adjacency
st_rook <- function(a, b = a, ...) st_relate(a, b, pattern = "F***1****", ...)
grid_A <- st_rook(sf_grid, sparse = TRUE)
A_ss <- as(grid_A, "sparseMatrix")
A_ss <- as(A_ss, "TsparseMatrix")
I_ss <- Matrix::sparseMatrix(
    i = 1:length(sf_grid),
    j = 1:length(sf_grid),
    x = rep(1, length(sf_grid))
)
I_ss <- as(I_ss, "TsparseMatrix")

# using igraph for sparse-matrix calculations
graph <- igraph::graph_from_adjacency_matrix(A_ss, weighted = TRUE)
f2 <- function(x, extra = NULL) {
    as.vector(x %*% A_ss)
}
rho_min <- 1 / Re(igraph::arpack(f2,
    sym = FALSE,
    options = list(n = nrow(A_ss), nev = 3, ncv = 8, which = "SR")
)$values[1])
rho_max <- 1 / Re(igraph::arpack(f2,
    sym = FALSE,
    options = list(n = nrow(A_ss), nev = 3, ncv = 8, which = "LR")
)$values[1])

# make ICAR precision
rho <- rho_max * 0.75
var_X1 <- 0.5^2
var_X2 <- 0.5^2
Q <- diag(rep(1, nrow(A_ss))) - rho * A_ss
X1 <- rmvnorm_prec(mu = rep(0, nrow(A_ss)), prec = Q / var_X1, n.sims = 1)
X2 <- rmvnorm_prec(mu = rep(0, nrow(A_ss)), prec = Q / var_X2, n.sims = 1)
lambda_s <- exp(1 + X1 + X2)

n_obs <- nrow(A_ss)
s_i <- sample(1:nrow(A_ss), size = n_obs, replace = FALSE)
c_i <- rpois(n_obs, lambda_s[s_i])

formula <- ~0
X_sk <- model.matrix(formula, data.frame("X1" = X1, "X2" = X2))

#--------------------------------------------------------------------------------
# fit in RTMB
#--------------------------------------------------------------------------------

gamma_k <- ifelse(ncol(X_sk) == 0, numeric(1), numeric(ncol(X_sk)))
par <- list(
    "beta0" = 0,
    "rho_prime" = 0,
    "log_sig" = 0,
    "gamma_k" = gamma_k,
    "omega_s" = rep(0, nrow(A_ss))
)

data <- list(
    "c_i" = c_i,
    "rho_bounds" = c(rho_min, rho_max),
    "s_i" = s_i,
    "X_sk" = rep(0, ncol(X_sk)),
    "I_ss" = I_ss,
    "A_ss" = A_ss
)

f <- function(par) {
    getAll(data, par, warn = FALSE)
    rho <- (1 - plogis(rho_prime)) * (rho_bounds[2] - rho_bounds[1]) + rho_bounds[1]
    jnll <- 0
    Q_ss <- (I_ss - rho * A_ss) / exp(2 * log_sig) # note this is scaled
    jnll <- jnll - dgmrf(omega_s, 0.0, Q_ss, TRUE)
    # control situations when X_sk is NULL, which I think is automatically done in TMB:
    if (is.null(dim(X_sk))) {
        lambda_s <- exp(beta0 + omega_s)
    } else {
        lambda_s <- exp(beta0 + omega_s + X_sk * gamma_k)
    }
    # Pr(data|fixed and random coefficients):
    for (i in 1:length(c_i)) {
        jnll <- jnll - dpois(c_i[i], lambda_s[s_i[i]], TRUE)
    }
    sum_lambda <- sum(lambda_s)
    ADREPORT(sum_lambda)
    REPORT(rho)
    jnll
}

obj <- MakeADFun(f, par,
    random = "omega_s",
    map = list("gamma_k" = factor(NA)) # map this off to match book
)

opt <- nlminb(obj$par, obj$fn, obj$gr)
opt # --> 620.538 is book solution if set.seed(1)
sdr <- sdreport(obj)
sdr
