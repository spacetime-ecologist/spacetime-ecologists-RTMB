##################################
##################################
##################################
##################################
##################################
# ! WARNING - not yet finished
# ! Confusing - inner hessian failure
# ! with ln_D_st...
##################################
##################################
##################################
##################################
##################################

library(RTMB)
library(sf)
library(Matrix)
library(raster)
library(elevatr)
library(rnaturalearth)
library(ggplot2)
library(gridExtra)

#--------------------------------------------------------------------------------
# read in data and format for analysis
#--------------------------------------------------------------------------------

# get spatial domain
sf_states <- ne_states(c("United States of America", "Canada"), return = "sf")
sf_states <- sf_states[pmatch(c(
    "Brit", "Alas", "Yukon",
    "Washington", "Ore", "Calif"
), sf_states$name_en), ]

p <- st_polygon(list(matrix(c(
    -180, 0, -180, 80, -50,
    80, -50, 0, -180, 0
), byrow = TRUE, ncol = 2)))
sf_states <- st_intersection(st_sfc(p,
    crs = "+proj=longlat +datum=WGS84"
), sf_states)
sf_states <- st_union(sf_states)
sf_states <- st_transform(sf_states, crs = st_crs("+init=epsg:3338 +units=km"))

# get coastline
sf_coast <- ne_coastline(scale = 50, return = "sf")
sf_coast <- st_transform(sf_coast, crs = st_crs("+init=epsg:3338 +units=km"))

# keep mainland Alaska
sf_states <- st_cast(sf_states, "POLYGON")
which_max <- which.max(st_area(sf_states))
sf_states <- sf_states[which_max]

# simplify geometry somewhat
sf_states <- st_simplify(sf_states, dTolerance = 10)

# load copNDVI (saved from rasterdiv)
copNDVI <- raster("data/NDVI.tif")

samples <- st_read("data/samples_3520.csv",
    options = c("X_POSSIBLE_NAMES=X", "Y_POSSIBLE_NAMES=Y")
)
st_crs(samples) <- "+proj=longlat +datum=WGS84"
samples <- st_transform(samples, crs = st_crs("+init=epsg:3338 +units=km"))
samples$Year <- as.numeric(samples$Year)
samples$SpeciesTotal <- as.numeric(samples$SpeciesTotal)
cellsize <- 200

df_grid <- st_read("data/df_grid.shp")
sf_grid <- st_geometry(st_read("data/sf_grid.shp"))

# bin into grids
samples <- st_intersection(samples, sf_grid)
grid_i <- as.integer(st_intersects(samples, sf_grid))

# get adjacency
st_rook <- function(a, b = a, ...) st_relate(a, b, pattern = "F***1****", ...)
grid_A <- st_rook(sf_grid, sparse = TRUE)
A_ss <- as(grid_A, "sparseMatrix")
A_ss <- as(A_ss, "TsparseMatrix")
I_ss <- Matrix::sparseMatrix(
    i = 1:length(sf_grid),
    j = 1:length(sf_grid), x = rep(1, length(sf_grid))
)
I_ss <- as(I_ss, "TsparseMatrix")
At_zz <- cbind(attr(A_ss, "i"), attr(A_ss, "j"))
colsumA_s <- colSums(A_ss)

#--------------------------------------------------------------------------------
# estimate it with RTMB
#--------------------------------------------------------------------------------

preference_formula <- ~ 0 + poly(elevatn, 2, raw = TRUE) +
    poly(NDVI, 2, raw = TRUE) + poly(dst_t_c, 2, raw = TRUE)
X_sz <- model.matrix(preference_formula, data = df_grid)

data <- list(
    "CTMC_version" = 1,
    "DeltaD" = cellsize,
    "n_t" = diff(range(samples$Year)) + 1,
    "n_s" = nrow(X_sz),
    "c_i" = samples$SpeciesTotal,
    "s_i" = grid_i - 1,
    "t_i" = samples$Year - min(samples$Year),
    "X_sz" = X_sz,
    "colsumA_s" = colsumA_s,
    "At_zz" = At_zz + 1,
    "rho_bounds" = 1 / range(eigen(as.matrix(grid_A))$values),
    "I_ss" = I_ss,
    "A_ss" = A_ss
)
par <- list(
    "rho_prime" = 0,
    "ln_sigmaO" = 0,
    "ln_sigmaB" = 0,
    "ln_D" = log(80^2 / (2 * 2 * 5)),
    "gamma_z" = rep(0, ncol(data$X_sz)),
    "beta_t" = rep(0, data$n_t),
    "ln_D_st" = matrix(0, nrow = data$n_s, ncol = data$n_t)
)

make_M <- function(CTMC_version, n_g, DeltaD, At_zz, ln_D, h_g, colsumA_g) {
    n_z <- nrow(At_zz)
    D <- exp(ln_D)
    # Mrate_gg <- AD(Matrix(0, nrow = n_g, ncol = n_g)) # sparse
    Mrate_gg <- matrix(0, nrow = n_g, ncol = n_g) # dense matrix
    # standard approach
    if (CTMC_version == 0) {
        # ! TODO
    }
    # log-space to ensure Metzler matrix
    if (CTMC_version != 0) {
        # combined taxis and diffusion
        Mrate_gg[At_zz] <- Mrate_gg[At_zz] +
            D / DeltaD^2 * exp((h_g[At_zz[, 2]] - h_g[At_zz[, 1]]) / DeltaD)
        row_sums <- rowSums(Mrate_gg)
        diag(Mrate_gg) <- diag(Mrate_gg) - as.vector(row_sums)
    }
    Mrate_gg
}

f <- function(par) {
    getAll(data, par, warn = FALSE)
    "[<-" <- ADoverload("[<-")
    rho <- plogis(rho_prime) * (rho_bounds[2] - rho_bounds[1]) + rho_bounds[1]
    ln_Dhat_st <- matrix(0, nrow = n_s, ncol = n_t)
    jnll <- 0
    Q_ss <- (I_ss - rho * A_ss) / exp(2 * ln_sigmaO)
    # assemble movement
    h_s <- X_sz %*% gamma_z
    Mrate_ss <- make_M(CTMC_version, n_s, DeltaD, At_zz, ln_D, h_s, colsumA_g)
    M_ss <- expm(Mrate_ss)
    ln_Dhat_st[, 1] <- beta_t[1]
    # initial density
    jnll <- jnll - dgmrf(x = ln_D_st[, 1], mu = ln_Dhat_st[, 1], Q = Q_ss, log = TRUE)
    # project density forward
    tmp_s <- numeric(n_s)
    for (t in 2:n_t) {
        tmp_s <- exp(ln_D_st[, t - 1])
        tmp_s <- as.vector(t(tmp_s) %*% M_ss)
        ln_D_st[, t] <- beta_t[t] + log(tmp_s)
        jnll <- jnll - dgmrf(x = ln_D_st[, t], mu = ln_Dhat_st[, t], Q = Q_ss, log = TRUE)
    }
    # project density w/o process errors
    proj_st <- AD(Matrix(0, n_s, n_t))
    proj_st[, 1] <- ln_D_st[, 1]
    for (t in 2:n_t) {
        proj_st[, t] <- t(proj_st[, t - 1]) %*% M_ss
    }
    # Pr(random coefficients)
    for (t in 2:n_t) {
        jnll <- jnll - dnorm(beta_t[t], beta_t[t - 1], exp(ln_sigmaB), TRUE)
    }
    # Pr(data|fixed and random coefs)
    for (i in 1:length(c_i)) {
        jnll <- jnll - dpois(c_i[i], exp(ln_D_st[s_i[i] + 1, t_i[i] + 1]), TRUE)
    }
    REPORT(jnll)
    jnll
}

f(par)
TapeConfig(atomic = "disable")
TMB::config(tmbad.sparse_hessian_compress = TRUE) # ??
map <- list(ln_D = factor(NA))
obj <- MakeADFun(f, par, random = c("beta_t", "ln_D_st"), map = map)

image(obj$env$spHess(random = TRUE)) # ! something very wrong with inner hessian ln_D_st

# opt <- nlminb(obj$par, obj$fn, obj$gr) # --> book solves to 8833.474
