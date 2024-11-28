library(terra)
library(sf)
library(Matrix)
library(rnaturalearth)
library(ggplot2)
library(marginaleffects)
library(RTMB)

#--------------------------------------------------------------------------------
# read in data and format for analysis
#--------------------------------------------------------------------------------

# load data
bathy_terra <- readRDS("data/ai_bathy_3km.Rds")
likelihood_terra <- readRDS("data/likelihood_3km.Rds")

# get land layer
sf_states <- ne_states("united states of america", return = "sf")
sf_states <- sf_states[pmatch("Alas", sf_states$name_en), ]
sf_states <- st_transform(sf_states,
    crs = st_crs("+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
)
sf_states <- st_geometry(sf_states)

# change resolution
fact <- 1
bathy_terra <- aggregate(bathy_terra, fact = fact)
likelihood_terra <- as.array(aggregate(likelihood_terra, fact = fact))

# convert bathymetry
bathy_sf <- st_as_sf(as.polygons(bathy_terra, trunc = FALSE, dissolve = FALSE))
bathy_sf$northings <- st_coordinates(st_centroid(bathy_sf))[, "Y"]
bathy_sf$eastings <- st_coordinates(st_centroid(bathy_sf))[, "X"]

# format likelihood
L_gt <- apply(likelihood_terra, MARGIN = 3, FUN = function(mat) {
    as.vector(t(mat))
})
sf_L_gt <- st_sf(st_geometry(bathy_sf), L_gt)

# get adjacency matrix using Raster
A_gg <- adjacent(bathy_terra, cells = 1:prod(dim(bathy_terra)), pairs = TRUE)
A_gg <- Matrix::sparseMatrix(i = A_gg[, 1], j = A_gg[, 2], x = rep(1, nrow(A_gg)))
A_gg <- as(A_gg, "TsparseMatrix")

# drop geometry and likelihood with depth <1 m
which_exclude <- which(bathy_sf$ai_bathy_fill <= 1)
sf_L_gt <- sf_L_gt[-which_exclude, ]
bathy_sf <- bathy_sf[-which_exclude, ]
A_gg <- A_gg[-which_exclude, -which_exclude]
bathy_sf$ai_bathy_fill <- bathy_sf$ai_bathy_fill / 1000

# assemble inputs
colsumA_g <- colSums(A_gg)
I_gg <- Matrix::sparseMatrix(
    i = 1:nrow(bathy_sf),
    j = 1:nrow(bathy_sf), x = rep(1, nrow(bathy_sf))
)
I_gg <- as(I_gg, "TsparseMatrix")
D_gg <- Matrix::sparseMatrix(
    i = 1:nrow(bathy_sf),
    j = 1:nrow(bathy_sf), x = colsumA_g
)
D_gg <- as(D_gg, "TsparseMatrix")
At_zz <- cbind(attr(A_gg, "i"), attr(A_gg, "j"))

# Make covariates
preference_formula <- ~ 0 + poly(ai_bathy_fill, 2)
X_gz <- model.matrix(preference_formula, data = bathy_sf)

#--------------------------------------------------------------------------------
# estimate using RTMB
#--------------------------------------------------------------------------------

set.seed(1) # --> ensure reproducibility for re starting values
par <- list(
    "ln_D" = 1,
    "gamma_z" = 0.01 * rnorm(ncol(X_gz))
)

data <- list(
    "CTMC_version" = 1, # 0=diffusion-taxis;  1=logspace diffusion-taxis
    "expm_version" = 0, # 0=uniformization;  1=series
    "Nmax" = ifelse(mean(res(bathy_terra)) == 3000, 500, 250),
    "DeltaD" = mean(res(bathy_terra)) / 1000,
    "colsumA_g" = colsumA_g,
    "L_gt" = as.matrix(st_drop_geometry(sf_L_gt)),
    "X_gz" = X_gz,
    "At_zz" = At_zz + 1 # for 1-based indexing
)

make_M <- function(CTMC_version, n_g, DeltaD, At_zz, ln_D, h_g, colsumA_g) {
    n_z <- nrow(At_zz)
    D <- exp(ln_D)
    Mrate_gg <- AD(Matrix(0, nrow = n_g, ncol = n_g))
    ones <- matrix(1, ncol = 1, nrow = n_g)
    # standard approach
    if (CTMC_version == 0) {

    }
    # log-space to ensure Metzler matrix
    if (CTMC_version != 0) {
        # combined taxis and diffusion
        Mrate_gg[At_zz] <- Mrate_gg[At_zz] +
            D / DeltaD^2 * exp((h_g[At_zz[, 2]] - h_g[At_zz[, 1]]) / DeltaD)
        row_sums <- Mrate_gg %*% ones
        diag(Mrate_gg) <- diag(Mrate_gg) - as.vector(row_sums)
    }
    Mrate_gg
}

f <- function(par) {
    getAll(data, par, warn = FALSE)
    jnll <- 0
    n_g <- nrow(L_gt)
    n_t <- ncol(L_gt)
    forward_prob_gt <- forward_pred_gt <- forward_gt <-
        backward_prob_gt <- backward_pred_gt <-
        backward_gt <- matrix(0, n_g, n_t)

    # calculate movement matrix
    h_g <- X_gz %*% gamma_z
    Mrate_gg <- make_M(CTMC_version, n_g, DeltaD, At_zz, ln_D, h_g, colsumA_g)
    str(Mrate_gg)
    REPORT(Mrate_gg)
    jnll
}

f(par) # works
# build and optimize object
AD(Matrix(0, nrow = 5, ncol = 5))
obj <- MakeADFun(f, par)

head(obj$report()$`Mrate_gg`)
# opt <- nlminb(start = obj$par, obj = obj$fn, gr = obj$gr)
# opt # 161.4145 is book solution

head(Mrate_gg)

#  # diag_g <- diag(Mrate_gg)
rho <- 0
if (length(diag_g) > 0) {
    M <- diag_g[1]
    for (i in 1:length(diag_g)) {
        # M <- min(M, diag_g[i]) # ! NOTE MINIMUM BREAKS RTMB
    }
    rho <- rho - M
}
A_gg <- Mrate_gg
diag(A_gg) <- rho
A_prime_gg <- t(A_gg)

# ever so slight difference Mrate_gg --> rounding errors?
# explain AD(Matrix)
