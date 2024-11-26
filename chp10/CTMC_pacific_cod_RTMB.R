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

bathy_terra <- readRDS("data/ai_bathy_3km.Rds")
likelihood_terra <- readRDS("data/likelihood_3km.Rds")

# Get land layer
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

# Convert bathymetry
bathy_sf <- st_as_sf(as.polygons(bathy_terra, trunc = FALSE, dissolve = FALSE))
bathy_sf$northings <- st_coordinates(st_centroid(bathy_sf))[, "Y"]
bathy_sf$eastings <- st_coordinates(st_centroid(bathy_sf))[, "X"]

# Format likelihood
L_gt <- apply(likelihood_terra, MARGIN = 3, FUN = function(mat) {
    as.vector(t(mat))
})
sf_L_gt <- st_sf(st_geometry(bathy_sf), L_gt)

# Get adjacency matrix using Raster
A_gg <- adjacent(bathy_terra, cells = 1:prod(dim(bathy_terra)), pairs = TRUE)
A_gg <- Matrix::sparseMatrix(i = A_gg[, 1], j = A_gg[, 2], x = rep(1, nrow(A_gg)))
A_gg <- as(A_gg, "TsparseMatrix")

# Drop geometry and likelihood with depth <1 m
which_exclude <- which(bathy_sf$ai_bathy_fill <= 1)
sf_L_gt <- sf_L_gt[-which_exclude, ]
bathy_sf <- bathy_sf[-which_exclude, ]
A_gg <- A_gg[-which_exclude, -which_exclude]
bathy_sf$ai_bathy_fill <- bathy_sf$ai_bathy_fill / 1000

# Assemble inputs
colsumA_g <- colSums(A_gg)
I_gg <- Matrix::sparseMatrix(i = 1:nrow(bathy_sf), j = 1:nrow(bathy_sf), x = rep(1, nrow(bathy_sf)))
I_gg <- as(I_gg, "TsparseMatrix")
D_gg <- Matrix::sparseMatrix(i = 1:nrow(bathy_sf), j = 1:nrow(bathy_sf), x = colsumA_g)
D_gg <- as(D_gg, "TsparseMatrix")
At_zz <- cbind(attr(A_gg, "i"), attr(A_gg, "j"))

# Make covariates
preference_formula <- ~ 0 + poly(ai_bathy_fill, 2)
X_gz <- model.matrix(preference_formula, data = bathy_sf)

#--------------------------------------------------------------------------------
# estimate using RTMB
#--------------------------------------------------------------------------------

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
    "At_zz" = At_zz
)
make_M <- function(CTMC_version, n_g, DeltaD, At_zz, ln_D, h_g, solsumA_g) {
  n_z <- nrow(At_zz)
  D <- exp(ln_D)
  if(CTMB_version==){
    
  }
}
f <- function(par) {
    getAll(data, par, warn = FALSE)
    jnll <- 0
    n_g <- nrow(L_gt)
    n_t <- ncol(L_gt)
    forward_prob_gt <- forward_pred_gt <- forward_gt <-
        backward_prob_gt <- backward_pred_gt <- backward_gt <- matrix(0, n_g, n_t)

    # calculate movement matrix
    h_g <- X_gz %*% gamma_z
    ...
    ...
    ...
}

# Build and optimize object
obj <- MakeADFun(data = Data, parameters = Params)
opt <- nlminb(start = obj$par, obj = obj$fn, gr = obj$gr)
opt # 161.4145 is book solution
