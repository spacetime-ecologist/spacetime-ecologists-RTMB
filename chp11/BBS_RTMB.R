library(sf)
library(ape)
library(rnaturalearth)
library(phylobase)
library(phylosignal)
library(viridisLite)
library(fmesher)
library(RTMB)

#--------------------------------------------------------------------------------
# load data, manipulate it for analysis
#--------------------------------------------------------------------------------

DF <- read.csv(file = "Top20_Samples.csv")
trait_set <- read.csv("Top20_Traits.csv")

# sf_df
sf_grid <- st_geometry(st_read("sf_grid.shp"))

# sf_grid
sf_DF <- read.csv("sf_DF.csv")
sf_DF <- st_as_sf(sf_DF,
    coords = c("X", "Y"),
    crs = st_crs("+proj=longlat +datum=WGS84")
)
sf_DF$Genus_species <- factor(sf_DF$Genus_species)

# df_grid
df_grid <- read.csv("df_grid.csv")
df_grid <- st_as_sf(df_grid,
    coords = c("X", "Y"),
    crs = st_crs("+proj=longlat +datum=WGS84")
)
taxa <- levels(sf_DF$Genus_species)

# create mesh
loc_DF <- st_coordinates(sf_DF)
loc_grid <- st_coordinates(st_centroid(sf_grid))
mesh <- fm_mesh_2d(loc_grid, refine = TRUE, cutoff = 1)

# create matrices in INLA
spde <- fm_fem(mesh, order = 2)
# create projection matrix from vertices to samples
A_is <- fm_evaluator(mesh, loc = loc_DF)$proj$A
# create projection matrix from vertices to grid
A_gs <- fm_evaluator(mesh, loc = loc_grid)$proj$A

# read tree
tree <- read.tree("Top20_Tree.tre")
tree$node.label <- c("root", seq_len(Nnode(tree) - 1))

# extract path from species to root
root <- Ntip(tree) + 1 # index for root
paths <- sapply(match(taxa, tree$tip.label),
    FUN = nodepath, phy = tree, to = root
)
names(paths) <- taxa

# convert to triplet form i=species; j=edge; x=edge_length
i <- unlist(lapply(seq_along(paths), FUN = function(i) {
    rep(i, length(paths[[i]]))
}))
j <- unlist(paths)
edge <- match(j, tree$edge[, 2])
x <- tree$edge.length[edge]

# remove root and renumber
ijx <- na.omit(cbind(i, j, x))
ijx[, 1:2] <- ifelse(ijx[, 1:2] > Ntip(tree), ijx[, 1:2] - 1, ijx[, 1:2])

# convert to dense design matrix
PhyloDesign_gc <- matrix(0, nrow = Nedge(tree), ncol = length(taxa))
PhyloDesign_gc[ijx[, 2:1]] <- ijx[, 3]

#--------------------------------------------------------------------------------
# alternatives to run
#--------------------------------------------------------------------------------

version <- c("Residual" = TRUE, "Phylo" = TRUE, "Traits" = TRUE, "Habitat" = TRUE)

# factors
if (version["Residual"]) {
    Lform_jc <- diag(1:nlevels(sf_DF$Genus_species))
} else {
    Lform_jc <- matrix(0, nrow = 0, ncol = length(taxa))
}

# phylogenetic design matrix
if (version["Phylo"]) {
    C_gc <- PhyloDesign_gc / mean(colSums(PhyloDesign_gc))
} else {
    C_gc <- matrix(0, nrow = 0, ncol = length(taxa))
}

# traits
if (version["Traits"] == TRUE) {
    T_ch <- model.matrix(
        ~ 0 + log(Mass) +
            factor(Primary.Lifestyle) +
            log(Hand.Wing.Index),
        data = trait_set[match(taxa, trait_set$Genus_species), ]
    )
    T_hc <- t(scale(T_ch))
} else {
    T_hc <- matrix(0, nrow = 0, ncol = length(taxa))
}

if (version["Habitat"] == TRUE) {
    habitat_formula <- ~ 0 + poly(log_elevation_km, 2, raw = TRUE) +
        poly(scale_NDVI, 2, raw = TRUE) +
        poly(log_pop_dens, 2, raw = TRUE)
    X_gk <- model.matrix(habitat_formula, data = df_grid)
    X_ik <- model.matrix(habitat_formula, data = sf_DF)
} else {
    X_gk <- model.matrix(~0, data = df_grid)
}

#--------------------------------------------------------------------------------
# fit in RTMB
#--------------------------------------------------------------------------------

rmatrix <- function(nrow, ncol, ...) {
    matrix(data = rnorm(n = nrow * ncol, ...), nrow = nrow, ncol = ncol)
}

set.seed(1) # ensure reproducibility of start values
par <- list(
    "ln_kappa" = log(1),
    "ln_sigmaM" = log(1),
    "ln_sigmaC" = rep(0, ifelse(nrow(C_gc) > 0, 1, 0)),
    "ln_sigmaT_h" = rep(0, nrow(T_hc)),
    "alpha_c" = rep(0, length(taxa)),
    "lambda_jc" = ifelse(Lform_jc == 0, 0, 1) *
        rmatrix(sd = 0.1, nrow = nrow(Lform_jc), ncol = ncol(Lform_jc)),
    "beta_kc" = rmatrix(sd = 0.1, nrow = ncol(X_ik), ncol = length(taxa)),
    "gamma_sh" = rmatrix(sd = 0.1, nrow = nrow(spde$c0), ncol = nrow(T_hc)),
    "delta_sg" = rmatrix(sd = 0.1, nrow = nrow(spde$c0), ncol = nrow(C_gc)),
    "omega_sj" = rmatrix(sd = 0.1, nrow = nrow(spde$c0), ncol = nrow(Lform_jc))
)

data <- list(
    "n_i" = sf_DF$SpeciesTotal,
    "c_i" = as.numeric(sf_DF$Genus_species),
    "C_gc" = C_gc,
    "T_hc" = T_hc,
    "X_gk" = X_gk,
    "X_ik" = X_ik,
    "M0" = spde$c0,
    "M1" = spde$g1,
    "M2" = spde$g2,
    "A_is" = A_is,
    "A_gs" = A_gs
)

f <- function(par) {
    getAll(data, par, warn = FALSE)
    jnll <- 0
    tau <- 1 / sqrt(4 * pi * exp(2 * ln_kappa))
    Q <- tau^2 * (exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2)
    # assemble predictions at data
    beta_ic <- X_ik %*% beta_kc
    gamma_ic <- A_is %*% gamma_sh %*% T_hc
    delta_ic <- A_is %*% delta_sg %*% C_gc
    omega_ic <- A_is %*% omega_sj %*% lambda_jc
    p_ic <- beta_ic + gamma_ic + delta_ic + omega_ic
    # assemble predictions at extrapolation grid
    beta_gc <- X_gk %*% beta_kc
    gamma_gc <- A_gs %*% gamma_sh %*% T_hc
    delta_gc <- A_gs %*% delta_sg %*% C_gc
    omega_gc <- A_gs %*% omega_sj %*% lambda_jc
    p_gc <- beta_gc + gamma_gc + delta_gc + omega_gc
    # Pr(phylogeny)
    for (g in 1:ncol(delta_sg)) {
        jnll <- jnll - dgmrf(
            delta_sg[, g], 0, Q,
            scale = exp(ln_sigmaC[1]),
            log = TRUE
        )
    }
    # Pr(traits)
    for (h in 1:ncol(gamma_sh)) {
        jnll <- jnll - dgmrf(
            gamma_sh[, h], 0, Q,
            scale = exp(ln_sigmaT_h[h]),
            log = TRUE
        )
    }
    # Pr(factors)
    for (j in 1:ncol(omega_sj)) {
        jnll <- jnll - dgmrf(omega_sj[, j], 0, Q, log = TRUE)
    }
    # Pr(data|random effects)
    for (i in 1:length(n_i)) {
        mu <- exp(alpha_c[c_i[i]] + p_ic[i, c_i[i]])
        jnll <- jnll - dnbinom2(n_i[i], mu, mu * (1 + exp(2 * ln_sigmaM)), TRUE)
    }
    REPORT(Q) # for variance calculations
    REPORT(beta_gc)
    REPORT(gamma_gc)
    REPORT(delta_gc)
    REPORT(omega_gc)
    REPORT(p_gc)
    ADREPORT(beta_kc) # get covariance even when REML
    jnll
}

map <- list("lambda_jc" = factor(ifelse(Lform_jc == 0, NA, Lform_jc)))
random <- c("gamma_sh", "delta_sg", "omega_sj")
# use restricted maximum likelihood (REML):
random <- c(random, "alpha_c", "beta_kc")

TapeConfig(atomic = "disable") # use if RTMB hangs doing matrix algebra
obj <- MakeADFun(f, par, random = random, map = map)
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt
# not run because it takes a while to run:
# sdr <- sdreport(obj, getReportCovariance = FALSE)
TapeConfig(atomic = "enable") # set it back
