library(sf)
library(RTMB)
library(tmbstan)

#--------------------------------------------------------------------------------
# load data, manipulate it for analysis
#--------------------------------------------------------------------------------

vismba <- readRDS("data/vismba.rds")
samples <- data.frame("x" = vismba$gx, "y" = vismba$gy, "agb" = vismba$agb)
samples <- st_as_sf(samples, coords = c("x", "y"))

grid <- st_make_grid(st_bbox(c(xmin = 0, xmax = 1000, ymin = 0, ymax = 500)),
    cellsize = c(250, 125)
)

grid_i <- st_intersects(samples, grid)

count_i <- tapply(samples$agb, INDEX = factor(unlist(grid_i),
    levels = 1:length(grid)
), FUN = length)

data <- data.frame(st_coordinates(st_centroid(grid)),
    "count" = ifelse(is.na(count_i), 0, count_i)
)

grid_sf <- st_sf(grid, count = count_i)

#--------------------------------------------------------------------------------
# fit with RTMB
#--------------------------------------------------------------------------------

data <- list("y_i" = data$count)
par <- list(
    "eps_i" = rep(0, length(data$y_i)),
    "ln_mu" = 0,
    "ln_sd" = 0
)

f <- function(par) {
    getAll(data, par, warn = FALSE)
    y_i <- OBS(y_i)
    jnll <- 0
    # Pr(random coefficients)
    jnll <- jnll - sum(dnorm(eps_i, ln_mu, exp(ln_sd), TRUE))
    # Pr(data | fixed + random effects):
    yhat_i <- exp(eps_i)
    jnll <- jnll - sum(dpois(y_i, yhat_i, TRUE))
    yhat_sum <- sum(yhat_i)
    REPORT(jnll)
    REPORT(yhat_i)
    REPORT(yhat_sum)
    ADREPORT(yhat_i)
    ADREPORT(yhat_sum)
    jnll
}

obj <- MakeADFun(f, par, random = "eps_i")
obj$fn()
obj$gr()
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt # book solution 63.20921

sdr <- sdreport(obj, getJointPrecision = TRUE)
parhat <- summary(sdr, "fixed")
sdr
parhat

#--------------------------------------------------------------------------------
# fit with tmbstan
#--------------------------------------------------------------------------------

fit <- tmbstan(obj, chains = 1)
plot(fit, pars = setdiff(names(fit), "lp__"))
