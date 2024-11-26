rmvnorm_prec <-
    function(mu, # estimated fixed and random effects
             prec, # estimated joint precision
             n.sims) {
        require(Matrix)
        # Simulate values
        z0 <- matrix(rnorm(length(mu) * n.sims), ncol = n.sims)
        # Q = t(P) * L * t(L) * P
        L <- Cholesky(prec, super = TRUE)
        # Calcualte t(P) * solve(t(L)) * z0 in two steps
        z <- solve(L, z0, system = "Lt") # z = Lt^-1 * z
        z <- solve(L, z, system = "Pt") # z = Pt    * z
        mu + as.matrix(z)
    }
