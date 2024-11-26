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