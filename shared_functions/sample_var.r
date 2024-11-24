# function to generate samples
sample_var <- function(obj, # TMB object
                       var_name, # name of REPORTed variable
                       mu, # Estimated mean
                       prec, # Estimated precision
                       n_samples = 500, # Number of samples
                       fun1 = function(x, ...) x, # optional samples transform
                       fun2 = sd, # calculate fun2 across samples
                       ...) {
    require(abind)
    # Take samples using precision matrix
    u_zr <- rmvnorm_prec(mu = mu, prec = prec, n.sims = n_samples)
    # Calculate REPORTed variable for each sample
    for (rI in 1:n_samples) {
        Var <- obj$report(par = u_zr[, rI])[[var_name]]
        # Transform REPORTed variable using fun1
        Var <- fun1(Var, ...)
        # Bind samples together into array
        if (is.vector(Var)) Var <- as.array(Var)
        if (rI == 1) Var_zr <- Var
        if (rI >= 2) {
            Var_zr <- abind(Var_zr, Var, along = length(dim(Var)) + 1)
        }
    }
    # summarize across samples using fun2
    out <- apply(Var_zr, MARGIN = 1:(length(dim(Var_zr)) - 1), FUN = fun2)
    # Return value
    return(out)
}
