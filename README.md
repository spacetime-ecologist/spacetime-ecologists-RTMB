# Convert Thorson and Kristensen 2024 TMB examples to RTMB

## Issues
            
## Chapter 5 

    - `spatial_GLMM_RTMB.r` 
    - joint distribution (ctl = 2) version of this model hangs at `makeADFun()`

## Chapter 7 

    - `ICAR_covariate_simulation_RTMB.r` does not solve to same jnll as example in book, so may be a bug on my end
    - `Integrated_model.R` TMB model from book does not converge so not yet attempted in RTMB

## Chapter 9 

    - `sea_ice_RTMB.r` hangs when calling `makeADFun()`

## Chapter 10

    - 'CTMC_eagles.R` example from book doesn't run for me, errors out on line 56

## Chapter 11 

    - `BBS.R` example from book doesn't run for me, errors out on line 42