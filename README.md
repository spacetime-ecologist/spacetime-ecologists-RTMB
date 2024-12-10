## Thorson and Kristensen 2024 in RTMB

RTMB code translations of TMB code associated with textbook:

Thorson, J., and Kristensen, K. 2024. Spatio-Temporal Models for Ecologists. In 1st edition. Chapman and Hall/CRC, Boca Raton, FL.

[https://github.com/James-Thorson/Spatio-temporal-models-for-ecologists](https://github.com/James-Thorson/Spatio-temporal-models-for-ecologists)

All RTMB examples coded by Christopher Cahill (Quantitative Fisheries Center at Michigan State University) with assistance from Kasper Kristensen, Anders Nielsen, and Christoffer Albertsen (DTU).  

See also RTMB:

[https://github.com/kaskr/RTMB](https://github.com/kaskr/RTMB)

## Run RTMB code for a single chapter

Enter a given chapter's folder (i.e., set the working directory or open a project in a chapter folder) and run the scripts ending with '_RTMB.R'.  

## Run all RTMB models in the book (Linux)

All working RTMB models can be executed by entering the ST_for_ecologists_RTMB folder and typing

```shell
make -j 4 -l 2
```
which runs in approximately 17 minutes on a laptop with 64 GB RAM and 20 cores.  Each RTMB model's fitted negative log likelihood can then be inspected via

```shell
cat RTMB.res
```
Note that all packages required to run the individual chapter scripts must be installed prior to using the `Makefile` 

## Known issues

- chapter 10 models do not yet work (in progress)
- chapter 11 `BBS.R` from book does not run, errors out on line 42

