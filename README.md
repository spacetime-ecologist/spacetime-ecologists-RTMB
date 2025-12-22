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

All working models can be executed by entering the spacetime-ecologists-RTMB folder and typing

```shell
make   
```
Each model's fitted negative log likelihood is stored and can be inspected via

```shell
cat RTMB.res
```
Note that all packages required to run the individual chapter scripts must be installed prior to using the `Makefile`, which can be installed by sourcing the `install.R` file

## Funding acknowledgements

Cahill's effort was funded by supporting partners of the Quantitative Fisheries Center, 
which includes Michigan State University, Michigan Department of Natural Resources (DNR), 
the Great Lakes Fishery Commission, and multiple Council of Lake Committee Agencies including the Michigan DNR, Ohio DNR, Minnesota DNR, 
Illinois DNR, Wisconsin DNR, Pennsylvania Fish and Game Commission, 
New York State Department of Environmental Conservation, 
Ontario Ministry of Natural Resources, and the Great Lakes Indian Fish and Wildlife Commission. 
