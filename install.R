# install all packages needed for the book
pkgs <- readLines("packages")
installed <- rownames(installed.packages())
install.packages(setdiff(pkgs, installed))

