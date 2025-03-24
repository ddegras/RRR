# R package RRR: Reduced Rank Ridge Regression

The R package `RRR` provides functions to: 
-   Perform multivariate regression: reduced rank regression (RR), multivariate ridge regression, partial least squares regression (PLS), and principal component regression (PCR).
-   Simulate multivariate regression models: RR, PLS, PCR.
-   Select tuning parameters (ridge parameter and/or regression rank) by cross-validation.

It depends on the package `RSpectra`. To install the package: 
```
library(devtools)
# Uncomment the line below if needed
# install.packages("RSpectra")
install_github("ddegras/RRR")
```
