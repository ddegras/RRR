\name{kcv}
\alias{kcv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
K-Fold Cross-Validation 
}
\description{
Select tuning parameters for multivariate regression methods by k-fold cross-validation.  
}

\usage{
kcv(x, y, lambda = NULL, r = NULL, nfolds = 5, method = c("rrr", "rr", "ridge", "pls", "pcr"), seed = NULL)
}

\arguments{
  \item{x}{
predictor matrix (cases-by-predictors)
}
  \item{y}{
response matrix (cases-by-responses)
}
  \item{lambda}{
smoothing parameters to use in ridge and reduced-rank ridge regression
}
  \item{r}{
ranks to use in all methods but ridge regression
}
  \item{nfolds}{
number of cross-validation folds
}
  \item{method}{
regression methods to use (can be abbreviated)
}
  \item{seed}{
optional seed for random number generator. Useful for reproducibility of results
}
}
\details{
\code{rrr} = reduced rank ridge regression, \code{rr} = reduced rank regression
}


\value{
A list with components 
\item{error}{List of average prediction error across folds. Each list component corresponds to one regression method and contains a matrix of prediction errors for each combination of tuning parameters.}
\item{pars}{Best tuning parameters for each regression method.}
\item{lambda}{\code{lambda} argument}
\item{r}{\code{r} argument}
}

\references{
\url{https://en.wikipedia.org/wiki/Cross-validation_(statistics)#k-fold_cross-validation}
}




\examples{
## Simulate data
sim <- simulate.rr.model(n = 100, p = 50, q = 20, rho = .5, sigma2 = .25)
x <- sim$x
y <- sim$y
b <- sim$b # regression matrix
qr(b)$rank # true rank = 10

## Specify grid of tuning parameters
lambda <- c(0, .01, .1, 1, 10, 100)
r <- c(1, 5, 10, 15, 20) 

## Perform 5-fold cross-validation 
cross <- kcv(x, y, lambda, r, seed = 22)
cross$pars # best parameters
# RRR and RR recover the true regression rank 
# The ranks in PLS and PCR have a different meaning as they pertain to different data models
# They should not be expected to match the true rank in a reduced-rank regression model
}
