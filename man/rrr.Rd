\name{rrr}
\alias{rrr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reduced Rank Ridge Regression}

\description{
Perform reduced rank ridge regression
}
\usage{
rrr(x, y, lambda = NULL, r = NULL, nlambda = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Predictor matrix (cases-by-predictors)
}
  \item{y}{
Response matrix (cases-by-responses)
}
  \item{lambda}{
Smoothing parameters used in ridge regression stage
}
  \item{r}{
Ranks used in reduced rank regression stage
}
  \item{nlambda}{
Number of values in \code{lambda}. Will be ignored if \code{lambda} is specified explicitly. 
}
}

\details{
See Mukherjee and Zhu (2011) for details. 

This function covers reduced rank regression (\code{\link{rr}}) and ridge regression (\code{\link{ridge}}) as special cases (\eqn{\lambda=0} and \eqn{r} equal to smallest dimension of \code{x} and \code{y}, respectively). 

}

\value{
A list with components
\item{coef}{regression coefficients (predictor-by-response-by-lambda value-by-rank)}
\item{fitted}{fitted responses (case-by-response-by-lambda value-by-rank)}
\item{lambda}{\code{lambda} argument}
\item{r}{effective values of \code{r}}
}

\references{
Ashin Mukherjee and Ji Zhu. \emph{Reduced rank ridge regression and its kernel extensions.} Statistical analysis and data mining, 4(6):612–622, 2011.
}


\seealso{
\code{\link{rr}}, \code{\link{ridge}}, \code{\link{simulate.rr.model}}
}

\examples{
## Simulate data
sim <- simulate.rr.model(n = 100, p = 50, q = 20, rho = .5, sigma2 = .25)
x <- sim$x
y <- sim$y
b <- sim$b # regression matrix
qr(b)$rank # true rank = 10

## Fit RRR regression
lambda <- c(0, .01, .1, 1, 10, 100)
r <- c(1, 5, 10, 15, 20) 
rrrfit <- rrr(x, y, lambda, r)

## Evaluate mean squared error
mse <- apply(rrrfit$coef - as.vector(b), 3:4, function(x) mean(x^2))
rownames(mse) <- lambda
colnames(mse) <- r
mse
}

