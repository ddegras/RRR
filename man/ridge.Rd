\name{ridge}
\alias{ridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Ridge Regression
}
\description{
Perform multivariate ridge regression with common tuning parameter for all responses
}
\usage{
ridge(x, y, lambda = NULL, nlambda = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Predictor matrix (cases-by-predictors)}
  \item{y}{
Response matrix (cases-by-responses)}
  \item{lambda}{
Tuning parameters to use in ridge regression}
  \item{nlambda}{
Number of tuning parameters to use if \code{lambda} is not specified
}
}

\details{
If \code{lambda} is not specified, the function will generate a grid of values spanning a wide dynamic range based on the eigenvalues of \eqn{X'X}. 
}

\value{
A list with components
\item{coef}{regression coefficients (predictor-by-response-by-lambda value-by-rank)}
\item{fitted}{fitted responses (case-by-response-by-lambda value-by-rank)}
\item{lambda}{\code{lambda} argument}
}

\references{
\url{https://en.wikipedia.org/wiki/Ridge_regression}
}



\seealso{
\code{\link{rrr}}
}
\examples{
sim <- simulate.rr.model(n = 10, p = 5, q = 5, rho = .5, sigma2 = .1, model = 2)
ridgefit <- ridge(sim$x, sim$y)
mse <- apply(ridgefit$coef - as.vector(sim$b), 3, function(x) mean(x^2))
plot(ridgefit$lambda, mse, xlab = "lambda", log = "xy")
}

