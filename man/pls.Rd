\name{pls}
\alias{pls}

\title{Partial Least Squares Regression}

\description{
Performs partial least squares regression (PLSR)
}

\usage{
pls(x, y, r = NULL)
}

\arguments{
  \item{x}{
Predictor matrix (cases-by-predictors)
}
  \item{y}{
Response matrix (cases-by-responses)
}
  \item{r}{
Ranks to use in PLSR
}
}


\value{
A list with components
\item{coef}{Regression coefficients (predictors-by-responses-by-rank)}
\item{fitted}{Regression coefficients (cases-by-responses-by-rank)}
\item{r}{Effective ranks used}
}

\references{
\url{https://en.wikipedia.org/wiki/Partial_least_squares_regression} 
}


\seealso{
\code{\link{simulate.pls.model}} 
}

\examples{
n <- 10
p <- q <- 5
r <- 2
sigma2 <- c(1, 1, .05, .05)
sim <- simulate.pls.model(n, p, q, r, sigma2)
plsfit <- pls(sim$x, sim$y, r)
mean((plsfit$coef - sim$b)^2) # mean squared error
}

