\name{rr}
\alias{rr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reduced Rank Regression}

\description{
Perform reduced rank regression}

\usage{
rr(x, y, r = NULL)
}

\arguments{
  \item{x}{
Predictor matrix (cases-by-predictors)
}
  \item{y}{
Response matrix (cases-by-responses)
}
  \item{r}{
Ranks used in reduced rank regression stage
}
}



\value{
A list with components
\item{coef}{regression coefficients (predictor-by-response-by-lambda value-by-rank)}
\item{fitted}{fitted responses (case-by-response-by-lambda value-by-rank)}
\item{r}{effective values of \code{r}}
}

\references{
Gregory C. Reinsel, Raja P. Velu, and Kun Chen. \emph{Multivariate Reduced-Rank Regression. Theory, Methods and Applications.} Springer, 2022. 
}


\seealso{
\code{\link{rrr}}, \code{\link{simulate.rr.model}} 
}

\examples{
## Simulate data
sim <- simulate.rr.model(n = 100, p = 50, q = 20, rho = .5, sigma2 = .25)
x <- sim$x
y <- sim$y
b <- sim$b # regression matrix
qr(b)$rank # true rank = 10

## Fit RR regression
r <- c(1, 5, 10, 15, 20) 
rrfit <- rrr(x, y, r)

## Evaluate mean squared error
mse <- apply(rrfit$coef - as.vector(b), 3, function(x) mean(x^2))
plot(r, mse)
}


