\name{cv.ridge}
\alias{cv.ridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-Validation for Univariate Ridge Estimator
}
\description{
Selects the shrinkage parameter of the ridge estimator by cross-validation and generalized cross-validation. 
}
\usage{
cv.ridge(object, y, x, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object created by a call to \code{\link{ridge}}. 
}
  \item{y}{
The vector or matrix of responses used to create \code{object}
}
  \item{x}{
The matrix of predictors used to create \code{object}  
}
  \item{lambda}{
The vector of shrinkage parameters used to create \code{object}
}
}

\details{
The function \code{cv.ridge} can work in one of two ways: 1) call the function \code{\link{ridge}} first (with argument \code{retsvd = TRUE}) and then call \code{cv.ridge} on the result or 2) directly do ridge estimation and CV/GCV in one step. In the first case, only the arguments \code{object} and \code{y} need to be specified in \code{cv.ridge}; the other arguments can be omitted. In the second case, only the arguments \code{y}, \code{x}, and \code{lambda} need to be specified. 

The argument \code{y} can be a vector or a matrix, in which case CV/GCV will be calculated separately for each column. 
}
\value{
A list with components 
\item{lambda.best}{Matrix of best \code{lambda} value for each response variable (rows) and criterion (CV/GCV in columns).}
\item{lambda.best.idx}{Matrix of indices associated to \code{lambda.best} and \code{lambda}.}
\item{lambda}{Vector of shrinkage parameters used in ridge regression.}
\item{cv}{Matrix of CV scores for each response variable (rows) and criterion (CV/GCV, columns).}
\item{gcv}{Matrix of GCV scores for each response variable (rows) and criterion (CV/GCV, columns).}
}
\references{
Golub, G., Wahba, G. and Heat, C. (1979). Generalized Cross Validation as a Method for Choosing a Good Ridge Parameter. \emph{Technometrics}. 21, 215–223. \doi{10.2307/1268518}.
}


\seealso{
\code{\link{ridge}}
}
\examples{
dat <- simulate.rr.model(100, 50, 20, 10, .9, .5)
lambda <- 10^seq(-3, -1, by = .05)
fit <- ridge(dat$x, dat$y, lambda)
lambda.hat <- cv.ridge(fit, dat$y)
# Equivalent to the two lines above:
# lambda.hat <- cv.ridge(x = dat$x, y = dat$y, lambda = lambda) 
matplot(lambda, cbind(lambda.hat$cv[1,], lambda.hat$gcv[1,]), type = "l",
	xlab = "lambda", ylab = "CV", log = "x")
title("CV and GCV Plots for Y1")
legend("top", c("CV", "GCV"), lty = 1:2, col = 1:2, cex = .9)
abline(v = lambda.hat$lambda.best[1,], lty = 2, col = "grey", lwd = .75)
}
