\name{simulate.rr.model}
\alias{simulate.rr.model}

\title{
Multivariate Regression Models of Mukherjee and Zhu (2011)}

\description{
Simulate the models used in Mukherjee and Zhu (2011)
}

\usage{
simulate.rr.model(n, p, q, r, rho, sigma, center = TRUE)
}


\arguments{
  \item{n}{
number of cases}
  \item{p}{
number of predictor variables
}
  \item{q}{
number of response variables 
}
\item{r}{
	rank of regression matrix
}
  \item{rho}{
correlation among predictor variables}
  \item{sigma}{
error standard deviation
}
\item{center}{Should the predictors and responses be centered? Default = TRUE}
}

\details{
Simulates the regression model \deqn{Y = XB + E}
where \eqn{Y} is a response matrix of size n-by-q, 
\eqn{X} is a response matrix of size n-by-p, 
\eqn{B} is a regression matrix of size p-by-q,
and \eqn{E} is an error matrix. 
Refer to Mukherjee and Zhu (2011) for details.
}

\value{
A list with components
\item{x}{Predictor matrix (n-by-p)}
\item{y}{Response matrix (n-by-q)}
\item{b}{Regression coefficient matrix (p-by-q)}
}

\references{
Ashin Mukherjee and Ji Zhu. \emph{Reduced rank ridge regression and its kernel extensions.} Statistical analysis and data mining, 4(6):612–622, 2011.}


\seealso{
\code{\link{simulate.pls.model}}
}

\examples{
sim <- simulate.rr.model(n = 100, p = 50, q = 20, r = 10, rho = .5, sigma = .5)
x <- sim$x
y <- sim$y
b <- sim$b
}
