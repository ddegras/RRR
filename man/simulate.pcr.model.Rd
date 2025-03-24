\name{simulate.pcr.model}
\alias{simulate.pcr.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Principal Component Regression Model
}
\description{
Simulate PCR model }
\usage{
simulate.pcr.model(n, p, q, r, sigmax = NULL, sigmay = 0.01, center = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of cases}
  \item{p}{
Number of predictor variables}
  \item{q}{
Number of response variables}
  \item{r}{
Number of principal components generating the response
}
  \item{sigmax}{
Covariance matrix of the predictors
}
  \item{sigmay}{
Standard deviation of error term in response
}
\item{center}{Should the data be centered?}
}


\details{
If \code{sigmax} is not specified, it is generated randomly and scaled as a correlation matrix. 

A n-by-p matrix \eqn{S} of PC scores is created by: 1) randomly generating an n-by-p matrix from the N(0,1) distribution, 2) multiplying each column this matrix by the square root of the corresponding eigenvalue of \eqn{\Sigma_X} (eigenvalues are sorted in descending order). 
After that, the predictor matrix  obtains as \eqn{X=SV'} where \eqn{V} is a matrix of eigenvectors of \eqn{\Sigma_X} (principal directions, \eqn{V'V=I_p}).

An r-by-q regression matrix \eqn{\Gamma} (in PC coordinates) is generated from the N(0,1) distribution. The response matrix is obtained as \eqn{Y=S_r\Gamma + E}, where the matrix \eqn{S_r} contains the first \eqn{r} columns of \eqn{S} and the error matrix \eqn{E} follows the \eqn{N(0,\sigma_Y^2)} distribution. 

The p-by-q regression matrix in the original \eqn{x} coordinates can be recovered as \eqn{B=V\Gamma}. The regression model equivalently expresses as \eqn{Y=XB+E}. 

Note: the data matrices \eqn{X} and \eqn{Y} are typically centered in PCR. However, we do not scale them as this would break the regression model (the parameters \eqn{\Sigma_X}, \eqn{B}, and \eqn{\sigma_Y} would change). In this way, the user can generate data according to the PCR model, calculate the PCR estimate \eqn{\hat{B}}, and compare it to the true \eqn{B}. 
}

\value{
A list with components
\item{x}{Predictor matrix (n-by-p)}
\item{y}{Response matrix (n-by-q)}
\item{b}{regression matrix (p-by-q)}

}



\seealso{
\code{\link{simulate.pls.model}}, \code{\link{simulate.rr.model}}
}

\references{ 
\url{https://en.wikipedia.org/wiki/Principal_component_regression}
}

\examples{
## Generate data
sim <- simulate.pcr.model(n=10, p=5, q=2, r=2)
## Run PCR
pcrfit <- pcr(sim$x, sim$y, 2)
## Compare the estimated regression coefficients to their true values
pcrfit$coef
sim$b
}


