\name{pcr}
\alias{pcr}

\title{
Principal Component Regression}

\description{
Regression based on the largest principal components of the predictor matrix
}

\usage{
pcr(x, y, r = NULL)
}

\arguments{
  \item{x}{
Predictor matrix (cases-by-predictors)
}
  \item{y}{
Response matrix (cases-by-responses)
}
  \item{r}{
Numbers of principal components used in regression (vector)
}
}

\value{
A list with components 
\item{coef}{Regression coefficients (predictors-by-responses-by-numbers of PCs)}
\item{fitted}{Regression coefficients (cases-by-responses-by-numbers of PCs)}
\item{r}{Effective numbers of principal components used}
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
pcrfit$coeffs
sim$b
}
