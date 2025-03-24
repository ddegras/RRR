\name{simulate.pls.model}
\alias{simulate.pls.model}

\title{Partial Least Squares Regression Model}

\description{
Simulate a Partial Least Squares Regression (PLSR) Model
}

\usage{
simulate.pls.model(n, p, q, r, sigma2 = NULL, center = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of cases
}
  \item{p}{
Number of predictors
}
  \item{q}{
Number of responses
}
  \item{r}{
Number of factors. Should be less than or equal to \code{p} and \code{q}. 
}
  \item{sigma2}{Scales of model parameters. Should have length 1 or 4. Default = 1. 
}
\item{center}{Should the data be centered?}

}
\details{
The PLSR model has the form \deqn{X = TP'+E, \quad Y = TQ' + F}
where \eqn{T} is an n-by-r matrix of common factors to \eqn{X} and \eqn{Y}, 
\eqn{P} is an p-by-r matrix of x-loadings, \eqn{Q} is an q-by-r matrix of y-loadings, 
\eqn{E} and \eqn{F} are measurement error matrices for \eqn{X} and \eqn{Y}, respectively. 

The factor matrix \eqn{T} is generated from the N(0,1). The matrices \eqn{P}, \eqn{Q}, \eqn{E}, and \eqn{F} are independently generated from normal distributions with mean zero and variances given in \code{sigma2}. 

}

\value{
A list with components
\item{x}{Predictor matrix (n-by-p)}
\item{y}{Response matrix (n-by-q)}
\item{b}{Regression matrix (p-by-q)}
}


\references{
\url{https://en.wikipedia.org/wiki/Partial_least_squares_regression}
}


\seealso{\code{\link{pls}}
}

\examples{
simulate.pls.model(10, 5, 5, 2)
}
