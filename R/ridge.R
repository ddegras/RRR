ridge <- function(x, y, lambda = NULL, nlambda = 20)
{
if (!is.matrix(x)) x <- as.matrix(x)
if (!is.matrix(y)) y <- as.matrix(y)
stopifnot(nrow(x) == nrow(y))
stopifnot(is.null(lambda) || all(lambda >= 0))
if (!is.null(lambda)) nlambda <- length(lambda)	
n <- nrow(x)
p <- ncol(x)
q <- ncol(y)

## Trivial case
if (all(x == 0)) {
	return(list(coeffs = matrix(0, p, q), fitted = matrix(0, n, q), r = 1))
}

coeffs <- array(dim = c(p, q, nlambda))
fitted <- array(dim = c(n, q, nlambda))

svdx <- svd(x)
d <- svdx$d
d2 <- d^2
uty <- crossprod(svdx$u, y)
if (is.null(lambda)) {
	lambda <- 10^seq(from = log10(d2[1]), 
		to = log10(d2[1]) - 8, len = nlambda)
}

for (l in 1:nlambda) {
	scl <- d / (d2 + lambda[l])
	scl[is.nan(scl)] <- 0
	coeffs[,,l] <- svdx$v %*% (uty * scl)
	fitted[,,l] <- svdx$u %*% (uty * (d*scl))
}

if (nlambda == 1) {
	dim(coeffs) <- c(p, q)
	dim(fitted) <- c(n, q)
} else {
	dimnames(coeffs) <- list(pred = NULL, resp = NULL, lambda = lambda)
	dimnames(fitted) <- list(case = NULL, resp = NULL, lambda = lambda)	
}


list(coef = coeffs, fitted = fitted, lambda = lambda)	
}
