ridge <- function(x, y, lambda = NULL, nlambda = 20, retsvd = TRUE)
{
if (!is.matrix(x)) x <- as.matrix(x)
if (!is.matrix(y)) y <- as.matrix(y)
stopifnot(nrow(x) == nrow(y))
stopifnot(is.null(lambda) || all(lambda >= 0))
if (!is.null(lambda)) nlambda <- length(lambda)	
n <- nrow(x)
p <- ncol(x)
q <- ncol(y)

coeffs <- array(dim = c(p, q, nlambda))
fitted <- array(dim = c(n, q, nlambda))
dimnames(coeffs) <- list(pred = NULL, resp = NULL, lambda = NULL)
dimnames(fitted) <- list(case = NULL, resp = NULL, lambda = NULL)	

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

if (any(c(p, q, nlambda) == 1))
	coeffs <- drop(coeffs)
if (any(c(n, q, nlambda) == 1))
	fitted <- drop(fitted)

list(coef = coeffs, fitted = fitted, lambda = lambda, 
	svdx = if (retsvd) svdx else NULL)	
}
