rrr <- function(x, y, lambda = NULL, r = NULL, nlambda = 20)
{
if (!is.matrix(x)) x <- as.matrix(x)
if (!is.matrix(y)) y <- as.matrix(y)
stopifnot(nrow(x) == nrow(y))
stopifnot(is.null(lambda) || all(lambda >= 0))
if (!is.null(lambda)) nlambda <- length(lambda)	
n <- nrow(x)
p <- ncol(x)
q <- ncol(y)
rankmax <- min(n,p,q)
if (is.null(r)) {
	r <- 1:rankmax
} else {
	r <- as.integer(r)
	stopifnot(all(r > 0))
	r[r > rankmax] <- rankmax
	r <- unique(r)
}	
nr <- length(r)
rmax <- max(r)

## Trivial case
if (all(x == 0)) {
	return(list(coeffs = matrix(0, p, q), fitted = matrix(0, n, q), r = 1))
}

coeffs <- array(dim = c(p, q, nlambda, nr))
fitted <- array(dim = c(n, q, nlambda, nr))

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
	bridge <- svdx$v %*% (scl * uty)
	yridge <- svdx$u %*% ((d * scl) * uty)
	yrext <- if (lambda[l] == 0) yridge else rbind(yridge, bridge)
	v <- tryCatch(
				suppressWarnings(svds(yrext, rmax, nu = 0)$v),
				error = function(e) svd(yrext, 0, rmax)$v)
	for (j in 1:nr) {
		coeffs[,,l,j] <- tcrossprod(bridge %*% v[,1:r[j]], v[,1:r[j]])
		fitted[,,l,j] <- tcrossprod(yridge %*% v[,1:r[j]], v[,1:r[j]])
	}	
}

## Reduce output dimension if single r value
keepdim <- c(TRUE, TRUE, nlambda > 1, nr > 1)
if (any(!keepdim)) {
	dim(coeffs) <- dim(coeffs)[keepdim]
	dim(fitted) <- dim(fitted)[keepdim]
}
dimnames(coeffs) <- list(pred = NULL, resp = NULL, 
	lambda = lambda, r = r)[keepdim]
dimnames(fitted) <- list(case = NULL, resp = NULL, 
	lambda = lambda, r = r)[keepdim]


list(coef = coeffs, fitted = fitted, lambda = lambda, r = r)	
}
