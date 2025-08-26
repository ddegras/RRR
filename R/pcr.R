pcr <-
function(x, y, r = NULL)
{
if (!is.matrix(x)) x <- as.matrix(x)
if (!is.matrix(y)) y <- as.matrix(y)
stopifnot(nrow(x) == nrow(y))
stopifnot(is.null(r) || all(r >= 1))
n <- nrow(x)
p <- ncol(x)
q <- ncol(y)
rankmax <- min(n, p)
if (is.null(r)) {
	r <- 1:rankmax
} else {
	r <- as.integer(r)
	r[r > rankmax] <- rankmax 
	r <- unique(r)
}
nr <- length(r)
rmax <- max(r)

## Trivial case
if (all(x == 0)) {
	return(list(coeffs = matrix(0, p, q), fitted = matrix(0, n, q), r = 1))
}

## Carry out SVD	
svdx <- if (rmax <= rankmax / 2) {
	tryCatch(RSpectra::svds(x, rmax), error = function(e) NULL) 
	} else NULL
if (is.null(svdx)) svdx <- svd(x, rmax, rmax)
rankx <- sum(svdx$d >= svdx$d[1] * 1e-8) # ensure good conditioning of x
	
## PC Regression
coeffs <- array(dim = c(p, q, nr), 
	dimnames = list(pred = NULL, resp = NULL, r = r))
fitted <- array(dim = c(n, q, nr),
	dimnames = list(case = NULL, resp = NULL, r = r))
for (i in 1:nr) {
	ri <- min(r[i], rankx)
	uty <- crossprod(svdx$u[,1:ri], y)
	coeffs[,,i] <- svdx$v[,1:ri] %*% (uty / svdx$d[1:ri])
	fitted[,,i] <- svdx$u[,1:ri] %*% uty
}
# if (nr == 1) { 
	# dim(coeffs) <- c(p, q)
	# dim(fitted) <- c(n, q)
# } else {
	# dimnames(coeffs) <- list(pred = NULL, resp = NULL, r = r)
	# dimnames(fitted) <- list(case = NULL, resp = NULL, r = r)	
# }

list(coef = coeffs, fitted = fitted, r = r)	
}
