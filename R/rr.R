rr <- function(x, y, r = NULL)
{
if (!is.matrix(x)) x <- as.matrix(x)
if (!is.matrix(y)) y <- as.matrix(y)
stopifnot(nrow(x) == nrow(y))
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

## OLS fit
olsfit <- suppressWarnings(lsfit(x, y, intercept = FALSE))
bols <- olsfit$coef
yols <- y - olsfit$residuals

# svdx <- if (maxsvd(x)
# keep <- (svdx$d >= max(svdx$d[1]*1e-8, 1e-14))
# if (!all(keep)) {
	# svdx$u <- svdx$u[,keep]
	# svdx$d <- svdx$d[keep]
	# svdx$v <- svdx$v[,keep]
# }
# uty <- crossprod(svdx$u,y)
# bols <- svdx$v %*% (uty / svdx$d) # regression coeffs
# yols <- svdx$u %*% uty # fitted values
if (all(r == rankmax)) 
	return(list(coeffs = bols, fitted = yols, r = r))

## Projections on lower-rank spaces
coeffs <- array(dim = c(p, q, nr))
fitted <- array(dim = c(n, q, nr))
svdyols <- tryCatch(supressWarnings(svds(yols, k = rmax, nu = 0)), 
	error = function(e) NULL)
if (is.null(svdyols)) svdyols <- svd(yols, nu = 0, nv = rmax)
for (i in 1:nr) {
	if (r[i] == rankmax) {
		coeffs[,,i] <- bols
		fitted[,,i] <- yols
	} else {
		v <- svdyols$v[,1:r[i]]
		coeffs[,,i] <- tcrossprod(bols %*% v, v)
		fitted[,,i] <- tcrossprod(yols %*% v, v)
	}	
}

## Reduce output dimension if single r value
if (nr == 1) {
	dim(coeffs) <- c(p, q)
	dim(fitted) <- c(n, q)
} else {
	dimnames(coeffs) <- list(pred = NULL, resp = NULL, r = r)
	dimnames(fitted) <- list(case = NULL, resp = NULL, r = r)	
}

list(coef = coeffs, fitted = fitted, r = r)
}
