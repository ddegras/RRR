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

## OLS fit: minimum-norm solution
olsfit <- NULL
if (n >= p) {
	olsfit <- suppressWarnings(lsfit(x, y, intercept = FALSE))
	if (olsfit$qr$rank == p) {
		bols <- olsfit$coef
		pivot <- olsfit$qr$pivot
		if (!identical(pivot, 1:p)) bols[pivot,] <- bols
		yols <- y - olsfit$residuals
	} else {
		olsfit <- NULL
	}
}
if (is.null(olsfit)) {
	svdx <- svd(x)
	pos <- (svdx$d > max(svdx$d[1] * 1e-8, 1e-14))
	if (any(pos)) {
		uty <- crossprod(svdx$u[,pos], y)
		bols <- svdx$v[,pos] %*% (uty / svdx$d[pos])
		yols <- svdx$u[,pos] %*% uty
	} else {
		bols <- rep(0, p)
		yols <- matrix(0, n, q)
	}
}
if ((nr == 1) && (r == rankmax)) 
	return(list(coef = bols, fitted = yols, r = r))

## Projections on lower-rank spaces
coeffs <- array(, c(p, q, nr), list(pred = NULL, resp = NULL, r = r))
fitted <- array(, c(n, q, nr), list(case = NULL, resp = NULL, r = r))
	
svdyols <- tryCatch(suppressWarnings(svds(yols, k = rmax, nu = 0)), 
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
# if (nr == 1) {
	# dim(coeffs) <- c(p, q)
	# dim(fitted) <- c(n, q)
# } else {
	# dimnames(coeffs) <- list(pred = NULL, resp = NULL, r = r)
	# dimnames(fitted) <- list(case = NULL, resp = NULL, r = r)	
# }

list(coef = coeffs, fitted = fitted, r = r)
}
