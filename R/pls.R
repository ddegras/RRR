pls <- function(x, y, r = NULL)
{
if (!is.matrix(x)) x <- as.matrix(x)
if (!is.matrix(y)) y <- as.matrix(y)
stopifnot(nrow(x) == nrow(y))
n <- nrow(x)
p <- ncol(x)
q <- ncol(y)
rankmax <- min(dim(x))
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

pp <- w <- matrix(0, ncol(x), rmax)
cc <- matrix(0, ncol(y), rmax)
b <- numeric(rmax)
coeffs <- array(dim = c(p, q, nr))
fitted <- array(dim = c(n, q, nr))
y0 <- y
count <- 0
for (i in 1:rmax) {
	svdxty <- tryCatch(svds(crossprod(x, y), 1), 
		error = function(e) svd(crossprod(x, y), 1, 1))
	if (svdxty$d[1] <= 1e-14) {
		fitted[,,r >= i] <- as.vector(y0 - y)
		break
	}
	w[,i] <- svdxty$u # x weights 
	tt <- x %*% w[,i] # x scores
	tt <- tt / sqrt(sum(tt^2)) # rescale
	pp[,i] <- crossprod(x, tt) # regression coeff for raw x vs x scores
	xhat <- tcrossprod(tt, pp[,i]) # fitted x
	cc[,i] <- svdxty$v # y weights
	u <- y %*% cc[,i] # y scores
	b[i] <- sum(tt * u) # slope in regression of x score vs y score 
	yhat <- tcrossprod(b[i] * tt, cc[,i]) # fitted y
	x <- x - xhat # deflate x
	y <- y - yhat # deflate y
	if (i %in% r) {
		count <- count + 1
		fitted[,,count] <- y0 - y
	}
}	

## Calculate PLS regression coeffs and effective rank
for (i in 1:nr) {
	ri <- r[i] 
	ptw <- crossprod(pp[,1:ri], w[,1:ri])
	btc <- b[1:ri] * t(cc[,1:ri])
	bb <- tryCatch(w[,1:ri] %*% solve(ptw, btc), 
		error = function(e) NULL)
	if (is.null(bb)) {
		svdptw <- svd(ptw)
		nz <- (svdptw$d >= min(1e-14, svdptw$d[1] * 1e-8))
		svdptw$d[nz] <- 1 / svdptw$d[nz]
		svdptw$d[!nz] <- 0
		invptw <- svdptw$v %*% (svdptw$d * t(svdptw$u))
		bb <- w[,1:ri] %*% invptw %*% btc
	}
	coeffs[,,i] <- bb
}

if (nr == 1) {
	dim(coeffs) <- c(p, q)
	dim(fitted) <- c(n, q) 
} else {
	dimnames(coeffs) <- list(pred = NULL, resp = NULL, r = r)
	dimnames(fitted) <- list(case = NULL, resp = NULL, r = r)	
	
}

list(coef = coeffs, fitted = fitted, r = r)	
}
