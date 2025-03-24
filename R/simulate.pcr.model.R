simulate.pcr.model <- function(n, p, q, r, sigmax = NULL, sigmay = .01, center = TRUE)
{
r <- min(p, q, r)
if (is.null(sigmax)) {
	sigmax <- crossprod(matrix(rnorm(p^2), p, p)) 
	sigmax <- cov2cor(sigmax)
} else {
	stopifnot((nrow(sigmax) == p) && ncol(sigmax) == p)
}
eig <- eigen(sigmax, TRUE)
loadings <- eig$vectors
scores <- matrix(rnorm(n*p), n, p) %*% diag(eig$values)
x <- tcrossprod(scores, loadings)
gamma <- matrix(rnorm(q*r), r, q) # regression coefficients in PC basis
b <- loadings[,1:r] %*% gamma # regression coefficients in original x coordinates
e <- matrix(rnorm(n*q, sd = sigmay), n, q)
y <- scores[,1:r] %*% gamma + e # same as y = x %*% b + e
if (center) {
	x <- x - rep(colMeans(x), each = n)
	y <- y - rep(colMeans(y), each = n)
}
list(x = x, y = y, b = b)		
}
