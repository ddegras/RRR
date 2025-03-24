simulate.pls.model <- function(n, p, q, r, sigma2 = NULL, center = TRUE)
{
r <- min(p, q, r)
if (is.null(sigma2)) {
	sigma <- rep(1, 4)
} else {
	stopifnot(all(sigma2 > 0))
	sigma <- rep_len(sqrt(sigma2), 4)
}
tt <- matrix(rnorm(n*r), n, r) # common factors
pp <- matrix(rnorm(p*r, sd = sigma[1]), p, r) # x factor loadings 
qq <- matrix(rnorm(q*r, sd = sigma[2]), q, r) # y factor loadings
ee <- matrix(rnorm(n*p, sd = sigma[3]), n, p) # noise for x
ff <- matrix(rnorm(n*q, sd = sigma[4]), n, q) # noise for y
x <- tcrossprod(tt, pp) + ee
y <- tcrossprod(tt, qq) + ff
b <- pp %*% solve(crossprod(pp), t(qq))
if (center) {
	x <- sweep(x, 2, colMeans(x), "-")
	y <- sweep(y, 2, colMeans(y), "-")
}

list(x = x, y = y, b = b)
}
