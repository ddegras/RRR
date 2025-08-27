simulate.rr.model <- function(n, p, q, r, rho, sigma, center = TRUE) 
{
    stopifnot(abs(rho) < 1)
    n <- as.integer(n)
    p <- as.integer(p)
    q <- as.integer(q)
    r <- as.integer(r)
    stopifnot(min(n, p, q, r) > 0)
    stopifnot(r <= min(p, q))
    sigma <- rep_len(sigma, p)
    stopifnot(all(sigma >= 0))
    
    alpha <- sqrt(1 - rho^2)
    x <- matrix(rnorm(n * p, sd = alpha), p, n)
    x[1,] <- x[1,] / alpha
    x <- filter(x, filter = rho, method = "recursive")
    x <- t(x)
    b <- matrix(rnorm(p * r), p, r)
	score <- matrix(rnorm(q * r), r, q)
	e <- matrix(rnorm(n * q, sd = sigma), n, q, byrow = TRUE)
    y <- (x %*% b) %*% score + e
    b <- b %*% score
    if (center) {
        x <- sweep(x, 2, colMeans(x), "-")
        y <- sweep(y, 2, colMeans(y), "-")
    }
    list(x = x, y = y, b = b)
}


# simulate.rr.model <- function(n, p, q, rho, sigma2, model = 1, center = TRUE)
# {
# stopifnot(model %in% (1:3))
# stopifnot(abs(rho) < 1)

# ## Regression coefficients
# b <- matrix(rnorm(p*q), p, q)
# svdb <- svd(b)
# if (model == 1) {
	# half <- floor(min(p,q)/2)
	# svdb$d[1:half] <- 2
	# svdb$d[-(1:half)] <- 0
# } else if (model == 2) {
	# svdb$d <- rep(1, min(p,q))
# } else { # model 3
	# svdb$d <- rep(0, min(p,q))
	# svdb$d[1] <- 5
# }
# b <- svdb$u %*% (svdb$d * t(svdb$v))

# ## Predictors
# x <- matrix(
	# data = c(rnorm(n), rnorm(n*(p-1), sd = sqrt(1 - rho^2))),
	# p, n, byrow = TRUE)
# x <- filter(x, filter = rho, method = "recursive")
# x <- t(x)

# ## Responses
# y <- x %*% b + rnorm(n*q, sd = sqrt(sigma2))

# if (center) {
	# x <- sweep(x, 2, colMeans(x), "-")
	# y <- sweep(y, 2, colMeans(y), "-")
# }

# list(x = x, y = y, b = b)
# }
