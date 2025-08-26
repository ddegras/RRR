cv.ridge <- function(object, y, x, lambda) {
if (!is.matrix(y)) 
	y <- as.matrix(y)
if (missing(object)) 
	object <- ridge(x, y, lambda, retsvd = TRUE)
if (missing(lambda)) 
	lambda <- object$lambda
nlambda <- length(lambda)
n <- nrow(y)
q <- ncol(y)
cv <- gcv <- matrix(, q, nlambda)
dims <- dim(object$fitted)
if (length(dims) < 3) 
	dim(object$fitted) <- c(n, q, nlambda)
svdx <- object$svdx
d2 <- svdx$d^2
dut <- svdx$d * t(svdx$u)
for (j in 1:nlambda) {
	e <- y - object$fitted[,,j]
	h <- colSums(dut * (dut / (d2 + lambda[j])))
	cv[,j] <- colMeans((e/(1-h))^2)
	gcv[,j] <- colMeans(e^2) / (1 - mean(h))^2
}
cv.best <- apply(cv, 1, which.min) 
gcv.best <- apply(gcv, 1, which.min)
list(lambda.best = cbind(cv = lambda[cv.best], gcv = lambda[gcv.best]),
	lambda.best.idx = cbind(cv = cv.best, gcv = gcv.best),
	lambda = lambda, cv = cv, gcv = gcv)
}
