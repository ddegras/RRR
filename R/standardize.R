standardize <-
function(x) {
if (!is.matrix(x)) {
	x <- x - mean(x)
	nrm <- sqrt(sum(x^2))
	return(if (nrm == 0) x else x/nrm)
}
x <- x - matrix(colMeans(x), nrow(x), ncol(x), byrow = TRUE)
nrm <- sqrt(colMeans(x^2))
nrm[nrm == 0] <- 1
x / matrix(nrm, nrow(x), ncol(x), byrow = TRUE)
}
