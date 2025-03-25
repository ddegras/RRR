kcv <-
function(x, y, lambda = NULL, r = NULL, nfolds = 5, 
	method = c("rrr", "rr", "ridge", "pls", "pcr"), 
	seed = NULL)
{
stopifnot((!is.null(lambda)) || (!is.null(r)))
stopifnot(is.null(lambda) || all(lambda >= 0))
if (!is.null(r)) r <- as.integer(r) 
stopifnot(is.null(r) || all(r > 0))
if (!is.matrix(x)) x <- as.matrix(x)
if (!is.matrix(y)) y <- as.matrix(y)
stopifnot(nrow(x) == nrow(y))
method <- match.arg(method, several.ok = TRUE)
n <- nrow(y)
if (!is.null(seed)) set.seed(seed)
folds <- sample(rep_len(1:nfolds, n))

if (is.null(r)) {
	method <- setdiff(method, c("rr", "rrr", "pls", "pcr"))
} else {
	rpcr <- r
	rankmax <- min(dim(x))
	rpcr[r > rankmax] <- rankmax
	rpcr <- unique(rpcr)
	rankmax <- min(dim(x), dim(y))
	r[r > rankmax] <- rankmax
	r <- unique(r)
}
if (is.null(lambda)) {
	method <- setdiff(method, c("ridge", "rrr"))
}
if (length(method) == 0)
	stop("Make sure to specify 'lambda' and/or 'r' as required by 'method'.")
nlam <- length(lambda)
nr <- length(r)
nrpcr <- length(rpcr)
err <- list(rrr = array(dim = c(nfolds, nlam, nr)), 
	rr = matrix(, nfolds, nr), ridge = matrix(, nfolds, nlam), 
	pls = matrix(, nfolds, nr), pcr = matrix(, nfolds, nrpcr))
err <- err[method]
	
for (k in 1: nfolds) {
	train <- which(folds != k)
	holdout <- which(folds == k)
	
	xtrain <- x[train,]
	ytrain <- y[train,]
	xout <- x[holdout,]
	yout <- y[holdout,]
	
	## Fit method to training sample
	errfun <- function(b) sum((yout - xout %*% b)^2)
	if ("rrr" %in% method) {
		b <- rrr(xtrain, ytrain, lambda, r)$coef
		err[["rrr"]][k,,] <- if (nr == 1 && nlam == 1) {
			errfun(b) } else apply(b, 3:4, errfun)
	}
	if ("rr" %in% method) {
		b <- rr(xtrain, ytrain, r)$coef
		err[["rr"]][k,] <- if (nr == 1) {
			errfun(b) } else apply(b, 3, errfun)
	}
	if ("ridge" %in% method) {
		b <- ridge(xtrain, ytrain, lambda)$coef
		err[["ridge"]][k,] <- if (nlam == 1) {
			errfun(b) } else apply(b, 3, errfun)
	}
	if ("pls" %in% method) {
		b <- pls(xtrain, ytrain, r)$coef
		err[["pls"]][k,] <- if (nr == 1) {
			errfun(b) } else apply(b, 3, errfun)
	}
	if ("pcr" %in% method) {
		b <- pcr(xtrain, ytrain, rpcr)$coef
		err[["pcr"]][k,] <- if (nrpcr == 1) {
			errfun(b) } else apply(b, 3, errfun)
	}		
}

## Average prediction error across folds
ape <- lapply(err, function(x) colSums(x) / length(y))

## Best parameter(s)
pars <- list()
idx <- sapply(ape, which.min)
if ("rrr" %in% method) {
	idx1 <- arrayInd(idx["rrr"], c(nlam, nr))
	pars$rrr <- c(lambda = lambda[idx1[1]], r = r[idx1[2]])
}
if ("rr" %in% method) pars$rr <- c(r = r[idx["rr"]])
if ("ridge" %in% method) pars$ridge <- c(lambda = lambda[idx["ridge"]])
if ("pls" %in% method) pars$pls <- c(r = r[idx["pls"]])
if ("pcr" %in% method) pars$pcr <- c(r = rpcr[idx["pcr"]])

list(error = ape, pars = pars, lambda = lambda, r = r)
}
