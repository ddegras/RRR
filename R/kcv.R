kcv <- function(x, y, lambda = NULL, r = NULL, nfolds = 5, 
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
	
	xtrain <- x[train,,drop=FALSE]
	ytrain <- y[train,,drop=FALSE]
	xout <- x[holdout,,drop=FALSE]
	yout <- y[holdout,,drop=FALSE]
	
	## Fit method to training sample
	errfun <- function(b) sum((yout - xout %*% b)^2)
	if ("rrr" %in% method) {
		b <- rrr(xtrain, ytrain, lambda, r)$coef
		# ndims <- length(dim(b))
		# err[["rrr"]][k,,] <- if (ndims == 2) {
			# errfun(b) } else apply(b, 3:ndims, errfun)
		errk <- apply(b, 3:4, errfun)
		if (nrow(errk) < nr) {
			nrk <- nrow(errk)
			errk <- errk[c(1:nrk, rep(nrk, nr - nrk)),, drop = FALSE]
		}
		err[["rrr"]][k,,] <- errk
	}
	if ("rr" %in% method) {
		b <- rr(xtrain, ytrain, r)$coef
		errk <- apply(b, 3, errfun)
		if (length(errk) < nr) {
			nrk <- length(errk)
			errk <- c(errk, rep(errk[nrk], nr - nrk))
		}
		err[["rr"]][k,] <- errk
		# err[["rr"]][k,] <- if (nr == 1) {
			# errfun(b) } else apply(b, 3, errfun)
	}
	if ("ridge" %in% method) {
		b <- ridge(xtrain, ytrain, lambda)$coef
		# err[["ridge"]][k,] <- if (nlam == 1) {
			# errfun(b) } else apply(b, 3, errfun)
		errk <- apply(b, 3, errfun)
		if (length(errk) < nr) {
			nrk <- length(errk)
			errk <- c(errk, rep(errk[nrk], nr - nrk))
		}
		err[["ridge"]][k,] <- errk
	}
	if ("pls" %in% method) {
		b <- pls(xtrain, ytrain, r)$coef
		errk <- apply(b, 3, errfun)
		if (length(errk) < nr) {
			nrk <- length(errk)
			errk <- c(errk, rep(errk[nrk], nr - nrk))
		}
		err[["pls"]][k,] <- errk
		# err[["pls"]][k,] <- if (nr == 1) {
			# errfun(b) } else apply(b, 3, errfun)
	}
	if ("pcr" %in% method) {
		b <- pcr(xtrain, ytrain, rpcr)$coef
		errk <- apply(b, 3, errfun)
		if (length(errk) < nrpcr) {
			nrk <- length(errk)
			errk <- c(errk, rep(errk[nrk], nrpcr - nrk))
		}
		err[["pcr"]][k,] <- errk
		# err[["pcr"]][k,] <- if (nrpcr == 1) {
			# errfun(b) } else apply(b, 3, errfun)
	}		
}

## Average prediction error across folds
ape <- lapply(err, function(x) colSums(x) / length(y))

## Best parameter(s)
pars <- vector("list", length(method))
names(pars) <- method
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
