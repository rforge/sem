

robust_vcov <- function (sem.obj, adj.obj, data.obj, useFit=FALSE, useGinv=FALSE){
	if (missing(adj.obj) && missing(data.obj)) 
		stop("Need a data or sbchisq object")
	if (missing(adj.obj)) {
		adj.obj <- sbchisq(sem.obj, data.obj, useFit = useFit)
	}
	ncases <- sem.obj$N
	hes <- sem_hessian(adj.obj$w_mat, adj.obj$p_deriv_mat)
	info_m <- try(solve(hes), silent = TRUE)
	if (class(info_m) == "try-error" && useGinv == TRUE) {
		info_m <- ginv(hes)
		ginvFlag <- TRUE
	}
	acov <- info_m %*% (t(adj.obj$p_deriv_mat) %*% adj.obj$w_mat %*% 
			adj.obj$w_adf %*% adj.obj$w_mat %*% adj.obj$p_deriv_mat) %*% 
		info_m
	acov <- acov/(ncases - 1)
	rownames(acov) <- colnames(acov) <- colnames(adj.obj$p_deriv_mat)
	acov
}


sbchisq <- function (sem.obj, sem.data, adj = 1e-04, useFit = FALSE, useGinv = FALSE) 
{
	props <- sem.props(sem.obj)
	sem.prop <- props$chisq
	chisq <- props$chisq
	df <- props$df
	w_adf <- adf.wmat(sem.data)
	w_mat <- ml.wmat(sem.obj, useFit = useFit)
	p_deriv_mat <- delta_matrix(sem.obj)
	ginvFlag <- FALSE
	invMat <- try(solve(t(p_deriv_mat) %*% w_mat %*% p_deriv_mat), 
		silent = TRUE)
	if (class(invMat) == "try-error" && useGinv == TRUE) {
		invMat <- ginv(t(p_deriv_mat) %*% w_mat %*% p_deriv_mat)
		ginvFlag <- TRUE
	}
	res_u <- w_mat - (w_mat %*% p_deriv_mat %*% invMat %*% t(p_deriv_mat) %*% 
			w_mat)
	ug <- res_u %*% w_adf
	scale_stat <- sum(diag(ug))/df
	chisq.scaled <- chisq/scale_stat
	p.old <- pchisq(chisq, df, lower.tail=FALSE)
	p <- pchisq(chisq.scaled, df, lower.tail=FALSE)
	ret <- list(chisq = chisq, t = sem.obj$t, df = df, p.old = p.old, 
		c = scale_stat, chisq.scaled = chisq.scaled, p = p, w_mat = w_mat, 
		p_deriv_mat = p_deriv_mat, w_adf = w_adf, res_u = res_u, 
		N = length(sem.data[, 1]), ginvFlag = ginvFlag)
	class(ret) <- "adjchisq"
	return(ret)
}

adf.wmat <- function (raw_data) 
{
#	n <- length(raw_data[, 1])
#	n.col <- length(names(raw_data))
	names <- colnames(raw_data)
	n <- nrow(raw_data)
	n.col <- ncol(raw_data)
	nc.star <- n.col * (n.col + 1)/2
	nc2 <- n.col^2
#	i1 <- rep(1, n)
#	xbar <- apply(raw_data, 2, sum)/n
#	z <- raw_data - kronecker(xbar, i1)
	z <- scale(raw_data, center=TRUE, scale=FALSE)
	sc <- vector(nc.star, mode="list")
	outnames <- vector(nc.star, mode="character")
	ind <- combn(n.col + 1, 2)
	ind[2, ] <- ind[2, ] - 1
	for (q in 1:nc.star) {
		i <- ind[1, q]
		j <- ind[2, q]
#		a.name <- paste(names(raw_data)[i], names(raw_data)[j], 
#			sep = "_")
		outnames[q] <- paste(names[i], names[j], sep="_")
		sc[[q]] <- z[, i] * z[, j]
	}
	names(sc) <- outnames
	adf_mat <- var(data.frame(sc)) * (n - 1)/n
	return(adf_mat)
}

ml.wmat <- function (sem.obj, useFit = FALSE) 
{
	p <- nrow(sem.obj$C) # length(rownames(sem.obj$C))
	if (useFit) {
		An <- sem.obj$C
	}
	else {
		An <- sem.obj$S
	}
	Dp <- Ktrans(p)
	An.inv <- solve(An)
	w_mat <- 0.5 * t(Dp) %*% kronecker(An.inv, An.inv) %*% Dp
	rownames(w_mat) <- vech(matrix.names(sem.obj$C))
	colnames(w_mat) <- rownames(w_mat)
	return(w_mat)
}

#delta_matrix <- function (sem.object, adj = 1e-04) 
#{
#	p.star <- sem.object$n * (sem.object$n + 1)/2
#	par.idx <- which(sem.object$ram[, "parameter"] > 0)
#	nparams <- length(par.idx)
#	delta.mat <- matrix(0, nparams, p.star)
#	rownames(delta.mat) <- rep(NA, nparams)
#	colnames(delta.mat) <- vech(matrix.names(sem.object$C))
#	vars <- sem.object$var.names
#	C.vect <- vech(sem.object$C)
#	J <- sem.object$J
#	m <- sem.object$m
#	for (i in 1:nparams) {
#		A <- sem.object$A
#		P <- sem.object$P
#		from <- sem.object$ram[par.idx[i], 2]
#		to <- sem.object$ram[par.idx[i], 3]
#		path_type <- sem.object$ram[par.idx[i], 1]
#		if (path_type == 1) {
#			adjust <- abs(A[from, to])*adj
#			A[from, to] <- A[from, to] + adjust
#		}
#		else {
#			adjust <- abs(P[from, to])*adj
#			P[from, to] <- P[from, to] + adjust
#			P[to, from] <- P[to, from]
#		}
#		I.Ainv <- solve(diag(m) - A)
#		C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
#		C.vech <- vech(C)
#		delta.mat[i, ] <- (C.vech - C.vect)/adjust
#		rownames(delta.mat)[i] <- rownames(sem.object$ram)[par.idx[i]]
#		if (rownames(delta.mat)[i] == "") {
#			rownames(delta.mat)[i] <- paste(vars[sem.object$ram[par.idx[i], 
#						3]], vars[sem.object$ram[par.idx[i], 2]], sep = "-")
#		}
#	}
#	delta.mat <- t(delta.mat)
#	return(delta.mat)
#}

delta_matrix <- function (sem.object, adj = 1e-04) {
	p.star <- sem.object$n * (sem.object$n + 1)/2
	pars <- names(sem.object$coeff)
	nparms <- length(pars)
	delta.mat <- matrix(0, nparms, p.star)
	rownames(delta.mat) <- pars
	vars <- sem.object$var.names
	J <- sem.object$J
	m <- sem.object$m
	for (j in 1:nparms) {
		A <- sem.object$A
		P <- sem.object$P
		i <- which(rownames(sem.object$ram) == pars[j])
		from <- sem.object$ram[i, 2]
		to <- sem.object$ram[i, 3]
		path_type <- sem.object$ram[i, 1][1]
		if (path_type == 1){
			AA <- A[cbind(from, to)][1]
			adjust <- abs(AA) * adj
			A[cbind(from, to)] <- AA + adjust
		}
		else {
			PP <- P[cbind(to, from)][1]
			adjust <- PP * adj
			P[cbind(from, to)] <- P[cbind(to, from)] <- PP + adjust
		}
		I.Ainv <- solve(diag(m) - A)
		C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
		delta.mat[j, ] <- (vech(C) - vech(sem.object$C))/adjust
	}
	t(delta.mat)
}


aicc.adjchisq <- function (adj.obj) 
{
	t <- adj.obj$t
	N <- adj.obj$N
	ret <- adj.obj$chisq.scaled + 2 * t * (t + 1)/(N - t - 1)
	return(ret)
}

bic.adjchisq <- function (adj.obj) 
{
	ret <- adj.obj$chisq.scaled - adj.obj$df * log(adj.obj$N)
	return(ret)
}

sem.props <- function (object) 
{
	ret <- list()
	ret$N <- object$N
	N <- ret$N
	ret$n <- object$n
	n <- ret$n
	ret$t <- object$t
	t <- ret$t
	ret$n.fix <- object$n.fix
	n.fix <- ret$n.fix
	ret$df <- n * (n + 1)/2 - t - n.fix * (n.fix + 1)/2
	ret$chisq <- object$criterion * (N - (!object$raw))
	return(ret)
}

Ktrans <- function (num.vars) 
{
	D.mat <- duplication.matrix(num.vars)
	return(D.mat)
}

matrix.names <- function (mat, sep = "_") 
{
	rnames <- rownames(mat)
	cnames <- colnames(mat)
	for (i in 1:length(mat[, 1])) {
		for (j in 1:length(mat[, 1])) {
			mat[i, j] <- paste(rnames[i], cnames[j], sep = sep)
		}
	}
	return(mat)
}

sem_hessian <- function (w_mat, delta_mat) 
{
	ret <- t(delta_mat) %*% w_mat %*% delta_mat
	return(ret)
}
