sem <- function(model, ...){
	if (is.character(model)) class(model) <- 'semmod'
	UseMethod('sem', model)
}

sem.semmod <- function (model, S, N, data, raw=FALSE, obs.variables=rownames(S), fixed.x=NULL, formula= ~ ., robust=TRUE, debug=FALSE, ...){
    parse.path <- function(path) {                                           
        path.1 <- gsub('-', '', gsub(' ','', path))
        direction <- if (regexpr('<>', path.1) > 0) 2 
            else if (regexpr('<', path.1) > 0) -1
            else if (regexpr('>', path.1) > 0) 1
            else stop(paste('ill-formed path:', path))
        path.1 <- strsplit(path.1, '[<>]')[[1]]
        list(first=path.1[1], second=path.1[length(path.1)], direction=direction)
        }
	if (missing(S)){
		if (missing(data)) stop("S and data cannot both be missing")
		N.all <- nrow(data)
		data <- na.omit(data)
		data <- as.matrix(data)
		if (!is.numeric(data)) stop("data must be numeric")
		S <- if (raw) rawMoments(formula, data=data) else cov(data)
		N <- nrow(data)
		if (N < N.all) warning(N - N.all, " observations removed due to missingness")
	}
    if ((!is.matrix(model)) | ncol(model) != 3) stop ('model argument must be a 3-column matrix')
    startvalues <- as.numeric(model[,3])
    par.names <- model[,2]
    n.paths <- length(par.names)
    heads <- from <- to <- rep(0, n.paths)
    for (p in 1:n.paths){
        path <- parse.path(model[p,1])
        heads[p] <- abs(path$direction)
        to[p] <- path$second
        from[p] <- path$first
        if (path$direction == -1) {
            to[p] <- path$first
            from[p] <- path$second
            }
        }
    ram <- matrix(0, p, 5)
    all.vars <- unique(c(to, from))
    latent.vars <- setdiff(all.vars, obs.variables)
    not.used <- setdiff(obs.variables, all.vars)
    if (length(not.used) > 0){
        rownames(S) <- colnames(S) <- obs.variables
        obs.variables <- setdiff(obs.variables, not.used)
        S <- S[obs.variables, obs.variables]
        warning("The following observed variables are in the input covariance or raw-moment matrix ",
            "but do not appear in the model:\n",
            paste(not.used, collapse=", "), "\n")
        }
    vars <- c(obs.variables, latent.vars)
    pars <- na.omit(unique(par.names))
    ram[,1] <- heads
    ram[,2] <- apply(outer(vars, to, '=='), 2, which)
    ram[,3] <- apply(outer(vars, from, '=='), 2, which)   
    par.nos <- apply(outer(pars, par.names, '=='), 2, which)
    if (length(par.nos) > 0)
        ram[,4] <- unlist(lapply(par.nos, function(x) if (length(x) == 0) 0 else x))
    ram[,5]<- startvalues
    colnames(ram) <- c('heads', 'to', 'from', 'parameter', 'start')
    if (!is.null(fixed.x)) fixed.x <- apply(outer(vars, fixed.x, '=='), 2, which)
    n <- length(obs.variables)
    m <- length(all.vars)
    t <- length(pars)
    if (debug) {
        cat('\n observed variables:\n') 
        print(paste(paste(1:n,':', sep=''), obs.variables, sep=''))
        cat('\n')
        if (m > n){ 
            cat('\n latent variables:\n')
            print(paste(paste((n+1):m,':', sep=''), latent.vars, sep=''))
            cat('\n')
            }
        cat('\n parameters:\n') 
        print(paste(paste(1:t,':', sep=''), pars, sep=''))
        cat('\n\n RAM:\n')
        print(ram)
        }
	if (missing(data)) data <- NULL
    sem(ram, S=S, N=N, raw=raw, data=data, param.names=pars, var.names=vars, fixed.x=fixed.x,
		semmod=model, robust=robust, debug=debug, ...)
    }
     
sem.default <- function(model, S, N, data=NULL, raw=FALSE, param.names, 
    var.names, fixed.x=NULL, robust=TRUE, semmod=NULL, debug=FALSE,
    analytic.gradient=TRUE, warn=FALSE, maxiter=500, par.size=c('ones', 'startvalues'), 
    refit=TRUE, start.tol=1E-6, optimizer=optimizerNlm, objective=objectiveML, ...){
    ord <- function(x) 1 + apply(outer(unique(x), x, "<"), 2, sum)
    is.triangular <- function(X) {
        is.matrix(X) && (nrow(X) == ncol(X)) && 
            (all(0 == X[upper.tri(X)])) || (all(0 == X[lower.tri(X)]))
        } 
	ram <- model
    S <- unclass(S) # in case S is a rawmoment object
    if (nrow(S) > 1 && is.triangular(S)) S <- S + t(S) - diag(diag(S))
    if (!isSymmetric(S)) stop('S must be a square triangular or symmetric matrix')
	if (qr(S)$rank < ncol(S)) warning("S is numerically singular: expect problems")
	if (any(eigen(S, symmetric=TRUE, only.values=TRUE)$values <= 0)) 
		warning("S is not positive-definite: expect problems")
    if ((!is.matrix(ram)) | ncol(ram) != 5 | (!is.numeric(ram)))
        stop ('ram argument must be a 5-column numeric matrix')
    par.size <- if (missing(par.size)) {
        range <- range(diag(S))
        if (range[2]/range[1] > 100) 'startvalues' else 'ones'
        }
        else match.arg(par.size)
    n <- nrow(S)
    observed <- 1:n
    n.fix <- length(fixed.x)
    if (!is.null(fixed.x)){
        for (i in 1:n.fix){
            for (j in 1:i){
                ram <- rbind(ram, c(2, fixed.x[i], fixed.x[j], 
                    0, S[fixed.x[i], fixed.x[j]]))
                }
            }
        }
    m <- max(ram[,c(2,3)])
    missing.variances <- setdiff(1:m, ram[,2][ram[,2] == ram[,3]])
    if (length(missing.variances) > 0) warning(paste("The following variables have no variance or error-variance parameter (double-headed arrow):\n",
        paste(var.names[missing.variances], collapse=", "), 
        "\nThe model is almost surely misspecified; check also for missing covariances.\n"))
    t <- max(ram[,4])
    df <- n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2
    if (df < 0) stop(paste("The model has negative degrees of freedom =", df))
    J <- matrix(0, n, m)
    correct <- matrix(2, m, m)
    diag(correct) <- 1
    J[cbind(1:n, observed)]<-1
    par.posn <-  sapply(1:t, function(i) which(ram[,4] == i)[1])
    colnames(ram)<-c("heads", "to", "from", "parameter", "start value")
    rownames(ram)<-rep("",nrow(ram))
    if (length(param.names) > 0) rownames(ram)[par.posn]<-param.names
    fixed <- ram[,4] == 0
    sel.free <- ram[,4]
    sel.free[fixed] <- 1
    one.head <- ram[,1] == 1
    one.free <- which( (!fixed) & one.head )
    two.free <- which( (!fixed) & (!one.head) )
    arrows.1 <- ram[one.head, c(2,3), drop=FALSE]
    arrows.2 <- ram[!one.head, c(2,3), drop=FALSE]
    arrows.2t <- ram[!one.head, c(3,2), drop=FALSE]
    arrows.1.free <- ram[one.free,c(2,3), drop=FALSE]
    arrows.2.free <- ram[two.free,c(2,3), drop=FALSE]
    sel.free.1 <- sel.free[one.free]
    sel.free.2 <- sel.free[two.free]
    unique.free.1 <- unique(sel.free.1)
    unique.free.2 <- unique(sel.free.2)
    result <- list()
    result$var.names <- var.names
    result$ram <- ram
    rownames(S) <- colnames(S) <- var.names[observed]
    result$S <- S
    result$J <- J
    result$n.fix <- n.fix
    result$n <- n
    result$N <- N
    result$m <- m
    result$t <- t
    result$raw <- raw
	result$data <- data
	result$semmod <- semmod
    if (length(param.names)== 0){
        warning("there are no free parameters in the model")
        }
    else {
        start <- if (any(is.na(ram[,5][par.posn])))
            startvalues(S, ram, debug=debug, tol=start.tol)
            else ram[,5][par.posn]
      typsize <- if (par.size == 'startvalues') abs(start) else rep(1,t)
		model.description <- list(S=S, logdetS=log(det(S)), invS=solve(S), N=N, m=m, n=n, t=t, fixed=fixed, ram=ram, sel.free=sel.free, arrows.1=arrows.1, arrows.1.free=arrows.1.free,
			one.head=one.head, arrows.2t=arrows.2t, arrows.2=arrows.2, arrows.2.free=arrows.2.free, unique.free.1=unique.free.1, unique.free.2=unique.free.2,
			J=J, correct=correct, param.names=param.names, var.names=var.names, observed=observed, raw=raw)
		res <- optimizer(start=start, 
			objective=objective, gradient=analytic.gradient, maxiter=maxiter, debug=debug, par.size=par.size, 
			model.description=model.description, warn=warn, ...)
        ram[par.posn, 5] <- start
        par.code <- paste(var.names[ram[,2]], c('<---', '<-->')[ram[,1]],
        var.names[ram[,3]])
        result$coeff <- res$par
		result$vcov <- res$vcov
        result$par.posn <- par.posn
        result$convergence <- res$convergence
        result$iterations <- res$iterations
        result$raw <- raw
    result$criterion <-  res$criterion
	result$C <- res$C
	result$A <- res$A
	result$P <- res$P
	}
	cls <- gsub("\\.", "", deparse(substitute(objective)))
    class(result) <- c(cls, "sem")
	if (robust && !is.null(data) && inherits(result, "objectiveML")){
		result$adj.obj <- sbchisq(result, data)
		result$robust.vcov <- robust_vcov(result, adj.obj=result$adj.obj)
    }
	result
}
	
vcov.sem <- function(object, robust=FALSE, ...) {
	if (robust) object$robust.vcov else object$vcov
}

coef.sem <- function(object, ...){
	object$coeff
}