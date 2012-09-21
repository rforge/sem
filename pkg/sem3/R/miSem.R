# last modified 2012-09-20 by J. Fox

miSem <- function(model, ...){
    if (!require(mi)) stop("mi package missing")
    UseMethod("miSem")
}

miSem.semmod <- function(model, ..., data, formula = ~., raw=FALSE, fixed.x=NULL, objective=objectiveML,
                  n.imp=5, n.iter=30, seed=sample(1e6, 1), mi.args=list()){
    cls <- gsub("\\.", "", deparse(substitute(objective)))
    cls <- gsub("2", "", cls)
    cls <- c(cls, "sem")
    warn <- options(warn=-1)
    on.exit(options(warn))
    initial.fit <- sem(model, ..., data=data, formula=formula, raw=raw, fixed.x=fixed.x,
                       objective = if (raw) objectiveFIML else objective)
    options(warn)
    class(initial.fit) <- if (raw) c("objectiveFIML", "sem") else cls
    coefficients <- coefficients(initial.fit)
    coef.names <- names(coefficients)
    var.names <- initial.fit$var.names 
    ram <- initial.fit$ram 
    ram[coef.names, "start value"] <- coefficients 
    N <- nrow(data)
    if (!is.null(fixed.x)) fixed.x <- apply(outer(var.names, fixed.x, "=="), 2, which)
    mi.args$n.imp <- n.imp
    mi.args$n.iter <- n.iter
    mi.args$seed <- seed
    mi.args$object <- data
    mi.data <- do.call("mi", mi.args)
    fits <- vector(n.imp, mode="list")
    has.tcltk <- require("tcltk")
	if (has.tcltk) pb <- tkProgressBar("Fitting", "Imputation no.: ", 0, n.imp)
    for (i in 1:n.imp){
        if (has.tcltk) setTkProgressBar(pb, i, label=sprintf("Imputation no.: %d", i))
        data.i <- mi.data.frame(mi.data, m=i)
        data.i <- model.frame(formula, data=data.i)
        data.i <- model.matrix(formula, data=data.i)
        colnames(data.i)[colnames(data.i) == "(Intercept)"] <- "Intercept"
    	S <- if (raw) rawMoments(data.i) else {
			data.i <-  data.i[, colnames(data.i) != "Intercept"]
			cov(data.i)
		}
        fit <- sem(ram, S=S, N=N, data=data.i, raw=raw, param.names=coef.names, var.names=var.names, fixed.x=fixed.x,
                              optimizer=initial.fit$optimizer, objective=objective, ...)
        class(fit) <- cls   
	    fits[[i]] <- fit
    }
    if (has.tcltk) close(pb)
    result <- list(initial.fit=initial.fit, mi.fits=fits, imputations=mi.data, seed=seed)
    class(result) <- "miSem"
    result
}

print.miSem <- function(x, ...){
    coefs <- lapply(x$mi.fits, coef)
    ses <- lapply(x$mi.fits, function(x) sqrt(diag(vcov(x))))
    pooled <- mi.pooled(coefs, ses)
    table <- matrix(0, length(pooled$coefficients), 4)
    table[, 1] <- pooled$coefficients
    table[, 2] <- pooled$se
    table[, 3] <- table[, 1]/table[, 2]
    table[, 4] <- 2*pnorm(abs(table[, 3]), lower.tail=FALSE)
    rownames(table) <- names(pooled$coefficients)
    cat("\nCoefficients:\n")
    colnames(table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    printCoefmat(table, ...)
    invisible(x)
}

summary.miSem <- function(object, digits=max(3, getOption("digits") - 2), ...){
    coefs <- lapply(object$mi.fits, coef)
    table <- do.call("cbind", coefs)
    rownames(table) <- names(coefs[[1]])
    table <- cbind(table, rowMeans(table), coef(object$initial.fit))
    colnames(table) <- c(paste("Imputation", 1:length(coefs)), "Averaged", "Initial Fit")
    result <- list(object=object, mi.results=table, digits=digits)
    class(result) <- "summary.miSem"
    result
}

print.summary.miSem <- function(x, ...){
    cat("\nCoefficients by imputation:\n")
    print(x$mi.results, digits=x$digits, ...)
    print(x$object, digits=x$digits, ...)
    invisible(x)
}


