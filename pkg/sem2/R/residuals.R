# last modified 2011-08-07 by J. Fox

residuals.sem <- function(object, ...){
    object$S - object$C
    }
	
standardized.residuals <- function(...){
	.Deprecated("standardizedResiduals", package="sem")
	standardizedResiduals(...)
}

standardizedResiduals <- function(object, ...){
    UseMethod("standardizedResiduals")
    }

standardizedResiduals.sem <- function(object, ...){
    res <- residuals(object)
    s <- diag(object$S)
    res/sqrt(outer(s, s))
    }

normalized.residuals <- function(object, ...){
	.Deprecated("normalizedResiduals", package="sem")
	normalizedResiduals(...)
    }
	
normalizedResiduals <- function(object, ...){
	UseMethod("normalizedResiduals")
}
    
normalizedResiduals.sem <- function(object, ...){
    res <- residuals(object)
    N <- object$N - (!object$raw)
    C <- object$C
    c <- diag(C)
    res/sqrt((outer(c,c) + C^2)/N)
    }
