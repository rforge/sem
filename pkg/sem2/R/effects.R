# added 2011-10-31 by J. Fox

effects.sem <- function(object, ...) {
	A <- object$A
	m <- object$m
	I <- diag(m)
	endog <- classifyVariables(object$semmod)$endogenous  
	AA <- - A
	diag(AA) <- 1
	Total <- solve(AA) - I
	Indirect <-  Total - A
	result <- list(Total=Total[endog, ], Direct=A[endog, ], Indirect=Indirect[endog, ])
	class(result) <- "semeffects"
	result
}

print.semeffects <- function(x, digits=getOption("digits"), ...){
	cat("\n Total Effects\n")
	print(x$Total, digits=digits)
	cat("\n Direct Effects\n")
	print(x$Direct, digits=digits)
	cat("\n Indirect Effects\n")
	print(x$Indirect, digits=digits)
	invisible(x)
}
