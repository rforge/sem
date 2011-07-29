# methods to generate 
# various information criteria from
# sem or adjchisq objects
# as well as generate and AIC table

# last modified 2011-07-29 by J. Fox

logLik.sem <- function(object, ...){
	deviance(object)*-2	
	}

# generics

AICc <- function(object) UseMethod ("AICc", object)

CAIC <- function(object) UseMethod ("CAIC", object)


# methods for sem objects

AIC.sem<-function(object, ..., k) {
	deviance(object) + 2*object$t
}

# small sample second order corrected aic
AICc.sem<-function(object) {
	deviance(object) + 2*object$t*(object$t + 1)/(object$N - object$t - 1)
}

# Consistent Akaike Information Criterion
CAIC.sem<-function(object) {
	props<-sem.props(object)
	props$chisq - props$df*(1 + log(object$N))
}

BIC.sem<-function(object, ...) {
	deviance(object) + object$t*log(object$N)
}

# weights

aicW <- function(a.list, func=AICc){
	aiclist <- sapply(a.list, function(x) eval(func(x)), simplify=TRUE)			
	delta.i <- aiclist - min(aiclist)
	aicw <- exp(-0.5*delta.i)/sum(exp(-0.5*delta.i))
	return.matrix <- matrix(c(aiclist,delta.i, aicw), ncol=3)
	colnames(return.matrix) <- c("IC", "delta.i", "weight")
	rownames(return.matrix) <- 1:length(return.matrix[,1])
	return(return.matrix)
}
