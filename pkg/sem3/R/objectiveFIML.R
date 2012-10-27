# last modified 2012-01-06 by J. Fox
# Modified for Compiled code in C/C++ by Zhenghua Nie.

objectiveFIML <- function (gradient = TRUE, hessian=FALSE) 
{
    result <- list(objective = function(par, model.description) {
        with(model.description, {
						 
						 res <- CompiledObjective(par=par, model.description=model.description, objective="objectiveFIML", gradient=gradient, hessian=hessian) 
				
						 f <- res$f
				     C <- res$C
				     A <- res$A
				     P <- res$P

						 grad <- NULL
						 if(gradient)
								 grad <- res$gradient
						 hess <- NULL
						 if(hessian) 
								 hess <- res$hessian

						 attributes(f) <- list(C = C, A = A, P = P, gradient=grad, hessian=hess)
						 f
}
				)
}
		)
		class(result) <- "semObjective"
		result
}
