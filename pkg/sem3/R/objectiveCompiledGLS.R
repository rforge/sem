# last modified 2011-11-04 by J. Fox

objectiveCompiledGLS <- function (gradient = FALSE) 
{
    result <- list(objective = function(par, model.description) {
        with(model.description, {
						 
						 res <- CompiledObjective(par=par, model.description=model.description, objective="objectiveGLS") 
				
						 f <- res$f
				     C <- res$C
				     A <- res$A
				     P <- res$P

            attributes(f) <- list(C = C, A = A, P = P)
            f
        }
				)
    }
		)
    class(result) <- "semObjective"
    result
}
