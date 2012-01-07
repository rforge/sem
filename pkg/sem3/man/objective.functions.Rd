\name{objective.functions}

\alias{objective.functions}
\alias{objectiveML}
\alias{objectiveGLS}
\alias{objectiveML2}
\alias{objectiveGLS2}
\alias{msemObjectiveML}

\title{
sem Objective-Function Builders
}

\description{
These functions return objective functions suitable for use with \code{\link{optimizers}} called by \code{\link{sem}}. The user would not
normally call these functions directly, but rather supply one of them in the \code{objective} argument to
\code{sem}. Users may also write their own objective functions. \code{objectiveML} and \code{objectiveML2} are for multinormal maximum-likelihood
estimation; \code{objectiveGLS} and \code{objectiveGLS2} are for generalized least squares. \code{objectiveML} and \code{objectiveGLS} use
compiled code and are therefore substantially faster. \code{objectiveML2} and \code{objectiveGLS2} are provided primarily to illustrate
how to write \code{sem} objective functions in R. \code{msemObjectiveML} is for fitting multi-group models by multinormal maximum likelihood.
}

\usage{
objectiveML(gradient=TRUE, hessian=FALSE)
objectiveML2(gradient=TRUE)

objectiveGLS(gradient=FALSE)
objectiveGLS2(gradient=FALSE)

msemObjectiveML(gradient=TRUE)
}

\arguments{
  \item{gradient}{If \code{TRUE}, the object that's returned includes a function for computing an analytic gradient; there is at present no
  analytic gradient available for \code{objectiveGLS} or \code{objectiveGLS2}.}
  \item{hessian}{If \code{TRUE}, the objected returned includes a function to compute an analytic Hessian; only avaiable for \code{objectiveML}
  and not generally recommended.}
}

\value{
These functions return an object of class \code{"semObjective"}, with up to two elements:
\item{objective}{an objective function.}
\item{gradient}{a gradient function.}
}

\author{John Fox \email{jfox@mcmaster.ca}}

\references{
See \code{\link{sem}}.
}

\seealso{
\code{\link{sem}}, \code{\link{optimizers}}
}

\keyword{models}

