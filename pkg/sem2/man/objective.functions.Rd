\name{objective.functions}
\alias{objectiveML}
\alias{objectiveGLS}

\title{
sem Objective Function Builders
}

\description{
These functions return objective functions suitable for use with optimizers called by \code\link{{sem}}. The user would not
normally call these functions directly, but rather merely supply one of them in the \code{objective} argument to
\code{sem}.
}

\usage{
objectiveML(gradient = TRUE)

objectiveGLS(gradient = FALSE)
}

\arguments{
  \item{gradient}{If \code{TRUE}, the object that's returned includes a function for computing an analytic gradient; there is at present no
  analytic gradient available for \code{objectiveGLS}.}
}

\value{
These functions return an object of class \code{"semObjective"}, with up to two elements:
\item{objective}{an objective function.}
\item{gradient}{a gradient function.}
}

\author{John Fox \email{jfox@mcmaster.ca}

\seealso{
\code{\link{sem}}, \code{\link{optimizers}}
}

\keyword{models}

