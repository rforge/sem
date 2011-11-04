\name{residuals.sem}
\alias{residuals.sem}
\alias{standardizedResiduals}
\alias{standardizedResiduals.sem}
\alias{normalizedResiduals}
\alias{normalizedResiduals.objectiveML}
\alias{normalizedResiduals.objectiveGLS}

\title{Residual Covariances for a Structural Equation Model}
\description{
  These functions compute residual covariances, variance-standardized
  residual covariances, and normalized residual covariances
  for the observed variables in a structural-equation
  model fit by \code{sem}.
}
\usage{
\method{residuals}{sem}(object, ...)

\method{standardizedResiduals}{sem}(object, ...)

\method{normalizedResiduals}{objectiveML}(object, ...)

\method{normalizedResiduals}{objectiveGLS}(object, ...)
}

\arguments{
  \item{object}{an object inheriting from class \code{sem} returned by the \code{\link{sem}} function.}
  \item{...}{not for the user.}
}

\details{
  Residuals are defined as \eqn{S - C}, where \eqn{S} is the sample covariance matrix
  of the observed variables and \eqn{C} is the model-reproduced covariance matrix.
  The \emph{standardized} residual covariance for a pair of variables divides the
  residual covariance by the product of the sample standard deviations of the 
  two variables, \eqn{(s_{ij} - c_{ij})/(s_{ii}s_{jj})^{1/2}}. The \emph{normalized} residual
  is given by
  \deqn{\frac{s_{ij}-c_{ij}} {[(c_{ii}c_{jj}-c_{ij}^2)/N^{*}]^{1/2}}}{%
    (s[ij] - c[ij])/[(c[ii]c[ii] + c[ij]^2)/N*]^[1/2]}
  where \eqn{N^{*}} is the number of observations minus one if the model is fit to a
  covariance matrix, or the number of observations if it is fit to a raw moment matrix.
}

\value{
  Each function returns a matrix of residuals.
}

\references{
    Bollen, K. A. (1989) 
    \emph{Structural Equations With Latent Variables.} Wiley.
    }
    
\author{John Fox \email{jfox@mcmaster.ca}}


\seealso{\code{\link{sem}}}

\examples{
    \dontrun{
# ------------- assumes that Duncan, Haller and Portes peer-influences model
# -------------     has been fit and is in sem.dhp.1

residuals(sem.dhp.1)
normalizedResiduals(sem.dhp.1) 
standardizedResiduals(sem.dhp.1)
    }
}

\keyword{models}
