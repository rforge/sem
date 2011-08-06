\name{standardizedCoefficients}
\alias{standardizedCoefficients}
\alias{std.coef}

\title{Standardized Coefficients for Structural Equation Models}
\description{
  These functions calculate standardized regression coefficients
  for structural equation models. The function \code{std.coef} is
  simply an abbreviation for \code{standardized.coefficients}.
}

\usage{
standardizedCoefficients(object, digits = 5, oneheaded=TRUE, twoheaded=TRUE)

stdCoef(...)
}

\arguments{
  \item{object}{an object of class \code{sem} returned by the 
    \code{\link{sem}} function.}
  \item{digits}{number of digits for printed output.}
  \item{oneheaded}{standardize path coefficients? Default is \code{TRUE}.}
  \item{twoheaded}{standardize variances and covariances? Default is \code{TRUE}.}
  \item{...}{arguments to pass to \code{standardized.coefficients}.}
}

\value{
  Returns a data frame with the standardized coefficients, labelled 
  both by parameter names and by arrows in the path diagram for the
  model.
}

\references{
    Bollen, K. A. (1989) 
    \emph{Structural Equations With Latent Variables.} Wiley.
}

\author{John Fox \email{jfox@mcmaster.ca} and Adam Kramer}

\seealso{\code{\link{sem}}}

\examples{

# ------------- assumes that Duncan, Haller and Portes peer-influences model
# -------------     has been fit and is in sem.dhp.1     
\dontrun{
standardizedCoefficients(sem.dhp.1) 
	}
}

\keyword{models}
