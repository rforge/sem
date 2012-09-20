\name{miSem}
\alias{miSem}
\alias{miSem.semmod}
\alias{print.miSem}
\alias{summary.miSem}

\title{
Estimate a Structural Equation Model By Multiple Imputation
}

\description{
\code{miSem} uses the \code{\link[mi]{mi}} function in the \pkg{mi} package to generate multiple imputations of missing
data, fitting the specified model to each completed data set.
}

\usage{
miSem(model, ...)

\method{miSem}{semmod}(model, ..., data, formula = ~., raw = FALSE, 
    fixed.x = NULL, objective = objectiveML, n.imp = 5, n.iter = 30, 
    seed = sample(1e+06, 1), mi.args = list())

\method{print}{miSem}(x, ...)

\method{summary}{miSem}(object, digits=max(3, getOption("digits") - 2), ...)
}

\arguments{
  \item{model}{
    an SEM model-description object of class \code{semmod}, created by \code{\link{specifyEquations}}
    \code{\link{cfa}}, or \code{\link{specifyModel}}.
    }
  \item{\dots, formula, raw, fixed.x, objective}{
    arguments to be passed to \code{\link{sem}}.
    }
  \item{data}{
    an R data frame, presumably with some missing data (encoded as \code{NA}), containing the data for fitting
    the SEM, possibly along with other variables to use to obtain multiple imputations of missing values.
    }

  \item{n.imp}{
    number of imputations (default \code{5}).
    }

  \item{n.iter}{
    number of iterations for the multiple-imputation process (default \code{30}).
    }

  \item{seed}{
    seed for the random-number generator (default is an integer sampled from 1 to 1E6); stored in the resulting object.
    }

  \item{mi.args}{
    other arguments to be passed to \code{\link[mi]{mi}}.
    }
  \item{x, object}{
    an object of class \code{"miSem"}.
    }
  \item{digits}{
    for printing numbers.
  }
}

\value{
    \code{miSem} returns an object of class \code{"miSem"} with the following components:
    \item{initial.fit}{an \code{sem} model object produced using \code{\link{objectiveFIML}} if \code{raw=TRUE},
    or the objective function given by the \code{objective} argument otherwise.}
    \item{mi.fits}{a list of \code{sem} model objects, one for each imputed data set.}
    \item{imputation}{the object produced by \code{\link[mi]{mi}}, containing the imputed data sets.}
    \item{seed}{the seed used for the random number generator.}
}

\references{
    Yu-Sung Su, Andrew Gelman, Jennifer Hill, Masanao Yajima. (2011). 
    ``Multiple imputation with diagnostics (mi) in R: Opening windows into the black box.'' 
    \emph{Journal of Statistical Software} 45(2).
}
\author{
    John Fox \email{jfox@mcmaster.ca}
    }

\seealso{
\code{\link{sem}}, \code{\link[mi]{mi}}
}
\examples{
    \dontrun{
mod.cfa.tests <- cfa(raw=TRUE)
verbal: x1, x2, x3
math: y1, y2, y3

imps <- miSem(mod.cfa.tests, data=Tests, fixed.x="Intercept", raw=TRUE, seed=12345)    
summary(imps, digits=3)    
    }
}

\keyword{models}
