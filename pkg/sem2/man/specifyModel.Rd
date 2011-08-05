\name{specifyModel}
\alias{specifyModel}
\alias{print.semmod}
\title{Specify a Structural Equation Model}
\description{
  Create the RAM specification of a structural equation model.
}
\usage{
specifyModel(file="", exog.variances=FALSE, endog.variances=TRUE, covs)

classifyVariables(model)
removeRedundantPaths(model, warn=TRUE)

\method{print}{semmod}(x, ...)
}

\arguments{
  \item{file}{The (quoted) file from which to read the model specification,
    including the path to the file if it is not in the current directory. If
    \code{""} (the default), then the specification is read from the standard
    input stream, and is terminated by a blank line.}
  \item{exog.variances}{If \code{TRUE} (the default is \code{FALSE}), free variance
  	parameters are added for the exogenous variables that lack them.}
  \item{endog.variances}{If \code{TRUE} (the default), free error-variance parameters
    are added for the endogenous variables that lack them.}
  \item{covs}{optional: a character vector of one or more elements, with each element
  	giving a string of variable names, separated by commas. Variances and covariances
  	among all variables in each such string are added to the model.}
  \item{x, model}{An object of class \code{semmod}, as produced by
    \code{specifyModel}.}
  \item{warn}{print a warning if redundant paths are detected.}
  \item{...}{Ignored.}
}

\details{
  Each line of the RAM specification consists of three (unquoted) entries, 
  separated by commas:
  
    \describe{
      \item{1. Arrow specification:}{This is a simple formula, of the form
        \code{A -> B} or, equivalently, \code{B <- A} for a regression
        coefficient (i.e., a single-headed or directional arrow); 
        \code{A <-> A} for a variance or \code{A <-> B} for a covariance
        (i.e., a double-headed or bidirectional arrow). Here, \code{A} and
        \code{B} are variable names in the model. If a name does not correspond
        to an observed variable, then it is assumed to be a latent variable.
        Spaces can appear freely in an arrow specification, and
        there can be any number of hyphens in the arrows, including zero: Thus,
        e.g., \code{A->B}, \code{A --> B}, and \code{A>B} are all legitimate
        and equivalent.}
      \item{2. Parameter name:}{The name of the regression coefficient, variance,
        or covariance specified by the arrow. Assigning the same name to two or
        more arrows results in an equality constraint. Specifying the parameter name
        as \code{NA} produces a fixed parameter.}
      \item{3. Value:}{start value for a free parameter or value of a fixed parameter.
        If given as \code{NA}, \code{sem} will compute the start value.}
      }
      
  Lines may end in a comment following \code{\#}.
      
  See \code{\link{sem}} for further details on model specification.
}

\value{
  \code{specifyModel} and \code{removeRedundantPaths} return an object of class \code{semmod}, suitable as input for \code{\link{sem}}.
  \code{classifyVariables} returns a list with two character vectors: \code{endogenous}, containing the names of endogenous
  variables in the model; and \code{exogenous}, containing the names of exogenous variables.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{sem}}}

\examples{

# Note: These examples can't be run via example() because the default file
#  argument of specifyModel() requires that the model specification be entered
#  at the command prompt. The examples can be copied and run in the R console,
#  however.

    \dontrun{
model.dhp <- specifyModel()
    RParAsp  -> RGenAsp, gam11,  NA
    RIQ      -> RGenAsp, gam12,  NA
    RSES     -> RGenAsp, gam13,  NA
    FSES     -> RGenAsp, gam14,  NA
    RSES     -> FGenAsp, gam23,  NA
    FSES     -> FGenAsp, gam24,  NA
    FIQ      -> FGenAsp, gam25,  NA
    FParAsp  -> FGenAsp, gam26,  NA
    FGenAsp  -> RGenAsp, beta12, NA
    RGenAsp  -> FGenAsp, beta21, NA
    RGenAsp  -> ROccAsp,  NA,     1
    RGenAsp  -> REdAsp,  lam21,  NA
    FGenAsp  -> FOccAsp,  NA,     1
    FGenAsp  -> FEdAsp,  lam42,  NA
    RGenAsp <-> RGenAsp, ps11,   NA
    FGenAsp <-> FGenAsp, ps22,   NA
    RGenAsp <-> FGenAsp, ps12,   NA
    ROccAsp <-> ROccAsp, theta1, NA
    REdAsp  <-> REdAsp,  theta2, NA
    FOccAsp <-> FOccAsp, theta3, NA
    FEdAsp  <-> FEdAsp,  theta4, NA
    
model.dhp
    
# an equivalent specification, allowing specifyModel() to generate
#  variance parameters for endogenous variables:
 
model.dhp <- specifyModel()
RParAsp  -> RGenAsp, gam11,  NA
RIQ      -> RGenAsp, gam12,  NA
RSES     -> RGenAsp, gam13,  NA
FSES     -> RGenAsp, gam14,  NA
RSES     -> FGenAsp, gam23,  NA
FSES     -> FGenAsp, gam24,  NA
FIQ      -> FGenAsp, gam25,  NA
FParAsp  -> FGenAsp, gam26,  NA
FGenAsp  -> RGenAsp, beta12, NA
RGenAsp  -> FGenAsp, beta21, NA
RGenAsp  -> ROccAsp,  NA,     1
RGenAsp  -> REdAsp,  lam21,  NA
FGenAsp  -> FOccAsp,  NA,     1
FGenAsp  -> FEdAsp,  lam42,  NA
RGenAsp <-> FGenAsp, ps12,   NA

model.dhp

# Another equivalent specification, telling specifyModel to add paths for 
#   variances and covariance of RGenAsp and FGenAsp:
 
model.dhp <- specifyModel(covs="RGenAsp, FGenAsp")
RParAsp  -> RGenAsp, gam11,  NA
RIQ      -> RGenAsp, gam12,  NA
RSES     -> RGenAsp, gam13,  NA
FSES     -> RGenAsp, gam14,  NA
RSES     -> FGenAsp, gam23,  NA
FSES     -> FGenAsp, gam24,  NA
FIQ      -> FGenAsp, gam25,  NA
FParAsp  -> FGenAsp, gam26,  NA
FGenAsp  -> RGenAsp, beta12, NA
RGenAsp  -> FGenAsp, beta21, NA
RGenAsp  -> ROccAsp,  NA,     1
RGenAsp  -> REdAsp,  lam21,  NA
FGenAsp  -> FOccAsp,  NA,     1
FGenAsp  -> FEdAsp,  lam42,  NA

model.dhp

classifyVariables(model.dhp)
    }
}

\keyword{models}
