# last modified 2012-07-29 by J. Fox

## this is the straightforward approach summing over observations:
# objectiveFIML2 <- function(gradient=FALSE){
#   result <- list(
#     objective = function(par, model.description){
#       with(model.description, {
#         A <- P <- matrix(0, m, m)
#         val <- ifelse (fixed, ram[,5], par[sel.free])
#         A[arrows.1] <- val[one.head]
#         P[arrows.2t] <- P[arrows.2] <- val[!one.head]
#         I.Ainv <- solve(diag(m) - A)
#         C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
#         f <- 0
#         log.2pi <- log(2*pi)
#         for (i in 1:nrow(data)){
#           sel <- valid[i, ]
#           x <- data[i, sel]
#           CC <- C[sel, sel]
#           f <- f + log(det(CC)) + as.vector(x %*% solve(CC) %*% x) + log.2pi*(sum(sel))
#         }
#         f <- f/N
#         attributes(f) <- list(C=C, A=A, P=P)
#         f
#       })
#     }
#   )
#   class(result) <- "semObjective"
#   result
# }

## this avoids duplication by using the svd for the inverse and determinant but is actually slower in R
# objectiveFIML2 <- function(gradient=FALSE){
#   result <- list(
#     objective = function(par, model.description){
#       with(model.description, {
#         A <- P <- matrix(0, m, m)
#         val <- ifelse (fixed, ram[,5], par[sel.free])
#         A[arrows.1] <- val[one.head]
#         P[arrows.2t] <- P[arrows.2] <- val[!one.head]
#         I.Ainv <- solve(diag(m) - A)
#         C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
#         f <- 0
#         log.2pi <- log(2*pi)
#         for (i in 1:nrow(data)){
#           sel <- valid[i, ]
#           x <- data[i, sel]
#           CC <- C[sel, sel]
#           svdCC <-svd(CC)
#           f <- f + log(prod(svdCC$d)) + as.vector(x %*% svdCC$v %*% diag(1/svdCC$d) %*% t(svdCC$u) %*% x) + log.2pi*(sum(sel))
#         }
#         f <- f/N
#         attributes(f) <- list(C=C, A=A, P=P)
#         f
#       })
#     }
#   )
#   class(result) <- "semObjective"
#   result
# }

## this caches results for distinct valid-data patterns but still sums over observations:
# objectiveFIML2 <- function(gradient=FALSE){
#   result <- list(
#     objective = function(par, model.description){
#       with(model.description, {
#         A <- P <- matrix(0, m, m)
#         val <- ifelse (fixed, ram[,5], par[sel.free])
#         A[arrows.1] <- val[one.head]
#         P[arrows.2t] <- P[arrows.2] <- val[!one.head]
#         I.Ainv <- solve(diag(m) - A)
#         C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
#         f <- 0
#         log.2pi <- log(2*pi)
#         CC <- vector(length(unique.patterns), mode="list")
#         names(CC) <- unique.patterns
#         log.det.plus <- CC
#         for (i in 1:nrow(data)){
#           sel <- valid[i, ]
#           if (is.null(log.det.plus[[valid.pattern[i]]])){
#             log.det.plus[[valid.pattern[i]]] <- log(det(C[sel, sel])) + log.2pi*(sum(sel))
#             CC[[valid.pattern[i]]] <- solve(C[sel, sel])
#           }
#           x <- data[i, sel]
#           f <- f + as.vector(x %*% CC[[valid.pattern[i]]] %*% x) + log.det.plus[[valid.pattern[i]]]
#         }
#         f <- f/N
#         attributes(f) <- list(C=C, A=A, P=P)
#         f
#       })
#     }
#   )
#   class(result) <- "semObjective"
#   result
# }

## this sums over distinct valid data patterns rather than observations:
objectiveFIML2 <- function(gradient=FALSE){
  result <- list(
    objective = function(par, model.description){
      with(model.description, {
        A <- P <- matrix(0, m, m)
        val <- ifelse (fixed, ram[,5], par[sel.free])
        A[arrows.1] <- val[one.head]
        P[arrows.2t] <- P[arrows.2] <- val[!one.head]
        I.Ainv <- solve(diag(m) - A)
        C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
        f <- 0
        log.2pi <- log(2*pi)
        n.pat <- nrow(valid.data.patterns)
        for (i in 1:n.pat){
          sel <- valid.data.patterns[i, ]
          X <- data[pattern.number == i, sel, drop=FALSE]
          # note: in the matrix product X %*% solve(C[sel, sel]) %*% t(X))
          #       only the trace is required and so it isn't necessary to form the whole product
          # f <- f + sum(diag(X %*% solve(C[sel, sel]) %*% t(X))) + nrow(X)*(log.2pi + log(det(C[sel, sel])))
          CC <- solve(C[sel, sel])
          for (j in 1:nrow(X)){
            f <- f + as.vector(X[j, ] %*% CC %*% X[j, ])
          }
          f <- f + nrow(X)*(log.2pi + log(det(C[sel, sel])))
        }
        f <- f/N
        attributes(f) <- list(C=C, A=A, P=P)
        f
      })
    }
  )
  class(result) <- "semObjective"
  result
}

logLik.objectiveFIML <- function(object, saturated=FALSE, iterlim=1000, ...){
  logLikSaturated <- function(object, iterlim, ...){
    objective <- function(par){
      C <- matrix(0, n, n)
      C[tri] <- par
      C <- C + t(C) - diag(diag(C))
      f <- 0
      for (i in 1:n.pat){
        sel <- valid.data.patterns[i, ]
        X <- data[pattern.number == i, sel, drop=FALSE]
        f <- f + sum(diag(X %*% solve(C[sel, sel]) %*% t(X))) + nrow(X)*(log.2pi + log(det(C[sel, sel])))
      }
      f
    }
    data <- object$data
    valid <- !is.na(data)
    valid.pattern <- apply(valid, 1, function(row) paste(row, collapse="."))
    unique.patterns <- unique(valid.pattern)
    pattern.number <- apply(outer(valid.pattern, unique.patterns, `==`), 1, which)
    valid.data.patterns <- t(sapply(strsplit(unique.patterns, "\\."), as.logical))
    n.pat <- nrow(valid.data.patterns) 
    log.2pi <- log(2*pi)
    n <- ncol(data)
    N <- nrow(data)
    C <- object$C
    tri <- lower.tri(C, diag=TRUE)
    start <- C[tri]
    opt <- options(warn=-1)
    on.exit(options(opt))
    res <- nlm(objective, start, iterlim=iterlim)
    logL <- - res$minimum/2
    C <- matrix(0, n, n)
    C[tri] <- res$estimate
    C <- C + t(C) - diag(diag(C))
    list(logL=logL, C=C, code=res$code)
  }
  if (saturated) {
    res <- logLikSaturated(object, iterlim=iterlim)
    if (res$code > 3) warning("nlm return code = ", res$code)
    logL <- res$logL
    attr(logL, "C") <- res$C
    return(logL)
  }
  else return(- object$criterion*object$N/2)
}


residuals.objectiveFIML <- function(object, S, ...){
  if (missing(S)) S <- attr(logLik(object, saturated=TRUE), "C")
  S - object$C
}

normalizedResiduals.objectiveFIML <- function(object, S, ...){
  if (missing(S)) S <- attr(logLik(object, saturated=TRUE), "C")
  res <- residuals(object, S=S)
  N <- object$N - (!object$raw)
  C <- object$C
  c <- diag(C)
  res/sqrt((outer(c, c) + C^2)/N)
}

standardizedResiduals.objectiveFIML <- function(object, S, ...){
  if (missing(S)) S <- attr(logLik(object, saturated=TRUE), "C")
  res <- residuals(object, S=S)
  s <- diag(S)
  res/sqrt(outer(s, s))
}

deviance.objectiveFIML <- function(object, saturated.logLik, ...){
  if (missing(saturated.logLik)) saturated.logLik <- logLik(object, saturated=TRUE)
  2*(as.vector(saturated.logLik) - logLik(object))
}

AIC.objectiveFIML <- function(object, saturated.logLik, ..., k) {
  if (missing(saturated.logLik)) saturated.logLik <- logLik(object, saturated=TRUE)
  deviance(object, saturated.logLik) + 2*object$t
}

AICc.objectiveFIML <- function(object, saturated.logLik, ...) {
  if (missing(saturated.logLik)) saturated.logLik <- logLik(object, saturated=TRUE)
  deviance(object, saturated.logLik) + 2*object$t*(object$t + 1)/(object$N - object$t - 1)
}

CAIC.objectiveFIML <- function(object, saturated.logLik, ...) {
  props <- semProps(object)
  if (missing(saturated.logLik)) saturated.logLik <- logLik(object, saturated=TRUE)
  deviance(object, saturated.logLik) - props$df*(1 + log(object$N))
}

BIC.objectiveFIML <- function(object, saturated.logLik, ...) {
  if (missing(saturated.logLik)) saturated.logLik <- logLik(object, saturated=TRUE)
  deviance(object, saturated.logLik) + object$t*log(object$N)
}

summary.objectiveFIML <- function(object, digits=5, conf.level=.90, robust=FALSE, analytic.se=FALSE, saturated=FALSE, saturated.logLik, ...) {
  vcov <- vcov(object, robust=robust, analytic=analytic.se)
  if (any(is.na(vcov))) stop("coefficient covariances cannot be computed")
  if (missing(saturated.logLik)) saturated.logLik <- if (saturated) logLik(object, saturated=TRUE) else NULL
  S <- attr(saturated.logLik, "C")
  norm.res <- if (saturated) normalizedResiduals(object, S) else NA
  se <- sqrt(diag(vcov))
  z <- object$coeff/se
  n.fix <- object$n.fix
  n <- object$n
  t <- object$t
  C <- object$C
  N <- object$N
  df <- n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2
  if (saturated){
    invC <- solve(C)
    CSC <- invC %*% (S - C)
    CSC <- CSC %*% CSC
    CS <- invC %*% S
    CS <- CS %*% CS
    chisq <- 2*(saturated.logLik - logLik(object))
    logLik <- NULL
  }
  else {
    chisq <- NULL
    logLik <- logLik(object)
  }
  Rsq <- RMSEA.U <- RMSEA.L <- RMSEA <- NFI <- NNFI <- CFI <- AGFI <- NA
  RMSEA <- c(RMSEA, RMSEA.L, RMSEA.U, conf.level)
  if (!is.null(object$coeff)){
    var.names <- rownames(object$A)
    ram <- object$ram[object$par.posn, , drop=FALSE]
    par.code <- paste(var.names[ram[,2]], c('<---', '<-->')[ram[,1]],
                      var.names[ram[,3]])
    coeff <- data.frame(object$coeff, se, z, 2*pnorm(abs(z), lower.tail=FALSE), par.code)
    names(coeff) <- c("Estimate", "Std Error", "z value", "Pr(>|z|)", " ")
    row.names(coeff) <- names(object$coeff)
  }
  else coeff <- NULL
  if(saturated){
    AIC <- AIC(object, saturated.logLik=saturated.logLik)
    AICc <- AICc(object, saturated.logLik=saturated.logLik)
    BIC <- BIC(object, saturated.logLik=saturated.logLik)
    CAIC <- CAIC(object, saturated.logLik=saturated.logLik)
    SRMR <- sqrt(sum(standardizedResiduals(object, S=S)^2 * 
      upper.tri(diag(n), diag=TRUE))/(n*(n + 1)/2))
  }
  else AIC <- AICc <- BIC <- CAIC <- SRMR <- NA
#   if (robust) { 
#     chisq.adjusted <- object$adj.obj$chisq.scaled
#     chisqNull.adjusted <- chisqNull/object$adj.obj$c 
#     NFI.adjusted <- (chisqNull.adjusted - chisq)/chisqNull.adjusted
#     NNFI.adjusted <- (chisqNull.adjusted/dfNull - chisq.adjusted/df)/(chisqNull.adjusted/dfNull - 1)
#     L1 <- max(chisq.adjusted - df, 0)
#     L0 <- max(L1, chisqNull.adjusted - dfNull)
#     CFI.adjusted <- 1 - L1/L0
#   }
#   else{
    chisq.adjusted <- chisqNull.adjusted <- NFI.adjusted <- NNFI.adjusted <- CFI.adjusted <- NULL
#   }
  ans <- list(chisq=chisq, logLik=logLik, df=df, chisqNull=chisqNull, dfNull=NA,
              GFI=NULL, AGFI=AGFI, RMSEA=RMSEA, NFI=NFI, NNFI=NNFI, CFI=CFI, BIC=BIC, SRMR=SRMR, 
              AIC=AIC, AICc=AICc, CAIC=CAIC, Rsq=Rsq,
              chisq.adjusted=chisq.adjusted, chisqNull.adjusted=chisqNull.adjusted, NFI.adjusted=NFI.adjusted,
              NNFI.adjusted=NNFI.adjusted, CFI.adjusted=CFI.adjusted,
              norm.res=norm.res, coeff=coeff, digits=digits, 
              iterations=object$iterations, aliased=object$aliased, raw=object$raw,
              robust=robust, robust.vcov=object$robust.vcov, adj.obj=object$adj.obj)
  class(ans) <- "summary.objectiveML"
  ans
}

print.objectiveFIML <- function(x, saturated=FALSE, ...) {
  n <- x$n
  t <- x$t
  n.fix <- x$n.fix
  df <- n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2
  if (saturated){
    cat("\n Model Chisquare = ", 2*(logLik(x, saturated=TRUE) - logLik(x)), 
        "  Df = ", df, "\n\n")
  }
  else{
    cat("\n Model log-likelihood = ", logLik(x), "  Df = ", df, "\n\n")
  }
  if (!is.null(x$coef)){
    print(x$coeff)
    if (!is.na(x$iterations)) cat("\n Iterations = ", x$iterations, "\n")
    if (!is.null(x$aliased)) cat("\n Aliased parameters:", x$aliased, "\n")
  }
  invisible(x)
}
