linmodEst <- function(x, y)
{
  ## compute QR-decomposition of x
  qx <- qr(x)
  ## compute (x’x)^(-1) x’y
  coef <- solve.qr(qx, y)
  ## degrees of freedom and standard deviation of residuals
  df <- nrow(x)-ncol(x)
  sigma2 <- sum((y - x%*%coef)^2)/df
  ## compute sigma^2 * (x’x)^-1
  vcov <- sigma2 * chol2inv(qx$qr)
  colnames(vcov) <- rownames(vcov) <- colnames(x)
  list(coefficients = coef,
       vcov = vcov,
       sigma = sqrt(sigma2),
       df = df)
}

#' Linear Regression
#' 
#' Fit a linear regression model.
#' 
#' @aliases linmod linmod.default linmod.formula print.linmod predict.linmod
#' summary.linmod print.summary.linmod
#' @param x a numeric design matrix for the model.
#' @param y a numeric vector of responses.
#' @param formula a symbolic description of the model to be fit.
#' @param data an optional data frame containing the variables in the model.
#' @param object an object of class \code{"linmod"}, i.e., a fitted model.
#' @param \dots not used.
#' @return An object of class \code{logreg}, basically a list including
#' elements \item{coefficients}{ a named vector of coefficients } \item{vcov}{
#' covariance matrix of coefficients } \item{fitted.values}{ fitted values }
#' \item{residuals}{ residuals }
#' @author Friedrich Leisch
#' @keywords regression
#' @examples
#' 
#' data(cats, package="MASS")
#' mod1 <- linmod(Hwt~Bwt*Sex, data=cats)
#' mod1
#' summary(mod1)
#' 
#' @export
linmod <- function(x, ...) UseMethod("linmod")

#' @rdname linmod
#' @export

linmod.default <- function(x, y, ...)
{
  x <- as.matrix(x)
  y <- as.numeric(y)
  est <- linmodEst(x, y)
  est$fitted.values <- as.vector(x %*% est$coefficients)
  est$residuals <- y - est$fitted.values
  est$call <- match.call()
  class(est) <- "linmod"
  est
}

#' @rdname linmod
#' @export linmod.formula

linmod.formula <- function(formula, data=list(), ...)
{
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  est <- linmod.default(x, y, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}

#' @rdname linmod
#' @export

print.linmod <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

#' @rdname linmod
#' @export

summary.linmod <- function(object, ...)
{
  se <- sqrt(diag(object$vcov))
  tval <- coef(object) / se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               t.value = tval,
               p.value = 2*pt(-abs(tval), df=object$df))
  res <- list(call=object$call,
              coefficients=TAB)
  class(res) <- "summary.linmod"
  res
}

#' @rdname linmod
#' @export

print.summary.linmod <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
}
