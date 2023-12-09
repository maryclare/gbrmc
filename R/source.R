#' @name sampler
#'
#' @description Returns posterior samples for power penalized linear or
#' logistic regression
#'
#' @export
sampler <- function(y, X, sigma.sq, tau.sq, q = 2, print.iter = FALSE,
                    means, Vars, nu = NULL, samps, type = "linear",
                    offset = rep(0, length(y)), reset = FALSE,
                    gamma.means = rep(0, ncol(X)), gamma.sds = rep(1, ncol(X)),
                    max.iter.cd = 100, z.start = rep(0, ncol(X)),
                    joint = FALSE, pappr = 0,
                    a = NULL, B = NULL) {

  p <- ncol(X)
  Vars.orig <- Vars
  Corrs.orig <- cov2cor(Vars.orig)
  Vars.ei <- eigen(Vars/2 + t(Vars/2))
  Vars.rt <- Vars.ei$vectors%*%diag(ifelse(Vars.ei$values > 0,
                                           sqrt(Vars.ei$values), 0))%*%t(Vars.ei$vectors)

  sig <- sqrt(sigma.sq)
  tau <- sqrt(tau.sq)

  if (!joint) {
    theta <- rep(pi, p)
  } else {
    theta <- pi
  }

  betas <- gammas <- matrix(nrow = samps, ncol = p)
  thetas <- matrix(nrow = samps, ncol = length(theta))

  z <- z.start

  for (i in 1:samps) {
    if (print.iter) {cat("i=", i, "\n")}

    slice <- sliceztheta(z = z, nu = nu, theta = theta,
                         y = y, X = X, sig = sig, q = q, tau = tau,
                         means = means, Vars.rt = Vars.rt, type = type,
                         offset = offset,
                         gamma.means = gamma.means, gamma.sds = gamma.sds,
                         pappr = pappr, a = a, B = B)
    z <- slice$z
    theta <- slice$theta
    gammas[i, ] <- z
    thetas[i, ] <- theta
    betas[i, ] <- gammatobeta(gamma = gammas[i, ], means = means,
                              Vars.rt = Vars.rt,
                              gamma.means = gamma.means, gamma.sds = gamma.sds)
  }

  return(list("beta" = betas, "theta" = thetas,
              "gammas" = gammas))

}

#' @name likpriprop
#'
#' @description Returns the part of the sum of the likelihood, prior,
#' and proposal density that is proportional to the regression coefficients
#' as a function of the transformed gamma parameters and the mean and variance
#' needed to transform them back,
#' beta = Vars.rt(gamma.sds*gamma + gamma.means) + mean
#' Currently supports linear and logistic regression and allows for an offset
#' in the linear component of either
#'
#' @export
likpriprop <- function(y, gamma,
                       X, sig, q, tau,
                       means, Vars.rt, nu, type,
                       offset = rep(0, length(y)),
                       gamma.means, gamma.sds, pappr,
                       a, B) {
  beta <- gammatobeta(gamma = gamma, means = means, Vars.rt = Vars.rt,
                      gamma.means = gamma.means, gamma.sds = gamma.sds)
  Xbeta <- X%*%beta + offset
  if (type == "linear") {
    llik <- -sum((y - Xbeta)^2)/(2*sig^2)
  } else if (type == "logistic") {
    probs <- 1/(1 + exp(-Xbeta))
    eps <- 10^(-14)
    probs[probs > (1 - eps)] <- 1 - eps
    probs[probs < eps] <- eps
    llik <- sum(y*log(probs) + (1 - y)*log(1 - probs))
  }
  lprior <- -(gamma(3/q)/gamma(1/q))^(q/2)*sum(abs(beta/tau)^q)
  if (is.null(nu)) {
    lpro <-  sum(gamma^2)/2
  } else {
    lpro <- sum((nu + 1)*log(1 + (gamma^2)/((nu - 2)))/2)
  }
  logscale <- (1 - pappr)*(llik + lprior + lpro)
  if (pappr != 0) {
    logscale <- logscale + pappr*(-sum((beta - a)%*%solve(B)%*%(beta - a))/2)
  }
  exp(logscale)
}

#' @name gammatobeta
#'
#' @description Transforms gammas to betas
#'
#' @export
gammatobeta <- function(gamma, means, Vars.rt,
                        gamma.means, gamma.sds) {
  Vars.rt%*%(gamma.sds*gamma + gamma.means) + means
}

#' @name sliceztheta
#'
#' @description Performs the slice sampling step
#'
#' @export
sliceztheta <- function(z, nu, theta, y, X, sig, q, tau,
                        means, Vars.rt, type,
                        offset, gamma.means, gamma.sds, pappr,
                        a, B) {
  p <- length(z)
  delta <- z

  new <- ifelse(is.null(nu), 1, sqrt((nu - 2)/nu))*rnorm(p)/ifelse(is.null(nu), 1,
                                                                 sqrt(rgamma(length(nu), nu/2, nu/2)))
  v0 <- delta*sin(theta) + new*cos(theta)
  v1 <- delta*cos(theta) - new*sin(theta)

  theta.new <- theta
  z.new <- v0*sin(theta.new) + v1*cos(theta.new)
  for (j in 1:length(theta)) {
    l <- runif(1, 0, likpriprop(y = y, gamma = z.new,
                              X = X, sig = sig, q = q, tau = tau,
                              means = means,
                              Vars.rt = Vars.rt, nu = nu, type = type,
                              offset = offset, gamma.means = gamma.means,
                              gamma.sds = gamma.sds, pappr = pappr,
                              a = a, B = B))
    a <- 0
    b <- 2*pi

    theta.new[j] <- runif(1, a, b)
    z.new <- v0*sin(theta.new) + v1*cos(theta.new)

    while (likpriprop(y = y, gamma = z.new,
                    X = X, sig = sig, q = q, tau = tau,
                    means = means, Vars.rt = Vars.rt, nu = nu, type = type,
                    offset = offset, gamma.means = gamma.means,
                    gamma.sds = gamma.sds, pappr = pappr,
                    a = a, B = B) < l) {
    if (theta.new[j] < theta[j]) {
      a <- theta.new[j]
    } else {
      b <- theta.new[j]
    }
    theta.new[j] <- runif(1, a, b)
    z.new <- v0*sin(theta.new) + v1*cos(theta.new)
    }
  }
  theta <- theta.new
  z <- z.new
  return(list("theta" = theta.new, "z" = z.new))
}

