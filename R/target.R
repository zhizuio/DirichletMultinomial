###################################################
## various fitted quantities useful to calculate ##
## the necessary expressions                     ##
###################################################
fitfun <- function(par) { # the last element of par is 'phi'

  beta <- matrix(par[1:(p * (q - 1))], ncol = (q - 1), nrow = p, byrow = FALSE)

  # mean quantities
  mu_all <- tmp <- matrix(nrow = n, ncol = q)
  for (j in 1:(q - 1)) {
    tmp[, j] <- exp(X %*% beta[, j])
  }
  mu_all[, q] <- 1 / (1 + rowSums(tmp[, -q]))
  for (j in 1:(q - 1)) {
    mu_all[, j] <- tmp[, j] / (1 + rowSums(tmp[, -q]))
  }

  # dispersion quantities
  phi <- par[-(1:(p * (q - 1)))]

  list(
    beta = beta,
    mu_all = mu_all,
    phi = phi
  )
}

###################################
## log-likelihood function ########
###################################
loglikfun <- function(par, fit = NULL) {
  ## extract fitted quantities involved in the likelihood
  if (is.null(fit)) {
    fit <- fitfun(par)
  }
  with(fit, {
    ## compute log-likelihood
    ll <- rep(0, n)

    for (i in 1:n) {
      C12 <- 0
      for (j in 1:q) {
        if (Y[i, j] == 0) {
          C12 <- C12 + 0
        } else {
          for (k in 1:Y[i, j]) {
            C12 <- C12 + log(mu_all[i, j] * (1 - phi) + (k - 1) * phi)
          }
        }
      }
      mi <- sum(Y[i, ])
      C3 <- 0
      if (mi == 0) {
        C3 <- C3 + 0
      } else {
        for (k3 in 1:mi) {
          C3 <- C3 + log(1 - phi + (k3 - 1) * phi)
        }
      }

      ll[i] <- C12 - C3
    }

    -sum(ll)
  })
}
