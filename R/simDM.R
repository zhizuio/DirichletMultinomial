# =================================================
# Function to simulate Dirichlet-Multinomial data
#
# author: Zhi Zhao (zhiz@uio.no)
# date: 21-11-2024
# =================================================

simDM <- function(n, m, mu, phi) {
  # n: number of individuals
  # m: integer or vector of total number(s) of counts; if integer, every individual has the same total counts
  # mu: vector of mean/proportion parameters of all categories
  # phi: dispersion parameter (a scalar)

  if (length(m) == 1) {
    mm <- rep(m, n) # make every individual the same total counts
  } else {
    mm <- m
  }

  if (phi >= 1 || phi <= 0) {
    stop("Argument 'phi' has to be in (0, 1)!")
  }

  # Calculate the Dirichlet concentration parameters
  # alpha0 <- 1 / phi - 1
  # alpha <- mu * alpha0
  alpha <- mu * (1 - phi) / phi

  # Use 'rgamma(1, shape = alpha[j], rate = 1)'
  # to generate n-individual Dirichlet proportions
  g <- matrix(0, nrow = n, ncol = length(alpha))
  for (j in 1:length(alpha)) {
    g[, j] <- rgamma(n, shape = alpha[j], rate = 1)
  }
  if (any(rowSums(g) == 0)) {
    stop("The Dirichlet concentration parameters are all close to zero!")
  }
  mu_ij <- g / rowSums(g)

  Y <- matrix(0, nrow = n, ncol = length(alpha))
  for (i in 1:n) {
    Y[i, ] <- rmultinom(1, mm[i], mu_ij[i, ])
  }

  return(Y)
}
