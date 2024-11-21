# DirichletMultinomial

This is a repository for the work on Dirichlet-Multinomial models.

## Simulations

### Simulations without covariates

```r
rm(list=ls())
# Load the function to simulate Dirichlet-Multinomial distributed data
source("R/simDM.R")

#================================
# Simulate data without covariates
#================================
set.seed(123)
# Simulate two individuals with the same total counts 'mu' of each individual
Y1 <- simDM(n = 2, m = 100, mu = c(0.1, 0.5, 0.4), phi = 0.3)

# Simulate two individuals with different total counts 'mu'
Y2 <- simDM(n = 2, m = c(100, 120), mu = c(0.1, 0.5, 0.4), phi = 0.3)
```

### Simulations with covariates

We can simualte Dirichlet-Multinomial distributed responses with their mean modeled by covariates $\mathbf X$ and their dispersion modeled by covariates $\mathbf Z$ as follows.

$$
\begin{aligned}
    \boldsymbol y_{i \cdot} | \boldsymbol\mu,\phi  \sim \mathcal{DM}(\boldsymbol\mu,\phi ) 
    &= \binom{m}{(y_{i1},...,y_{iq})} \frac{\prod_{j=1}^q\prod_{k=1}^{y_{ik}}[\mu_j(1-\phi) + (k-1)\phi]}{\prod_{k=1}^m [1-\phi + (k-1)\phi]}, \\
    %
    g_\mu (\mu_{ij}) &= \mathbf X_i \boldsymbol\beta_j, \\
    %
    g_\phi (\phi_{i}) &= \mathbf Z_i \boldsymbol\gamma, 
\end{aligned}
$$

where $g_\mu(\cdot)$ and $g_\phi(\cdot)$ are logistic and log link functions, respectively.

```r
set.seed(123)
n <- 100
p <- 5
q <- 3

# simulate covariates X linked to DM's mean
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
betas <- matrix(runif(p * q, -1, 1), nrow = p, ncol = q)
mu_all <- eta <- matrix(nrow = n, ncol = q)
for (j in 1:(q - 1)) {
  eta[, j] <- exp(X %*% betas[, j])
}
mu_all[, q] <- 1 / (1 + rowSums(eta[, -q]))
for (j in 1:(q - 1)) {
  mu_all[, j] <- eta[, j] / (1 + rowSums(eta[, -q]))
}

# simulate covariates Z linked to DM's dispersion
set.seed(123)
Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
gammas <- runif(p, -1, 1)
phi_all <- exp(Z %*% gammas) ## log link
phi_all[phi_all > 1] <- 0.9

# phi_all <- 1 - exp(-exp( Z %*% gammas )) ## clog-log link
# phi_all[phi_all > 0.8] <- 0.8 + runif(sum(phi_all > 0.8), -0.1, 0.1) # It's better to be smaller

# simulate outcomes Y
Y <- matrix(nrow = n, ncol = q)
for (i in 1:n) {
  Y[i, ] <- simDM(n = 1, m = 100, mu = mu_all[i, ], phi = phi_all[i])
}

head(Y)
```

```
##      [,1] [,2] [,3]
## [1,]    8   12   80
## [2,]    5   71   24
## [3,]    0   93    7
## [4,]    1   89   10
## [5,]  100    0    0
## [6,]   18   27   55
```