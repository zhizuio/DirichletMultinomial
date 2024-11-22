# DirichletMultinomial

This is a repository for the work on Dirichlet-Multinomial models.

## Simulations

### Strategy to simulate data from the Dirichlet-Multinomial distribution

The Dirichlet-Multinomial distribution is as follows.

$$
    \boldsymbol y_{i \cdot} | \boldsymbol\mu,\phi  \sim \mathcal{DM}(\boldsymbol\mu,\phi ) 
    = \binom{m}{(y_{i1},...,y_{iq})} \frac{\prod_{j=1}^q\prod_{k=1}^{y_{ik}}[\mu_j(1-\phi) + (k-1)\phi]}{\prod_{k=1}^m [1-\phi + (k-1)\phi]}.
$$

We can use the following algorithm to simulate Dirichlet-Multinomial distributed data.

1. Set number of individuals $n$, number of total counts $m$ in every individual, proportions of $q$ categories $\boldsymbol\mu=(\mu_1,...,\mu_q)$, and dispersion parameter $\phi$.

2. Compute the Dirichlet's concentration parameters $\alpha_j=\mu_j \phi(1-\phi)$, $j=1,...,q$.

3. Use gamma distribution to generate Dirichlet's probabilities of one individual:
    * For $j$-th category in $i$-th individual, generate $c_{ij} \sim Gamma(\alpha_j, 1)$, $j=1,...,q$.
    * Generate $\mu_{ij} = c_{ij} / \sum_{l=1}^q c_{ij}$.
    * Repeat the two steps above for $i=1,...,n$ to generate data for $n$ individuals.

4. Generate each individual's counts $\boldsymbol y_{i\cdot} \sim Multinomial(m; \mu_{i1},...,\mu_{iq})$.
 

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

### Simulations with covariates linked to the mean

We can simualte Dirichlet-Multinomial distributed responses with their mean modeled by covariates $\mathbf X$.

$$
    g_\mu (\mu_{ij}) = \mathbf X_i \boldsymbol\beta_j,
$$

where $g_\mu(\cdot)$ is the logit link functions. 


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

# simulate outcomes Y
Y <- matrix(nrow = n, ncol = q)
for (i in 1:n) {
  Y[i, ] <- simDM(n = 1, m = 100, mu = mu_all[i, ], phi = 0.1)
}

head(Y)
```
```
##      [,1] [,2] [,3]
## [1,]    9   36   55
## [2,]    8   67   25
## [3,]    1   37   62
## [4,]   15   58   27
## [5,]   31   31   38
## [6,]    3   19   78
```

### Simulations with covariates linked to both mean and dispersion

We can simualte Dirichlet-Multinomial distributed responses with their mean modeled by covariates $\mathbf X$ and their dispersion modeled by covariates $\mathbf Z$ as follows.

$$
\begin{aligned}
    g_\mu (\mu_{ij}) &= \mathbf X_i \boldsymbol\beta_j, \\
    %
    g_\phi (\phi_{i}) &= \mathbf Z_i \boldsymbol\gamma, 
\end{aligned}
$$

where $g_\mu(\cdot)$ and $g_\phi(\cdot)$ are both logit link functions. 
Note that the dispersion $\phi_i \in (0,1)$ and $\phi_i$ near $1$ might cause the simulation stopped due to only zero values generated from the gamma distribution. 
Therefore, we will use a trick (from Euloge) to tune the intercept $\gamma_0$ and control the mean of $\phi_i$ as a fixed value e.g. 0.2.

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
# To be modified with logit link and control the prevalence with mean e.g. 0.2
set.seed(123)
p <- 1
#Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
Z <- matrix(rbinom(100, 0:1, prob=c(0.5, 0.5)), nrow = n) # 
gammas <- runif(p, -1, 1)

# simulate outcomes Y
Y <- matrix(nrow = n, ncol = q)
for (i in 1:n) {
  Y[i, ] <- simDM(n = 1, m = 100, mu = mu_all[i, ], phi = phi_all[i])
}

head(Y)
```

```
##      [,1] [,2] [,3]
## [1,]    0    0  100
## [2,]    0   26   74
## [3,]    0   68   32
## [4,]    2   30   68
## [5,]    4   51   45
## [6,]    0    1   99
```