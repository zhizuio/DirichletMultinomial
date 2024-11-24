
# Function to simulate Dirichlet-Multinomial data

DirichletMulti <- function(n, m, mu, phi){
  # n: number of individuals
  # m: integer or vector of total number(s) of counts; if integer, every individual has the same total counts
  # mu: vector of mean/proportion parameters of all categories
  # phi: dispersion parameter (a scalar)
  
  if (length(m) == 1) {
    mm <- rep(m, n)
  }else{mm = m}
  
  if (phi >= 1 || phi <= 0) {
    stop("Argument 'phi' has to be in (0, 1)!")
  }
  
  # Calculate the Dirichlet concentration parameters
  alpha0 <- 1 / phi - 1
  alpha <- alpha0 * mu
  
  # Use 'rgamma(1, shape = alpha[j], rate = 1)' 
  # to generate n-individual Dirichlet proportions
  Gam <- matrix(0, nrow = n, ncol = length(alpha))
  for(j in 1:length(alpha)) {
    Gam[, j] <- rgamma(n, shape = alpha[j], rate = 1)
  }
  mu_ij <- Gam / rowSums(Gam)
  
  Y <- matrix(0, nrow = n, ncol = length(alpha))
  for(i in 1:n) {
    Y[i, ] <- rmultinom(1, mm[i], mu_ij[i, ])
  }
  
  return( Y )
}


# Simulate data with covariates
set.seed(123)
n <- 100
p <- 2
q <- 3

X <- array(rnorm(n * p ), dim = c(n, p))
betas <- matrix(c(rep(1,p),rep(2,p*(q-1))), nrow = p, ncol = q)
mu_all <- tmp <- matrix(nrow = n, ncol = q)
for (j in 1:(q - 1)) {
  tmp[, j] <- exp(X %*% betas[, j])
}
mu_all[, q] <- 1 / (1 + rowSums(tmp[, -q]))
for (j in 1:(q - 1)) {
  mu_all[, j] <- tmp[, j] / (1 + rowSums(tmp[, -q]))
}


Y <- matrix(nrow = n, ncol = q)
for( i in 1:n){
  Y[i, ] <- DirichletMulti(n = 1, m = 20, mu = mu_all[i, ], phi = 0.25)
}

###################################################
## various fitted quantities useful to calculate ##
## the necessary expressions                     ##
###################################################
fitfun <- function(par) {#par的前p个是beta^1,
  
  beta <- matrix(par[1:(p*(q-1))],ncol=(q-1),nrow=p,byrow = FALSE)

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
  phi <- par[-(1:(p*(q-1)))]
  
  
  
  list(
    beta = beta,
    mu_all=mu_all,
    phi = phi
  )
}

###################################
## log-likelihood function ########
###################################
loglikfun <- function(par, fit = NULL) {
  ## extract fitted quantities involved in the likelihood
  if(is.null(fit)) {
    fit <- fitfun(par)
  }
  with(fit, {
    
    phi <- 1 / (1 + exp(-phi))
    ## compute log-likelihood
    ll <- rep(0,n)
    
    for (i in 1:n) {
      C12 <- 0
      for (j in 1:q) {
        
        if(Y[i,j]==0){
          C12=C12+0
          }else{
          for (k in 1:Y[i,j]) {
            C12 = C12 + log(mu_all[i,j] * (1-phi) + (k-1)*phi)
          }
        }
      
      }
      mi <- sum(Y[i,])
      C3 <- 0
      if(mi==0){
        C3 = C3 + 0
      }else{
        for (k3 in 1:mi) {
          C3 <- C3 +  log(1-phi +(k3-1)*phi)
        }
      }
     
      ll[i] <-  C12-C3
    }
    
    -sum(ll)
    
  })
}



opt <- nlminb(start = c(rep(0.1,p*(q-1)),0.1), 
              objective = loglikfun,
              lower = c(rep(-100,p*(q-1)),0),
              upper = c(rep(100,p*(q-1)),1)
             )
(par <- opt$par)

 

#############################################
## score function (gradient) by default    ## 
## or gradient contributions (sum = FALSE) ##
#############################################
# gradfun <- function(par, sum = TRUE, fit = NULL) {
#   ## extract fitted quantities involved in the score function
#   if(is.null(fit)) {
#     fit <- fitfun(par, deriv = 1L)
#   }
#   with(fit, {
#     mu <- mu
#     phi <- phi
#     eta <- eta
#     phi_eta <- phi_eta
#     d1mu <- d1mu
#     d1phi <- d1phi
#     mustar <- mustar
#     ## compute gradient contributions
#     res <- cbind(
#       phi * (ystar - mustar) * d1mu * weights * x,
#       (mu * (ystar - mustar) + log(1-y) - digamma((1-mu)*phi) + digamma(phi)) * d1phi * weights * z
#     )
#     if(sum) colSums(res) else res
#   })
# }

#########################################
## Fisher information (default)  ########
## or covariance matrix (inverse=TRUE) ##
#########################################


