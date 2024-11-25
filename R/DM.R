


rm(list=ls())

DM.sim.data <- function(n,m,beta,phi,Nsim=100,trace=TRUE,seed)
{
  set.seed(seed)
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
  
  p <- nrow(beta)
  q <- ncol(beta)+1
  X <- matrix(rnorm(n*p),n,p)#matrix(rnorm(n*p,sd=1/sqrt(n)),n,p)
  mu_all <- tmp <- matrix(nrow = n, ncol = q)
  for (j in 1:(q - 1)) {
    tmp[, j] <- exp(X %*% beta[, j])
  }
  mu_all[, q] <- 1 / (1 + rowSums(tmp[, -q]))
  for (j in 1:(q - 1)) {
    mu_all[, j] <- tmp[, j] / (1 + rowSums(tmp[, -q]))
  }
  
  
  
  data.sim <- as.list(numeric(Nsim))
  
  for (ii in 1:Nsim)
  {
    Y <- matrix(nrow = n, ncol = q)
    for( i in 1:n){
      Y[i, ] <- DirichletMulti(n = 1, m = 20, mu = mu_all[i, ], phi = 0.25)
    }
    data.sim[[ii]] <- list(id=ii,Y=Y)
  }
  attr(data.sim,"X") <- X
  attr(data.sim,"seed") <- seed
  data.sim
  
}

DM.sim.single <- function(datasingle,X,trace,initial_values)
{ # computes ML and BR estimates for one dataset
  
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
    # phi_eta <- as.vector(z %*% gamma + offset[[2L]])
    # phi <- phi_linkinv(phi_eta)
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
  if (trace) print(datasingle$id)
  #mle
  require(brglm2)
  Y <- datasingle$Y
  q <- ncol(Y)
  colnames(Y) <- paste("Y",1:q,sep="")
  p <- ncol(X)
  n <- nrow(X)
  
  opt <- nlminb(start = initial_values, 
                objective = loglikfun,
                lower = c(rep(-100,p*(q-1)),1e-3),
                upper = c(rep(100,p*(q-1)),1-1e-3)
  )
  
  mle <- opt$par
  # 使用 expand.grid 生成所有组合
  combinations <- expand.grid(1:p, 1:(q-1))
  names(mle) <- c(paste('beta',combinations$Var1,combinations$Var2,sep="" ),"phi")
  
  ########################
  ## multinomial regression
  ###########################
  
  library(MGLM)
  mnreg <- MGLMreg(formula = Y ~ X+0,  dist = "MN")
  mmle <- coef(mnreg)
  names(mmle) <- paste0(paste('beta',combinations$Var1,combinations$Var2,sep="" ),"MN")
  
  
  ####result
  c(mle,c(mmle))

}

DM.sim.all <- function(dataall,cores=1,trace=FALSE,initial_values)
{
  require("plyr")
  #require("brglm2")
  #require(doMC)
  #registerDoMC(cores)
  require("foreach")
  require("doParallel")
  require(Rcpp)
  require(RcppArmadillo)
  
  registerDoParallel(cores=cores)
  
  
  X <- attributes(dataall)$X
  
  out <- ldply(.data=dataall,.fun=DM.sim.single,
               X=X,trace=trace,initial_values=initial_values,.parallel=TRUE)
  out
}


##################
## various p 
## fix n=100 and q=3
##################

setwd("/Users/caizhu/Desktop/DM")
set.seed(123)
n <- 100
ps <- c(4,5)
q <- 3
m <- 20
phi <- 0.25



Nsim <- 1000

for (iii in 1:length(ps)) {
  p <- ps[iii]
  
  beta <- matrix(c(rep(-2,p*ceiling((q-1)/2)),rep(2,p*(q-1-ceiling((q-1)/2)) )), nrow = p, ncol = (q-1) )
  
  initial_values <- c(rep(0.1,p*(q-1)),0.1)
  datasim <- DM.sim.data(n=n,m=m,beta=beta,phi=phi,Nsim=Nsim,seed=123)
  
  system.time({result.DM <- DM.sim.all(dataall = datasim,cores = 5,
                                       trace = TRUE,
                                       initial_values=initial_values)})
  
  
  filename <- paste("DM_n",n,"_p",p,"_q",q,"_m",m,"_phi",phi,"_no_gamma_Nsim",Nsim,".RData",sep ="")
  save(result.DM=result.DM,datasim,file=filename)
  
}



##################
## various sample sizes n
## fix p=2 and q=3
##################

setwd("/Users/caizhu/Desktop/DM")
set.seed(123)
ns <- c(20,50,100,500,1000)
p <- 2
q <- 3
m <- 20
phi <- 0.25

beta <- matrix(c(rep(-2,p*ceiling((q-1)/2)),rep(2,p*(q-1-ceiling((q-1)/2)) )), nrow = p, ncol = (q-1) )

initial_values <- c(rep(0.1,p*(q-1)),0.1)

Nsim <- 1000

for (iii in 1:length(ns)) {
  n <- ns[iii]
  datasim <- DM.sim.data(n=n,m=m,beta=beta,phi=phi,Nsim=Nsim,seed=123)
  
  system.time({result.DM <- DM.sim.all(dataall = datasim,cores = 5,
                                       trace = TRUE,
                                       initial_values=initial_values)})
  
  
  filename <- paste("DM_n",n,"_p",p,"_q",q,"_m",m,"_phi",phi,"_no_gamma_Nsim",Nsim,".RData",sep ="")
  save(result.DM=result.DM,datasim,file=filename)
  
}


################
## result 
################

Nsim=1000
ns <- c(20,50,100,500,1000)
p <- 2
q <- 3
m <- 20
phi <- 0.25

n <- ns[3]

(filename <- paste("DM_n",n,"_p",p,"_q",q,"_m",m,"_phi",phi,"_no_gamma_Nsim",Nsim,".RData",sep ="") )
load(file=filename)

op <- par(mfrow=c(1,2))
library(wesanderson)
WAcol <- c("#F21A00","#78B7C5","#F2AD00")
combinations <- expand.grid(1:p,1:(q-1))
boxplot(result.DM[,1:(p*(q-1)+1)],col=WAcol[2],pch='.',axes=FALSE,
        xlab = "Coefficients",ylab="Estimates",
        main="DM Maximum likelihood")
box()
axis(2)
axis(1,at=1:5,labels = c(paste("Beta",combinations$Var1,combinations$Var2,sep=""),"Phi"))
segments(0.5,-2,2.5,-2,col="gold",lwd=3)
segments(2.5,2,4.5,2,col = "gold",lwd=3)
segments(4.5,0.25,5.5,0.25,col="gold",lwd=3)
legend("topleft",paste("n=",n,",","p=",p,",","q=",q))


boxplot(result.DM[,(p*(q-1)+2):ncol(result.DM)],col=WAcol[3],pch='.',axes=FALSE,
        xlab = "Coefficients",ylab="Estimates",
        main="MN Maximum likelihood")
box()
axis(2)
axis(1,at=1:4,labels = paste("Beta",combinations$Var1,combinations$Var2,sep=""))
segments(0.5,-2,2.5,-2,col="gold",lwd=3)
segments(2.5,2,4.5,2,col = "gold",lwd=3)


############################################################

windows()
par(mfrow=c(2,2))
library(wesanderson)

pchsymb <- 20
cexval <- 1
WAcol <- c("#F21A00","#78B7C5","#F2AD00")
#par(mfrow=c(2,3))
boxplot(result.DM$beta11,result.DM$beta11MN,col=c(WAcol[2],WAcol[3]),names = c("DM","MN"),main="Beta11")
abline(h=-2,col="gold",lwd=3)
legend("topright",paste("n=",n,",","p=",p,",","q=",q))

boxplot(result.DM$beta12,result.DM$beta12MN,
        col=c(WAcol[2],WAcol[3]),
        names = c("DM","MN"),
        main="Beta12")
abline(h=2,col="gold",lwd=3)

boxplot(result.DM$beta21,result.DM$beta21MN,
        col=c(WAcol[2],WAcol[3]),
        names = c("DM","MN"),
        main="Beta21")
abline(h=-2,col="gold",lwd=3)


boxplot(result.DM$beta22,result.DM$beta22MN,
        col=c(WAcol[2],WAcol[3]),
        names = c("DM","MN"),
        main="Beta22")
abline(h=2,col="gold",lwd=3)


windows()
boxplot(result.DM$phi,
        col=WAcol[2],
        names =c("DM"),
        main="Phi",ylim=c(0.1,0.3))
abline(h=0.25,col="gold",lwd=3)
legend("topright",paste("n=",n,",","p=",p,",","q=",q))
