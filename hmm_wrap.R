source("hmm.R")
source("hmm_rmulti.R")

hmm_auto <- function(x, k, epsilon=0.01, maxloop=1000, ltr=FALSE, rmulti=TRUE){
  
  B <- 1:nrow(x)
  K <- length(k)
  
  emer <- matrix(0,K,length(B),byrow=T)
  for(i in 1:K){
    emer[i,] <- abs(rnorm(length(B)))
    emer[i,] <- emer[i,]/sum(emer[i,])
  }
  
  if(ltr){
    a <- matrix(numeric((K+1)*(K+1)),K+1,K+1)
    a[1,2] <- 1
    for(i in 2:K){
      a[i,i] <- 0.5
      a[i,i+1] <- 0.5
    }
    a[K+1,K+1] <- 1
  }else{
    a <- diag(K+1)
    a[1,1] <- 0
    a[,-1] <- 1/K
  }
  
  if(rmulti) hmm_rmulti(x, B, k, emer, a, epsilon, maxloop)
  else hmm(x, B, k, emer, a, epsilon, maxloop)
}
