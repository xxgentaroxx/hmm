source("hmm.R")
source("hmm_rmulti.R")

hmm_auto <- function(x, k, epsilon=0.01, maxloop=1000, ltr=FALSE, rmulti=TRUE, maxtry=10){
  
  B <- 1:nrow(x)
  K <- length(k)
  
  hmm_res <- list()
  hmm_res$logp <- -Inf
  
  for(iter in 1:maxtry){
    print(paste("----- Loop", iter, "-----"))
    
    emer <- matrix(0,K,length(B),byrow=T)
    for(i in 1:K){
      emer[i,] <- abs(rnorm(length(B)))
      emer[i,] <- emer[i,]/sum(emer[i,])
    }
    
    if(ltr){
      a <- matrix(numeric((K+1)*(K+1)),K+1,K+1)
      a[1,2] <- 1
      for(i in 2:K){
        a[i,i:(i+1)] <- c(MCMCpack::rdirichlet(1,1:2))
      }
      a[K+1,K+1] <- 1
    }else{
      a <- diag(K+1)
      a[1,1] <- 0
      a[,-1] <- c(MCMCpack::rdirichlet(1,1:K))
    }
    
    if(rmulti){
      hmm_try <- try(hmm_rmulti(x, B, k, emer, a, epsilon, maxloop))
      if(class(hmm_try)=="try-error") next
      else if(hmm_try$logp > hmm_res$logp) hmm_res <- hmm_try
    }else{
      hmm_try <- try(hmm(x, B, k, emer, a, epsilon, maxloop))
      if(class(hmm_try)=="try-error") next
      else if(hmm_try$logp > hmm_res$logp) hmm_res <- hmm_try
    }
  }
  
  hmm_res
}
