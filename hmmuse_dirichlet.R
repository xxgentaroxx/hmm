library(tidyverse)
source("hmm_rmulti.R")

# Create data -----
B <- 1:20
k <- 1:3
K <- length(k) #l=1:K
L <- 500

a_ans <- diag(K+1)
emer_ans <- matrix(0,K,length(B),byrow=T)

a_ans[1,1] <- 0
a_ans[1,-1] <- 1/K
for(i in 1:K){
  a_ans[i+1,-1] <- MCMCpack::rdirichlet(1,sample(exp(k),K))
  emer_ans[i,] <- MCMCpack::rdirichlet(1,sample(exp(B),length(B)))
}

xdata <- create_observation_rmulti(B, k, emer_ans, a_ans, L, sample(1000:10000,L,replace=TRUE))

# Initialization -----
# Dual way
# a_init <- diag(K+1)
# a_init[1,1] <- 0
# a_init[,-1] <- 1/K

# Left-to-Right
a_init <- matrix(numeric((K+1)*(K+1)),K+1,K+1)
a_init[1,2] <- 1
for(i in 2:K){
  a_init[i,i] <- 0.5
  a_init[i,i+1] <- 0.5
}
a_init[K+1,K+1] <- 1

emer_init <- matrix(0,K,length(B),byrow=T)
for(i in 1:K){
  emer_init[i,] <- abs(rnorm(length(B)))
  emer_init[i,] <- emer_init[i,]/sum(emer_init[i,])
}

# HMM estimation -----
hmm_estim <- hmm_rmulti(xdata$x, B, k, emer_init, a_init)
