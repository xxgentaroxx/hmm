source("hmm.R")

# Create observation data -----
B <- 1:6
k <- 1:2
K <- length(k) #l=1:K

a_ans <- diag(K+1)
a_ans[1,1] <- 0
a_ans[1,2] <- 0.5; a_ans[1,3] <- 1-a_ans[1,2]
a_ans[2,3] <- 0.05; a_ans[2,2] <- 1-a_ans[2,3]
a_ans[3,2] <- 0.1; a_ans[3,3] <- 1-a_ans[3,2]

emer_ans <- matrix(0,K,length(B),byrow=T)
emer_ans[1,] <- rep(1/6,length(B))
emer_ans[2,] <- c(rep(1/10,5),1/2)
emer_ans <- emer_ans/rowSums(emer_ans)

xdata <- create_observation(B, k, emer_ans, a_ans, 30000)

# Initialization -----

a_init <- diag(K+1)
a_init[1,1] <- 0
a_init[1,2] <- 0.5; a_init[1,3] <- 1-a_init[1,2]
a_init[2,3] <- 0.5; a_init[2,2] <- 1-a_init[2,3]
a_init[3,2] <- 0.5; a_init[3,3] <- 1-a_init[3,2]

emer_init <- matrix(0,K,length(B),byrow=T)
emer_init[1,] <- c(5,3,2,8,1,1)/20
emer_init[2,] <- c(1,2,3,8,1,5)/20

# HMM estimation -----

hmm_estim <- hmm(xdata$x, B, k, emer_init, a_init)

mean(hmm_estim$pi_star==xdata$state)
