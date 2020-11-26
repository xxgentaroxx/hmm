library(tidyverse)
source("hmm_wrap.R")

# Create data -----
B <- 1:10
k <- 1:3
K <- length(k) #l=1:K
L <- 300

# Dual way
a_ans <- diag(K+1)*(K-1)+1
a_ans[1,1] <- 0
a_ans[1,-1] <- 1/K
for(i in 1:K){
  a_ans[i+1,-1] <- c(MCMCpack::rdirichlet(1,exp(a_ans[i+1,-1])))
}

# Left-to-right
a_ans_ltr <- matrix(numeric((K+1)*(K+1)),K+1,K+1)
a_ans_ltr[1,2] <- 1
for(i in 2:K){
  a_ans_ltr[i,i:(i+1)] <- c(MCMCpack::rdirichlet(1,c(sample(50:80,1),1)))
}
a_ans_ltr[K+1,K+1] <- 1

emer_ans <- matrix(0,K,length(B),byrow=T)
for(i in 1:K){
  emer_ans[i,] <- MCMCpack::rdirichlet(1,sample(exp(seq(1,length(B),length=2*length(B))),length(B)))
}

xdata <- create_observation_rmulti(B, k, emer_ans, a_ans, L, sample(1000:10000,L,replace=TRUE))
xdata_ltr <- create_observation_rmulti(B, k, emer_ans, a_ans_ltr, L, sample(1000:10000,L,replace=TRUE))

# HMM estimation -----
hmm_estim <- hmm_auto(xdata$x, k, epsilon=0.01, maxloop=1000, rmulti=TRUE, maxtry=20)
hmm_estim_ltr <- hmm_auto(xdata_ltr$x, k, epsilon=0.01, maxloop=1000, ltr=TRUE, rmulti=TRUE, maxtry=20)

# View -----

# Dual
datadf <- as.data.frame(xdata$x) %>% mutate(res=row_number()) %>% 
  pivot_longer(cols=-res,names_to="trial",names_prefix="V",values_to="value") %>% 
  mutate(trial=as.numeric(trial),res=as.character(res)) %>% 
  group_by(trial) %>% mutate(relvalue=value/sum(value)) %>% ungroup
datadf2 <- datadf %>% distinct(trial) %>% 
  mutate(state=factor(hmm_estim$pi_star)) %>% mutate(xend=lead(trial)) %>% 
  mutate(xstart=trial-0.5,xend=replace_na(xend,max(trial)+1)-0.5)

datadf %>% 
  ggplot()+
  geom_rect(data=datadf2, aes(ymin=-Inf,ymax=Inf,xmin=xstart,xmax=xend,
                              fill=state,color=NULL),alpha=0.3)+
  geom_line(size=1,aes(trial,relvalue,color=res))

order_ans <- order(table(xdata$state))
order_estim <- order(table(hmm_estim$pi_star))
pitable <- table(xdata$state,hmm_estim$pi_star)[order_ans,order_estim]
atable <- data.frame(ans=c(a_ans[-1,-1]),j=rep(order_ans,length(k)),k=rep(order_ans,each=length(k))) %>% 
  full_join(data.frame(estim=c(hmm_estim$a[-1,-1]),j=rep(order_estim,length(k)),k=rep(order_estim,each=length(k))))
emertable <- data.frame(ans=c(emer_ans),res=rep(B,each=length(k)),k=rep(order_ans,length(B))) %>% 
  full_join(data.frame(estim=c(hmm_estim$e),res=rep(B,each=length(k)),k=rep(order_estim,length(B))))

ggplot(atable, aes(ans,estim,color=factor(k))) + geom_point()
ggplot(emertable, aes(ans,estim,color=factor(k))) + geom_point()

# Left-to-right
datadf_ltr <- as.data.frame(xdata_ltr$x) %>% mutate(res=row_number()) %>% 
  pivot_longer(cols=-res,names_to="trial",names_prefix="V",values_to="value") %>% 
  mutate(trial=as.numeric(trial),res=as.character(res)) %>% 
  group_by(trial) %>% mutate(relvalue=value/sum(value)) %>% ungroup
datadf2_ltr <- datadf_ltr %>% distinct(trial) %>% 
  mutate(state=factor(hmm_estim_ltr$pi_star)) %>% mutate(xend=lead(trial)) %>% 
  mutate(xstart=trial-0.5,xend=replace_na(xend,max(trial)+1)-0.5)

datadf_ltr %>% 
  ggplot()+
  geom_rect(data=datadf2_ltr, aes(ymin=-Inf,ymax=Inf,xmin=xstart,xmax=xend,
                              fill=state,color=NULL),alpha=0.3)+
  geom_line(size=1,aes(trial,relvalue,color=res))

order_ans_ltr <- order(table(xdata_ltr$state))
order_estim_ltr <- order(table(hmm_estim_ltr$pi_star))
pitable_ltr <- table(xdata_ltr$state,hmm_estim_ltr$pi_star)[order_ans_ltr,order_estim_ltr]
atable_ltr <- data.frame(ans=c(a_ans_ltr[-1,-1]),j=rep(order_ans_ltr,length(k)),k=rep(order_ans_ltr,each=length(k))) %>% 
  full_join(data.frame(estim=c(hmm_estim_ltr$a[-1,-1]),j=rep(order_estim_ltr,length(k)),k=rep(order_estim_ltr,each=length(k))))
emertable_ltr <- data.frame(ans=c(emer_ans),res=rep(B,each=length(k)),k=rep(order_ans_ltr,length(B))) %>% 
  full_join(data.frame(estim=c(hmm_estim_ltr$e),res=rep(B,each=length(k)),k=rep(order_estim_ltr,length(B))))

ggplot(atable_ltr, aes(ans,estim,color=factor(k))) + geom_point()
ggplot(emertable_ltr, aes(ans,estim,color=factor(k))) + geom_point()
