source("hmm_wrap.R")

# Create observation data -----
B <- 1:6
k <- 1:2
K <- length(k) #l=1:K
L <- 300

a_ans_d <- diag(K+1)
a_ans_d[1,1] <- 0
a_ans_d[1,2] <- 0.5; a_ans_d[1,3] <- 1-a_ans_d[1,2]
a_ans_d[2,3] <- 0.05; a_ans_d[2,2] <- 1-a_ans_d[2,3]
a_ans_d[3,2] <- 0.1; a_ans_d[3,3] <- 1-a_ans_d[3,2]

emer_ans_d <- matrix(0,K,length(B),byrow=T)
emer_ans_d[1,] <- rep(1/6,length(B))
emer_ans_d[2,] <- c(rep(1/10,5),1/2)
emer_ans_d <- emer_ans_d/rowSums(emer_ans_d)

xdata_d <- create_observation_rmulti(B, k, emer_ans_d, a_ans_d, L, sample(10:100,L,replace=TRUE))

# HMM estimation -----
hmm_estim_d <- hmm_auto(xdata_d$x, k)

mean(hmm_estim_d$pi_star==xdata_d$state)

# View -----

datadf_d <- as.data.frame(xdata_d$x) %>% mutate(res=row_number()) %>% 
  pivot_longer(cols=-res,names_to="trial",names_prefix="V",values_to="value") %>% 
  mutate(trial=as.numeric(trial),res=as.character(res)) %>% 
  group_by(trial) %>% mutate(relvalue=value/sum(value)) %>% ungroup
datadf2_d <- datadf_d %>% distinct(trial) %>% 
  mutate(state=factor(hmm_estim_d$pi_star)) %>% mutate(xend=lead(trial)) %>% 
  mutate(xstart=trial-0.5,xend=replace_na(xend,max(trial)+1)-0.5)

ggplot(datadf_d)+
  geom_rect(data=datadf2_d, aes(ymin=-Inf,ymax=Inf,xmin=xstart,xmax=xend,
                          fill=state,color=NULL),alpha=0.3)+
  geom_line(size=1,aes(trial,relvalue,color=res))

order_ans_d <- order(table(xdata_d$state))
order_estim_d <- order(table(hmm_estim_d$pi_star))
pitable_d <- table(xdata_d$state,hmm_estim_d$pi_star)[order_ans_d,order_estim_d]
atable_d <- data.frame(ans=c(a_ans_d[-1,-1]),j=rep(order_ans_d,length(k)),k=rep(order_ans_d,each=length(k))) %>% 
  full_join(data.frame(estim=c(hmm_estim_d$a[-1,-1]),j=rep(order_estim_d,length(k)),k=rep(order_estim_d,each=length(k))))
emertable_d <- data.frame(ans=c(emer_ans_d),res=rep(B,each=length(k)),k=rep(order_ans_d,length(B))) %>% 
  full_join(data.frame(estim=c(hmm_estim_d$e),res=rep(B,each=length(k)),k=rep(order_estim_d,length(B))))

ggplot(atable_d, aes(ans,estim,color=factor(k))) + geom_point()
ggplot(emertable_d, aes(ans,estim,color=factor(k))) + geom_point()
