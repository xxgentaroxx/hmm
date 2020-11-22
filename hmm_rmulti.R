create_observation_rmulti <- function(B,k,emer,a,L,n){
  x <- matrix(nrow=length(B), ncol=L)
  if(length(n)==1) n <- rep(n, L)
  state <- numeric(L)
  state[1] <- sample(k,1,prob=a[1,-1])
  x[,1] <- rmultinom(1, size = n[1], prob = emer[state[1],])
  for(i in 2:L){
    state[i] <- sample(k,1,prob=a[state[i-1]+1,-1])
    x[,i] <- rmultinom(1, size = n[i], prob = emer[state[i],])
  }
  return(list(
    "x"=x,
    "state"=state
  ))
}

logsumexp <- function(pk){
  if(sum(is.nan(pk-max(pk)))>0)return(-Inf)
  log(sum(exp(pk-max(pk))))+max(pk)
}

forward_rmulti <- function(x,B,k,emer,a){
  K <- length(k) #l=i:K
  L <- ncol(x) #i=1:L
  f <- matrix(numeric((K+1)*L),K+1,L,byrow=T)
  f[1,-1] <- -Inf
  a <- log(a)
  f[-1,1] <- -Inf
  
  for(i in 2:L){
    for(l in 1:K){
      f[l+1,i] <- dmultinom(x=x[,i],prob=emer[l,],log=TRUE)+logsumexp(f[,i-1]+a[,l+1])
    }
  }
  
  return(list(
    "f"=f,
    "p"=logsumexp(f[,L])
  ))
}

backward_rmulti <- function(x,B,k,emer,a){
  K <- length(k) #l=i:K
  L <- ncol(x) #i=1:L
  
  b <- matrix(numeric((K+1)*L),K+1,L,byrow=T)
  b[-1,L] <- log(0.5)
  a <- log(a)
  
  for(i in (L-1):1){
    for(l in K:1){
      b[l+1,i] <- logsumexp(unlist(lapply(1:K,function(j)a[l+1,j+1]+
                                            dmultinom(x=x[,i+1],prob=emer[j,],log=TRUE)+b[j+1,i+1])))
    }
  }
  
  return(list(
    "b"=b,
    "p"=logsumexp(unlist(lapply(1:K,function(l)a[1,l+1]+
                                  dmultinom(x=x[,1],prob=emer[l,],log=TRUE)+b[l+1,1])))
  ))
}

bw_expectation_rmulti <- function(x,B,k,emer,a,fw,bw,p){
  K <- length(k) #l=1:K
  L <- ncol(x) #i=1:L
  a[-1,-1] <- log(a[-1,-1])
  
  E <- mapply(function(b)apply(fw+bw+rep(log(x[b,]),each=K)-rep(log(colSums(x)),each=K),1,logsumexp),B)-p
  # E <- mapply(function(b)apply((fw+bw)[,x[b,]==1],1,logsumexp),B)-p
  A <- matrix(mapply(function(from,to)logsumexp(fw[from,-L]+a[from+1,to+1]+
              mapply(function(xi)dmultinom(x=xi,prob=emer[to,],log=TRUE),
                     xi=lapply(2:L,function(i)(x[,i])))+bw[to,-1]),
                     from=rep(1:K,each=K),to=rep(1:K,K)),K,K,byrow=T)-p
  emer <- exp(E-apply(E,1,logsumexp))
  a[-1,-1] <- exp(A-apply(A,1,logsumexp))
  
  return(list("e"=emer,"a"=a))
}

baum_welch_rmulti <- function(x,B,k,emer,a,epsilon=0.01,maxloop=1000){
  forward_res <- forward_rmulti(x,B,k,emer,a)
  backward_res <- backward_rmulti(x,B,k,emer,a)
  
  fw <- forward_res$f[-1,]
  bw <- backward_res$b[-1,]
  p <- backward_res$p
  
  for(loop in 1:maxloop){
    bw_res <- bw_expectation_rmulti(x,B,k,emer,a,fw,bw,p)
    old_p <- p
    
    a <- bw_res$a
    emer <- bw_res$e
    
    forward_res <- forward_rmulti(x,B,k,emer,a)
    backward_res <- backward_rmulti(x,B,k,emer,a)
    
    fw <- forward_res$f[-1,]
    bw <- backward_res$b[-1,]
    p <- backward_res$p
    
    print(paste(loop,":",p-old_p))
    if(p-old_p<epsilon) break
  }
  
  return(list("e"=emer, "a"=a, "p"=p))
}

post_decode_rmulti <- function(x,B,k,emer,a){
  forward_res <- forward_rmulti(x,B,k,emer,a)
  backward_res <- backward_rmulti(x,B,k,emer,a)
  
  fw <- forward_res$f[-1,]
  bw <- backward_res$b[-1,]
  p <- backward_res$p
  
  return(list("pi_star"=apply(fw+bw-p, 2, which.max), "p"=p))
}

hmm_rmulti <- function(x,B,k,emer,a,epsilon=0.01,maxloop=1000){
  bw_res <- baum_welch_rmulti(x,B,k,emer,a,epsilon,maxloop)
  pd_res <- post_decode_rmulti(x,B,k,bw_res$e,bw_res$a)
  
  return(list("e"=bw_res$e, "a"=bw_res$a, 
              "pi_star"=pd_res$pi_star, "logp"=pd_res$p))
}
